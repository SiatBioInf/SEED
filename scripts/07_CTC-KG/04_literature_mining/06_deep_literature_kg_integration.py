"""
深度文献→知识图谱集成管道 (v3)
==============================
从182+篇CTC文献中提取多维度关系，与现有KG交叉验证，
创建多模态证据链。

新功能 vs v2 (llm_extract_mechanism.py):
  1. 创建Publication节点（证据溯源）
  2. 提取CancerType上下文（知道关系发生在哪个癌种）
  3. 提取CTC亚群信息（epithelial/mesenchymal/hybrid）
  4. 提取检测方法（CellSearch/微流控/RT-PCR等）
  5. 提取药物敏感/耐药信息
  6. 提取生物标志物证据（诊断/预后/预测）
  7. KG交叉验证：基因是否在scRNA数据中表达？是否在通路中？是否有已知药物？
  8. 新型边：BIOMARKER_FOR, DETECTED_BY, EPITHELIAL_MESENCHYMAL_STATE,
     DRUG_RESISTANCE_IN, DRUG_SENSITIVITY_IN, CO_OCCURS_WITH
  9. 证据链：文献边 ↔ scRNA边 ↔ 通路边 形成推理闭环

用法:
  python deep_literature_kg_integration.py [--mode full|extract|import|crossref]
"""

import json, time, sys, os, re
from collections import defaultdict

sys.path.insert(0, r'D:\科研\CTC知识图谱\v2')
sys.stdout.reconfigure(encoding='utf-8')

from openai import OpenAI
from neo4j import GraphDatabase

# ============================================================
# Configuration
# ============================================================
API_KEY = os.environ["SEED_LLM_API_KEY"]
BASE_URL = "https://uni-api.cstcloud.cn/v1"
MODEL = "deepseek-v4-flash"
NEO4J_URI = "bolt://127.0.0.1:7687"
NEO4J_AUTH = (os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"])

DATA_DIR = r"D:\科研\CTC知识图谱\v2\data\processed"
ABSTRACTS_FILE = os.path.join(DATA_DIR, "abstracts.json")
EXTRACTED_DEEP_FILE = os.path.join(DATA_DIR, "extracted_relations_deep.json")
CROSSREF_RESULTS_FILE = os.path.join(DATA_DIR, "crossref_results.json")

# ============================================================
# CTC知识图谱中的癌种名称映射
# ============================================================
CANCER_TYPE_MAP = {
    "breast cancer": "Breast Cancer", "breast": "Breast Cancer",
    "lung cancer": "Lung Cancer", "nsclc": "NSCLC", "sclc": "SCLC",
    "colorectal cancer": "Colorectal Cancer", "colon cancer": "Colorectal Cancer",
    "colorectal": "Colorectal Cancer", "crc": "Colorectal Cancer",
    "prostate cancer": "Prostate Cancer", "prostate": "Prostate Cancer",
    "pancreatic cancer": "Pancreatic Cancer", "pdac": "Pancreatic Cancer",
    "pancreatic": "Pancreatic Cancer",
    "gastric cancer": "Gastric Cancer", "stomach cancer": "Gastric Cancer",
    "gastric": "Gastric Cancer", "gc": "Gastric Cancer",
    "hepatocellular carcinoma": "Liver Cancer", "hcc": "Liver Cancer",
    "liver cancer": "Liver Cancer",
    "ovarian cancer": "Ovarian Cancer", "ovarian": "Ovarian Cancer",
    "melanoma": "Melanoma",
    "renal cell carcinoma": "Renal Cancer", "kidney cancer": "Renal Cancer",
    "rcc": "Renal Cancer",
    "bladder cancer": "Bladder Cancer", "bladder": "Bladder Cancer",
    "endometrial cancer": "Endometrial Cancer", "endometrial": "Endometrial Cancer",
    "osteosarcoma": "Osteosarcoma",
}

# ============================================================
# Neo4j连接
# ============================================================
driver = GraphDatabase.driver(NEO4J_URI, auth=NEO4J_AUTH)
client = OpenAI(api_key=API_KEY, base_url=BASE_URL)


# ============================================================
# Phase 0: 创建Publication节点
# ============================================================
def create_publication_nodes(abstracts):
    """为每篇文献创建Publication节点"""
    print("\n" + "=" * 60)
    print("Phase 0: Creating Publication nodes")
    print("=" * 60)

    created = 0
    with driver.session() as s:
        for pmid, data in abstracts.items():
            title = (data.get('title', '') or '')[:500]
            abstract = (data.get('abstract', '') or '')[:5000]
            # Extract year from PMID range (approximate) or title
            year = data.get('year', None)
            if not year:
                # Estimate year from PMID (rough heuristic)
                try:
                    pmid_int = int(pmid)
                    if pmid_int < 20000000: year = 2002
                    elif pmid_int < 23000000: year = 2010
                    elif pmid_int < 27000000: year = 2014
                    elif pmid_int < 32000000: year = 2018
                    elif pmid_int < 37000000: year = 2022
                    else: year = 2025
                except:
                    year = 2020

            try:
                s.run("""
                    MERGE (pub:Publication {pmid: $pmid})
                    SET pub.title = $title,
                        pub.abstract = $abstract,
                        pub.year = $year
                """, pmid=pmid, title=title, abstract=abstract, year=year)
                created += 1
            except Exception as e:
                print(f"  Error PMID {pmid}: {e}")

        count = s.run("MATCH (p:Publication) RETURN count(p) AS n").single()
        print(f"  Created/updated {created} Publication nodes")
        print(f"  Total Publication nodes: {count['n']}")
    return created


# ============================================================
# Phase 1: 增强LLM提取
# ============================================================
SYSTEM_PROMPT_V3 = """You are an expert CTC biologist and knowledge graph curator.

Extract ALL of the following from the PubMed abstract. Be thorough - extract every mention, not just the main finding.

Output EXACTLY this JSON structure (no markdown, no explanation):
{
  "cancer_types": ["standardized cancer name"],
  "ctc_subpopulations": [
    {"type": "epithelial|mesenchymal|hybrid|platelet_covered|cluster|stem_like|dormant", "evidence": "quote"}
  ],
  "genes": [
    {
      "symbol": "GENE_SYMBOL",
      "role": "driver|biomarker|target|resistance_marker|passenger|regulator",
      "evidence": "exact quote from abstract"
    }
  ],
  "detection_methods": [
    {"name": "CellSearch|RT-PCR|microfluidic|FACS|etc", "evidence": "quote"}
  ],
  "drugs": [
    {"name": "drug name", "relation": "resistance|sensitivity|targets", "gene": "target gene symbol", "evidence": "quote"}
  ],
  "biomarkers": [
    {"gene": "GENE", "cancer_type": "cancer name", "type": "diagnostic|prognostic|predictive|monitoring", "evidence": "quote", "direction": "favorable|unfavorable"}
  ],
  "gene_phenotype_relations": [
    {"gene": "GENE", "phenotype": "EMT_induction|MET_induction|metastasis|invasion|migration|CTC_cluster_formation|immune_evasion|drug_resistance|anoikis_resistance|apoptosis_resistance|stemness|proliferation|angiogenesis|metabolic_reprogramming|senescence|dormancy|extravasation|intravasation|platelet_adhesion|NETosis", "direction": "DRIVES|SUPPRESSES|PROMOTES|INHIBITS|ASSOCIATED_WITH", "evidence": "quote", "confidence": "high|medium|low"}
  ],
  "gene_gene_relations": [
    {"gene_a": "GENE", "gene_b": "GENE", "relation": "UPREGULATES|DOWNREGULATES|ACTIVATES|INHIBITS|PHOSPHORYLATES|INTERACTS|BINDS|COMPLEX_WITH|TRANSCRIPTIONALLY_REGULATES", "mechanism": "brief mechanism", "evidence": "quote", "confidence": "high|medium|low"}
  ],
  "gene_pathway_relations": [
    {"gene": "GENE", "pathway": "pathway name", "relation": "ACTIVATES|INHIBITS|MEMBER_OF|REGULATES", "evidence": "quote", "confidence": "high|medium|low"}
  ],
  "clinical_correlations": [
    {"finding": "description", "endpoint": "OS|PFS|DFS|ORR|metastasis_free_survival", "direction": "favorable|unfavorable", "evidence": "quote", "hazard_ratio": null}
  ],
  "mechanism_summary": "one sentence describing the core molecular mechanism"
}

Phenotype types (use EXACTLY these): EMT_induction, MET_induction, metastasis, invasion, migration, CTC_cluster_formation, immune_evasion, drug_resistance, anoikis_resistance, apoptosis_resistance, stemness, proliferation, angiogenesis, metabolic_reprogramming, senescence, dormancy, extravasation, intravasation, platelet_adhesion, NETosis, chromosomal_instability

Confidence rules:
- high: direct experimental evidence (KD/KO/inhibitor/overexpression with phenotype change)
- medium: correlation/association with statistical significance
- low: inferred, indirect, or speculative

IMPORTANT:
- Return only the JSON, no markdown code fences, no explanation before or after.
- Gene symbols must be standard HGNC symbols (uppercase for human genes).
- If a section has no findings, return empty array [].
"""


def extract_from_paper(title, abstract, pmid):
    """LLM提取一篇论文的多维关系"""
    user_msg = f"TITLE: {title}\n\nABSTRACT: {abstract[:4000]}"

    try:
        resp = client.chat.completions.create(
            model=MODEL,
            messages=[
                {"role": "system", "content": SYSTEM_PROMPT_V3},
                {"role": "user", "content": user_msg}
            ],
            temperature=0,
            max_tokens=2000,
            timeout=30  # 30 second timeout
        )

        content = resp.choices[0].message.content.strip()

        # Clean up response
        if '```json' in content:
            content = content.split('```json')[1].split('```')[0]
        elif '```' in content:
            content = content.split('```')[1].split('```')[0]
        if '{' in content:
            content = content[content.find('{'):content.rfind('}') + 1]

        result = json.loads(content)
        result['pmid'] = pmid
        result['title'] = title[:200]
        return result

    except json.JSONDecodeError as e:
        print(f"  JSON parse error: {str(e)[:80]}")
        return None
    except Exception as e:
        print(f"  API error: {str(e)[:100]}")
        return None


def run_extraction(abstracts, resume_from=None):
    """对所有论文运行增强提取"""
    print("\n" + "=" * 60)
    print("Phase 1: Enhanced LLM Extraction")
    print("=" * 60)

    results = []
    total = len(abstracts)
    success = 0
    skipped = 0

    items = list(abstracts.items())
    # Sort by PMID (desc) so newer papers processed first
    items.sort(key=lambda x: int(x[0]) if x[0].isdigit() else 0, reverse=True)

    for i, (pmid, data) in enumerate(items):
        title = data.get('title', '')[:500]
        abstract = data.get('abstract', '')

        if not abstract or len(abstract) < 50:
            skipped += 1
            continue

        # Skip if already processed (resume support)
        if resume_from and pmid in resume_from:
            continue

        result = extract_from_paper(title, abstract, pmid)

        if result:
            # Count extracted items
            n_genes = len(result.get('genes', []))
            n_phen = len(result.get('gene_phenotype_relations', []))
            n_gg = len(result.get('gene_gene_relations', []))
            n_gp = len(result.get('gene_pathway_relations', []))
            n_drugs = len(result.get('drugs', []))
            n_bio = len(result.get('biomarkers', []))
            n_ctc = len(result.get('ctc_subpopulations', []))

            results.append(result)
            success += 1

            parts = []
            if n_genes: parts.append(f'{n_genes}G')
            if n_phen: parts.append(f'{n_phen}Ph')
            if n_gg: parts.append(f'{n_gg}GG')
            if n_gp: parts.append(f'{n_gp}GP')
            if n_drugs: parts.append(f'{n_drugs}D')
            if n_bio: parts.append(f'{n_bio}Bio')
            if n_ctc: parts.append(f'{n_ctc}Sub')
            tag = '+'.join(parts) if parts else 'empty'
            print(f'[{i + 1}/{total}] PMID {pmid}: {tag}')
        else:
            print(f'[{i + 1}/{total}] PMID {pmid}: FAILED')

        # Rate limiting
        time.sleep(0.1)
        if (i + 1) % 50 == 0:
            print(f'  --- Progress: {i + 1}/{total}, {success} successful ---')
            # Save intermediate
            with open(EXTRACTED_DEEP_FILE, 'w', encoding='utf-8') as f:
                json.dump(results, f, ensure_ascii=False, indent=2)

    print(f'\nExtraction complete: {success}/{total} successful, {skipped} skipped (no abstract)')
    return results


# ============================================================
# Phase 2: KG交叉验证
# ============================================================
def cross_reference_with_kg(extracted_results):
    """将提取的基因与KG中已有数据交叉验证"""
    print("\n" + "=" * 60)
    print("Phase 2: Cross-referencing with Knowledge Graph")
    print("=" * 60)

    # Collect all unique genes from extraction
    all_genes = set()
    for r in extracted_results:
        for g in r.get('genes', []):
            symbol = g.get('symbol', '').strip().upper()
            if symbol and len(symbol) >= 2:
                all_genes.add(symbol)
        for rel in r.get('gene_phenotype_relations', []):
            symbol = rel.get('gene', '').strip().upper()
            if symbol and len(symbol) >= 2:
                all_genes.add(symbol)
        for rel in r.get('gene_gene_relations', []):
            for key in ['gene_a', 'gene_b']:
                symbol = rel.get(key, '').strip().upper()
                if symbol and len(symbol) >= 2:
                    all_genes.add(symbol)
        for rel in r.get('gene_pathway_relations', []):
            symbol = rel.get('gene', '').strip().upper()
            if symbol and len(symbol) >= 2:
                all_genes.add(symbol)

    print(f"  Unique genes extracted: {len(all_genes)}")

    # Batch query KG for each gene
    gene_kg_info = {}
    with driver.session() as s:
        for gene in all_genes:
            info = {
                'in_kg': False,
                'expressed_in_cancers': [],
                'pathways': [],
                'ppi_partners': [],
                'drugs_targeting': [],
                'prognostic_in': [],
                'top_specificity_score': 0
            }

            # Check if gene exists
            node = s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=gene).single()
            if not node:
                gene_kg_info[gene] = info
                continue

            info['in_kg'] = True

            # Get CTC expression data
            expr_results = s.run("""
                MATCH (g:Gene {symbol: $s})-[r:CTC_EXPRESSED_IN]->(c:CancerType)
                RETURN c.name AS cancer, r.specificity_score AS score, r.tier AS tier
                ORDER BY r.specificity_score DESC
            """, s=gene)
            for rec in expr_results:
                info['expressed_in_cancers'].append({
                    'cancer': rec['cancer'],
                    'score': rec['score'],
                    'tier': rec['tier']
                })
                if rec['score'] and rec['score'] > info['top_specificity_score']:
                    info['top_specificity_score'] = rec['score']

            # Get pathway memberships
            pw_results = s.run("""
                MATCH (g:Gene {symbol: $s})-[:MEMBER_OF_PATHWAY]->(p:Pathway)
                RETURN p.name AS pathway
            """, s=gene)
            info['pathways'] = [rec['pathway'] for rec in pw_results]

            # Get PPI partners
            ppi_results = s.run("""
                MATCH (g:Gene {symbol: $s})-[r:INTERACTS_WITH]-(p:Gene)
                RETURN p.symbol AS partner, r.score AS score
                ORDER BY r.score DESC LIMIT 10
            """, s=gene)
            info['ppi_partners'] = [{'partner': rec['partner'], 'score': rec['score']}
                                     for rec in ppi_results]

            # Get drugs targeting
            drug_results = s.run("""
                MATCH (d:Drug)-[:TARGETS]->(g:Gene {symbol: $s})
                RETURN d.name AS drug
            """, s=gene)
            info['drugs_targeting'] = [rec['drug'] for rec in drug_results]

            # Get prognostic data
            prog_results = s.run("""
                MATCH (g:Gene {symbol: $s})-[:PROGNOSTIC_IN]->(co:ClinicalOutcome)
                RETURN co.endpoint AS endpoint, co.hr AS hr
            """, s=gene)
            info['prognostic_in'] = [{'endpoint': rec['endpoint'], 'hr': rec['hr']}
                                      for rec in prog_results]

            gene_kg_info[gene] = info

    # Statistics
    in_kg_count = sum(1 for v in gene_kg_info.values() if v['in_kg'])
    with_expr = sum(1 for v in gene_kg_info.values() if v['expressed_in_cancers'])
    with_pathway = sum(1 for v in gene_kg_info.values() if v['pathways'])
    with_drug = sum(1 for v in gene_kg_info.values() if v['drugs_targeting'])
    with_ppi = sum(1 for v in gene_kg_info.values() if v['ppi_partners'])
    with_prog = sum(1 for v in gene_kg_info.values() if v['prognostic_in'])

    print(f"  In KG: {in_kg_count}/{len(all_genes)} ({100 * in_kg_count // max(1, len(all_genes))}%)")
    print(f"  With scRNA expression: {with_expr}")
    print(f"  With pathway membership: {with_pathway}")
    print(f"  With drug targeting: {with_drug}")
    print(f"  With PPI partners: {with_ppi}")
    print(f"  With prognostic data: {with_prog}")

    return gene_kg_info


# ============================================================
# Phase 3: 导入Neo4j（新边类型 + 交叉验证证据）
# ============================================================
def import_to_neo4j(extracted_results, gene_kg_info):
    """将提取结果导入Neo4j，创建多类型边+PUBLISHED_IN溯源"""
    print("\n" + "=" * 60)
    print("Phase 3: Importing to Neo4j")
    print("=" * 60)

    stats = defaultdict(int)

    with driver.session() as s:
        for paper in extracted_results:
            pmid = paper['pmid']
            cancer_types = paper.get('cancer_types', [])

            # --- 3a. Link to Publication node (evidence traceability) ---
            # Done at the end for each edge below

            # --- 3b. Import CTC Subpopulations ---
            for subpop in paper.get('ctc_subpopulations', []):
                stype = subpop.get('type', '').lower()
                evidence = subpop.get('evidence', '')[:300]
                if not stype: continue

                # Map to standard names
                type_map = {
                    'epithelial': 'Epithelial CTC',
                    'mesenchymal': 'Mesenchymal CTC',
                    'hybrid': 'Hybrid E/M CTC',
                    'platelet_covered': 'Platelet-Covered CTC',
                    'cluster': 'CTC Cluster',
                    'stem_like': 'Stem-like CTC',
                    'dormant': 'Dormant CTC',
                }
                display_name = type_map.get(stype, stype.replace('_', ' ').title() + ' CTC')

                s.run("MERGE (csp:CTC_Subpopulation {name: $name})", name=display_name)
                stats['ctc_subpopulations_created'] += 1

                # If genes mentioned in this paper, link them to subpopulation
                for g in paper.get('genes', []):
                    gene_symbol = g.get('symbol', '').strip().upper()
                    if not gene_symbol or len(gene_symbol) < 2: continue
                    check = s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=gene_symbol).single()
                    if not check: continue

                    s.run("""
                        MATCH (g:Gene {symbol: $gene})
                        MATCH (csp:CTC_Subpopulation {name: $subpop})
                        MERGE (g)-[r:EPITHELIAL_MESENCHYMAL_STATE {pmid: $pmid}]->(csp)
                        SET r.evidence = $evidence,
                            r.evidence_type = 'literature_llm'
                    """, gene=gene_symbol, subpop=display_name,
                         evidence=evidence, pmid=pmid)
                    stats['epithelial_mesenchymal_state_edges'] += 1

                    # Link to Publication
                    s.run("""
                        MATCH (g:Gene {symbol: $gene})
                        MATCH (csp:CTC_Subpopulation {name: $subpop})
                        MATCH (pub:Publication {pmid: $pmid})
                        MATCH (g)-[r:EPITHELIAL_MESENCHYMAL_STATE]->(csp)
                        MERGE (r)-[:PUBLISHED_IN]->(pub)
                    """, gene=gene_symbol, subpop=display_name, pmid=pmid)

            # --- 3c. Import Detection Methods ---
            for method in paper.get('detection_methods', []):
                mname = method.get('name', '').strip()
                evidence = method.get('evidence', '')[:300]
                if not mname: continue

                # Standardize method names
                method_std = {
                    'cellsearch': 'CellSearch', 'rt-pcr': 'RT-PCR', 'qrt-pcr': 'RT-PCR',
                    'qpcr': 'RT-PCR', 'microfluidic': 'Microfluidic Chip',
                    'facs': 'FACS', 'flow cytometry': 'FACS',
                    'immunomagnetic': 'Immunomagnetic Enrichment',
                    'epispot': 'EPISPOT', 'filtration': 'Microfiltration',
                    'iscope': 'iFISH/SET', 'fish': 'FISH',
                    'ddpcr': 'ddPCR', 'ngs': 'NGS', 'rna-seq': 'RNA-seq',
                    'scrnaseq': 'scRNA-seq', 'single cell rna': 'scRNA-seq',
                    'mass spectrometry': 'Mass Spectrometry',
                    'western blot': 'Western Blot', 'elisa': 'ELISA',
                }

                for key, std in method_std.items():
                    if key in mname.lower():
                        mname = std
                        break

                s.run("MERGE (dm:DetectionMethod {name: $name})", name=mname)
                s.run("""
                    MATCH (dm:DetectionMethod {name: $name})
                    MATCH (pub:Publication {pmid: $pmid})
                    MERGE (pub)-[:PUBLISHED_IN]->(dm)
                """, name=mname, pmid=pmid)
                stats['detection_methods_created'] += 1

                # Link genes to detection methods
                for g in paper.get('genes', []):
                    gene_symbol = g.get('symbol', '').strip().upper()
                    if not gene_symbol or len(gene_symbol) < 2: continue
                    check = s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=gene_symbol).single()
                    if not check: continue

                    s.run("""
                        MATCH (g:Gene {symbol: $gene})
                        MATCH (dm:DetectionMethod {name: $method})
                        MERGE (g)-[r:DETECTED_BY {pmid: $pmid}]->(dm)
                        SET r.evidence = $evidence,
                            r.evidence_type = 'literature_llm'
                    """, gene=gene_symbol, method=mname, evidence=evidence, pmid=pmid)
                    stats['detected_by_edges'] += 1

                    # Link edge to Publication
                    s.run("""
                        MATCH (g:Gene {symbol: $gene})-[r:DETECTED_BY]->(dm:DetectionMethod {name: $method})
                        MATCH (pub:Publication {pmid: $pmid})
                        MERGE (r)-[:PUBLISHED_IN]->(pub)
                    """, gene=gene_symbol, method=mname, pmid=pmid)

            # --- 3d. Import Biomarkers ---
            for bio in paper.get('biomarkers', []):
                gene_symbol = bio.get('gene', '').strip().upper()
                cancer_type_raw = bio.get('cancer_type', '').strip()
                bio_type = bio.get('type', 'prognostic')
                evidence = bio.get('evidence', '')[:300]
                direction = bio.get('direction', '')

                if not gene_symbol or len(gene_symbol) < 2: continue
                check = s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=gene_symbol).single()
                if not check: continue

                # Map cancer type to KG standard name
                cancer_std = CANCER_TYPE_MAP.get(cancer_type_raw.lower(), cancer_type_raw)

                # Try matching cancer in KG
                cancer_node = s.run("""
                    MATCH (c:CancerType)
                    WHERE toLower(c.name) CONTAINS toLower($name)
                       OR toLower($name) CONTAINS toLower(c.name)
                    RETURN c.name AS name LIMIT 1
                """, name=cancer_std).single()

                if cancer_node:
                    cancer_name = cancer_node['name']
                    conf_map = {
                        'diagnostic': 0.7,
                        'prognostic': 0.8,
                        'predictive': 0.9,
                        'monitoring': 0.7
                    }
                    conf = conf_map.get(bio_type, 0.7)

                    s.run("""
                        MATCH (g:Gene {symbol: $gene})
                        MATCH (c:CancerType {name: $cancer})
                        MERGE (g)-[r:BIOMARKER_FOR {pmid: $pmid}]->(c)
                        SET r.biomarker_type = $btype,
                            r.direction = $direction,
                            r.confidence = $conf,
                            r.evidence = $evidence,
                            r.evidence_type = 'literature_llm'
                    """, gene=gene_symbol, cancer=cancer_name, btype=bio_type,
                         direction=direction, conf=conf, evidence=evidence, pmid=pmid)
                    stats['biomarker_for_edges'] += 1

                    # Link to Publication
                    s.run("""
                        MATCH (g:Gene {symbol: $gene})-[r:BIOMARKER_FOR]->(c:CancerType {name: $cancer})
                        MATCH (pub:Publication {pmid: $pmid})
                        MERGE (r)-[:PUBLISHED_IN]->(pub)
                    """, gene=gene_symbol, cancer=cancer_name, pmid=pmid)

            # --- 3e. Import Drug relations ---
            for drug in paper.get('drugs', []):
                drug_name = drug.get('name', '').strip()
                relation = drug.get('relation', '').strip()
                target_gene = drug.get('gene', '').strip().upper()
                evidence = drug.get('evidence', '')[:300]

                if not drug_name or not target_gene: continue
                check = s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=target_gene).single()
                if not check: continue

                # Create or match Drug node
                s.run("MERGE (d:Drug {name: $name})", name=drug_name)

                if relation in ('resistance', 'confers_resistance'):
                    s.run("""
                        MATCH (g:Gene {symbol: $gene})
                        MATCH (d:Drug {name: $drug})
                        MERGE (g)-[r:DRUG_RESISTANCE_IN {pmid: $pmid}]->(d)
                        SET r.evidence = $evidence,
                            r.evidence_type = 'literature_llm'
                    """, gene=target_gene, drug=drug_name, evidence=evidence, pmid=pmid)
                    stats['drug_resistance_edges'] += 1

                    s.run("""
                        MATCH (g:Gene {symbol: $gene})-[r:DRUG_RESISTANCE_IN]->(d:Drug {name: $drug})
                        MATCH (pub:Publication {pmid: $pmid})
                        MERGE (r)-[:PUBLISHED_IN]->(pub)
                    """, gene=target_gene, drug=drug_name, pmid=pmid)

                elif relation in ('sensitivity', 'sensitizes'):
                    s.run("""
                        MATCH (g:Gene {symbol: $gene})
                        MATCH (d:Drug {name: $drug})
                        MERGE (g)-[r:DRUG_SENSITIVITY_IN {pmid: $pmid}]->(d)
                        SET r.evidence = $evidence,
                            r.evidence_type = 'literature_llm'
                    """, gene=target_gene, drug=drug_name, evidence=evidence, pmid=pmid)
                    stats['drug_sensitivity_edges'] += 1

                    s.run("""
                        MATCH (g:Gene {symbol: $gene})-[r:DRUG_SENSITIVITY_IN]->(d:Drug {name: $drug})
                        MATCH (pub:Publication {pmid: $pmid})
                        MERGE (r)-[:PUBLISHED_IN]->(pub)
                    """, gene=target_gene, drug=drug_name, pmid=pmid)

                elif relation in ('targets', 'inhibits'):
                    s.run("""
                        MATCH (g:Gene {symbol: $gene})
                        MATCH (d:Drug {name: $drug})
                        MERGE (d)-[r:TARGETS {pmid: $pmid}]->(g)
                        SET r.evidence = $evidence,
                            r.evidence_type = 'literature_llm'
                    """, gene=target_gene, drug=drug_name, evidence=evidence, pmid=pmid)
                    stats['drug_targets_edges'] += 1

            # --- 3f. Import Gene-Phenotype relations (enhanced REGULATES) ---
            for rel in paper.get('gene_phenotype_relations', []):
                gene_symbol = rel.get('gene', '').strip().upper()
                phenotype = rel.get('phenotype', '')
                direction = rel.get('direction', 'associated_with')
                evidence = rel.get('evidence', '')[:300]
                conf_map = {'high': 0.9, 'medium': 0.6, 'low': 0.3}
                conf = conf_map.get(rel.get('confidence', 'medium'), 0.5)

                if not gene_symbol or not phenotype or len(gene_symbol) < 2: continue
                check = s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=gene_symbol).single()
                if not check: continue

                display_name = phenotype.replace('_', ' ').title()
                s.run("MERGE (ph:CTC_Phenotype {name: $name})", name=display_name)

                # Add KG cross-reference data if available
                kg_data = gene_kg_info.get(gene_symbol, {})
                scRNA_validated = len(kg_data.get('expressed_in_cancers', [])) > 0
                pathway_validated = len(kg_data.get('pathways', [])) > 0

                s.run("""
                    MATCH (g:Gene {symbol: $gene})
                    MATCH (ph:CTC_Phenotype {name: $phenotype})
                    MERGE (g)-[r:REGULATES {pmid: $pmid}]->(ph)
                    SET r.direction = $direction,
                        r.mechanism = $evidence,
                        r.confidence = $conf,
                        r.evidence_type = 'literature_llm',
                        r.scrna_validated = $scrna_val,
                        r.pathway_validated = $pw_val,
                        r.pmid = $pmid
                """, gene=gene_symbol, phenotype=display_name, direction=direction,
                     evidence=evidence, conf=conf, pmid=pmid,
                     scrna_val=scRNA_validated, pw_val=pathway_validated)
                stats['regulates_edges'] += 1

                # Link to Publication
                s.run("""
                    MATCH (g:Gene {symbol: $gene})-[r:REGULATES]->(ph:CTC_Phenotype {name: $phenotype})
                    MATCH (pub:Publication {pmid: $pmid})
                    MERGE (r)-[:PUBLISHED_IN]->(pub)
                """, gene=gene_symbol, phenotype=display_name, pmid=pmid)

            # --- 3g. Import Gene-Gene relations ---
            for rel in paper.get('gene_gene_relations', []):
                gene_a = rel.get('gene_a', '').strip().upper()
                gene_b = rel.get('gene_b', '').strip().upper()
                relation = rel.get('relation', 'associated_with')
                mechanism = rel.get('mechanism', '')[:300]
                conf_map = {'high': 0.9, 'medium': 0.6, 'low': 0.3}
                conf = conf_map.get(rel.get('confidence', 'medium'), 0.5)

                if not gene_a or not gene_b or len(gene_a) < 2 or len(gene_b) < 2: continue
                check_a = s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=gene_a).single()
                check_b = s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=gene_b).single()
                if not check_a or not check_b: continue

                # Check if STRING PPI exists
                ppi_exists = s.run("""
                    MATCH (a:Gene {symbol: $a})-[r:INTERACTS_WITH]-(b:Gene {symbol: $b})
                    RETURN r.score AS score LIMIT 1
                """, a=gene_a, b=gene_b).single()
                string_score = ppi_exists['score'] if ppi_exists else None

                s.run("""
                    MATCH (a:Gene {symbol: $gene_a})
                    MATCH (b:Gene {symbol: $gene_b})
                    MERGE (a)-[r:REGULATES {pmid: $pmid}]->(b)
                    SET r.relation_type = $relation,
                        r.mechanism = $mechanism,
                        r.confidence = $conf,
                        r.evidence_type = 'literature_llm',
                        r.string_validated_score = $string_score,
                        r.pmid = $pmid
                """, gene_a=gene_a, gene_b=gene_b, relation=relation,
                     mechanism=mechanism, conf=conf, pmid=pmid,
                     string_score=string_score)
                stats['gene_gene_edges'] += 1

                # Also create CO_OCCURS_WITH for literature co-mention
                s.run("""
                    MATCH (a:Gene {symbol: $gene_a})
                    MATCH (b:Gene {symbol: $gene_b})
                    MERGE (a)-[r:CO_OCCURS_WITH {pmid: $pmid}]-(b)
                    SET r.relation_type = $relation,
                        r.evidence_type = 'literature_co_occurrence'
                """, gene_a=gene_a, gene_b=gene_b, relation=relation, pmid=pmid)
                stats['co_occurs_with_edges'] += 1

                # Link to Publication
                s.run("""
                    MATCH (a:Gene {symbol: $gene_a})-[r:REGULATES]->(b:Gene {symbol: $gene_b})
                    MATCH (pub:Publication {pmid: $pmid})
                    MERGE (r)-[:PUBLISHED_IN]->(pub)
                """, gene_a=gene_a, gene_b=gene_b, pmid=pmid)

            # --- 3h. Import Gene-Pathway relations ---
            for rel in paper.get('gene_pathway_relations', []):
                gene_symbol = rel.get('gene', '').strip().upper()
                pathway = rel.get('pathway', '').strip()
                relation = rel.get('relation', 'associated_with')
                evidence = rel.get('evidence', '')[:300]
                conf_map = {'high': 0.9, 'medium': 0.6, 'low': 0.3}
                conf = conf_map.get(rel.get('confidence', 'medium'), 0.5)

                if not gene_symbol or not pathway or len(gene_symbol) < 2: continue
                check = s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=gene_symbol).single()
                if not check: continue

                # Try match pathway in KG
                pw_result = s.run("""
                    MATCH (p:Pathway)
                    WHERE toLower(p.name) CONTAINS toLower($pathway)
                       OR toLower($pathway) CONTAINS toLower(p.name)
                    RETURN p.name AS name LIMIT 1
                """, pathway=pathway).single()

                if pw_result:
                    pw_name = pw_result['name']
                    s.run("""
                        MATCH (g:Gene {symbol: $gene})
                        MATCH (p:Pathway {name: $pw_name})
                        MERGE (g)-[r:REGULATES {pmid: $pmid}]->(p)
                        SET r.relation_type = $relation,
                            r.mechanism = $evidence,
                            r.confidence = $conf,
                            r.evidence_type = 'literature_llm',
                            r.pmid = $pmid
                    """, gene=gene_symbol, pw_name=pw_name, relation=relation,
                         evidence=evidence, conf=conf, pmid=pmid)
                    stats['gene_pathway_edges'] += 1

                    s.run("""
                        MATCH (g:Gene {symbol: $gene})-[r:REGULATES]->(p:Pathway {name: $pw_name})
                        MATCH (pub:Publication {pmid: $pmid})
                        MERGE (r)-[:PUBLISHED_IN]->(pub)
                    """, gene=gene_symbol, pw_name=pw_name, pmid=pmid)

    # --- Summary ---
    print(f"\nImport summary:")
    for key, val in sorted(stats.items()):
        print(f"  {key}: {val}")

    return stats


# ============================================================
# Phase 4: 证据链分析报告
# ============================================================
def generate_evidence_chains():
    """分析KG中的多模态证据链"""
    print("\n" + "=" * 60)
    print("Phase 4: Evidence Chain Analysis")
    print("=" * 60)

    with driver.session() as s:
        # 查找有文献+scRNA双重证据的基因
        dual_evidence = s.run("""
            MATCH (g:Gene)-[r_lit:REGULATES]->(ph:CTC_Phenotype)
            MATCH (g)-[r_expr:CTC_EXPRESSED_IN]->(c:CancerType)
            WHERE r_lit.evidence_type = 'literature_llm'
            RETURN DISTINCT g.symbol AS gene, ph.name AS phenotype,
                   r_lit.direction AS direction, c.name AS cancer,
                   r_expr.specificity_score AS specificity
            ORDER BY r_expr.specificity_score DESC LIMIT 20
        """)

        chains = []
        for rec in dual_evidence:
            chains.append({
                'gene': rec['gene'],
                'phenotype': rec['phenotype'],
                'direction': rec['direction'],
                'cancer': rec['cancer'],
                'specificity': rec['specificity'],
            })
            print(f"  {rec['gene']} → [{rec['direction']}] → {rec['phenotype']} "
                  f"(expressed in {rec['cancer']}, spec={rec['specificity']:.2f})")

        # 查找文献+STRING双重验证的基因关系
        ppi_lit = s.run("""
            MATCH (a:Gene)-[r_lit:REGULATES]->(b:Gene)
            MATCH (a)-[r_ppi:INTERACTS_WITH]-(b)
            WHERE r_lit.evidence_type = 'literature_llm'
            RETURN a.symbol AS gene_a, b.symbol AS gene_b,
                   r_lit.relation_type AS lit_rel, r_ppi.score AS ppi_score
            LIMIT 20
        """)

        ppi_chains = []
        for rec in ppi_lit:
            ppi_chains.append({
                'gene_a': rec['gene_a'],
                'gene_b': rec['gene_b'],
                'lit_rel': rec['lit_rel'],
                'ppi_score': rec['ppi_score'],
            })
            print(f"  {rec['gene_a']} → [{rec['lit_rel']}] → {rec['gene_b']} "
                  f"(STRING score={rec['ppi_score']})")

        # 查找文献+通路+scRNA三模态证据链
        triple = s.run("""
            MATCH (g:Gene)-[r_lit:REGULATES]->(ph:CTC_Phenotype)
            MATCH (g)-[r_expr:CTC_EXPRESSED_IN]->(c:CancerType)
            MATCH (g)-[:MEMBER_OF_PATHWAY]->(p:Pathway)
            WHERE r_lit.evidence_type = 'literature_llm'
            RETURN DISTINCT g.symbol AS gene, ph.name AS phenotype,
                   c.name AS cancer, p.name AS pathway
            LIMIT 20
        """)

        triple_chains = []
        for rec in triple:
            triple_chains.append({
                'gene': rec['gene'],
                'phenotype': rec['phenotype'],
                'cancer': rec['cancer'],
                'pathway': rec['pathway'],
            })
            print(f"  {rec['gene']} ∈ {rec['pathway']} → [{rec['phenotype']}] "
                  f"in {rec['cancer']} [TRIPLE EVIDENCE]")

    return {
        'dual_evidence_chains': chains,
        'ppi_literature_chains': ppi_chains,
        'triple_evidence_chains': triple_chains,
    }


# ============================================================
# Main
# ============================================================
def main(mode='full'):
    print("=" * 60)
    print("CTC-KG v3: Deep Literature Integration Pipeline")
    print("=" * 60)

    # Load abstracts
    with open(ABSTRACTS_FILE, 'r', encoding='utf-8') as f:
        abstracts = json.load(f)
    print(f"Loaded {len(abstracts)} abstracts")

    if mode in ('full', 'pubnodes'):
        create_publication_nodes(abstracts)
        if mode == 'pubnodes':
            driver.close()
            return

    if mode in ('full', 'extract'):
        # Check for resume
        resume_set = set()
        if os.path.exists(EXTRACTED_DEEP_FILE):
            with open(EXTRACTED_DEEP_FILE, 'r', encoding='utf-8') as f:
                existing = json.load(f)
            resume_set = {r['pmid'] for r in existing}
            print(f"Found {len(existing)} previously extracted papers (will skip)")

        extracted = run_extraction(abstracts, resume_from=resume_set)

        # Save final results
        if os.path.exists(EXTRACTED_DEEP_FILE):
            with open(EXTRACTED_DEEP_FILE, 'r', encoding='utf-8') as f:
                existing = json.load(f)
            existing_pmids = {r['pmid'] for r in existing}
            for r in extracted:
                if r['pmid'] not in existing_pmids:
                    existing.append(r)
            extracted = existing

        with open(EXTRACTED_DEEP_FILE, 'w', encoding='utf-8') as f:
            json.dump(extracted, f, ensure_ascii=False, indent=2)
        print(f"\nSaved {len(extracted)} extracted results to {EXTRACTED_DEEP_FILE}")
    else:
        with open(EXTRACTED_DEEP_FILE, 'r', encoding='utf-8') as f:
            extracted = json.load(f)
        print(f"Loaded {len(extracted)} previously extracted results")

    if mode in ('full', 'crossref'):
        gene_kg_info = cross_reference_with_kg(extracted)
        with open(CROSSREF_RESULTS_FILE, 'w', encoding='utf-8') as f:
            json.dump(gene_kg_info, f, ensure_ascii=False, indent=2)
    else:
        if os.path.exists(CROSSREF_RESULTS_FILE):
            with open(CROSSREF_RESULTS_FILE, 'r', encoding='utf-8') as f:
                gene_kg_info = json.load(f)
        else:
            gene_kg_info = {}

    if mode in ('full', 'import'):
        import_to_neo4j(extracted, gene_kg_info)

    if mode in ('full', 'report'):
        chains = generate_evidence_chains()
        with open(os.path.join(DATA_DIR, 'evidence_chains.json'), 'w', encoding='utf-8') as f:
            json.dump(chains, f, ensure_ascii=False, indent=2)

    # Final statistics
    print("\n" + "=" * 60)
    print("Final KG State:")
    print("=" * 60)
    with driver.session() as s:
        edge_types = ['CTC_EXPRESSED_IN', 'INTERACTS_WITH', 'MEMBER_OF_PATHWAY',
                      'ACTIVE_IN', 'TARGETS', 'PROGNOSTIC_IN', 'REGULATES',
                      'BIOMARKER_FOR', 'DETECTED_BY', 'EPITHELIAL_MESENCHYMAL_STATE',
                      'DRUG_RESISTANCE_IN', 'DRUG_SENSITIVITY_IN', 'CO_OCCURS_WITH']
        for et in edge_types:
            try:
                r = s.run(f"MATCH ()-[r:{et}]->() RETURN count(r) AS n").single()
                if r: print(f"  {et}: {r['n']}")
            except:
                pass

        node_types = ['Gene', 'CancerType', 'Pathway', 'Drug', 'ClinicalOutcome',
                      'CTC_Phenotype', 'Publication', 'CTC_Subpopulation',
                      'DetectionMethod', 'Mutation', 'miRNA']
        for nt in node_types:
            try:
                r = s.run(f"MATCH (n:{nt}) RETURN count(n) AS n").single()
                if r and r['n'] > 0: print(f"  {nt}: {r['n']}")
            except:
                pass

    driver.close()
    print("\nDone! Deep literature integration complete.")


if __name__ == '__main__':
    mode = sys.argv[2] if len(sys.argv) > 2 and sys.argv[1] == '--mode' else 'full'
    main(mode)
