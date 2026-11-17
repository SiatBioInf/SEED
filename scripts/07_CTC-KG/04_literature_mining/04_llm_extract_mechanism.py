"""
LLM提取文献机制关系 → 导入Neo4j

用CSTCloud API从182篇CTC文献摘要中提取：
  - 基因→基因 调控关系 (UPREGULATES, INHIBITS, etc.)
  - 基因→通路 关系
  - 基因→表型 关系

输出直接导入Neo4j作为 REGULATES 边。
"""
import os

import json, time, sys
sys.path.insert(0, r'D:\科研\CTC知识图谱\v2')
sys.stdout.reconfigure(encoding='utf-8')
from openai import OpenAI
from neo4j import GraphDatabase

# API配置
client = OpenAI(
    api_key=os.environ["SEED_LLM_API_KEY"],
    base_url="https://uni-api.cstcloud.cn/v1"
)

# Neo4j
d = GraphDatabase.driver('bolt://127.0.0.1:7687', auth=(os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"]))

# 读取已下载的摘要
with open(r'D:\科研\CTC知识图谱\v2\data\processed\abstracts.json', 'r', encoding='utf-8') as f:
    abstracts = json.load(f)
print(f'Loaded {len(abstracts)} abstracts')

# LLM提取Prompt
SYSTEM_PROMPT = """You are an expert in CTC (Circulating Tumor Cell) biology and cancer molecular mechanisms.

Extract MECHANISTIC relationships from the PubMed abstract below.
Rules:
1. Only extract relationships with clear CAUSAL or REGULATORY evidence.
   - "TWIST1 knockdown reduced invasion by 60%" = valid
   - "TWIST1 was detected by RT-PCR" = skip (no mechanism)
2. Identify the direction: a gene UPREGULATES, DOWNREGULATES, ACTIVATES, INHIBITS, DRIVES, or SUPPRESSES something.
3. Output ONLY valid JSON. No explanations.

Output format:
{
  "genes_mentioned": ["GENE1", "GENE2"],
  "pathways_mentioned": ["pathway name"],
  "phenotypes": [
    {"gene": "GENE1", "phenotype": "EMT_induction", "direction": "DRIVES",
     "evidence": "exact quote from abstract", "confidence": "high"}
  ],
  "gene_gene_relations": [
    {"gene_a": "GENE1", "gene_b": "GENE2", "relation": "UPREGULATES",
     "mechanism": "brief mechanism", "evidence": "exact quote", "confidence": "high"}
  ],
  "gene_pathway_relations": [
    {"gene": "GENE1", "pathway": "Wnt signaling", "relation": "ACTIVATES",
     "evidence": "exact quote", "confidence": "high"}
  ]
}

Phenotype types (use EXACTLY these): EMT_induction, MET_induction, metastasis, invasion, migration, CTC_cluster_formation, immune_evasion, drug_resistance, anoikis_resistance, apoptosis_resistance, stemness, proliferation, angiogenesis, metabolic_reprogramming

Confidence: high (direct experimental: KD/KO/inhibitor), medium (correlation/association), low (inferred/indirect)

If no mechanistic relationships, return: {"genes_mentioned": [], "pathways_mentioned": [], "phenotypes": [], "gene_gene_relations": [], "gene_pathway_relations": []}
"""

all_phenotypes = []
all_gene_gene = []
all_gene_pathway = []
total = len(abstracts)
success = 0
failed = 0

# Pre-filter: only process papers likely to contain mechanism info
MECHANISM_KEYWORDS = [
    'knockdown', 'overexpression', 'knockout', 'silencing',
    'activates', 'inhibits', 'upregulates', 'downregulates',
    'phosphorylates', 'regulates', 'mediates', 'suppresses',
    'promotes', 'enhances', 'reduces', 'induces', 'drives',
    'resistance', 'sensitive', 'signaling', 'pathway',
    'interacts', 'binds', 'complex', 'translocation',
    'EMT', 'epithelial-mesenchymal', 'MET',
    'migration', 'invasion', 'metastasis',
    'apoptosis', 'survival', 'proliferation',
    'stemness', 'differentiation',
]
SKIP_KEYWORDS = [
    'detected by', 'was detected', 'were detected', 'measured by',
    'expression of', 'expression level', 'expressed in',
    'CellSearch', 'enrichment', 'isolation',
]

filtered = {}
skipped = 0
for pmid, data in abstracts.items():
    abstract_lower = (data.get('abstract', '') + ' ' + data.get('title', '')).lower()

    # Skip detection-only papers
    has_skip = any(kw in abstract_lower for kw in SKIP_KEYWORDS)
    has_kw = any(kw in abstract_lower for kw in MECHANISM_KEYWORDS)

    if has_kw and not has_skip:
        filtered[pmid] = data
    else:
        skipped += 1

test_items = list(filtered.items())
print(f'Pre-filter: {len(test_items)} papers with mechanism content (skipped {skipped} detection-only papers)')

for i, (pmid, data) in enumerate(test_items):
    title = data.get('title', '')[:500]
    abstract = data.get('abstract', '')[:3000]

    if not abstract:
        continue

    user_msg = f"TITLE: {title}\n\nABSTRACT: {abstract}"

    try:
        resp = client.chat.completions.create(
            model="S1-Base-Ultra",
            messages=[
                {"role": "system", "content": SYSTEM_PROMPT},
                {"role": "user", "content": user_msg}
            ],
            temperature=0,
            max_tokens=1000
        )

        content = resp.choices[0].message.content.strip()
        # Remove markdown code fences if present
        if '```json' in content:
            content = content.split('```json')[1].split('```')[0]
        elif content.startswith('```'):
            content = content.split('```')[1].split('```')[0]
        # Remove any non-JSON prefix/suffix
        if '{' in content:
            content = content[content.find('{'):content.rfind('}')+1]

        result = json.loads(content)

        n_phen = len(result.get('phenotypes', []))
        n_gg = len(result.get('gene_gene_relations', []))
        n_gp = len(result.get('gene_pathway_relations', []))

        if n_phen + n_gg + n_gp > 0:
            for p in result.get('phenotypes', []):
                p['pmid'] = pmid
                all_phenotypes.append(p)
            for r in result.get('gene_gene_relations', []):
                r['pmid'] = pmid
                all_gene_gene.append(r)
            for r in result.get('gene_pathway_relations', []):
                r['pmid'] = pmid
                all_gene_pathway.append(r)

        success += 1
        tag = f'{n_phen}P + {n_gg}GG + {n_gp}GP' if n_phen + n_gg + n_gp > 0 else 'no relations'
        print(f'[{i+1}/{total}] PMID {pmid}: {tag}')

    except json.JSONDecodeError:
        failed += 1
        print(f'[{i+1}/{total}] PMID {pmid}: JSON parse error')
    except Exception as e:
        failed += 1
        print(f'[{i+1}/{total}] PMID {pmid}: ERROR - {str(e)[:80]}')

    # Rate limiting
    if i % 10 == 0:
        print(f'  Progress: {i+1}/{total}', flush=True)
    time.sleep(0.1)

print(f'\nDone. Success: {success}/{total}, Failed: {failed}')
print(f'Extracted: {len(all_phenotypes)} phenotypes, {len(all_gene_gene)} gene-gene, {len(all_gene_pathway)} gene-pathway')

# Save extracted relations
with open(r'D:\科研\CTC知识图谱\v2\data\processed\extracted_relations.json', 'w', encoding='utf-8') as f:
    json.dump({
        'phenotypes': all_phenotypes,
        'gene_gene': all_gene_gene,
        'gene_pathway': all_gene_pathway
    }, f, ensure_ascii=False, indent=2)

# ============================================================
# Import to Neo4j
# ============================================================
print('\nImporting to Neo4j...')
imported = 0

with d.session() as s:
    # Import gene-phenotype relations
    for p in all_phenotypes:
        gene = p['gene'].upper()
        phenotype = p['phenotype']
        direction = p.get('direction', 'associated_with')
        evidence = p.get('evidence', '')[:200]
        conf_map = {'high': 0.9, 'medium': 0.6, 'low': 0.3}
        conf = conf_map.get(p.get('confidence', 'medium'), 0.5)
        pmid = p['pmid']

        display = phenotype.replace('_', ' ').title()

        s.run("MERGE (ph:CTC_Phenotype {name: $name})", name=display)

        check = s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=gene).single()
        if not check:
            continue

        s.run("""
            MATCH (g:Gene {symbol: $gene})
            MATCH (ph:CTC_Phenotype {name: $phenotype})
            MERGE (g)-[r:REGULATES {pmid: $pmid}]->(ph)
            SET r.direction = $direction,
                r.mechanism = $evidence,
                r.confidence = $conf,
                r.evidence_type = 'literature_llm',
                r.pmid = $pmid
        """, gene=gene, phenotype=display, direction=direction,
             evidence=evidence, conf=conf, pmid=pmid)
        imported += 1

    # Import gene-gene relations
    for r in all_gene_gene:
        gene_a = r['gene_a'].upper()
        gene_b = r['gene_b'].upper()
        relation = r.get('relation', 'associated_with')
        mechanism = r.get('mechanism', '')[:200]
        conf_map = {'high': 0.9, 'medium': 0.6, 'low': 0.3}
        conf = conf_map.get(r.get('confidence', 'medium'), 0.5)
        pmid = r['pmid']

        check_a = s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=gene_a).single()
        check_b = s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=gene_b).single()
        if not check_a or not check_b:
            continue

        s.run("""
            MATCH (a:Gene {symbol: $gene_a})
            MATCH (b:Gene {symbol: $gene_b})
            MERGE (a)-[r:REGULATES {pmid: $pmid}]->(b)
            SET r.direction = $relation,
                r.mechanism = $mechanism,
                r.confidence = $conf,
                r.evidence_type = 'literature_llm',
                r.pmid = $pmid
        """, gene_a=gene_a, gene_b=gene_b, relation=relation,
             mechanism=mechanism, conf=conf, pmid=pmid)
        imported += 1

    # Import gene-pathway relations
    for r in all_gene_pathway:
        gene = r['gene'].upper()
        pathway = r['pathway']
        relation = r.get('relation', 'associated_with')
        evidence = r.get('evidence', '')[:200]
        conf_map = {'high': 0.9, 'medium': 0.6, 'low': 0.3}
        conf = conf_map.get(r.get('confidence', 'medium'), 0.5)
        pmid = r['pmid']

        check = s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=gene).single()
        if not check:
            continue

        # Match pathway by name (fuzzy)
        result = s.run("""
            MATCH (p:Pathway)
            WHERE toLower(p.name) CONTAINS toLower($pathway)
               OR toLower($pathway) CONTAINS toLower(p.name)
            RETURN p LIMIT 1
        """, pathway=pathway)
        pw = result.single()
        if not pw:
            continue

        s.run("""
            MATCH (g:Gene {symbol: $gene})
            MATCH (p:Pathway {name: $pw_name})
            MERGE (g)-[r:REGULATES {pmid: $pmid}]->(p)
            SET r.direction = $relation,
                r.mechanism = $evidence,
                r.confidence = $conf,
                r.evidence_type = 'literature_llm',
                r.pmid = $pmid
        """, gene=gene, pw_name=pw['p'].get('name'), relation=relation,
             evidence=evidence, conf=conf, pmid=pmid)
        imported += 1

print(f'Imported {imported} REGULATES edges')

# Verify
with d.session() as s:
    r = s.run("MATCH ()-[r:REGULATES]->() RETURN count(r) AS n").single()
    print(f'Total REGULATES edges: {r["n"]}')

    r = s.run("MATCH (ph:CTC_Phenotype) RETURN count(ph) AS n").single()
    print(f'CTC_Phenotype nodes: {r["n"]}')

d.close()
print('Done.')
