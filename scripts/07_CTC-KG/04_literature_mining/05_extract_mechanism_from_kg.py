"""
Extract Mechanism Relations from Old KG Gene Descriptions

做什么：
  旧KG（v1）的基因description字段里已经包含了机制信息。
  比如 TWIST1 的描述写 "TWIST1 overexpression was sufficient to induce mesenchymal morphology"
  直接用关键词匹配 + 规则提取出 Gene -> CTC_Phenotype 关系。

优势：
  - 零API成本
  - 基于人工校验过的数据（旧KG经过三轮人工验证）
  - 立即可用
"""
import os

import json, re, sys
sys.path.insert(0, r'D:\科研\CTC知识图谱\v2')
sys.stdout.reconfigure(encoding='utf-8')
from neo4j import GraphDatabase
from collections import defaultdict

d = GraphDatabase.driver('bolt://127.0.0.1:7687', auth=(os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"]))

# Load old KG
with open(r'D:\科研\CTC知识图谱\neo4j_query_table_data_2026-6-29.json', 'r', encoding='utf-8') as f:
    old_kg = json.load(f)

# Build gene description map from old KG
gene_descs = {}
for item in old_kg:
    node = item.get('m') or item.get('n')
    if node and 'Gene' in node.get('labels', []):
        props = node.get('properties', {})
        name = props.get('name', '')
        desc = props.get('description', '')
        if name and desc:
            gene_descs[name] = desc

print(f'Genes with descriptions in old KG: {len(gene_descs)}')

# ============================================================
# Keyword-based phenotype extraction
# ============================================================

PHENOTYPE_PATTERNS = {
    'EMT_induction': [
        r'EMT', r'epithelial.mesenchymal.transition', r'mesenchymal.transition',
        r'promote.*EMT', r'induce.*EMT', r'drive.*EMT',
        r'upregulat.*VIM|upregulat.*vimentin|upregulat.*CDH2|upregulat.*N-cadherin',
        r'downregulat.*CDH1|downregulat.*E-cadherin',
    ],
    'MET_induction': [
        r'MET\b', r'mesenchymal.epithelial.transition',
        r'upregulat.*CDH1|upregulat.*E-cadherin',
    ],
    'metastasis': [
        r'metastas', r'distant.spread', r'invasion',
        r'promote.*metasta', r'enhance.*metasta',
    ],
    'drug_resistance': [
        r'drug.resistance', r'chemoresistance', r'chemotherapy.resistance',
        r'resistance to', r'resistant to',
        r'EGFR.TKI.resistance', r'trastuzumab.resistance',
        r'paclitaxel.resistance', r'cisplatin.resistance', r'doxorubicin.resistance',
        r'ABCC1|MRP1',
    ],
    'immune_evasion': [
        r'immune.evasion', r'immunosuppress', r'immune.escape',
        r'NK.cell', r'T.cell', r'MDSC', r'macrophage',
        r'PD-L1|CD274', r'CTLA4',
        r'immun.*microenvironment',
    ],
    'CTC_cluster_formation': [
        r'CTC.cluster', r'cluster.formation', r'homotypic',
        r'ICAM1', r'PLXNB2',
    ],
    'anoikis_resistance': [
        r'anoikis', r'suspension.survival', r'detachment.*survival',
    ],
    'apoptosis_resistance': [
        r'apoptosis', r'anti.apoptotic', r'pro.apoptotic',
        r'BCL2|BCL-XL|MCL1|survivin',
        r'caspase',
    ],
    'stemness': [
        r'stemness', r'cancer.stem.cell', r'stem.like',
        r'CD44|CD133|PROM1|ALDH', r'self.renewal',
        r'Nanog|OCT4|SOX2|KLF4',
    ],
    'proliferation': [
        r'proliferation', r'cell.cycle', r'growth',
        r'Ki.67|MKI67',
    ],
    'angiogenesis': [
        r'angiogenesis', r'VEGF', r'vascular',
        r'endothelial',
    ],
    'metabolic_reprogramming': [
        r'glycolysis', r'metabolic.reprogramming', r'Warburg',
        r'glucose.metabolism', r'lipid.metabolism',
        r'HIF1A|HK2|LDHA|PKM2',
    ],
}

# Build gene → phenotype mapping
gene_phenotypes = defaultdict(list)

for gene, desc in gene_descs.items():
    desc_lower = desc.lower()

    for phenotype, patterns in PHENOTYPE_PATTERNS.items():
        matched = False
        for pat in patterns:
            if re.search(pat, desc_lower, re.IGNORECASE):
                matched = True
                break
        if matched:
            gene_phenotypes[gene].append(phenotype)

print(f'Genes with phenotype associations: {len(gene_phenotypes)}')

# Show distribution
phenotype_counts = defaultdict(int)
for phenotypes in gene_phenotypes.values():
    for p in phenotypes:
        phenotype_counts[p] += 1
print('\nPhenotype distribution:')
for p, c in sorted(phenotype_counts.items(), key=lambda x: -x[1]):
    print(f'  {p}: {c} genes')

# ============================================================
# Import into Neo4j
# ============================================================
print('\nImporting to Neo4j...')
imported = 0

with d.session() as s:
    for gene, phenotypes in gene_phenotypes.items():
        # Check gene exists in Neo4j
        check = s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=gene).single()
        if not check:
            continue

        for phenotype in phenotypes:
            # Create phenotype node (friendly name)
            display_name = phenotype.replace('_', ' ').title()

            s.run("""
                MERGE (ph:CTC_Phenotype {name: $name})
                ON CREATE SET ph.description = $desc
            """, name=display_name, desc=f'Extracted from old KG gene descriptions')

            # Create REGULATES edge
            try:
                s.run("""
                    MATCH (g:Gene {symbol: $gene})
                    MATCH (ph:CTC_Phenotype {name: $name})
                    MERGE (g)-[r:REGULATES]->(ph)
                    SET r.direction = 'associated_with',
                        r.evidence_type = 'literature_derived',
                        r.confidence = 0.5,
                        r.source = 'old_kg_description'
                """, gene=gene, name=display_name)
                imported += 1
            except Exception as e:
                pass

print(f'Imported {imported} gene-phenotype associations')

# Verify
with d.session() as s:
    r = s.run("MATCH ()-[r:REGULATES]->() RETURN count(r) AS n").single()
    print(f'REGULATES edges: {r["n"]}')

    r = s.run("MATCH (ph:CTC_Phenotype) RETURN count(ph) AS n").single()
    print(f'CTC_Phenotype nodes: {r["n"]}')

    # Top phenotype associations
    r = s.run("""
        MATCH (ph:CTC_Phenotype)<-[r:REGULATES]-(g:Gene)
        RETURN ph.name AS phenotype, count(g) AS n_genes
        ORDER BY n_genes DESC
    """)
    print('\nPhenotypes with most genes:')
    for row in r:
        print(f'  {row["phenotype"]:25s} {row["n_genes"]} genes')

d.close()
print('Done.')
