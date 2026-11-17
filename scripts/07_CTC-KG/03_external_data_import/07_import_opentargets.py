"""
Import gene-disease/phenotype associations from Open Targets Platform.

Open Targets (https://platform.opentargets.org) 整合了：
  - GWAS Catalog
  - ClinVar
  - Gene2Phenotype
  - 文献挖掘
  - 药物靶点数据

每个基因-疾病关联都有 evidence score (0-1)，比手写的description可靠得多。

用法：提取关键CTC基因与癌症相关表型的关联，导入Neo4j作为 REGULATES 边。
"""
import os

import sys, requests, time
sys.path.insert(0, r'D:\科研\CTC知识图谱\v2')
sys.stdout.reconfigure(encoding='utf-8')
from neo4j import GraphDatabase

d = GraphDatabase.driver('bolt://127.0.0.1:7687', auth=(os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"]))

# CTC-relevant phenotypes/traits from Open Targets (EFO IDs)
CTC_PHENOTYPES = {
    'neoplasm': 'EFO_0000616',
    'carcinoma': 'EFO_0000237',
    'metastasis': 'EFO_0009713',
    'cell_migration': 'GO_0016477',
    'epithelial_to_mesenchymal_transition': 'GO_0001837',
    'drug_resistance': 'EFO_0002618',
    'angiogenesis': 'GO_0001525',
    'apoptotic_process': 'GO_0006915',
    'cell_proliferation': 'GO_0008283',
    'immune_response': 'GO_0006955',
}

# CTC-relevant disease terms
CTC_DISEASES = [
    ('breast_carcinoma', 'Orphanet_180257'),
    ('lung_carcinoma', 'EFO_0001071'),
    ('colorectal_carcinoma', 'EFO_0005848'),
    ('prostate_carcinoma', 'EFO_0001663'),
    ('pancreatic_carcinoma', 'EFO_0002618'),
    ('melanoma', 'EFO_0000756'),
    ('liver_carcinoma', 'EFO_0002642'),
    ('gastric_carcinoma', 'EFO_0003909'),
    ('ovarian_carcinoma', 'EFO_0003896'),
    ('head_and_neck_carcinoma', 'EFO_0003867'),
]

# Get genes with high specificity (Tier 1) from our KG
with d.session() as s:
    result = s.run("""
        MATCH (g:Gene)-[e:CTC_EXPRESSED_IN]->(c:CancerType)
        WHERE e.confidence_tier = 1
        RETURN DISTINCT g.symbol AS symbol
        LIMIT 200
    """)
    top_genes = [r['symbol'] for r in result]
    print(f'Top Tier1 genes to query: {len(top_genes)}')

# Query Open Targets for each gene-disease pair
# Use their GraphQL API
api_url = 'https://api.platform.opentargets.org/api/v4/graphql'

# Batch genes into groups of 50 to keep queries manageable
associations = []

for i in range(0, min(len(top_genes), 100), 20):
    batch = top_genes[i:i+20]
    gene_query = '", "'.join(batch)
    disease_query = '"), ("'.join(f'{name}' for name, _ in CTC_DISEASES)

    query = """
    query GeneDiseaseAssociations {
      diseases(efoIds: ["Orphanet_180257", "EFO_0001071", "EFO_0000756"]) {
        id
        name
        associatedTargets(rows: 500) {
          rows {
            target {
              id
              approvedSymbol
            }
            score
            datatypeScores {
              id
              score
            }
          }
        }
      }
    }
    """

    # Simplified: query one disease at a time
    for disease_name, disease_id in CTC_DISEASES[:5]:  # Limit to 5 diseases for speed
        simple_query = f"""
        {{
          disease(efoId: "{disease_id}") {{
            id
            name
            associatedTargets(rows: 500) {{
              rows {{
                target {{
                  approvedSymbol
                }}
                score
              }}
            }}
          }}
        }}
        """

        try:
            resp = requests.post(api_url, json={'query': simple_query}, timeout=60)
            if resp.status_code == 200:
                data = resp.json()
                disease_data = data.get('data', {}).get('disease', {})
                disease_name_ot = disease_data.get('name', disease_name)

                rows = disease_data.get('associatedTargets', {}).get('rows', [])
                for row in rows:
                    symbol = row.get('target', {}).get('approvedSymbol', '')
                    score = row.get('score', 0)

                    if symbol and score > 0.3 and symbol in top_genes:
                        associations.append({
                            'gene': symbol,
                            'disease': disease_name_ot,
                            'score': score
                        })

                print(f'  {disease_name}: {len(rows)} targets, {len([r for r in rows if r.get("target",{}).get("approvedSymbol","") in top_genes])} match our genes')

        except Exception as e:
            print(f'  {disease_name}: ERROR - {e}')

        time.sleep(0.5)

print(f'\nTotal associations found: {len(associations)}')

# Map disease names to CTC phenotypes
DISEASE_TO_PHENOTYPE = {
    'breast carcinoma': 'Metastasis',
    'breast cancer': 'Metastasis',
    'lung carcinoma': 'Metastasis',
    'lung cancer': 'Metastasis',
    'colorectal carcinoma': 'Metastasis',
    'prostate carcinoma': 'Metastasis',
    'pancreatic carcinoma': 'Metastasis',
    'melanoma': 'Metastasis',
    'liver carcinoma': 'Metastasis',
    'gastric carcinoma': 'Metastasis',
    'ovarian carcinoma': 'Metastasis',
    'head and neck carcinoma': 'Metastasis',
}

# Import to Neo4j as REGULATES edges to CTC_Phenotype
if associations:
    imported = 0
    with d.session() as s:
        for assoc in associations:
            gene = assoc['gene']
            disease = assoc['disease']
            score = assoc['score']

            # Map to phenotype
            phenotype = DISEASE_TO_PHENOTYPE.get(
                disease.lower(),
                disease.replace(' carcinoma', '').title()
            )

            # Create phenotype node and edge
            s.run("""
                MERGE (ph:CTC_Phenotype {name: $name})
            """, name=phenotype)

            try:
                s.run("""
                    MATCH (g:Gene {symbol: $gene})
                    MATCH (ph:CTC_Phenotype {name: $phenotype})
                    MERGE (g)-[r:REGULATES]->(ph)
                    SET r.direction = 'associated_with',
                        r.confidence = $score,
                        r.evidence_type = 'genetic_association',
                        r.source = 'OpenTargets'
                """, gene=gene, phenotype=phenotype, score=score)
                imported += 1
            except:
                pass

    print(f'Imported {imported} REGULATES edges from Open Targets')

# Verify
with d.session() as s:
    r = s.run("MATCH ()-[r:REGULATES]->() RETURN count(r) AS n").single()
    print(f'\nREGULATES edges: {r["n"]}')

    r = s.run("MATCH (ph:CTC_Phenotype) RETURN count(ph) AS n").single()
    print(f'CTC_Phenotype nodes: {r["n"]}')

d.close()
print('Done.')
