"""
Import drug-target interactions from DGIdb via GraphQL API.
DGIdb aggregates DrugBank, TTD, PharmGKB, ChEMBL, CIViC, etc.
"""
import os
import sys, requests, time, json
sys.path.insert(0, r'D:\科研\CTC知识图谱\v2')
from neo4j import GraphDatabase

driver = GraphDatabase.driver('bolt://127.0.0.1:7687', auth=(os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"]))

# ============================================================
# Step 1: Get CTC-relevant genes
# ============================================================
print("Step 1: Getting CTC-relevant genes from KG")
with driver.session() as s:
    result = s.run("""
        MATCH (g:Gene)-[e:CTC_EXPRESSED_IN]->(c:CancerType)
        WHERE e.confidence_tier IN [1, 2]
        RETURN DISTINCT g.symbol AS symbol,
               max(e.specificity_score) AS max_spec,
               collect(DISTINCT c.name) AS cancers
        ORDER BY max_spec DESC
    """)
    genes_data = [dict(r) for r in result]

# Get all Tier1+2 genes
all_tier_genes = [g['symbol'] for g in genes_data[:2000]]
# Also add known CTC hub genes that might not have high specificity
hub_genes = ['EPCAM', 'VIM', 'CDH1', 'CDH2', 'TWIST1', 'SNAI1', 'ZEB1', 'ZEB2',
             'ERBB2', 'EGFR', 'MET', 'CD44', 'KRT19', 'ALDH1A1', 'PROM1',
             'PIK3CA', 'AKT1', 'MTOR', 'NFKB1', 'STAT3', 'CXCR4', 'CXCL12',
             'ABCC1', 'ABCG2', 'CD274', 'PDCD1', 'CTLA4', 'VEGFA', 'HIF1A',
             'AR', 'ESR1', 'BRCA1', 'BRCA2', 'TP53', 'MYC', 'KRAS', 'PTEN',
             'CCND1', 'CDK4', 'CDK6', 'PARP1', 'BCL2', 'MKI67', 'KLR4',
             'HOXC10', 'CXCL5', 'PLXNB2', 'ICAM1', 'CDH17', 'CHI3L1',
             'PLAU', 'IL6', 'TNF', 'NOTCH1', 'JAK2', 'TGFB1', 'WNT3A']
all_genes = list(set(all_tier_genes + hub_genes))
print(f"  Total genes to query DGIdb: {len(all_genes)}")

# ============================================================
# Step 2: Query DGIdb GraphQL in batches
# ============================================================
print("\nStep 2: Querying DGIdb GraphQL API")

DGIDB_GQL = 'https://dgidb.org/api/graphql'

def build_query(gene_list):
    """Build GraphQL query for a batch of genes"""
    names = '", "'.join(gene_list)
    return f"""
    {{
      genes(names: ["{names}"]) {{
        nodes {{
          name
          interactions {{
            drug {{ name approved antiNeoplastic immunotherapy }}
            interactionTypes {{ type }}
            interactionScore
            sources {{ sourceDbName }}
          }}
        }}
      }}
    }}
    """

def _process_result(data, interactions_list):
    """Extract interactions from DGIdb GraphQL response"""
    count = 0
    if 'data' not in data or not data['data'].get('genes'):
        return 0
    for gene_node in data['data']['genes']['nodes']:
        gene_name = gene_node['name']
        for inter in gene_node.get('interactions', []):
            drug = inter.get('drug', {})
            drug_name = drug.get('name', '')
            if not drug_name:
                continue
            types = [t['type'] for t in inter.get('interactionTypes', [])]
            sources = [s['sourceDbName'] for s in inter.get('sources', [])]
            score = inter.get('interactionScore', 0)
            interactions_list.append({
                'gene': gene_name,
                'drug': drug_name,
                'types': '|'.join(set(types)) if types else 'unknown',
                'score': score,
                'sources': '|'.join(sources[:5]),
                'approved': drug.get('approved', False),
                'antineoplastic': drug.get('antiNeoplastic', False),
                'immunotherapy': drug.get('immunotherapy', False),
            })
            count += 1
    return count

all_interactions = []
batch_size = 40  # Keep batch small to avoid query size limits
total_batches = (len(all_genes) - 1) // batch_size + 1

for i in range(0, len(all_genes), batch_size):
    batch = all_genes[i:i+batch_size]
    bi = i // batch_size + 1

    try:
        query = build_query(batch)
        resp = requests.post(DGIDB_GQL,
            json={'query': query},
            headers={'Content-Type': 'application/json'},
            timeout=120)

        if resp.status_code != 200:
            print(f"  Batch {bi}/{total_batches}: HTTP {resp.status_code}")
            continue

        data = resp.json()
        if 'errors' in data:
            err_msg = str(data['errors'][0].get('message', ''))[:120]
            print(f"  Batch {bi}/{total_batches}: GraphQL error - {err_msg}")
            # Try smaller batch
            if len(batch) > 10:
                for j in range(i, min(i+len(batch), len(all_genes)), 10):
                    small_batch = all_genes[j:j+10]
                    try:
                        query2 = build_query(small_batch)
                        resp2 = requests.post(DGIDB_GQL,
                            json={'query': query2},
                            headers={'Content-Type': 'application/json'},
                            timeout=120)
                        if resp2.status_code == 200:
                            data2 = resp2.json()
                            batch_count = _process_result(data2, all_interactions)
                    except:
                        pass
                    time.sleep(0.2)
            continue

        batch_count = _process_result(data, all_interactions)
        print(f"  Batch {bi}/{total_batches}: {len(batch)} genes -> {batch_count} interactions (total: {len(all_interactions)})")

    except Exception as e:
        print(f"  Batch {bi}/{total_batches}: ERROR {str(e)[:100]}")

    time.sleep(0.25)

print(f"\n  Total DGIdb interactions: {len(all_interactions)}")

# Deduplicate
seen = set()
unique_inter = []
for inter in all_interactions:
    key = (inter['gene'].upper(), inter['drug'].upper())
    if key not in seen:
        seen.add(key)
        unique_inter.append(inter)
    else:
        # Merge: keep higher score
        for existing in unique_inter:
            if (existing['gene'].upper(), existing['drug'].upper()) == key:
                existing['score'] = max(existing['score'], inter['score'])
                break
print(f"  Unique interactions: {len(unique_inter)}")

# ============================================================
# Step 3: Import to Neo4j
# ============================================================
print(f"\nStep 3: Importing {len(unique_inter)} drug-target edges to Neo4j")

with driver.session() as s:
    # Remove old DGIdb data
    s.run("MATCH ()-[r:TARGETS]->() WHERE r.source = 'DGIdb' DELETE r")
    s.run("MATCH (d:Drug) WHERE d.source = 'DGIdb' DETACH DELETE d")

    imported = 0
    skipped = 0
    for inter in unique_inter:
        try:
            result = s.run("""
                MERGE (d:Drug {name: $drug})
                  ON CREATE SET d.source = 'DGIdb',
                                d.approved = $approved,
                                d.antineoplastic = $antineoplastic,
                                d.immunotherapy = $immunotherapy
                  ON MATCH SET d.approved = CASE WHEN $approved THEN true ELSE d.approved END,
                              d.antineoplastic = CASE WHEN $antineoplastic THEN true ELSE d.antineoplastic END,
                              d.immunotherapy = CASE WHEN $immunotherapy THEN true ELSE d.immunotherapy END
                WITH d
                MATCH (g:Gene {symbol: $gene})
                MERGE (d)-[r:TARGETS]->(g)
                SET r.action = $types,
                    r.interaction_score = $score,
                    r.source = 'DGIdb',
                    r.sources_detail = $sources
            """, drug=inter['drug'], gene=inter['gene'],
                 types=inter['types'], score=inter['score'],
                 sources=inter['sources'],
                 approved=inter['approved'],
                 antineoplastic=inter['antineoplastic'],
                 immunotherapy=inter['immunotherapy'])
            imported += 1
        except Exception as e:
            skipped += 1

    print(f"  Imported: {imported}, Skipped: {skipped}")

# ============================================================
# Step 4: Verification
# ============================================================
print("\n" + "=" * 60)
print("VERIFICATION")

with driver.session() as s:
    r = s.run("MATCH (d:Drug) RETURN count(d) AS n").single()
    print(f"  Total Drug nodes: {r['n']}")

    r = s.run("MATCH ()-[r:TARGETS]->() RETURN count(r) AS n").single()
    print(f"  Total TARGETS edges: {r['n']}")

    r = s.run("""
        MATCH ()-[r:TARGETS]->()
        RETURN r.source AS source, count(r) AS n
        ORDER BY n DESC
    """)
    print("  By source:")
    for row in r:
        print(f"    {row['source']}: {row['n']}")

    # Approved antineoplastic drugs targeting CTC genes
    r = s.run("""
        MATCH (d:Drug {approved: true, antineoplastic: true})-[t:TARGETS]->(g:Gene)
        WHERE t.interaction_score > 0.3
        RETURN d.name AS drug, collect(DISTINCT g.symbol)[..5] AS targets,
               max(t.interaction_score) AS score
        ORDER BY score DESC
        LIMIT 25
    """)
    print("\n  Top approved antineoplastic drugs targeting CTC genes:")
    for row in r:
        print(f"    {row['drug']:30s} -> {', '.join(row['targets'])}")

    # Best covered CTC genes
    r = s.run("""
        MATCH (d:Drug)-[t:TARGETS]->(g:Gene)
        MATCH (g)-[e:CTC_EXPRESSED_IN]->(c:CancerType)
        WHERE e.confidence_tier = 1
        WITH g, collect(DISTINCT d.name) AS drugs, count(DISTINCT d) AS n_drugs,
             max(e.specificity_score) AS spec
        WHERE n_drugs > 5
        RETURN g.symbol AS gene, n_drugs, spec, drugs[..5] AS examples
        ORDER BY n_drugs DESC
        LIMIT 20
    """)
    print("\n  CTC genes with most drugs:")
    for row in r:
        print(f"    {row['gene']:12s} {row['n_drugs']:3d} drugs (spec={row['spec']:.2f}): {', '.join(row['examples'])}")

driver.close()
print("\nDone!")
