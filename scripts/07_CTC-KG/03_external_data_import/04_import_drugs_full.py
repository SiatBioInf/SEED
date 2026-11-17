"""
Full drug-target import via DGIdb + OpenTargets APIs.
No download required. Free APIs.

DGIdb aggregates from 30+ sources: DrugBank, TTD, PharmGKB, ChemblInteractions, etc.
OpenTargets provides gene-disease association scores.
"""
import os
import sys, requests, time, json
sys.path.insert(0, r'D:\科研\CTC知识图谱\v2')
from neo4j import GraphDatabase

NEO4J_URI = "bolt://127.0.0.1:7687"
NEO4J_AUTH = (os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"])
driver = GraphDatabase.driver(NEO4J_URI, auth=NEO4J_AUTH)

# ============================================================
# Step 1: Get CTC-relevant genes from KG (Tier 1 + Tier 2)
# ============================================================
print("=" * 60)
print("Step 1: Getting CTC-relevant genes from KG")
with driver.session() as s:
    # Get genes with high-confidence CTC expression
    result = s.run("""
        MATCH (g:Gene)-[e:CTC_EXPRESSED_IN]->(c:CancerType)
        WHERE e.confidence_tier IN [1, 2]
        RETURN DISTINCT g.symbol AS symbol,
               max(e.specificity_score) AS max_spec,
               collect(DISTINCT c.name) AS cancers
        ORDER BY max_spec DESC
    """)
    genes_data = [dict(r) for r in result]
    print(f"  Tier1+2 genes: {len(genes_data)}")

# Split into: priority genes (top 500 by specificity) + all others
priority_genes = [g['symbol'] for g in genes_data[:500]]
all_tier_genes = [g['symbol'] for g in genes_data[:2000]]
print(f"  Priority genes: {len(priority_genes)} (top specificity)")
print(f"  Total to query: {len(all_tier_genes)}")

# ============================================================
# Step 2: DGIdb - Drug-Gene Interactions
# ============================================================
print("\n" + "=" * 60)
print("Step 2: Querying DGIdb for drug-gene interactions")

DGIDB_API = "https://dgidb.org/api/v2"

def query_dgidb_batch(genes):
    """Query DGIdb for a batch of genes. Returns list of drug-gene interactions."""
    interactions = []
    gene_list = ','.join(genes)

    try:
        # DGIdb interactions endpoint
        resp = requests.post(
            f"{DGIDB_API}/interactions",
            json={"genes": gene_list},
            headers={"Content-Type": "application/json"},
            timeout=120
        )
        if resp.status_code == 200:
            data = resp.json()
            for match in data.get('matchedTerms', []):
                gene_name = match.get('geneName', '')
                for interaction in match.get('interactions', []):
                    drug_name = interaction.get('drugName', '')
                    interaction_type = interaction.get('interactionType', '')
                    sources = interaction.get('sources', [])
                    score = interaction.get('interactionScore', 0)

                    if drug_name and gene_name:
                        interactions.append({
                            'gene': gene_name,
                            'drug': drug_name,
                            'type': interaction_type,
                            'score': score,
                            'sources': '|'.join(sources[:3])  # top 3 sources
                        })
        else:
            print(f"  DGIdb API error: {resp.status_code}")
    except Exception as e:
        print(f"  DGIdb error: {e}")

    return interactions

# Query DGIdb in batches of 100 genes
all_dgidb = []
batch_size = 100
for i in range(0, len(all_tier_genes), batch_size):
    batch = all_tier_genes[i:i+batch_size]
    interactions = query_dgidb_batch(batch)
    all_dgidb.extend(interactions)
    print(f"  Batch {i//batch_size + 1}/{(len(all_tier_genes)-1)//batch_size + 1}: "
          f"{len(batch)} genes -> {len(interactions)} interactions (total: {len(all_dgidb)})")
    time.sleep(0.3)

print(f"\n  Total DGIdb interactions: {len(all_dgidb)}")

# Deduplicate
seen = set()
unique_dgidb = []
for inter in all_dgidb:
    key = (inter['gene'], inter['drug'])
    if key not in seen:
        seen.add(key)
        unique_dgidb.append(inter)
print(f"  Unique interactions: {len(unique_dgidb)}")

# ============================================================
# Step 3: OpenTargets - Gene-Disease Associations
# ============================================================
print("\n" + "=" * 60)
print("Step 3: Querying OpenTargets for gene-disease associations")

OT_API = "https://api.platform.opentargets.org/api/v4/graphql"

# CTC-relevant diseases
DISEASES = [
    ("breast carcinoma", "EFO_0003869"),
    ("lung carcinoma", "EFO_0001071"),
    ("colorectal carcinoma", "EFO_0005848"),
    ("prostate carcinoma", "EFO_0001663"),
    ("pancreatic carcinoma", "EFO_0002618"),
    ("melanoma", "EFO_0000756"),
    ("liver carcinoma", "EFO_0002642"),
    ("gastric carcinoma", "EFO_0003909"),
    ("ovarian carcinoma", "EFO_0003896"),
    ("head and neck carcinoma", "EFO_0003867"),
]

ot_associations = []

for disease_name, efo_id in DISEASES:
    query = """
    query targetAssociations($efoId: String!) {
      disease(efoId: $efoId) {
        id
        name
        associatedTargets(page: {size: 500}) {
          rows {
            target {
              approvedSymbol
            }
            score
          }
        }
      }
    }
    """

    try:
        resp = requests.post(OT_API,
            json={'query': query, 'variables': {'efoId': efo_id}},
            timeout=60)

        if resp.status_code == 200:
            data = resp.json()
            disease_data = data.get('data', {}).get('disease', {})
            if not disease_data:
                print(f"  {disease_name}: no data")
                continue

            rows = disease_data.get('associatedTargets', {}).get('rows', [])
            matched = 0
            for row in rows:
                symbol = row.get('target', {}).get('approvedSymbol', '')
                score = row.get('score', 0)

                if symbol and score > 0.2 and symbol in all_tier_genes:
                    ot_associations.append({
                        'gene': symbol,
                        'disease': disease_name,
                        'score': score
                    })
                    matched += 1

            print(f"  {disease_name}: {len(rows)} targets, {matched} match CTC genes")
        else:
            print(f"  {disease_name}: HTTP {resp.status_code}")

    except Exception as e:
        print(f"  {disease_name}: ERROR - {e}")

    time.sleep(0.3)

print(f"\n  Total OpenTargets associations: {len(ot_associations)}")

# ============================================================
# Step 4: Import to Neo4j
# ============================================================
print("\n" + "=" * 60)
print("Step 4: Importing to Neo4j")

with driver.session() as s:
    # Clean old DGIdb/OpenTargets data
    s.run("MATCH ()-[r:TARGETS]->() WHERE r.source IN ['DGIdb', 'OpenTargets'] DELETE r")
    s.run("MATCH (d:Drug) WHERE d.source IN ['DGIdb', 'OpenTargets'] DETACH DELETE d")

    # Import DGIdb drug-target edges
    dgidb_imported = 0
    for inter in unique_dgidb:
        try:
            s.run("""
                MERGE (d:Drug {name: $drug})
                  ON CREATE SET d.source = 'DGIdb'
                WITH d
                MATCH (g:Gene {symbol: $gene})
                MERGE (d)-[r:TARGETS]->(g)
                SET r.action = $type,
                    r.source = 'DGIdb',
                    r.interaction_score = $score,
                    r.sources_detail = $sources
            """, drug=inter['drug'], gene=inter['gene'],
                 type=inter['type'], score=inter['score'],
                 sources=inter['sources'])
            dgidb_imported += 1
        except Exception as e:
            pass

    print(f"  DGIdb TARGETS edges imported: {dgidb_imported}")

    # Import OpenTargets gene-disease associations as REGULATES edges
    ot_imported = 0
    for assoc in ot_associations:
        try:
            # Map disease name to phenotype
            disease = assoc['disease']
            phenotype = disease.replace(' carcinoma', '').title()
            # Standardize
            pheno_map = {
                'Breast': 'Metastasis',
                'Lung': 'Metastasis',
                'Colorectal': 'Metastasis',
                'Prostate': 'Metastasis',
                'Pancreatic': 'Metastasis',
                'Melanoma': 'Metastasis',
                'Liver': 'Metastasis',
                'Gastric': 'Metastasis',
                'Ovarian': 'Metastasis',
                'Head And Neck': 'Metastasis',
            }
            phenotype = pheno_map.get(phenotype, phenotype)

            s.run("""
                MERGE (ph:CTC_Phenotype {name: $name})
                WITH ph
                MATCH (g:Gene {symbol: $gene})
                MERGE (g)-[r:REGULATES]->(ph)
                SET r.direction = 'associated_with',
                    r.confidence = $score,
                    r.evidence_type = 'genetic_association',
                    r.source = 'OpenTargets'
            """, name=phenotype, gene=assoc['gene'], score=assoc['score'])
            ot_imported += 1
        except:
            pass

    print(f"  OpenTargets REGULATES edges imported: {ot_imported}")

# ============================================================
# Step 5: Verification
# ============================================================
print("\n" + "=" * 60)
print("Step 5: Verification")

with driver.session() as s:
    # Count drugs
    r = s.run("MATCH (d:Drug) RETURN count(d) AS n").single()
    print(f"  Total Drug nodes: {r['n']}")

    # Count TARGETS edges
    r = s.run("MATCH ()-[r:TARGETS]->() RETURN count(r) AS n").single()
    print(f"  Total TARGETS edges: {r['n']}")

    # Count by source
    r = s.run("""
        MATCH ()-[r:TARGETS]->()
        RETURN r.source AS source, count(r) AS n
        ORDER BY n DESC
    """)
    print("  TARGETS by source:")
    for row in r:
        print(f"    {row['source']}: {row['n']}")

    # Count REGULATES
    r = s.run("MATCH ()-[r:REGULATES]->() RETURN count(r) AS n").single()
    print(f"  Total REGULATES edges: {r['n']}")

    # Test query: drugs targeting breast cancer CTC genes
    print("\n  Test: Top drugs for Breast Cancer CTC genes (Tier 1):")
    r = s.run("""
        MATCH (g:Gene)-[e:CTC_EXPRESSED_IN]->(c:CancerType {name: 'Breast Cancer'})
        WHERE e.confidence_tier = 1
        MATCH (d:Drug)-[t:TARGETS]->(g)
        RETURN g.symbol AS gene, e.specificity_score AS spec,
               collect(DISTINCT d.name)[..3] AS drugs
        ORDER BY spec DESC
        LIMIT 10
    """)
    for row in r:
        drugs_str = ', '.join(row['drugs'])
        print(f"    {row['gene']:12s} (spec={row['spec']:.2f}): {drugs_str}")

driver.close()
print("\nDone!")
