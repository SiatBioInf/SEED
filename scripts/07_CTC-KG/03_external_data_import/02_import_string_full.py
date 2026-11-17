import os
import sys, time, requests
sys.stdout.reconfigure(encoding='utf-8')
sys.path.insert(0, r'D:\科研\CTC知识图谱\v2')
from neo4j import GraphDatabase

NEO4J = ('bolt://127.0.0.1:7687', 'neo4j', os.environ["NEO4J_PASSWORD"])
STRING_BASE = 'https://string-db.org/api/json'
SCORE_THRESHOLD = 700  # 高置信度（之前900太严了）
BATCH_SIZE = 1500      # 每批基因数

d = GraphDatabase.driver(NEO4J[0], auth=(NEO4J[1], NEO4J[2]))

# Get all genes
with d.session() as s:
    genes = [dict(r)['g.symbol'] for r in s.run('MATCH (g:Gene) WHERE g.symbol IS NOT NULL RETURN g.symbol')]
print(f'Total genes in KG: {len(genes)}')

# Split into batches
batches = [genes[i:i+BATCH_SIZE] for i in range(0, len(genes), BATCH_SIZE)]
print(f'Batches: {len(batches)} (up to {BATCH_SIZE} genes each)')

all_edges = []
total_mapped = 0

for bi, batch in enumerate(batches):
    print(f'\nBatch {bi+1}/{len(batches)} ({len(batch)} genes)...', end=' ')

    # Step 1: Map genes to STRING IDs
    try:
        resp = requests.post(f'{STRING_BASE}/get_string_ids',
            data={'identifiers': '\r'.join(batch), 'species': 9606, 'limit': 1},
            timeout=60)
        id_map = {}
        for item in resp.json():
            sid = item.get('stringId', '')
            query = item.get('queryItem', '')
            if sid:
                id_map[sid] = query
        mapped = len(id_map)
        total_mapped += mapped
        print(f'{mapped} mapped to STRING IDs', end='')
    except Exception as e:
        print(f'ERROR: {e}')
        time.sleep(2)
        continue

    time.sleep(0.8)

    # Step 2: Get PPI network
    if len(id_map) < 2:
        print(' -> skip (too few genes)')
        continue

    try:
        string_ids = list(id_map.keys())
        resp = requests.post(f'{STRING_BASE}/network',
            data={
                'identifiers': '\r'.join(string_ids),
                'species': 9606,
                'required_score': SCORE_THRESHOLD,
                'network_type': 'physical',
                'caller_identity': 'ctc_kg_v2'
            },
            timeout=120)

        batch_edges = 0
        for item in resp.json():
            # Use preferredName directly (avoids mapping through stringId)
            gene_a = item.get('preferredName_A', '')
            gene_b = item.get('preferredName_B', '')
            # score is 0-1 float in STRING API, convert to 0-1000 int
            raw_score = item.get('score', 0)
            if gene_a and gene_b and gene_a != gene_b and raw_score > 0:
                all_edges.append({
                    'gene_a': gene_a,
                    'gene_b': gene_b,
                    'score': int(raw_score * 1000),
                    'experiments': int(item.get('escore', 0)),
                    'database': int(item.get('dscore', 0)),
                    'textmining': int(item.get('tscore', 0))
                })
                batch_edges += 1

        print(f' -> {batch_edges} edges')
    except Exception as e:
        print(f'ERROR: {e}')

    time.sleep(1)

print(f'\n\nTotal: {len(all_edges)} PPI edges from {total_mapped} mapped genes')

if len(all_edges) == 0:
    print('No edges found. Check STRING API.')
    d.close()
    exit()

# Import to Neo4j
print('\nImporting to Neo4j...')
imported = 0
with d.session() as s:
    for i in range(0, len(all_edges), 500):
        batch = all_edges[i:i+500]
        s.run("""
            UNWIND $edges AS edge
            MATCH (a:Gene {symbol: edge.gene_a})
            MATCH (b:Gene {symbol: edge.gene_b})
            MERGE (a)-[r:INTERACTS_WITH]->(b)
            SET r.score = edge.score,
                r.experiments = edge.experiments,
                r.database = edge.database,
                r.textmining = edge.textmining,
                r.source = 'STRING_v12.5'
        """, edges=batch)
        imported += len(batch)
    print(f'Imported {imported} edges')

    # Delete test edges from broken previous runs
    s.run("MATCH ()-[r:INTERACTS_WITH]->() WHERE r.source = 'test' DELETE r")

# Verify
with d.session() as s:
    r = s.run('MATCH ()-[r:INTERACTS_WITH]->() RETURN count(r) AS n, avg(r.score) AS avg_score, max(r.score) AS max_score').single()
    print('\nVerification: PPI edges=' + str(r['n']) + ', avg_score=' + str(r['avg_score']) + ', max_score=' + str(r['max_score']))

    # Test: VIM interactions
    r = s.run("""
        MATCH (g:Gene {symbol: 'VIM'})-[r:INTERACTS_WITH]-(g2:Gene)
        RETURN g2.symbol AS interactor, r.score AS score
        ORDER BY r.score DESC LIMIT 10
    """)
    print('VIM interactors:')
    for row in r:
        print(f'  {row["interactor"]:10s} score={row["score"]}')

d.close()
print('Done.')
