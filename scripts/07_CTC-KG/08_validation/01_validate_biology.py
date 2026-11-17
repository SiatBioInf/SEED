import os
import sys
sys.stdout.reconfigure(encoding='utf-8')
from neo4j import GraphDatabase

d = GraphDatabase.driver('bolt://127.0.0.1:7687', auth=(os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"]))

def query(title, cypher, **params):
    print(f'\n=== {title} ===')
    with d.session() as s:
        for r in s.run(cypher, **params):
            print('  ' + ' | '.join(f'{k}={v}' for k,v in r.items()))

with d.session() as s:
    # 1. 经典CTC标志物验证
    print('='*60)
    print('VALIDATION 1: Known CTC markers across cancers')
    print('='*60)
    markers = ['EPCAM', 'VIM', 'CDH1', 'CDH2', 'KRT19', 'PROM1', 'TWIST1', 'SNAI1']
    for gene in markers:
        r = s.run("""
            MATCH (g:Gene {symbol: $gene})-[e:CTC_EXPRESSED_IN]->(c:CancerType)
            WHERE e.confidence_tier <= 2
            RETURN c.name AS cancer, round(e.mean_percentile,3) AS pct,
                   round(e.pct_delta,3) AS delta, e.confidence_tier AS tier
            ORDER BY e.mean_percentile DESC
        """, gene=gene)
        rows = list(r)
        if rows:
            cancers_str = ', '.join(f"{r['cancer']}({r['pct']})" for r in rows[:4])
            print(f'  {gene:8s}: found in {len(rows)} cancers [{cancers_str}]')
        else:
            print(f'  {gene:8s}: NOT FOUND in any cancer (Tier1+2)')

    # 2. KRT19 (CK19) -- 最重要的上皮CTC标志物之一
    print('\n' + '='*60)
    print('VALIDATION 2: KRT19 (CK19, key epithelial CTC marker)')
    print('='*60)
    r = s.run("""
        MATCH (g:Gene {symbol: 'KRT19'})-[e:CTC_EXPRESSED_IN]->(c:CancerType)
        RETURN c.name AS cancer, round(e.mean_percentile,3) AS pct,
               round(e.pct_delta,3) AS delta, e.detection_rate AS det,
               e.confidence_tier AS tier
        ORDER BY e.mean_percentile DESC
    """)
    for row in r:
        print(f'  {row["cancer"]:40s} pct={row["pct"]:.3f} delta={row["delta"]:.3f} det={row["det"]:.2f} tier={row["tier"]}')

    # 3. 关键通路检查
    print('\n' + '='*60)
    print('VALIDATION 3: Key CTC pathways across cancers')
    print('='*60)
    key_pathways = ['EPITHELIAL_MESENCHYMAL', 'TGF_BETA', 'PI3K_AKT', 'HYPOXIA', 'GLYCOLYSIS']
    for pw in key_pathways:
        r = s.run("""
            MATCH (p:Pathway)-[a:ACTIVE_IN]->(c:CancerType)
            WHERE p.name CONTAINS $pw AND a.confidence_tier <= 2
            RETURN p.name AS pathway, c.name AS cancer,
                   round(a.pct_delta,4) AS delta, a.confidence_tier AS tier
            ORDER BY a.pct_delta DESC LIMIT 5
        """, pw=pw)
        rows = list(r)
        if rows:
            tops = ', '.join(f'{r["cancer"]}({r["delta"]:.3f})' for r in rows[:3])
            print(f'  {pw}: {tops}')

    # 4. 检查异常值：效应量极大但检出率极低的基因
    print('\n' + '='*60)
    print('VALIDATION 4: Suspicious genes (high delta but low detection)')
    print('='*60)
    r = s.run("""
        MATCH (g:Gene)-[e:CTC_EXPRESSED_IN]->(c:CancerType)
        WHERE e.confidence_tier = 1 AND e.detection_rate < 0.1
        RETURN g.symbol AS gene, c.name AS cancer, e.pct_delta AS delta,
               e.detection_rate AS det, c.n_ctc_cells AS n
        ORDER BY e.pct_delta DESC LIMIT 15
    """)
    rows = list(r)
    if rows:
        for row in rows:
            print(f'  {row["gene"]:12s} in {row["cancer"]:30s} delta={row["delta"]:.3f} det={row["det"]:.2f} n_cells={row["n"]}')
    else:
        print('  None found -- all Tier1 genes have good detection rates (good!)')

    # 5. 药物-靶点验证
    print('\n' + '='*60)
    print('VALIDATION 5: Drug targets -- are they expressed in expected cancers?')
    print('='*60)
    r = s.run("""
        MATCH (d:Drug)-[:TARGETS]->(g:Gene)-[e:CTC_EXPRESSED_IN]->(c:CancerType)
        WHERE e.confidence_tier <= 2
        RETURN d.name AS drug, g.symbol AS target, c.name AS cancer,
               round(e.specificity_score,2) AS spec
        ORDER BY d.name, e.specificity_score DESC
    """)
    for row in r:
        print(f'  {row["drug"]:20s} -> {row["target"]:8s} (expressed in {row["cancer"]:30s} spec={row["spec"]:.2f})')

d.close()
