import os
import sys
sys.stdout.reconfigure(encoding='utf-8')
from neo4j import GraphDatabase

d = GraphDatabase.driver('bolt://127.0.0.1:7687', auth=(os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"]))

with d.session() as s:
    # 获取所有癌种
    cancers = [r['name'] for r in s.run('MATCH (c:CancerType) RETURN c.name AS name ORDER BY name')]

    for cn in cancers:
        print('=' * 60)
        print('  ' + cn)
        print('=' * 60)

        # Tier 1+2 基因数
        r = s.run("""
            MATCH (g:Gene)-[e:CTC_EXPRESSED_IN]->(c:CancerType {name: $cn})
            WHERE e.confidence_tier <= 2
            RETURN count(g) AS n
        """, cn=cn).single()
        tier12_count = r['n']

        # 总基因数
        r = s.run("""
            MATCH (g:Gene)-[e:CTC_EXPRESSED_IN]->(c:CancerType {name: $cn})
            RETURN count(g) AS n
        """, cn=cn).single()
        total_count = r['n']

        print(f'  Tier1+2 genes: {tier12_count}  |  Total: {total_count}')

        # Top 10 by specificity_score
        r = s.run("""
            MATCH (g:Gene)-[e:CTC_EXPRESSED_IN]->(c:CancerType {name: $cn})
            WHERE e.confidence_tier <= 2
            RETURN g.symbol AS gene, e.pct_delta AS delta,
                   round(e.cohens_d, 2) AS d, round(e.specificity_score, 2) AS spec
            ORDER BY e.specificity_score DESC LIMIT 10
        """, cn=cn)
        top_genes = list(r)
        if top_genes:
            print(f'  Top genes:')
            for row in top_genes:
                print(f'    {row["gene"]:12s}  delta={row["delta"]:.3f}  d={row["d"]:.1f}  spec={row["spec"]:.2f}')

        # 靶向药
        r = s.run("""
            MATCH (g:Gene)-[e:CTC_EXPRESSED_IN]->(c:CancerType {name: $cn})
            WHERE e.confidence_tier <= 2
            MATCH (d:Drug)-[:TARGETS]->(g)
            RETURN g.symbol AS gene, e.specificity_score AS spec,
                   collect(DISTINCT d.name) AS drugs
            ORDER BY e.specificity_score DESC LIMIT 10
        """, cn=cn)
        drugs = list(r)
        if drugs:
            print(f'  Druggable genes:')
            for row in drugs:
                drug_list = ', '.join(row['drugs'])
                print(f'    {row["gene"]:12s} (spec={row["spec"]:.2f}): {drug_list}')

        # 最活跃的通路
        r = s.run("""
            MATCH (p:Pathway)-[a:ACTIVE_IN]->(c:CancerType {name: $cn})
            WHERE a.confidence_tier <= 2
            RETURN p.name AS pathway, round(a.pct_delta, 4) AS delta
            ORDER BY a.pct_delta DESC LIMIT 8
        """, cn=cn)
        pathways = list(r)
        if pathways:
            print(f'  Top pathways:')
            for row in pathways:
                pname = row['pathway'].replace('HALLMARK_','')
                print(f'    {pname:45s}  delta={row["delta"]:.4f}')

        print()

d.close()
