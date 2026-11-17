import os
import json, sys
sys.stdout.reconfigure(encoding='utf-8')
from neo4j import GraphDatabase

d = GraphDatabase.driver('bolt://127.0.0.1:7687', auth=(os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"]))

with open(r'D:\科研\CTC知识图谱\v2\data\processed\extracted_relations.json', 'r', encoding='utf-8') as f:
    data = json.load(f)

print('Phenotypes:', len(data['phenotypes']))
print('Gene-gene:', len(data['gene_gene']))
print('Gene-pathway:', len(data['gene_pathway']))

# Show samples
for p in data['phenotypes'][:5]:
    print('  ' + str(p.get('gene')) + ' -> ' + str(p.get('phenotype')) + ' [' + str(p.get('pmid')) + ']')

imported = 0
conf_map = {'high': 0.9, 'medium': 0.6, 'low': 0.3}

with d.session() as s:
    for p in data['phenotypes']:
        gene = (p.get('gene') or '').upper()
        phenotype = (p.get('phenotype') or '').replace('_', ' ').title()
        if not gene or not phenotype:
            continue
        pmid = p.get('pmid', '')
        direction = p.get('direction', 'associated_with')
        evidence = (p.get('evidence') or '')[:200]
        conf = conf_map.get(p.get('confidence', 'medium'), 0.5)

        s.run("MERGE (ph:CTC_Phenotype {name: $name})", name=phenotype)

        check = s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=gene).single()
        if not check:
            continue

        s.run("""
            MATCH (g:Gene {symbol: $gene})
            MATCH (ph:CTC_Phenotype {name: $name})
            MERGE (g)-[r:REGULATES {pmid: $pmid}]->(ph)
            SET r.direction = $dir, r.mechanism = $ev, r.confidence = $conf,
                r.evidence_type = 'literature_llm', r.pmid = $pmid
        """, gene=gene, name=phenotype, dir=direction, ev=evidence, conf=conf, pmid=pmid)
        imported += 1

    for r in data['gene_gene']:
        a = (r.get('gene_a') or '').upper()
        b = (r.get('gene_b') or '').upper()
        if not a or not b:
            continue
        relation = r.get('relation', 'associated_with')
        mechanism = (r.get('mechanism') or '')[:200]
        conf = conf_map.get(r.get('confidence', 'medium'), 0.5)
        pmid = r.get('pmid', '')

        if not s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=a).single():
            continue
        if not s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=b).single():
            continue

        s.run("""
            MATCH (a:Gene {symbol: $a}), (b:Gene {symbol: $b})
            MERGE (a)-[r:REGULATES {pmid: $pmid}]->(b)
            SET r.direction = $rel, r.mechanism = $mech, r.confidence = $conf,
                r.evidence_type = 'literature_llm', r.pmid = $pmid
        """, a=a, b=b, rel=relation, mech=mechanism, conf=conf, pmid=pmid)
        imported += 1

    for r in data['gene_pathway']:
        gene = (r.get('gene') or '').upper()
        pathway = r.get('pathway') or ''
        if not gene or not pathway:
            continue
        pmid = r.get('pmid', '')
        relation = r.get('relation', 'associated_with')
        evidence = (r.get('evidence') or '')[:200]
        conf = conf_map.get(r.get('confidence', 'medium'), 0.5)

        if not s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=gene).single():
            continue
        pw = s.run("MATCH (p:Pathway) WHERE toLower(p.name) CONTAINS toLower($pw) RETURN p LIMIT 1", pw=pathway).single()
        if not pw:
            continue

        s.run("""
            MATCH (g:Gene {symbol: $gene}), (p:Pathway {name: $pn})
            MERGE (g)-[r:REGULATES {pmid: $pmid}]->(p)
            SET r.direction = $rel, r.mechanism = $ev, r.confidence = $conf,
                r.evidence_type = 'literature_llm', r.pmid = $pmid
        """, gene=gene, pn=pw['p'].get('name'), rel=relation, ev=evidence, conf=conf, pmid=pmid)
        imported += 1

print('Imported', imported, 'REGULATES edges')

# Verify
r = s.run("MATCH ()-[r:REGULATES]->() RETURN count(r) AS n").single()
print('Total REGULATES:', r['n'])
r = s.run("MATCH (ph:CTC_Phenotype) RETURN count(ph) AS n").single()
print('Phenotype nodes:', r['n'])

d.close()
print('Done.')
