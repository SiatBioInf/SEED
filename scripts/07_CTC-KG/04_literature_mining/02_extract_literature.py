"""
Literature Mechanism Extraction

做什么：
  从182篇CTC文献中，用LLM提取基因之间的调控关系、基因-表型关系。
  比如从"TWIST1 knockdown reduced CTC invasion"提取出：
    TWIST1 → DRIVES → CTC invasion (PMID:33143160)

输出：
  每条关系对应Neo4j中的一条 REGULATES 边。

分两步：
  Step 1: 获取182篇文献的摘要（Bio.Entrez API）
  Step 2: LLM提取结构化关系 → JSON → Neo4j
"""

import json, csv, time, os
from pathlib import Path

DATA_DIR = r'D:\科研\CTC知识图谱\v2\data'

# ============================================================
# Step 1: 获取所有182篇文献的摘要
# ============================================================
def fetch_abstracts(pmids, output_file):
    """
    从PubMed获取摘要。
    需要先安装: pip install biopython

    如果没有biopython，也可以跳过摘要获取，直接用旧KG中已有的结论数据。
    """
    from Bio import Entrez
    Entrez.email = "ctc_kg@research.org"  # NCBI要求提供邮箱

    abstracts = {}
    batch_size = 50

    print(f'Fetching abstracts for {len(pmids)} PMIDs...')
    pmid_list = sorted(pmids)

    for i in range(0, len(pmid_list), batch_size):
        batch = pmid_list[i:i+batch_size]
        try:
            handle = Entrez.efetch(db="pubmed", id=','.join(batch),
                                   rettype="xml", retmode="xml")
            records = Entrez.read(handle)
            handle.close()

            for article in records.get('PubmedArticle', []):
                medline = article.get('MedlineCitation', {})
                pmid = str(medline.get('PMID', ''))
                article_info = medline.get('Article', {})
                title = article_info.get('ArticleTitle', '')
                abstract_parts = article_info.get('Abstract', {}).get('AbstractText', [])

                # Abstract can be a list of strings or list of StringElement objects
                if isinstance(abstract_parts, list):
                    abstract = ' '.join(
                        str(p) if isinstance(p, str) else str(p)
                        for p in abstract_parts
                    )
                else:
                    abstract = str(abstract_parts)

                abstracts[pmid] = {'title': str(title), 'abstract': abstract}

            print(f'  Batch {i//batch_size + 1}: {len(batch)} PMIDs')
            time.sleep(0.5)

        except Exception as e:
            print(f'  Batch {i//batch_size + 1} ERROR: {e}')
            # Fallback: try individual PMIDs
            for pmid in batch:
                try:
                    handle = Entrez.efetch(db="pubmed", id=pmid,
                                           rettype="xml", retmode="xml")
                    records = Entrez.read(handle)
                    handle.close()
                    for article in records.get('PubmedArticle', []):
                        medline = article.get('MedlineCitation', {})
                        article_info = medline.get('Article', {})
                        title = article_info.get('ArticleTitle', '')
                        abstract_parts = article_info.get('Abstract', {}).get('AbstractText', [])
                        if isinstance(abstract_parts, list):
                            abstract = ' '.join(str(p) for p in abstract_parts)
                        else:
                            abstract = str(abstract_parts)
                        abstracts[str(medline.get('PMID', ''))] = {
                            'title': str(title), 'abstract': abstract
                        }
                    time.sleep(0.3)
                except Exception as e2:
                    print(f'    PMID {pmid}: {e2}')

    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(abstracts, f, ensure_ascii=False, indent=2)

    print(f'Saved {len(abstracts)} abstracts to {output_file}')
    return abstracts


# ============================================================
# Step 2: LLM提取结构化关系
# ============================================================

EXTRACTION_PROMPT = """You are an expert in CTC (Circulating Tumor Cell) biology and molecular oncology.

Extract MECHANISTIC relationships from the following PubMed abstract.
IMPORTANT: Only extract relationships where there is a clear causal or regulatory connection.
Do NOT extract genes that are merely "detected" or "measured" without functional context.

For each relationship, provide:
1. subject: the gene/protein that acts
2. predicate: the relationship type:
   - UPREGULATES / DOWNREGULATES (gene → gene expression change)
   - ACTIVATES / INHIBITS (gene → pathway or phenotype)
   - DRIVES / SUPPRESSES (gene → CTC phenotype)
   - PHOSPHORYLATES (gene → gene, post-translational)
2. object: the gene, pathway, or phenotype being affected
3. object_type: Gene, Pathway, or Phenotype
4. evidence_sentence: the EXACT sentence from the abstract supporting this claim
5. confidence: high/medium/low based on evidence strength
   - high: direct experimental evidence (knockdown, inhibitor, etc.)
   - medium: correlation or association evidence
   - low: indirect or inferred evidence

Phenotypes must be one of:
  EMT_induction, MET_induction, metastasis, invasion, migration,
  CTC_cluster_formation, immune_evasion, drug_resistance,
  anoikis_resistance, apoptosis_resistance, stemness,
  proliferation, angiogenesis, metabolic_reprogramming

Pathway names should be standard (e.g., "Wnt signaling", "PI3K/AKT/mTOR", "TGF-beta signaling")

Output ONLY valid JSON, no explanation:
{
  "pmid": "XXXXX",
  "relations": [
    {
      "subject": "TWIST1",
      "predicate": "DRIVES",
      "object": "EMT_induction",
      "object_type": "Phenotype",
      "evidence_sentence": "TWIST1 overexpression was sufficient to induce mesenchymal morphology",
      "confidence": "high"
    }
  ]
}

If no mechanistic relationships found, return {"pmid": "XXXXX", "relations": []}

ABSTRACT:
Title: {title}
Abstract: {abstract}
"""


def extract_relations_with_llm(abstracts_file, output_file, llm_client):
    """
    用LLM逐篇提取关系。

    llm_client: 一个函数，接受(prompt)返回(str)的JSON响应
    """
    with open(abstracts_file, 'r', encoding='utf-8') as f:
        abstracts = json.load(f)

    all_relations = []
    total = len(abstracts)

    for i, (pmid, data) in enumerate(abstracts.items()):
        title = data.get('title', '')
        abstract = data.get('abstract', '')

        if not abstract:
            print(f'  [{i+1}/{total}] PMID {pmid}: no abstract, skip')
            continue

        prompt = EXTRACTION_PROMPT.format(title=title, abstract=abstract[:3000])

        try:
            response = llm_client(prompt)
            result = json.loads(response)
            relations = result.get('relations', [])

            if relations:
                # Add PMID to each relation
                for rel in relations:
                    rel['pmid'] = pmid
                all_relations.extend(relations)
                print(f'  [{i+1}/{total}] PMID {pmid}: {len(relations)} relations')
            else:
                print(f'  [{i+1}/{total}] PMID {pmid}: 0 relations')

        except Exception as e:
            print(f'  [{i+1}/{total}] PMID {pmid}: ERROR - {e}')

        # Rate limiting
        time.sleep(0.5)

    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(all_relations, f, ensure_ascii=False, indent=2)

    print(f'\nTotal extracted relations: {len(all_relations)}')
    print(f'Saved to {output_file}')
    return all_relations


# ============================================================
# Step 3: 导入关系到Neo4j
# ============================================================
def import_relations_to_neo4j(relations_file):
    """将LLM提取的关系导入Neo4j作为REGULATES边"""
    from neo4j import GraphDatabase
    d = GraphDatabase.driver('bolt://127.0.0.1:7687', auth=(os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"]))

    with open(relations_file, 'r', encoding='utf-8') as f:
        relations = json.load(f)

    imported = 0
    with d.session() as s:
        for rel in relations:
            subject = rel['subject'].upper()
            predicate = rel['predicate']
            obj = rel['object']
            obj_type = rel.get('object_type', 'Phenotype')
            pmid = rel['pmid']
            evidence = rel.get('evidence_sentence', '')
            confidence = rel.get('confidence', 'medium')

            # Map confidence to float
            conf_map = {'high': 0.9, 'medium': 0.6, 'low': 0.3}
            conf_score = conf_map.get(confidence, 0.5)

            # For Phenotype target: create CTC_Phenotype node
            if obj_type == 'Phenotype':
                s.run("""
                    MERGE (ph:CTC_Phenotype {name: $name})
                """, name=obj)

            # Create REGULATES edge based on object type
            try:
                if obj_type == 'Gene':
                    s.run("""
                        MATCH (g1:Gene {symbol: $subject})
                        MATCH (g2:Gene {symbol: $object})
                        MERGE (g1)-[r:REGULATES {pmid: $pmid}]->(g2)
                        SET r.direction = $predicate,
                            r.mechanism = $evidence,
                            r.confidence = $conf,
                            r.pmid = $pmid
                    """, subject=subject, object=obj, predicate=predicate,
                         evidence=evidence[:200], conf=conf_score, pmid=pmid)

                elif obj_type == 'Phenotype':
                    s.run("""
                        MATCH (g1:Gene {symbol: $subject})
                        MATCH (ph:CTC_Phenotype {name: $object})
                        MERGE (g1)-[r:REGULATES {pmid: $pmid}]->(ph)
                        SET r.direction = $predicate,
                            r.mechanism = $evidence,
                            r.confidence = $conf,
                            r.pmid = $pmid
                    """, subject=subject, object=obj, predicate=predicate,
                         evidence=evidence[:200], conf=conf_score, pmid=pmid)

                elif obj_type == 'Pathway':
                    # Match pathway by name (fuzzy)
                    s.run("""
                        MATCH (g1:Gene {symbol: $subject})
                        MATCH (p:Pathway)
                        WHERE toLower(p.name) CONTAINS toLower($object)
                           OR toLower($object) CONTAINS toLower(p.name)
                        MERGE (g1)-[r:REGULATES {pmid: $pmid}]->(p)
                        SET r.direction = $predicate,
                            r.mechanism = $evidence,
                            r.confidence = $conf,
                            r.pmid = $pmid
                    """, subject=subject, object=obj, predicate=predicate,
                         evidence=evidence[:200], conf=conf_score, pmid=pmid)

                imported += 1
            except Exception as e:
                print(f'  Failed: {subject} -> {obj}: {e}')

    print(f'Imported {imported} REGULATES edges')
    d.close()


# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    import sys

    # Get PMIDs from old KG
    kg_file = r'D:\科研\CTC知识图谱\neo4j_query_table_data_2026-6-29.json'
    with open(kg_file, 'r', encoding='utf-8') as f:
        old_kg = json.load(f)

    pmids = set()
    for d in old_kg:
        if d['m'] is not None and 'PMID' in d['m']['labels']:
            pmids.add(d['m']['properties']['pmid'])
        if d['n'] is not None and 'PMID' in d['n']['labels']:
            pmids.add(d['n']['properties']['pmid'])

    print(f'Found {len(pmids)} unique PMIDs')

    abstracts_file = os.path.join(DATA_DIR, 'processed', 'abstracts.json')
    relations_file = os.path.join(DATA_DIR, 'processed', 'relations.json')

    if len(sys.argv) > 1 and sys.argv[1] == 'fetch':
        # Step 1: Fetch abstracts
        fetch_abstracts(pmids, abstracts_file)

    elif len(sys.argv) > 1 and sys.argv[1] == 'extract':
        # Step 2: LLM extraction
        # Define your LLM client here
        def llm_client(prompt):
            # TODO: Replace with actual API call
            # Example with OpenAI:
            # import openai
            # resp = openai.ChatCompletion.create(
            #     model="gpt-4", messages=[{"role":"user","content":prompt}],
            #     temperature=0)
            # return resp.choices[0].message.content
            raise NotImplementedError("Configure your LLM API client")

        extract_relations_with_llm(abstracts_file, relations_file, llm_client)

    elif len(sys.argv) > 1 and sys.argv[1] == 'import':
        # Step 3: Import to Neo4j
        import_relations_to_neo4j(relations_file)

    else:
        print("Usage:")
        print("  python extract_literature.py fetch    — Get abstracts from PubMed")
        print("  python extract_literature.py extract  — LLM extraction (needs API key)")
        print("  python extract_literature.py import   — Import relations to Neo4j")
