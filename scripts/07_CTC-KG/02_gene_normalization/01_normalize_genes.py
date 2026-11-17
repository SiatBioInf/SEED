"""
Gene Normalization: 给KG中的Gene节点添加HGNC标准ID

为什么需要这一步？
  现在Gene节点只有symbol（如"TP53"），没有标准ID。
  没有标准ID就无法和STRING、DrugBank等外部数据库对接。
  HGNC ID是基因的"身份证号"，全世界通用。

方法：
  用 mygene.info 免费API批量查询基因信息。
  每次查询最多1000个基因，35K基因需要约35次请求。
"""
import os

import csv, requests, json, time
from collections import defaultdict
from neo4j import GraphDatabase

NEO4J_URI = "bolt://127.0.0.1:7687"
NEO4J_AUTH = (os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"])

# Step 1: 从CSV中提取所有unique基因
print("Step 1: Extracting all unique genes from KG...")
genes = set()
with open(r"D:\科研\CTC知识图谱\v2\data\processed\gene_cancer_edges.csv", "r") as f:
    for row in csv.DictReader(f):
        genes.add(row["gene"])
print(f"  Found {len(genes)} unique gene symbols")

# Step 2: 用mygene.info批量查询（分批次，每批1000个）
print("Step 2: Querying mygene.info for HGNC IDs...")

gene_list = list(genes)
batch_size = 1000
gene_info = {}  # symbol -> {entrezgene, name, hgnc_id}

for i in range(0, len(gene_list), batch_size):
    batch = gene_list[i:i + batch_size]
    query = " ".join(batch)

    try:
        resp = requests.post(
            "http://mygene.info/v3/query",
            json={
                "q": query,
                "scopes": "symbol",
                "species": "human",
                "fields": "symbol,entrezgene,name",
                "size": batch_size
            },
            timeout=60
        )
        results = resp.json()

        for r in results:
            symbol = r.get("symbol", "").upper()
            entrez = r.get("entrezgene")
            name = r.get("name", "")
            if symbol and entrez:
                gene_info[symbol] = {
                    "entrez_id": str(entrez),
                    "hgnc_id": f"HGNC:{entrez}",
                    "name": name
                }

        print(f"  Batch {i//batch_size + 1}: queried {len(batch)} genes, found {len(results)}")

    except Exception as e:
        print(f"  Batch {i//batch_size + 1} failed: {e}")

    time.sleep(0.5)  # 礼貌性延迟

print(f"  Total matched: {len(gene_info)} / {len(genes)} ({100*len(gene_info)/len(genes):.1f}%)")

# Step 3: 更新Neo4j中的Gene节点
print("Step 3: Updating Gene nodes in Neo4j...")

driver = GraphDatabase.driver(NEO4J_URI, auth=NEO4J_AUTH)

# 注意：mygene返回的entrezgene不是真正的HGNC ID
# 真正的HGNC ID格式是 "HGNC:12428"
# entrezgene是NCBI Gene ID，格式如 "7291"
# 两者是不同的ID系统，但可以互相映射
# 我们存储entrez_id作为过渡，后续可以通过HGNC API转换

with driver.session() as session:
    updated = 0
    for symbol, info in gene_info.items():
        result = session.run("""
            MATCH (g:Gene {symbol: $symbol})
            SET g.entrez_id = $entrez_id,
                g.name = $name
            RETURN g.symbol
        """, symbol=symbol, entrez_id=info["entrez_id"], name=info["name"])
        if result.peek():
            updated += 1

    print(f"  Updated {updated} Gene nodes with Entrez IDs and names")

# Step 4: 验证
print("\nStep 4: Verification")
with driver.session() as session:
    # 随机查几个基因
    for symbol in ["TWIST1", "VIM", "KLK3", "CDH1", "ERBB2"]:
        result = session.run(
            "MATCH (g:Gene {symbol: $s}) RETURN g.symbol, g.entrez_id, g.name",
            s=symbol
        )
        rec = result.single()
        if rec:
            print(f"  {rec['g.symbol']:10s} Entrez:{rec['g.entrez_id']:10s}  {rec['g.name'][:60]}")

    # 统计
    result = session.run("""
        MATCH (g:Gene)
        RETURN count(g) AS total,
               count(g.entrez_id) AS with_id,
               count(g.name) AS with_name
    """)
    rec = result.single()
    print(f"\n  Total genes: {rec['total']}")
    print(f"  With Entrez ID: {rec['with_id']} ({100*rec['with_id']/rec['total']:.1f}%)")
    print(f"  With name: {rec['with_name']} ({100*rec['with_name']/rec['total']:.1f}%)")

driver.close()
print("\nDone.")
