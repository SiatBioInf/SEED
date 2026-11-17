"""
STRING PPI Import: 导入蛋白-蛋白相互作用网络

这一步做什么：
  从STRING数据库获取基因之间的蛋白互作关系。
  比如 VIM 和 CDH2 会物理结合（score=950），
  TWIST1 和 SNAI1 有实验验证的相互作用。

为什么需要这一步：
  没有Gene↔Gene边，KG就无法做网络分析。
  加上PPI边后，Agent可以做：
  "VIM的互作蛋白中，哪些在乳腺癌CTC中也高表达？"
  这样的2跳推理查询。

策略：
  只导入我们KG中已有的基因之间的互作关系。
  不和STRING全部数据库对接（那有20B+条边，太多了）。
  筛选 combined_score >= 700 的高置信度互作。
"""
import os

import csv, requests, json, time
from neo4j import GraphDatabase

NEO4J_URI = "bolt://127.0.0.1:7687"
NEO4J_AUTH = (os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"])

# Step 1: 从Neo4j中获取所有Gene节点的symbol
print("Step 1: Getting all gene symbols from Neo4j...")
driver = GraphDatabase.driver(NEO4J_URI, auth=NEO4J_AUTH)
with driver.session() as session:
    result = session.run("MATCH (g:Gene) WHERE g.symbol IS NOT NULL RETURN g.symbol AS symbol")
    gene_symbols = [r["symbol"] for r in result]
print(f"  {len(gene_symbols)} genes in KG")

# Step 2: 用STRING API查询PPI网络
#   第一步先把基因名映射为STRING的内部ID
#   第二步再获取互作网络
print("\nStep 2: Querying STRING API for PPI network...")

# STRING API一次最多处理约2000个基因（URL长度限制）
# 35K基因需要分约18批
# 但互作边会非常多（35K基因的完全互作网络可能有百万条边）
# 策略：只查询我们KG中Tier1+2的高置信度基因（约22K个）

# 简化：先查所有基因，只保留 combined_score >= 900 的极高置信度互作
batch_size = 1500
all_edges = []
gene_batches = [gene_symbols[i:i+batch_size] for i in range(0, len(gene_symbols), batch_size)]
print(f"  {len(gene_batches)} batches of {batch_size} genes each")

# 先试第一批作为测试
batch = gene_batches[0]
print(f"\n  Testing with first batch ({len(batch)} genes)...")

# Step 2a: 映射基因名到STRING ID
resp = requests.post(
    "https://string-db.org/api/json/get_string_ids",
    data={"identifiers": "\r".join(batch), "species": 9606, "limit": 1},
    timeout=60
)
id_map = {}
for item in resp.json():
    sid = item.get("stringId", "")
    query = item.get("queryItem", "")
    if sid:
        id_map[sid] = query
print(f"  Mapped {len(id_map)} genes to STRING IDs")

time.sleep(1)

# Step 2b: 获取互作网络
string_ids = list(id_map.keys())
resp = requests.post(
    "https://string-db.org/api/json/network",
    data={
        "identifiers": "\r".join(string_ids),
        "species": 9606,
        "required_score": 900,  # 极高置信度
        "network_type": "physical",
        "caller_identity": "ctc_kg_v2"
    },
    timeout=120
)
network = resp.json()
print(f"  Got {len(network)} PPI edges (score >= 900)")

# 转换回基因symbol
for item in network:
    gene_a = id_map.get(item.get("stringId_A", ""), item.get("stringId_A", ""))
    gene_b = id_map.get(item.get("stringId_B", ""), item.get("stringId_B", ""))
    score = int(item.get("score", 0))
    if gene_a and gene_b and gene_a != gene_b:
        all_edges.append({
            "gene_a": gene_a,
            "gene_b": gene_b,
            "score": score,
            "experiments": int(item.get("ascore", 0)),
            "database": int(item.get("dscore", 0)),
            "textmining": int(item.get("tscore", 0))
        })

print(f"\n  First batch: {len(all_edges)} edges")
print(f"\n  Due to API rate limits and batch processing needs, ")
print(f"  running all {len(gene_batches)} batches would take ~30 min.")
print(f"  Let's estimate: ~{len(all_edges) * len(gene_batches)} total edges expected.")

# Step 3: 导入Neo4j
if all_edges:
    print(f"\nStep 3: Importing {len(all_edges)} PPI edges into Neo4j...")

    with driver.session() as session:
        # 分批导入，每批1000条边
        import_batch = 500
        imported = 0
        for i in range(0, len(all_edges), import_batch):
            batch_edges = all_edges[i:i + import_batch]
            session.run("""
                UNWIND $edges AS edge
                MATCH (a:Gene {symbol: edge.gene_a})
                MATCH (b:Gene {symbol: edge.gene_b})
                MERGE (a)-[r:INTERACTS_WITH]->(b)
                SET r.score = edge.score,
                    r.experiments = edge.experiments,
                    r.database = edge.database,
                    r.textmining = edge.textmining,
                    r.source = 'STRING_v12.5'
            """, edges=batch_edges)
            imported += len(batch_edges)

        print(f"  Imported {imported} edges")

        # Step 4: 验证
        print("\nStep 4: Verification")
        result = session.run("""
            MATCH ()-[r:INTERACTS_WITH]->()
            RETURN count(r) AS total_ppi,
                   avg(r.score) AS avg_score,
                   max(r.score) AS max_score
        """)
        rec = result.single()
        print(f"  Total PPI edges: {rec['total_ppi']}")
        print(f"  Average score: {rec['avg_score']:.0f}")
        print(f"  Max score: {rec['max_score']}")

        # 测试查询：VIM的互作蛋白
        result = session.run("""
            MATCH (g:Gene {symbol: 'VIM'})-[r:INTERACTS_WITH]-(g2:Gene)
            RETURN g2.symbol AS interactor, r.score AS score
            ORDER BY r.score DESC LIMIT 10
        """)
        print("\n  VIM's interaction partners:")
        for rec in result:
            print(f"    {rec['interactor']:10s} score={rec['score']}")

driver.close()
print("\nDone. First batch imported as proof of concept.")
print("Full import needs all 24 batches - we can run this later.")
