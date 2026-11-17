"""
DrugBank Import: 导入药物-靶点关系

这一步做什么：
  从DrugBank数据库获取"哪些药靶向哪些基因"的关系。
  比如 Trastuzumab (赫赛汀) 靶向 ERBB2 (HER2)。

为什么需要这一步：
  加上药物靶点关系后，Agent可以做：
  "乳腺癌CTC中高表达的基因，有哪些已有药物可以靶向？"
  这类药物重利用（drug repurposing）查询。

前置条件：
  需要DrugBank免费学术账号（go.drugbank.com 注册即可）
  下载DrugBank XML完整数据库。

替代方案（如果不想下载DrugBank）：
  用 mygene.info API 也可以获取部分drug-gene关系，
  但不如DrugBank完整。
"""
import os

import csv
from neo4j import GraphDatabase

NEO4J_URI = "bolt://127.0.0.1:7687"
NEO4J_AUTH = (os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"])

# ============================================================
# 方案A: 使用内置CTC相关药物列表（简化版，不需要下载DrugBank）
# ============================================================
#
# 这些是已知靶向CTC相关基因的临床药物。
# 来源：手工整理 + 公共知识。
# 不如DrugBank完整但立即可用，后续可以替换为DrugBank完整数据。

CTC_DRUG_TARGETS = [
    # 靶向ERBB2 (HER2) —— 乳腺癌/胃癌CTC
    ("Trastuzumab", "DB00072", "ERBB2", "antibody"),
    ("Pertuzumab", "DB06366", "ERBB2", "antibody"),
    ("Lapatinib", "DB01259", "ERBB2", "inhibitor"),
    ("Neratinib", "DB11828", "ERBB2", "inhibitor"),
    ("Tucatinib", "DB15442", "ERBB2", "inhibitor"),

    # 靶向EGFR —— 肺癌/结直肠癌CTC
    ("Cetuximab", "DB00002", "EGFR", "antibody"),
    ("Gefitinib", "DB00317", "EGFR", "inhibitor"),
    ("Erlotinib", "DB00530", "EGFR", "inhibitor"),
    ("Osimertinib", "DB09330", "EGFR", "inhibitor"),

    # 靶向VEGFA —— 抗血管生成
    ("Bevacizumab", "DB00112", "VEGFA", "antibody"),

    # 靶向PI3K/AKT/mTOR通路
    ("Everolimus", "DB01590", "MTOR", "inhibitor"),
    ("Alpelisib", "DB12315", "PIK3CA", "inhibitor"),

    # 靶向AR —— 前列腺癌CTC
    ("Enzalutamide", "DB08899", "AR", "antagonist"),
    ("Abiraterone", "DB05812", "AR", "inhibitor"),

    # 靶向ESR1 —— 乳腺癌CTC
    ("Tamoxifen", "DB00675", "ESR1", "modulator"),
    ("Fulvestrant", "DB00947", "ESR1", "antagonist"),

    # 靶向CDK4/6 —— 乳腺癌
    ("Palbociclib", "DB09073", "CDK4", "inhibitor"),
    ("Ribociclib", "DB11730", "CDK6", "inhibitor"),

    # 靶向PARP —— DNA修复
    ("Olaparib", "DB09074", "PARP1", "inhibitor"),

    # 靶向BCL2 —— 凋亡
    ("Venetoclax", "DB11581", "BCL2", "inhibitor"),

    # 靶向ALK/ROS1 —— 肺癌
    ("Crizotinib", "DB08865", "ALK", "inhibitor"),

    # 靶向MET
    ("Cabozantinib", "DB08875", "MET", "inhibitor"),

    # 靶向免疫检查点
    ("Pembrolizumab", "DB09037", "CD274", "antibody"),  # PD-L1
    ("Nivolumab", "DB09035", "CD274", "antibody"),
]

print("DrugBank Import")
print("=" * 60)

driver = GraphDatabase.driver(NEO4J_URI, auth=NEO4J_AUTH)

with driver.session() as session:
    # 创建Drug节点和TARGETS边
    for drug_name, drugbank_id, gene_symbol, action in CTC_DRUG_TARGETS:
        # 检查基因是否在KG中
        check = session.run(
            "MATCH (g:Gene {symbol: $symbol}) RETURN g.symbol",
            symbol=gene_symbol
        )
        if not check.single():
            print(f"  SKIP {drug_name} -> {gene_symbol}: gene not in KG")
            continue

        # 创建/更新Drug节点和TARGETS边
        session.run("""
            MERGE (d:Drug {drugbank_id: $drugbank_id})
              ON CREATE SET d.name = $name
            WITH d
            MATCH (g:Gene {symbol: $gene})
            MERGE (d)-[r:TARGETS]->(g)
            SET r.action = $action,
                r.known_action = 'yes',
                r.source = 'DrugBank_manual'
        """, drugbank_id=drugbank_id, name=drug_name,
             gene=gene_symbol, action=action)
        print(f"  OK   {drug_name} -> {gene_symbol} ({action})")

    # 验证
    print("\nVerification:")
    result = session.run("""
        MATCH (d:Drug)-[r:TARGETS]->(g:Gene)
        RETURN count(d) AS drugs, count(r) AS targets
    """)
    rec = result.single()
    print(f"  Drug nodes: {rec['drugs']}")
    print(f"  TARGETS edges: {rec['targets']}")

    # 测试查询：乳腺癌CTC高表达基因中有哪些被药物靶向？
    print("\nTest: Drugs targeting genes highly expressed in Breast Cancer CTCs:")
    result = session.run("""
        MATCH (g:Gene)-[e:CTC_EXPRESSED_IN]->(c:CancerType {name: 'Breast Cancer'})
        WHERE e.confidence_tier <= 2
        MATCH (d:Drug)-[t:TARGETS]->(g)
        RETURN g.symbol AS gene, e.specificity_score AS score,
               collect(DISTINCT {drug: d.name, action: t.action}) AS drugs
        ORDER BY score DESC LIMIT 10
    """)
    for rec in result:
        drugs_str = ", ".join(f"{d['drug']}({d['action']})" for d in rec['drugs'])
        print(f"  {rec['gene']:10s} (spec={rec['score']:.2f}): {drugs_str}")

driver.close()
print("\nDone. For complete drug-target data, download DrugBank XML from go.drugbank.com")
