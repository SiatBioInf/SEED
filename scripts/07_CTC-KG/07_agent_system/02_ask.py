"""
CTC-KG 问答入口

用法：在终端里直接问问题
  python ask.py "乳腺癌CTC有哪些高表达基因？"
  python ask.py "VIM参与什么通路？"
  python ask.py "有什么药可以靶向前列腺癌CTC的基因？"
  python ask.py "TWIST1在CTC中有什么证据？"
"""
import os

import sys
sys.path.insert(0, r'D:\科研\CTC知识图谱\v2')
sys.stdout.reconfigure(encoding='utf-8')
from src.agent.spark_agent import CTCSparkAgent

if len(sys.argv) < 2:
    print("用法: python ask.py \"你的问题\"")
    print()
    print("示例:")
    print('  python ask.py "乳腺癌CTC有哪些高表达基因？"')
    print('  python ask.py "VIM参与什么通路？"')
    print('  python ask.py "有什么药可以靶向前列腺癌CTC的基因？"')
    print('  python ask.py "TWIST1在CTC中有什么证据？"')
    print('  python ask.py "肺癌CTC中EMT通路相关的基因有哪些？"')
    sys.exit(0)

question = ' '.join(sys.argv[1:])

# 接S1-Base-Ultra做真正的自然语言理解
from openai import OpenAI
llm = OpenAI(
    api_key=os.environ["SEED_LLM_API_KEY"],
    base_url="https://uni-api.cstcloud.cn/v1"
)

def llm_plan(question):
    """用LLM理解问题，生成查询计划"""
    schema_desc = """
    你是一个Neo4j Cypher专家。数据库schema如下：

    节点:
    - Gene {symbol, name}  (35,679个)
    - CancerType {name}  (12个: Breast Cancer, Non-Small Cell Lung Cancer, Colorectal Cancer, Prostate Cancer, Pancreatic Ductal Adenocarcinoma, Melanoma, Hepatocellular Carcinoma, Gastric Cancer, Ovarian Cancer, Head and Neck Squamous Cell Carcinoma, Nasopharyngeal Carcinoma, Small Cell Lung Cancer)
    - Pathway {name, source_db}  (47个, MSigDB Hallmark通路)
    - Drug {name, drugbank_id}  (22个)
    - ClinicalOutcome {outcome_id, outcome_type, hr, p_value}
    - CTC_Phenotype {name}

    边:
    - (Gene)-[CTC_EXPRESSED_IN {mean_percentile, pct_delta, specificity_score, confidence_tier}]->(CancerType)
    - (Gene)-[MEMBER_OF_PATHWAY]->(Pathway)
    - (Gene)-[INTERACTS_WITH {score}]->(Gene)  (STRING PPI, score 0-1000)
    - (Drug)-[TARGETS]->(Gene)
    - (Gene)-[PROGNOSTIC_IN {hr, p_value}]->(ClinicalOutcome)
    - (Pathway)-[ACTIVE_IN {pct_delta}]->(CancerType)
    - (Gene)-[REGULATES {direction}]->(CTC_Phenotype)

    生成一个查询计划JSON。只输出JSON：
    {"steps": [{"name": "步骤名", "cypher": "Cypher查询", "params": {"param": "value"}}]}

    要求：
    1. 用 c.name = $cancer 精确匹配癌种名
    2. 基因用 g.symbol = $gene
    3. LIMIT不超过30
    4. 只用 MATCH/RETURN，不用 CREATE/DELETE/SET
    5. 不确定癌种名时用 CONTAINS
    """
    resp = llm.chat.completions.create(
        model="S1-Base-Ultra",
        messages=[
            {"role": "system", "content": schema_desc},
            {"role": "user", "content": question}
        ],
        temperature=0, max_tokens=2000
    )
    content = resp.choices[0].message.content.strip()
    # 提取JSON
    if '```json' in content:
        content = content.split('```json')[1].split('```')[0]
    elif '```' in content:
        content = content.split('```')[1].split('```')[0]
    if '{' in content:
        content = content[content.find('{'):content.rfind('}')+1]
    return json.loads(content)

agent = CTCSparkAgent(llm_client=llm, llm_plan_fn=llm_plan)

print('='*60)
print('Q:', question)
print('='*60)

result = agent.ask(question)

print()
print('A:', result['answer'])
print()
print('Query path:', result['query_path'])
print('Confidence:', result['confidence'])

# Show detailed evidence
print()
print('--- Detailed Evidence ---')
for step_name, rows in result['evidence'].items():
    if rows:
        print(f'\n[{step_name}] ({len(rows)} results)')
        for row in rows[:5]:
            # Show first 3 keys of each row
            keys = list(row.keys())[:4]
            vals = {k: str(row[k])[:80] for k in keys}
            print(f'  {vals}')

agent.close()
