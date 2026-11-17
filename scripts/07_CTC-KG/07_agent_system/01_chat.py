"""
CTC-KG 对话助手 v2

策略：不让LLM凭空写Cypher（容易出错），而是让LLM选模板+填参数。
模板是预定义的正确Cypher，LLM只负责理解问题和填参数。
"""
import os

import sys, json
sys.path.insert(0, r'D:\科研\CTC知识图谱\v2')
sys.stdout.reconfigure(encoding='utf-8')
try: import readline
except ImportError: pass
from openai import OpenAI
from neo4j import GraphDatabase

llm = OpenAI(api_key=os.environ["SEED_LLM_API_KEY"],
             base_url="https://uni-api.cstcloud.cn/v1")
db = GraphDatabase.driver('bolt://127.0.0.1:7687', auth=(os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"]))

# ============================================================
# 预定义模板（已验证，不会出错）
# ============================================================
TEMPLATES = {
    "genes_in_cancer": {
        "desc": "某癌种CTC中高表达基因（含功能分类）",
        "params": ["cancer"],
        "cypher": """
            MATCH (g:Gene)-[e:CTC_EXPRESSED_IN]->(c:CancerType)
            WHERE toLower(c.name) CONTAINS toLower($cancer) AND e.confidence_tier <= 2
            RETURN g.symbol AS gene, g.protein_class AS function,
                   e.specificity_score AS specificity,
                   round(e.cohens_d, 2) AS effect_size, e.confidence_tier AS tier
            ORDER BY e.specificity_score DESC LIMIT 20
        """
    },
    "gene_network": {
        "desc": "基因的完整网络（功能+通路+互作+表达+药物+预后+表型）",
        "params": ["gene"],
        "cypher": """
            MATCH (g:Gene {symbol: $gene})
            OPTIONAL MATCH (g)-[:MEMBER_OF_PATHWAY]->(p:Pathway)
            OPTIONAL MATCH (g)-[ppi:INTERACTS_WITH]-(g2:Gene)
            WHERE ppi.score > 700
            OPTIONAL MATCH (g)-[e:CTC_EXPRESSED_IN]->(ct:CancerType)
            WHERE e.confidence_tier <= 2
            OPTIONAL MATCH (d:Drug)-[:TARGETS]->(g)
            OPTIONAL MATCH (g)-[prog:PROGNOSTIC_IN]->(co:ClinicalOutcome)
            OPTIONAL MATCH (g)-[reg:REGULATES]->(target)
            RETURN g.symbol AS gene, g.name AS gene_name,
                   g.protein_class AS function, g.function_summary AS summary,
                   collect(DISTINCT p.name) AS pathways,
                   collect(DISTINCT {gene: g2.symbol, score: ppi.score}) AS interactors,
                   collect(DISTINCT {cancer: ct.name, specificity: e.specificity_score}) AS expressed_in,
                   collect(DISTINCT d.name) AS drugs,
                   collect(DISTINCT {outcome: co.outcome_type, HR: prog.hr, p_value: prog.p_value}) AS prognosis
        """
    },
    "drugs_for_cancer": {
        "desc": "靶向某癌种CTC高表达基因的药物",
        "params": ["cancer"],
        "cypher": """
            MATCH (g:Gene)-[e:CTC_EXPRESSED_IN]->(c:CancerType)
            WHERE toLower(c.name) CONTAINS toLower($cancer) AND e.confidence_tier <= 2
            WITH g, e ORDER BY e.specificity_score DESC LIMIT 20
            MATCH (d:Drug)-[:TARGETS]->(g)
            RETURN g.symbol AS gene, round(e.specificity_score, 3) AS ctc_specificity,
                   collect(DISTINCT d.name) AS drugs
        """
    },
    "pathways_in_cancer": {
        "desc": "某癌种CTC中活跃的通路",
        "params": ["cancer"],
        "cypher": """
            MATCH (p:Pathway)-[a:ACTIVE_IN]->(c:CancerType)
            WHERE toLower(c.name) CONTAINS toLower($cancer) AND a.confidence_tier <= 2
            RETURN p.name AS pathway, round(a.pct_delta, 4) AS activity
            ORDER BY a.pct_delta DESC LIMIT 20
        """
    },
    "pan_cancer_genes": {
        "desc": "跨多个癌种高表达的泛癌基因",
        "params": ["min_cancers"],
        "cypher": """
            MATCH (g:Gene)-[e:CTC_EXPRESSED_IN]->(c:CancerType)
            WHERE e.confidence_tier <= 2 AND e.mean_percentile > 0.7
            WITH g.symbol AS gene, count(c) AS n_cancers, avg(e.mean_percentile) AS avg_pct
            WHERE n_cancers >= $min_cancers
            RETURN gene, n_cancers, round(avg_pct, 3) AS avg_percentile
            ORDER BY n_cancers DESC LIMIT 30
        """
    },
    "gene_specificity": {
        "desc": "某基因在各癌种中的表达特异性",
        "params": ["gene"],
        "cypher": """
            MATCH (g:Gene {symbol: $gene})-[e:CTC_EXPRESSED_IN]->(c:CancerType)
            RETURN c.name AS cancer, round(e.mean_percentile, 3) AS mean_pct,
                   round(e.pct_delta, 3) AS delta, e.confidence_tier AS tier
            ORDER BY e.mean_percentile DESC
        """
    },
    "prognosis_by_cancer": {
        "desc": "某癌种中与预后显著相关的基因",
        "params": ["cancer"],
        "cypher": """
            MATCH (g:Gene)-[r:PROGNOSTIC_IN]->(co:ClinicalOutcome)
            WHERE r.p_value < 0.05 AND toLower(co.cancer_type) CONTAINS toLower($cancer)
            RETURN g.symbol AS gene, round(r.hr, 2) AS HR,
                   round(r.p_value, 4) AS p_value, co.cancer_type AS cancer
            ORDER BY r.p_value ASC LIMIT 20
        """
    },
}

# ============================================================
# LLM Router：选模板 + 填参数
# ============================================================
ROUTER_PROMPT = """You are a query router for a CTC knowledge graph. Choose the best template and fill parameters.

Available templates:
{templates}

Cancer names: Breast Cancer, Non-Small Cell Lung Cancer, Colorectal Cancer, Prostate Cancer, Pancreatic Ductal Adenocarcinoma, Melanoma, Hepatocellular Carcinoma, Gastric Cancer, Ovarian Cancer, Head and Neck Squamous Cell Carcinoma, Nasopharyngeal Carcinoma, Small Cell Lung Cancer

Output ONLY valid JSON:
{{"template": "template_name", "params": {{"param": "value"}}}}"""

def route_question(question):
    """LLM选择模板和参数"""
    template_desc = "\n".join(f"- {k}: {v['desc']} (params: {v['params']})" for k, v in TEMPLATES.items())

    resp = llm.chat.completions.create(
        model="S1-Base-Ultra",
        messages=[
            {"role": "system", "content": ROUTER_PROMPT.format(templates=template_desc)},
            {"role": "user", "content": question}
        ],
        temperature=0, max_tokens=300
    )
    content = resp.choices[0].message.content.strip()
    if '```' in content:
        content = content.split('```')[1].split('```')[0]
    try:
        return json.loads(content)
    except json.JSONDecodeError:
        if '{' in content:
            return json.loads(content[content.find('{'):content.rfind('}')+1])
        return None

def execute_query(template_name, params):
    """执行模板查询"""
    t = TEMPLATES.get(template_name)
    if not t:
        return None
    with db.session() as s:
        result = s.run(t['cypher'], **params)
        return [dict(r) for r in result]

# ============================================================
# 对话循环
# ============================================================
print("=" * 60)
print("  CTC-KG 对话助手 (S1-Base-Ultra + 模板)")
print("  输入问题，输入 quit 退出")
print("=" * 60)

while True:
    try:
        question = input("\nYou: ").strip()
    except (EOFError, KeyboardInterrupt):
        break

    if question.lower() in ('quit', 'exit', 'q'):
        print("再见！")
        break
    if not question:
        continue

    # Step 1: LLM选模板+填参数
    plan = route_question(question)
    if not plan:
        print("Agent: 抱歉，没理解你的问题。试试问'乳腺癌有哪些高表达基因'这样的问题？")
        continue

    tname = plan.get('template')
    params = plan.get('params', {})
    tdesc = TEMPLATES.get(tname, {}).get('desc', tname)

    # Step 2: 执行
    print(f"(查询: {tdesc})", end=" ", flush=True)
    results = execute_query(tname, params)

    if results is None:
        print("模板不存在")
        continue

    if not results:
        # 没结果，尝试换个模板
        print("(没有直接结果)", end=" ")
        # Try gene_network if a gene was asked
        if 'gene' in params:
            results = execute_query('gene_network', params)
        elif 'cancer' in params:
            results = execute_query('genes_in_cancer', params)

    if not results:
        print("Agent: 没有找到相关数据。")
        continue

    # Step 3: LLM生成自然语言答案
    results_str = json.dumps(results[:15], ensure_ascii=False, indent=2)
    resp = llm.chat.completions.create(
        model="S1-Base-Ultra",
        messages=[
            {"role": "system", "content": "你是CTC研究助手。根据查询结果用自然语言回答用户问题，引用具体数值。用中文回答。简洁。"},
            {"role": "user", "content": f"问题: {question}\n\n查询结果:\n{results_str}\n\n请回答:"}
        ],
        temperature=0, max_tokens=800
    )
    print("Agent:", resp.choices[0].message.content)

db.close()
