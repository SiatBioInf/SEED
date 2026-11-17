"""
TCGA Survival Import: 导入基因-预后关系

做什么：
  从已有的TCGA生存分析结果中导入"基因高/低表达 → 生存预后"的关系。
  例如：HOXC10高表达在结直肠癌中 HR=1.47 (风险增加), P=0.041。

为什么需要：
  这是"KG发现 → TCGA验证 → CTC回验"三层验证体系的第二层。
  让Agent能回答："这个基因高表达对预后有什么影响？"
"""
import os

import sys, csv
sys.path.insert(0, r'D:\科研\CTC知识图谱\v2')
sys.stdout.reconfigure(encoding='utf-8')
from neo4j import GraphDatabase

d = GraphDatabase.driver('bolt://127.0.0.1:7687', auth=(os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"]))

# TCGA缩写 → 癌种全名映射
TCGA_MAP = {
    'BRCA': 'Breast Cancer',
    'LUAD': 'Non-Small Cell Lung Cancer',
    'COADREAD': 'Colorectal Cancer',
    'LIHC': 'Hepatocellular Carcinoma',
    'PAAD': 'Pancreatic Ductal Adenocarcinoma',
    'STAD': 'Gastric Cancer',
    'HNSC': 'Head and Neck Squamous Cell Carcinoma',
    'SKCM': 'Melanoma',
    'OV': 'Ovarian Cancer',
    'BLCA': 'Bladder Cancer',
    'PRAD': 'Prostate Cancer',
}

# Read TCGA data
tcga_file = r'D:\ClaudeCode记录\CTC-KG论文\three_module_tcga_all.csv'
rows = []
with open(tcga_file, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        rows.append(row)

print(f'TCGA rows: {len(rows)}')

# Filter to significant results and genes in our KG
imported = 0
with d.session() as s:
    for row in rows:
        gene = row['gene']
        cox_p = float(row['cox_p'])
        hr = float(row['HR'])

        # Only significant results (p < 0.05)
        if cox_p >= 0.05:
            continue

        # Map cancer name
        tcga_abbr = row['cancer_tcga']
        cancer_name = TCGA_MAP.get(tcga_abbr, tcga_abbr)

        # Check gene exists in KG
        check = s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=gene).single()
        if not check:
            continue

        # Check cancer exists
        check = s.run("MATCH (c:CancerType {name: $n}) RETURN c", n=cancer_name).single()
        if not check:
            continue

        # Create ClinicalOutcome node and PROGNOSTIC_IN edge
        outcome_id = 'TCGA_' + tcga_abbr + '_' + gene + '_OS'
        hr_lower = float(row.get('HR_lower', 0) or 0)
        hr_upper = float(row.get('HR_upper', 0) or 0)
        logrank_p = float(row.get('logrank_p', 0) or 0)
        n_patients = int(row.get('n', 0) or 0)

        s.run("""
            MERGE (co:ClinicalOutcome {outcome_id: $oid})
              ON CREATE SET co.outcome_type = 'Overall Survival',
                            co.cancer_type = $cancer,
                            co.cohort = 'TCGA',
                            co.n_patients = $n,
                            co.hr = $hr,
                            co.hr_ci_lower = $hr_low,
                            co.hr_ci_upper = $hr_high,
                            co.cox_p = $cox_p,
                            co.logrank_p = $logrank_p
            WITH co
            MATCH (g:Gene {symbol: $gene})
            MERGE (g)-[r:PROGNOSTIC_IN]->(co)
            SET r.hr = $hr,
                r.hr_ci_lower = $hr_low,
                r.hr_ci_upper = $hr_high,
                r.p_value = $cox_p,
                r.effect_direction = CASE WHEN $hr > 1 THEN 'unfavorable' ELSE 'favorable' END,
                r.data_source = 'TCGA',
                r.cancer_type = $cancer
        """, oid=outcome_id, gene=gene, cancer=cancer_name,
             hr=hr, hr_low=hr_lower, hr_high=hr_upper,
             cox_p=cox_p, logrank_p=logrank_p, n=n_patients)

        imported += 1

print(f'Imported {imported} prognostic relationships')

# Verify
with d.session() as s:
    r = s.run("MATCH ()-[r:PROGNOSTIC_IN]->() RETURN count(r) AS n").single()
    print(f'PROGNOSTIC_IN edges: {r["n"]}')

    r = s.run("MATCH (co:ClinicalOutcome) RETURN count(co) AS n").single()
    print(f'ClinicalOutcome nodes: {r["n"]}')

    # Top unfavorable genes
    print('\nTop unfavorable prognostic genes (HR > 2):')
    r = s.run("""
        MATCH (g:Gene)-[r:PROGNOSTIC_IN]->(co:ClinicalOutcome)
        WHERE r.hr > 2.0
        RETURN g.symbol AS gene, co.cancer_type AS cancer,
               round(r.hr, 2) AS HR, round(r.p_value, 4) AS p
        ORDER BY r.hr DESC LIMIT 10
    """)
    for row in r:
        print(f'  {row["gene"]:10s} {row["cancer"]:35s} HR={row["HR"]:.2f} p={row["p"]:.4f}')

d.close()
print('Done.')
