"""
给KG中的Gene节点添加功能注释

数据来源: mygene.info (NCBI Gene + UniProt + GO)
每个基因加上: protein_class(功能分类), summary(一句话描述), go_terms(GO功能)

为什么有意义：
  Agent不只告诉你"RPLP0P2效应量2.39"，
  而是"RPLP0P2是一个核糖体蛋白，效应量2.39"——
  研究者一眼就知道这是管家基因不是真正的靶点。
"""
import os

import requests, time, json, sys
sys.path.insert(0, r'D:\科研\CTC知识图谱\v2')
sys.stdout.reconfigure(encoding='utf-8')
from neo4j import GraphDatabase

d = GraphDatabase.driver('bolt://127.0.0.1:7687', auth=(os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"]))

# ============================================================
# 功能分类规则（基于GO和基因名）
# ============================================================
def classify_gene(symbol, info):
    """根据GO term + 基因描述判断基因的功能类别"""
    name = (info.get('name') or '').lower()
    summary = (info.get('summary') or '').lower()
    go_data = info.get('go')
    if isinstance(go_data, dict):
        go_bp_raw = go_data.get('BP', [])
        go_mf_raw = go_data.get('MF', [])
        go_bp = [g.get('term','') if isinstance(g, dict) else str(g) for g in (go_bp_raw if isinstance(go_bp_raw, list) else [go_bp_raw])]
        go_mf = [g.get('term','') if isinstance(g, dict) else str(g) for g in (go_mf_raw if isinstance(go_mf_raw, list) else [go_mf_raw])]
    else:
        go_bp = []
        go_mf = []
    all_text = name + ' ' + summary + ' ' + ' '.join(go_bp) + ' ' + ' '.join(go_mf)

    categories = []

    # 受体
    if any(kw in all_text for kw in ['receptor', 'gpcr', 'tyrosine kinase receptor']):
        if 'tyrosine kinase' in all_text or 'egf' in all_text or symbol in ['EGFR','ERBB2','ERBB3','ERBB4','MET','KDR','FLT1','FLT4','PDGFRA','PDGFRB','IGF1R','FGFR','KIT','RET','ALK','ROS1','AXL','EPHA','EPHB','NTRK']:
            categories.append('Receptor Tyrosine Kinase')
        elif 'g protein' in all_text or 'gpcr' in all_text:
            categories.append('GPCR')
        elif 'nuclear receptor' in all_text or 'steroid receptor' in all_text:
            categories.append('Nuclear Receptor')
        elif 'immune' in all_text or symbol.startswith('CD') or symbol in ['CD3E','CD3D','CD4','CD8A','CD8B','CD19','CD274','CTLA4','PDCD1']:
            categories.append('Immune Receptor')
        else:
            categories.append('Receptor')

    # 转录因子
    if any(kw in all_text for kw in ['transcription factor', 'transcriptional', 'homeobox', 'zinc finger', 'forkhead', 'helix-loop-helix', 'bhlh', 'leucine zipper', 'nuclear factor']):
        categories.append('Transcription Factor')

    # 激酶
    if 'kinase' in all_text and 'tyrosine kinase' not in ' '.join(categories) and 'receptor tyrosine' not in ' '.join(categories):
        categories.append('Kinase')

    # 细胞粘附
    if any(kw in all_text for kw in ['cadherin', 'cell adhesion', 'integrin', 'claudin', 'occludin', 'junctional', 'tight junction', 'adherens', 'desmosome']):
        categories.append('Cell Adhesion')

    # 细胞骨架
    if any(kw in all_text for kw in ['intermediate filament', 'actin', 'cytoskeleton', 'microtubule', 'tubulin', 'keratin', 'vimentin', 'myosin', 'tropomyosin']):
        categories.append('Cytoskeleton')

    # 分泌蛋白/蛋白酶
    if any(kw in all_text for kw in ['secreted', 'extracellular', 'serine protease', 'kallikrein', 'metalloproteinase', 'mmp', 'trypsin', 'chymotrypsin', 'growth factor', 'cytokine', 'chemokine', 'interleukin']):
        if any(kw in all_text for kw in ['growth factor']):
            categories.append('Growth Factor')
        elif any(kw in all_text for kw in ['cytokine', 'chemokine', 'interleukin']):
            categories.append('Cytokine/Chemokine')
        else:
            categories.append('Secreted Protein')

    # 免疫相关
    if any(kw in all_text for kw in ['immune', 'immunoglobulin', 'antibody', 't cell', 'b cell', 'nk cell', 'macrophage', 'neutrophil', 'lymphocyte', 'mhc', 'hla', 'antigen', 'complement']):
        categories.append('Immune Related')

    # 代谢酶
    if any(kw in all_text for kw in ['dehydrogenase', 'synthase', 'reductase', 'oxidase', 'transferase', 'hydrolase', 'isomerase', 'lyase', 'ligase', 'metabolism', 'metabolic', 'glycolysis', 'oxidative phosphorylation', 'tca']):
        categories.append('Metabolic Enzyme')

    # 转运蛋白
    if any(kw in all_text for kw in ['transporter', 'channel', 'pump', 'abc transporter', 'solute carrier', 'ion channel', 'aquaporin']):
        categories.append('Transporter')

    # 染色质/表观遗传
    if any(kw in all_text for kw in ['chromatin', 'histone', 'methyltransferase', 'acetyltransferase', 'epigenetic', 'nucleosome', 'dna methyl', 'hdac']):
        categories.append('Epigenetic/Chromatin')

    # 核糖体
    if any(kw in all_text for kw in ['ribosomal', 'ribosome', 'translation', 'translational']) or symbol.startswith(('RPS','RPL')):
        categories.append('Ribosomal')

    # DNA修复
    if any(kw in all_text for kw in ['dna repair', 'dna damage', 'mismatch repair', 'base excision', 'nucleotide excision', 'homologous recombination']):
        categories.append('DNA Repair')

    # 细胞周期
    if any(kw in all_text for kw in ['cell cycle', 'cyclin', 'cdk', 'mitotic', 'mitosis', 'checkpoint', 'centrosome', 'spindle']):
        categories.append('Cell Cycle')

    # 凋亡
    if any(kw in all_text for kw in ['apoptosis', 'apoptotic', 'caspase', 'bcl', 'death receptor']):
        categories.append('Apoptosis')

    # 泛素化/降解
    if any(kw in all_text for kw in ['ubiquitin', 'proteasome', 'e3 ligase', 'e2 conjugating', 'degradation']):
        categories.append('Protein Degradation')

    # 线粒体
    if symbol.startswith('MT-') or any(kw in all_text for kw in ['mitochondrial', 'mitochondria']):
        categories.append('Mitochondrial')

    # 没有匹配到 → Other
    if not categories:
        categories.append('Other/Unknown')

    return categories


# ============================================================
# 批量查询 + 更新
# ============================================================
print("Step 1: Getting gene symbols from Neo4j...")
with d.session() as s:
    genes = [r['symbol'] for r in s.run("MATCH (g:Gene) RETURN g.symbol AS symbol")]
print(f"  {len(genes)} genes")

# 分批查询 mygene.info（每批1000个，每批间隔0.5s）
batch_size = 1000
updated = 0
batches = [genes[i:i+batch_size] for i in range(0, len(genes), batch_size)]

print(f"Step 2: Querying mygene.info in {len(batches)} batches...")
all_gene_info = {}

for bi, batch in enumerate(batches):
    try:
        resp = requests.post('http://mygene.info/v3/query', json={
            'q': ' '.join(batch), 'scopes': 'symbol', 'species': 'human',
            'fields': 'symbol,name,summary,go.BP,go.MF,uniprot',
            'size': batch_size
        }, timeout=120)
        for item in resp.json():
            sym = item.get('symbol', '').upper()
            if sym:
                all_gene_info[sym] = item

        print(f"  Batch {bi+1}/{len(batches)}: {len(all_gene_info)} genes so far")
    except Exception as e:
        print(f"  Batch {bi+1} ERROR: {e}")

    time.sleep(0.5)

print(f"\n  Got info for {len(all_gene_info)} / {len(genes)} genes")

# ============================================================
# Step 3: 更新Neo4j
# ============================================================
print("\nStep 3: Updating Gene nodes (batch mode)...")

# Prepare batch data
batch_data = []
for symbol, info in all_gene_info.items():
    name = info.get('name', '')
    summary = (info.get('summary') or '')[:300]
    categories = classify_gene(symbol, info)
    go_data = info.get('go')
    if isinstance(go_data, dict):
        go_bp_raw = go_data.get('BP', [])
        go_bp = [g.get('term','') if isinstance(g, dict) else str(g) for g in (go_bp_raw if isinstance(go_bp_raw, list) else [go_bp_raw])]
    else:
        go_bp = []
    go_terms = list(dict.fromkeys(go_bp))[:8]
    uniprot_id = (info.get('uniprot') or {}).get('Swiss-Prot', '') if isinstance(info.get('uniprot'), dict) else ''

    batch_data.append({
        'symbol': symbol, 'name': name, 'summary': summary,
        'categories': categories, 'go_terms': go_terms,
        'uniprot_id': uniprot_id
    })

# Batch update (500 per transaction)
batch_size = 500
with d.session() as s:
    for i in range(0, len(batch_data), batch_size):
        batch = batch_data[i:i+batch_size]
        s.run("""
            UNWIND $rows AS row
            MATCH (g:Gene {symbol: row.symbol})
            SET g.name = row.name,
                g.function_summary = row.summary,
                g.protein_class = row.categories,
                g.go_biological_process = row.go_terms,
                g.uniprot_id = row.uniprot_id
        """, rows=batch)
        updated += len(batch)
        if (i // batch_size) % 10 == 0:
            print(f"  Updated {updated}/{len(batch_data)} genes...")

print(f"\n  Total updated: {updated} genes")

# ============================================================
# Step 4: 验证
# ============================================================
print("\nStep 4: Verification")
with d.session() as s:
    for gene in ['TWIST1','VIM','KLK3','ERBB2','EPCAM','EGFR','CDH1','SNAI1','RPLP0P2','KRT19']:
        r = s.run("""
            MATCH (g:Gene {symbol: $s})
            RETURN g.symbol AS symbol, g.protein_class AS class, g.name AS name,
                   g.function_summary AS summary
        """, s=gene).single()
        if r:
            cls = ', '.join(r['class']) if r['class'] else 'N/A'
            print(f"  {r['symbol']:10s} [{cls}]")
            print(f"    {r['name']}")

    # 功能分类统计
    r = s.run("""
        MATCH (g:Gene) WHERE g.protein_class IS NOT NULL
        UNWIND g.protein_class AS cls
        RETURN cls, count(g) AS n
        ORDER BY n DESC LIMIT 15
    """)
    print("\n  Gene functional categories:")
    for row in r:
        print(f"    {row['cls']:30s} {row['n']} genes")

d.close()
print("\nDone.")
