"""Quick stats on R script output files"""
import csv, os
from collections import Counter, defaultdict

def load_csv(path):
    with open(path, 'r') as f:
        return list(csv.DictReader(f))

edges = load_csv(r'D:\科研\CTC知识图谱\v2\data\processed\gene_cancer_edges.csv')
print('=== Gene-Cancer Edges ===')
print(f'Total: {len(edges)}')

tiers = Counter(int(r['confidence_tier']) for r in edges)
for t in [1,2,3]:
    pct = 100 * tiers[t] / len(edges)
    print(f'  Tier {t}: {tiers[t]:5d} ({pct:.1f}%)')

ct_genes = defaultdict(lambda: {1:0, 2:0, 3:0})
ct_cells = {}
for r in edges:
    ct_genes[r['cancer']][int(r['confidence_tier'])] += 1
    ct_cells[r['cancer']] = int(r['n_cells'])

print('\nPer cancer (Tier 1+2 genes):')
for cn in sorted(ct_genes.keys()):
    t12 = ct_genes[cn][1] + ct_genes[cn][2]
    nc = ct_cells.get(cn, 0)
    flag = ' ** <50 cells, unreliable **' if nc < 50 else ''
    print(f'  {cn:6s}: {t12:4d} genes (n={nc} cells){flag}')

print('\nTop gene per cancer by effect size (Cohen d):')
for cn in sorted(ct_genes.keys()):
    cn_edges = [r for r in edges if r['cancer'] == cn and int(r['confidence_tier']) == 1]
    if cn_edges:
        best = max(cn_edges, key=lambda r: float(r['cohens_d'] or 0))
        d = float(best['cohens_d'])
        delta = float(best['pct_delta'])
        print(f'  {cn:6s}: {best["gene"]:12s} delta={delta:.4f} d={d:.2f}')

# Pathway edges
pw_path = r'D:\科研\CTC知识图谱\v2\data\processed\pathway_cancer_edges.csv'
if os.path.exists(pw_path):
    pw = load_csv(pw_path)
    print(f'\n=== Pathway-Cancer Edges: {len(pw)} ===')
    pw_tiers = Counter(int(r['confidence_tier']) for r in pw)
    for t in sorted(pw_tiers):
        print(f'  Tier {t}: {pw_tiers[t]}')
    if pw:
        pw_sorted = sorted(pw, key=lambda r: float(r['pct_delta'] or 0), reverse=True)
        print('Top 10 pathways:')
        for r in pw_sorted[:10]:
            pname = r['pathway'][:50]
            delta = float(r['pct_delta'])
            print(f'  {pname:50s} {r["cancer"]:6s} delta={delta:.4f}')

gp_path = r'D:\科研\CTC知识图谱\v2\data\processed\gene_pathway_edges.csv'
if os.path.exists(gp_path):
    gp = load_csv(gp_path)
    genes_set = {r['gene'] for r in gp}
    pws_set = {r['pathway'] for r in gp}
    print(f'\n=== Gene-Pathway Edges: {len(gp)} ===')
    print(f'Unique genes: {len(genes_set)}, Unique pathways: {len(pws_set)}')
