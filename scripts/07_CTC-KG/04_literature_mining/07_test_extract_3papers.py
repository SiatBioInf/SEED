"""Quick test: extract 3 papers to verify pipeline works"""
import os
import json, time, sys
sys.path.insert(0, r'D:\科研\CTC知识图谱\v2')
sys.stdout.reconfigure(encoding='utf-8')
from openai import OpenAI

client = OpenAI(
    api_key=os.environ["SEED_LLM_API_KEY"],
    base_url='https://uni-api.cstcloud.cn/v1'
)

# Load first 3 recent papers
with open(r'D:\科研\CTC知识图谱\v2\data\processed\abstracts.json', 'r', encoding='utf-8') as f:
    abstracts = json.load(f)

pmids = sorted(abstracts.keys(), reverse=True)[:3]

# Slim prompt for testing
SYSTEM_PROMPT = """You are a CTC biologist. Extract structured information from this abstract.
Return ONLY a JSON object with these keys:
- cancer_types: list of cancer type names mentioned
- genes: list of {symbol, role (driver/biomarker/target/resistance_marker), evidence}
- gene_phenotype_relations: list of {gene, phenotype, direction (DRIVES/SUPPRESSES/ASSOCIATED_WITH), evidence, confidence (high/medium/low)}
- gene_gene_relations: list of {gene_a, gene_b, relation, mechanism, evidence, confidence}
- gene_pathway_relations: list of {gene, pathway, relation, evidence, confidence}
- biomarkers: list of {gene, cancer_type, type (diagnostic/prognostic/predictive), evidence}
- drugs: list of {name, relation (resistance/sensitivity/targets), gene, evidence}
- ctc_subpopulations: list of {type (epithelial/mesenchymal/hybrid/platelet_covered/cluster/stem_like), evidence}
Phenotype types: EMT_induction, MET_induction, metastasis, invasion, migration, CTC_cluster_formation, immune_evasion, drug_resistance, anoikis_resistance, apoptosis_resistance, stemness, proliferation, angiogenesis, metabolic_reprogramming, senescence, dormancy
Confidence: high (KD/KO/inhibitor data), medium (correlation), low (inferred)
No markdown, no explanation, ONLY the JSON object."""

for pmid in pmids:
    data = abstracts[pmid]
    title = data.get('title', '')[:500]
    abstract = data.get('abstract', '')[:3000]
    if not abstract: continue

    t0 = time.time()
    try:
        resp = client.chat.completions.create(
            model='S1-Base-Ultra',
            messages=[
                {'role': 'system', 'content': SYSTEM_PROMPT},
                {'role': 'user', 'content': f'TITLE: {title}\n\nABSTRACT: {abstract}'}
            ],
            temperature=0,
            max_tokens=1500,
            timeout=45
        )
        content = resp.choices[0].message.content.strip()
        if '```json' in content:
            content = content.split('```json')[1].split('```')[0]
        if '{' in content:
            content = content[content.find('{'):content.rfind('}')+1]
        result = json.loads(content)
        n_genes = len(result.get('genes', []))
        n_phen = len(result.get('gene_phenotype_relations', []))
        n_gg = len(result.get('gene_gene_relations', []))
        n_bio = len(result.get('biomarkers', []))
        n_ctc = len(result.get('ctc_subpopulations', []))
        elapsed = time.time()-t0
        print(f'PMID {pmid}: {n_genes}G + {n_phen}Ph + {n_gg}GG + {n_bio}Bio + {n_ctc}Sub ({elapsed:.1f}s)')
        if result.get('cancer_types'):
            print(f'  Cancer: {result["cancer_types"]}')
        if result.get('genes'):
            print(f'  Genes: {[(g.get("symbol","?"), g.get("role","?")) for g in result["genes"][:5]]}')
    except json.JSONDecodeError as e:
        print(f'PMID {pmid}: JSON error - {str(e)[:60]} ({time.time()-t0:.1f}s)')
    except Exception as e:
        print(f'PMID {pmid}: {type(e).__name__}: {str(e)[:100]} ({time.time()-t0:.1f}s)')
