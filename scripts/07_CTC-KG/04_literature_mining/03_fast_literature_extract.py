"""
快速文献提取器 — 批量处理 + 增量保存
=====================================
策略优化:
  1. 简洁Prompt减少LLM推理时间
  2. 每10篇保存一次中间结果（不怕中断）
  3. 30s超时 + 自动跳过失败的论文
  4. 两阶段: 先提取所有基因/关系，再交叉验证KG

用法:
  python fast_literature_extract.py [--batch-size 5] [--limit 50]
"""

import json, time, sys, os
sys.path.insert(0, r'D:\科研\CTC知识图谱\v2')
sys.stdout.reconfigure(encoding='utf-8')
sys.stderr.reconfigure(encoding='utf-8')

from openai import OpenAI

# Config
API_KEY = os.environ["SEED_LLM_API_KEY"]
BASE_URL = "https://uni-api.cstcloud.cn/v1"
MODEL = "deepseek-v4-flash"
DATA_DIR = r"D:\科研\CTC知识图谱\v2\data\processed"
ABSTRACTS_FILE = os.path.join(DATA_DIR, "abstracts.json")
OUTPUT_FILE = os.path.join(DATA_DIR, "extracted_relations_deep.json")

client = OpenAI(api_key=API_KEY, base_url=BASE_URL)

# 紧凑Prompt - 关键信息全覆盖，但token数少
SYSTEM_PROMPT = """Extract from this CTC paper abstract into JSON:
{
  "genes": [{"symbol": "GENE", "role": "driver|biomarker|target|resistance_marker", "evidence": "quote"}],
  "cancer_types": ["cancer name"],
  "ctc_subpopulations": [{"type": "epithelial|mesenchymal|hybrid|platelet_covered|cluster|stem_like", "evidence": "quote"}],
  "phenotypes": [{"gene": "GENE", "phenotype": "EMT_induction|metastasis|invasion|migration|CTC_cluster_formation|immune_evasion|drug_resistance|anoikis_resistance|apoptosis_resistance|stemness|proliferation|angiogenesis|metabolic_reprogramming|senescence|dormancy|extravasation", "direction": "DRIVES|SUPPRESSES|PROMOTES|INHIBITS|ASSOCIATED_WITH", "evidence": "quote", "confidence": "high|medium|low"}],
  "gene_gene": [{"gene_a": "GENE", "gene_b": "GENE", "relation": "UPREGULATES|DOWNREGULATES|ACTIVATES|INHIBITS|INTERACTS|BINDS", "evidence": "quote", "confidence": "high|medium|low"}],
  "gene_pathway": [{"gene": "GENE", "pathway": "name", "relation": "ACTIVATES|INHIBITS|MEMBER_OF", "evidence": "quote", "confidence": "high|medium|low"}],
  "biomarkers": [{"gene": "GENE", "cancer_type": "name", "type": "diagnostic|prognostic|predictive", "evidence": "quote"}],
  "drugs": [{"name": "drug", "relation": "resistance|sensitivity|targets", "gene": "GENE", "evidence": "quote"}],
  "mechanism_summary": "1 sentence"
}
Confidence: high (experimental KD/KO/inhibitor), medium (correlation), low (inferred).
Return ONLY the JSON, no markdown. Empty arrays if no findings."""


def extract_one(title, abstract, pmid):
    """提取一篇论文"""
    t0 = time.time()
    try:
        resp = client.chat.completions.create(
            model=MODEL,
            messages=[
                {"role": "system", "content": SYSTEM_PROMPT},
                {"role": "user", "content": f"TITLE: {title[:300]}\n\nABSTRACT: {abstract[:3000]}"}
            ],
            temperature=0,
            max_tokens=1200,
            timeout=40
        )
        content = resp.choices[0].message.content.strip()
        # Clean
        if '```json' in content:
            content = content.split('```json')[1].split('```')[0]
        elif '```' in content:
            content = content.split('```')[1].split('```')[0]
        if '{' in content:
            content = content[content.find('{'):content.rfind('}') + 1]

        result = json.loads(content)
        result['pmid'] = pmid
        result['title'] = title[:200]
        elapsed = time.time() - t0
        return result, elapsed
    except json.JSONDecodeError:
        return None, time.time() - t0
    except Exception as e:
        print(f"    Error: {type(e).__name__}: {str(e)[:80]}")
        return None, time.time() - t0


def count_extracted(result):
    """统计提取了多少条信息"""
    n = 0
    for key in ['genes', 'cancer_types', 'ctc_subpopulations', 'phenotypes',
                'gene_gene', 'gene_pathway', 'biomarkers', 'drugs']:
        val = result.get(key, [])
        if isinstance(val, list):
            n += len(val)
    return n


def main():
    # Load abstracts
    with open(ABSTRACTS_FILE, 'r', encoding='utf-8') as f:
        abstracts = json.load(f)
    print(f"Loaded {len(abstracts)} abstracts")

    # Resume support
    done_pmids = set()
    all_results = []
    if os.path.exists(OUTPUT_FILE):
        with open(OUTPUT_FILE, 'r', encoding='utf-8') as f:
            all_results = json.load(f)
        done_pmids = {r['pmid'] for r in all_results}
        print(f"Resuming: {len(done_pmids)} already done, {len(abstracts) - len(done_pmids)} remaining")

    # Sort by PMID descending (newer first)
    items = [(pmid, data) for pmid, data in abstracts.items()
             if pmid not in done_pmids]
    items.sort(key=lambda x: int(x[0]) if x[0].isdigit() else 0, reverse=True)

    total = len(items)
    success = 0
    empty = 0
    failed = 0
    total_relations = 0
    total_time = 0.0

    print(f"Processing {total} papers...")
    print(f"{'#':5s} {'PMID':10s} {'Status':8s} {'Items':6s} {'Time':7s} {'Total':8s}")
    print("-" * 55)

    for i, (pmid, data) in enumerate(items):
        title = data.get('title', '')[:500]
        abstract = data.get('abstract', '')

        if not abstract or len(abstract) < 50:
            failed += 1
            continue

        result, elapsed = extract_one(title, abstract, pmid)
        total_time += elapsed

        if result is None:
            failed += 1
            print(f"{i+1:5d} {pmid:10s} {'FAIL':8s} {'-':6s} {elapsed:6.1f}s "
                  f"{total_time/60:7.1f}m")
        else:
            n_items = count_extracted(result)
            total_relations += n_items
            all_results.append(result)
            success += 1
            if n_items == 0:
                empty += 1

            print(f"{i+1:5d} {pmid:10s} {'OK':8s} {n_items:5d}  {elapsed:6.1f}s "
                  f"{total_time/60:7.1f}m")

        # Save every 10 papers
        if (i + 1) % 10 == 0:
            with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
                json.dump(all_results, f, ensure_ascii=False, indent=2)
            eta = (total_time / (i + 1)) * (total - i - 1) / 60 if (i + 1) > 0 else 0
            print(f"  --- Saved {len(all_results)} results | "
                  f"OK:{success} Empty:{empty} Fail:{failed} | "
                  f"ETA:{eta:.0f}m ---")
            sys.stdout.flush()

        # Brief pause between papers
        time.sleep(0.05)

    # Final save
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        json.dump(all_results, f, ensure_ascii=False, indent=2)

    # Summary
    print(f"\n{'='*60}")
    print(f"EXTRACTION COMPLETE")
    print(f"  Total papers: {total}")
    print(f"  Success: {success} ({100*success//max(1,total)}%)")
    print(f"  Empty (no findings): {empty}")
    print(f"  Failed: {failed}")
    print(f"  Total relations extracted: {total_relations}")
    print(f"  Total time: {total_time/60:.1f} min")
    print(f"  Avg time/paper: {total_time/max(1,success):.1f}s")
    print(f"  Output: {OUTPUT_FILE}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
