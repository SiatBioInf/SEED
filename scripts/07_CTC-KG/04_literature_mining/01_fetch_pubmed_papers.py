"""
PubMed批量抓取CTC相关文献
=========================
搜索策略：分多个主题维度搜索，覆盖CTC研究的各个方向
去重后合并到abstracts.json
"""

import json, time, urllib.request, urllib.parse, xml.etree.ElementTree as ET
import sys, os
sys.stdout.reconfigure(encoding='utf-8')

DATA_DIR = r"D:\科研\CTC知识图谱\v2\data\processed"
ABSTRACTS_FILE = os.path.join(DATA_DIR, "abstracts.json")

# ============================================================
# PubMed E-utilities API
# ============================================================
PUBMED_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def search_pubmed(query, max_results=50):
    """Search PubMed and return list of PMIDs"""
    params = urllib.parse.urlencode({
        'db': 'pubmed',
        'term': query,
        'retmax': max_results,
        'retmode': 'json',
        'sort': 'relevance',
    })
    url = f"{PUBMED_BASE}/esearch.fcgi?{params}"
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = json.loads(resp.read())
        return data.get('esearchresult', {}).get('idlist', [])
    except Exception as e:
        print(f"  Search error: {e}")
        return []


def fetch_abstracts(pmids, batch_size=50):
    """Fetch abstracts for a list of PMIDs"""
    results = {}
    for i in range(0, len(pmids), batch_size):
        batch = pmids[i:i + batch_size]
        params = urllib.parse.urlencode({
            'db': 'pubmed',
            'id': ','.join(batch),
            'rettype': 'xml',
            'retmode': 'xml',
        })
        url = f"{PUBMED_BASE}/efetch.fcgi?{params}"
        try:
            with urllib.request.urlopen(url, timeout=60) as resp:
                root = ET.fromstring(resp.read())
            for article in root.findall('.//PubmedArticle'):
                pmid = article.findtext('.//PMID', '')
                title = article.findtext('.//ArticleTitle', '')
                abstract_parts = article.findall('.//AbstractText')
                if abstract_parts:
                    abstract = ' '.join(
                        e.text or '' for e in abstract_parts if e.text
                    )
                else:
                    abstract = article.findtext('.//Abstract/AbstractText', '')
                if pmid and (title or abstract):
                    results[pmid] = {
                        'title': title[:500],
                        'abstract': abstract[:5000],
                    }
        except Exception as e:
            print(f"  Fetch error for batch {i}: {e}")
        time.sleep(0.5)  # Rate limiting
    return results


# ============================================================
# 多维度搜索策略
# ============================================================
SEARCH_QUERIES = [
    # 1. CTC mechanism & EMT
    ('CTC mechanism EMT signaling',
     '("circulating tumor cells"[MeSH Terms] OR "circulating tumor cell"[Title/Abstract]) '
     'AND (mechanism OR signaling OR pathway) AND (EMT OR "epithelial mesenchymal") '
     'AND 2024:2026[dp]', 60),

    # 2. CTC drug resistance
    ('CTC drug resistance',
     '("circulating tumor cells"[MeSH Terms] OR "circulating tumor cell"[Title/Abstract]) '
     'AND ("drug resistance" OR chemoresistance OR "therapy resistance") '
     'AND 2024:2026[dp]', 60),

    # 3. CTC immunotherapy & immune evasion
    ('CTC immunotherapy PD-L1',
     '("circulating tumor cells"[MeSH Terms] OR "circulating tumor cell"[Title/Abstract]) '
     'AND (immunotherapy OR "immune checkpoint" OR "PD-L1" OR "immune evasion") '
     'AND 2024:2026[dp]', 60),

    # 4. CTC single-cell omics
    ('CTC single-cell genomics',
     '("circulating tumor cells"[MeSH Terms] OR "circulating tumor cell"[Title/Abstract]) '
     'AND ("single-cell" OR scRNA-seq OR scATAC-seq OR "single cell RNA" '
     'OR "single cell sequencing" OR transcriptomics) '
     'AND 2023:2026[dp]', 60),

    # 5. CTC clusters & metastasis
    ('CTC cluster metastasis',
     '("circulating tumor cells"[MeSH Terms] OR "circulating tumor cell"[Title/Abstract]) '
     'AND (cluster OR "CTC cluster" OR "tumor microemboli") '
     'AND (metastasis OR colonization OR extravasation) '
     'AND 2020:2026[dp]', 60),

    # 6. CTC epigenetics
    ('CTC epigenetics methylation',
     '("circulating tumor cells"[MeSH Terms] OR "circulating tumor cell"[Title/Abstract]) '
     'AND (epigenetic OR methylation OR "histone modification" OR "noncoding RNA" OR lncRNA OR miRNA) '
     'AND 2020:2026[dp]', 60),

    # 7. CTC metabolism
    ('CTC metabolism metabolic',
     '("circulating tumor cells"[MeSH Terms] OR "circulating tumor cell"[Title/Abstract]) '
     'AND (metabolism OR metabolic OR glycolysis OR "oxidative phosphorylation" OR ferroptosis) '
     'AND 2020:2026[dp]', 50),

    # 8. CTC liquid biopsy clinical applications
    ('CTC liquid biopsy clinical prognosis',
     '("circulating tumor cells"[MeSH Terms] OR "circulating tumor cell"[Title/Abstract]) '
     'AND ("liquid biopsy" OR prognostic OR "clinical outcome" OR "survival analysis") '
     'AND 2023:2026[dp]', 50),

    # 9. CTC microenvironment & platelet/immune interactions
    ('CTC platelet neutrophil immune interaction',
     '("circulating tumor cells"[MeSH Terms] OR "circulating tumor cell"[Title/Abstract]) '
     'AND (platelet OR neutrophil OR macrophage OR NK cell OR "tumor microenvironment") '
     'AND 2020:2026[dp]', 60),

    # 10. CTC organotropism & dormancy
    ('CTC dormancy organotropism bone marrow',
     '("circulating tumor cells"[MeSH Terms] OR "circulating tumor cell"[Title/Abstract]) '
     'AND (dormancy OR "disseminated tumor cell" OR organotropism OR "bone marrow" '
     'OR "metastatic niche") AND 2020:2026[dp]', 50),
]


def main():
    # Load existing abstracts
    existing = {}
    if os.path.exists(ABSTRACTS_FILE):
        with open(ABSTRACTS_FILE, 'r', encoding='utf-8') as f:
            existing = json.load(f)
    print(f"Existing abstracts: {len(existing)}")

    all_new_pmids = set()

    for query_name, query, max_results in SEARCH_QUERIES:
        print(f"\nSearching: {query_name}...")
        pmids = search_pubmed(query, max_results)
        new_pmids = [p for p in pmids if p not in existing]
        all_new_pmids.update(new_pmids)
        print(f"  Found {len(pmids)} results, {len(new_pmids)} new")

        # Be nice to PubMed API
        time.sleep(0.5)

    print(f"\nTotal new unique PMIDs to fetch: {len(all_new_pmids)}")

    if not all_new_pmids:
        print("No new papers to add.")
        return

    # Fetch abstracts
    print("Fetching abstracts...")
    pmid_list = sorted(all_new_pmids)
    new_abstracts = fetch_abstracts(pmid_list, batch_size=50)
    print(f"  Successfully fetched: {len(new_abstracts)}")

    # Merge
    merged = dict(existing)
    merged.update(new_abstracts)

    # Save
    with open(ABSTRACTS_FILE, 'w', encoding='utf-8') as f:
        json.dump(merged, f, ensure_ascii=False, indent=2)

    print(f"\nSaved {len(merged)} abstracts ({len(new_abstracts)} new)")
    print("Done!")


if __name__ == '__main__':
    main()
