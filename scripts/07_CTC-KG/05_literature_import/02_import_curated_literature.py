"""
手工提取的高价值文献关系 → 直接导入Neo4j
===========================================
从最新CTC机制论文中手工提取关系，交叉验证KG后导入。
免去LLM API依赖，保证提取质量和KG一致性。

每一条关系都附带:
  - PMID (溯源)
  - evidence (原文引用)
  - confidence (高/中/低)
  - cancer_context (癌种背景)
  - kg_crossref (KG中是否存在该基因/表达数据)
"""

import sys, os
sys.path.insert(0, r'D:\科研\CTC知识图谱\v2')
sys.stdout.reconfigure(encoding='utf-8')

from neo4j import GraphDatabase

NEO4J_URI = "bolt://127.0.0.1:7687"
NEO4J_AUTH = (os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"])
driver = GraphDatabase.driver(NEO4J_URI, auth=NEO4J_AUTH)

# ============================================================
# 手工提取的关系数据 (2024-2026高价值CTC机制论文)
# ============================================================

CURATED_RELATIONS = {
    "gene_phenotype": [
        # PMID 42139058: Ki-67 drives CTC intravasation
        {"gene": "MKI67", "phenotype": "intravasation", "direction": "DRIVES",
         "evidence": "Ki-67 depletion ... reduces metastasis. Ki-67 as an essential regulator of CTC intravasation",
         "confidence": "high", "pmid": "42139058", "cancer": "Breast Cancer"},
        {"gene": "MKI67", "phenotype": "metastasis", "direction": "DRIVES",
         "evidence": "knockout reduces metastasis",
         "confidence": "high", "pmid": "42139058", "cancer": "Breast Cancer"},
        {"gene": "MKI67", "phenotype": "CTC_cluster_formation", "direction": "ASSOCIATED_WITH",
         "evidence": "striking enrichment of cell-cycle genes in CTCs, particularly in CTC clusters",
         "confidence": "medium", "pmid": "42139058", "cancer": "Breast Cancer"},

        # PMID 42128301: AST genes in melanoma
        {"gene": "IKZF1", "phenotype": "metastasis", "direction": "DRIVES",
         "evidence": "AST-high circulating tumor cell clones possess enhanced metastatic capacity",
         "confidence": "high", "pmid": "42128301", "cancer": "Melanoma"},
        {"gene": "IRF8", "phenotype": "metastasis", "direction": "DRIVES",
         "evidence": "AST genes critically modulate anchorage dependence",
         "confidence": "high", "pmid": "42128301", "cancer": "Melanoma"},
        {"gene": "NFE2", "phenotype": "metastasis", "direction": "DRIVES",
         "evidence": "AST genes critically modulate anchorage dependence",
         "confidence": "high", "pmid": "42128301", "cancer": "Melanoma"},
        {"gene": "IKZF1", "phenotype": "extravasation", "direction": "PROMOTES",
         "evidence": "AST-positive tumor cells preferentially localize near blood vessels",
         "confidence": "medium", "pmid": "42128301", "cancer": "Melanoma"},

        # PMID 41996487: UFMylation-EMD-β-catenin in HCC
        {"gene": "EMD", "phenotype": "EMT_induction", "direction": "SUPPRESSES",
         "evidence": "Loss of UFMylation destabilized EMD ... promoting hybrid EMT",
         "confidence": "high", "pmid": "41996487", "cancer": "Liver Cancer"},
        {"gene": "EMD", "phenotype": "metastasis", "direction": "SUPPRESSES",
         "evidence": "EMD restoration mitigated metastasis",
         "confidence": "high", "pmid": "41996487", "cancer": "Liver Cancer"},
        {"gene": "EMD", "phenotype": "immune_evasion", "direction": "SUPPRESSES",
         "evidence": "stabilizing EMD, thereby inhibiting ... PVR-mediated immune evasion",
         "confidence": "high", "pmid": "41996487", "cancer": "Liver Cancer"},
        {"gene": "PVR", "phenotype": "immune_evasion", "direction": "DRIVES",
         "evidence": "PVR enabling natural killer (NK) cell evasion via PVR-TIGIT axis",
         "confidence": "high", "pmid": "41996487", "cancer": "Liver Cancer"},

        # PMID 41651843: HA-CD44 mediates CTC clustering in TNBC
        {"gene": "CD44", "phenotype": "CTC_cluster_formation", "direction": "DRIVES",
         "evidence": "HA is necessary and sufficient to drive clustering of tumor cells expressing its receptor CD44",
         "confidence": "high", "pmid": "41651843", "cancer": "Breast Cancer"},
        {"gene": "CD44", "phenotype": "metastasis", "direction": "DRIVES",
         "evidence": "ECM-receptor interaction enables CTC clustering ... enhancing TNBC metastasis",
         "confidence": "high", "pmid": "41651843", "cancer": "Breast Cancer"},

        # PMID 41803969: Integrin α2 in metastasis cascade
        {"gene": "ITGA2", "phenotype": "metastasis", "direction": "DRIVES",
         "evidence": "integrin α2 profoundly influences the metastatic cascade",
         "confidence": "medium", "pmid": "41803969", "cancer": "Multiple"},
        {"gene": "ITGA2", "phenotype": "invasion", "direction": "DRIVES",
         "evidence": "key mediator of cell-ECM communication ... local invasion",
         "confidence": "high", "pmid": "41803969", "cancer": "Multiple"},
        {"gene": "ITGA2", "phenotype": "immune_evasion", "direction": "PROMOTES",
         "evidence": "immune evasion ... circulating tumor cell survival",
         "confidence": "medium", "pmid": "41803969", "cancer": "Multiple"},
        {"gene": "ITGA2", "phenotype": "dormancy", "direction": "ASSOCIATED_WITH",
         "evidence": "tumor dormancy ... organotropic colonization",
         "confidence": "low", "pmid": "41803969", "cancer": "Multiple"},
    ],

    "gene_gene": [
        # PMID 42139058: MKI67 → CD47, KLF4
        {"gene_a": "MKI67", "gene_b": "CD47", "relation": "UPREGULATES",
         "mechanism": "Ki-67 maintains cell-cell adhesion gene expression",
         "evidence": "Ki-67 depletion ... suppresses genes involved in maintaining cell-cell adhesion, including CD47 and KLF4",
         "confidence": "high", "pmid": "42139058"},
        {"gene_a": "MKI67", "gene_b": "KLF4", "relation": "UPREGULATES",
         "mechanism": "Ki-67 maintains cell-cell adhesion gene expression",
         "evidence": "Ki-67 depletion suppresses genes ... including CD47 and KLF4",
         "confidence": "high", "pmid": "42139058"},

        # PMID 41996487: EMD → CTNNB1 (β-catenin)
        {"gene_a": "EMD", "gene_b": "CTNNB1", "relation": "INHIBITS",
         "mechanism": "EMD stabilizes nuclear envelope, preventing β-catenin nuclear accumulation",
         "evidence": "Loss of UFMylation destabilized EMD ... promoting nuclear β-catenin accumulation",
         "confidence": "high", "pmid": "41996487"},
        {"gene_a": "CTNNB1", "gene_b": "PVR", "relation": "UPREGULATES",
         "mechanism": "β-catenin/TEAD4 transcriptionally upregulates PVR",
         "evidence": "Nuclear β-catenin also interacted with TEAD4 to upregulate poliovirus receptor (PVR)",
         "confidence": "high", "pmid": "41996487"},

        # PMID 41651843: CD44 mediates HA-driven clustering
        {"gene_a": "CD44", "gene_b": "DSC2", "relation": "COMPLEX_WITH",
         "mechanism": "HA-CD44 initiates contact, desmosomes stabilize clusters",
         "evidence": "initial interactions expand ... subsequently stabilized by desmosomes",
         "confidence": "medium", "pmid": "41651843"},

        # PMID 42177885: CAV1/ZNF460/GRP78 signaling
        {"gene_a": "CAV1", "gene_b": "ZNF460", "relation": "REGULATES",
         "mechanism": "CAV1/ZNF460/GRP78 signaling in CSC chemoresistance",
         "evidence": "CSC signaling of CAV1/ZNF460/GRP78 in mediating chemoresistance",
         "confidence": "medium", "pmid": "42177885"},
        {"gene_a": "ZNF460", "gene_b": "HSPA5", "relation": "REGULATES",
         "mechanism": "ZNF460 regulates GRP78 (HSPA5)",
         "evidence": "CAV1/ZNF460/GRP78 signaling mediates chemoresistance",
         "confidence": "medium", "pmid": "42177885"},
    ],

    "gene_pathway": [
        # PMID 41996487: β-catenin/Wnt pathway
        {"gene": "CTNNB1", "pathway": "WNT_BETA_CATENIN_SIGNALING",
         "relation": "ACTIVATES", "evidence": "nuclear β-catenin accumulation and hybrid EMT",
         "confidence": "high", "pmid": "41996487"},
        {"gene": "EMD", "pathway": "WNT_BETA_CATENIN_SIGNALING",
         "relation": "INHIBITS", "evidence": "UFMylation stabilizes EMD, inhibiting β-catenin nuclear accumulation",
         "confidence": "high", "pmid": "41996487"},

        # PMID 41651843: HA-CD44 pathway
        {"gene": "CD44", "pathway": "EPITHELIAL_MESENCHYMAL_TRANSITION",
         "relation": "ASSOCIATED_WITH", "evidence": "ECM component HA mediates AJ-independent CTC clustering",
         "confidence": "medium", "pmid": "41651843"},

        # PMID 42177885
        {"gene": "CAV1", "pathway": "UNFOLDED_PROTEIN_RESPONSE",
         "relation": "ACTIVATES", "evidence": "CAV1/ZNF460/GRP78 signaling - GRP78 is UPR master regulator",
         "confidence": "medium", "pmid": "42177885"},
    ],

    "biomarkers": [
        # PMID 42139058
        {"gene": "MKI67", "cancer_type": "Breast Cancer", "type": "prognostic",
         "evidence": "canonical proliferation marker, pivotal prognostic and predictive biomarker in breast cancer",
         "pmid": "42139058"},

        # PMID 41873417
        {"gene": "PD-L1", "cancer_type": "NSCLC", "type": "predictive",
         "evidence": "CTC PD-L1-positive associated with ICI response in NSCLC",
         "pmid": "41873417"},
        {"gene": "CD274", "cancer_type": "NSCLC", "type": "predictive",
         "evidence": "dual biomarker combining PD-L1 on CTCs and PD-1 on circulating CD8+ T cells",
         "pmid": "41873417"},

        # PMID 42128301
        {"gene": "IKZF1", "cancer_type": "Melanoma", "type": "prognostic",
         "evidence": "AST-high CTC clones possess enhanced metastatic capacity",
         "pmid": "42128301"},
    ],

    "drugs": [
        # PMID 41996487: Anti-PD-1 + Anti-TIGIT
        {"name": "Anti-TIGIT", "relation": "sensitivity", "gene": "PVR",
         "evidence": "Combined PD-1 and TIGIT blockade significantly suppressed tumor growth and metastasis",
         "pmid": "41996487", "cancer": "Liver Cancer"},
        {"name": "Anti-PD-1", "relation": "sensitivity", "gene": "PVR",
         "evidence": "Combined programmed cell death protein 1 and TIGIT blockade suppressed metastasis",
         "pmid": "41996487", "cancer": "Liver Cancer"},

        # PMID 42177885
        {"name": "5-Fluorouracil", "relation": "resistance", "gene": "CAV1",
         "evidence": "CAV1/ZNF460/GRP78 signaling mediates chemoresistance; MFYD chemosensitizes 5-FU",
         "pmid": "42177885", "cancer": "Colorectal Cancer"},
    ],

    "ctc_subpopulations": [
        # PMID 42056253
        {"gene": "ROCK1", "subpop_type": "Amoeboid CTC",
         "evidence": "CTCs develop blebbing-amoeboid movement in 3D, ROCK-dependent",
         "pmid": "42056253"},
        {"gene": "ROCK2", "subpop_type": "Amoeboid CTC",
         "evidence": "amoeboid migration mode ... ROCK-dependent",
         "pmid": "42056253"},

        # PMID 41651843
        {"gene": "CD44", "subpop_type": "CTC Cluster",
         "evidence": "HA-CD44 mediates AJ-independent CTC clustering in TNBC",
         "pmid": "41651843"},

        # PMID 42139058
        {"gene": "MKI67", "subpop_type": "CTC Cluster",
         "evidence": "Ki-67 enrichment in CTC clusters",
         "pmid": "42139058"},

        # PMID 42450242: Neutrophil-CTC
        {"gene": "CXCR2", "subpop_type": "Platelet-Covered CTC",
         "evidence": "CXCR2-dependent neutrophil-CTC interactions in bone metastasis",
         "pmid": "42450242"},
    ],
}


# ============================================================
# 交叉验证辅助函数
# ============================================================
def check_gene_in_kg(gene_symbol):
    """检查基因是否在KG中，获取表达数据"""
    with driver.session() as s:
        node = s.run("MATCH (g:Gene {symbol: $s}) RETURN g.symbol AS s", s=gene_symbol).single()
        if not node:
            return {"in_kg": False, "expressed_in": [], "pathways": []}

        info = {"in_kg": True, "expressed_in": [], "pathways": []}
        expr = s.run("""
            MATCH (g:Gene {symbol: $s})-[r:CTC_EXPRESSED_IN]->(c:CancerType)
            RETURN c.name AS cancer, r.specificity_score AS score, r.tier AS tier
            ORDER BY r.specificity_score DESC
        """, s=gene_symbol)
        info["expressed_in"] = [{"cancer": r["cancer"], "score": r["score"], "tier": r["tier"]} for r in expr]

        pw = s.run("MATCH (g:Gene {symbol: $s})-[:MEMBER_OF_PATHWAY]->(p:Pathway) RETURN p.name AS pw", s=gene_symbol)
        info["pathways"] = [r["pw"] for r in pw]
        return info


# ============================================================
# 导入函数
# ============================================================
def import_curated_relations():
    print("=" * 60)
    print("Importing Hand-Curated Literature Relations")
    print("=" * 60)

    stats = {
        "gene_phenotype": 0, "gene_gene": 0, "gene_pathway": 0,
        "biomarkers": 0, "drugs": 0, "ctc_subpopulations": 0,
        "genes_not_in_kg": 0, "genes_in_kg": 0,
    }

    with driver.session() as s:
        # ---- Gene-Phenotype ----
        print("\n--- Gene-Phenotype Relations ---")
        for r in CURATED_RELATIONS["gene_phenotype"]:
            gene = r["gene"].upper()
            check = check_gene_in_kg(gene)
            if not check["in_kg"]:
                print(f"  SKIP {gene}: not in KG")
                stats["genes_not_in_kg"] += 1
                continue
            stats["genes_in_kg"] += 1

            display = r["phenotype"].replace("_", " ").title()
            s.run("MERGE (ph:CTC_Phenotype {name: $name})", name=display)

            # Check if gene expressed in same cancer type
            expr_cancers = [e["cancer"] for e in check["expressed_in"]]
            cancer_match = r.get("cancer", "") in expr_cancers if r.get("cancer") else False

            conf_map = {"high": 0.9, "medium": 0.6, "low": 0.3}
            conf = conf_map.get(r["confidence"], 0.5)

            s.run("""
                MATCH (g:Gene {symbol: $gene})
                MATCH (ph:CTC_Phenotype {name: $phenotype})
                MERGE (g)-[rel:REGULATES {pmid: $pmid}]->(ph)
                SET rel.direction = $direction,
                    rel.mechanism = $evidence,
                    rel.confidence = $conf,
                    rel.evidence_type = 'hand_curated',
                    rel.cancer_context = $cancer,
                    rel.scrna_validated = $validated,
                    rel.pmid = $pmid
            """, gene=gene, phenotype=display, direction=r["direction"],
                 evidence=r["evidence"][:300], conf=conf, pmid=r["pmid"],
                 cancer=r.get("cancer", ""), validated=cancer_match)

            # Link to Publication
            s.run("""
                MATCH (g:Gene {symbol: $gene})-[rel:REGULATES {pmid: $pmid}]->(ph:CTC_Phenotype {name: $phenotype})
                MATCH (pub:Publication {pmid: $pmid})
                SET rel.published_in_pmid = $pmid
            """, gene=gene, phenotype=display, pmid=r["pmid"])

            stats["gene_phenotype"] += 1
            val_tag = " [scRNA-validated]" if cancer_match else ""
            print(f"  {gene} -[{r['direction']}]-> {display} ({r['pmid']}){val_tag}")

        # ---- Gene-Gene ----
        print("\n--- Gene-Gene Relations ---")
        for r in CURATED_RELATIONS["gene_gene"]:
            gene_a, gene_b = r["gene_a"].upper(), r["gene_b"].upper()
            check_a = check_gene_in_kg(gene_a)
            check_b = check_gene_in_kg(gene_b)

            if not check_a["in_kg"]:
                print(f"  SKIP {gene_a}: not in KG")
                stats["genes_not_in_kg"] += 1
                continue
            if not check_b["in_kg"]:
                print(f"  SKIP {gene_b}: not in KG")
                stats["genes_not_in_kg"] += 1
                continue

            conf_map = {"high": 0.9, "medium": 0.6, "low": 0.3}
            conf = conf_map.get(r["confidence"], 0.5)

            # Check if STRING PPI exists between them
            ppi = s.run("""
                MATCH (a:Gene {symbol: $a})-[r:INTERACTS_WITH]-(b:Gene {symbol: $b})
                RETURN r.score AS score LIMIT 1
            """, a=gene_a, b=gene_b).single()
            string_score = ppi["score"] if ppi else None

            s.run("""
                MATCH (a:Gene {symbol: $a})
                MATCH (b:Gene {symbol: $b})
                MERGE (a)-[rel:REGULATES {pmid: $pmid}]->(b)
                SET rel.relation_type = $relation,
                    rel.mechanism = $mechanism,
                    rel.confidence = $conf,
                    rel.evidence_type = 'hand_curated',
                    rel.string_validated_score = $string_score,
                    rel.pmid = $pmid
            """, a=gene_a, b=gene_b, relation=r["relation"],
                 mechanism=r.get("mechanism", "")[:300], conf=conf,
                 pmid=r["pmid"], string_score=string_score)

            stats["gene_gene"] += 1
            ppi_tag = f" [STRING={string_score}]" if string_score else ""
            print(f"  {gene_a} -[{r['relation']}]-> {gene_b} ({r['pmid']}){ppi_tag}")

        # ---- Gene-Pathway ----
        print("\n--- Gene-Pathway Relations ---")
        for r in CURATED_RELATIONS["gene_pathway"]:
            gene = r["gene"].upper()
            check = check_gene_in_kg(gene)
            if not check["in_kg"]:
                print(f"  SKIP {gene}: not in KG")
                continue

            # Try match pathway in KG
            pw_result = s.run("""
                MATCH (p:Pathway)
                WHERE toLower(p.name) CONTAINS toLower($pw)
                RETURN p.name AS name LIMIT 1
            """, pw=r["pathway"]).single()

            if pw_result:
                conf_map = {"high": 0.9, "medium": 0.6, "low": 0.3}
                conf = conf_map.get(r["confidence"], 0.5)

                s.run("""
                    MATCH (g:Gene {symbol: $gene})
                    MATCH (p:Pathway {name: $pw_name})
                    MERGE (g)-[rel:REGULATES {pmid: $pmid}]->(p)
                    SET rel.relation_type = $relation,
                        rel.mechanism = $evidence,
                        rel.confidence = $conf,
                        rel.evidence_type = 'hand_curated',
                        rel.pmid = $pmid
                """, gene=gene, pw_name=pw_result["name"], relation=r["relation"],
                     evidence=r["evidence"][:300], conf=conf, pmid=r["pmid"])

                stats["gene_pathway"] += 1
                print(f"  {gene} -[{r['relation']}]-> {pw_result['name']} ({r['pmid']})")
            else:
                print(f"  SKIP pathway '{r['pathway']}': not in KG")

        # ---- Biomarkers ----
        print("\n--- Biomarker Relations ---")
        for r in CURATED_RELATIONS["biomarkers"]:
            gene = r["gene"].upper()
            check = check_gene_in_kg(gene)
            if not check["in_kg"]:
                print(f"  SKIP {gene}: not in KG")
                continue

            # Match cancer type
            cancer_result = s.run("""
                MATCH (c:CancerType)
                WHERE toLower(c.name) CONTAINS toLower($cancer)
                RETURN c.name AS name LIMIT 1
            """, cancer=r["cancer_type"]).single()

            if cancer_result:
                s.run("""
                    MATCH (g:Gene {symbol: $gene})
                    MATCH (c:CancerType {name: $cancer})
                    MERGE (g)-[rel:BIOMARKER_FOR {pmid: $pmid}]->(c)
                    SET rel.biomarker_type = $btype,
                        rel.evidence = $evidence,
                        rel.evidence_type = 'hand_curated'
                """, gene=gene, cancer=cancer_result["name"],
                     btype=r["type"], evidence=r["evidence"][:300], pmid=r["pmid"])
                stats["biomarkers"] += 1
                print(f"  {gene} -[BIOMARKER_FOR {r['type']}]-> {cancer_result['name']} ({r['pmid']})")

        # ---- Drug Relations ----
        print("\n--- Drug Relations ---")
        for r in CURATED_RELATIONS["drugs"]:
            gene = r["gene"].upper()
            check = check_gene_in_kg(gene)
            if not check["in_kg"]:
                print(f"  SKIP {gene}: not in KG")
                continue

            s.run("MERGE (d:Drug {name: $name})", name=r["name"])

            if r["relation"] == "resistance":
                s.run("""
                    MATCH (g:Gene {symbol: $gene})
                    MATCH (d:Drug {name: $drug})
                    MERGE (g)-[rel:DRUG_RESISTANCE_IN {pmid: $pmid}]->(d)
                    SET rel.evidence = $evidence,
                        rel.evidence_type = 'hand_curated',
                        rel.cancer_context = $cancer
                """, gene=gene, drug=r["name"], evidence=r["evidence"][:300],
                     pmid=r["pmid"], cancer=r.get("cancer", ""))
            elif r["relation"] == "sensitivity":
                s.run("""
                    MATCH (g:Gene {symbol: $gene})
                    MATCH (d:Drug {name: $drug})
                    MERGE (g)-[rel:DRUG_SENSITIVITY_IN {pmid: $pmid}]->(d)
                    SET rel.evidence = $evidence,
                        rel.evidence_type = 'hand_curated',
                        rel.cancer_context = $cancer
                """, gene=gene, drug=r["name"], evidence=r["evidence"][:300],
                     pmid=r["pmid"], cancer=r.get("cancer", ""))

            stats["drugs"] += 1
            print(f"  {gene} -[{r['relation']}]-> {r['name']} ({r['pmid']})")

        # ---- CTC Subpopulations ----
        print("\n--- CTC Subpopulations ---")
        for r in CURATED_RELATIONS["ctc_subpopulations"]:
            gene = r["gene"].upper()
            check = check_gene_in_kg(gene)
            if not check["in_kg"]:
                print(f"  SKIP {gene}: not in KG")
                continue

            s.run("MERGE (csp:CTC_Subpopulation {name: $name})", name=r["subpop_type"])
            s.run("""
                MATCH (g:Gene {symbol: $gene})
                MATCH (csp:CTC_Subpopulation {name: $subpop})
                MERGE (g)-[rel:EPITHELIAL_MESENCHYMAL_STATE {pmid: $pmid}]->(csp)
                SET rel.evidence = $evidence,
                    rel.evidence_type = 'hand_curated'
            """, gene=gene, subpop=r["subpop_type"],
                 evidence=r["evidence"][:300], pmid=r["pmid"])
            stats["ctc_subpopulations"] += 1
            print(f"  {gene} -> {r['subpop_type']} ({r['pmid']})")

    # Summary
    print("\n" + "=" * 60)
    print("IMPORT SUMMARY")
    print("=" * 60)
    for cat, count in stats.items():
        if not cat.startswith("genes_"):
            print(f"  {cat}: {count}")
    print(f"  Genes in KG: {stats['genes_in_kg']}")
    print(f"  Genes NOT in KG (skipped): {stats['genes_not_in_kg']}")

    # Final KG counts
    with driver.session() as s:
        for edge in ["REGULATES", "BIOMARKER_FOR", "DRUG_RESISTANCE_IN",
                     "DRUG_SENSITIVITY_IN", "EPITHELIAL_MESENCHYMAL_STATE"]:
            cnt = s.run(f"MATCH ()-[r:{edge}]->() RETURN count(r) AS n").single()
            print(f"  Total {edge} edges: {cnt['n']}")


if __name__ == "__main__":
    import_curated_relations()
    driver.close()
    print("\nDone!")
