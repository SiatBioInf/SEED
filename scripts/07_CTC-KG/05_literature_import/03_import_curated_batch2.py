"""
Batch 2: 手工提取 - PAK2/DSG2/PERK/RBC/Desmosome CTC机制
"""

import sys, os
sys.path.insert(0, r'D:\科研\CTC知识图谱\v2')
sys.stdout.reconfigure(encoding='utf-8')
from neo4j import GraphDatabase

driver = GraphDatabase.driver('bolt://127.0.0.1:7687', auth=(os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"]))

BATCH2 = {
    "gene_phenotype": [
        # PMID 41423632: PAK2 promotes CTC cluster formation
        {"gene":"PAK2","phenotype":"CTC_cluster_formation","direction":"DRIVES",
         "evidence":"PAK2 promotes CTC cluster formation by increasing E-cadherin phosphorylation at Ser840",
         "confidence":"high","pmid":"41423632","cancer":"Breast Cancer"},
        {"gene":"PAK2","phenotype":"metastasis","direction":"DRIVES",
         "evidence":"Pharmacologic inhibition of PAK2 ... reduced metastasis and decreased CTC cluster numbers in vivo",
         "confidence":"high","pmid":"41423632","cancer":"Breast Cancer"},

        # PMID 41381723: HIF1A-DSG2 desmosome junction in CTC clusters
        {"gene":"HIF1A","phenotype":"CTC_cluster_formation","direction":"DRIVES",
         "evidence":"HIF-1a positively controls DSG2-mediated desmosome junctions in CTC clusters",
         "confidence":"high","pmid":"41381723","cancer":"Breast Cancer"},
        {"gene":"DSG2","phenotype":"CTC_cluster_formation","direction":"DRIVES",
         "evidence":"DSG2 mediates conversion between desmosome and adhesion junctions in CTC clusters",
         "confidence":"high","pmid":"41381723","cancer":"Breast Cancer"},
        {"gene":"DSG2","phenotype":"metastasis","direction":"DRIVES",
         "evidence":"hirudin ... disruption of CTC clusters rather than directly inhibiting cell migration",
         "confidence":"high","pmid":"41381723","cancer":"Breast Cancer"},

        # PMID 41212905: PERK/MAT2A/PDGFB axis in CTC cluster survival
        {"gene":"EIF2AK3","phenotype":"CTC_cluster_formation","direction":"DRIVES",
         "evidence":"PERK was activated in CTC clusters and maintained SAM availability, facilitating adaptive survival",
         "confidence":"high","pmid":"41212905","cancer":"Multiple"},
        {"gene":"PDGFB","phenotype":"CTC_cluster_formation","direction":"DRIVES",
         "evidence":"PDGFB functioned as a shared survival signal triggering PI3K/AKT pathway",
         "confidence":"high","pmid":"41212905","cancer":"Multiple"},
        {"gene":"EIF2AK3","phenotype":"metastasis","direction":"DRIVES",
         "evidence":"Inhibition of PERK and PDGFRbeta profoundly impaired survival signaling and suppressed metastatic dissemination",
         "confidence":"high","pmid":"41212905","cancer":"Multiple"},

        # PMID 41558870: Platelet-adherent CTCs
        {"gene":"PF4","phenotype":"EMT_induction","direction":"DRIVES",
         "evidence":"CTCs exhibited high expression of platelet-related genes such as PF4 and PPBP, suggesting platelet-driven EMT activation",
         "confidence":"medium","pmid":"41558870","cancer":"Pancreatic Cancer"},
        {"gene":"PPBP","phenotype":"EMT_induction","direction":"DRIVES",
         "evidence":"high expression of platelet-related genes such as PF4 and PPBP, suggesting platelet-driven EMT activation",
         "confidence":"medium","pmid":"41558870","cancer":"Pancreatic Cancer"},

        # PMID 41275332: G-CSF drives pre-metastatic niche
        {"gene":"CSF3","phenotype":"metastasis","direction":"DRIVES",
         "evidence":"tumor-derived G-CSF modulates the pre-metastatic niche through activating VEGFA-KDR signaling axis",
         "confidence":"high","pmid":"41275332","cancer":"Breast Cancer"},
        {"gene":"CSF3","phenotype":"angiogenesis","direction":"DRIVES",
         "evidence":"G-CSF ... promoting angiogenesis and PMN development",
         "confidence":"high","pmid":"41275332","cancer":"Breast Cancer"},
    ],

    "gene_gene": [
        # PMID 41423632: PAK2 → CDH1
        {"gene_a":"PAK2","gene_b":"CDH1","relation":"PHOSPHORYLATES",
         "mechanism":"PAK2 phosphorylates E-cadherin at Ser840 increasing cell-cell adhesion",
         "evidence":"PAK2 promotes CTC cluster formation by increasing E-cadherin phosphorylation at Ser840",
         "confidence":"high","pmid":"41423632"},

        # PMID 41381723: HIF1A → DSG2
        {"gene_a":"HIF1A","gene_b":"DSG2","relation":"UPREGULATES",
         "mechanism":"HIF-1a transcriptionally activates DSG2 expression in hypoxic CTC clusters",
         "evidence":"HIF-1a positively controls DSG2-mediated desmosome junctions",
         "confidence":"high","pmid":"41381723"},

        # PMID 41212905: PERK→ATF4→MAT2A cascade
        {"gene_a":"EIF2AK3","gene_b":"ATF4","relation":"UPREGULATES",
         "mechanism":"PERK mediates ATF4 upregulation as part of UPR",
         "evidence":"PERK mediated the upregulation of ATF4, which enhanced MAT2A expression",
         "confidence":"high","pmid":"41212905"},
        {"gene_a":"ATF4","gene_b":"MAT2A","relation":"UPREGULATES",
         "mechanism":"ATF4 transcriptionally activates MAT2A for SAM synthesis",
         "evidence":"ATF4 ... enhanced MAT2A expression, contributing to SAM synthesis",
         "confidence":"high","pmid":"41212905"},
        {"gene_a":"PDGFB","gene_b":"PDGFRB","relation":"ACTIVATES",
         "mechanism":"PDGFB binds PDGFRbeta triggering PI3K/AKT survival signaling",
         "evidence":"PDGFB ... triggering PI3K/AKT pathway via PDGFRbeta",
         "confidence":"high","pmid":"41212905"},

        # PMID 41275332: G-CSF → KDR
        {"gene_a":"CSF3","gene_b":"KDR","relation":"UPREGULATES",
         "mechanism":"G-CSF activates VEGFA-KDR signaling axis in lung endothelium",
         "evidence":"G-CSF ... activating the VEGFA-KDR signaling axis in the lung",
         "confidence":"high","pmid":"41275332"},
    ],

    "gene_pathway": [
        {"gene":"EIF2AK3","pathway":"UNFOLDED_PROTEIN_RESPONSE","relation":"ACTIVATES",
         "evidence":"PERK, unfolded protein response mediator, activated in CTC clusters","confidence":"high","pmid":"41212905"},
        {"gene":"PDGFB","pathway":"PI3K_AKT_MTOR_SIGNALING","relation":"ACTIVATES",
         "evidence":"PDGFB triggering PI3K/AKT pathway via PDGFRbeta","confidence":"high","pmid":"41212905"},
        {"gene":"PAK2","pathway":"EPITHELIAL_MESENCHYMAL_TRANSITION","relation":"ASSOCIATED_WITH",
         "evidence":"PAK2 plays critical role in cytoskeletal remodeling, advanced tumor progression","confidence":"medium","pmid":"41423632"},
        {"gene":"CSF3","pathway":"ANGIOGENESIS","relation":"ACTIVATES",
         "evidence":"G-CSF-KDR axis promoting angiogenesis and PMN development","confidence":"high","pmid":"41275332"},
    ],

    "biomarkers": [
        {"gene":"PAK2","cancer_type":"Breast Cancer","type":"prognostic",
         "evidence":"PAK2 elevated expression associated with poor patient prognosis","pmid":"41423632"},
        {"gene":"ASCL1","cancer_type":"SCLC","type":"diagnostic",
         "evidence":"CyTOF identified CTC populations of SCLC subtypes driven by ASCL1, NeuroD1, POU2F3","pmid":"41446220"},
    ],

    "drugs": [
        {"name":"FRAX597","relation":"targets","gene":"PAK2",
         "evidence":"FRAX597 impaired CTC cluster formation, suppressed tumor growth, reduced metastasis","pmid":"41423632","cancer":"Breast Cancer"},
        {"name":"Hirudin","relation":"targets","gene":"HIF1A",
         "evidence":"Hirudin inhibits HIF-1a-controlled DSG2-mediated desmosome junctions, disintegrating CTC clusters","pmid":"41381723","cancer":"Breast Cancer"},
        {"name":"Cabazitaxel","relation":"resistance","gene":"CIN",
         "evidence":"CTC chromosomal instability (CTC-CIN) can predict taxane resistance in mCRPC","pmid":"41186991","cancer":"Prostate Cancer"},
    ],

    "ctc_subpopulations": [
        {"gene":"PAK2","subpop_type":"CTC Cluster","evidence":"PAK2 promotes CTC cluster formation via E-cadherin phosphorylation","pmid":"41423632"},
        {"gene":"DSG2","subpop_type":"CTC Cluster","evidence":"DSG2 mediates desmosome junctions maintaining CTC cluster stability","pmid":"41381723"},
        {"gene":"HIF1A","subpop_type":"CTC Cluster","evidence":"HIF-1a controls DSG2-mediated desmosome junctions in CTC clusters","pmid":"41381723"},
        {"gene":"EIF2AK3","subpop_type":"CTC Cluster","evidence":"PERK activated in CTC clusters, supports adaptive survival in bloodstream","pmid":"41212905"},
        {"gene":"PF4","subpop_type":"Platelet-Covered CTC","evidence":"CTCs exhibited high expression of platelet-related genes PF4 and PPBP","pmid":"41558870"},
        {"gene":"PPBP","subpop_type":"Platelet-Covered CTC","evidence":"PF4 and PPBP, suggesting platelet-driven EMT activation","pmid":"41558870"},
        {"gene":"ASCL1","subpop_type":"Stem-like CTC","evidence":"CyTOF identified SCLC CTC subtypes driven by ASCL1 transcription factor","pmid":"41446220"},
    ],
}


def check_gene(gene):
    with driver.session() as s:
        node = s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=gene.upper()).single()
        if node:
            expr = s.run("MATCH (g:Gene {symbol: $s})-[r:CTC_EXPRESSED_IN]->(c:CancerType) RETURN c.name AS c, r.specificity_score AS sc ORDER BY sc DESC", s=gene.upper())
            return True, [{"cancer":r["c"],"score":r["sc"]} for r in expr]
        return False, []

def import_batch():
    print("="*60)
    print("Batch 2: Importing CTC Cluster & Survival Mechanisms")
    print("="*60)
    stats = {"gp":0,"gg":0,"gpw":0,"bio":0,"drug":0,"ctc":0,"skip":0}

    with driver.session() as s:
        # Gene-Phenotype
        print("\n--- Gene-Phenotype ---")
        for r in BATCH2["gene_phenotype"]:
            gene = r["gene"].upper()
            ok, expr = check_gene(gene)
            if not ok:
                print(f"  SKIP {gene}: not in KG"); stats["skip"]+=1; continue
            display = r["phenotype"].replace("_"," ").title()
            s.run("MERGE (ph:CTC_Phenotype {name: $n})", n=display)
            cm = {"high":0.9,"medium":0.6,"low":0.3}
            s.run("""MATCH (g:Gene {symbol:$g}) MATCH (ph:CTC_Phenotype {name:$ph})
                MERGE (g)-[rel:REGULATES {pmid:$pm}]->(ph)
                SET rel.direction=$d, rel.mechanism=$e, rel.confidence=$c,
                    rel.evidence_type='hand_curated', rel.cancer_context=$ca, rel.pmid=$pm
            """, g=gene, ph=display, d=r["direction"], e=r["evidence"][:300],
                 c=cm.get(r["confidence"],0.5), pm=r["pmid"], ca=r.get("cancer",""))
            stats["gp"]+=1
            print(f"  {gene} -[{r['direction']}]-> {display} ({r['pmid']})")

        # Gene-Gene
        print("\n--- Gene-Gene ---")
        for r in BATCH2["gene_gene"]:
            ga,gb = r["gene_a"].upper(), r["gene_b"].upper()
            oa,_ = check_gene(ga); ob,_ = check_gene(gb)
            if not oa or not ob: print(f"  SKIP {ga}->{gb}"); stats["skip"]+=1; continue
            cm = {"high":0.9,"medium":0.6,"low":0.3}
            ppi = s.run("MATCH (a:Gene {symbol:$a})-[r:INTERACTS_WITH]-(b:Gene {symbol:$b}) RETURN r.score AS sc LIMIT 1", a=ga,b=gb).single()
            ss = ppi["sc"] if ppi else None
            s.run("""MATCH (a:Gene {symbol:$a}) MATCH (b:Gene {symbol:$b})
                MERGE (a)-[rel:REGULATES {pmid:$pm}]->(b)
                SET rel.relation_type=$rt, rel.mechanism=$m, rel.confidence=$c,
                    rel.evidence_type='hand_curated', rel.string_validated_score=$ss, rel.pmid=$pm
            """, a=ga,b=gb,rt=r["relation"],m=r.get("mechanism","")[:300],
                 c=cm.get(r["confidence"],0.5),pm=r["pmid"],ss=ss)
            stats["gg"]+=1
            tag = f" [STRING={ss}]" if ss else ""
            print(f"  {ga} -[{r['relation']}]-> {gb} ({r['pmid']}){tag}")

        # Gene-Pathway
        print("\n--- Gene-Pathway ---")
        for r in BATCH2["gene_pathway"]:
            gene = r["gene"].upper()
            ok,_ = check_gene(gene)
            if not ok: print(f"  SKIP {gene}"); continue
            pw = s.run("""MATCH (p:Pathway) WHERE toLower(p.name) CONTAINS toLower($pw) RETURN p.name AS n LIMIT 1""", pw=r["pathway"]).single()
            if pw:
                cm = {"high":0.9,"medium":0.6,"low":0.3}
                s.run("""MATCH (g:Gene {symbol:$g}) MATCH (p:Pathway {name:$pn})
                    MERGE (g)-[rel:REGULATES {pmid:$pm}]->(p)
                    SET rel.relation_type=$rt, rel.mechanism=$e, rel.confidence=$c,
                        rel.evidence_type='hand_curated', rel.pmid=$pm
                """, g=gene, pn=pw["n"], rt=r["relation"], e=r["evidence"][:300],
                     c=cm.get(r["confidence"],0.5), pm=r["pmid"])
                stats["gpw"]+=1
                print(f"  {gene} -[{r['relation']}]-> {pw['n']} ({r['pmid']})")
            else:
                print(f"  SKIP pathway '{r['pathway']}': not in KG")

        # Biomarkers
        print("\n--- Biomarkers ---")
        for r in BATCH2["biomarkers"]:
            gene = r["gene"].upper()
            ok,_ = check_gene(gene)
            if not ok: print(f"  SKIP {gene}"); continue
            cn = s.run("""MATCH (c:CancerType) WHERE toLower(c.name) CONTAINS toLower($cn) RETURN c.name AS n LIMIT 1""", cn=r["cancer_type"]).single()
            if cn:
                s.run("""MATCH (g:Gene {symbol:$g}) MATCH (c:CancerType {name:$cn})
                    MERGE (g)-[rel:BIOMARKER_FOR {pmid:$pm}]->(c)
                    SET rel.biomarker_type=$bt, rel.evidence=$e, rel.evidence_type='hand_curated'
                """, g=gene, cn=cn["n"], bt=r["type"], e=r["evidence"][:300], pm=r["pmid"])
                stats["bio"]+=1
                print(f"  {gene} -[BIOMARKER_FOR {r['type']}]-> {cn['n']} ({r['pmid']})")

        # Drugs
        print("\n--- Drugs ---")
        for r in BATCH2["drugs"]:
            gene = r["gene"].upper()
            ok,_ = check_gene(gene) if gene not in ["CIN"] else (True,[])
            if not ok: print(f"  SKIP {gene}"); continue
            s.run("MERGE (d:Drug {name:$n})", n=r["name"])
            edge_type = {"resistance":"DRUG_RESISTANCE_IN","sensitivity":"DRUG_SENSITIVITY_IN","targets":"TARGETS"}.get(r["relation"],"TARGETS")
            if edge_type == "TARGETS":
                s.run("""MATCH (d:Drug {name:$dn}) MATCH (g:Gene {symbol:$g})
                    MERGE (d)-[rel:TARGETS {pmid:$pm}]->(g)
                    SET rel.evidence=$e, rel.evidence_type='hand_curated', rel.cancer_context=$ca
                """, dn=r["name"],g=gene,e=r["evidence"][:300],pm=r["pmid"],ca=r.get("cancer",""))
            else:
                s.run(f"""MATCH (g:Gene {{symbol:$g}}) MATCH (d:Drug {{name:$dn}})
                    MERGE (g)-[rel:{edge_type} {{pmid:$pm}}]->(d)
                    SET rel.evidence=$e, rel.evidence_type='hand_curated', rel.cancer_context=$ca
                """, g=gene,dn=r["name"],e=r["evidence"][:300],pm=r["pmid"],ca=r.get("cancer",""))
            stats["drug"]+=1
            print(f"  {gene} -[{r['relation']}]-> {r['name']} ({r['pmid']})")

        # CTC Subpopulations
        print("\n--- CTC Subpopulations ---")
        for r in BATCH2["ctc_subpopulations"]:
            gene = r["gene"].upper()
            ok,_ = check_gene(gene)
            if not ok: print(f"  SKIP {gene}"); continue
            s.run("MERGE (csp:CTC_Subpopulation {name:$n})", n=r["subpop_type"])
            s.run("""MATCH (g:Gene {symbol:$g}) MATCH (csp:CTC_Subpopulation {name:$sp})
                MERGE (g)-[rel:EPITHELIAL_MESENCHYMAL_STATE {pmid:$pm}]->(csp)
                SET rel.evidence=$e, rel.evidence_type='hand_curated'
            """, g=gene,sp=r["subpop_type"],e=r["evidence"][:300],pm=r["pmid"])
            stats["ctc"]+=1
            print(f"  {gene} -> {r['subpop_type']} ({r['pmid']})")

    print(f"\nImport: gp={stats['gp']} gg={stats['gg']} gpw={stats['gpw']} bio={stats['bio']} drug={stats['drug']} ctc={stats['ctc']} skip={stats['skip']}")
    with driver.session() as s:
        for e in ["REGULATES","BIOMARKER_FOR","TARGETS","DRUG_RESISTANCE_IN","DRUG_SENSITIVITY_IN","EPITHELIAL_MESENCHYMAL_STATE"]:
            c = s.run(f"MATCH ()-[r:{e}]->() RETURN count(r) AS n").single()
            print(f"  Total {e}: {c['n']}")

if __name__ == "__main__":
    import_batch()
    driver.close()
    print("Done!")
