"""Batch 3: TAC1-TACR1 metastasis + Ovarian CTC genomics + mTOR in CTCs"""
import sys, os
sys.path.insert(0, r'D:\科研\CTC知识图谱\v2')
sys.stdout.reconfigure(encoding='utf-8')
from neo4j import GraphDatabase
driver = GraphDatabase.driver('bolt://127.0.0.1:7687', auth=(os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"]))

BATCH3 = {
    "gene_phenotype": [
        {"gene":"TWIST1","phenotype":"metastasis","direction":"DRIVES",
         "evidence":"TWIST1+ TAC1+ tumor subpopulation exhibited enhanced secretion of substance P driving metastasis",
         "confidence":"high","pmid":"41985871","cancer":"Neuroblastoma"},
        {"gene":"TAC1","phenotype":"metastasis","direction":"DRIVES",
         "evidence":"TAC1 overexpression was sufficient to accelerate tumorigenesis, CTC generation and metastatic progression",
         "confidence":"high","pmid":"41985871","cancer":"Neuroblastoma"},
        {"gene":"TAC1","phenotype":"angiogenesis","direction":"ASSOCIATED_WITH",
         "evidence":"TAC1 signaled to endothelial cells with TACR1 expression, leading to endothelial senescence and disrupted vascular integrity",
         "confidence":"medium","pmid":"41985871","cancer":"Neuroblastoma"},
        {"gene":"IL12A","phenotype":"immune_evasion","direction":"SUPPRESSES",
         "evidence":"EMT CTC showed upregulation of IL12A ... linked to the patient's beneficial prognosis",
         "confidence":"medium","pmid":"41367869","cancer":"Ovarian Cancer"},
        {"gene":"VEGFA","phenotype":"metastasis","direction":"DRIVES",
         "evidence":"epithelial CTC expressed high levels of VEGFA",
         "confidence":"medium","pmid":"41367869","cancer":"Ovarian Cancer"},
    ],
    "gene_gene": [
        {"gene_a":"TWIST1","gene_b":"TAC1","relation":"UPREGULATES",
         "mechanism":"TWIST1 transcription factor regulates TAC1 (substance P) expression",
         "evidence":"TAC1 expression was regulated by the transcription factor TWIST1",
         "confidence":"high","pmid":"41985871"},
        {"gene_a":"TAC1","gene_b":"TACR1","relation":"ACTIVATES",
         "mechanism":"Substance P (TAC1) paracrine signaling to TACR1+ endothelial cells",
         "evidence":"TAC1 signaled to tumor-associated endothelial cells with elevated TACR1 expression",
         "confidence":"high","pmid":"41985871"},
    ],
    "gene_pathway": [
        {"gene":"MTOR","pathway":"PI3K_AKT_MTOR_SIGNALING","relation":"ACTIVATES",
         "evidence":"EMT CTC showed upregulation of mTOR pathway (RPTOR, RICTOR, MTOR)","confidence":"medium","pmid":"41367869"},
        {"gene":"AKT1","pathway":"PI3K_AKT_MTOR_SIGNALING","relation":"ACTIVATES",
         "evidence":"AKT1 upregulated, putatively indicating mTOR pathway activation","confidence":"medium","pmid":"41367869"},
    ],
    "biomarkers": [
        {"gene":"FOLR1","cancer_type":"Ovarian Cancer","type":"diagnostic",
         "evidence":"FRalpha as HGSOC-specific CTC marker; significant reduction post-treatment","pmid":"41437733"},
        {"gene":"TAC1","cancer_type":"Neuroblastoma","type":"prognostic",
         "evidence":"TAC1 overexpression accelerated CTC generation and metastatic progression","pmid":"41985871"},
    ],
    "ctc_subpopulations": [
        {"gene":"IL12A","subpop_type":"Hybrid E/M CTC","evidence":"EMT CTC showed upregulation of IL12A, mTOR pathway genes","pmid":"41367869"},
        {"gene":"VEGFA","subpop_type":"Epithelial CTC","evidence":"epithelial CTC expressed high levels of VEGFA","pmid":"41367869"},
        {"gene":"TAC1","subpop_type":"Hybrid E/M CTC","evidence":"C2 subtype marked by MYCN amplification and dual TWIST1/TAC1 expression","pmid":"41985871"},
    ],
}

def check(gene):
    with driver.session() as s:
        n = s.run("MATCH (g:Gene {symbol:$s}) RETURN g", s=gene.upper()).single()
        return n is not None

def run():
    print("Batch 3: TAC1-TACR1 + Ovarian CTC")
    stats = {"gp":0,"gg":0,"gpw":0,"bio":0,"ctc":0,"skip":0}
    with driver.session() as s:
        for r in BATCH3["gene_phenotype"]:
            g = r["gene"].upper()
            if not check(g): print(f"  SKIP {g}"); stats["skip"]+=1; continue
            d = r["phenotype"].replace("_"," ").title()
            s.run("MERGE (ph:CTC_Phenotype {name:$n})", n=d)
            cm = {"high":0.9,"medium":0.6,"low":0.3}
            s.run("""MATCH (g:Gene {symbol:$g}) MATCH (ph:CTC_Phenotype {name:$ph})
                MERGE (g)-[rel:REGULATES {pmid:$pm}]->(ph)
                SET rel.direction=$d, rel.mechanism=$e, rel.confidence=$c,
                    rel.evidence_type='hand_curated', rel.cancer_context=$ca, rel.pmid=$pm
            """, g=g, ph=d, d=r["direction"], e=r["evidence"][:300],
                 c=cm.get(r["confidence"],0.5), pm=r["pmid"], ca=r.get("cancer",""))
            stats["gp"]+=1; print(f"  {g} -[{r['direction']}]-> {d} ({r['pmid']})")

        for r in BATCH3["gene_gene"]:
            ga,gb = r["gene_a"].upper(), r["gene_b"].upper()
            if not check(ga) or not check(gb): print(f"  SKIP {ga}->{gb}"); stats["skip"]+=1; continue
            cm = {"high":0.9,"medium":0.6,"low":0.3}
            s.run("""MATCH (a:Gene {symbol:$a}) MATCH (b:Gene {symbol:$b})
                MERGE (a)-[rel:REGULATES {pmid:$pm}]->(b)
                SET rel.relation_type=$rt, rel.mechanism=$m, rel.confidence=$c,
                    rel.evidence_type='hand_curated', rel.pmid=$pm
            """, a=ga,b=gb,rt=r["relation"],m=r.get("mechanism","")[:300],
                 c=cm.get(r["confidence"],0.5),pm=r["pmid"])
            stats["gg"]+=1; print(f"  {ga} -[{r['relation']}]-> {gb} ({r['pmid']})")

        for r in BATCH3["gene_pathway"]:
            g = r["gene"].upper()
            if not check(g): print(f"  SKIP {g}"); continue
            pw = s.run("""MATCH (p:Pathway) WHERE toLower(p.name) CONTAINS toLower($pw) RETURN p.name AS n LIMIT 1""", pw=r["pathway"]).single()
            if pw:
                cm = {"high":0.9,"medium":0.6,"low":0.3}
                s.run("""MATCH (g:Gene {symbol:$g}) MATCH (p:Pathway {name:$pn})
                    MERGE (g)-[rel:REGULATES {pmid:$pm}]->(p)
                    SET rel.relation_type=$rt, rel.mechanism=$e, rel.confidence=$c,
                        rel.evidence_type='hand_curated', rel.pmid=$pm
                """, g=g, pn=pw["n"], rt=r["relation"], e=r["evidence"][:300],
                     c=cm.get(r["confidence"],0.5), pm=r["pmid"])
                stats["gpw"]+=1; print(f"  {g} -[{r['relation']}]-> {pw['n']} ({r['pmid']})")

        for r in BATCH3["biomarkers"]:
            g = r["gene"].upper()
            if not check(g): print(f"  SKIP {g}"); continue
            cn = s.run("""MATCH (c:CancerType) WHERE toLower(c.name) CONTAINS toLower($cn) RETURN c.name AS n LIMIT 1""", cn=r["cancer_type"]).single()
            if cn:
                s.run("""MATCH (g:Gene {symbol:$g}) MATCH (c:CancerType {name:$cn})
                    MERGE (g)-[rel:BIOMARKER_FOR {pmid:$pm}]->(c)
                    SET rel.biomarker_type=$bt, rel.evidence=$e, rel.evidence_type='hand_curated'
                """, g=g, cn=cn["n"], bt=r["type"], e=r["evidence"][:300], pm=r["pmid"])
                stats["bio"]+=1; print(f"  {g} -[BIOMARKER_FOR {r['type']}]-> {cn['n']} ({r['pmid']})")

        for r in BATCH3["ctc_subpopulations"]:
            g = r["gene"].upper()
            if not check(g): print(f"  SKIP {g}"); continue
            s.run("MERGE (csp:CTC_Subpopulation {name:$n})", n=r["subpop_type"])
            s.run("""MATCH (g:Gene {symbol:$g}) MATCH (csp:CTC_Subpopulation {name:$sp})
                MERGE (g)-[rel:EPITHELIAL_MESENCHYMAL_STATE {pmid:$pm}]->(csp)
                SET rel.evidence=$e, rel.evidence_type='hand_curated'
            """, g=g, sp=r["subpop_type"], e=r["evidence"][:300], pm=r["pmid"])
            stats["ctc"]+=1; print(f"  {g} -> {r['subpop_type']} ({r['pmid']})")

    print(f"\nBatch3: gp={stats['gp']} gg={stats['gg']} gpw={stats['gpw']} bio={stats['bio']} ctc={stats['ctc']} skip={stats['skip']}")
    with driver.session() as s:
        for e in ["REGULATES","BIOMARKER_FOR","EPITHELIAL_MESENCHYMAL_STATE"]:
            c = s.run(f"MATCH ()-[r:{e}]->() RETURN count(r) AS n").single()
            print(f"  Total {e}: {c['n']}")

if __name__ == "__main__":
    run()
    driver.close()
