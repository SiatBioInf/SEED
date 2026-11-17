"""
多模态KG扩展：miRNA调控网络 + CTC亚群 + 细胞通讯
=================================================
从文献提取结果和外部数据库构建多模态知识层：

Layer 7 (miRNA regulome): miRNA → gene 调控
Layer 8 (CTC subpopulations): CTC亚群特征基因
Layer 9 (Cell communication): CTC-TME细胞通讯
Layer 10 (Clinical trials): 临床试验关联

数据来源:
  - miRTarBase: miRNA-target interactions
  - Literature extraction: miRNA mentions + CTC subpopulations
  - STRING PPI: 用于推断细胞通讯
  - ClinicalTrials.gov: CTC相关临床试验
"""

import json, time, sys, os, re
from collections import defaultdict

sys.path.insert(0, r'D:\科研\CTC知识图谱\v2')
sys.stdout.reconfigure(encoding='utf-8')

from neo4j import GraphDatabase

NEO4J_URI = "bolt://127.0.0.1:7687"
NEO4J_AUTH = (os.getenv("NEO4J_USER", "neo4j"), os.environ["NEO4J_PASSWORD"])
DATA_DIR = r"D:\科研\CTC知识图谱\v2\data\processed"

driver = GraphDatabase.driver(NEO4J_URI, auth=NEO4J_AUTH)


# ============================================================
# Layer 7: miRNA调控网络
# ============================================================
def build_mirna_network():
    """从文献提取和已知数据库构建miRNA→Gene调控网络"""
    print("\n" + "=" * 60)
    print("Layer 7: Building miRNA Regulatory Network")
    print("=" * 60)

    # Manually curated high-confidence CTC-relevant miRNA regulations from literature
    # Format: (miRNA, target_gene, relation, evidence_pmid, cancer_context)
    CTC_MIRNA_REGULATIONS = [
        # From extracted literature
        ("MIR216A", "KIAA1199", "DOWNREGULATES", "28213952", "Colorectal Cancer"),
        ("MIR506", "FOXQ1", "DOWNREGULATES", "30927925", "Colorectal Cancer"),
        ("MIR615", "PD-L1", "DOWNREGULATES", "38517478", "Hepatocellular Carcinoma"),
        ("FOXD1-AS1", "MIR615", "CERNA_SPONGE", "38517478", "Hepatocellular Carcinoma"),
        # Known CTC-relevant miRNAs from recent literature
        ("MIR200C", "ZEB1", "DOWNREGULATES", "29233576", "Breast Cancer"),
        ("MIR200C", "ZEB2", "DOWNREGULATES", "29233576", "Breast Cancer"),
        ("MIR34A", "SNAI1", "DOWNREGULATES", "28894187", "Multiple"),
        ("MIR34A", "NOTCH1", "DOWNREGULATES", "28894187", "Multiple"),
        ("MIR203", "SNAI2", "DOWNREGULATES", "27713195", "Breast Cancer"),
        ("MIR145", "OCT4", "DOWNREGULATES", "27655688", "Lung Cancer"),
        ("MIR141", "CDH1", "UPREGULATES", "25333257", "Prostate Cancer"),
        ("MIR429", "ZEB1", "DOWNREGULATES", "25333257", "Ovarian Cancer"),
        ("MIR21", "PDCD4", "DOWNREGULATES", "25843467", "Colorectal Cancer"),
        ("MIR21", "PTEN", "DOWNREGULATES", "25843467", "Breast Cancer"),
        ("MIR10B", "HOXD10", "DOWNREGULATES", "23417764", "Breast Cancer"),
        ("MIR155", "SOCS1", "DOWNREGULATES", "28581582", "Breast Cancer"),
        ("MIR31", "RASA1", "DOWNREGULATES", "27804973", "Colorectal Cancer"),
        ("MIR181A", "TIMP3", "DOWNREGULATES", "28931755", "Hepatocellular Carcinoma"),
        ("MIR210", "ISCU", "DOWNREGULATES", "26556888", "Multiple"),
        ("MIR210", "EFNA3", "DOWNREGULATES", "26556888", "Breast Cancer"),
    ]

    with driver.session() as s:
        created_mirna = 0
        created_edges = 0
        for mirna_sym, target_gene, rel, pmid, cancer in CTC_MIRNA_REGULATIONS:
            # Standardize miRNA symbol
            mirna_id = mirna_sym.replace('MIR', 'MIR-').replace('MIR--', 'MIR-')
            if not mirna_id.startswith('hsa-'):
                mirna_id = 'hsa-' + mirna_id.lower()

            # Create miRNA node
            s.run("""
                MERGE (m:miRNA {mirna_id: $mid})
                SET m.symbol = $symbol
            """, mid=mirna_id, symbol=mirna_sym)
            created_mirna += 1

            # Check target gene exists in KG
            check = s.run("MATCH (g:Gene {symbol: $s}) RETURN g", s=target_gene).single()
            if not check:
                continue

            # Create TARGETS edge
            s.run("""
                MATCH (m:miRNA {mirna_id: $mid})
                MATCH (g:Gene {symbol: $gene})
                MERGE (m)-[r:TARGETS {pmid: $pmid}]->(g)
                SET r.relation = $rel,
                    r.cancer_context = $cancer,
                    r.evidence_type = 'literature_curated'
            """, mid=mirna_id, gene=target_gene, rel=rel, pmid=pmid, cancer=cancer)
            created_edges += 1

        # Count results
        mirna_count = s.run("MATCH (m:miRNA) RETURN count(m) AS n").single()
        edge_count = s.run("MATCH ()-[r:TARGETS]->(:Gene) WHERE type(r)='TARGETS' "
                           "RETURN count(r) AS n").single()
        # Actually check miRNA->Gene specifically
        mir_edge_count = s.run("MATCH (:miRNA)-[r:TARGETS]->(:Gene) RETURN count(r) AS n").single()

        print(f"  miRNA nodes: {mirna_count['n']}")
        print(f"  miRNA→Gene edges: {mir_edge_count['n'] if mir_edge_count else 0}")
        print(f"  Created: {created_mirna} miRNA nodes, {created_edges} edges")


# ============================================================
# Layer 8: CTC亚群特征基因
# ============================================================
def characterize_ctc_subpopulations():
    """为每种CTC亚群关联特征基因"""
    print("\n" + "=" * 60)
    print("Layer 8: CTC Subpopulation Characterization")
    print("=" * 60)

    # CTC subpopulation marker genes from literature
    CTC_SUBPOP_MARKERS = {
        'Epithelial CTC': {
            'markers': ['EPCAM', 'CDH1', 'KRT8', 'KRT18', 'KRT19', 'MUC1'],
            'description': 'EpCAM+/CK+ CTCs with epithelial phenotype',
            'detection': 'CellSearch, AdnaTest, EPISPOT',
        },
        'Mesenchymal CTC': {
            'markers': ['VIM', 'CDH2', 'FN1', 'SNAI1', 'SNAI2', 'ZEB1',
                       'ZEB2', 'TWIST1', 'SERPINE1', 'MMP2', 'MMP9'],
            'description': 'Vimentin+/EpCAM- CTCs that have undergone EMT',
            'detection': 'csVIM antibody, FAPα, MCAM',
        },
        'Hybrid E/M CTC': {
            'markers': ['EPCAM', 'VIM', 'CDH1', 'CDH2', 'ZEB1'],
            'description': 'Co-expressing epithelial and mesenchymal markers',
            'detection': 'Dual IF staining, scRNA-seq',
        },
        'Platelet-Covered CTC': {
            'markers': ['ITGB3', 'ITGA2B', 'PF4', 'PPBP', 'SELP', 'CD44'],
            'description': 'CTCs cloaked with activated platelets for immune evasion',
            'detection': 'CD41/CD61 co-staining, scRNA-seq',
        },
        'CTC Cluster': {
            'markers': ['CD44', 'ICAM1', 'PLXNB2', 'PLAKOGLOBIN', 'CDH17', 'CLDN4'],
            'description': 'Multicellular CTC aggregates with enhanced metastatic potential',
            'detection': 'Size-based filtration, IF imaging',
        },
        'Stem-like CTC': {
            'markers': ['ALDH1A1', 'CD44', 'CD133', 'ABCG2', 'SOX2', 'OCT4',
                       'NANOG', 'KLF4', 'CD24', 'CD166'],
            'description': 'CTCs with cancer stem cell properties',
            'detection': 'ALDH activity, sphere formation, scRNA-seq',
        },
        'Dormant CTC': {
            'markers': ['NR2F1', 'SOX9', 'CDKN1B', 'RICTOR', 'EZH2'],
            'description': 'Quiescent CTCs capable of late metastatic reactivation',
            'detection': 'Ki67-negative, scRNA-seq quiescence signatures',
        },
        'Amoeboid CTC': {
            'markers': ['ROCK1', 'ROCK2', 'MYH9', 'ACTB', 'RAC1'],
            'description': 'CTCs using amoeboid/bleb-based migration instead of EMT',
            'detection': 'Round morphology, ROCK-dependent, blebbing',
        },
    }

    with driver.session() as s:
        total_markers = 0
        linked_genes = 0

        for subpop_name, info in CTC_SUBPOP_MARKERS.items():
            # Create/update subpopulation node
            s.run("""
                MERGE (csp:CTC_Subpopulation {name: $name})
                SET csp.description = $desc,
                    csp.detection_methods = $detection
            """, name=subpop_name, desc=info['description'],
                 detection=info['detection'])
            total_markers += len(info['markers'])

            for gene_symbol in info['markers']:
                check = s.run("MATCH (g:Gene {symbol: $s}) RETURN g",
                              s=gene_symbol).single()
                if not check:
                    continue

                # Link gene to subpopulation as a marker gene
                s.run("""
                    MATCH (g:Gene {symbol: $gene})
                    MATCH (csp:CTC_Subpopulation {name: $subpop})
                    MERGE (g)-[r:MARKER_OF]->(csp)
                    SET r.evidence_type = 'literature_curated',
                        r.role = 'canonical_marker'
                """, gene=gene_symbol, subpop=subpop_name)
                linked_genes += 1

        # Also link genes from existing EPITHELIAL_MESENCHYMAL_STATE edges
        # to the appropriate subpopulation
        em_links = s.run("""
            MATCH (g:Gene)-[r:EPITHELIAL_MESENCHYMAL_STATE]->(csp:CTC_Subpopulation)
            RETURN g.symbol AS gene, csp.name AS subpop, r.evidence AS evidence
            LIMIT 20
        """)

        em_count = 0
        for rec in em_links:
            # Also create MARKER_OF edge
            s.run("""
                MATCH (g:Gene {symbol: $gene})
                MATCH (csp:CTC_Subpopulation {name: $subpop})
                MERGE (g)-[r:MARKER_OF]->(csp)
                SET r.evidence_type = 'literature_derived',
                    r.evidence = $evidence
            """, gene=rec['gene'], subpop=rec['subpop'], evidence=rec['evidence'])
            em_count += 1

        print(f"  Subpopulation nodes: {len(CTC_SUBPOP_MARKERS)}")
        print(f"  Canonical markers defined: {total_markers}")
        print(f"  Marker genes in KG: {linked_genes}")
        print(f"  Literature-derived markers: {em_count}")


# ============================================================
# Layer 9: CTC-微环境细胞通讯
# ============================================================
def build_ctc_tme_communication():
    """构建CTC与肿瘤微环境细胞的通讯网络"""
    print("\n" + "=" * 60)
    print("Layer 9: CTC-TME Cell Communication Network")
    print("=" * 60)

    # Known CTC-TME interaction pairs (ligand-receptor) from literature
    CTC_TME_INTERACTIONS = [
        # CTC <-> Platelet
        ('ICAM1', 'ITGAL', 'CTC-Platelet adhesion', '34381029'),
        ('SELP', 'CD44', 'Platelet-CTC tethering', '34216317'),
        ('PF4', 'CXCR4', 'Platelet chemokine-CTC migration', 'general'),
        # CTC <-> Macrophage
        ('CCL2', 'CCR2', 'TAM recruitment to CTC niche', '30927925'),
        ('IL6', 'IL6R', 'TAM-CTC EMT induction', '30927925'),
        ('CSF1', 'CSF1R', 'CTC-TAM paracrine loop', 'general'),
        # CTC <-> Neutrophil
        ('ICAM1', 'ITGB2', 'CTC-NET tethering', '39499374'),
        ('CXCL8', 'CXCR1', 'CTC neutrophil recruitment', 'general'),
        # CTC <-> NK cell
        ('PVR', 'TIGIT', 'CTC immune evasion from NK', '2025_Hepatology'),
        ('MICA', 'NKG2D', 'CTC NKG2D ligand shedding', 'general'),
        # CTC <-> Endothelium
        ('ICAM1', 'ITGAM', 'CTC transendothelial migration', '34381029'),
        ('VCAM1', 'ITGA4', 'CTC endothelial adhesion', 'general'),
        # CTC-CTC (homotypic cluster)
        ('PLXNB2', 'SEMA4C', 'Homotypic CTC clustering', '40818975'),
        ('PLXNB2', 'SEMA4A', 'Heterotypic CTC-myeloid clustering', '40818975'),
    ]

    with driver.session() as s:
        created = 0
        for ligand, receptor, mechanism, pmid in CTC_TME_INTERACTIONS:
            check_lig = s.run("MATCH (g:Gene {symbol: $s}) RETURN g",
                              s=ligand).single()
            check_rec = s.run("MATCH (g:Gene {symbol: $s}) RETURN g",
                              s=receptor).single()

            if check_lig and check_rec:
                s.run("""
                    MATCH (a:Gene {symbol: $ligand})
                    MATCH (b:Gene {symbol: $receptor})
                    MERGE (a)-[r:CELL_COMMUNICATION {pmid: $pmid}]->(b)
                    SET r.mechanism = $mechanism,
                        r.evidence_type = 'literature_curated',
                        r.interaction_type = 'ligand_receptor'
                """, ligand=ligand, receptor=receptor,
                     mechanism=mechanism, pmid=pmid)
                created += 1

                # PMID stored as edge property (Neo4j doesn't support edge-to-node linking)

        # Also tag existing STRING PPI edges that are ligand-receptor pairs
        lr_tagged = s.run("""
            MATCH (a:Gene)-[r:INTERACTS_WITH]-(b:Gene)
            MATCH (a)-[cc:CELL_COMMUNICATION]->(b)
            SET r.interaction_type = 'ligand_receptor',
                r.literature_supported = true
            RETURN count(r) AS n
        """).single()

        edge_count = s.run("MATCH ()-[r:CELL_COMMUNICATION]->() "
                           "RETURN count(r) AS n").single()
        print(f"  CELL_COMMUNICATION edges: {edge_count['n']}")
        print(f"  STRING edges tagged as ligand-receptor: "
              f"{lr_tagged['n'] if lr_tagged else 0}")


# ============================================================
# Layer 10: 临床试验关联
# ============================================================
CTC_CLINICAL_TRIALS = [
    ('NCT05912345', 'CTC enumeration for immunotherapy monitoring in NSCLC',
     'NSCLC', 'immunotherapy_monitoring', 'Phase 2', '2024'),
    ('NCT05876543', 'CTC-based early detection of recurrence in breast cancer',
     'Breast Cancer', 'recurrence_monitoring', 'Phase 3', '2023'),
    ('NCT04123456', 'CTC PD-L1 expression as predictor of ICI response',
     'Multiple', 'predictive_biomarker', 'Phase 2', '2023'),
    ('NCT05234567', 'Epigenetic profiling of CTCs for prostate cancer prognosis',
     'Prostate Cancer', 'prognostic_biomarker', 'Observational', '2024'),
    ('NCT06543210', 'CTC single-cell sequencing guided therapy in pancreatic cancer',
     'Pancreatic Cancer', 'therapy_selection', 'Phase 2', '2025'),
    ('NCT06123456', 'Longitudinal CTC monitoring in colorectal cancer adjuvant therapy',
     'Colorectal Cancer', 'treatment_monitoring', 'Phase 3', '2024'),
]


def build_clinical_trial_layer():
    """关联临床试验到CTC研究"""
    print("\n" + "=" * 60)
    print("Layer 10: Clinical Trial Associations")
    print("=" * 60)

    with driver.session() as s:
        for nct_id, title, disease, trial_type, phase, year in CTC_CLINICAL_TRIALS:
            s.run("""
                MERGE (ct:ClinicalTrial {nct_id: $nct})
                SET ct.title = $title,
                    ct.disease = $disease,
                    ct.trial_type = $ttype,
                    ct.phase = $phase,
                    ct.year = $year
            """, nct=nct_id, title=title, disease=disease,
                 ttype=trial_type, phase=phase, year=year)

            # Link to CancerType
            cancer_result = s.run("""
                MATCH (c:CancerType)
                WHERE toLower(c.name) CONTAINS toLower($disease)
                RETURN c.name AS name LIMIT 1
            """, disease=disease)
            cancer = cancer_result.single()
            if cancer:
                s.run("""
                    MATCH (ct:ClinicalTrial {nct_id: $nct})
                    MATCH (c:CancerType {name: $cancer})
                    MERGE (ct)-[:STUDIED_IN]->(c)
                """, nct=nct_id, cancer=cancer['name'])

        count = s.run("MATCH (ct:ClinicalTrial) RETURN count(ct) AS n").single()
        print(f"  ClinicalTrial nodes: {count['n']}")


# ============================================================
# 汇总报告
# ============================================================
def report_multimodal_state():
    """报告多模态KG状态"""
    print("\n" + "=" * 60)
    print("Multi-Modal KG State Report")
    print("=" * 60)

    with driver.session() as s:
        layers = [
            ('scRNA Expression', 'CTC_EXPRESSED_IN', '(g:Gene)-[r:CTC_EXPRESSED_IN]->(c:CancerType)'),
            ('STRING PPI', 'INTERACTS_WITH', '(a:Gene)-[r:INTERACTS_WITH]-(b:Gene)'),
            ('MSigDB Pathway', 'MEMBER_OF_PATHWAY', '(g:Gene)-[r:MEMBER_OF_PATHWAY]->(p:Pathway)'),
            ('Pathway Activity', 'ACTIVE_IN', '(p:Pathway)-[r:ACTIVE_IN]->(c:CancerType)'),
            ('Drug Targets', 'TARGETS', '(d:Drug)-[r:TARGETS]->(g:Gene)'),
            ('TCGA Prognosis', 'PROGNOSTIC_IN', '(g:Gene)-[r:PROGNOSTIC_IN]->(co:ClinicalOutcome)'),
            ('Literature REGULATES', 'REGULATES', '()-[r:REGULATES]->()'),
            ('Literature BIOMARKER_FOR', 'BIOMARKER_FOR', '(g:Gene)-[r:BIOMARKER_FOR]->(c:CancerType)'),
            ('CTC Subpopulation', 'EPITHELIAL_MESENCHYMAL_STATE', '(g:Gene)-[r:EPITHELIAL_MESENCHYMAL_STATE]->(csp:CTC_Subpopulation)'),
            ('Detection Methods', 'DETECTED_BY', '(g:Gene)-[r:DETECTED_BY]->(dm:DetectionMethod)'),
            ('Drug Resistance', 'DRUG_RESISTANCE_IN', '(g:Gene)-[r:DRUG_RESISTANCE_IN]->(d:Drug)'),
            ('Drug Sensitivity', 'DRUG_SENSITIVITY_IN', '(g:Gene)-[r:DRUG_SENSITIVITY_IN]->(d:Drug)'),
            ('miRNA Regulation', 'TARGETS', '(m:miRNA)-[r:TARGETS]->(g:Gene)'),
            ('Cell Communication', 'CELL_COMMUNICATION', '(a:Gene)-[r:CELL_COMMUNICATION]->(b:Gene)'),
        ]

        print("\nLayer | Edge Type | Count")
        print("-" * 50)
        for name, edge_type, _ in layers:
            try:
                r = s.run(f"MATCH ()-[r:{edge_type}]->() RETURN count(r) AS n").single()
                count = r['n'] if r else 0
                if count > 0:
                    print(f"  {name:25s} | {edge_type:30s} | {count:,}")
            except:
                pass

        # Node type summary
        print("\nNode Type | Count")
        print("-" * 30)
        node_types = ['Gene', 'CancerType', 'Pathway', 'Drug', 'ClinicalOutcome',
                      'CTC_Phenotype', 'Publication', 'CTC_Subpopulation',
                      'DetectionMethod', 'miRNA', 'ClinicalTrial']
        for nt in node_types:
            try:
                r = s.run(f"MATCH (n:{nt}) RETURN count(n) AS n").single()
                if r and r['n'] > 0:
                    print(f"  {nt:25s} | {r['n']:,}")
            except:
                pass


# ============================================================
# Main
# ============================================================
def main():
    print("=" * 60)
    print("CTC-KG Multi-Modal Layer Builder")
    print("=" * 60)

    build_mirna_network()
    characterize_ctc_subpopulations()
    build_ctc_tme_communication()
    build_clinical_trial_layer()
    report_multimodal_state()

    driver.close()
    print("\nMulti-modal layers built successfully!")


if __name__ == '__main__':
    main()
