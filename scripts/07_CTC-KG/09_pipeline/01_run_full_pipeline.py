"""
一键运行完整深度文献→KG集成管道
===============================
Step 1: 对全部论文做增强LLM提取
Step 2: KG交叉验证
Step 3: 导入Neo4j
Step 4: 生成证据链报告
Step 5: 多模态层构建

用法:
  python 09_pipeline/01_run_full_pipeline.py [--resume]
"""

import subprocess, sys, os, json

SCRIPTS_DIR = r"D:\科研\CTC知识图谱\v2\scripts"
DATA_DIR = r"D:\科研\CTC知识图谱\v2\data\processed"


def run_step(name, script, args=None):
    """Run a pipeline step"""
    print(f"\n{'='*60}")
    print(f"STEP: {name}")
    print(f"{'='*60}")
    cmd = [sys.executable, os.path.join(SCRIPTS_DIR, script)]
    if args:
        cmd.extend(args)
    result = subprocess.run(cmd, capture_output=False, text=True)
    if result.returncode != 0:
        print(f"WARNING: {name} exited with code {result.returncode}")
    return result.returncode


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--resume', action='store_true', help='Resume from saved state')
    parser.add_argument('--skip-extract', action='store_true', help='Skip extraction (use existing)')
    args = parser.parse_args()

    print("=" * 60)
    print("CTC-KG v3: Full Deep Literature → KG Pipeline")
    print("=" * 60)

    # Load abstracts
    abstracts_file = os.path.join(DATA_DIR, "abstracts.json")
    with open(abstracts_file, 'r', encoding='utf-8') as f:
        abstracts = json.load(f)
    print(f"Total papers in library: {len(abstracts)}")

    # Step 1: Extract
    if not args.skip_extract:
        run_step("Extraction", r"04_literature_mining\06_deep_literature_kg_integration.py", ["--mode", "extract"])
    else:
        print("Skipping extraction (--skip-extract)")

    # Step 2+3: Cross-reference + Import
    # These are combined in --mode full's import step
    run_step("Cross-reference + Import", r"04_literature_mining\06_deep_literature_kg_integration.py", ["--mode", "import"])

    # Step 4: Report
    run_step("Evidence Chain Report", r"04_literature_mining\06_deep_literature_kg_integration.py", ["--mode", "report"])

    # Step 5: Multi-modal layers
    run_step("Multi-Modal Layers", r"06_multimodal_layers\01_build_multimodal_layers.py")

    print("\n" + "=" * 60)
    print("PIPELINE COMPLETE")
    print("=" * 60)


if __name__ == '__main__':
    main()
