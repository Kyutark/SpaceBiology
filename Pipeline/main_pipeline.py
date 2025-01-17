
import subprocess
from pathlib import Path

# 현재 스크립트 디렉토리
script_dir = Path(__file__).resolve().parent

# 스크립트 경로 설정
scripts = [
    "modules/dataset_downloader.py",
    "modules/rna_seq_processor.py",
    "modules/cds_downloader.py",
    "modules/genecount_modifier.py",
    "modules/gene_clustering.py",
    "modules/deg_analyzer.py"
]

# 각 스크립트 순차적으로 실행
for script in scripts:
    try:
        print(f"Running {script}...")
        script_path = script_dir / script
        subprocess.run(["python3", script_path], check=True, cwd=script_path.parent)
        print(f"{script} completed successfully.\n")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running {script}: {e}")
        break  # 오류 발생 시, 파이프라인 중단
