import os
import subprocess
import shutil

# 경로 설정
result_dir = "../results"
orthofinder_input_dir = "../results/orthofinder_input"
orthofinder_output_dir = "../results/orthofinder_output"

# 기존 디렉토리 삭제
if os.path.exists(orthofinder_output_dir):
    shutil.rmtree(orthofinder_output_dir)

# 입력 디렉토리 생성 (기존 디렉토리 삭제)
if os.path.exists(orthofinder_input_dir):
    shutil.rmtree(orthofinder_input_dir)
os.makedirs(orthofinder_input_dir, exist_ok=True)

# _cds.fasta 파일 수집 및 복사
fasta_files = [f for f in os.listdir(result_dir) if f.endswith("_cds.fasta")]
for fasta_file in fasta_files:
    full_path = os.path.join(result_dir, fasta_file)
    shutil.copy(full_path, orthofinder_input_dir)

# OrthoFinder 실행 (MMseqs 사용)
try:
    print("Running OrthoFinder with MMseqs...")
    subprocess.run(
        [
            "orthofinder",
            "-S", "mmseqs",  # MMseqs 사용 설정
            "-f", orthofinder_input_dir,
            "-o", orthofinder_output_dir
        ],
        check=True
    )
    print(f"OrthoFinder completed with MMseqs! Results are stored in {orthofinder_output_dir}.")
except subprocess.CalledProcessError as e:
    print(f"An error occurred during OrthoFinder execution with MMseqs: {e}")

