import subprocess
import pandas as pd

def run_command(command):
    """ 주어진 커맨드를 실행하고, 출력을 반환합니다. """
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        print("Error executing command: ", command)
        print(stderr.decode())
        return None
    return stdout.decode()

# 데이터 로드
data = pd.read_excel('dataset.xlsx', index_col=0)

# 각 SRR ID에 대해 처리 실행
for index, row in data.iterrows():
    srr_id = row['SRR']
    species = row['species']
    condition = row['condition']

    # SRR ID를 FASTQ로 변환
    run_command(f"fastq-dump --split-files {srr_id}")

    # Trim Galore 실행
    run_command(f"trim_galore --illumina {srr_id}_1.fastq")

    # BWA-MEM2 실행 (참조 게놈 필요)
    run_command(f"bwa mem -t 8 ref_{species}.fa {srr_id}_1_trimmed.fq > {srr_id}.sam")

    # FeatureCounts 실행
    run_command(f"featureCounts -a ref_{species}.gtf -o counts_{srr_id}.txt {srr_id}.sam")

    # EdgeR 실행 (R 스크립트 필요)
    run_command(f"Rscript edgeR_analysis.R counts_{srr_id}.txt results_{srr_id}.csv")

print("DEG 분석이 완료되었습니다.")
