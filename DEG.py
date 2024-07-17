import subprocess
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri

# R 패키지 로드
edgeR = importr('edgeR')
pandas2ri.activate()

def run_command(command):
    """ 주어진 커맨드를 실행하고, 출력을 반환합니다. """
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        print("Error executing command: ", command)
        print(stderr.decode())
        return None
    return stdout.decode()

def run_edgeR_analysis(counts_file, output_file):
    """ edgeR를 사용한 DEG 분석을 수행하는 함수 """
    ro.r('''
    library(edgeR)
    function(counts_file, output_file) {
      counts <- read.delim(counts_file, row.names = 1)
      group <- factor(c(rep("NG", times=3), rep("MG", times=3)))
      dge <- DGEList(counts = counts, group = group)
      dge <- calcNormFactors(dge)
      design <- model.matrix(~ group)
      dge <- estimateDisp(dge, design)
      fit <- glmQLFit(dge, design)
      results <- glmQLFTest(fit)
      write.csv(as.data.frame(results$table), file = output_file)
    }
    ''')(counts_file, output_file)

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

    # EdgeR 분석
    run_edgeR_analysis(f"counts_{srr_id}.txt", f"results_{srr_id}.csv")

print("DEG 분석이 완료되었습니다.")
