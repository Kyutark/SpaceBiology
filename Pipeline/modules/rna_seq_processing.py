from pathlib import Path
import os
import pandas as pd

def process_rna_seq():
    # 현재 스크립트 파일의 디렉토리를 기준으로 경로 설정
    script_dir = Path(__file__).resolve().parent
    mdsh_csv = script_dir.parent / "data" / "mdsh.csv"

    if not mdsh_csv.exists():
        raise FileNotFoundError(f"File {mdsh_csv} does not exist!")

    # CSV 파일 로드
    mdsh_df = pd.read_csv(mdsh_csv)
    runs = mdsh_df['run'].dropna().tolist()

    # 결과 디렉토리 생성
    result_dir = script_dir.parent / "results"
    os.makedirs(result_dir, exist_ok=True)

    # 각 run에 대해 RNA-seq 처리 실행
    for run in runs:
        ref = mdsh_df.loc[mdsh_df['run'] == run, 'ref'].values[0]
        fq1 = script_dir.parent / "data" / f"{run}_1.fastq"
        fq2 = script_dir.parent / "data" / f"{run}_2.fastq" if (script_dir.parent / "data" / f"{run}_2.fastq").exists() else None
        fna_file = script_dir.parent / "data" / f"{ref}.fna"
        gff_file = script_dir.parent / "data" / f"{ref}.gff"

        print(f"Processing RNA-seq for run: {run} with reference: {ref}")

        # 품질 제어
        if fq2:
            os.system(f"fastp -i {fq1} -I {fq2} -o {result_dir}/{run}_filtered_1.fastq -O {result_dir}/{run}_filtered_2.fastq")
        else:
            os.system(f"fastp -i {fq1} -o {result_dir}/{run}_filtered.fastq")

        # Bowtie2 인덱스 생성 및 정렬
        os.system(f"bowtie2-build {fna_file} {result_dir}/{run}_index")
        if fq2:
            os.system(f"bowtie2 -x {result_dir}/{run}_index -1 {result_dir}/{run}_filtered_1.fastq -2 {result_dir}/{run}_filtered_2.fastq -S {result_dir}/{run}.sam")
        else:
            os.system(f"bowtie2 -x {result_dir}/{run}_index -U {result_dir}/{run}_filtered.fastq -S {result_dir}/{run}.sam")

        # FeatureCounts를 통한 카운팅 - gff 파일에서 'ID' 속성 사용
        if fq2:
            # Paired-end인 경우 -p 옵션 추가
            os.system(f"featureCounts -a {gff_file} -o {result_dir}/{run}_counts.txt -t CDS -g ID -p {result_dir}/{run}.sam")
        else:
            # Single-end인 경우 -p 옵션 없이 실행
            os.system(f"featureCounts -a {gff_file} -o {result_dir}/{run}_counts.txt -t CDS -g ID {result_dir}/{run}.sam")

if __name__ == "__main__":
    process_rna_seq()

