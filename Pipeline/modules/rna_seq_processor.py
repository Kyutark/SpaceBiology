from pathlib import Path
import os
import pandas as pd

def process_rna_seq():
    # 현재 스크립트 파일의 디렉토리를 기준으로 경로 설정
    script_dir = Path(__file__).resolve().parent  # /modules
    data_dir = script_dir.parent / "data"  # /data
    result_dir = script_dir.parent / "results"  # /results

    # CSV 파일 로드
    mdsh_csv = data_dir / "mdsh.csv"
    if not mdsh_csv.exists():
        raise FileNotFoundError(f"File {mdsh_csv} does not exist!")

    mdsh_df = pd.read_csv(mdsh_csv)
    runs = mdsh_df['run'].dropna().tolist()

    # 결과 디렉토리 생성
    os.makedirs(result_dir, exist_ok=True)

    # 각 run에 대해 RNA-seq 처리 실행
    for run in runs:
        # 관련 파일 경로 설정
        ref = mdsh_df.loc[mdsh_df['run'] == run, 'ref'].values[0]
        fq1 = data_dir / f"{run}_1.fastq.gz"
        fq2 = data_dir / f"{run}_2.fastq.gz" if (data_dir / f"{run}_2.fastq.gz").exists() else None
        fna_file = data_dir / f"{ref}.fna"
        gff_file = data_dir / f"{ref}.gff"

        if not fq1.exists():
            print(f"Error: FASTQ file {fq1} does not exist for run {run}")
            continue

        print(f"Processing RNA-seq for run: {run} with reference: {ref}")

        # 품질 제어 (fastp로 압축된 파일 처리)
        if fq2:  # Paired-end
            os.system(f"fastp -i {fq1} -I {fq2} -o {result_dir}/{run}_filtered_1.fastq -O {result_dir}/{run}_filtered_2.fastq")
            # pigz로 병렬 압축
            os.system(f"pigz -p 24 {result_dir}/{run}_filtered_1.fastq")
            os.system(f"pigz -p 24 {result_dir}/{run}_filtered_2.fastq")
        else:  # Single-end
            os.system(f"fastp -i {fq1} -o {result_dir}/{run}_filtered.fastq")
            # pigz로 병렬 압축
            os.system(f"pigz -p 24 {result_dir}/{run}_filtered.fastq")

        # Bowtie2 인덱스 생성
        os.system(f"bowtie2-build {fna_file} {result_dir}/{run}_index")

        # Bowtie2 정렬 후 바로 BAM 파일 생성 및 정렬
        sorted_bam_file = result_dir / f"{run}_sorted.bam"
        if fq2:  # Paired-end
            os.system(f"bowtie2 -p 24 -x {result_dir}/{run}_index -1 {result_dir}/{run}_filtered_1.fastq.gz -2 {result_dir}/{run}_filtered_2.fastq.gz | samtools view -bS | samtools sort -o {sorted_bam_file}")
        else:  # Single-end
            os.system(f"bowtie2 -p 24 -x {result_dir}/{run}_index -U {result_dir}/{run}_filtered.fastq.gz | samtools view -bS | samtools sort -o {sorted_bam_file}")

        # BAM 파일 인덱싱 (선택 사항)
        os.system(f"samtools index {sorted_bam_file}")

        # FeatureCounts를 통한 카운팅 - gff 파일에서 'ID' 속성 사용
        if fq2:  # Paired-end
            os.system(f"featureCounts -a {gff_file} -o {result_dir}/{run}_counts.txt -t CDS -g ID -p {sorted_bam_file}")
        else:  # Single-end
            os.system(f"featureCounts -a {gff_file} -o {result_dir}/{run}_counts.txt -t CDS -g ID {sorted_bam_file}")

if __name__ == "__main__":
    process_rna_seq()

