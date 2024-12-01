from pathlib import Path
import os
import glob
import pandas as pd

def download_data(mdsh_csv=None):
    # /data 디렉토리 및 mdsh.csv 경로 설정
    script_dir = Path(__file__).resolve().parent  # /modules
    if mdsh_csv is None:
        mdsh_csv = script_dir.parent / "data" / "mdsh.csv"
    data_dir = script_dir.parent / "data"

    mdsh_df = pd.read_csv(mdsh_csv)
    unique_refs = mdsh_df['ref'].drop_duplicates().dropna().tolist()

    # GFF3와 FNA 파일 다운로드
    for ref in unique_refs:
        print(f"Downloading files for {ref}...")
        os.system(f"datasets download genome accession {ref} --include gff3,genome --filename {data_dir}/{ref}.zip")
        os.system(f"unzip -o {data_dir}/{ref}.zip -d {data_dir}")
        gff_files = glob.glob(f"{data_dir}/ncbi_dataset/data/{ref}/*.gff")
        fna_files = glob.glob(f"{data_dir}/ncbi_dataset/data/{ref}/*.fna")

        if gff_files:
            os.rename(gff_files[0], f"{data_dir}/{ref}.gff")
        else:
            print(f"Warning: GFF file not found for {ref}")

        if fna_files:
            os.rename(fna_files[0], f"{data_dir}/{ref}.fna")
        else:
            print(f"Warning: FNA file not found for {ref}")

    # FASTQ 파일 다운로드 및 pigz로 압축
    runs = mdsh_df['run'].dropna().tolist()
    for run in runs:
        if run.startswith('SRR'):  # SRA ID인 경우
            print(f"Downloading FASTQ files for {run}...")
            os.system(f"fastq-dump --outdir {data_dir} --split-files {run}")

            # 다운로드된 FASTQ 파일을 pigz로 압축
            fastq_files = glob.glob(f"{data_dir}/{run}_*.fastq")  # 다운로드된 파일 목록
            for fq in fastq_files:
                print(f"Compressing {fq} with pigz...")
                os.system(f"pigz -p 24 {fq}")  # pigz로 24개의 CPU 코어를 사용하여 압축
                print(f"Compressed {fq} to {fq}.gz")

if __name__ == "__main__":
    download_data()


