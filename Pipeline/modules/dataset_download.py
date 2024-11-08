from pathlib import Path
import os
import glob
import pandas as pd

def download_data(mdsh_csv=None):
    # mdsh_csv 경로 설정
    if mdsh_csv is None:
        script_dir = Path(__file__).resolve().parent
        mdsh_csv = script_dir.parent / "data" / "mdsh.csv"

    mdsh_df = pd.read_csv(mdsh_csv)
    unique_refs = mdsh_df['ref'].drop_duplicates().dropna().tolist()
    data_dir = script_dir.parent / "data"  # data 디렉토리 절대 경로 지정

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

    # Run 열에 있는 파일 다운로드 또는 복사
    runs = mdsh_df['run'].dropna().tolist()
    for run in runs:
        if run.startswith('SRR'):  # SRA ID인 경우 다운로드
            print(f"Downloading FASTQ files for {run}...")
            os.system(f"fastq-dump --outdir {data_dir} --split-files {run}")
            if os.path.isfile(f"{data_dir}/{run}.fastq"):
                os.rename(f"{data_dir}/{run}.fastq", f"{data_dir}/{run}_1.fastq")

if __name__ == "__main__":
    download_data()
