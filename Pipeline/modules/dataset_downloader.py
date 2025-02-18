import os
import glob
import pandas as pd
import shutil
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def generate_cds_faa(ref, data_dir):
    """ Generate a CDS protein sequence file (.faa) from GFF and FNA files if missing. """
    fna_file = data_dir / f"{ref}.fna"
    gff_file = data_dir / f"{ref}.gff"
    faa_output = data_dir / f"{ref}_translated_genome.cds.faa"

    if not fna_file.exists() or not gff_file.exists():
        print(f"Skipping {ref}: Required files (FNA/GFF) missing.")
        return

    genome_records = {record.id: record.seq for record in SeqIO.parse(fna_file, "fasta")}
    protein_records = []

    with open(gff_file, "r") as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue

            start = int(parts[3]) - 1  # GFF is 1-based, Python uses 0-based
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            protein_id = None
            for attr in attributes.split(";"):
                if attr.startswith("protein_id"):
                    protein_id = attr.split("=")[1]
                    break

            if protein_id and parts[0] in genome_records:
                cds_seq = genome_records[parts[0]][start:end]
                if strand == "-":
                    cds_seq = cds_seq.reverse_complement()
                protein_seq = cds_seq.translate(to_stop=True)
                protein_records.append(SeqRecord(Seq(str(protein_seq)), id=protein_id, description=""))

    if protein_records:
        SeqIO.write(protein_records, faa_output, "fasta")
        print(f"Generated {faa_output}")
    else:
        print(f"No CDS sequences found for {ref}.")

def download_data(mdsh_csv=None):
    """ Download genome data and generate missing CDS FAA files. """
    script_dir = Path(__file__).resolve().parent
    if mdsh_csv is None:
        mdsh_csv = script_dir.parent / "data" / "mdsh.csv"
    data_dir = script_dir.parent / "data"

    mdsh_df = pd.read_csv(mdsh_csv)
    unique_refs = mdsh_df['ref'].drop_duplicates().dropna().tolist()

    for ref in unique_refs:
        print(f"Downloading files for {ref}...")
        os.system(f"datasets download genome accession {ref} --include gff3,genome,cds --filename {data_dir}/{ref}.zip")

        # 압축 해제 시 바로 data/로 저장
        os.system(f"unzip -o {data_dir}/{ref}.zip -d {data_dir}")

        # 불필요한 압축 파일 삭제
        os.remove(f"{data_dir}/{ref}.zip")

        # FAA 파일 확인 후 생성 여부 결정
        faa_files = glob.glob(f"{data_dir}/*.faa")
        if not faa_files:
            print(f"No FAA file found for {ref}, generating from GFF and FNA...")
            generate_cds_faa(ref, data_dir)

    # Download FASTQ files and compress with pigz
    runs = mdsh_df['run'].dropna().tolist()
    for run in runs:
        if run.startswith('SRR'):
            print(f"Downloading FASTQ files for {run}...")
            os.system(f"fastq-dump --outdir {data_dir} --split-files {run}")

            fastq_files = glob.glob(f"{data_dir}/{run}_*.fastq")
            for fq in fastq_files:
                print(f"Compressing {fq} with pigz...")
                os.system(f"pigz -p 24 {fq}")

if __name__ == "__main__":
    download_data()
