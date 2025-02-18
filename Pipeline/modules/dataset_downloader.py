import subprocess
from pathlib import Path
import pandas as pd
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

    genome_records = SeqIO.to_dict(SeqIO.parse(fna_file, "fasta"))
    protein_records = []

    with gff_file.open() as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue

            start = int(parts[3]) - 1  # GFF is 1-based, Python uses 0-based
            end = int(parts[4])
            strand = parts[6]
            attributes = {k: v for k, v in (attr.split("=") for attr in parts[8].split(";"))}

            protein_id = attributes.get("protein_id")
            if protein_id and parts[0] in genome_records:
                cds_seq = genome_records[parts[0]].seq[start:end]
                if strand == "-":
                    cds_seq = cds_seq.reverse_complement()
                protein_records.append(SeqRecord(cds_seq.translate(to_stop=True), id=protein_id, description=""))

    if protein_records:
        SeqIO.write(protein_records, faa_output, "fasta")
        print(f"Generated {faa_output}")
    else:
        print(f"No CDS sequences found for {ref}.")

def download_data(mdsh_csv=None):
    """ Download genome data and generate missing CDS FAA files. """
    script_dir = Path(__file__).resolve().parent
    mdsh_csv = mdsh_csv or script_dir.parent / "data" / "mdsh.csv"
    data_dir = script_dir.parent / "data"

    mdsh_df = pd.read_csv(mdsh_csv)
    unique_refs = mdsh_df['ref'].drop_duplicates().dropna().tolist()

    for ref in unique_refs:
        print(f"Downloading files for {ref}...")
        subprocess.run(["datasets", "download", "genome", "accession", ref, "--include", "gff3,genome,cds", "--filename", str(data_dir / f"{ref}.zip")])

        # Extract files directly into the data/ directory
        subprocess.run(["unzip", "-o", str(data_dir / f"{ref}.zip"), "-d", str(data_dir)])
        (data_dir / f"{ref}.zip").unlink()  # Remove unnecessary compressed files

        # Check for FAA file and generate if missing
        if not list(data_dir.glob("*.faa")):
            print(f"No FAA file found for {ref}, generating from GFF and FNA...")
            generate_cds_faa(ref, data_dir)

    # Download FASTQ files and compress with pigz
    runs = mdsh_df['run'].dropna().tolist()
    for run in runs:
        if run.startswith('SRR'):
            print(f"Downloading FASTQ files for {run}...")
            subprocess.run(["fastq-dump", "--outdir", str(data_dir), "--split-files", run])

            for fq in data_dir.glob(f"{run}_*.fastq"):
                print(f"Compressing {fq} with pigz...")
                subprocess.run(["pigz", "-p", "24", str(fq)])

if __name__ == "__main__":
    download_data()
