from pathlib import Path
import subprocess
import pandas as pd

def process_rna_seq():
    """ Process RNA-seq data using fastp, bowtie2, samtools, and featureCounts. """
    # Set paths
    script_dir = Path(__file__).resolve().parent  # /module
    data_dir = script_dir.parent / "data"  # /data
    result_dir = script_dir.parent / "result"  # /result

    # Load metadata
    mdsh_csv = data_dir / "mdsh.csv"
    if not mdsh_csv.exists():
        raise FileNotFoundError(f"File {mdsh_csv} does not exist!")

    mdsh_df = pd.read_csv(mdsh_csv)
    runs = mdsh_df['run'].dropna().tolist()

    # Create result directory
    result_dir.mkdir(parents=True, exist_ok=True)

    # Process each RNA-seq run
    for run in runs:
        # Define file paths
        ref = mdsh_df.loc[mdsh_df['run'] == run, 'ref'].values[0]
        fq1 = data_dir / f"{run}_1.fastq.gz"
        fq2 = data_dir / f"{run}_2.fastq.gz" if (data_dir / f"{run}_2.fastq.gz").exists() else None
        fna_file = data_dir / f"{ref}.fna"
        gff_file = data_dir / f"{ref}.gff"

        # Check if required files exist
        if not fq1.exists():
            print(f"Error: FASTQ file {fq1} does not exist for run {run}")
            continue

        print(f"Processing RNA-seq for run: {run} with reference: {ref}")

        # Quality control with fastp
        if fq2:  # Paired-end
            subprocess.run(["fastp", "-i", fq1, "-I", fq2, "-o", result_dir / f"{run}_filtered_1.fastq.gz", "-O", result_dir / f"{run}_filtered_2.fastq.gz"], check=True)
        else:  # Single-end
            subprocess.run(["fastp", "-i", fq1, "-o", result_dir / f"{run}_filtered.fastq.gz"], check=True)

        # Build Bowtie2 index
        subprocess.run(["bowtie2-build", fna_file, result_dir / f"{run}_index"], check=True)

        # Align reads with Bowtie2 and convert to sorted BAM
        sorted_bam_file = result_dir / f"{run}_sorted.bam"
        if fq2:  # Paired-end
            subprocess.run(
                f"bowtie2 -p 24 -x {result_dir}/{run}_index -1 {result_dir}/{run}_filtered_1.fastq.gz -2 {result_dir}/{run}_filtered_2.fastq.gz | "
                f"samtools view -bS | samtools sort -o {sorted_bam_file}",
                shell=True,
                check=True
            )
        else:  # Single-end
            subprocess.run(
                f"bowtie2 -p 24 -x {result_dir}/{run}_index -U {result_dir}/{run}_filtered.fastq.gz | "
                f"samtools view -bS | samtools sort -o {sorted_bam_file}",
                shell=True,
                check=True
            )

        # FeatureCounts for gene expression quantification
        counts_file = result_dir / f"{run}_counts.txt"
        if fq2:  # Paired-end
            subprocess.run(["featureCounts", "-a", gff_file, "-o", counts_file, "-t", "CDS", "-g", "ID", "-p", sorted_bam_file], check=True)
        else:  # Single-end
            subprocess.run(["featureCounts", "-a", gff_file, "-o", counts_file, "-t", "CDS", "-g", "ID", sorted_bam_file], check=True)

if __name__ == "__main__":
    process_rna_seq()
