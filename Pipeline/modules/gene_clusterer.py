import os
import subprocess
import pandas as pd
import logging
import shutil

# Set up logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def merge_fasta_files(input_dir, output_file):
    """Merge all _cds.fasta files into a single FASTA file."""
    fasta_files = [f for f in os.listdir(input_dir) if f.endswith("_cds.fasta")]
    if not fasta_files:
        raise FileNotFoundError("_cds.fasta 파일을 찾을 수 없습니다. 병합할 파일이 필요합니다.")

    logging.info(f"Found {len(fasta_files)} FASTA files. Merging into {output_file}...")
    with open(output_file, "w") as outfile:
        for fasta_file in fasta_files:
            with open(os.path.join(input_dir, fasta_file), "r") as infile:
                outfile.write(infile.read())

def run_mmseqs_clustering(merged_fasta, output_dir):
    """Run MMseqs2 to perform clustering."""
    seq_db = os.path.join(output_dir, "seqDB")
    cluster_res = os.path.join(output_dir, "clusterRes")
    tmp_dir = os.path.join(output_dir, "tmp")
    tsv_file = os.path.join(output_dir, "clusterRes_named.tsv")

    # Remove any existing files/directories
    for path in [seq_db, cluster_res, tmp_dir, tsv_file]:
        if os.path.exists(path):
            if os.path.isfile(path):
                os.remove(path)
            else:
                shutil.rmtree(path)

    # Run MMseqs2 commands
    subprocess.run(["mmseqs", "createdb", merged_fasta, seq_db], check=True)
    subprocess.run(["mmseqs", "cluster", seq_db, cluster_res, tmp_dir, "--min-seq-id", "0.3", "-c", "0.8"], check=True)
    subprocess.run(["mmseqs", "createtsv", seq_db, seq_db, cluster_res, tsv_file], check=True)

    logging.info(f"Clustering completed. TSV file saved at {tsv_file}")
    return tsv_file

def relabel_tabular_files(input_dir, tsv_file, relabelled_prefix="relabelled_"):
    """Re-label `.tabular` files using MMseqs2 clustering results."""
    mmseqs_df = pd.read_csv(tsv_file, sep="\t", header=None, names=["Representative", "Member"])
    mapping = dict(zip(mmseqs_df["Member"], mmseqs_df["Representative"]))

    tabular_files = [f for f in os.listdir(input_dir) if f.endswith(".tabular") and not f.startswith(relabelled_prefix)]
    for tab_file in tabular_files:
        tabular_path = os.path.join(input_dir, tab_file)
        tabular_df = pd.read_csv(tabular_path, sep="\t", header=None, names=["CDS_ID", "Gene_Count"])

        # Map CDS_ID to Representative
        tabular_df["Representative"] = tabular_df["CDS_ID"].map(mapping)
        tabular_df = tabular_df.dropna(subset=["Representative"])  # Drop rows with no mapping

        # Aggregate counts by Representative
        aggregated_df = tabular_df.groupby("Representative")["Gene_Count"].sum().reset_index()

        # Save relabelled file
        relabelled_file = os.path.join(input_dir, relabelled_prefix + tab_file)
        aggregated_df.to_csv(relabelled_file, sep="\t", header=False, index=False)
        logging.info(f"Processed {tab_file}: saved relabelled file as {relabelled_file}")

def main():
    input_dir = "/home/synbio/Kyutark/analysis_pipeline/results"
    output_dir = "/home/synbio/Kyutark/analysis_pipeline/results/mmseqs_output"
    merged_fasta = os.path.join(output_dir, "merged.fasta")
    relabelled_prefix = "relabelled_"

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: Merge FASTA files
    merge_fasta_files(input_dir, merged_fasta)

    # Step 2: Run MMseqs2 clustering
    tsv_file = run_mmseqs_clustering(merged_fasta, output_dir)

    # Step 3: Re-label .tabular files
    relabel_tabular_files(input_dir, tsv_file, relabelled_prefix)

if __name__ == "__main__":
    main()

