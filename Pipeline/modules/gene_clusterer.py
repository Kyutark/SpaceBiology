import re
import subprocess
import shutil
from pathlib import Path
from Bio import SeqIO

def extract_protein_ids_from_gff(gff_file):
    """Extracts mapping of 'ID=cds-XXXXX;' to 'XXXXX' from a GFF file."""
    protein_id_map = {}

    with open(gff_file, "r") as file:
        for line in file:
            match = re.search(r"ID=cds-([\w\.]+?);", line)  # Stop at `;`
            if match:
                protein_id = match.group(1)  # Extracted protein ID
                protein_id_map[protein_id] = protein_id  # Store mapping

    return protein_id_map

def rename_fasta_headers(input_dir, gff_file, output_faa):
    """Renames FAA headers based on protein IDs extracted from GFF."""
    input_dir = Path(input_dir)
    faa_files = list(input_dir.glob("*.faa"))
    if not faa_files:
        raise FileNotFoundError("No .faa files found in the directory.")

    # Extract protein ID mapping from GFF
    protein_id_map = extract_protein_ids_from_gff(gff_file)

    with output_faa.open("w") as outfile:
        for faa_file in faa_files:
            for record in SeqIO.parse(faa_file, "fasta"):
                protein_id = protein_id_map.get(record.id)  # Match FASTA ID to GFF protein ID
                
                if protein_id:  # Only rename if a match is found
                    record.id = protein_id
                    record.description = ""  # Remove extra description
                    SeqIO.write(record, outfile, "fasta")

    print(f"Renamed protein sequences saved to {output_faa}")

def run_mmseqs_clustering(seq_db, output_dir):
    """Runs MMseqs2 clustering and retains only the final .tsv result."""
    cluster_res = output_dir / "clusterRes"
    tmp_dir = output_dir / "tmp"
    tsv_file = output_dir / "clusterRes_named.tsv"

    # Run MMseqs2 clustering
    subprocess.run(["mmseqs", "cluster", seq_db, cluster_res, tmp_dir, "--min-seq-id", "0.3", "-c", "0.8"], check=True)
    subprocess.run(["mmseqs", "createtsv", seq_db, seq_db, cluster_res, tsv_file], check=True)

    print(f"Clustering completed. TSV file saved at {tsv_file}")

    # Cleanup: Remove all temporary MMseqs2 files except for the final .tsv
    for path in [seq_db, cluster_res, tmp_dir]:
        if path.exists():
            shutil.rmtree(path)

    return tsv_file

def main():
    # Set up paths
    script_dir = Path(__file__).resolve().parent
    data_dir = script_dir.parent / "data"  # ../data
    result_dir = script_dir.parent / "result" / "mmseqs_output"  # ../result/mmseqs_output

    # Ensure output directory exists
    result_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Extract proteins from GFF and rename FAA headers
    renamed_faa_files = []
    
    for gff_file in data_dir.glob("*.gff"):  # Process each GFF file
        ref_id = gff_file.stem
        faa_output = result_dir / f"{ref_id}_renamed.faa"
        
        rename_fasta_headers(data_dir, gff_file, faa_output)
        renamed_faa_files.append(faa_output)

    # Step 2: Merge all renamed FAA files into one for clustering
    merged_fasta = result_dir / "merged_renamed_proteins.faa"
    
    with merged_fasta.open("w") as merged_file:
        for faa_file in renamed_faa_files:
            with faa_file.open("r") as single_faa:
                merged_file.write(single_faa.read())

    # Step 3: Run MMseqs2 clustering on the merged FAA file
    seq_db = result_dir / "seqDB"
    subprocess.run(["mmseqs", "createdb", merged_fasta, seq_db], check=True)

    run_mmseqs_clustering(seq_db, result_dir)

if __name__ == "__main__":
    main()
