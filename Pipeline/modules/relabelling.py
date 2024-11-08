import os
import glob
import pandas as pd


def load_genegroups(filepath):
    """Load gene group mappings from genegroups_mmseqs.txt."""
    genegroups = pd.read_csv(filepath, sep="\t", header=None, names=["genegroup", "geneID"])
    genegroups["geneID"] = genegroups["geneID"].str.split("-").str[0]  # Remove subid part (e.g., -1)
    return genegroups.set_index("geneID")["genegroup"].to_dict()


def relabel_and_aggregate_tabular_files(genegroups, results_dir="../results"):
    """Relabel gene IDs in .tabular files, aggregate counts for duplicate genegroups, and save relabelled files."""
    tabular_files = glob.glob(os.path.join(results_dir, "*.tabular"))

    for file_path in tabular_files:
        df = pd.read_csv(file_path, sep="\t")

        # Remove subid part from gene IDs in the .tabular file for matching
        df["CleanedGeneID"] = df["Geneid"].str.split("-").str[0]

        # Map gene IDs to genegroups and remove unmatched rows
        df["genegroup"] = df["CleanedGeneID"].map(genegroups)
        df = df.dropna(subset=["genegroup"])  # Drop rows without matching genegroup
        df = df.drop(columns=["CleanedGeneID"])  # Drop temporary column

        # Replace Geneid with genegroup
        df["Geneid"] = df["genegroup"]
        df = df.drop(columns=["genegroup"])  # Drop genegroup column after replacing

        # Aggregate rows by genegroup, summing the counts for duplicates
        df = df.groupby("Geneid", as_index=False).sum()

        # Save the relabelled file
        base_name = os.path.basename(file_path)
        output_file = os.path.join(results_dir, f"relabelled_{base_name}")
        df.to_csv(output_file, sep="\t", index=False)


if __name__ == "__main__":
    # Load gene group mappings
    genegroups_file = "../results/genegroups_mmseqs.txt"
    genegroups = load_genegroups(genegroups_file)

    # Relabel and aggregate .tabular files
    relabel_and_aggregate_tabular_files(genegroups)
