import pandas as pd
import io
from pathlib import Path

# Define result directory based on script location
script_dir = Path(__file__).resolve().parent
result_dir = script_dir.parent / "results"

# Find all _counts.txt files in the result directory
counts_files = list(result_dir.glob("*_counts.txt"))

def process_counts_file(file_path, output_dir):
    """Process a featureCounts output file and save it in tabular format."""
    
    # Read the file and remove comment lines
    with file_path.open("r") as file:
        data_lines = [line for line in file.read().splitlines() if not line.startswith("#")]
    
    # Convert to DataFrame
    data_str = "\n".join(data_lines)
    df = pd.read_csv(io.StringIO(data_str), sep="\t")

    # Select only the first and last columns
    df = df.iloc[:, [0, -1]].copy()  # Ensure DataFrame is not a view

    # Remove 'cds-' prefix from the first column
    df.iloc[:, 0] = df.iloc[:, 0].str.replace('cds-', '', regex=False)

    # Generate new filename
    new_name = file_path.name.replace("_counts.txt", ".tabular")
    output_path = output_dir / new_name

    # Save to a new file
    df.to_csv(output_path, sep="\t", index=False, header=False)

# Process all found counts files
for counts_file in counts_files:
    process_counts_file(counts_file, result_dir)

print("Processing completed.")
