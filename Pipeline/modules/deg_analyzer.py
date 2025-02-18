import subprocess
import re
from pathlib import Path
import pandas as pd

# Define paths relative to the script location
script_dir = Path(__file__).resolve().parent
metadata_path = script_dir.parent / "data" / "mdsh.csv"
results_dir = script_dir.parent / "results"
error_log_path = results_dir / "error_log.txt"

def sanitize_filename(name):
    """Sanitize filenames by replacing special characters with underscores."""
    return re.sub(r'[<>:"/\\|?*]', '_', name)

def load_metadata():
    """Load metadata from the CSV file."""
    metadata = pd.read_csv(metadata_path)
    experiments = metadata['experiment'].unique()
    return metadata, experiments

def prepare_experiment_files(metadata, experiment):
    """Retrieve control and treatment files for a given experiment."""
    exp_data = metadata[metadata['experiment'] == experiment]

    # Extract run IDs for each condition
    control_runs = exp_data[exp_data['condition'] == "control"]['run'].tolist()
    treatment_runs = exp_data[exp_data['condition'] == "treatment"]['run'].tolist()

    # Construct file paths (using Path objects for better compatibility)
    control_files = [results_dir / f"relabelled_{run}.tabular" for run in control_runs]
    treatment_files = [results_dir / f"relabelled_{run}.tabular" for run in treatment_runs]

    return control_files, treatment_files

def save_r_script(control_files, treatment_files, output_file):
    """Generate an R script to perform logFC and p-value calculations."""
    r_script_path = script_dir / "deg_analysis.R"  # Ensure the R script is saved in the module directory

    r_script = f"""
    library(edgeR)

    # Define file paths
    control_files <- c({', '.join([f'"{file}"' for file in control_files])})
    treatment_files <- c({', '.join([f'"{file}"' for file in treatment_files])})

    # Function to read count files
    read_counts_file <- function(filepath) {{
        df <- read.table(filepath, header = FALSE, sep = "\\t")
        colnames(df) <- c("geneID", "count")
        return(df)
    }}

    # Read count data
    control_data <- lapply(control_files, read_counts_file)
    treatment_data <- lapply(treatment_files, read_counts_file)

    # Merge all count data in a single step
    all_counts <- Reduce(function(x, y) merge(x, y, by = "geneID", all = TRUE), c(control_data, treatment_data))
    all_counts[is.na(all_counts)] <- 0  # Replace NA values with 0

    # Set row names and remove geneID column
    rownames(all_counts) <- all_counts$geneID
    all_counts <- all_counts[, -1]

    # Define experimental groups
    group <- factor(c(rep("control", length(control_files)), rep("treatment", length(treatment_files))))
    dge <- DGEList(counts = all_counts, group = group)
    dge <- calcNormFactors(dge)

    # Model fitting and differential expression analysis
    design <- model.matrix(~group)
    dge <- estimateDisp(dge, design)
    fit <- glmFit(dge, design)
    lrt <- glmLRT(fit)

    # Save logFC and P-value results
    deg_results <- topTags(lrt, n = Inf)$table  # `n = Inf` retrieves all results
    write.table(deg_results[, c("logFC", "PValue")], file = "{output_file}", sep = "\\t", row.names = TRUE, quote = FALSE)
    """

    with r_script_path.open("w") as f:
        f.write(r_script)

def run_deg_analysis():
    """Perform DEG analysis for each experiment."""
    metadata, experiments = load_metadata()
    for experiment in experiments:
        safe_experiment = sanitize_filename(experiment)
        try:
            # Prepare experiment files
            control_files, treatment_files = prepare_experiment_files(metadata, experiment)
            output_file = results_dir / f"DEG_results_{safe_experiment}.txt"

            # Generate and execute the R script
            save_r_script(control_files, treatment_files, output_file)
            subprocess.run(["Rscript", str(script_dir / "deg_analysis.R")], check=True)
            print(f"DEG analysis for experiment {experiment} completed. Results saved to {output_file}")
        except Exception as e:
            error_message = f"Error in experiment {experiment}: {e}\n"
            print(error_message)
            with error_log_path.open("a") as log_file:
                log_file.write(error_message)

if __name__ == "__main__":
    run_deg_analysis()
