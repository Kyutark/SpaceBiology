import os
import pandas as pd
import subprocess

# 파일 경로 설정
metadata_path = "../data/mdsh.csv"
results_dir = "../results/"
relabelled_files_pattern = "relabelled_*.tabular"


def load_metadata():
    """ 메타데이터 로드 및 조건 검증 """
    metadata = pd.read_csv(metadata_path)
    experiments = metadata['experiment'].unique()
    return metadata, experiments


def prepare_experiment_files(metadata, experiment):
    """ 실험에 따른 control과 treatment 파일 분리 및 유효성 검사 """
    exp_data = metadata[metadata['experiment'] == experiment]
    conditions = exp_data['condition'].unique()

    # 조건이 두 개로 나뉘는지 확인
    if len(conditions) != 2:
        raise ValueError(f"Experiment {experiment} has {len(conditions)} conditions, expected 2.")

    # 각 조건별 run 추출
    control_runs = exp_data[exp_data['condition'] == conditions[0]]['run'].tolist()
    treatment_runs = exp_data[exp_data['condition'] == conditions[1]]['run'].tolist()

    # 각 run에 대한 파일 경로 설정
    control_files = [os.path.join(results_dir, f"relabelled_{run}.tabular") for run in control_runs]
    treatment_files = [os.path.join(results_dir, f"relabelled_{run}.tabular") for run in treatment_runs]

    return control_files, treatment_files


def save_r_script(control_files, treatment_files, output_file):
    """ R 스크립트를 작성하여 logFC 및 p-value 계산 """
    r_script = f"""
    library(edgeR)
    control_files <- c({', '.join([f'"{file}"' for file in control_files])})
    treatment_files <- c({', '.join([f'"{file}"' for file in treatment_files])})

    read_counts_file <- function(filepath) {{
        df <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        colnames(df) <- c("geneID", "count")
        return(df)
    }}

    control_data <- lapply(control_files, read_counts_file)
    treatment_data <- lapply(treatment_files, read_counts_file)

    control_counts <- Reduce(function(x, y) merge(x, y, by = "geneID", all = TRUE), control_data)
    treatment_counts <- Reduce(function(x, y) merge(x, y, by = "geneID", all = TRUE), treatment_data)

    control_counts[is.na(control_counts)] <- 0
    treatment_counts[is.na(treatment_counts)] <- 0

    all_counts <- merge(control_counts, treatment_counts, by = "geneID")
    rownames(all_counts) <- all_counts$geneID
    all_counts <- all_counts[, -1]

    group <- factor(c(rep("control", length(control_files)), rep("treatment", length(treatment_files))))
    dge <- DGEList(counts = all_counts, group = group)
    dge <- calcNormFactors(dge)

    design <- model.matrix(~group)
    dge <- estimateDisp(dge, design)
    fit <- glmFit(dge, design)
    lrt <- glmLRT(fit)

    deg_results <- topTags(lrt, n = nrow(dge))$table
    write.table(deg_results[, c("logFC", "PValue")], file = "{output_file}", sep = "\t", row.names = TRUE, quote = FALSE)
    """

    with open("deg_analysis.R", "w") as f:
        f.write(r_script)


def run_deg_analysis():
    """ 실험 별로 DEG 분석 수행 """
    metadata, experiments = load_metadata()
    for experiment in experiments:
        try:
            # 파일 분리 및 조건 확인
            control_files, treatment_files = prepare_experiment_files(metadata, experiment)
            output_file = os.path.join(results_dir, f"DEG_results_{experiment}.txt")

            # R 스크립트 생성 및 실행
            save_r_script(control_files, treatment_files, output_file)
            subprocess.run(["Rscript", "deg_analysis.R"], check=True)
            print(f"DEG analysis for experiment {experiment} completed. Results saved to {output_file}")
        except ValueError as e:
            print(e)


if __name__ == "__main__":
    run_deg_analysis()
