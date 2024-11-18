import os
import pandas as pd

# 경로 설정
result_dir = "../results"
orthofinder_output_dir = "../results/orthofinder_output"
relabelled_prefix = "relabelled_"

# Orthogroups.tsv 파일 경로
orthogroups_path = None
for root, dirs, files in os.walk(orthofinder_output_dir):
    if "Orthogroups.tsv" in files:
        orthogroups_path = os.path.join(root, "Orthogroups.tsv")
        break

if not orthogroups_path:
    raise FileNotFoundError("Orthogroups.tsv 파일을 찾을 수 없습니다.")

# Orthogroups.tsv 파일 로드 및 처리
orthogroups_df = pd.read_csv(orthogroups_path, sep="\t", index_col=0)

# OGID와 CDS_ID 매핑을 딕셔너리로 변환
ogid_mapping = {}
for og, row_data in orthogroups_df.iterrows():
    organism_cds_ids = row_data.dropna().tolist()
    for cds_id in organism_cds_ids:
        ogid_mapping[cds_id] = og

# 모든 .tabular 파일에 대해 처리
tabular_files = [f for f in os.listdir(result_dir) if f.endswith(".tabular") and not f.startswith(relabelled_prefix)]

for tab_file in tabular_files:
    # .tabular 파일 읽기
    tabular_path = os.path.join(result_dir, tab_file)
    tabular_df = pd.read_csv(tabular_path, sep="\t", header=None, names=["CDS_ID", "Gene_Count"])

    # OGID 대체 및 필터링
    unique_counts = {}

    for idx, row in tabular_df.iterrows():
        cds_id = row["CDS_ID"]
        count = row["Gene_Count"]

        # OGID 찾기
        ogid = ogid_mapping.get(cds_id, None)

        # OGID가 없는 경우 해당 행 제외
        if ogid is None:
            continue

        # OGID가 있는 경우
        if ogid in unique_counts:
            unique_counts[ogid] += count
        else:
            unique_counts[ogid] = count

    # 결과를 데이터프레임으로 정리
    relabelled_df = pd.DataFrame(list(unique_counts.items()), columns=["OGID", "Gene_Count"])

    # relabelled 파일로 저장
    relabelled_filename = relabelled_prefix + tab_file
    relabelled_path = os.path.join(result_dir, relabelled_filename)
    relabelled_df.to_csv(relabelled_path, sep="\t", header=False, index=False)

    print(f"{tab_file} 처리 완료: {relabelled_filename}으로 저장되었습니다.")

