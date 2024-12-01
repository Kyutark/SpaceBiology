import os
import pandas as pd
import io  # StringIO 사용을 위해 import

# 결과 디렉토리 경로
result_dir = "../results"

# 결과 디렉토리에서 _counts.txt로 끝나는 파일 검색
counts_files = [f for f in os.listdir(result_dir) if f.endswith("_counts.txt")]

# 파일 처리 함수
def process_counts_file(file_path, output_dir):
    # 파일 읽기
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # 주석 제거 (#로 시작하는 줄)
    data_lines = [line for line in lines if not line.startswith("#")]
    
    # 데이터프레임으로 변환
    data_str = "\n".join(data_lines)
    df = pd.read_csv(io.StringIO(data_str), sep="\t")
    
    # 필요한 열만 선택 (0열과 마지막 열)
    df = df.iloc[:, [0, -1]]
    
    # 0열의 'cds-' 부분 제거
    df.iloc[:, 0] = df.iloc[:, 0].str.replace('cds-', '', regex=False)
    
    # 파일 이름 처리
    base_name = os.path.basename(file_path)
    new_name = base_name.replace("_counts.txt", ".tabular")
    output_path = os.path.join(output_dir, new_name)
    
    # 새로운 파일로 저장
    df.to_csv(output_path, sep="\t", index=False, header=False)

# 결과 디렉토리에 저장
output_dir = result_dir
for counts_file in counts_files:
    file_path = os.path.join(result_dir, counts_file)
    process_counts_file(file_path, output_dir)

print("처리가 완료되었습니다.")

