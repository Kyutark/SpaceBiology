import os
import pandas as pd


def process_rna_seq():
    # mdsh.csv에서 run과 ref 정보를 읽어옴
    mdsh_df = pd.read_csv("../data/mdsh.csv")
    runs = mdsh_df['run'].dropna().tolist()

    # 결과 디렉토리 생성
    os.makedirs("../results", exist_ok=True)

    # 각 run에 대해 RNA-seq 처리 실행
    for run in runs:
        ref = mdsh_df.loc[mdsh_df['run'] == run, 'ref'].values[0]
        fq1 = f"../data/{run}_1.fastq"
        fq2 = f"../data/{run}_2.fastq" if os.path.isfile(f"../data/{run}_2.fastq") else None
        fna_file = f"../data/{ref}.fna"
        gff_file = f"../data/{ref}.gff"

        print(f"Processing RNA-seq for run: {run} with reference: {ref}")

        # 품질 제어
        if fq2:
            os.system(f"fastp -i {fq1} -I {fq2} -o ../results/{run}_filtered_1.fastq -O ../results/{run}_filtered_2.fastq")
        else:
            os.system(f"fastp -i {fq1} -o ../results/{run}_filtered.fastq")

        # Bowtie2 인덱스 생성 및 정렬
        os.system(f"bowtie2-build {fna_file} ../results/{run}_index")
        if fq2:
            os.system(
                f"bowtie2 -x ../results/{run}_index -1 ../results/{run}_filtered_1.fastq -2 ../results/{run}_filtered_2.fastq -S ../results/{run}.sam")
        else:
            os.system(f"bowtie2 -x ../results/{run}_index -U ../results/{run}_filtered.fastq -S ../results/{run}.sam")

        # FeatureCounts를 통한 카운팅 - gff 파일에서 'ID' 속성 사용
        # FeatureCounts 실행
        if fq2:
            # Paired-end인 경우 -p 옵션 추가
            os.system(
                f"featureCounts -a {gff_file} -o ../results/{run}_counts.txt -t CDS -g ID -p ../results/{run}.sam")
        else:
            # Single-end인 경우 -p 옵션 없이 실행
            os.system(f"featureCounts -a {gff_file} -o ../results/{run}_counts.txt -t CDS -g ID ../results/{run}.sam")


if __name__ == "__main__":
    process_rna_seq()
