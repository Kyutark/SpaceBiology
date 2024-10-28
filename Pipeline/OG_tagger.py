import pandas as pd
import os


def OG_tagger(mmseqs_file, tabular_dir):
    mmseqs_df = pd.read_csv(mmseqs_file, sep='\t', header=None, names=['geneID', 'OGID'])

    # 중복 geneID 제거 (첫 번째 값만 사용)
    mmseqs_df = mmseqs_df.drop_duplicates(subset=['geneID'], keep='first')

    # 디렉토리 내의 모든 .tabular 파일 처리
    for filename in os.listdir(tabular_dir):
        if filename.endswith('.tabular'):
            srr_file = os.path.join(tabular_dir, filename)
            srr_df = pd.read_csv(srr_file, sep='\t', header=None)
            srr_df['OGID'] = srr_df[0].map(mmseqs_df.set_index('geneID')['OGID'])
            srr_df = srr_df.dropna(subset=['OGID'])

            # 중복 OGID 처리
            srr_df = srr_df.groupby('OGID').agg({1: 'sum'}).reset_index()

            output_file = os.path.join(tabular_dir, f"OG_{filename}")
            srr_df.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    mmseqs_file = "clusterRes30.tsv"
    tabular_dir = "."  # .tabular 파일이 있는 디렉토리
    OG_tagger(mmseqs_file, tabular_dir)
