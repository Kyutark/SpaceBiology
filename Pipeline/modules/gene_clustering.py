import os
import glob
import pandas as pd
from Bio import Entrez, SeqIO

Entrez.email = "kyutkim01@korea.ac.kr"

os.makedirs("../results/mmseqs2", exist_ok=True)

def pre_clustering():
    """ 모든 *_counts.txt 파일에서 Geneid와 마지막 열을 추출하고 .tabular 파일로 저장 """
    count_files = glob.glob("../results/*_counts.txt")
    for count_file in count_files:
        base_name = os.path.basename(count_file).replace("_counts.txt", "")
        output_file = f"../results/{base_name}.tabular"
        df = pd.read_csv(count_file, sep="\t", comment="#")
        df = df[['Geneid', df.columns[-1]]]  # Geneid와 마지막 열만 남김
        df['Geneid'] = df['Geneid'].str.replace("^cds-", "", regex=True)  # 'cds-' 접두사 제거
        df.columns = ['Geneid', base_name]  # 마지막 열의 이름을 SRR ID로 변경
        df.to_csv(output_file, sep="\t", index=False)

def gather_unique_gene_ids():
    """ 모든 .tabular 파일에서 GeneID를 수집하여 중복 없이 정리한 후 저장 """
    tabular_files = glob.glob("../results/*.tabular")
    gene_ids = set()

    for file in tabular_files:
        df = pd.read_csv(file, sep="\t")
        # gene ID에서 - 이후 부분 제거
        cleaned_ids = [gene_id.split('-')[0] for gene_id in df['Geneid'].unique()]
        gene_ids.update(cleaned_ids)

    with open("../results/mmseqs2/unique_gene_ids.txt", "w") as f:
        f.write("\n".join(gene_ids))
    return gene_ids

def download_sequences(gene_ids, output_fasta="../results/mmseqs2/all_genes.fasta"):
    with open(output_fasta, "w") as output_handle:
        for gene_id in gene_ids:
            try:
                handle = Entrez.efetch(db="protein", id=gene_id, rettype="fasta", retmode="text")
                seq_record = SeqIO.read(handle, "fasta")
                SeqIO.write(seq_record, output_handle, "fasta")
                handle.close()
                print(f"Downloaded {gene_id}")
            except Exception as e:
                print(f"Error downloading {gene_id}: {e}")

def run_mmseqs_clustering(input_file="../results/mmseqs2/all_genes.fasta", output_dir="../results/mmseqs2/mmseqs_clustering"):
    """ MMseqs2를 이용하여 gene clustering 수행 """
    # output 디렉토리 생성
    os.makedirs(output_dir, exist_ok=True)

    # MMseqs2 clustering pipeline
    os.system(f"mmseqs createdb {input_file} {output_dir}/genesDB")
    os.system(f"mmseqs cluster {output_dir}/genesDB {output_dir}/genesDB_clu tmp --min-seq-id 0.3 -c 0.8")
    os.system(
        f"mmseqs createtsv {output_dir}/genesDB {output_dir}/genesDB {output_dir}/genesDB_clu ../results/genegroups_mmseqs.txt")

def cluster_genes():
    """ 전체 클러스터링 파이프라인 실행 """
    pre_clustering()
    gene_ids = gather_unique_gene_ids()
    download_sequences(gene_ids)
    run_mmseqs_clustering()

if __name__ == "__main__":
    cluster_genes()