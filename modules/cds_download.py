import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

# 경로 설정
data_dir = "../data"
results_dir = "../results"
mdsh_file = os.path.join(data_dir, "mdsh.csv")

# 결과 디렉토리 생성
os.makedirs(results_dir, exist_ok=True)

# 메타데이터 로드 및 organism- ref 매핑
metadata = pd.read_csv(mdsh_file)
organism_ref_map = metadata[['organism', 'ref']].drop_duplicates().set_index('ref')['organism'].to_dict()

# 각 ref에 대한 fasta 파일 생성
for ref, organism in organism_ref_map.items():
    fna_file = os.path.join(data_dir, f"{ref}.fna")
    gff_file = os.path.join(data_dir, f"{ref}.gff")
    output_fasta = os.path.join(results_dir, f"{organism.replace(' ', '_')}_cds.fasta")

    # fna 파일에서 모든 레코드를 가져오기
    genome_records = {record.id: record.seq for record in SeqIO.parse(fna_file, "fasta")}

    # GFF 파일 읽고 CDS 추출
    protein_records = []
    with open(gff_file, "r") as gff_file:
        for line in gff_file:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if parts[2] == "CDS":
                start = int(parts[3]) - 1  # GFF는 1-based, Python은 0-based index
                end = int(parts[4])
                strand = parts[6]
                attributes = parts[8]

                # protein_id 추출
                protein_id = None
                for attr in attributes.split(";"):
                    if attr.startswith("protein_id"):
                        protein_id = attr.split("=")[1]
                        break

                if protein_id:
                    # 해당 CDS가 포함된 시퀀스 ID 추출
                    seq_id = parts[0]
                    cds_seq = genome_records[seq_id][start:end]
                    
                    # 염기서열 변환 (역상보 서열 포함)
                    if strand == "-":
                        cds_seq = cds_seq.reverse_complement()
                    protein_seq = cds_seq.translate(to_stop=True)

                    # FASTA 레코드 생성
                    protein_record = SeqRecord(Seq(str(protein_seq)), id=protein_id, description="")
                    protein_records.append(protein_record)

    # FASTA 파일로 저장
    SeqIO.write(protein_records, output_fasta, "fasta")
    print(f"Protein sequences for {organism} saved to {output_fasta}")

