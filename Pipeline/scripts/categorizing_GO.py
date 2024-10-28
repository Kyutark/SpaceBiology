import pandas as pd
import requests
from xml.etree import ElementTree as ET
from Bio import Entrez
import time

# NCBI에서 정보를 가져오기 위해 이메일 설정
Entrez.email = "your.email@example.com"

# 엑셀 파일 로드
file_path = "C:\\Users\\Kyutark Kim\\PycharmProjects\\re-labelling2\\validOG_list.xlsx"
df = pd.read_excel(file_path)

# 0열의 단백질 ID 추출 (0행은 인덱스이므로 무시)
protein_ids = df.iloc[:, 0].tolist()


# RefSeq ID를 UniProt Accession으로 변환
def get_uniprot_accession(refseq_id):
    url = "https://rest.uniprot.org/idmapping/run"
    params = {
        'from': 'RefSeq_Protein',
        'to': 'UniProtKB',
        'ids': refseq_id
    }
    response = requests.post(url, data=params)
    response.raise_for_status()
    job_id = response.json()['jobId']

    # 결과를 가져오기 위해 일정 시간 대기 후 반복 확인
    result_url = f"https://rest.uniprot.org/idmapping/results/{job_id}"
    while True:
        result_response = requests.get(result_url)
        if result_response.status_code == 200:
            result_data = result_response.json()
            if 'results' in result_data and len(result_data['results']) > 0:
                return result_data['results'][0]['to']
        elif result_response.status_code == 404:
            time.sleep(5)
        else:
            result_response.raise_for_status()

    return None


# 단백질 ID로부터 GO terms 가져오기 함수 정의
def get_go_terms_uniprot(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
    response = requests.get(url)
    response.raise_for_status()

    root = ET.fromstring(response.content)

    # GO terms 추출
    go_terms = []
    for entry in root.findall(".//{http://uniprot.org/uniprot}dbReference[@type='GO']"):
        for prop in entry.findall(".//{http://uniprot.org/uniprot}property[@type='term']"):
            go_terms.append(prop.attrib['value'])

    return go_terms


# 각 단백질 ID에 대해 UniProt Accession으로 변환 후 GO terms 조회 및 출력
for protein_id in protein_ids:
    try:
        uniprot_id = get_uniprot_accession(protein_id)
        if uniprot_id:
            go_terms = get_go_terms_uniprot(uniprot_id)
            print(f"Protein ID: {protein_id} (UniProt ID: {uniprot_id}), GO terms: {go_terms}")
        else:
            print(f"Protein ID: {protein_id} - UniProt ID를 찾을 수 없습니다.")
    except requests.exceptions.HTTPError as e:
        print(f"HTTP Error for Protein ID: {protein_id} - {str(e)}")
