#read GTF 
#필요한 attribute 분리
#protein-coding gene만 filtering

import pandas as pd 

#gtf 파일을 읽어서 dataframe으로 저장 
def read_gtf(gtf_path):

    #progress tracking 
    print("[1/3] reading gtf file")

    columns = [
        "chromosome",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attributes",
    ]

    gtf = pd.read_csv(
        gtf_path,
        sep="\t", #column 구분 tab으로 
        comment="#", #hashtag으로 시작하는 줄은 무시
        names=columns, #datagram column 이름 지정
    )
    print(f"GTF loaded: {len(gtf):,} rows")

    return gtf

#GTF의 attributes column에서 
#gene_id, gene_type, gene_name, transcript_id 추출
def extract_attribute(gtf):

    gtf=gtf.copy() #원본 dataframe은 수정 X 

    gtf["gene_id"] = (
        gtf["attributes"]
        .str.extract(r'gene_id "([^"]+)"') #기호가 뭔지 알아보기
    )

    gtf["gene_type"] = (
        gtf["attributes"]
        .str.extract(r'gene_type "([^"]+)"')
    )

    gtf["gene_name"] = (
        gtf["attributes"]
        .str.extract(r'gene_name "([^"]+)"')
    )

    gtf["transcript_id"] = (
        gtf["attributes"]
        .str.extract(r'transcript_id "([^"]+)"')
    )

    return gtf

#protein coding gene에 속하는 annotation만 반환
def filter_protein_coding(gtf):
    protein_coding_gtf = gtf[
        gtf["gene_type"] == "protein_coding"
    ].copy()

    protein_coding_gtf = protein_coding_gtf.reset_index(
        drop=True
    )

    return protein_coding_gtf

def load_annotation(gtf_path):

    gtf = read_gtf(gtf_path)
    gtf = extract_attribute(gtf)
    gtf = filter_protein_coding(gtf)

    return gtf