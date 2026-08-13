#실행용 파일
#read -> extract -> filter

from annotation import (
    read_gtf,
    extract_attribute,
    filter_protein_coding
)

gtf_path = (
    r"C:\Users\User\Desktop\biglab\gc_project\data\annotation\gencode.v50.annotation.gtf"
)

#gtf 읽어서 df로 저장
gtf = read_gtf(gtf_path) 
#gtf attributes 열에서 
# gene_id/gene_type/gene_name/transcript_id 뽑아서 새로운 col로 추가 
gtf = extract_attribute(gtf) 
protein_coding_gtf = filter_protein_coding(gtf)

print("\n=== Feature counts ===")

print(
    protein_coding_gtf["feature"]
    .value_counts()
)