#go analysis 
import gseapy as gp
from top_bottom_run import top100, bottom100

gene_sets = [
    "GO_Biological_Process_2025",
    "GO_Molecular_Function_2025",
    "GO_Cellular_Component_2025",
]

# Top 100 분석 
print("=== Top 100 GO Enrichment ===")
# enricher로 분석 끝난 결과를 top_result에 저장 
# gene_name column만 가져오고 / python list로 바꿔줌
top_result = gp.enrichr(
    gene_list=top100["gene_name"].tolist(),
    gene_sets=gene_sets,
    organism="human",
    outdir=None #따로 파일로 저장 X 
).results

# Bottom 100
print("=== Bottom 100 GO Enrichment ===")
bottom_result = gp.enrichr(
    gene_list=bottom100["gene_name"].tolist(),
    gene_sets=gene_sets,
    organism="human",
    outdir=None
).results

# 유의한 결과만 선택: p-value < 0.05
print("=== Filtering significant GO terms based on P-val ===")
top_go = top_result[
    top_result["Adjusted P-value"] < 0.05 #이따 0.05로 다시 돌리기
].copy()

bottom_go = bottom_result[
    bottom_result["Adjusted P-value"] < 0.05
].copy()

# BP / MF / CC 이름 간단하게 변경
go_type = {
    "GO_Biological_Process_2025": "BP",
    "GO_Molecular_Function_2025": "MF",
    "GO_Cellular_Component_2025": "CC",
}
#바뀐 이름으로 저장
top_go["Type"] = top_go["Gene_set"].replace(go_type)
bottom_go["Type"] = bottom_go["Gene_set"].replace(go_type)

# 필요한 column만 출력
top_go = top_go[
    ["Type", "Term", "Adjusted P-value"]
]

bottom_go = bottom_go[
    ["Type", "Term", "Adjusted P-value"]
]

print("=== GO Enrichment Completed ===")
print("\n=== Top 100 GO ===")
print(top_go.to_string(index=False))

print("\n=== Bottom 100 GO ===")
print(bottom_go.to_string(index=False))