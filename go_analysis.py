import gseapy as gp

def enrich_go(gene_df):

    gene_sets = [
        "GO_Biological_Process_2025",
        "GO_Molecular_Function_2025",
        "GO_Cellular_Component_2025"
    ]

    result = gp.enrichr(
        gene_list=gene_df["gene_name"].dropna().tolist(), #NaN (결측값) 제거
        gene_sets=gene_sets,
        outdir=None #enrichr 결과를 별도 파일로 저장 X 
    ).results

    return result[
        result["Adjusted P-value"] < 0.05
    ]