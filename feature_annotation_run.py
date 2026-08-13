from annotation import load_annotation
from feature_annotation import (
    get_regions,
    union_intervals,
    calculate_introns,
    get_utr_regions,
    calculate_intergenic,
)
gtf_path=(r"C:\Users\User\Desktop\biglab\gc_project\data\annotation\gencode.v50.annotation.gtf")

# annotation 불러오기
gtf = load_annotation(gtf_path)

print("\nFeature annotation started...")

# TP53 선택
gene_df = gtf[gtf["gene_name"] == "TP53"]

# gene / CDS / exon 좌표
gene_regions, cds_regions, exon_regions = get_regions(gene_df)

# 겹치는 영역 합치기
cds_union = union_intervals(cds_regions)
exon_union = union_intervals(exon_regions)

# intron 계산
introns = calculate_introns(exon_union)

# UTR 계산
utr5, utr3 = get_utr_regions(gene_df)

# 결과 출력
print("\n=== TP53 ===")
print("Gene:", gene_regions)
print("CDS:", cds_union)
print("Exon:", exon_union)
print("Intron:", introns)
print("5'UTR:", utr5)
print("3'UTR:", utr3)

print("\nFeature annotation completed!")