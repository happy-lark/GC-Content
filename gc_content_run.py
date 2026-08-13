from main import gtf, genome
from gc_content import feature_gc, calculate_gene_gc
from feature_annotation import (
    get_regions,
    union_intervals,
    calculate_introns,
    get_utr_regions,
)

print("\nGC Content calculation started...")

"""
[ Test용 : TP53 (범위 제한)] 
gene_df = gtf[gtf["gene_name"] == "TP53"]

# chromosome 확인
chromosome = gene_df["chromosome"].iloc[0]

# feature 좌표 구하기
gene_regions, cds_regions, exon_regions = get_regions(gene_df)

cds_union = union_intervals(cds_regions)
exon_union = union_intervals(exon_regions)

introns = calculate_introns(exon_union)

utr5, utr3 = get_utr_regions(gene_df)

# GC content 계산
print("\n=== TP53 GC Content ===")

print("Gene:", feature_gc(genome, chromosome, gene_regions))
print("CDS:", feature_gc(genome, chromosome, cds_union))
print("Exon:", feature_gc(genome, chromosome, exon_union))
print("Intron:", feature_gc(genome, chromosome, introns))
print("5'UTR:", feature_gc(genome, chromosome, utr5))
print("3'UTR:", feature_gc(genome, chromosome, utr3))

print("\nGC Content calculation completed!")
""" 
#전체 protein-coding gene 
gene_gc_df = calculate_gene_gc(gtf, genome)

print("\n=== Gene GC Content ===")
print(gene_gc_df.head()) #확인용으로 앞의 5행만 출력햇음

print("\nGC Content calculation completed!")
