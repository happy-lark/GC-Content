#GC Content 기준으로 total/region별 상하위 gene 100개 추출

from main import gtf, genome 
from gc_content import calculate_gene_gc, calculate_region_gc

# total
total_gc_df=calculate_gene_gc(gtf,genome)

total_bottom100=total_gc_df.nsmallest(100,"GC_content")
total_top100=total_gc_df.nlargest(100,"GC_content")

# region 
region_gc_df=calculate_region_gc(gtf,genome)

regions = [
    "CDS_GC",
    "Exon_GC",
    "Intron_GC",
    "5UTR_GC",
    "3UTR_GC"
]

region_top100 = {}
region_bottom100 = {}

for region in regions:
    region_top100[region] = region_gc_df.nlargest(100, region)
    region_bottom100[region] = region_gc_df.nsmallest(100, region)

#결과 출력 
print("\n ===Total Bottom 100 Genes===")
print(total_bottom100.to_string(index=False))

print("\n ===Total Top 100 Genes===")
print(total_top100.to_string(index=False))

for region in regions:
    print(f"\n=== {region} Bottom 100 Genes ===")
    print(region_bottom100[region].to_string(index=False))

    print(f"\n=== {region} Top 100 Genes ===")
    print(region_top100[region].to_string(index=False))

print("\nCompleted")