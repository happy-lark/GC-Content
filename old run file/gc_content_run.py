#전체 GC Content와 region별 GC Content 계산 및 출력 (head로 5개까지)
import pandas as pd
#모든 column을 생략 없이 출력
pd.set_options("display.max_columns",None)
from main import gtf, genome
from gc_content import (
    calculate_gene_gc,
    calculate_region_gc)

print("\n GC Content calculation started!!")

#전체 protein-coding gene의 gc
print("===Total GC Content===")
total_gc_df = calculate_gene_gc(gtf, genome)
print(total_gc_df.head()) #확인용으로 앞의 5행만 출력햇음

#region별 gc
print("\n=== Region GC Content ===")
region_gc_df = calculate_region_gc(gtf, genome)
print(region_gc_df.head()) 

print("\nGC Content calculation completed!")
