#GC Content 기준으로 상하위 100개

from main import gtf, genome 
from gc_content import calculate_gene_gc 

gene_gc_df=calculate_gene_gc(gtf,genome)
bottom100=gene_gc_df.nsmallest(100,"GC_content")
top100=gene_gc_df.nlargest(100,"GC_content")

print("\n Top/Bottom 100 gene selection start")

print("\n== Bottom 100 Genes==")
print(bottom100)
print("\n== Top 100 Genes==")
print(top100)
print("\nC Completed")