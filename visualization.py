import matplotlib.pyplot as plt
import numpy as np

#region별 GC content 분포: boxplot
def plot_region_gc(region_gc_df):
    regions = [
        "CDS_GC",
        "Exon_GC",
        "Intron_GC",
        "5UTR_GC",
        "3UTR_GC"
    ]
    plt.boxplot(
        [region_gc_df[region].dropna() for region in regions],
        tick_labels=regions
    )
    plt.ylabel("GC Content (%)")
    plt.title("GC Content by Region")
    plt.tight_layout()
    plt.show()

# 전체 gene GC Content 분포: Histogram
def plot_total_gc(total_gc_df):
    plt.hist(total_gc_df["GC_content"])
    plt.xlabel("GC Content (%)")
    plt.ylabel("Number of Genes")
    plt.title("Total Gene GC Content")
    plt.tight_layout()
    plt.show()

#GO Enrichment 결과: horizontal bar plot  
def plot_go(go_df, title):
    significance = -np.log10(go_df["Adjusted P-value"])
    plt.barh(
        go_df["Term"],
        significance
    )
    plt.xlabel("-log10(Adjusted P-value)")
    plt.title(title)
    plt.tight_layout()
    plt.show()