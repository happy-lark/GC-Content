import matplotlib.pyplot as plt
import numpy as np

# feature별 GC content 분포: boxplot
def plot_feature_gc(feature_gc_df, features, feature_names):
    plt.boxplot(
        [feature_gc_df[feature].dropna() for feature in features],
        tick_labels=[feature_names[feature] for feature in features]
    )

    plt.ylabel("GC Content (%)")
    plt.title("Feature-level GC Content")
    plt.xticks(rotation=30)
    plt.tight_layout()
    plt.show()

# 전체 gene GC Content 분포: Histogram
def plot_gene_gc(gene_gc_df):

    plt.hist(gene_gc_df["GC_content"])

    plt.xlabel("GC Content (%)")
    plt.ylabel("Number of Genes")
    plt.title("Gene-level GC Content")
    plt.tight_layout()
    plt.show()

# feature별 median GC Content
def plot_feature_median(median_gc):

    plt.bar(
        median_gc.index,
        median_gc.values
    )

    plt.ylabel("Median GC Content (%)")
    plt.title("Feature-level GC Content Comparison")
    plt.xticks(rotation=30)
    plt.tight_layout()
    plt.show()

# GO Enrichment 결과: horizontal bar plot
def plot_go(go_df, title):

    significance = -np.log10(
        go_df["Adjusted P-value"]
    )

    plt.barh(
        go_df["Term"],
        significance
    )

    plt.xlabel("-log10(Adjusted P-value)")
    plt.title(title)
    plt.tight_layout()
    plt.show()