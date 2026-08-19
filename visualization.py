import matplotlib.pyplot as plt
import numpy as np

# chromosome별 GC Content 
# bar graph
def plot_chromosome_gc(chromosome_gc_df):

    x = np.arange(len(chromosome_gc_df)) #chr 개수만큼 x축 위치 
    width = 0.4

    #GC with N
    plt.bar(
        x - width / 2, #bar끼리 겹치지 않도록 조금 왼쪽으로 이동 
        chromosome_gc_df["GC_with_N"], #y축 (GC Content 값)
        width,
        label="GC with N"
    )

    #GC without N
    plt.bar(
        x + width / 2,
        chromosome_gc_df["GC_without_N"], #조금 오른쪽으로 이동
        width,
        label="GC without N"
    )

    #x축 눈금 설정 (실제 chr 이름)
    plt.xticks(
        x,
        chromosome_gc_df["Chromosome"],
        rotation=45
    )

    plt.ylabel("GC Content (%)")
    plt.title("Chromosome-level GC Content")
    plt.legend() # 범례/label 표시 
    plt.tight_layout() # 여백 자동 조정 
    plt.show() # 그래프 출력 

# feature별 GC content 분포: boxplot
def plot_feature_gc(feature_gc_df, features, feature_names):

    plt.boxplot(
        [feature_gc_df[feature].dropna() for feature in features],
        tick_labels=[feature_names[feature] for feature in features] #x축 이름 설정
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


# feature별 median GC Content: bar plot
def plot_feature_median(median_gc):

    plt.bar(
        median_gc.index, #x축에 표시할 feature 이름 
        median_gc.values #각 feature의 median gc content
    )

    plt.ylabel("Median GC Content (%)")
    plt.title("Feature-level GC Content Comparison")
    plt.xticks(rotation=30)
    plt.tight_layout()
    plt.show()


# 상위 10개 gene 선택 / gene별 gc content를 horizontal bar plot으로 
def plot_gene_group(gene_df, column, title):

    genes = gene_df.head(10)

    plt.barh(
        genes["gene_name"],
        genes[column]
    )

    plt.xlabel("GC Content (%)")
    plt.title(title)
    plt.tight_layout()
    plt.show()

# GO Enrichment 결과: horizontal bar plot
def plot_go(go_df, title):
    plt.barh(go_df["Term"], -np.log10(go_df["Adjusted P-value"]))
    plt.xlabel("-log10(Adjusted P-value)")
    plt.title(title)
    plt.tight_layout()
    plt.show()