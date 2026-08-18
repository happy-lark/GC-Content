# 최종 실행 파일
from main import gtf, genome
from gc_content import (
    calculate_chromosome_gc,
    calculate_gene_gc,
    calculate_region_gc
)
from go_analysis import enrich_go
from visualization import (
    plot_chromosome_gc,
    plot_gene_gc,
    plot_feature_gc,
    plot_feature_median,
    plot_gene_group,
    plot_go
)

class GCAnalysis:

    def __init__(self, gtf, genome):
        self.gtf = gtf
        self.genome = genome

        self.feature_names = {
            "CDS_GC": "CDS",
            "Exon_GC": "Exon",
            "Intron_GC": "Intron",
            "5UTR_GC": "5' UTR",
            "3UTR_GC": "3' UTR",
            "Upstream_GC": "Upstream Intergenic",
            "Downstream_GC": "Downstream Intergenic"
        }

        self.features = list(self.feature_names)

    # GC Content 계산
    def calculate_gc(self):

        print("\nCalculating Chromosome-level GC Content...")
        self.chromosome_gc = calculate_chromosome_gc(self.genome)

        print("Calculating Gene-level GC Content...")
        self.gene_gc = calculate_gene_gc(self.gtf, self.genome)

        print("Calculating Feature-level GC Content...")
        self.feature_gc = calculate_region_gc(self.gtf, self.genome)

        print("GC Content calculation completed!")

    # Feature별 median GC 계산
    def get_feature_median(self):

        median_gc = self.feature_gc[self.features].median()
        median_gc.index = self.feature_names.values()

        return median_gc

    # Feature별 median GC 출력
    def show_feature_median(self):

        print("\n=== Feature-level Median GC Content ===")
        print(self.get_feature_median().to_string())

    # Top / Bottom 100
    def get_top_bottom(self, feature="Gene"):

        if feature == "Gene":
            df = self.gene_gc
            column = "GC_content"
        else:
            df = self.feature_gc
            column = feature

        return (
            df.nlargest(100, column),
            df.nsmallest(100, column)
        )

    # Top / Bottom 출력
    def show_top_bottom(self, feature="Gene"):

        top, bottom = self.get_top_bottom(feature)
        name = self.feature_names.get(feature, "Gene")

        print(f"\n=== {name} Top 100 ===")
        print(top.to_string(index=False))

        print(f"\n=== {name} Bottom 100 ===")
        print(bottom.to_string(index=False))

    # GO 분석
    def show_go(self):

        feature = input(f"Gene or feature {self.features}: ")

        if feature.lower() == "gene":
            feature = "Gene"

        if feature != "Gene" and feature not in self.features:
            print("Invalid feature")
            return

        group = input("Top or Bottom: ").lower()

        if group not in ["top", "bottom"]:
            print("Invalid group")
            return

        top, bottom = self.get_top_bottom(feature)
        genes = top if group == "top" else bottom

        print("\nGO enrichment started...")
        result = enrich_go(genes)

        if result.empty:
            print(
                f"{group.title()} 100: "
                "No significant GO terms after Adjusted P-value filtering"
            )
            return

        print(result.to_string(index=False))

        name = self.feature_names.get(feature, "Gene")
        plot_go(result, f"{name} {group.title()} 100 GO")

    # 메뉴
    def run(self):

        self.calculate_gc()

        while True:

            print("""
1. Chromosome-level GC Content
2. Gene-level GC Content
3. Feature-level GC Content
4. Feature-level Median GC Content
5. Gene-level Top/Bottom 100
6. Feature-level Top/Bottom 100
7. GO Enrichment
8. Visualization
0. Exit
""")

            choice = input("Select option: ")

            if choice == "1":
                print(self.chromosome_gc.to_string(index=False))

            elif choice == "2":
                print(self.gene_gc.to_string(index=False))

            elif choice == "3":
                display_df = self.feature_gc.rename(
                    columns=self.feature_names
                )
                print(display_df.to_string(index=False))

            elif choice == "4":
                self.show_feature_median()

            elif choice == "5":
                self.show_top_bottom()

            elif choice == "6":
                feature = input(f"Select feature {self.features}: ")

                if feature in self.features:
                    self.show_top_bottom(feature)
                else:
                    print("Invalid feature")

            elif choice == "7":
                self.show_go()

            elif choice == "8":

                graph = input(
                    "1. Chromosome-level GC\n"
                    "2. Gene-level GC\n"
                    "3. Feature-level GC\n"
                    "4. Feature-level Median GC\n"
                    "5. Gene-level Top/Bottom\n"
                    "6. Feature-level Top/Bottom\n"
                    "Select: "
                )

                if graph == "1":
                    plot_chromosome_gc(
                        self.chromosome_gc
                    )

                elif graph == "2":
                    plot_gene_gc(
                        self.gene_gc
                    )

                elif graph == "3":
                    plot_feature_gc(
                        self.feature_gc,
                        self.features,
                        self.feature_names
                    )

                elif graph == "4":
                    plot_feature_median(
                        self.get_feature_median()
                    )

                elif graph == "5":
                    top, bottom = self.get_top_bottom()

                    plot_gene_group(
                        top,
                        "GC_content",
                        "Gene-level Top 10"
                    )

                    plot_gene_group(
                        bottom,
                        "GC_content",
                        "Gene-level Bottom 10"
                    )

                elif graph == "6":
                    feature = input(
                        f"Select feature {self.features}: "
                    )

                    if feature in self.features:
                        top, bottom = self.get_top_bottom(feature)
                        name = self.feature_names[feature]

                        plot_gene_group(
                            top,
                            feature,
                            f"{name} Top 10"
                        )

                        plot_gene_group(
                            bottom,
                            feature,
                            f"{name} Bottom 10"
                        )

                    else:
                        print("Invalid feature")

                else:
                    print("Invalid option")

            elif choice == "0":
                print("Analysis completed!")
                break

            else:
                print("Invalid option")


analysis = GCAnalysis(gtf, genome)
analysis.run()