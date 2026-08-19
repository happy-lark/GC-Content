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

        #입력 받은 GTF, genome 데이터를 객체 내부 저장
        self.gtf = gtf
        self.genome = genome

        #dictionary
        self.feature_names = {
            "CDS_GC": "CDS",
            "Exon_GC": "Exon",
            "Intron_GC": "Intron",
            "5UTR_GC": "5' UTR",
            "3UTR_GC": "3' UTR",
            "Upstream_GC": "Upstream Intergenic",
            "Downstream_GC": "Downstream Intergenic"
        }
        #dictionary key만 리스트로 저장 
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

    # Feature 선택
    def select_feature(self):
        #self.features에 있는 feature들을 하나씩 반복하면서 번호 붙임 
        # enumerate(): 리스트 반복하면서 값/번호를 같이 반환
        for i, feature in enumerate(self.features, start=1):
            print(f"{i}. {self.feature_names[feature]}")

        #사용자에게 번호를 입력받아 choice에 저장
        choice = input("Select feature: ")

        if not choice.isdigit():
            print("Invalid feature")
            return None

        index = int(choice) - 1

        #입력한 번호가 실제 feature 안에 있는지 확인 
        if index not in range(len(self.features)):
            print("Invalid feature")
            return None
        
        #선택한 feature의 실제 column 이름 반환 
        return self.features[index]

    # Feature별 median GC content 계산
    def get_feature_median(self):
        median_gc = self.feature_gc[self.features].median()
        median_gc.index = self.feature_names.values()
        return median_gc

    # Feature별 median GC 출력
    def show_feature_median(self):
        print("\n=== Feature-level Median GC Content ===")
        print(self.get_feature_median().to_string())

    # Top / Bottom 100
    def get_top_bottom(self, feature):
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
    def show_top_bottom(self, feature):
        top, bottom = self.get_top_bottom(feature)
        name = self.feature_names.get(feature, "Gene")

        print(f"\n=== {name} Top 100 ===")
        print(top.to_string(index=False))

        print(f"\n=== {name} Bottom 100 ===")
        print(bottom.to_string(index=False))

    # GO 분석
    def show_go(self):
        level = input("1. Gene-level\n2. Feature-level\nSelect: ")

        if level == "1":
            feature = "Gene"
        elif level == "2":
            feature = self.select_feature()
            if feature is None:
                return
        else:
            print("Invalid option")
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
                f"{group.title()} 100: No significant GO terms after Adjusted P-value filtering")
            return

        print(result.to_string(index=False))

        name = self.feature_names.get(feature, "Gene")
        plot_go(result, f"{name} {group.title()} 100 GO")

    # Visualization
    def show_visualization(self):
        print("""
1. Chromosome-level GC
2. Gene-level GC
3. Feature-level GC
4. Feature-level Median GC
5. Gene-level Top/Bottom
6. Feature-level Top/Bottom
""")
        graph = input("Select: ")

        if graph == "1":
            plot_chromosome_gc(self.chromosome_gc)
        elif graph == "2":
            plot_gene_gc(self.gene_gc)
        elif graph == "3":
            plot_feature_gc(
                self.feature_gc,
                self.features,
                self.feature_names
            )
        elif graph == "4":
            plot_feature_median(self.get_feature_median())
        elif graph==["5","6"]:
            feature="Gene" if graph=="5" else self.select.feature() 
            if feature is None:
                return 
            top, bottom=self.get_top_bottom(feature)

            column="GC_content" if feature=="Gene" else feature
            name=self.feature_names.get(feature,"Gene")

            plot_gene_group(top,column,f"{name} Top 100")
            plot_gene_group(bottom,column, f"{name} Bottom 100")
        else:
            print("Invalid option")

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
                feature = self.select_feature()
                if feature:
                    self.show_top_bottom(feature)
            elif choice == "7":
                self.show_go()
            elif choice == "8":
                self.show_visualization()
            elif choice == "0":
                print("Analysis completed!")
                break
            else:
                print("Invalid option")
analysis = GCAnalysis(gtf, genome)
analysis.run()

