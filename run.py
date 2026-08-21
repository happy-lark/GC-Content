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

    def calculate_gc(self):
        print("\nCalculating Chromosome level GC Content")
        self.chromosome_gc = calculate_chromosome_gc(self.genome)
        print("Calculating Gene level GC Content")
        self.gene_gc = calculate_gene_gc(self.gtf,self.genome)
        print("Calculating Feature level GC Content")
        self.feature_gc = calculate_region_gc(self.gtf,self.genomed)
        print("GC Content calculation completed")

    def select_feature(self):
        print("1. CDS")
        print("2. Exon")
        print("3. Intron")
        print("4. 5' UTR")
        print("5. 3' UTR")
        print("6. Upstream Intergenic")
        print("7. Downstream Intergenic")

        choice=input("Select Feature: ")

        if choice == "1":
            return "CDS_GC"
        elif choice == "2":
            return "Exon_GC"
        elif choice == "3":
            return "Intron_GC"
        elif choice == "4":
            return "5UTR_GC"
        elif choice == "5":
            return "3UTR_GC"
        elif choice == "6":
            return "Upstream_GC"
        elif choice == "7":
            return "Downstream_GC"
        else:
            print("Invalid feature")
            return None

    def get_feature_median(self):
        #feature GC column만 선택 
        feature_data=self.feature_gc[self.features]
        #각 feature의 median 계산
        median_gc=feature_data.median()
        #column이름 변경 (eg. CDS_GC -> CDS)
        median_gc.index=self.feature_names.values()
        return median_gc

    def show_feature_median(self):
        print("\n=== Feature level Median GC Content ===")
        print(self.get_feature_median().to_string()) #문자열로

    def get_top_bottom(self, feature):
        if feature == "Gene":
            data = self.gene_gc
            gc_column = "GC_content"
        else:
            data = self.feature_gc
            gc_column = feature

        top = data.nlargest(100, gc_column)
        bottom = data.nsmallest(100, gc_column)

        return top, bottom

    # Top / Bottom 출력
    def show_top_bottom(self, feature):
        top, bottom = self.get_top_bottom(feature)

        if feature=="Gene":
            name="Gene"
        else:
            name=self.feature_names[feature]

        print(f"\n=== {name} Top 100 ===")
        print(top.to_string(index=False))

        print(f"\n=== {name} Bottom 100 ===")
        print(bottom.to_string(index=False)) 

    # GO 분석
    def show_go(self):
        level = input(
            "1. Gene level\n"
            "2. Feature level\n"
            "Select: "
        )

        #Gene/Feature 선택 
        if level == "1":
            feature = "Gene"
        elif level == "2":
            feature = self.select_feature()
            if feature is None:
                return
        else:
            print("Invalid option")
            return

        # Top / Bottom 선택
        group=input("Top of Bottom: ").lower() 

        top, bottom =self.get_top_bottom(feature)

        if group=="top":
            genes=top
        elif group=="bottom":
            genes=bottom
        else:
            print("Invalid group")
            return

        print("\nGO enrichment started.")

        result = enrich_go(genes) 

        # 유의한 GO term이 없는 경우
        if result.empty:
            print("No significant GO terms after Adjusted p-value filtering")
            return
      
        print(result.to_string(index=False))

        # 그래프 title용 이름
        if feature=="Gene":
            name="Gene"
        else:
            name=self.feature_names[feature]

        plot_go(result,
                f"{name} {group.title()} 100 GO")

    # Visualization
    def show_visualization(self):

        print("""
1. Chromosome level GC
2. Gene level GC 
3. Feature level GC 
4. Feature level Median GC 
5. Gene level Top/Bottom
6. Feature level Top/Bottom
""")
        graph = input("Select: ") 

        if graph == "1":
            plot_chromosome_gc(self.chromosome_gc)

        elif graph == "2":
            plot_gene_gc(self.gene_gc)

        elif graph == "3":
            plot_feature_gc(self.feature_gc, self.features, self.feature_names) 

        elif graph == "4":
            plot_feature_median(self.get_feature_median()) 

        elif graph == "5":
            top, bottom = (self.get_top_bottom("Gene")) 
            plot_gene_group(
                top, #그래프에 사용할 데이터
                "GC_content", #그래프에 사용할 gc column
                "Gene Top 100" #그래프 제목
                ) 
            plot_gene_group(bottom,"GC_content","Gene Bottom 100")

        elif graph == "6":
            feature = (self.select_feature())

            if feature is None:
                return

            top, bottom = (self.get_top_bottom(feature))

            if feature=="Gene":
                name="Gene"

            else:
                name=self.feature_names[feature]

            plot_gene_group(top,feature,f"{name} Top 100")
            plot_gene_group(bottom,feature,f"{name} Bottom 100")

        else:
            print("Invalid Option")

    # 메인 메뉴
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
                display_df=(self.feature_gc.rename(columns=self.feature_names))
                print(display_df.to_string(inedex=False))

            elif choice == "4":
                self.show_feature_median()

            elif choice == "5":
                self.show_top_bottom("Gene")

            elif choice == "6":
                feature=(self.select_feature()) 
                if feature:
                    self.show_top_bottom(feature) 

            elif choice == "7":
                self.show_go()

            elif choice == "8":
                self.show_visualization()

            elif choice == "0":
                print("Analysis completed")
                break 
            
            else: 
                print("Invlalid Option. Try again")

analysis = GCAnalysis(gtf, genome) 
analysis.run() 

