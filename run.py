#최종 실행 파일
from main import gtf, genome
from gc_content import calculate_gene_gc, calculate_region_gc
from go_analysis import enrich_go
from visualization import plot_total_gc, plot_region_gc, plot_go

class GCAnalysis:

    def __init__(self, gtf, genome):
        self.gtf = gtf
        self.genome = genome
        self.regions = [
            "CDS_GC",
            "Exon_GC",
            "Intron_GC",
            "5UTR_GC",
            "3UTR_GC",
            "Upstream_GC",
            "Downstream_GC"
        ]

    # GC Content 계산
    def calculate_gc(self):
        print("\nCalculating GC Content...")

        self.total_gc = calculate_gene_gc(
            self.gtf, self.genome
        )
        self.region_gc = calculate_region_gc(
            self.gtf, self.genome
        )

        print("GC Content completed!")


    # Top / Bottom 100
    def get_top_bottom(self, region="Total"):

        if region == "Total":
            df = self.total_gc
            column = "GC_content"

        else:
            df = self.region_gc
            column = region

        return (
            df.nlargest(100, column),
            df.nsmallest(100, column)
        )


    # Top / Bottom 출력
    def show_top_bottom(self, region="Total"):

        top, bottom = self.get_top_bottom(region)

        print(f"\n=== {region} Top 100 ===")
        print(top.to_string(index=False))

        print(f"\n=== {region} Bottom 100 ===")
        print(bottom.to_string(index=False))


    # GO 분석
    def show_go(self):

        region = input("Total or region: ")
        group = input("Top or Bottom: ").lower()

        if region != "Total" and region not in self.regions:
            print("Invalid region")
            return

        if group not in ["top", "bottom"]:
            print("Invalid group")
            return

        top, bottom = self.get_top_bottom(region)
        genes = top if group == "top" else bottom

        print("\nGO enrichment started...")
        result = enrich_go(genes)

        if result.empty:
            print("No significant GO terms")
            return

        print(result.to_string(index=False))

        plot_go(
            result,
            f"{region} {group.title()} 100 GO"
        )


    # 메뉴
    def run(self):

        self.calculate_gc()

        while True:

            print("""
1. Total GC Content
2. Region GC Content
3. Total Top/Bottom 100
4. Region Top/Bottom 100
5. GO Enrichment
6. Visualization
0. Exit
""")

            choice = input("Select option: ")


            if choice == "1":
                print(self.total_gc.to_string(index=False))


            elif choice == "2":
                print(self.region_gc.to_string(index=False))


            elif choice == "3":
                self.show_top_bottom()


            elif choice == "4":
                region = input(f"Select region {self.regions}: ")

                if region in self.regions:
                    self.show_top_bottom(region)
                else:
                    print("Invalid region")


            elif choice == "5":
                self.show_go()


            elif choice == "6":
                graph = input(
                    "1. Total GC\n"
                    "2. Region GC\n"
                    "Select: "
                )

                if graph == "1":
                    plot_total_gc(self.total_gc)

                elif graph == "2":
                    plot_region_gc(self.region_gc)


            elif choice == "0":
                break


            else:
                print("Invalid option")

analysis = GCAnalysis(gtf, genome)
analysis.run()