#GC Content 구하기
 
from fasta import read_fasta
from feature_annotation import get_regions
import pandas as pd

#dna 문자열의 gc 계산
def gc_content(sequence,include_n):
    sequence=sequence.upper()
    gc_count=sequence.count("G")+sequence.count("C")

    if include_n:
        total=len(sequence)
    else:
        total=sum(sequence.count(base) for base in "ATGC")
    return gc_count/total*100

# feature 좌표로 fasta에서 dna 추출 (좌표 -> 실제 dna sequence)
def get_sequence(genome,chromosome,regions):
    sequence=""
    for start,end in regions:
        sequence+=genome[chromosome][start-1:end]
        return sequence

#위 두 함수를 연결해서 feature gc 계산 
def feature_gc(genome, chromosome, regions):
    sequence = get_sequence(genome, chromosome, regions)
    return {
        "GC_with_N": gc_content(sequence, True),
        "GC_without_N": gc_content(sequence, False)
    }

#gtf: gene의 위치정보가 들어있는 dataframe
#genome: chromosome별 실제 dna sequence가 들어 잇는 dictionary

def calculate_gene_gc(gtf, genome):
    results = []
    for gene_id, gene_df in gtf.groupby("gene_id"):
        chromosome = gene_df["chromosome"].iloc[0]
        gene_name = gene_df["gene_name"].iloc[0]
        gene_regions, _, _ = get_regions(gene_df)

        gc = feature_gc(genome, chromosome, gene_regions)

        results.append([
            gene_id,
            gene_name,
            chromosome,
            gc["GC_without_N"]
        ])

    return pd.DataFrame(
        results,
        columns=["gene_id", "gene_name", "chromosome", "GC_content"]
    )