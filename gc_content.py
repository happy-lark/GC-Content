# GC Content 구하기 (total + region별로)
# region별 gc content에 progress tracking 코드 잇으니까 최종본에서는 제거하기

import pandas as pd

from feature_annotation import (
    get_regions,
    union_intervals,
    calculate_introns,
    get_utr_regions,
    calculate_intergenic
)

# dna 문자열의 gc 계산
def gc_content(sequence, include_n):
    sequence = sequence.upper()
    gc_count = sequence.count("G") + sequence.count("C")

    if include_n:
        total = len(sequence)
    else:
        total = sum(sequence.count(base) for base in "ATGC")

    return gc_count / total * 100 if total > 0 else None


# feature 좌표로 fasta에서 dna 추출 (좌표 -> 실제 dna sequence)
def get_sequence(genome, chromosome, regions):
    sequence = ""

    for start, end in regions:
        sequence += genome[chromosome][start - 1:end]

    return sequence


# 위 두 함수를 연결해서 feature gc 계산
def feature_gc(genome, chromosome, regions):
    sequence = get_sequence(genome, chromosome, regions)

    return {
        "GC_with_N": gc_content(sequence, True),
        "GC_without_N": gc_content(sequence, False)
    }

#=====GC Cotent 계산
#1. chromosome별 GC Content 계산
# 일반 GC / noN GC로 구분
def calculate_chromosome_gc(genome):
    results = []
    
    chromosomes = [
        f"chr{i}" for i in range(1, 23)
    ] + ["chrX", "chrY"]

    for chromosome in chromosomes:
        sequence = genome[chromosome]
        results.append([
            chromosome,
            len(sequence),
            gc_content(sequence, True),
            gc_content(sequence, False)
        ])

    return pd.DataFrame(
        results,
        columns=[
            "Chromosome",
            "Length",
            "GC_with_N",
            "GC_without_N"
        ]
    )

# 2. gene level: 각 protein-coding gene 전체 영역의 GC
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
        columns=[
            "gene_id",
            "gene_name",
            "chromosome",
            "GC_content"
        ]
    )

#3. feature-level
# CDS / exon / intron / 5'UTR / 3'UTR / upstream / downstream
def calculate_region_gc(gtf, genome):
    results = []
    # gene annotation만 사용해서 chromosome별 intergenic 계산
    genes = gtf[gtf["feature"] == "gene"]
    intergenic = {}
    for chromosome, chromosome_genes in genes.groupby("chromosome"):
        for gene_id, upstream, downstream in calculate_intergenic(chromosome_genes):
            intergenic[gene_id] = (upstream, downstream)

    grouped = gtf.groupby("gene_id")
    total = gtf["gene_id"].nunique()

    for i, (gene_id, gene_df) in enumerate(grouped, start=1):

        chromosome = gene_df["chromosome"].iloc[0]
        gene_name = gene_df["gene_name"].iloc[0]

        _, cds, exon = get_regions(gene_df)

        cds = union_intervals(cds)
        exon = union_intervals(exon)

        intron = calculate_introns(exon)
        utr5, utr3 = get_utr_regions(gene_df)

        upstream, downstream = intergenic.get(gene_id, (None, None))

        results.append([
            gene_id,
            gene_name,
            chromosome,
            feature_gc(genome, chromosome, cds)["GC_without_N"],
            feature_gc(genome, chromosome, exon)["GC_without_N"],
            feature_gc(genome, chromosome, intron)["GC_without_N"],
            feature_gc(genome, chromosome, utr5)["GC_without_N"],
            feature_gc(genome, chromosome, utr3)["GC_without_N"],
            feature_gc(genome, chromosome, [upstream] if upstream else [])["GC_without_N"],
            feature_gc(genome, chromosome, [downstream] if downstream else [])["GC_without_N"]
        ])

        # progress tracking
        if i % 100 == 0 or i == total:
            print(f"Region GC: {i}/{total} genes completed")

    return pd.DataFrame(
        results,
        columns=[
            "gene_id",
            "gene_name",
            "chromosome",
            "CDS_GC",
            "Exon_GC",
            "Intron_GC",
            "5UTR_GC",
            "3UTR_GC",
            "Upstream_GC",
            "Downstream_GC"
        ]
    )

