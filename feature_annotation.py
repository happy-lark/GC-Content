#GTF 좌표를 이용해서 CDS / exon / intron / UTR / intergenic 구간 정의 

# 1. gene, CDS, exon 좌표 추출
def get_regions(gene_df):

    gene_regions = (
        gene_df[gene_df["feature"] == "gene"]
        [["start", "end"]].values.tolist()
    )

    cds_regions = (
        gene_df[gene_df["feature"] == "CDS"]
        [["start", "end"]].values.tolist()
    )

    exon_regions = (
        gene_df[gene_df["feature"] == "exon"]
        [["start", "end"]].values.tolist()
    )

    return gene_regions, cds_regions, exon_regions

# 2. 겹쳐지거나 이어진 구간 CDS, exon 병합
def union_intervals(regions):

    #좌표가 하나도 없을때 
    if len(regions) == 0:
        return []

    #좌표를 start 기준으로 정렬 (작은 순서대로)
    regions = sorted(regions) 

    merged = [regions[0]]

    for start, end in regions[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:
            merged[-1][1] = max(last_end, end)
        else:
            merged.append([start, end])

    return merged

#3. exon union 사이를 intron으로 계산
def calculate_introns(exons):
    introns=[]
    #exon 2개 비교
    for i in range(len(exons)-1):
        intron_start=exons[i][1]+1
        intron_end=exons[i+1][0]-1

        if intron_start<=intron_end:
            introns.append([intron_start,intron_end])

    return introns

# 4. transcript에서 CDS가 아닌 부분을 찾아서 5'UTR / 3'UTR 계산 
def get_utr_regions(gene_df):

    utr5 = []
    utr3 = []

    transcripts = gene_df[gene_df["feature"] == "transcript"]

    for _, t in transcripts.iterrows():
        cds = gene_df[
            (gene_df["feature"] == "CDS")
            & (gene_df["transcript_id"] == t["transcript_id"])
        ]

        if cds.empty:
            continue

        cds_start = cds["start"].min()
        cds_end = cds["end"].max()

        left = [t["start"], cds_start - 1]
        right = [cds_end + 1, t["end"]]

        if t["strand"] == "+":
            if left[0] <= left[1]:
                utr5.append(left)
            if right[0] <= right[1]:
                utr3.append(right)

        else:
            if left[0] <= left[1]:
                utr3.append(left)
            if right[0] <= right[1]:
                utr5.append(right)

    return union_intervals(utr5), union_intervals(utr3)

# 5. gene과 gene 사이 intergenic 구간 계산
#strand에 따라 upstream/downstream 구분
def calculate_intergenic(genes):

    results = []

    genes = genes.sort_values("start").reset_index(drop=True)

    for i in range(len(genes)):

        gene = genes.iloc[i]

        left = None
        right = None

        if i > 0:
            left = [genes.iloc[i - 1]["end"] + 1, gene["start"] - 1]

        if i < len(genes) - 1:
            right = [gene["end"] + 1, genes.iloc[i + 1]["start"] - 1]

        if gene["strand"] == "+":
            upstream, downstream = left, right
        else:
            upstream, downstream = right, left

        results.append([
            gene["gene_id"],
            upstream,
            downstream
        ])

    return results