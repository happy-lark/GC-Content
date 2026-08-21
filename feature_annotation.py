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

    #입력 좌표가 없으면 빈 리스트 반환 
    if not regions:
        return []

    regions = sorted(regions) 

    #첫번째 구간을 병합 기준으로 설정 
    merged = [regions[0]] 

    #두번째 구간부터 하나씩 이전 구간과 비교 
    for current in regions[1:]:
        current_start, current_end = current 
        previous = merged [-1] 

        if current_start <= previous[1] + 1:
            previous[1]=max(previous[1],current_end)
        else:
            merged.append(current)

    return merged

#3. 병합된 exon 사이의 빈 구간을 intron 좌표로 계산 후 리스트로 반환
def calculate_introns(exons):
    introns=[]
    
    for i in range(len(exons)-1):
        intron_start=exons[i][1]+1
        intron_end=exons[i+1][0]-1

        if intron_start<=intron_end:
            introns.append([intron_start,intron_end])

    return introns

# 4. transcript에서 CDS가 아닌 양쪽 영역을 UTR 후보로 계산하고
#   strand에 따라 5'UTR / 3' UTR로 구분한 후 병합해서 좌표를 리스트로 반환 
def get_utr_regions(gene_df):

    utr5 = []
    utr3 = []

    transcripts = gene_df[gene_df["feature"] == "transcript"]

    #iterrows: transcripts DataFrame 한행씩 반복
    #gene_df에서 CDS인 행만 선택 && 현재 반복중인 transcript와transcript_id 동일한 행만 선택
    #현재 transcript에 속하는 cds 행들만 gene_df에서 추출해서 cds에 저장
    for _, transcript in transcripts.iterrows():
        cds = gene_df[
            (gene_df["feature"] == "CDS")
            & (gene_df["transcript_id"] == transcript["transcript_id"])
        ]

        #skip transcripts without CDS
        if cds.empty:
            continue

        # 현재 transcript의 CDS 범위
        cds_start = cds["start"].min()
        cds_end = cds["end"].max()

        #left, right 범위
        left = [transcript["start"], cds_start - 1]
        right = [cds_end + 1, transcript["end"]]
        #left[0]:start // left[1]:end
        if transcript["strand"] == "+":
            if left[0] <= left[1]: #실제로 존재하는 유효한 좌표일때만 저장
                utr5.append(left)
            if right[0] <= right[1]:
                utr3.append(right)
        #transcript strand -  
        else:
            if left[0] <= left[1]:
                utr3.append(left)
            if right[0] <= right[1]:
                utr5.append(right)
    #5'UTR, 3'UTR 좌표 병합하고 결과 반환     
    return union_intervals(utr5), union_intervals(utr3)

# 5. gene과 gene 사이 intergenic 구간 계산
def calculate_intergenic(genes):

    results = []

    genes = genes.sort_values("start").reset_index(drop=True)

    for i in range(len(genes)):

        current_gene=genes.iloc[i]

        left_region=None
        right_region=None

        #이전 gene이 있을 때
        if i>0:
            previous_gene=genes.iloc[i-1]

            left_start=previous_gene["end"]+1
            left_end=current_gene["start"]-1

            left_region=[left_start, left_end]

        #담 gene이 있을 때
        if i< len(genes)-1: #왜 
            next_gene=genes.iloc[i+1]

            right_start=current_gene["end"]+1
            right_end=next_gene["start"]-1

            right_region=[right_start,right_end]

        #strand에 따라 upstream/downstream 결정
        if current_gene["strand"]=="+":
            upstream=left_region
            downstream=right_region
        else:
            upstream=right_region
            downstream=left_region

        results.append([
            current_gene["gene_id"], upstream, downstream])

    return results
