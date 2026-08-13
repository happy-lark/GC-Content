#run 파일들 중복되는 부분,,,, 나중에 확인 끝나면 run 파일 앖애기
from annotation import load_annotation
from fasta import read_fasta

# file path
gtf_path = r"data\annotation\gencode.v50.annotation.gtf"
fasta_path = r"data\genome\GRCh38.p14.genome.fa"

# load GTF, FASTA
gtf = load_annotation(gtf_path)
genome = read_fasta(fasta_path)
print("\nData loading completed!")