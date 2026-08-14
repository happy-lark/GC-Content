# 실행용 파일
from fasta import read_fasta

fasta_path = (
    r"C:\Users\User\Desktop\biglab\gc_project\data\genome\GRCh38.p14.genome.fa"
)
genome = read_fasta(fasta_path)

print("\nSequence IDs:")
print(genome.keys())

print("\nchr1 length:")
print(len(genome["chr1"]))

