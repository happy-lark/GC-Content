#FASTA 파일 읽어서 dictionary로 저장
#chromosome/sequence ID ->key
# DNA sequence -> value
# SeqIO: 생물학적 sequence 파일 읽어주는 기능 
# len(genome): dictionary에 들어있는 key개수

from Bio import SeqIO

def read_fasta(fasta_path): #fasta_path: fasta 파일 경로
    print("Reading FASTA file...") #progress tracking
    genome = {} #dictionary 
    with open(fasta_path, "r") as file: 
        for record in SeqIO.parse(file, "fasta"): #SeqIO가 파일을 fasta 형식으로 해석 > fasta안에 잇는 sequence를 하나씩 가져옴
            genome[record.id] = str(record.seq).upper() #chromosome을 dictionary에 저장 
    print(f"FASTA loaded: {len(genome)} sequences") #fasta 다 읽고 몇개의 sequence가 저장됏는지 출력

    return genome
