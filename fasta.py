#FASTA 파일 읽어서 dictionary로 저장

from Bio import SeqIO

def read_fasta(fasta_path):

    print("Reading FASTA file...") #progress tracking

    genome = {} #dictionary로 저장

    with open(fasta_path, "r") as file:
        for record in SeqIO.parse(file, "fasta"): #file을 fasta형식으로 해석 
            genome[record.id] = str(record.seq).upper() #chromosome을 dictionary에 저장

    print(f"FASTA loaded: {len(genome)} sequences")

    return genome

