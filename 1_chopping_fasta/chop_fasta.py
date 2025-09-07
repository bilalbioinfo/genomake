#!/usr/bin/env python3


import sys
from Bio import SeqIO

def fasta_to_fastq(fasta_file, read_len, step, output_fastq):
    with open(output_fastq, 'w') as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = str(record.seq)
            for i in range(0, len(seq) - read_len + 1, step):
                read_seq = seq[i:i + read_len]
                qual = "I" * read_len  # Dummy Q40
                read_id = f"{record.id}_{i+1}-{i+read_len}"
                out.write(f"@{read_id}\n{read_seq}\n+\n{qual}\n")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python chop_fasta.py <input.fasta> <read_len> <step> <output.fq>")
        sys.exit(1)
    fasta_file = sys.argv[1]
    read_len = int(sys.argv[2])
    step = int(sys.argv[3])
    output_fastq = sys.argv[4]

    fasta_to_fastq(fasta_file, read_len, step, output_fastq)
