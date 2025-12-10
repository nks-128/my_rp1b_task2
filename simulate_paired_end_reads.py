#!/usr/bin/env python3

# This script simulates paired-end reads from the mutated genome.
# It generates:
#   - R1.fastq  (forward reads)
#   - R2.fastq  (reverse-complement reads)
# Inputs:
#   mutated.fasta   – mutated genome
#   coverage        – e.g. 30x
#   read_length     – e.g. 100 bp
# Reads are sampled randomly, with realistic insert sizes (100–300 bp),
# and have perfect quality for testing variant callers.

import sys
import random

random.seed(123) # reproducible simulation

def read_fasta(path):
    header = None
    seq_parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:]   
            else:
                # Ensure sequence is uppercase (A,C,G,T)
                seq_parts.append(line.upper())
    return header, "".join(seq_parts)

# Function: revcomp
# Returns the reverse complement of a DNA sequence.
# Used to generate the R2 read, which comes from the opposite strand.
def revcomp(seq):
    comp = str.maketrans("ACGT", "TGCA")
    return seq.translate(comp)[::-1]

# Function: simulate_paired_end
#   - calculates how many read pairs are needed to reach the target coverage
#   - randomly selects start positions in the genome
#   - extracts R1 (forward) and R2 (reverse-complement) reads
#   - writes the reads in FASTQ format

def simulate_paired_end(seq, coverage, read_len, out_R1, out_R2):

    genome_len = len(seq)

    # Number of read pairs needed:
    # coverage = (2 * read_len * num_pairs) / genome_len
    num_pairs = int((coverage * genome_len) / (2 * read_len))

    print(f"Genome length: {genome_len}")
    print(f"Coverage: {coverage}x")
    print(f"Read length: {read_len} bp (paired-end)")
    print(f"Number of read pairs: {num_pairs}")

    with open(out_R1, "w") as f1, open(out_R2, "w") as f2:

        for i in range(num_pairs):

            # Choose a random starting point to ensure there is enough
            # space for both reads and the insert size
            start = random.randint(0, genome_len - (2 * read_len))

            # R1 comes directly from the forward strand
            read1 = seq[start:start + read_len]

            # Insert size: distance between R1 start and R2 start
            insert_size = read_len + random.randint(100, 300)
            r2_start = start + insert_size

            # Ensure R2 is fully inside the genome
            if r2_start + read_len >= genome_len:
                continue

            # R2 is reverse-complemented
            read2 = revcomp(seq[r2_start:r2_start + read_len])

            # Perfect quality string (“I” = ASCII 73 ≈ Q40)
            qual = "I" * read_len

            # Write FASTQ entries for both reads
            f1.write(f"@read{i}/1\n{read1}\n+\n{qual}\n")
            f2.write(f"@read{i}/2\n{read2}\n+\n{qual}\n")

# This section reads command-line arguments and runs the
# simulator with the provided inputs.
if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: simulate_paired_end_reads.py <mutated.fasta> <coverage> <read_length> <R1.fastq> <R2.fastq>")
        sys.exit(1)

    fasta = sys.argv[1]
    coverage = float(sys.argv[2])
    read_len = int(sys.argv[3])
    out_R1 = sys.argv[4]
    out_R2 = sys.argv[5]

    header, seq = read_fasta(fasta)
    simulate_paired_end(seq, coverage, read_len, out_R1, out_R2)
