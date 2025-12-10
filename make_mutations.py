#!/usr/bin/env python3

# This script introduces random mutations into a reference genome:
# - 300 SNPs (single base changes)
# - 20 INDELs (small insertions/deletions of 1–10 bp)
# It outputs:
#   1. A mutated genome FASTA
#   2. A truth table listing each mutation (for scoring the variant caller)

import random
import sys

# Fix the random seed so the results are reproducible every time
random.seed(42)

# Function: read_fasta
# Reads FASTA file and returns:
#   - header ("NC_037282.1 Plasmodium falciparum 3D7 genome assembly, chromosome: 11")
#   - full sequence as one long string
def read_fasta(path):
    header = None
    seq_parts = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # Save the header without the ">"
                header = line[1:]
            else:
                # Add sequence line to a growing list
                seq_parts.append(line)

    # Join all sequence lines into one continuous string
    return header, "".join(seq_parts)

# Function: write_fasta
# Writes a header + sequence to an output FASTA file
# Wraps sequence at 60 bases per line (standard FASTA format)
def write_fasta(header, seq, outpath):
    with open(outpath, "w") as out:
        out.write(f">{header}\n")
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60] + "\n")

# Function: random_base_not
# Returns a base that is not the original (e.g., A → C/T/G)
def random_base_not(b):
    bases = ["A","C","G","T"]
    if b in bases:
        bases.remove(b)   # remove the original base
    return random.choice(bases)

# Function: make_mutations
#   - randomly selects positions for SNPs (ensuring no duplicates)
#   - randomly selects positions + lengths for insertions/deletions
#   - applies these mutations to create a mutated genome sequence
def make_mutations(seq, num_snps=300, num_indels=20):
    L = len(seq)
    used = set()  # stores positions already mutated (avoid double mutation)

  
    # 1. Generate random SNPs
 
    snps = []
    while len(snps) < num_snps:
        pos = random.randint(1, L)

        # avoid mutating same position twice
        if pos in used:
            continue

        ref = seq[pos-1]
        if ref == "N":  # skip unknown bases
            continue

        alt = random_base_not(ref)
        used.add(pos)
        snps.append((pos, ref, alt))

    
    # 2. Generate random INDELs
   
    indels = []
    while len(indels) < num_indels:
        pos = random.randint(1, L)

        if pos in used:
            continue

        used.add(pos)
        indel_type = random.choice(["INS", "DEL"])
        length = random.randint(1, 10)  # size 1–10bp
        indels.append((pos, indel_type, length))

    
    # 3. Combine SNPs + INDELs to a single list
    
    mutations = []

    # Format SNPs into dictionaries
    for pos, ref, alt in snps:
        mutations.append({"type":"SNP", "pos":pos, "ref":ref, "alt":alt})

    # For INDELs:
    bases = ["A","C","G","T"]
    for pos, typ, length in indels:

        if typ == "INS":
            # Insert random bases
            alt = "".join(random.choice(bases) for _ in range(length))
            ref = "-"                # truth format uses "-" to represent insertion
        else:
            # Delete bases from the reference genome
            ref = seq[pos-1:pos-1+length]
            alt = "-"

        mutations.append({"type":typ, "pos":pos, "ref":ref, "alt":alt})

    
    # 4. Apply mutations to the original sequence
    #   Applied from right → left (reverse sorted by position),
    #   ensuring INDELs do not shift the positions of later mutations.
   
    mut_seq = seq

    for m in sorted(mutations, key=lambda x: x["pos"], reverse=True):
        i = m["pos"] - 1  # convert to 0-based index

        if m["type"] == "SNP":
            mut_seq = mut_seq[:i] + m["alt"] + mut_seq[i+1:]

        elif m["type"] == "INS":
            mut_seq = mut_seq[:i] + m["alt"] + mut_seq[i:]

        elif m["type"] == "DEL":
            mut_seq = mut_seq[:i] + mut_seq[i+len(m["ref"]):]

    return mut_seq, mutations


# Function: write_truth
# Writes a tab-separated truth table:
# chrom   pos   ref   alt   type
# (This file is used to compare with the VCF later) 
def write_truth(chrom, muts, out):
    with open(out, "w") as f:
        f.write("chrom\tpos\tref\talt\ttype\n")
        for m in sorted(muts, key=lambda x: x["pos"]):
            f.write(f"{chrom}\t{m['pos']}\t{m['ref']}\t{m['alt']}\t{m['type']}\n")


# This section runs the script from the command line.
# Usage:
#   make_mutations.py <input.fa> <mutated.fa> <truth.tsv>
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: make_mutations.py <input.fa> <output.fa> <truth.tsv>")
        sys.exit(1)

    # 1. Read the reference genome
    chrom, seq = read_fasta(sys.argv[1])

    # 2. Make mutations
    mutseq, muts = make_mutations(seq)

    # 3. Write mutated genome
    write_fasta(chrom, mutseq, sys.argv[2])

    # 4. Write truth file listing all mutations
    write_truth(chrom, muts, sys.argv[3])
