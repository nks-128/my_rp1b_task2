#!/usr/bin/env python3

# This script compares called variants (VCF) against a known truth set
# generated during simulation. It reports precision and recall for:
#   - SNPs (exact match on chrom, pos, ref, alt)
#   - INDELs (match on chrom with a ±2 bp position window)
# Inputs:
#   truth.tsv      – list of simulated mutations (ground truth)
#   calls.vcf.gz   – VCF produced by bcftools/snippy
# Outputs:
#   Printed summary of TP, FP, FN, precision, recall
# The purpose of this script is to evaluate the accuracy of variant
# callers using the known simulated dataset.

import sys
import gzip

# Normalize chromosome names
def clean_chrom(name):
    return name.split()[0]


# Load truth mutations from truth.tsv
# SNPs stored as (chrom, pos, ref, alt)
# INDELs stored as (chrom, pos) only
def load_truth(tsv_path):
    truth_snp = set()
    truth_indel = []

    with open(tsv_path) as f:
        next(f)  # skip header

        for line in f:
            chrom, pos, ref, alt, t = line.strip().split("\t")
            chrom = clean_chrom(chrom)
            pos = int(pos)

            if t == "SNP":
                truth_snp.add((chrom, pos, ref, alt))
            else:
                # For INDEL validation only the position is checked
                truth_indel.append((chrom, pos))

    return truth_snp, truth_indel

# Load called variants from a VCF or VCF.gz
# SNPs stored as exact tuples
# INDELs stored by chrom and position only
def load_calls(vcf_path):
    snp = set()
    indel = []

    def open_vcf(x):
        return gzip.open(x, "rt") if x.endswith(".gz") else open(x)

    with open_vcf(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.strip().split("\t")
            chrom = clean_chrom(cols[0])
            pos = int(cols[1])
            ref = cols[3]
            alts = cols[4].split(",")

            # Handle multi-allelic sites
            for alt in alts:
                if len(ref) == 1 and len(alt) == 1:
                    snp.add((chrom, pos, ref, alt))
                else:
                    indel.append((chrom, pos))

    return snp, indel


# SNP evaluation: exact match on chrom/pos/ref/alt
def eval_snp(truth, calls):
    tp = len(truth & calls)
    fp = len(calls - truth)
    fn = len(truth - calls)

    prec = tp / (tp + fp) if (tp + fp) else 0
    rec  = tp / (tp + fn) if (tp + fn) else 0

    return tp, fp, fn, prec, rec

# INDEL evaluation:
# A match counts if chrom is the same and |pos_true - pos_call| ≤ 2
def eval_indel(truth, calls):
    tp = 0
    fn = 0
    fp = 0

    # True positives and false negatives
    for c, p in truth:
        match = any((cc == c) and (abs(pp - p) <= 2) for cc, pp in calls)
        if match:
            tp += 1
        else:
            fn += 1

    # False positives
    for cc, pp in calls:
        match = any((cc == c) and (abs(pp - p) <= 2) for c, p in truth)
        if not match:
            fp += 1

    prec = tp / (tp + fp) if (tp + fp) else 0
    rec  = tp / (tp + fn) if (tp + fn) else 0

    return tp, fp, fn, prec, rec

# Run evaluator from command line
def main():
    if len(sys.argv) != 3:
        print("Usage: evaluate_variants.py truth.tsv calls.vcf.gz")
        sys.exit(1)

    truth_file = sys.argv[1]
    vcf_file   = sys.argv[2]

    truth_snp, truth_indel = load_truth(truth_file)
    call_snp,  call_indel  = load_calls(vcf_file)

    # SNP results
    tp, fp, fn, prec, rec = eval_snp(truth_snp, call_snp)
    print("=== SNPs ===")
    print(f"TP={tp} FP={fp} FN={fn}")
    print(f"Precision={prec:.3f} Recall={rec:.3f}\n")

    # INDEL results
    tp, fp, fn, prec, rec = eval_indel(truth_indel, call_indel)
    print("=== INDELs ===")
    print(f"TP={tp} FP={fp} FN={fn}")
    print(f"Precision={prec:.3f} Recall={rec:.3f}")
    
# Run main() only when script is called from the command line
if __name__ == "__main__":
    main()
