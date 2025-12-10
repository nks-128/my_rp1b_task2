#!/usr/bin/env bash

# This script runs a complete variant-calling pipeline using:
#   - minimap2  (read mapping)
#   - samtools  (BAM sorting/indexing)
#   - bcftools (variant calling)
#   - snippy   (second variant caller)
# Inputs:
#   reference.fasta
#   R1.fastq
#   R2.fastq
#   prefix     (name for output files)
# Outputs:
#   prefix.bam
#   prefix.bcftools.vcf.gz
#   prefix.snippy.vcf.gz
# The script automatically builds the minimap2 index if missing.

set -euo pipefail   # Safe-mode: exit on errors, undefined vars, pipe failures

# Check correct number of args
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <ref.fasta> <R1.fastq> <R2.fastq> <prefix>" >&2
    exit 1
fi

# Assign input parameters
REF=$1
R1=$2
R2=$3
PREFIX=$4

echo "[INFO] Reference: $REF"
echo "[INFO] Reads: $R1 / $R2"
echo "[INFO] Prefix: $PREFIX"

# 1. Build minimap2 index (.mmi) if it does not already exist
MMI="${REF%.fasta}.mmi"

if [ ! -f "$MMI" ]; then
    echo "[INFO] Building minimap2 index: $MMI"
    minimap2 -d "$MMI" "$REF"
else
    echo "[INFO] Using existing index: $MMI"
fi

# 2. Map reads with minimap2 and sort BAM with samtools
echo "[INFO] Mapping reads with minimap2 and sorting BAM..."

minimap2 -ax sr "$MMI" "$R1" "$R2" | \
    samtools sort -o "${PREFIX}.bam"

samtools index "${PREFIX}.bam"

# 3. Variant calling with bcftools (mpileup + call)
echo "[INFO] Calling variants with bcftools..."

bcftools mpileup -Ou -f "$REF" "${PREFIX}.bam" | \
    bcftools call -mv -Oz -o "${PREFIX}.bcftools.vcf.gz"

bcftools index "${PREFIX}.bcftools.vcf.gz"

# 4. Variant calling with snippy
echo "[INFO] Running snippy..."

SNPDIR="${PREFIX}.snippy"

snippy \
    --outdir "$SNPDIR" \
    --ref "$REF" \
    --R1 "$R1" \
    --R2 "$R2" \
    --cpus 4

# Confirm snippy generated a VCF
if [ ! -f "${SNPDIR}/snps.vcf" ]; then
    echo "[ERROR] Snippy did not produce snps.vcf in $SNPDIR" >&2
    exit 1
fi

# 5. Compress and index snippy VCF
echo "[INFO] Compressing snippy VCF..."
bgzip -c "${SNPDIR}/snps.vcf" > "${PREFIX}.snippy.vcf.gz"
bcftools index "${PREFIX}.snippy.vcf.gz"

echo "[INFO] Done. Generated:"
echo "  - ${PREFIX}.bcftools.vcf.gz"
echo "  - ${PREFIX}.snippy.vcf.gz"

