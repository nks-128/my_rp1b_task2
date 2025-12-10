#!/usr/bin/env python3

import sys
import gzip

# Open VCF files (whether they are .vcf or .vcf.gz)
def open_any(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")   # read gzipped VCF
    return open(path, "r")             # read normal VCF

# Parse a VCF and store variants in a dictionary
# Each variant is keyed by:
#   (chrom, pos, ref, alt)
# For each variant record:
#   - which caller detected it
#   - its QUAL score (per caller)
def parse_vcf(path, caller_name):
    variants = {}
    headers = []

    with open_any(path) as f:
        for line in f:
            # Store header lines separately, keep them for output
            if line.startswith("#"):
                headers.append(line)
                continue

            # Split variant columns
            cols = line.strip().split("\t")
            chrom, pos, _id, ref, alt, qual, flt, info = cols[:8]
            pos = int(pos)

            key = (chrom, pos, ref, alt)

            # Convert QUAL to float if possible
            try:
                qval = float(qual) if qual != "." else None
            except ValueError:
                qval = None

            # If this variant has not been seen before, create a record
            if key not in variants:
                variants[key] = {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "callers": set(),   # e.g. {"bcftools", "snippy"}
                    "quals": {},        # per-caller QUAL scores
                }

            # Add caller information
            variants[key]["callers"].add(caller_name)
            variants[key]["quals"][caller_name] = qval

    return headers, variants

# Merge two VCFs (bcftools + snippy)
# Output: a combined VCF where each variant has:
#   - CALLERS = which callers detected it
#   - SCORE   = 2 if both callers, 1 if only one
def main():
    if len(sys.argv) != 4:
        print("Usage: combine_vcfs.py <bcftools.vcf.gz> <snippy.vcf.gz> <out.vcf>", file=sys.stderr)
        sys.exit(1)

    bcftools_vcf = sys.argv[1]
    snippy_vcf   = sys.argv[2]
    out_vcf      = sys.argv[3]

    # Parse both input VCFs
    bc_headers, bc_vars = parse_vcf(bcftools_vcf, "bcftools")
    sn_headers, sn_vars = parse_vcf(snippy_vcf, "snippy")

    # Start with bcftools variants
    merged = bc_vars

    # Add snippy variants into the merged dictionary
    # If variant exists in both → combine caller info
    for key, val in sn_vars.items():
        if key in merged:
            merged[key]["callers"].update(val["callers"])
            merged[key]["quals"].update(val["quals"])
        else:
            merged[key] = val

    # Build output headers
    # Use bcftools headers and add new INFO fields for:
    #   CALLERS=...
    #   SCORE=...
    new_headers = []
    for line in bc_headers:
        new_headers.append(line)

    new_headers.append('##INFO=<ID=CALLERS,Number=.,Type=String,Description="Variant callers supporting this record">\n')
    new_headers.append('##INFO=<ID=SCORE,Number=1,Type=Integer,Description="Support score: 2=both callers, 1=single caller">\n')

    # Ensure the final column header line exists
    if not any(line.startswith("#CHROM") for line in new_headers):
        new_headers.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    # Write final combined VCF
    with open(out_vcf, "w") as out:
        # Write all header lines
        for line in new_headers:
            out.write(line)

        # Write variants sorted by chromosome and position
        for key in sorted(merged.keys(), key=lambda x: (x[0], x[1])):
            v = merged[key]

            chrom = v["chrom"]
            pos   = v["pos"]
            ref   = v["ref"]
            alt   = v["alt"]

            callers = sorted(v["callers"])
            score   = len(callers)   # 1 or 2

            # Use the maximum QUAL across caller
            qual_values = [q for q in v["quals"].values() if q is not None]
            qual_str = f"{max(qual_values):.2f}" if qual_values else "."

            info = f"CALLERS={','.join(callers)};SCORE={score}"

            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual_str}\t.\t{info}\n")

# Run main() only when script is called from the command line
if __name__ == "__main__":
    main()
