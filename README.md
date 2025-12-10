# RP1B: TASK 2
This repository contains my implementation of Task 2, which involves:  
1. Generating simulated genomic mutations and validating variant calling accuracy.    
2. Building a small SNP/indel-calling pipeline using minimap2, bcftools, and Snippy.
#### Output files location: `jovyan:/shared/team/people/nella/task_2`
---

# Part 1 - Validation Using Simulated Data

## 1. Generating Mutations
I used (`make_mutations.py`) that introduces:
- 300 SNPs
- 20 small indels (1–10 bp)

I applied this to:
- `NC_037282.1.fasta` (*Plasmodium falciparum* reference genome)

The script outputs:
- A mutated FASTA  
- A truth table of mutations (`NC_037282.1.truth.tsv`)

---

## 2. Simulating Illumina Reads
I used `simulate_paired_end_reads.py` to generate:
- 30× depth  
- 100 bp perfect paired-end reads  

Outputs:
- `R1.fq`
- `R2.fq`

---

## 3. Mapping with minimap2 and Calling Variants (bcftools)
Using:   
```bash
minimap2 -ax sr ref.fasta R1.fq R2.fq > sample.sam   
samtools sort -o sample.bam sample.sam   
samtools index sample.bam   
bcftools mpileup -Ou -f ref.fasta sample.bam | bcftools call -mv -Oz -o sample.vcf.gz
```
---

## 4. Comparing Expected vs Observed Variants

I evaluated results using `evaluate_variants.py`.

### Final accuracy metrics:
#### (TP = True Positives, FP = False Positives, FN = False Negatives)
```
=== SNPs ===   
TP = 300   
FP = 0   
FN = 0   
Precision = 1.000   
Recall = 1.000   

=== INDELs ===   
TP = 18   
FP = 4   
FN = 2   
Precision = 0.818   
Recall = 0.900   
```
### Explanation
- SNP calling was perfect (0 FP, 0 FN)  
- Indel accuracy was lower because short-read data struggled around repetitive regions, making alignment ambiguous and causing callers to misplace gaps or confuse true indels with sequencing errors

---

# Part 2 - Variant Calling Pipeline

The pipeline I implemented does the following:

1. Takes a reference genome together with FASTQ reads as input  
2. Maps reads using minimap2  
3. Calls variants with:  
   - bcftools  
   - Snippy  
4. For the simulated *Plasmodium falciparum* dataset, combines the bcftools and Snippy calls using `combine_vcfs.py`  
5. Runs on:  
   - Simulated *Plasmodium falciparum* data (with a known truth set)  
   - Real *E. coli* data (no truth set available)  
6. Assigns confidence scores **for the simulated dataset only** using `combine_vcfs.py`, which marks variants based on caller agreement

---

## 1. Calling Variants on Simulated Data

### Variant counts

|    Caller    | SNP/Indel Count |
|--------------|-----------------|
|   bcftools   |       320       |
|    snippy    |       318       |
| combined VCF |       330       |

### Why does the combined VCF have more?
Because it merges union of all positions found by either caller.

### Is the combined VCF better?

Yes, partly.

- Combined VCF recovered all true positives from both callers  
- It included a small number of extra positions (expected when merging callsets)  
- For SNPs, both callers were accurate; combining added redundancy and increased recall 
- For indels, combining improved sensitivity

---

## 2. Running Pipeline on Real *E. coli* Data

Inputs:
- `SRR25083113_1.fastq.gz`
- `SRR25083113_2.fastq.gz`
- `EcoliK12-MG1655.fasta`

### bcftools result:
71,402 variants

### Snippy result:
55,911 variants

### Interpretation
- Real data contained real sequencing errors, natural variation, and mapping artefacts  
- Snippy was more conservative (filtered low-quality sites)  
- bcftools was more permissive and therefore called more variants

---

## 3. Assigning a Trust Score to Variants

To assign confidence levels to variants in the simulated dataset, I used the script `combine_vcfs.py`, which merged the VCF outputs from bcftools and Snippy.

For each variant, the script examined which callers detected it and assigned a trust score:
- Score = 2 → variant detected by both bcftools and Snippy (high confidence)
- Score = 1 → variant detected by only one caller (lower confidence)

The script also added an INFO field listing which callers detected each variant (e.g., `CALLERS=bcftools,snippy`).

---

## 4. Inspecting Variants with samtools tview

I inspected high and low confidence sites using: 
```bash
samtools tview ecoli_real.sorted.bam EcoliK12-MG1655.fasta
```
### Findings:

### High-confidence variants
- Strong, consistent pileup support
- Reads agreed on the alternate allele
- Clean alignments with minimal clipping
- Reasonably uniform depth

### Low-confidence variants
- Often located near repeats or homopolymers
- Mixed allele signals (not all reads agreed)
- Soft-clipped reads indicating alignment uncertainty
- Abnormally low or high depth, suggested possible sequencing or mapping artefacts

### Conclusion:   
tview qualitatively supported the scoring system; high-confidence variants appeared well supported by the reads, whereas low-confidence variants often displayed patterns typical of false positives or ambiguous alignments.

---

# Discussion - What Worked / What Didn’t

### What worked:
- SNP detection on simulated data achieved perfect accuracy (100% precision & recall)
- Pipeline successfully ran both callers
- Combined VCF improved recall
- Snippy ran correctly after resolving the samtools compatibility issue
- Real *E. coli* data was processed successfully, producing a plausible set of variant calls

### What didn’t work:
- Indel calling was less accurate than SNP calling, as expected, due to alignment ambiguity and short-read limitations
- Snippy initially crashed when used with samtools >1.20, due to a known fixmate-related bug (resolved by creating a new environment)
- Accuracy could not be evaluated on real data because no truth set exists, and the dataset contains noise typical of real sequencing

---

# Final Conclusions

- The pipeline worked well on simulated data and detected SNPs accurately   
- Combining caller outputs slightly improved sensitivity on simulated data
- Real data variant interpretation required additional filtering and visual inspection  
- tview supported the confidence-scoring approach

---

# Files in This Repository

`make_mutations.py`   
`simulate_paired_end_reads.py`   
`evaluate_variants.py`   
`combine_vcfs.py`   
`run_pipeline.sh`   
`README.md`
