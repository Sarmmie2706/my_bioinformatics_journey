# Mini Project 1 — NGS Quality Control Pipeline
**Date:** April 2026  
**Dataset:** SRR062634 (subset: 20,000 reads)  
**Source:** ENA (downloaded directly as fastq.gz)  
**Tools Used:** FastQC, fastp, MultiQC

---

## 1. Objective
Perform quality control on raw Illumina sequencing reads from a public 
dataset, assess quality metrics, apply trimming/filtering, and compare 
results before and after processing.

---

## 2. Dataset
- Accession: SRR062634
- Downloaded from: ENA FTP
- Subset used: 20,000 reads (80,000 lines) per file
- Files: SRR062634_1_subset.fastq.gz, SRR062634_2_subset.fastq.gz

---

## 3. Raw Read Quality (FastQC)
| Metric | R1 | R2 |
|--------|----|----|
| Total sequences | 20,000 | 20,000 |
| Sequence length | 100 | 100 |
| Adapter content | ~0% | ~0% |


**Key observations:**
- The per base sequence quality failed before the trimming, with a considerable amount of the later read
positions being of poor quality. Both sequences had a GC% of 40% getting a "Warn" status.
- Read 2's per base sequence content also got a "Warn" status. All the other checks passed.


---

## 4. Trimming & Filtering with fastp
**Settings used:**
- Qualified quality threshold: Q20
- Unqualified percent limit: 40%
- Minimum length: 36bp
- 3' tail cutting: window size 4, mean quality 20
- Adapter detection: automatic (paired-end)

**Results:**
| Metric | Before | After |
|--------|--------|-------|
| Total reads | 20,000 | 19,214 |
| Total bases | ~2 Mbp | ~1.7–1.8 Mbp |
| Reads removed | — | 786 (~3.9%) |
| Q30 bases | 79.20% | 85.80% |
| GC content | 40.37% | 39.33% |
---

## 5. Post-Trimming Quality (FastQC)
| Metric | R1 | R2 |
|--------|----|----|
| Sequence length | 36-100 | 36-100 |
| Total sequences | 19,214 | 19,214 |
| Adapter content | ~0% | ~0% |

**Key observations:**
- Per base sequence quality improved, particularly at 3' end
- Adapter content reduced further after trimming
- While everyother metric improved, the GC content reduced by 1.04%

---

## 6. MultiQC Summary
In addition to the results shown by the individual fastqc and fastp results. The following values were seen from running multiqc:
| Metric | R1 | R2 | R1_trimmed | R2_trimmed |
|--------|----|----|------------|------------|
| % Duplication | 0.2% | 0.2% | 0.1% | 0.2% |


---

## 7. Conclusions
- Raw data quality can be considered acceptable with Q30 rate of 79%
- Adapter contamination was minimal in this dataset
- fastp filtering improved Q30 to 85% and removed ~3.9% of reads
- Data is now ready for alignment to the reference genome

---

## 8. What I Learned
- In this project, I learned the metrics that sequences DNA are judged based on, the accepted values or ranges,
and what could cause deviations from these ranges, and what these deviations could mean biologically.
- I also learned the tools used to correct errors seen to ensure the good quality of sequences for downstream analysis
- Phred scores are log-scale so Q30 = 99.9% accuracy, not just "30% good"
- Near-zero adapter content doesn't always mean no trimming needed,  quality trimming is still valuable

---

## Pipeline Script
See: `scripts/pipeline.sh`
