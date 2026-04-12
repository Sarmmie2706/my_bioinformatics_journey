# Mini Project 2 — Read Alignment Pipeline
**Date:** April 2026  
**Dataset:** SRR062634 (subset: 20,000 reads). It's the same one used for the QC analysis. Check repo for details  
**Reference:** hg38 chr17 (UCSC)  
**Tools Used:** BWA-MEM, SAMtools, IGV

---

## 1. Objective
To align quality-controlled Illumina reads to the human reference genome 
(chr17), process the resulting BAM file, and visualise alignments at 
the TP53 locus using IGV.

---

## 2. Reference Genome
- Assembly: hg38 (UCSC naming, chr prefix)
- Chromosome used: chr17 (83,257,441 bp)
- Indexed with: samtools faidx, bwa index, samtools dict

---

## 3. Alignment — BWA-MEM
**Command used:** BWA-MEM with read group tags, and then piped through 
namesort → fixmate → coordinate sort → markdup

**Read group metadata:**
- ID: SRR062634
- Sample: SRR062634
- Platform: ILLUMINA
- Library: lib1

**Note:** The markdup had to be arun a second time because the initial markdup attempt failed with "no ms score tag" error. 
This was fixed by running samtools fixmate before markdup, which adds the 
required mate score tags.

---

## 4. Alignment Statistics (flagstat)

| Metric | Value |
|--------|-------|
| Total reads | 38,756 |
| Primary reads | 38,428 |
| Supplementary | 328 |
| Duplicates marked | 11 (0.028%) |
| Mapped reads | 10,944 (28.24%) |
| Properly paired | 5,920 (15.41%) |
| Singletons | 1,448 (3.77%) |

**Interpretation:**
The mapping rate of 28.24% is expected and not a quality concern. 
It is that low because the FASTQ files contain reads from the entire human genome, but the 
reference used is chr17 only, for space reasons.So reads that come from other chromosomes 
(26,364 reads as shown by idxstats) will have nowhere to map and are reported as 
unmapped. With a full hg38 reference, the mapping rate is expected to be >95%.

The same issue is seen in the properly paired rate of 15.41%.If both mates 
don't map to chr17, then they can't be considered properly paired.  


---

## 5. idxstats

| Chromosome | Length | Mapped Reads | Unmapped |
|------------|--------|--------------|---------|
| chr17 | 83,257,441 | 10,944 | 1,448 |
| * (other chrs) | 0 | 0 | 26,364 |

---

## 6. Duplicate Marking
Just 11 duplicates were marked, representing 0.028% of total reads
For a 20k read subset, duplication is expected to be minimal.

---

## 7. TP53 Region
- Coordinates: chr17:7,661,779-7,687,538
- Reads mapping to this region: 8
- Mean coverage at locus: 1.061%

---

## 8. IGV Visualisation
*To be completed tomorrow after reviewing alignment visualization.*

---

## 9. Conclusions
- The alignment was generally successful against chr17 reference
- Low mapping rate (28%) is understandable for the single-chromosome reference that was used
- Data is ready for variant calling

---

## 10. What I Learned
- I learnt the tools used for aligning read against a reference genome.
- I also learnt how to check metrics that qualify the alignment, as well as how to filter based on these metrics and extract certain regions
- The correct samtools markdup order is namesort → fixmate → sort → markdup
- Mapping rate depends heavily on whether your reference matches your data's origin

---

## Pipeline Script
See: `scripts/alignment_pipeline.sh`
