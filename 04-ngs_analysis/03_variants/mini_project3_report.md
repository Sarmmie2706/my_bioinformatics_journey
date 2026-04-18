# Capstone Project — Somatic Variant Detection and Annotation
**Programme:** 3-Month Cancer Bioinformatics Learning Schedule  
**Phase:** Phase 3 — NGS Analysis & Variant Calling  
**Date:** April 2026  
**Author:** Samuel Aladegbaiye

---

## 1. Introduction

### 1.1 Background
Cancer is fundamentally a disease of the genome. Somatic mutations are
acquired changes seen in the DNA tumour cells but not normal cells. In addition to
initiating tumour formation and progression, they also affect response to treatment. Identifying
these mutations from sequencing data is therefore important in cancer genomics research.

### 1.2 TP53 in Cancer
TP53 (tumour protein p53) is the most commonly mutated gene in human cancer,
altered in approximately 50% of all tumours. It encodes p53, a protein that regulates cell cycle, DNA repair, apoptosis, and stress responses.
Loss-of-function mutations in TP53 allow damaged cells to proliferate unchecked, driving cancer development.
TP53 mutations are common in lung, breast, colorectal, ovarian, and head & neck cancers.

### 1.3 Aims
1. Perform end-to-end quality control on raw Illumina sequencing data
2. Align reads to the human reference genome (hg38)
3. Call germline and somatic variants using GATK best practices tools
4. Annotate and interpret variants in the context of cancer biology
5. Document all steps, findings, and limitations transparently

---

## 2. Dataset

| Parameter | Details |
|-----------|---------|
| Sample | NA12878 (Genome in a Bottle reference sample) |
| Accession | NA12878_20k_hg38 |
| Source | GATK Test Data (42basepairs mirror) |
| Data type | Whole genome sequencing (WGS), Illumina paired-end |
| Reference genome | hg38 (UCSC, chr prefix naming) |
| Note | NA12878 is a germline reference sample, not a true tumour sample. Somatic calling was performed for learning purpose.

The whole analysis began with a 20k subset of SRR062634_1.fastq.gz, but was discontinued after
variant calling returned zero variants for both germline and somatic mutations. The details of the first variant are in the table below. 


**Additional dataset (Week 1 QC):**
| Parameter | Details |
|-----------|---------|
| Sample | SRR062634 |
| Source | ENA (downloaded directly as fastq.gz) |
| Subset used | 20,000 reads |

---

## 3. Software & Versions
| Tool | Version | Purpose |
|------|---------|---------|
| FastQC | FastQC v0.12.1 | Raw read quality assessment |
| fastp | fastp 1.3.2 | Adapter trimming and quality filtering |
| BWA-MEM | Version: 0.7.19-r1273 | Read alignment |
| SAMtools | samtools 1.22.1 | BAM processing and statistics |
| GATK4 | (GATK) v4.6.2.0 | Variant calling |
| bcftools | bcftools 1.22 | VCF manipulation and filtering |

---

## 4. Methods

### 4.1 Pipeline Overview
Raw FASTQ
↓
Quality Control (FastQC → fastp → FastQC → MultiQC)
↓
Alignment (BWA-MEM → namesort → fixmate → sort → markdup)
↓
BAM Processing (index → flagstat → idxstats → stats)
↓
Variant Calling (HaplotypeCaller + Mutect2)
↓
Variant Filtering (FilterMutectCalls / manual bcftools)
↓
Annotation (VEP / manual COSMIC + ClinVar lookup)
↓
Interpretation

### 4.2 Quality Control
- Raw reads assessed with FastQC
- Adapter trimming and quality filtering performed with fastp
- Settings: Q20 quality threshold, minimum length 36bp, automatic
  adapter detection, 3' tail cutting with window size 4, mean quality 20
- Post-trimming FastQC run to confirm improvement
- Reports aggregated with MultiQC

### 4.3 Alignment
- Reference: hg38 chr17 (SRR062634)
- Aligner: BWA-MEM with read group tags
- Post-alignment processing: namesort → samtools fixmate → coordinate
  sort → samtools markdup
- Initial markdup attempt failed with "no ms score tag" error.
  This was fixed by running samtools fixmate before markdup to add required
  mate score tags.

### 4.4 Variant Calling
- Germline variants: GATK HaplotypeCaller restricted to TP53 locus
  (chr17:7,661,779-7,687,538) on SRR062634
- Somatic variants: GATK Mutect2 tumor-only mode on NA12878, genome-wide
- BQSR (Base Quality Score Recalibration) was not performed because
  of resource constraints (requires large dbSNP/known sites files). The commands required to run
  it are in the script file.
- The variant calling on SRR062634 returned zero variants. This is due to the subsetting of the
  data and using chr17 alone as the reference.
- A new file (NA12878.bam) was then used for variant calling, and this file returned variants for
  both the germline and somatic mutations.

### 4.5 Variant Filtering
- FilterMutectCalls was run on the somatic mutation but it produced an empty output VCF.
- As a result of this, manual filtering had to be done using bcftools with thresholds:
  - FORMAT/DP ≥ 5 (minimum read depth)
  - FORMAT/AF ≥ 0.05 (minimum allele frequency 5%)
  - INFO/TLOD ≥ 6.3 (GATK default tumor log odds threshold)

### 4.6 Variant Annotation
The VEP results was queried on the Ensembl website directly because the offline version required downloading 
caches with very large size (~10GB), and the online version had issues connecting with Ensembl's servers. The 
query was run by uploading the variant calling result and using the default settings on the site. 

---

## 5. Results
### 5.1 Variant Calling

#### Germline — HaplotypeCaller (NA12878, genome-wide)

| Metric | Value |
|--------|-------|
| Region | Whole genome |
| Variants called | 1590 |
| Interpretation | Low/zero calls expected due to insufficient depth |

#### Somatic — Mutect2 (NA12878, genome-wide)

| Stage | Variant Count |
|-------|--------------|
| Raw Mutect2 output | 13004 |
| After FilterMutectCalls | 0 (tool produced empty output) |
| After manual bcftools filtering | 86 |
The massive drop in number of variant calls after filtering is explained by the
points listed the subsection 5.2

**Manual filter thresholds applied:**
- FORMAT/DP[0:0] ≥ 5
- FORMAT/AF[0:0] ≥ 0.05
- INFO/TLOD ≥ 6.3

### 5.2 VCF Quality Metrics

| Metric | Value |
|--------|-------|
| Somatic variants (filtered) | 86 |
| SNPs | 76 |
| Indels | 11 |
| Ts/Tv ratio | 0.90 / 0.93 |

**Ts/Tv Interpretation:**
The Ts/Tv ratio of ~0.9 is significantly below the expected 2.0-2.1
for human WGS data. This indicates a high false positive rate in the
callset. This can be attributed to the fact that:
- BQSR was not performed which means that there would be uncorrected base quality scores
- There was no matched normal sample meaning there would be germline variants contaminating somatic calls
- NA12878 is a germline reference sample, not a true tumour
- Manual filtering used instead of FilterMutectCalls

### 5.3 Variant Annotation
| Gene | HGVSc | HGVSp | rsID | SIFT | PolyPhen | ClinVar | COSMIC | gnomAD |
|------|-------|-------|------|------|----------|---------|--------|--------|
| MT-CYB | c.85G>A | p.Ala29Thr | rs199795644 | tolerated (0.16) | benign (0) | Different alt allele at same position (m.14831G>C, VCV000693773.1) — not a direct match | Not found | Not found |
| MT-CYB | c.580A>G | p.Thr194Ala | rs2853508 | tolerated (0.21) | benign (0.009) | Not found | Not found | Not found |

### 5.4 Database Lookup — MT-CYB Variants

Both MT-CYB variants were queried against COSMIC, ClinVar, and 
gnomAD using rsID and genomic position. The results are interpreted
and explained in the Discussion section.

---

## 6. Discussion

### 6.1 Key Findings
The complete somatic variant calling pipeline was successfully 
implemented from raw FASTQ through to annotated variants. However, 
no clinically significant somatic cancer mutations were identified. 
This outcome is expected and informative given the nature of the 
dataset used.

Out of the 86 variants that passed manual quality filters, VEP annotation 
showed that the majority fell in intergenic regions or 
mitochondrial genes. Only two variants carried MODERATE impact 
predictions and both were missense variants in MT-CYB (p.Ala29Thr and 
p.Thr194Ala). No HIGH impact variants were identified, and no 
variants were found in canonical nuclear cancer driver genes such 
as TP53, BRCA1, BRCA2, or KRAS.

### 6.2 Database Lookup Summary
Systematic querying of COSMIC, ClinVar, and gnomAD confirmed that 
neither MT-CYB variant represents a known somatic cancer mutation.

**COSMIC:** Neither variant was found at the queried positions.
Other MT-CYB mutations of the same substitution type (G>A, A>G)
exist in the database but at different positions. This confirms these
are not known recurrent somatic cancer mutations.

**ClinVar:** The first variant (c.85G>A, rs199795644) returned a
ClinVar entry at the same genomic position (m.14831) but with a
different alternate allele (G>C vs G>A) which is not a direct match.
The second variant (c.580A>G, rs2853508) was not found in ClinVar.
Therefore, there was no way to assign Direct ClinVAr association to
either of the variants.

**gnomAD:** Neither variant was retrieved by rsID search. This could be because of
the complexity of mitochondrial variant representation in
population databases, where heteroplasmy levels and mitochondrial
copy number make allele frequency calculations complicated, instead of
indicating that the variants are absent from the population.

**Overall interpretation:** Both variants are predicted benign by
SIFT and PolyPhen. The absence from COSMIC, combined with benign
in silico predictions, supports the conclusion that these are
natural mitochondrial heteroplasmy in a healthy individual (NA12878)
and not somatic cancer driver mutations. Their presence in the
somatic callset is consistent with a well known limitation of
tumor-only Mutect2 calling without a matched normal sample.

### 6.3 Why These Results Are Expected
The absence of real somatic cancer mutations is as a result of 
the dataset used. NA12878 is a well-characterised healthy 
reference individual from the Genome in a Bottle consortium and not 
a cancer patient. Mutect2 in tumor-only mode, without a matched 
normal sample, cannot distinguish true somatic mutations from:

- Natural mitochondrial heteroplasmy (the MT-CYB variants that were seen)
- Germline variants present in the healthy individual
- Low-quality sequencing artifacts at low-depth positions

The low Ts/Tv ratio of ~0.9 (expected: ~2.0 for WGS) independently 
confirms a high false positive rate in the raw callset, consistent 
with these limitations.

### 6.4 Limitations of This Analysis

| Limitation | Impact on Results | Real-world Solution |
|------------|-------------------|-------------------|
| NA12878 is a healthy germline sample | No true somatic mutations exist to find | Use matched tumour-normal pairs from a cancer patient |
| No matched normal sample | Germline variants and MT heteroplasmy called as somatic | Sequence matched blood/adjacent normal tissue |
| No BQSR performed | Uncorrected base quality errors inflate false positive rate | Download GATK resource bundle, run BQSR before calling |
| FilterMutectCalls produced empty output | Manual bcftools filtering used instead | Investigate GATK logs, ensure stats file generated correctly |
| 20k read subset for QC phase | Near-zero coverage on chr17, no variants called | Use full dataset at clinical sequencing depth (30-60x) |
| Single chromosome reference (chr17) for SRR062634 | 28% mapping rate, low proper-pairing rate | Use full hg38 reference |
| VEP run via web interface, not command line | Cannot be scripted or automated | Download 20GB hg38 VEP cache for offline use |
| No panel of normals (PoN) for Mutect2 | Higher false positive somatic calls | Build PoN from multiple normal samples |

### 6.5 What a Real Clinical Pipeline Would Look Like
A production cancer genomics pipeline would include:

- Matched tumour-normal WGS at 60x (tumour) and 30x (normal)
- Full hg38 reference with BQSR using dbSNP and known indel sites
- Panel of Normals built from multiple healthy individuals
- Mutect2 in tumour-normal mode for clean somatic subtraction
- Automated VEP annotation with local cache including COSMIC, 
  ClinVar, and gnomAD
- Mutational signature analysis (SigProfiler) to identify 
  mutagenic processes
- Variant interpretation by a trained clinical molecular biologist
- Structured clinical report with actionability tiers (OncoKB)

Despite these limitations, this analysis successfully demonstrated 
every step of the somatic variant calling workflow, from raw data 
through to biological interpretation, and developed a clear 
understanding of why each step exists and what can go wrong when 
best practices are not fully followed. 

---

## 7. Conclusions
- A complete somatic variant calling pipeline was implemented from
  raw FASTQ to annotated variants
- Quality control with fastp improved Q30 rate from 79% to 85%
- Alignment to chr17 gave expected low mapping rate due to
  single-chromosome reference
- No variants were called on the learning dataset (SRR062634) due
  to insufficient depth — consistent with expectations
- 86 variants passed manual filtering on NA12878 genome-wide data
- Low Ts/Tv ratio (0.9) reflects the limitations of tumor-only calling
  without BQSR on a germline sample
- All pipeline limitations are understood and documented with
  real-world solutions identified

---

## 8. What I Learned

- Getting zero variants on the first data revealed the importance of sequencing depth in NGS analysis
- RUnning a somatic variants call without a matched normal will not yield optimal results, as the comparison is between normal and tumour cells
- With most of the pipeline running correctly but showing suboptimal results, this analysis shows the need to use the right data for analysis

---

## Pipeline Scripts
- Week 1 QC: `../01_qc/scripts/pipeline.sh`
- Week 2 Alignment: `../02_alignment/scripts/alignment_pipeline.sh`
- Week 3 Variants: `scripts/variant_pipeline.sh`

