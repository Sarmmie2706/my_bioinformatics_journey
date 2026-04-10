#!/usr/bin/bash

# All the file opening commands are commented out to prevent multiple tabs opening at once. The commands are present in the open_reports.sh to be run when review is to be done. 
# Activate the environment
conda activate bioinfo

# Change to your qc directory
cd projects/bioinformatics_portfolio/04-ngs_analysis/01_qc/

# Download and subset your data
wget -P data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062634/SRR062634_1.fastq.gz
wget -P data/ ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR062/SRR062634/SRR062634_2.fastq.gz
zcat data/SRR062634_1.fastq.gz | head -80000 | gzip > data/SRR062634_1_subset.fastq.gz
zcat data/SRR062634_1.fastq.gz | head -80000 | gzip > data/SRR062634_2_subset.fastq.gz

# Run fastqc on the subset and check the results
fastqc data/SRR062634_1_subset.fastq.gz data/SRR062634_2_subset.fastq.gz -o results/fastqc_raw/
# xdg-open results/fastqc_raw/SRR062634_1_subset_fastqc.html
# xdg-open results/fastqc_raw/SRR062634_2_subset_fastqc.html

# Create a folder for the fastp results, run the fastp and output the results into it
mkdir -p results/fastp/

fastp \
    --in1 data/SRR062634_1.fastq \
    --in2 data/SRR062634_2.fastq \
    --out1 results/fastp/SRR062634_1_clean.fastq.gz \
    --out2 results/fastp/SRR062634_2_clean.fastq.gz \
    --json results/fastp/SRR062634_fastp.json \
    --html results/fastp/SRR062634_fastp.html \
    --thread 4 \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --unqualified_percent_limit 40 \
    --length_required 36 \
    --cut_tail \
    --cut_tail_window_size 4 \
    --cut_tail_mean_quality 20

# Check the fastp results 
# xdg-open results/fastp/SRR062634_fastp.html

# Run fastqc on the trimmed results and check the results
fastqc results/fastp/SRR062634_1_clean.fastq.gz results/fastp/SRR062634_2_clean.fastq.gz -o results/fastqc_trimmed/
# xdg-open results/fastqc_trimmed/SRR062634_1_clean_fastqc.html
# xdg-open results/fastqc_trimmed/SRR062634_2_clean_fastqc.html

# Run multiqc on the fastqc and fastp results
multiqc results/ -o results/multiqc/ -n week1_multiqc_report
# xdg-open results/multiqc/week1_multiqc_report.html
