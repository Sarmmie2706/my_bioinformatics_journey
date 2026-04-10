#!/usr/bin/bash

# Check raw FastQC results
xdg-open results/fastqc_raw/SRR062634_1_subset_fastqc.html
xdg-open results/fastqc_raw/SRR062634_2_subset_fastqc.html

# Check fastp report
xdg-open results/fastp/SRR062634_fastp.html

# Check trimmed FastQC results
xdg-open results/fastqc_trimmed/SRR062634_1_clean_fastqc.html
xdg-open results/fastqc_trimmed/SRR062634_2_clean_fastqc.html

# Check MultiQC report
xdg-open results/multiqc/week1_multiqc_report.html
