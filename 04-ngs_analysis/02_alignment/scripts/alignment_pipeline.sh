#!/usr/bin/bash

set -eu

# Activate environment and change directories
conda activate bioinfo
cd projects/bioinformatics_portfolio/04-ngs_analysis/02_alignment

# Create directories
mkdir -p reference scripts results/{aligned,sorted,stats} logs

# Download and unzip chr 17 of the reference genome
wget -P reference/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr17.fa.gz
gunzip reference/chr17.fa.gz
head -20 reference/chr17.fa
wc -l reference/chr17.fa

# Crate indexfiles for the unzipped reference
samtools faidx reference/chr17.fa
bwa index reference/chr17.fa
samtools dict reference/chr17.fa -o reference/chr17.dict

# View the newly created index files
cat reference/chr17.fa.fai

bwa mem \
    -t 2 \
    -R "@RG\tID:SRR062634\tSM:SRR062634\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
    reference/chr17.fa \
    ../week1_qc/results/fastp/SRR062634_1_clean.fastq.gz \
    ../week1_qc/results/fastp/SRR062634_2_clean.fastq.gz \
    > results/aligned/SRR062634.sam

# Convert the sam to a bam file and delete the sam file to save space
samtools view -bS     results/aligned/SRR062634.sam     -o results/aligned/SRR062634.bam

# View first 10 alignments (human readable)
samtools view results/aligned/SRR062634.bam | head -10

# View with header
samtools view -H results/aligned/SRR062634.bam

# View only mapped reads
samtools view -F 4 results/aligned/SRR062634.bam | head -10

# View first 10 unmapped reads
samtools view -f 4 results/aligned/SRR062634.bam | head -10
# -F is used for exclusion while -f is used for inclusion

samtools sort \
    results/aligned/SRR062634.bam \
    -o results/sorted/SRR062634_sorted.bam \
    -@ 2

# The alignment and sorting can be done in the same step and skipping the intermediate sam file and unsorted bam file
# bwa mem -t 4 \
#    -R "@RG\tID:SRR062634\tSM:SRR062634\tPL:ILLUMINA\tLB:lib1\tPU:unit1" \
#    reference/chr17.fa \
#    ../week1_qc/results/fastp/SRR062634_1_clean.fastq.gz \
#    ../week1_qc/results/fastp/SRR062634_2_clean.fastq.gz \
#    | samtools sort -@ 4 -o results/sorted/SRR062634_sorted.bam

# Create a markdup (mark duplicates) bam file and index it
samtools sort -n -@ 2 results/sorted/SRR062634_sorted.bam \
    | samtools fixmate -m - - \
    | samtools sort -@ 2 \
    | samtools markdup -@ 2 - results/sorted/SRR062634_markdup.bam
samtools index results/sorted/SRR062634_markdup.bam

# Check the statistics of the markdup bam file
samtools flagstat results/sorted/SRR062634_markdup.bam     > results/stats/SRR062634_flagstat.txt
cat results/stats/SRR062634_flagstat.txt

samtools idxstats results/sorted/SRR062634_markdup.bam \
    > results/stats/SRR062634_idxstats.txt
cat results/stats/SRR062634_idxstats.txt

# Checking average depth across chr17
samtools depth results/sorted/SRR062634_markdup.bam \
    | awk '{sum+=$3; count++} END {print "Mean depth:", sum/count}' 

# Coverage at a specific region 
samtools depth -r chr17:7661779-7687538 \
    results/sorted/SRR062634_markdup.bam \
    | head -20

# Keep only properly paired reads (FLAG 2)
samtools view -b -f 2 \
    results/sorted/SRR062634_markdup.bam \
    -o results/sorted/SRR062634_proper_pairs.bam

# Remove duplicates and unmapped reads
# -F 1804 excludes: unmapped(4) + mate unmapped(8) + not primary(256) + duplicate(1024) + supplementary(512)
# 4+8+256+1024+512 = 1804
samtools view -b -F 1804 -f 2 \
    results/sorted/SRR062634_markdup.bam \
    -o results/sorted/SRR062634_filtered.bam

# Filter by mapping quality (keep MAPQ >= 20)
samtools view -b -q 20 \
    results/sorted/SRR062634_markdup.bam \
    -o results/sorted/SRR062634_mapq20.bam

# Extract only reads mapping to TP53 locus
samtools view -b \
    results/sorted/SRR062634_markdup.bam \
    chr17:7661779-7687538 \
    -o results/sorted/SRR062634_tp53.bam

# Indexing and calculating statistics of extracted part
samtools index results/sorted/SRR062634_tp53.bam
samtools flagstat results/sorted/SRR062634_tp53.bam

# Mean depth coverage
samtools depth results/sorted/SRR062634_tp53.bam     | awk '{sum+=3; count++} END {print Mean depth:, sum/count}'

# Reads that mapped to the region
samtools view -c results/sorted/SRR062634_markdup.bam chr17:7661779-7687538

# Checking full comprehensive QC statistics of the alignment
samtools stats results/sorted/SRR062634_markdup.bam \
    > results/stats/SRR062634_stats.txt

# Summary numbers only
grep "^SN" results/stats/SRR062634_stats.txt | head -30

# Open igv using the command igv to see the alignment visually. 
igv
# Note: This can also be donre directly from the destop app or the web app. It only requires loading the reference genome, alignment file and their indexes
