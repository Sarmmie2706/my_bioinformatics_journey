#!/usr/bin/bash

set -eu

# Run the HaplotypeCaller that detects germline mutation. The tp53 region on the chromosome was selected
gatk HaplotypeCaller \
    -R ../02_alignment/reference/chr17.fa \
    -I ../02_alignment/results/sorted/SRR062634_markdup.bam \
    -O results/variants/germline/SRR062634_germline.vcf.gz \
    --sample-name SRR062634 \
    -L chr17:7661779-7687538

# Count the number of variants
bcftools view -H results/variants/germline/SRR062634_germline.vcf.gz | wc -l

# Run the HaplotypeCaller that detects germline mutation. The tp53 region on the chromosome was also selected
gatk Mutect2 \
    -R ../02_alignment/reference/chr17.fa \
    -I ../o2_alignment/results/sorted/SRR062634_markdup.bam \
    -tumor SRR062634 \
    -O results/variants/somatic/SRR062634_somatic.vcf.gz \
    -L chr17:7661779-7687538 \
    --af-of-alleles-not-in-resource 0.0000025

# Count the number of variants
bcftools view -H results/variants/somatic/SRR062634_somatic.vcf.gz | wc -l
# Both results showed zero variants. So a different dataset was downloaded along with the full hg38 genome to be used as reference

# Reference genome downloaded, unzipped, and indexed
wget -P data/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip data/hg38.fa.gz 
samtools faidx data/hg38.fa 
samtools dict data/hg38.fa -o data/hg38.dict

# Dataset to run variant calling on
wget -P data/ https://42basepairs.com/download/s3/gatk-test-data/wgs_bam/NA12878_20k_hg38/NA12878.bam
wget -P data/ https://42basepairs.com/download/s3/gatk-test-data/wgs_bam/NA12878_20k_hg38/NA12878.bai

# Run the Haplotype Caller on the new dataset
gatk HaplotypeCaller \
    -R data/hg38.fa \
    -I data/NA12878.bam \
    -O results/variants/NA12878_germline.vcf.gz

# Count the number of germline variants called
bcftools view -H results/variants/NA12878_germline.vcf.gz | wc -l

# Run the Mutect2 caller for somatic mutation and then filter using FilterMutectCalls
gatk Mutect2 \
     -R data/hg38.fa \
     -I data/NA12878.bam \
     -tumor NA12878 \
     -O results/variants/NA12878_somatic.vcf.gz \
     --af-of-alleles-not-in-resource 0.0000025

gatk FilterMutectCalls \
     -R data/hg38.fa \
     -V results/variants/NA12878_somatic.vcf.gz \
     -O results/variants/NA12878_somatic_filtered.vcf.gz

# Count the number of somatic variations, unfiltered and filtered
bcftools view -H results/variants/NA12878_somatic.vcf.gz | wc -l
bcftools view -H results/variants/NA12878_somatic_filtered.vcf.gz | wc -l
# None of the variants passed the Mutect2 calls so manual filtering will be done using other metrics

# Filtering based on read depth (DP), allele frequency (AF), and Tumour log odds (TLOD) and view file
bcftools view \
    -i 'FORMAT/DP[0:0]>=5 && FORMAT/AF[0:0]>=0.05 && INFO/TLOD>=6.3' \
     results/variants/NA12878_somatic.vcf.gz \
    -o results/variants/NA12878_somatic_filtered.vcf.gz -Oz

bcftools view -H results/variants/NA12878_somatic_filtered.vcf.gz | head 

# Extract important columns and save to a table and view contents
bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%AF]\n' \
     results/variants/NA12878_somatic_filtered.vcf.gz \
    -o results/variants/NA12878_somatic_variants_table.tsv

head -20 results/variants/NA12878_somatic_variants_table.tsv 

# Create header file and concatenate to the saved variant table
echo -e "CHROM\tPOS\tREF\tALT\tDP\tAF" > results/variants/header.txt
cat results/variants/header.txt results/variants/NA12878_somatic_variants_table.tsv > results/variants/NA12878_somatic_variant_table.tsv 
head -20 results/variants/NA12878_somatic_variant_table.tsv 

# Remove variant table without header file to prevent confusion
rm results/variants/NA12878_somatic_variants_table.tsv 

# SAve stats and view contents using grep
bcftools stats results/variants/NA12878_somatic_filtered.vcf.gz > results/variants/NA12878_somatic_stats.txt
grep "^SN" results/variants/NA12878_somatic_stats.txt 
grep "^TSTV" results/variants/NA12878_somatic_stats.txt 

# Normalise and view normalised data
bcftools norm -f data/hg38.fa -m -any results/variants/NA12878_somatic_filtered.vcf.gz -o results/variants/NA12878_somatic_normalised.vcf.gz -Oz
bcftools view -H results/variants/NA12878_somatic_normalised.vcf.gz | wc -l





