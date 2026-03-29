library(tidyverse)
library(patchwork)

set.seed(99)
n_mutations <- 200

gene_pool <- c(
  rep("TP53",   25), rep("KRAS",   20), rep("BRCA1",  15),
  rep("EGFR",   14), rep("MYC",    12), rep("PTEN",   10),
  rep("RB1",     9), rep("APC",     8), rep("BRAF",    7),
  rep("PIK3CA",  6), rep("CDKN2A",  5), rep("ERBB2",   4),
  rep("MDM2",    4), rep("VEGFA",   3), rep("CDH1",    3),
  rep("OTHER",  55)
)

mutations <- data.frame(
  mutation_id   = paste0("MUT", 1:n_mutations),
  gene          = sample(gene_pool, n_mutations, replace=FALSE),
  sample_id     = paste0("S", sample(1:40, n_mutations, replace=TRUE)),
  mutation_type = sample(c("SNP","Indel","CNV"),
                         n_mutations, replace=TRUE,
                         prob=c(0.60, 0.25, 0.15)),
  chromosome    = sample(1:22, n_mutations, replace=TRUE),
  vaf           = c(
    rbeta(120, 2, 3),       
    rbeta(50,  1.5, 4),      
    rbeta(30,  3, 2)       
  ),
  cancer_type   = sample(c("Lung","Breast","Colorectal"),
                         n_mutations, replace=TRUE,
                         prob=c(0.40, 0.35, 0.25))
)
gene_pool
mutations
# Make VAF realistic — scale to 0-1
mutations$vaf <- round(mutations$vaf, 3)

head(mutations)
str(mutations)
summary(mutations)

cat(sprintf("%s %d", "Total number of mutations:", nrow(mutations)))
cat(sprintf("%s %d", "Total number of genes affected:", length(unique(mutations$gene))))
cat(sprintf("%s %d", "Total number of samples:", length(unique(mutations$sample_id))))

table(mutations$mutation_type)
prop.table(table(mutations$mutation_type)) * 100

mean(mutations$vaf)
median(mutations$vaf)

mutations %>% 
  summarise(
    mean_vaf = mean(vaf),
    .by = mutation_type
  ) %>% 
  print()































