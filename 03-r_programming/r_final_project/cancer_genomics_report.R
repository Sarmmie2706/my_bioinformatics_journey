library(tidyverse)
library(GenomicRanges)
library(org.Hs.eg.db)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(patchwork)
library(corrplot)

set.seed(42)

# Clinical data 
clinical <- data.frame(
  patient_id  = paste0("P", 1:60),
  age         = round(runif(60, 30, 80)),
  sex         = sample(c("Male","Female"), 60, replace=TRUE),
  cancer_type = sample(c("Lung","Breast","Colorectal"),
                       60, replace=TRUE,
                       prob=c(0.40, 0.35, 0.25)),
  stage       = sample(c("I","II","III","IV"),
                       60, replace=TRUE,
                       prob=c(0.15, 0.30, 0.35, 0.20)),
  treatment   = sample(c("Chemo","Targeted","Immunotherapy"),
                       60, replace=TRUE),
  stringsAsFactors = FALSE
)

clinical$survival_time <- c(
  rexp(30, rate=0.035),   # tumour patients — shorter survival
  rexp(30, rate=0.020)    # normal patients — longer survival
)
clinical$survival_time <- round(clinical$survival_time, 1)
clinical$status        <- rbinom(60, 1, 0.65)
clinical$group         <- rep(c("Tumour","Normal"), each=30)

# Expression data — 15 cancer genes, 60 patients 
gene_names <- c("TP53","BRCA1","KRAS","EGFR","MYC",
                "PTEN","RB1","APC","BRAF","PIK3CA",
                "CDKN2A","ERBB2","MDM2","VEGFA","CDH1")

# First 30 patients are tumour, last 30 are normal
expr_matrix <- matrix(
  c(
    rnbinom(15*30, mu=800,  size=1),  # tumour
    rnbinom(15*30, mu=400,  size=1)   # normal
  ),
  nrow = 15,
  ncol = 60
)
rownames(expr_matrix) <- gene_names
colnames(expr_matrix) <- clinical$patient_id

# Add signal — specific genes upregulated in tumour
expr_matrix["TP53",   1:30] <- expr_matrix["TP53",   1:30] + 500
expr_matrix["MYC",    1:30] <- expr_matrix["MYC",    1:30] + 400
expr_matrix["KRAS",   1:30] <- expr_matrix["KRAS",   1:30] + 350
expr_matrix["EGFR",   1:30] <- expr_matrix["EGFR",   1:30] + 300
expr_matrix["PIK3CA", 1:30] <- expr_matrix["PIK3CA", 1:30] + 250

# Downregulated in tumour
expr_matrix["BRCA1",  1:30] <- pmax(expr_matrix["BRCA1", 1:30] - 300, 0)
expr_matrix["PTEN",   1:30] <- pmax(expr_matrix["PTEN",  1:30] - 250, 0)
expr_matrix["CDH1",   1:30] <- pmax(expr_matrix["CDH1",  1:30] - 200, 0)

col_data <- DataFrame(clinical)
col_data$disease_state <- rep(c("Tumour", "Normal"), each = 30)

# Create annotation table
anno <- select(
  org.Hs.eg.db,
  keys = gene_names,
  columns = c("ENTREZID","GENENAME", "CHR"),
  keytype = "SYMBOL"
)

row_data <- DataFrame(
  gene = gene_names,
  full_name = anno$GENENAME,
  chromosome = anno$CHR
)

# Build the summarized experiment object
se = SummarizedExperiment(
  assays = list(counts = expr_matrix),
  colData = col_data,
  rowData = row_data
)
assay(se)
colData(se)
rowData(se)

# Add a second assay of log2 transformed counts 
assay(se, "log2counts") <- log2(assay(se, "counts") + 1)

# Subset for tumour patients
tumour_se <- se[, colData(se)$disease_state == "Tumour"]

# Mutation data 
n_muts <- 150

mutation_gene_pool <- c(
  rep("TP53",   30), rep("KRAS",   22), rep("BRCA1",  18),
  rep("EGFR",   15), rep("MYC",    12), rep("PTEN",   10),
  rep("RB1",     8), rep("APC",     8), rep("BRAF",    7),
  rep("PIK3CA",  6), rep("CDKN2A",  5), rep("ERBB2",   4),
  rep("MDM2",    3), rep("VEGFA",   2)
)

mutations <- data.frame(
  mutation_id   = paste0("MUT", 1:n_muts),
  patient_id    = sample(clinical$patient_id, n_muts, replace=TRUE),
  gene          = sample(mutation_gene_pool, n_muts, replace=FALSE),
  mutation_type = sample(c("SNP","Indel","CNV"),
                         n_muts, replace=TRUE,
                         prob=c(0.60, 0.25, 0.15)),
  chromosome    = sample(1:22, n_muts, replace=TRUE),
  vaf           = round(rbeta(n_muts, 2, 3), 3),
  stringsAsFactors = FALSE
)

# Check for NA values
colSums(is.na(clinical))
colSums(is.na(mutations))

# Join clinical metadata onto the mutations table so each mutation
# row also has the patient's age, cancer type, stage, and group
mutations_metadata <- left_join(mutations, clinical, by = "patient_id")

# Generate each patient's mutation burden, add to the clinical metadata and
# replace NAs with zero, if any
mutation_burden <- mutations_metadata |> 
  summarise(
    mut_burden = n(),
    .by = patient_id
  )
clinical <- left_join(clinical, mutation_burden, by = "patient_id")

clinical <- clinical |> 
  mutate(
    mut_burden = replace_na(mut_burden, 0)
  )

# Convert the expression matrix to a dataframe in order to pivot to the
# long format and then join to the clinical metadata
expr_df <- as.data.frame(expr_matrix)
expr_df <- expr_df |> 
  rownames_to_column(var = "genes")
expr_df_long <- pivot_longer(
  expr_df,
  cols = starts_with("P"),
  names_to = "patient_id",
  values_to = "expression"
)
clinical_expr_long <- left_join(expr_df_long, clinical, by = "patient_id")

# Join the long expression data to make a full annotation table
full_anno_table <- left_join(expr_df_long, anno, by = c("genes" = "SYMBOL"))
print(full_anno_table)

# Separate the long expression data based on group to check for normality
tumour_expr <- clinical_expr_long |> 
  filter(group == "Tumour") |> 
  mutate(
    log2_expr = log2(expression + 1)
  )
normal_expr <- clinical_expr_long |> 
  filter(group == "Normal") |> 
  mutate(
    log2_expr = log2(expression + 1)
  )
shapiro.test(tumour_expr$log2_expr)
shapiro.test(normal_expr$log2_expr)

clinical_expr_long <- clinical_expr_long|> 
  mutate(
    log2_expr = log2(expression + 1)
  )

# Using Wilcoxon test as data is not normally distributed
wilcox_p.values <- clinical_expr_long %>%
  group_by(genes) %>%
  summarise(
    p_value = wilcox.test(expression ~ group)$p.value
  )
padj <- p.adjust(wilcox_p.values$p_value, method = "BH")

# Calculate each gene's tumour and normal expression to generate log2FC
mean_tumour <- tumour_expr |> 
  summarise(
    mean_tumour = round(mean(expression)),
    .by = genes
  )
mean_normal <- normal_expr |> 
  summarise(
    mean_tumour = round(mean(expression)),
    .by = genes
  )
log2FC <- log2(mean_tumour$mean_tumour/mean_normal$mean_tumour)

#Save the important results to a dataframe
results_df <- data.frame(
  gene = mean_normal$genes,
  log2FC = round(log2FC, 3),
  pvalue = wilcox_p.values$p_value,
  padj = padj
)
results_df

#Filter the significant genes
sig_results_df <- results_df |> 
  filter(padj < 0.05)

mean_mut_burden <- mean(clinical$mut_burden)
mean_df <- data.frame(
  x = 50,  
  y = mean_mut_burden + 0.2,
  label = paste("Mean =", round(mean_mut_burden, 2))
)

# Bar chart of mutation  burden per patient
p1 <- ggplot(clinical, aes(x = fct_reorder(patient_id, -mut_burden), y = mut_burden, fill = cancer_type)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = mean_mut_burden,
             linetype = "dashed",
             colour = "red") +
  geom_text(data = mean_df, 
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE) +
  labs(
    title = "Mutation Burden in Different Patients",
    x = "Patients",
    y = "Mutation Burden"
  ) +
  facet_wrap(~ group) +
  theme(axis.text.x = element_blank(),    
        axis.ticks.x = element_blank()) 
ggsave(filename = "fig1_mutation_burden.png", plot = p1,
       width = 9, height = 5, dpi = 300)

# Add a status to the results dataframe based on log2FC value
results_df <- results_df |> 
  mutate(
    status = case_when(
      padj < 0.05 & log2FC > 0 ~ "Up",
      padj < 0.05 & log2FC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )

# Volcano plot showing log2FC and -log10(pvalue)
p2 <- ggplot(results_df, aes(log2FC, -log10(pvalue), colour = status)) +
  geom_point() +
  geom_text(data = sig_results_df, aes(x = log2FC, y = -log10(pvalue), label = gene),
             inherit.aes = FALSE,
             vjust = -1) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "red") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  labs(
    title = "Volcano Plot Showing DEGs with Significant P-Values",
    x = "Fold Change",
    y = "P-Values (Negative Logarithm)"
  )
ggsave(filename = "fig2_volcano.png", plot = p2,
       width = 7, height = 7, dpi = 300)

# Principal Component Analysis
log2_expr_matrix <- log2(expr_matrix + 1)
pca <- prcomp(t(log2_expr_matrix), center = TRUE, scale. = TRUE)
variance_exp <- round(((pca$sdev)^2 / sum((pca$sdev) ^ 2) * 100), 2)
pca$var_exp <- variance_exp


pca_df <- cbind(clinical, pca$x)

# Principal Components Coloured Based on Group
p3 <- ggplot(pca_df, aes(PC1, PC2, colour = group)) +
  geom_point() +
  labs(
    x = paste0("PC1 Variance explained: ", variance_exp[1], "%"),
    y = paste0("PC2 Variance explained: ", variance_exp[2], "%")
  )

# PCs Coloured Based on Cancer Type
p4 <- ggplot(pca_df, aes(PC1, PC2, colour = cancer_type)) +
  geom_point() +
  labs(
    x = paste0("PC1 Variance explained: ", variance_exp[1], "%"),
    y = paste0("PC2 Variance explained: ", variance_exp[2], "%")
  )

scree_df <- data.frame(
  PC       = 1:length(variance_exp),
  variance = variance_exp
)

#Scree plots of first 10 PCs
scree_df <- scree_df[1:10, ]
p6 <- ggplot(scree_df, aes(x = PC, y = variance)) +
  geom_bar(stat = "identity", fill = "blue") + 
  geom_line(colour = "red") +
  geom_point(colour = "red") +
  labs(
    title = "Variance Explained by First 10 PCs",
    x = "PCs",
    y = "Variance Explained (%)"
  ) +
  scale_x_continuous(breaks = seq(0, 10, by = 1))

p5 <- (p3 + p4) +
    plot_layout(
    axis_titles = "collect",
    guides = "collect"
  )
        
p7 <- p5 / p6 + plot_annotation(
  title = "PC1 vs PC2 based on Group and Cancer Type"
)
ggsave(filename = "fig3_pca.png", plot = p7,
       width = 9, height = 7, dpi = 300)

# Survival analysis
# Create survival object and fit KM curves by group and cancer type
surv_obj <- Surv(clinical$survival_time, clinical$status)
km_fit_group <- survfit(surv_obj ~ group, data =  clinical)
summary(km_fit_group)

# KM plot showing  pvalue, confidence intervals and risk table
p8 <- ggsurvplot(km_fit_group,
           data = clinical,
           title = "Change in Probability of Survival with Time",
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE)
g_to_save <- survminer:::.build_ggsurvplot(p8)
ggsave("fig4_survival.png", g_to_save, width = 10, height = 7, dpi = 300)

survdiff(surv_obj ~ group, data=clinical)

km_fit_cancer <- survfit(surv_obj ~ cancer_type, data =  clinical)
summary(km_fit_cancer)

ggsurvplot(km_fit_cancer,
           data = clinical,
           pval = TRUE,
           conf.int = TRUE,
           risk.table = TRUE)
survdiff(surv_obj ~ cancer_type, data=clinical)

# Fit a cox proportional hazards model
cox_model <- coxph(surv_obj ~ group + age + stage + cancer_type,
                   data = clinical)
summary(cox_model)

# Extract hazard ratios and p_values
hr_table <- summary(cox_model)$coefficients
hr_values <- exp(coef(cox_model))  
pvals <- hr_table[, "Pr(>|z|)"]  # P-values

# Run cox.zph on the cox model
ph_test <- cox.zph(cox_model)
ggcoxzph(ph_test)

# Correlation Analysis
cor_matrix <- cor(t(expr_matrix[, 1:30]), method = "spearman")
# Open PNG device to save 
png("fig5_correlation.png", width = 1200, height = 1200)
corrplot(cor_matrix, 
         method = "color", 
         type = "upper",
         addCoef.col = "black",
         number.cex = 1.5,
         tl.cex = 1.5)
title(main = "Cancer Genes Correlation Across 30 Tumour Samples",
               cex.main = 3.0, line = 1.5)
dev.off()

# Find strongest correlated gene pair
spearman_test <- cor.test(expr_matrix["EGFR", 1:30], 
                          expr_matrix["VEGFA", 1:30], 
                          method = "spearman")
print(spearman_test)

# Genomic ranges
# Create GRanges object
cancer_genes_gr <- GRanges(
  seqnames = c("chr17","chr17","chr12","chr7","chr8",
               "chr10","chr13","chr5","chr7","chr3",
               "chr9","chr17","chr12","chr6","chr16"),
  ranges   = IRanges(
    start = c(7668402,  41196312, 25398284, 55086725,
              127735434, 89692905, 48303747, 112707498,
              140719327, 179148114, 21967751, 37844393,
              114705880, 43770209,  68771195),
    end   = c(7687538,  41277500, 25403854, 55275031,
              127742951, 89728532, 49056122, 112846239,
              140924929, 179240093, 21994330, 37886679,
              114803509, 43786478,  68862835)
  ),
  strand = c("-","+","-","+","+",
             "-","-","+","+","+",
             "+","+","-","+","+"),
  gene   = gene_names
)

# Somatic variants
somatic_variants <- data.frame(
  CHROM  = c("chr17","chr17","chr12","chr7","chr8",
             "chr10","chr17","chr12","chr7","chr17",
             "chr3", "chr9", "chr17","chr6", "chr16"),
  POS    = c(7674220,  41200000, 25400000, 55200000,
             127738000, 89700000, 7680000,  25401000,
             55150000,  41250000, 179150000,21970000,
             37845000,  43772000, 68773000),
  REF    = c("A","C","G","T","A","C","G","T","A","C",
             "G","T","A","C","G"),
  ALT    = c("G","T","A","C","G","T","A","G","T","A",
             "T","C","G","T","A"),
  QUAL   = c(250,180,320,195,410,288,175,290,340,220,
             310,185,260,195,305),
  FILTER = c("PASS","PASS","PASS","PASS","PASS",
             "PASS","PASS","PASS","PASS","PASS",
             "PASS","PASS","PASS","PASS","PASS"),
  AF     = c(0.35,0.48,0.29,0.42,0.52,
             0.38,0.22,0.61,0.44,0.33,
             0.55,0.41,0.37,0.29,0.48),
  DP     = c(120,195,150,180,200,
             210,180,175,230,190,
             165,220,140,195,175),
  stringsAsFactors = FALSE
)

# Convert somatic variants to GRanges object but give an end position first
somatic_variants$END <- somatic_variants$POS
somatic_gr <- makeGRangesFromDataFrame(somatic_variants, 
                               keep.extra.columns=TRUE,
                               start.field="POS", 
                               end.field="END")

# Find overlapping genes
hits <- findOverlaps(somatic_gr, cancer_genes_gr)

# Result table for overlapping genes
overlaps_df <- data.frame(
  gene_name = cancer_genes_gr$gene[subjectHits(hits)],
  variant_pos = start(somatic_gr)[queryHits(hits)],
  ref = somatic_gr$REF[queryHits(hits)],
  alt = somatic_gr$ALT[queryHits(hits)],
  af = somatic_gr$AF[queryHits(hits)],
  dp = somatic_gr$DP[queryHits(hits)]
)
overlaps_df

# Which gene has most variants?
table(overlaps_df$gene)

# Export overlapping variants as BED file
overlapping_variants <- somatic_gr[queryHits(hits)]

# Distance of each variant to the nearest cancer gene
distance_to_nearest <- distanceToNearest(somatic_gr,cancer_genes_gr)

# Write overlapping variants to a bed file
bed_df <- data.frame(
  chrom  = as.character(seqnames(overlapping_variants)),
  start  = start(overlapping_variants) - 1,  # 0-based
  end    = end(overlapping_variants),
  name   = cancer_genes_gr$gene[subjectHits(hits)],
  score  = overlapping_variants$AF,
  strand = "."
)

write.table(bed_df, "somatic_overlaps.bed",
            sep       = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote     = FALSE)

# Logistic Regression
# Create a binary group to use to fit model
clinical <- clinical |> 
  mutate(
    group_binary = case_when(group == "Tumour" ~ 1,
                              TRUE ~ 0)
  )
glm_fit <- glm(group_binary ~ age + stage + cancer_type,
               family = "binomial",
               data = clinical)
summary(glm_fit)

# Extract odds ratio and McFadden's pseudo R
odds_ratio <- exp(coef(glm_fit))
mcfadden_r2 <- 1 - (glm_fit$deviance / glm_fit$null.deviance)
odds_ratio

# Derive predicted probabilities and plot against age
clinical$predicted_prob <- predict(glm_fit, type = "response")
ggplot(clinical, aes(x=age, y=predicted_prob, colour = stage)) +
  geom_line(colour="red", linewidth=1) +
  geom_point(aes(y=group_binary), alpha=0.3) +
  labs(title="How Age Affects Disease Status",
       x="Age", y="Predicted Probability") +
  theme_bw()

# Final Summary Table
mutations_per_gene <- mutations |> 
  summarise(
    n_mutations = n(),
    .by = gene
  )
results_df <- left_join(results_df, mutations_per_gene, by = "gene")
results_df <- left_join(results_df, anno, by = c("gene" = "SYMBOL"))
results_df <- results_df |> 
  dplyr::select(gene, GENENAME, CHR, ENTREZID, log2FC, pvalue, padj, status, n_mutations)


