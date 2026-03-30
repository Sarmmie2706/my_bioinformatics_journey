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

# Mutation Type Distribution
type_counts <- mutations %>% 
  filter(gene != "OTHER") %>% 
  summarise(
    gene_count = n(),
    .by = gene
  )
ggplot(type_counts, aes(x = gene_count, y = reorder(gene, gene_count))) +
geom_bar(stat = "identity", aes(fill = gene_count)) +
scale_fill_gradient(low="lightblue", high="darkblue") +
geom_text(aes(label = gene_count), hjust = 1, colour = "white") +
labs(
  title = "Most Mutated Genes",
  x = "Gene Count",
  y = "Genes"
) +
theme_classic()
ggsave(filename = "fig1_top_genes.png", dpi = 300, width = 9, height = 5)

# Mutation Type Distribution
mutation_counts <- mutations %>% 
  count(mutation_type, name = "counts") %>% 
  mutate(
    pct = (counts / sum(counts)) * 100
  )

p1 <- ggplot(mutation_counts, aes(mutation_type, counts, fill = mutation_type)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none") +
  labs(
    x = "Mutation Type",
    y = "Count"
  ) +
  coord_flip()

p2 <- ggplot(mutation_counts, aes(x="", y=pct, fill=mutation_type)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y") +
  theme_void()
p1 + p2 + plot_annotation(
  title = "Most Common Mutation Types and their Distributions"
  ) +
  plot_layout(guides = "collect")
ggsave(filename = "fig2_mutation_types.png", dpi = 300, height = 5, width = 8)

# VAF distribution by mutation type
ggplot(mutations, aes(x = mutation_type, y = vaf)) +
  geom_violin(aes(fill = mutation_type)) +
  geom_boxplot(width = 0.1) +
  geom_jitter(width = 0.1, alpha = 0.3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "red") +
  labs(
    title = "VAF Distribution by Mutation Type",
    x = "Mutation Type",
    y = "VAF"
  ) +
  theme_classic()
ggsave(filename = "fig3_vaf_distribution.png", dpi = 300, height = 5, width = 8)

# Chromosome Distribution
ggplot(mutations, aes(x = factor(chromosome, levels = 1:22), fill = cancer_type)) +
  geom_bar(position = "fill") +
  labs(
    title = "Chromosome Distribution of Different Cancer Types",
    x = "Chromosomes",
    y = "Proportion"
  ) +
  theme_classic()
ggsave(filename = "fig4_chromosome_distribution.png", dpi = 300, height = 5, width = 9)

########## STATISTICAL TESTS #############
# Do SNPs have higher VAF than Indels?
SNPs <- mutations %>% 
  filter(mutation_type == "SNP")
shapiro.test(SNPs$vaf)

indels <- mutations %>% 
  filter(mutation_type == "Indel")
shapiro.test(indels$vaf)

ggplot(data.frame(x = indels$vaf), aes(sample = x)) +
  stat_qq() +
  stat_qq_line(colour = "red")

wilcox_results <- wilcox.test(SNPs$vaf, indels$vaf)
# Using the wilcox test, as both data are not normally distributed,
# the test shows a W value of 2,366 and a p-value of 0.8572,
# meaning that any observed differences in the mean is most likely sue to chance,
# and hence, we fail to reject the null hypothesis.

# Is mutation type independent of cancer type?
mutation_cancer_table <- table(mutations$mutation_type, mutations$cancer_type)
chisq_results <- chisq.test(mutation_cancer_table)
chisq_results
chisq_results$residuals
chisq_results$expected
# At an X-squared value of 2.0911, df of 4 and p-value of 0.719,
# There is no significant association between the variables.

# Does VAF differ across all three mutation groups?
aov_results <- aov(vaf ~ mutation_type, data = mutations)
TukeyHSD(aov_results)
anova_pvalue <- summary(aov_results)[[1]][["Pr(>F)"]][1]
# the aov results showed a p-value of 0.0276. Tukey's test showed 
# that the pair with a significant difference in means is SNP-CNV.

# Multiple Testing Correction
all_pvals <- c(wilcox_results$p.value, chisq_results$p.value, anova_pvalue)
p.adjust(all_pvals, method="BH")

# After correction, the anova p.value increased from 0.0275 to 0.0827, now making it 
# non-significant. The conclusions of the other two tests remained non-significant, as 
# the p_values were still > 0.05

sample_data <- mutations %>% 
  summarise(
    mut_count = n(),
    .by = sample_id
  ) %>% 
  mutate(
    age = round(mut_count * 4 + rnorm(40, 30, 5), 0)
  )

lm_results <- lm(mut_count ~ age, data = sample_data)
lm_rsquared <- round(summary(lm_results)$r.squared * 100, 2)
summary(lm_results)$coefficients

ggplot(sample_data, aes(x = age, y = mut_count)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(
    title = "How Much Can Age Predict Mutation Burden?",
    subtitle = paste0("This model explains ", lm_rsquared, "% of the variation."),
    x = "Age",
    y = "Mutation Count"
  )
ggsave(filename = "fig5_mutation_vs_age.png", dpi = 300, height = 5, width = 8)
# Age significantly predicts mutation burden (β = 0.2205, p < 0.001),
# and explains a large proportion of the variance (R² = 0.8209)

par(mfrow = c(2,2))
plot(lm_results)

mutations
mutation_matrix <- mutations %>% 
  filter(gene != "OTHER") %>% 
  count(sample_id, gene, name="mut_count") %>% 
  pivot_wider(names_from=gene, 
              values_from=mut_count, 
              values_fill=0) %>% 
  column_to_rownames("sample_id")

pca <- prcomp(mutation_matrix, scale. = TRUE)
pca$x
pca$sdev
pca$rotation

var_explained <- (pca$sdev^2) / sum(pca$sdev^2) * 100

# Scree Plot
scree_df <- data.frame(
  PC      = 1:length(var_explained),
  variance = var_explained
)

p1 <- ggplot(scree_df, aes(x=PC, y=variance)) +
  geom_bar(stat="identity", fill="#2E75B6") +
  geom_line(colour="red") +
  geom_point(colour="red") +
  labs(title="Scree Plot", 
       x="Principal Component",
       y="Variance Explained (%)") +
  theme_bw()

p2 <- ggplot(data.frame(pca$x), aes(x = PC1, y = PC2)) +
  geom_point() +
  labs(
    title = "PC1 vs PC2",
    x = paste0("PC1 (", round(var_explained[1],1), "% variance)"),
    y = paste0("PC2 (", round(var_explained[2],1), "% variance)")
  ) +
  theme_classic()

p1 + p2 + plot_annotation(
  title = "Sample Mutation Profiles using PCA"
)
ggsave(
  filename = "fig6_pca_sample_mutation_profiles.png",
  dpi = 300, height = 5, width = 9
)
# Which genes contribute most to PC1?
loadings_pc1 <- pca$rotation[,"PC1"]
sort(abs(loadings_pc1), decreasing=TRUE)[1:10]
# These are the 10 genes driving the most variation on PC1







