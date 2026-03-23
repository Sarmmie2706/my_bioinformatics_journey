library(tidyverse)

cancer_db <- data.frame(
  name       = c("TP53","BRCA1","KRAS","EGFR","MYC","PTEN",
                 "RB1","APC","BRAF","PIK3CA","CDKN2A","ERBB2"),
  full_name  = c("Tumor Protein p53","Breast Cancer 1",
                 "KRAS Proto-oncogene","Epidermal Growth Factor Receptor",
                 "MYC Proto-oncogene","Phosphatase and Tensin Homolog",
                 "Retinoblastoma 1","APC Regulator of WNT Signaling",
                 "B-Raf Proto-oncogene","PI3-Kinase Catalytic Alpha",
                 "Cyclin Dependent Kinase Inhibitor 2A","Erb-B2 Receptor Tyrosine Kinase 2"),
  length     = c(1182,5592,567,3633,1600,1212,
                 2739,8538,2301,3207,1236,3770),
  chromosome = c(17,17,12,7,8,10,
                 13,5,7,3,9,17),
  gc_content = c(41.5,42.1,54.3,49.8,58.2,44.6,
                 42.8,47.3,51.2,52.6,55.1,48.9),
  type       = c("TSG","TSG","Oncogene","Oncogene","Oncogene","TSG",
                 "TSG","TSG","Oncogene","Oncogene","TSG","Oncogene"),
  stringsAsFactors = FALSE
)
cancer_db

lookup_gene <- function(gene) {
  if (!(gene %in% cancer_db$name)) {
    print("Gene not present in database")
  } else {
    i <- which(gene == cancer_db$name)
    cat("======", cancer_db$name[i], "======\n")
    cat(sprintf("%-13s: %s\n", "Full name: ", cancer_db$full_name[i]))
    cat(sprintf("%-13s: %d bp\n", "Length: ", cancer_db$length[i]))
    cat(sprintf("%-13s: %s\n", "Chromosome: ", cancer_db$chromosome[i]))
    cat(sprintf("%-13s: %.1f%%\n", "GC_Content: ", cancer_db$gc_content[i]))
    cat(sprintf("%-13s: %s\n", "Type: ", cancer_db$type[i]))
  }
}
lookup_gene("BRAF")
lookup_gene("NOGENE")

genes_on_chromosome <- function(chr) {
  if (!(chr %in% cancer_db$chromosome)) {
    print("Chromosome not present in database")
  } else {
    filtered_db <- cancer_db %>% 
      filter(chromosome == chr) %>% 
      arrange(desc(length))
    cat(sprintf("%d gene(s) found on chr %d\n", nrow(filtered_db), chr))
    print(filtered_db)
  }
}
genes_on_chromosome(5)
genes_on_chromosome(170)

summary_by_type <- function() {
  cancer_db %>% 
    group_by(type) %>% 
    summarise(
      count = n(),
      mean_length = mean(length),
      mean_gc_content = mean(gc_content)
    )
}
summary_by_type()
cancer_db
add_gene <- function(db, name, full_name, length, chromosome, gc_content, type) {
  if (name %in% db$name) {
    cat(sprintf("Warning: %s already exists in the database", name))
  } else {
    db <- db %>% 
      add_row(
        name = name,
        full_name = full_name,
        length = length,
        chromosome = chromosome,
        gc_content = gc_content,
        type = type
      )
    return(db)
  }
}
cancer_db <- add_gene(cancer_db, "VEGFA", "Vascular Endothelial Growth Factor A", 1472, 6, 50.3, "Oncogene")
cancer_db

classify_gene <- function(length) {
  if (length < 500) {
    return("Small")
  } else if (length >= 500 & length < 2000) {
    return("Medium")
  } else {
    return("Large")
  }
}

classify_all <- function() {
  invisible(sapply(1:nrow(cancer_db), function(i) {
   cat(sprintf("%-8s %-45s %s\n", 
      cancer_db$name[i],
      cancer_db$full_name[i],
      classify_gene(cancer_db$length[i])
    ))
  }))
}
classify_all()

write.csv(cancer_db, "cancer_genes.csv", row.names = FALSE)
cancer_db_reloaded <- read.csv("cancer_genes.csv", stringsAsFactors = FALSE)
View(cancer_db_reloaded)


