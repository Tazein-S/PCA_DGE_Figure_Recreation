library(edgeR)

# Load counts and define sample groups
counts <- read.table("merged_counts_matrix.txt", header=TRUE, row.names=1, sep="\t")
samples <- colnames(counts)
group <- ifelse(grepl("^CTC", samples), "CTC", 
                ifelse(grepl("^BLOOD", samples), "PB", "PT"))
coldata <- data.frame(sample = samples, group = group)

# Subset only CTC and PT samples
dge_counts <- counts[, coldata$group %in% c("CTC","PT")]
dge_coldata <- coldata[coldata$group %in% c("CTC","PT"), , drop=FALSE]

# edgeR differential expression
y <- DGEList(counts = dge_counts, group = dge_coldata$group)

# keep genes with at least 10 CPM in 2 samples
keep <- rowSums(cpm(y) > 10) >= 2
y <- y[keep, ]
y$samples$lib.size <- colSums(y$counts)

# normalize
y <- calcNormFactors(y)

# design matrix
design <- model.matrix(~ group, data = y$samples)

# dispersion
y <- estimateDisp(y, design)
# plotBCV(y)  # optional

# Fit model and test
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = 2)

# Extract all DEGs
deg_table <- topTags(lrt, n = Inf)$table
deg_table$ensembl_id <- rownames(deg_table)

# Significant DEGs
deg_sig <- deg_table[deg_table$FDR < 0.05, ]

### ==========================================================
### 🔥 Map Ensembl IDs → Gene Symbols
### ==========================================================

# Load GTF
gtf <- read.table("data/gencode.v25.annotation.gtf", 
                  sep="\t", quote="", comment.char="#", 
                  stringsAsFactors=FALSE)
gtf_genes <- gtf[gtf$V3 == "gene", ]

# Helper to extract attributes
extract_attr <- function(attr, key) {
  m <- regmatches(attr, regexpr(paste0(key, ' "([^"]+)"'), attr))
  if (length(m) > 0) {
    gsub(paste0(key, ' "([^"]+)"'), "\\1", m)
  } else {
    NA
  }
}

gtf_genes$gene_id_clean <- gsub("\\..*", "", 
                                sapply(gtf_genes$V9, extract_attr, key="gene_id"))
gtf_genes$gene_name <- sapply(gtf_genes$V9, extract_attr, key="gene_name")
gtf_genes <- gtf_genes[!duplicated(gtf_genes$gene_id_clean), ]

# Prepare deg_table for merge
res_df <- deg_table
res_df$ensembl_id_clean <- gsub("\\..*", "", res_df$ensembl_id)

# Merge gene symbols
res_df <- merge(
  res_df,
  gtf_genes[, c("gene_id_clean", "gene_name")],
  by.x = "ensembl_id_clean",
  by.y = "gene_id_clean",
  all.x = TRUE
)

# Use gene_name or fallback to Ensembl ID
res_df$gene_symbol <- ifelse(
  !is.na(res_df$gene_name) & res_df$gene_name != "",
  res_df$gene_name,
  res_df$ensembl_id
)

# Reorder columns
col_order <- c("gene_symbol", "ensembl_id", "ensembl_id_clean", "gene_name",
               setdiff(names(res_df), c("gene_symbol", "ensembl_id",
                                        "ensembl_id_clean", "gene_name")))
res_df <- res_df[, col_order]

# Optional: set rownames to gene symbols
#rownames(res_df) <- res_df$gene_symbol

### ==========================================================
### 🔥 Output results
### ==========================================================

write.table(res_df,
            file = "edgeR_all_DEGs_CTC_vs_PT_geneSymbols.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

deg_sig_full <- res_df[res_df$FDR < 0.05, ]
write.table(deg_sig_full,
            file = "edgeR_significant_DEGs_CTC_vs_PT_geneSymbols.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Quick summary
summary(decideTests(lrt))
