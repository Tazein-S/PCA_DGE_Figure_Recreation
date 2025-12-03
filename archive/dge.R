library(DESeq2)
library(dplyr)

# ===== Load counts and define sample groups =====
counts <- read.table("merged_counts_matrix.txt", header=TRUE, row.names=1, sep="\t")
samples <- colnames(counts)
group <- ifelse(grepl("^CTC", samples), "CTC", 
                ifelse(grepl("^BLOOD", samples), "PB", "PT"))
coldata <- data.frame(sample = samples, group = group)

# Subset only CTC and PT samples
dge_counts <- counts[, coldata$group %in% c("CTC","PT")]
dge_coldata <- coldata[coldata$group %in% c("CTC","PT"), , drop=FALSE]

# ===== DESeq2 differential expression =====
dds <- DESeqDataSetFromMatrix(countData = dge_counts, 
                              colData = dge_coldata, 
                              design = ~ group)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "CTC", "PT"))
res_df <- as.data.frame(res)
res_df$ensembl_id <- rownames(res_df)

# ===== Map Ensembl IDs to gene symbols =====
gtf <- read.table("data/gencode.v25.annotation.gtf", 
                  sep="\t", 
                  quote="", 
                  comment.char="#", 
                  stringsAsFactors = FALSE)
gtf_genes <- gtf[gtf$V3 == "gene", ]

extract_attr <- function(attr, key) {
  m <- regmatches(attr, regexpr(paste0(key, ' "([^"]+)"'), attr))
  if(length(m) > 0) gsub(paste0(key, ' "([^"]+)"'), "\\1", m) else NA
}

gtf_genes$gene_id_clean <- gsub("\\..*", "", 
                                sapply(gtf_genes$V9, extract_attr, key="gene_id"))
gtf_genes$gene_name <- sapply(gtf_genes$V9, extract_attr, key="gene_name")
gtf_genes <- gtf_genes[!duplicated(gtf_genes$gene_id_clean), ]

res_df$ensembl_id_clean <- gsub("\\..*", "", res_df$ensembl_id)

# Merge to add gene names
res_df <- merge(res_df, 
                gtf_genes[, c("gene_id_clean","gene_name")], 
                by.x = "ensembl_id_clean", 
                by.y = "gene_id_clean", 
                all.x = TRUE)

# Create gene_symbol
res_df$gene_symbol <- ifelse(!is.na(res_df$gene_name) & res_df$gene_name != "", 
                             res_df$gene_name, 
                             res_df$ensembl_id)

# Reorder columns
col_order <- c("ensembl_id", "ensembl_id_clean", "gene_symbol", "gene_name", 
               setdiff(names(res_df), c("ensembl_id", "ensembl_id_clean", 
                                        "gene_symbol", "gene_name")))
res_df <- res_df[, col_order]

# ===== Prepare for plotting =====
res_df <- res_df[!is.na(res_df$padj), ]
res_df$log10padj <- -log10(res_df$padj)

padj_threshold <- 0.05
lfc_threshold <- 1

res_df$regulation <- "NS"
res_df$regulation[res_df$padj < padj_threshold & 
                    res_df$log2FoldChange > lfc_threshold] <- "Up in CTC"
res_df$regulation[res_df$padj < padj_threshold & 
                    res_df$log2FoldChange < -lfc_threshold] <- "Up in PT"
res_df$regulation <- factor(res_df$regulation, 
                            levels = c("Up in CTC", "NS", "Up in PT"))

# ===== Save results for Python =====
write.table(res_df, 
            "DESeq2_results_for_python.tsv", 
            sep="\t", 
            quote=FALSE, 
            row.names=FALSE)

write.table(data.frame(top_genes_all = top_genes_all), 
            "top_genes_to_label.txt", 
            sep="\t", 
            quote=FALSE, 
            row.names=FALSE, 
            col.names=FALSE)