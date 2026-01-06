#===============================================================================
# BALLGOWN GENE-LEVEL DIFFERENTIAL EXPRESSION ANALYSIS
# Comparison: TUMOR vs NORMAL
# CORRECTED VERSION - Using q-value (FDR) for significance
#===============================================================================

#-------------------------------------------------------------------------------
# 1. SETUP AND DEPENDENCIES
#-------------------------------------------------------------------------------

library(ballgown)
library(genefilter)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(ggrepel)
library(biomaRt)
library(tidyr)

# Resolve function conflicts (biomaRt/AnnotationDbi masks dplyr functions)
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

set.seed(42)

# Create output directory
dir.create("results", showWarnings = FALSE)
dir.create("results/plots", showWarnings = FALSE)

#-------------------------------------------------------------------------------
# 2. DATA LOADING
#-------------------------------------------------------------------------------

# Set working directory to your Ballgown folder
# setwd("/path/to/your/Ballgown")

# Load metadata
metadata <- read.csv("metadata.csv", header = TRUE)

cat("Metadata columns:", colnames(metadata), "\n")
cat("Number of samples:", nrow(metadata), "\n")
cat("Conditions:\n")
print(table(metadata$condition))

# Create ballgown object
bg <- ballgown(
  dataDir = ".",
  samplePattern = "WT",
  pData = metadata,
  meas = "all"
)

cat("\nLoaded samples:\n")
print(sampleNames(bg))

save(bg, file = "results/bg_raw.rda")

#-------------------------------------------------------------------------------
# 3. FILTERING LOW-EXPRESSION GENES
#-------------------------------------------------------------------------------

bg_filt <- subset(bg, "rowVars(texpr(bg)) > 1", genomesubset = TRUE)

cat("\n--- Filtering Summary ---\n")
cat("Transcripts before filtering:", nrow(texpr(bg)), "\n")
cat("Transcripts after filtering:", nrow(texpr(bg_filt)), "\n")

save(bg_filt, file = "results/bg_filtered.rda")

#-------------------------------------------------------------------------------
# 4. GENE-LEVEL DIFFERENTIAL EXPRESSION ANALYSIS
#-------------------------------------------------------------------------------

results_genes <- stattest(
  bg_filt,
  feature = "gene",
  meas = "FPKM",
  covariate = "condition",
  getFC = TRUE
)

results_genes$log2FC <- log2(results_genes$fc)
results_genes$log2FC[is.infinite(results_genes$log2FC)] <- NA

#-------------------------------------------------------------------------------
# 5. GENE ANNOTATION (Ensembl ID to Gene Symbol)
#-------------------------------------------------------------------------------

cat("\n--- Annotating Genes with biomaRt ---\n")

results_genes$ensembl_id <- gsub("\\..*", "", results_genes$id)

ensembl <- tryCatch({
  useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
}, error = function(e) {
  useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "www")
})

annotations <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description",
                 "gene_biotype", "chromosome_name", "start_position",
                 "end_position", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = unique(results_genes$ensembl_id),
  mart = ensembl
)

annotations$description <- gsub(" \\[Source:.*\\]", "", annotations$description)
annotations <- annotations[!duplicated(annotations$ensembl_gene_id), ]

results_genes <- results_genes %>%
  left_join(annotations, by = c("ensembl_id" = "ensembl_gene_id")) %>%
  dplyr::rename(gene_symbol = external_gene_name, gene_name = description,
                biotype = gene_biotype, chr = chromosome_name,
                start = start_position, end = end_position, entrez_id = entrezgene_id)

results_genes$label <- ifelse(!is.na(results_genes$gene_symbol) & results_genes$gene_symbol != "",
                              results_genes$gene_symbol, results_genes$id)

#-------------------------------------------------------------------------------
# 6. ADD SIGNIFICANCE CATEGORIES (USING Q-VALUE!)
#-------------------------------------------------------------------------------

results_genes <- results_genes %>%
  mutate(
    significance = case_when(
      is.na(qval) ~ "Not tested",
      qval < 0.05 & abs(log2FC) >= 1 ~ "Significant (q<0.05, |log2FC|>=1)",
      qval < 0.05 ~ "Significant (q<0.05)",
      qval < 0.1 ~ "q<0.1",
      TRUE ~ "Not significant"
    ),
    direction = case_when(
      is.na(log2FC) ~ "NA",
      log2FC > 0 ~ "Upregulated in TUMOR",
      log2FC < 0 ~ "Downregulated in TUMOR",
      TRUE ~ "No change"
    )
  )

results_genes <- arrange(results_genes, qval)

results_genes <- results_genes %>%
  dplyr::select(id, ensembl_id, gene_symbol, gene_name, biotype, chr, start, end,
                entrez_id, fc, log2FC, pval, qval, significance, direction, label)

cat("\n--- DE Summary (q-value) ---\n")
cat("Total genes:", nrow(results_genes), "\n")
cat("Significant (q < 0.05):", sum(results_genes$qval < 0.05, na.rm = TRUE), "\n")
cat("Stringent (q < 0.05 & |log2FC| >= 1):", sum(results_genes$qval < 0.05 & abs(results_genes$log2FC) >= 1, na.rm = TRUE), "\n")

#-------------------------------------------------------------------------------
# 7. EXTRACT GENE EXPRESSION MATRIX
#-------------------------------------------------------------------------------

gene_expr <- as.data.frame(gexpr(bg_filt))
colnames(gene_expr) <- gsub("FPKM\\.", "", colnames(gene_expr))
gene_expr$ensembl_id <- gsub("\\..*", "", rownames(gene_expr))
gene_expr <- gene_expr %>%
  left_join(results_genes %>% dplyr::select(ensembl_id, gene_symbol) %>% distinct(), by = "ensembl_id")
gene_expr$label <- ifelse(!is.na(gene_expr$gene_symbol) & gene_expr$gene_symbol != "",
                          gene_expr$gene_symbol, gene_expr$ensembl_id)
sample_cols <- setdiff(colnames(gene_expr), c("ensembl_id", "gene_symbol", "label"))
gene_expr_log2 <- gene_expr
gene_expr_log2[, sample_cols] <- log2(gene_expr[, sample_cols] + 1)

#-------------------------------------------------------------------------------
# 8. SAVE RESULTS TABLES
#-------------------------------------------------------------------------------

write.csv(results_genes, "results/all_genes_DE_results.csv", row.names = FALSE)
write.csv(results_genes %>% dplyr::filter(qval < 0.05), "results/significant_genes_q05.csv", row.names = FALSE)
write.csv(results_genes %>% dplyr::filter(qval < 0.05 & abs(log2FC) >= 1), "results/significant_genes_stringent.csv", row.names = FALSE)
write.csv(results_genes %>% dplyr::filter(qval < 0.05 & log2FC > 0), "results/significant_upregulated.csv", row.names = FALSE)
write.csv(results_genes %>% dplyr::filter(qval < 0.05 & log2FC < 0), "results/significant_downregulated.csv", row.names = FALSE)
write.csv(results_genes %>% dplyr::filter(qval < 0.05 & log2FC >= 1), "results/significant_upregulated_stringent.csv", row.names = FALSE)
write.csv(results_genes %>% dplyr::filter(qval < 0.05 & log2FC <= -1), "results/significant_downregulated_stringent.csv", row.names = FALSE)
write.csv(gene_expr, "results/gene_expression_FPKM.csv", row.names = FALSE)
write.csv(gene_expr_log2, "results/gene_expression_log2FPKM.csv", row.names = FALSE)

#-------------------------------------------------------------------------------
# 9. VISUALIZATIONS
#-------------------------------------------------------------------------------

# Volcano plot (q-value)
volcano_plot <- EnhancedVolcano(results_genes, lab = results_genes$label, x = "log2FC", y = "qval",
                                 pCutoff = 0.05, FCcutoff = 1, title = "TUMOR vs NORMAL",
                                 subtitle = "Gene-level DE (FDR-adjusted)", pointSize = 2.0, labSize = 3.5,
                                 col = c("grey30", "forestgreen", "royalblue", "red2"), colAlpha = 0.6,
                                 drawConnectors = TRUE, max.overlaps = 25, ylab = bquote(~-Log[10] ~ italic(q-value)))
ggsave("results/plots/volcano_plot_qvalue.png", volcano_plot, width = 14, height = 12, dpi = 300)
ggsave("results/plots/volcano_plot_qvalue.pdf", volcano_plot, width = 14, height = 12)

# Additional plots... (MA, heatmap, PCA, etc.)
# See full script for complete visualization code

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to results/ directory\n")
save.image("results/ballgown_analysis.RData")
