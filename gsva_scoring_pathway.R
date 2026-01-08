#===============================================================================
# GSVA SAMPLE-LEVEL PATHWAY ACTIVITY SCORING
# Calculates pathway enrichment scores per sample
# Creates heatmaps showing pathway activity across samples
#===============================================================================

#-------------------------------------------------------------------------------
# 1. INSTALL AND LOAD PACKAGES
#-------------------------------------------------------------------------------

# Install packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages_bioc <- c("GSVA", "msigdbr", "limma", "org.Hs.eg.db")
packages_cran <- c("pheatmap", "RColorBrewer", "dplyr", "tidyr", "ggplot2")

for (pkg in packages_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Load libraries
library(GSVA)
library(msigdbr)
library(limma)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(ggplot2)

#-------------------------------------------------------------------------------
# 2. CONFIGURATION
#-------------------------------------------------------------------------------

# Set working directory (modify as needed)
# setwd("/path/to/your/Ballgown")

# Create output directory
dir.create("results/gsva", showWarnings = FALSE, recursive = TRUE)
dir.create("results/gsva/plots", showWarnings = FALSE)

# Analysis parameters
GSVA_METHOD <- "gsva"  # Options: "gsva", "ssgsea", "zscore", "plage"
MIN_GENES_PER_PATHWAY <- 10
MAX_GENES_PER_PATHWAY <- 500

#-------------------------------------------------------------------------------
# 3. LOAD EXPRESSION DATA
#-------------------------------------------------------------------------------

cat("=== Loading Expression Data ===\n\n")

# Option A: Load from previous Ballgown analysis
if (file.exists("results/ballgown_analysis.RData")) {
  load("results/ballgown_analysis.RData")
  cat("Loaded data from ballgown_analysis.RData\n")
} else if (file.exists("results/gene_expression_log2FPKM.csv")) {
  # Option B: Load from CSV
  gene_expr_log2 <- read.csv("results/gene_expression_log2FPKM.csv")
  metadata <- read.csv("metadata.csv")
  cat("Loaded data from CSV files\n")
} else {
  stop("No expression data found. Run ballgown_corrected.R first.")
}

# Get sample columns (exclude annotation columns)
sample_cols <- setdiff(colnames(gene_expr_log2), c("ensembl_id", "gene_symbol", "label"))

cat("Samples:", length(sample_cols), "\n")
cat("Genes:", nrow(gene_expr_log2), "\n\n")

#-------------------------------------------------------------------------------
# 4. PREPARE EXPRESSION MATRIX
#-------------------------------------------------------------------------------

cat("=== Preparing Expression Matrix ===\n\n")

# Create expression matrix with gene symbols as rownames
expr_matrix <- as.matrix(gene_expr_log2[, sample_cols])

# Use gene symbols as row names
if ("gene_symbol" %in% colnames(gene_expr_log2)) {
  rownames(expr_matrix) <- gene_expr_log2$gene_symbol
} else if ("label" %in% colnames(gene_expr_log2)) {
  rownames(expr_matrix) <- gene_expr_log2$label
}

# Remove rows with NA or empty gene symbols
valid_rows <- !is.na(rownames(expr_matrix)) & rownames(expr_matrix) != ""
expr_matrix <- expr_matrix[valid_rows, ]

# Remove duplicate gene symbols (keep highest expressed)
row_means <- rowMeans(expr_matrix, na.rm = TRUE)
expr_matrix <- expr_matrix[order(row_means, decreasing = TRUE), ]
expr_matrix <- expr_matrix[!duplicated(rownames(expr_matrix)), ]

# Remove rows with zero variance
row_vars <- apply(expr_matrix, 1, var, na.rm = TRUE)
expr_matrix <- expr_matrix[row_vars > 0, ]

cat("Final matrix:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n\n")

#-------------------------------------------------------------------------------
# 5. GET PATHWAY GENE SETS FROM MSigDB
#-------------------------------------------------------------------------------

cat("=== Loading Pathway Gene Sets ===\n\n")

# Get C2 collection (canonical pathways)
c2_pathways <- msigdbr(species = "Homo sapiens", category = "C2")

# KEGG Legacy pathways (classic KEGG)
kegg_pathways <- c2_pathways %>% filter(gs_subcat == "CP:KEGG_LEGACY")
kegg_list <- split(kegg_pathways$gene_symbol, kegg_pathways$gs_name)
cat("KEGG pathways:", length(kegg_list), "\n")

# Hallmark gene sets (highly curated, recommended)
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)
cat("Hallmark gene sets:", length(hallmark_list), "\n")

# GO Biological Process
c5_pathways <- msigdbr(species = "Homo sapiens", category = "C5")
go_bp <- c5_pathways %>% filter(gs_subcat == "GO:BP")
go_bp_list <- split(go_bp$gene_symbol, go_bp$gs_name)
cat("GO BP terms:", length(go_bp_list), "\n")

# Reactome pathways
reactome <- c2_pathways %>% filter(gs_subcat == "CP:REACTOME")
reactome_list <- split(reactome$gene_symbol, reactome$gs_name)
cat("Reactome pathways:", length(reactome_list), "\n\n")

#-------------------------------------------------------------------------------
# 6. FILTER GENE SETS BY SIZE
#-------------------------------------------------------------------------------

filter_genesets <- function(gset_list, min_size = 10, max_size = 500, expr_genes) {
  # Keep only genes that are in our expression matrix
  gset_list <- lapply(gset_list, function(x) intersect(x, expr_genes))
  
  # Filter by size
  sizes <- sapply(gset_list, length)
  gset_list <- gset_list[sizes >= min_size & sizes <= max_size]
  
  return(gset_list)
}

expr_genes <- rownames(expr_matrix)

kegg_list_filtered <- filter_genesets(kegg_list, MIN_GENES_PER_PATHWAY, MAX_GENES_PER_PATHWAY, expr_genes)
hallmark_list_filtered <- filter_genesets(hallmark_list, MIN_GENES_PER_PATHWAY, MAX_GENES_PER_PATHWAY, expr_genes)
go_bp_list_filtered <- filter_genesets(go_bp_list, MIN_GENES_PER_PATHWAY, MAX_GENES_PER_PATHWAY, expr_genes)
reactome_list_filtered <- filter_genesets(reactome_list, MIN_GENES_PER_PATHWAY, MAX_GENES_PER_PATHWAY, expr_genes)

cat("Filtered gene sets:\n")
cat("  KEGG:", length(kegg_list_filtered), "\n")
cat("  Hallmark:", length(hallmark_list_filtered), "\n")
cat("  GO BP:", length(go_bp_list_filtered), "\n")
cat("  Reactome:", length(reactome_list_filtered), "\n\n")

#-------------------------------------------------------------------------------
# 7. RUN GSVA
#-------------------------------------------------------------------------------

cat("=== Running GSVA ===\n\n")

run_gsva_analysis <- function(expr_mat, gset_list, method = "gsva", name = "pathway") {
  
  if (length(gset_list) == 0) {
    cat("No gene sets to analyze for", name, "\n")
    return(NULL)
  }
  
  cat("Running GSVA for", name, "with", length(gset_list), "gene sets...\n")
  
  # GSVA parameters
  gsva_param <- gsvaParam(
    exprData = expr_mat,
    geneSets = gset_list,
    kcdf = "Gaussian",  # For log-transformed continuous data
    minSize = MIN_GENES_PER_PATHWAY,
    maxSize = MAX_GENES_PER_PATHWAY
  )
  
  # Run GSVA
  gsva_results <- gsva(gsva_param, verbose = FALSE)
  
  cat("  Completed:", nrow(gsva_results), "pathways scored\n")
  
  return(gsva_results)
}

# Run GSVA for each gene set collection
gsva_kegg <- run_gsva_analysis(expr_matrix, kegg_list_filtered, GSVA_METHOD, "KEGG")
gsva_hallmark <- run_gsva_analysis(expr_matrix, hallmark_list_filtered, GSVA_METHOD, "Hallmark")
gsva_go_bp <- run_gsva_analysis(expr_matrix, go_bp_list_filtered, GSVA_METHOD, "GO_BP")
gsva_reactome <- run_gsva_analysis(expr_matrix, reactome_list_filtered, GSVA_METHOD, "Reactome")

cat("\n")

#-------------------------------------------------------------------------------
# 8. CREATE SAMPLE ANNOTATION
#-------------------------------------------------------------------------------

# Create annotation dataframe for heatmap
sample_annotation <- data.frame(
  Condition = metadata$condition,
  row.names = metadata$ids
)

# Match sample names if needed
if (!all(colnames(expr_matrix) %in% rownames(sample_annotation))) {
  # Try matching without .ballgown suffix
  rownames(sample_annotation) <- gsub("\\.ballgown$", "", metadata$ids)
}

# Define annotation colors
ann_colors <- list(
  Condition = c(TUMOR = "firebrick", NORMAL = "steelblue")
)

# If you have Treatment information, add it
if ("treatment" %in% colnames(metadata)) {
  sample_annotation$Treatment <- metadata$treatment
  ann_colors$Treatment <- c(Treated = "darkgreen", Untreated = "grey70")
}

#-------------------------------------------------------------------------------
# 9. DIFFERENTIAL PATHWAY ANALYSIS (TUMOR vs NORMAL)
#-------------------------------------------------------------------------------

cat("=== Differential Pathway Analysis ===\n\n")

run_pathway_diff <- function(gsva_mat, sample_anno, contrast_col = "Condition", 
                              group1 = "TUMOR", group2 = "NORMAL") {
  
  if (is.null(gsva_mat)) return(NULL)
  
  # Match samples
  common_samples <- intersect(colnames(gsva_mat), rownames(sample_anno))
  gsva_mat <- gsva_mat[, common_samples]
  sample_anno <- sample_anno[common_samples, , drop = FALSE]
  
  # Create design matrix
  group <- factor(sample_anno[[contrast_col]])
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  # Fit linear model
  fit <- lmFit(gsva_mat, design)
  
  # Create contrast
  contrast_formula <- paste0(group1, "-", group2)
  contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
  
  # Fit contrasts
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  # Get results
  results <- topTable(fit2, number = Inf, sort.by = "P")
  results$pathway <- rownames(results)
  results <- results[, c("pathway", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
  
  return(results)
}

# Run differential analysis
diff_kegg <- run_pathway_diff(gsva_kegg, sample_annotation)
diff_hallmark <- run_pathway_diff(gsva_hallmark, sample_annotation)
diff_go_bp <- run_pathway_diff(gsva_go_bp, sample_annotation)
diff_reactome <- run_pathway_diff(gsva_reactome, sample_annotation)

# Save differential pathway results
if (!is.null(diff_kegg)) {
  write.csv(diff_kegg, "results/gsva/GSVA_KEGG_differential.csv", row.names = FALSE)
  cat("Significant KEGG pathways (adj.P < 0.05):", sum(diff_kegg$adj.P.Val < 0.05), "\n")
}
if (!is.null(diff_hallmark)) {
  write.csv(diff_hallmark, "results/gsva/GSVA_Hallmark_differential.csv", row.names = FALSE)
  cat("Significant Hallmark pathways (adj.P < 0.05):", sum(diff_hallmark$adj.P.Val < 0.05), "\n")
}
if (!is.null(diff_go_bp)) {
  write.csv(diff_go_bp, "results/gsva/GSVA_GO_BP_differential.csv", row.names = FALSE)
  cat("Significant GO BP terms (adj.P < 0.05):", sum(diff_go_bp$adj.P.Val < 0.05), "\n")
}
if (!is.null(diff_reactome)) {
  write.csv(diff_reactome, "results/gsva/GSVA_Reactome_differential.csv", row.names = FALSE)
  cat("Significant Reactome pathways (adj.P < 0.05):", sum(diff_reactome$adj.P.Val < 0.05), "\n")
}

cat("\n")

#-------------------------------------------------------------------------------
# 10. CREATE HEATMAPS
#-------------------------------------------------------------------------------

cat("=== Creating Heatmaps ===\n\n")

create_gsva_heatmap <- function(gsva_mat, sample_anno, ann_colors, 
                                 diff_results = NULL, title = "GSVA Scores",
                                 filename = "gsva_heatmap.png",
                                 top_n = 50, show_significant = TRUE) {
  
  if (is.null(gsva_mat)) {
    cat("No data for", title, "\n")
    return(NULL)
  }
  
  # Match samples
  common_samples <- intersect(colnames(gsva_mat), rownames(sample_anno))
  gsva_mat <- gsva_mat[, common_samples]
  sample_anno_sub <- sample_anno[common_samples, , drop = FALSE]
  
  # Select pathways to show
  if (!is.null(diff_results) && show_significant) {
    # Show top significant pathways
    sig_pathways <- diff_results %>%
      filter(adj.P.Val < 0.05) %>%
      arrange(adj.P.Val) %>%
      head(top_n) %>%
      pull(pathway)
    
    if (length(sig_pathways) == 0) {
      # If no significant, show top by variance
      row_vars <- apply(gsva_mat, 1, var)
      sig_pathways <- names(sort(row_vars, decreasing = TRUE))[1:min(top_n, nrow(gsva_mat))]
    }
  } else {
    # Show top variable pathways
    row_vars <- apply(gsva_mat, 1, var)
    sig_pathways <- names(sort(row_vars, decreasing = TRUE))[1:min(top_n, nrow(gsva_mat))]
  }
  
  # Subset matrix
  plot_mat <- gsva_mat[sig_pathways, , drop = FALSE]
  
  # Clean pathway names for display
  rownames(plot_mat) <- gsub("^KEGG_|^HALLMARK_|^GOBP_|^REACTOME_", "", rownames(plot_mat))
  rownames(plot_mat) <- gsub("_", " ", rownames(plot_mat))
  rownames(plot_mat) <- substr(rownames(plot_mat), 1, 50)  # Truncate long names
  
  # Create heatmap
  png(paste0("results/gsva/plots/", filename), width = 14, height = 12, units = "in", res = 300)
  pheatmap(
    plot_mat,
    annotation_col = sample_anno_sub,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 10,
    main = title,
    border_color = NA,
    treeheight_row = 30,
    treeheight_col = 30
  )
  dev.off()
  
  # Also save PDF
  pdf(paste0("results/gsva/plots/", gsub(".png$", ".pdf", filename)), width = 14, height = 12)
  pheatmap(
    plot_mat,
    annotation_col = sample_anno_sub,
    annotation_colors = ann_colors,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 8,
    fontsize_col = 10,
    main = title,
    border_color = NA
  )
  dev.off()
  
  cat("Created:", filename, "\n")
}

# Create heatmaps for each gene set collection
create_gsva_heatmap(gsva_kegg, sample_annotation, ann_colors, diff_kegg,
                    "KEGG Pathway Activity (TUMOR vs NORMAL)", "GSVA_KEGG_heatmap.png")

create_gsva_heatmap(gsva_hallmark, sample_annotation, ann_colors, diff_hallmark,
                    "Hallmark Pathway Activity (TUMOR vs NORMAL)", "GSVA_Hallmark_heatmap.png")

create_gsva_heatmap(gsva_go_bp, sample_annotation, ann_colors, diff_go_bp,
                    "GO Biological Process Activity (TUMOR vs NORMAL)", "GSVA_GO_BP_heatmap.png", top_n = 40)

create_gsva_heatmap(gsva_reactome, sample_annotation, ann_colors, diff_reactome,
                    "Reactome Pathway Activity (TUMOR vs NORMAL)", "GSVA_Reactome_heatmap.png")

#-------------------------------------------------------------------------------
# 11. CREATE SEPARATED UP/DOWN HEATMAPS
#-------------------------------------------------------------------------------

cat("\n=== Creating Up/Down Separated Heatmaps ===\n\n")

create_updown_heatmap <- function(gsva_mat, sample_anno, ann_colors, diff_results,
                                   title_prefix = "KEGG", filename_prefix = "GSVA_KEGG",
                                   top_n = 30) {
  
  if (is.null(gsva_mat) || is.null(diff_results)) return(NULL)
  
  # Match samples
  common_samples <- intersect(colnames(gsva_mat), rownames(sample_anno))
  gsva_mat <- gsva_mat[, common_samples]
  sample_anno_sub <- sample_anno[common_samples, , drop = FALSE]
  
  # Get upregulated pathways (higher in TUMOR)
  up_pathways <- diff_results %>%
    filter(adj.P.Val < 0.05, logFC > 0) %>%
    arrange(adj.P.Val) %>%
    head(top_n) %>%
    pull(pathway)
  
  # Get downregulated pathways (lower in TUMOR)
  down_pathways <- diff_results %>%
    filter(adj.P.Val < 0.05, logFC < 0) %>%
    arrange(adj.P.Val) %>%
    head(top_n) %>%
    pull(pathway)
  
  # Upregulated heatmap
  if (length(up_pathways) > 0) {
    plot_mat_up <- gsva_mat[up_pathways, , drop = FALSE]
    rownames(plot_mat_up) <- gsub("^KEGG_|^HALLMARK_|^GOBP_|^REACTOME_", "", rownames(plot_mat_up))
    rownames(plot_mat_up) <- gsub("_", " ", rownames(plot_mat_up))
    rownames(plot_mat_up) <- substr(rownames(plot_mat_up), 1, 50)
    
    png(paste0("results/gsva/plots/", filename_prefix, "_upregulated_heatmap.png"), 
        width = 14, height = max(6, length(up_pathways) * 0.3), units = "in", res = 300)
    pheatmap(
      plot_mat_up,
      annotation_col = sample_anno_sub,
      annotation_colors = ann_colors,
      color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 9,
      fontsize_col = 10,
      main = paste(title_prefix, "Pathways - Upregulated in TUMOR"),
      border_color = NA
    )
    dev.off()
    cat("Created:", filename_prefix, "_upregulated_heatmap.png (", length(up_pathways), "pathways)\n")
  }
  
  # Downregulated heatmap
  if (length(down_pathways) > 0) {
    plot_mat_down <- gsva_mat[down_pathways, , drop = FALSE]
    rownames(plot_mat_down) <- gsub("^KEGG_|^HALLMARK_|^GOBP_|^REACTOME_", "", rownames(plot_mat_down))
    rownames(plot_mat_down) <- gsub("_", " ", rownames(plot_mat_down))
    rownames(plot_mat_down) <- substr(rownames(plot_mat_down), 1, 50)
    
    png(paste0("results/gsva/plots/", filename_prefix, "_downregulated_heatmap.png"), 
        width = 14, height = max(6, length(down_pathways) * 0.3), units = "in", res = 300)
    pheatmap(
      plot_mat_down,
      annotation_col = sample_anno_sub,
      annotation_colors = ann_colors,
      color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 9,
      fontsize_col = 10,
      main = paste(title_prefix, "Pathways - Downregulated in TUMOR"),
      border_color = NA
    )
    dev.off()
    cat("Created:", filename_prefix, "_downregulated_heatmap.png (", length(down_pathways), "pathways)\n")
  }
}

# Create up/down heatmaps
create_updown_heatmap(gsva_kegg, sample_annotation, ann_colors, diff_kegg, "KEGG", "GSVA_KEGG")
create_updown_heatmap(gsva_hallmark, sample_annotation, ann_colors, diff_hallmark, "Hallmark", "GSVA_Hallmark")
create_updown_heatmap(gsva_reactome, sample_annotation, ann_colors, diff_reactome, "Reactome", "GSVA_Reactome")

#-------------------------------------------------------------------------------
# 12. SAVE GSVA SCORES
#-------------------------------------------------------------------------------

cat("\n=== Saving GSVA Scores ===\n\n")

save_gsva_scores <- function(gsva_mat, filename) {
  if (is.null(gsva_mat)) return(NULL)
  
  scores_df <- as.data.frame(gsva_mat)
  scores_df$pathway <- rownames(scores_df)
  scores_df <- scores_df[, c("pathway", setdiff(colnames(scores_df), "pathway"))]
  
  write.csv(scores_df, filename, row.names = FALSE)
  cat("Saved:", filename, "\n")
}

save_gsva_scores(gsva_kegg, "results/gsva/GSVA_KEGG_scores.csv")
save_gsva_scores(gsva_hallmark, "results/gsva/GSVA_Hallmark_scores.csv")
save_gsva_scores(gsva_go_bp, "results/gsva/GSVA_GO_BP_scores.csv")
save_gsva_scores(gsva_reactome, "results/gsva/GSVA_Reactome_scores.csv")

#-------------------------------------------------------------------------------
# 13. BARPLOT OF TOP DIFFERENTIAL PATHWAYS
#-------------------------------------------------------------------------------

cat("\n=== Creating Differential Pathway Barplots ===\n\n")

create_diff_barplot <- function(diff_results, title = "KEGG", filename = "diff_barplot.png", top_n = 20) {
  
  if (is.null(diff_results)) return(NULL)
  
  # Get top pathways by significance
  top_up <- diff_results %>%
    filter(adj.P.Val < 0.05, logFC > 0) %>%
    arrange(desc(logFC)) %>%
    head(top_n)
  
  top_down <- diff_results %>%
    filter(adj.P.Val < 0.05, logFC < 0) %>%
    arrange(logFC) %>%
    head(top_n)
  
  plot_data <- bind_rows(
    top_up %>% mutate(direction = "Upregulated in TUMOR"),
    top_down %>% mutate(direction = "Downregulated in TUMOR")
  )
  
  if (nrow(plot_data) == 0) {
    cat("No significant pathways for", title, "\n")
    return(NULL)
  }
  
  # Clean names
  plot_data$pathway_clean <- gsub("^KEGG_|^HALLMARK_|^GOBP_|^REACTOME_", "", plot_data$pathway)
  plot_data$pathway_clean <- gsub("_", " ", plot_data$pathway_clean)
  plot_data$pathway_clean <- substr(plot_data$pathway_clean, 1, 40)
  
  p <- ggplot(plot_data, aes(x = reorder(pathway_clean, logFC), y = logFC, fill = direction)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("Upregulated in TUMOR" = "firebrick", 
                                  "Downregulated in TUMOR" = "steelblue")) +
    labs(
      title = paste("Differential", title, "Pathway Activity"),
      subtitle = "TUMOR vs NORMAL (adj.P < 0.05)",
      x = "",
      y = "GSVA Score Difference (logFC)",
      fill = ""
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      axis.text.y = element_text(size = 8)
    )
  
  ggsave(paste0("results/gsva/plots/", filename), p, width = 12, height = 10, dpi = 300)
  cat("Created:", filename, "\n")
}

create_diff_barplot(diff_kegg, "KEGG", "GSVA_KEGG_differential_barplot.png")
create_diff_barplot(diff_hallmark, "Hallmark", "GSVA_Hallmark_differential_barplot.png")
create_diff_barplot(diff_reactome, "Reactome", "GSVA_Reactome_differential_barplot.png")

#===============================================================================
# 14. SUMMARY
#===============================================================================

cat("\n")
cat("================================================================\n")
cat("GSVA ANALYSIS COMPLETE\n")
cat("================================================================\n")

cat("\n--- Output Files ---\n")
cat("\nGSVA Scores (per sample):\n")
cat("  results/gsva/GSVA_KEGG_scores.csv\n")
cat("  results/gsva/GSVA_Hallmark_scores.csv\n")
cat("  results/gsva/GSVA_GO_BP_scores.csv\n")
cat("  results/gsva/GSVA_Reactome_scores.csv\n")

cat("\nDifferential Pathway Analysis:\n")
cat("  results/gsva/GSVA_KEGG_differential.csv\n")
cat("  results/gsva/GSVA_Hallmark_differential.csv\n")
cat("  results/gsva/GSVA_GO_BP_differential.csv\n")
cat("  results/gsva/GSVA_Reactome_differential.csv\n")

cat("\nHeatmaps:\n")
cat("  results/gsva/plots/GSVA_KEGG_heatmap.png\n")
cat("  results/gsva/plots/GSVA_Hallmark_heatmap.png\n")
cat("  results/gsva/plots/GSVA_GO_BP_heatmap.png\n")
cat("  results/gsva/plots/GSVA_Reactome_heatmap.png\n")
cat("  results/gsva/plots/GSVA_KEGG_upregulated_heatmap.png\n")
cat("  results/gsva/plots/GSVA_KEGG_downregulated_heatmap.png\n")

cat("\nBarplots:\n")
cat("  results/gsva/plots/GSVA_KEGG_differential_barplot.png\n")
cat("  results/gsva/plots/GSVA_Hallmark_differential_barplot.png\n")
cat("  results/gsva/plots/GSVA_Reactome_differential_barplot.png\n")

# Save session info
sink("results/gsva/session_info.txt")
cat("GSVA analysis completed:", as.character(Sys.time()), "\n\n")
cat("Method:", GSVA_METHOD, "\n")
cat("Min genes per pathway:", MIN_GENES_PER_PATHWAY, "\n")
cat("Max genes per pathway:", MAX_GENES_PER_PATHWAY, "\n\n")
sessionInfo()
sink()

cat("\nSession info saved to results/gsva/session_info.txt\n")
cat("\nDone!\n")
