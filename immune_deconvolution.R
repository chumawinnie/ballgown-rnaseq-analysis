#===============================================================================
# IMMUNE CELL DECONVOLUTION FROM RNA-SEQ DATA
# Using immunedeconv package with xCell method
#===============================================================================

#-------------------------------------------------------------------------------
# 1. INSTALL REQUIRED PACKAGES
#-------------------------------------------------------------------------------

# Install remotes if not available
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install immunedeconv from GitHub
if (!requireNamespace("immunedeconv", quietly = TRUE)) {
  remotes::install_github("omnideconv/immunedeconv")
}

# Install other required packages
packages <- c("ggplot2", "viridis", "RColorBrewer", "reshape2", "dplyr", "tidyr")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

#-------------------------------------------------------------------------------
# 2. LOAD LIBRARIES
#-------------------------------------------------------------------------------

library(immunedeconv)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(tidyr)

#-------------------------------------------------------------------------------
# 3. LOAD GENE EXPRESSION DATA
#-------------------------------------------------------------------------------

# Load gene counts from nf-core/rnaseq output
# Replace with your actual file path
input_file <- "salmon.merged.gene_counts.tsv"

cat("Loading gene expression data from:", input_file, "\n")

gene_expression <- read.table(input_file, header = TRUE, sep = "\t")

cat("Loaded", nrow(gene_expression), "genes from", ncol(gene_expression) - 1, "samples\n")

# Check the first few rows
head(gene_expression)

#-------------------------------------------------------------------------------
# 4. PREPARE EXPRESSION MATRIX
#-------------------------------------------------------------------------------

# Convert to numeric matrix (excluding the gene name column)
gene_expression_matrix <- as.matrix(gene_expression[, -1])

# Ensure row names are gene names
rownames(gene_expression_matrix) <- gene_expression$gene_name

# Remove rows with NA gene names
gene_expression_matrix <- gene_expression_matrix[!is.na(rownames(gene_expression_matrix)), ]

# Remove duplicate gene names (keep first occurrence)
gene_expression_matrix <- gene_expression_matrix[!duplicated(rownames(gene_expression_matrix)), ]

cat("Expression matrix dimensions:", nrow(gene_expression_matrix), "genes x", 
    ncol(gene_expression_matrix), "samples\n")

#-------------------------------------------------------------------------------
# 5. PERFORM IMMUNE CELL DECONVOLUTION
#-------------------------------------------------------------------------------

cat("\nPerforming immune cell deconvolution using xCell...\n")

# Available methods: xcell, mcp_counter, epic, quantiseq, timer, cibersort, cibersort_abs
# Note: CIBERSORT requires registration at cibersort.stanford.edu

result_deconv <- deconvolute(gene_expression_matrix, method = "xcell")

cat("Deconvolution complete!\n")
cat("Cell types identified:", nrow(result_deconv), "\n")

# View results
print(head(result_deconv))

#-------------------------------------------------------------------------------
# 6. SAVE RESULTS
#-------------------------------------------------------------------------------

# Create output directory
dir.create("results/deconvolution", showWarnings = FALSE, recursive = TRUE)

# Write results to CSV
write.csv(result_deconv, "results/deconvolution/deconvolution_results_xcell.csv", row.names = FALSE)

# Transpose for sample-centric view
result_transposed <- result_deconv %>%
  pivot_longer(cols = -cell_type, names_to = "sample", values_to = "score") %>%
  pivot_wider(names_from = cell_type, values_from = score)

write.csv(result_transposed, "results/deconvolution/deconvolution_by_sample.csv", row.names = FALSE)

cat("Results saved to results/deconvolution/\n")

#-------------------------------------------------------------------------------
# 7. VISUALIZATION
#-------------------------------------------------------------------------------

# Convert to long format for ggplot2
result_melted <- melt(result_deconv, id.vars = "cell_type", 
                      variable.name = "sample", value.name = "score")

# Define custom color palette
num_cell_types <- length(unique(result_melted$cell_type))

# Use a good categorical palette
if (num_cell_types <= 12) {
  colors <- brewer.pal(num_cell_types, "Set3")
} else if (num_cell_types <= 20) {
  colors <- c(brewer.pal(12, "Set3"), brewer.pal(8, "Set2"))
} else {
  colors <- colorRampPalette(brewer.pal(12, "Set3"))(num_cell_types)
}

#-------------------------------------------------------------------------------
# 7.1 Stacked Bar Plot
#-------------------------------------------------------------------------------

p_stacked <- ggplot(result_melted, aes(x = sample, y = score, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7) +
  theme_minimal() +
  scale_fill_manual(values = colors) +
  labs(
    title = "Immune Cell Deconvolution Results (xCell)",
    x = "Sample",
    y = "Score",
    fill = "Cell Type"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("results/deconvolution/deconvolution_stacked_barplot.png", p_stacked, 
       width = 14, height = 10, dpi = 300)

#-------------------------------------------------------------------------------
# 7.2 Heatmap
#-------------------------------------------------------------------------------

# Prepare matrix for heatmap
heatmap_matrix <- as.matrix(result_deconv[, -1])
rownames(heatmap_matrix) <- result_deconv$cell_type

# Filter to show only cell types with some signal
row_sums <- rowSums(heatmap_matrix)
heatmap_matrix_filtered <- heatmap_matrix[row_sums > 0, ]

# Scale by row for better visualization
heatmap_scaled <- t(scale(t(heatmap_matrix_filtered)))
heatmap_scaled[is.na(heatmap_scaled)] <- 0

# Create heatmap using pheatmap
library(pheatmap)

png("results/deconvolution/deconvolution_heatmap.png", width = 12, height = 14, 
    units = "in", res = 300)
pheatmap(
  heatmap_scaled,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  main = "Immune Cell Deconvolution Heatmap (xCell)",
  border_color = NA
)
dev.off()

#-------------------------------------------------------------------------------
# 7.3 Dot Plot (selected cell types)
#-------------------------------------------------------------------------------

# Select major immune cell types
major_cell_types <- c(
  "T cells", "CD8+ T cells", "CD4+ T cells", "Tregs",
  "B cells", "NK cells", "Macrophages", "Monocytes",
  "Dendritic cells", "Neutrophils", "Mast cells"
)

result_major <- result_melted %>%
  filter(cell_type %in% major_cell_types)

if (nrow(result_major) > 0) {
  p_dot <- ggplot(result_major, aes(x = sample, y = cell_type)) +
    geom_point(aes(size = score, color = score)) +
    scale_color_viridis(option = "plasma") +
    scale_size_continuous(range = c(1, 8)) +
    theme_minimal() +
    labs(
      title = "Major Immune Cell Types",
      x = "Sample",
      y = "Cell Type",
      size = "Score",
      color = "Score"
    ) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave("results/deconvolution/deconvolution_dotplot.png", p_dot, 
         width = 12, height = 8, dpi = 300)
}

#-------------------------------------------------------------------------------
# 7.4 Box Plot by Cell Type
#-------------------------------------------------------------------------------

# Show distribution of each cell type across samples
p_box <- ggplot(result_melted, aes(x = reorder(cell_type, -score, FUN = median), 
                                   y = score, fill = cell_type)) +
  geom_boxplot(show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  labs(
    title = "Distribution of Immune Cell Scores",
    x = "Cell Type",
    y = "Score"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 8)
  )

ggsave("results/deconvolution/deconvolution_boxplot.png", p_box, 
       width = 10, height = 12, dpi = 300)

#===============================================================================
# 8. ALTERNATIVE DECONVOLUTION METHODS
#===============================================================================

# Uncomment to try other methods:

# # MCP-counter (no external dependencies)
# result_mcp <- deconvolute(gene_expression_matrix, method = "mcp_counter")
# write.csv(result_mcp, "results/deconvolution/deconvolution_mcp_counter.csv", row.names = FALSE)

# # EPIC
# result_epic <- deconvolute(gene_expression_matrix, method = "epic")
# write.csv(result_epic, "results/deconvolution/deconvolution_epic.csv", row.names = FALSE)

# # quanTIseq
# result_quantiseq <- deconvolute(gene_expression_matrix, method = "quantiseq")
# write.csv(result_quantiseq, "results/deconvolution/deconvolution_quantiseq.csv", row.names = FALSE)

#===============================================================================
# 9. SESSION INFO
#===============================================================================

sink("results/deconvolution/session_info.txt")
cat("Immune cell deconvolution completed:", as.character(Sys.time()), "\n\n")
cat("Method: xCell\n")
cat("Input file:", input_file, "\n")
cat("Genes:", nrow(gene_expression_matrix), "\n")
cat("Samples:", ncol(gene_expression_matrix), "\n\n")
sessionInfo()
sink()

cat("\n")
cat("================================================================\n")
cat("DECONVOLUTION COMPLETE\n")
cat("================================================================\n")
cat("\nOutput files:\n")
cat("  results/deconvolution/deconvolution_results_xcell.csv\n")
cat("  results/deconvolution/deconvolution_by_sample.csv\n")
cat("  results/deconvolution/deconvolution_stacked_barplot.png\n")
cat("  results/deconvolution/deconvolution_heatmap.png\n")
cat("  results/deconvolution/deconvolution_dotplot.png\n")
cat("  results/deconvolution/deconvolution_boxplot.png\n")
cat("\nDone!\n")
