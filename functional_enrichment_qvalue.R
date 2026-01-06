#===============================================================================
# FUNCTIONAL ENRICHMENT ANALYSIS
# Gene Ontology (topGO) and KEGG Pathway Analysis
# CORRECTED VERSION - Using Q-VALUE (FDR) for significance
#===============================================================================

#-------------------------------------------------------------------------------
# 1. SETUP AND DEPENDENCIES
#-------------------------------------------------------------------------------

# Install packages if needed (uncomment if required)
# BiocManager::install(c("topGO", "clusterProfiler", "org.Hs.eg.db", "DOSE", "enrichplot", "pathview"))
# install.packages(c("ggplot2", "dplyr", "tidyr", "RColorBrewer"))

library(topGO)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(pathview)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Resolve conflicts
select <- dplyr::select
filter <- dplyr::filter

# Create output directories
dir.create("results/enrichment", showWarnings = FALSE)
dir.create("results/enrichment/GO", showWarnings = FALSE)
dir.create("results/enrichment/KEGG", showWarnings = FALSE)
dir.create("results/enrichment/plots", showWarnings = FALSE)

#-------------------------------------------------------------------------------
# 2. LOAD DIFFERENTIAL EXPRESSION RESULTS
#-------------------------------------------------------------------------------

# Load your DE results (from Ballgown analysis)
results_genes <- read.csv("results/all_genes_DE_results.csv")

cat("Loaded", nrow(results_genes), "genes\n")
cat("Columns:", colnames(results_genes), "\n")

#-------------------------------------------------------------------------------
# IMPORTANT: Define significance thresholds using Q-VALUE (FDR)
#-------------------------------------------------------------------------------
qval_cutoff <- 0.05      # FDR-adjusted p-value threshold
log2fc_cutoff <- 1       # |log2FC| >= 1 means 2-fold change

# Get significant genes using Q-VALUE
sig_genes <- results_genes %>%
 filter(qval < qval_cutoff)

sig_genes_stringent <- results_genes %>%
 filter(qval < qval_cutoff, abs(log2FC) >= log2fc_cutoff)

# Separate up and down regulated (using q-value)
sig_up <- results_genes %>%
 filter(qval < qval_cutoff, log2FC > 0)

sig_down <- results_genes %>%
 filter(qval < qval_cutoff, log2FC < 0)

# Stringent up/down
sig_up_stringent <- results_genes %>%
 filter(qval < qval_cutoff, log2FC >= log2fc_cutoff)

sig_down_stringent <- results_genes %>%
 filter(qval < qval_cutoff, log2FC <= -log2fc_cutoff)

cat("\n--- Significant Genes Summary (using q-value < ", qval_cutoff, ") ---\n")
cat("Total significant (q <", qval_cutoff, "):", nrow(sig_genes), "\n")
cat("Stringent (q <", qval_cutoff, ", |log2FC| >=", log2fc_cutoff, "):", nrow(sig_genes_stringent), "\n")
cat("Upregulated (q <", qval_cutoff, "):", nrow(sig_up), "\n")
cat("Downregulated (q <", qval_cutoff, "):", nrow(sig_down), "\n")
cat("Upregulated stringent (log2FC >=", log2fc_cutoff, "):", nrow(sig_up_stringent), "\n")
cat("Downregulated stringent (log2FC <= -", log2fc_cutoff, "):", nrow(sig_down_stringent), "\n")

#===============================================================================
# 3. GENE ONTOLOGY ANALYSIS WITH topGO
#===============================================================================

cat("\n")
cat("================================================================\n")
cat("GENE ONTOLOGY ENRICHMENT (topGO)\n")
cat("================================================================\n")

#-------------------------------------------------------------------------------
# 3.1 Prepare gene list for topGO
#-------------------------------------------------------------------------------

# Clean Ensembl IDs (remove version if present)
results_genes$ensembl_clean <- gsub("\\..*", "", results_genes$ensembl_id)

# Create named vector: gene scores (using Q-VALUES)
all_genes <- results_genes$qval
names(all_genes) <- results_genes$ensembl_clean

# Remove NAs
all_genes <- all_genes[!is.na(all_genes)]

cat("Total genes for GO analysis:", length(all_genes), "\n")

# Selection function (genes with q < 0.05)
topDiffGenes <- function(allScore) {
 return(allScore < qval_cutoff)
}

cat("Significant genes for topGO:", sum(topDiffGenes(all_genes)), "\n")

#-------------------------------------------------------------------------------
# 3.2 Run topGO for each ontology (BP, MF, CC)
#-------------------------------------------------------------------------------

run_topGO <- function(ontology, all_genes, gene_selection_fun) {
 
 cat("\nRunning topGO for", ontology, "...\n")
 
 # Create topGO object
 GOdata <- new("topGOdata",
               ontology = ontology,
               allGenes = all_genes,
               geneSel = gene_selection_fun,
               annot = annFUN.org,
               mapping = "org.Hs.eg.db",
               ID = "ENSEMBL")
 
 # Print summary
 cat("Genes in GO annotation:", length(genes(GOdata)), "\n")
 cat("Significant genes:", sum(gene_selection_fun(all_genes)), "\n")
 cat("GO terms available:", length(usedGO(GOdata)), "\n")
 
 # Run enrichment tests
 # Fisher's exact test
 result_fisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
 
 # Weight01 algorithm (recommended - accounts for GO hierarchy)
 result_weight <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
 
 # Kolmogorov-Smirnov test (uses all gene scores, not just sig/not sig)
 result_ks <- runTest(GOdata, algorithm = "classic", statistic = "ks")
 
 # Generate results table
 allRes <- GenTable(GOdata,
                    classicFisher = result_fisher,
                    weight01Fisher = result_weight,
                    classicKS = result_ks,
                    orderBy = "weight01Fisher",
                    ranksOf = "classicFisher",
                    topNodes = 50)
 
 # Add ontology column
 allRes$Ontology <- ontology
 
 # Convert p-values to numeric (they come as characters with "<" symbols)
 allRes$classicFisher <- as.numeric(gsub("< ", "", allRes$classicFisher))
 allRes$weight01Fisher <- as.numeric(gsub("< ", "", allRes$weight01Fisher))
 allRes$classicKS <- as.numeric(gsub("< ", "", allRes$classicKS))
 
 return(list(
   GOdata = GOdata,
   result_fisher = result_fisher,
   result_weight = result_weight,
   result_ks = result_ks,
   table = allRes
 ))
}

# Run for all three ontologies
GO_BP <- run_topGO("BP", all_genes, topDiffGenes)  # Biological Process
GO_MF <- run_topGO("MF", all_genes, topDiffGenes)  # Molecular Function
GO_CC <- run_topGO("CC", all_genes, topDiffGenes)  # Cellular Component

#-------------------------------------------------------------------------------
# 3.3 Save GO results
#-------------------------------------------------------------------------------

# Combine all results
GO_all_results <- bind_rows(
 GO_BP$table,
 GO_MF$table,
 GO_CC$table
)

# Save tables
write.csv(GO_BP$table, "results/enrichment/GO/GO_BP_results.csv", row.names = FALSE)
write.csv(GO_MF$table, "results/enrichment/GO/GO_MF_results.csv", row.names = FALSE)
write.csv(GO_CC$table, "results/enrichment/GO/GO_CC_results.csv", row.names = FALSE)
write.csv(GO_all_results, "results/enrichment/GO/GO_all_results.csv", row.names = FALSE)

cat("\nGO results saved to results/enrichment/GO/\n")

#-------------------------------------------------------------------------------
# 3.4 Visualize GO results
#-------------------------------------------------------------------------------

# Function to create GO bar plot
plot_GO_bar <- function(go_table, ontology_name, top_n = 15) {
 
 # Get top terms
 plot_data <- go_table %>%
   arrange(weight01Fisher) %>%
   head(top_n) %>%
   mutate(
     Term = factor(Term, levels = rev(Term)),
     neg_log10_pval = -log10(weight01Fisher),
     GeneRatio = Significant / Annotated
   )
 
 p <- ggplot(plot_data, aes(x = neg_log10_pval, y = Term)) +
   geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
   geom_text(aes(label = paste0(Significant, "/", Annotated)), 
             hjust = -0.1, size = 3) +
   labs(
     title = paste("GO Enrichment:", ontology_name),
     subtitle = paste("Based on q-value <", qval_cutoff),
     x = "-log10(p-value)",
     y = ""
   ) +
   theme_bw() +
   theme(
     plot.title = element_text(hjust = 0.5, face = "bold"),
     axis.text.y = element_text(size = 9)
   ) +
   xlim(0, max(plot_data$neg_log10_pval) * 1.2)
 
 return(p)
}

# Create bar plots
p_BP <- plot_GO_bar(GO_BP$table, "Biological Process")
p_MF <- plot_GO_bar(GO_MF$table, "Molecular Function")
p_CC <- plot_GO_bar(GO_CC$table, "Cellular Component")

ggsave("results/enrichment/plots/GO_BP_barplot.png", p_BP, width = 12, height = 8, dpi = 300)
ggsave("results/enrichment/plots/GO_MF_barplot.png", p_MF, width = 12, height = 8, dpi = 300)
ggsave("results/enrichment/plots/GO_CC_barplot.png", p_CC, width = 12, height = 8, dpi = 300)

# Combined dot plot
plot_GO_dotplot <- function(go_results, top_n = 10) {
 
 plot_data <- go_results %>%
   group_by(Ontology) %>%
   arrange(weight01Fisher) %>%
   slice_head(n = top_n) %>%
   ungroup() %>%
   mutate(
     Term = factor(Term, levels = rev(unique(Term))),
     neg_log10_pval = -log10(weight01Fisher),
     GeneRatio = Significant / Annotated
   )
 
 p <- ggplot(plot_data, aes(x = GeneRatio, y = Term)) +
   geom_point(aes(size = Significant, color = neg_log10_pval)) +
   scale_color_gradient(low = "blue", high = "red", name = "-log10(p-value)") +
   scale_size_continuous(name = "Gene Count", range = c(3, 10)) +
   facet_grid(Ontology ~ ., scales = "free_y", space = "free_y") +
   labs(
     title = "GO Enrichment Analysis",
     subtitle = paste("Significant genes: q-value <", qval_cutoff),
     x = "Gene Ratio (Significant/Annotated)",
     y = ""
   ) +
   theme_bw() +
   theme(
     plot.title = element_text(hjust = 0.5, face = "bold"),
     strip.text = element_text(face = "bold"),
     axis.text.y = element_text(size = 8)
   )
 
 return(p)
}

p_GO_combined <- plot_GO_dotplot(GO_all_results, top_n = 8)
ggsave("results/enrichment/plots/GO_combined_dotplot.png", p_GO_combined, width = 12, height = 14, dpi = 300)
ggsave("results/enrichment/plots/GO_combined_dotplot.pdf", p_GO_combined, width = 12, height = 14)

#-------------------------------------------------------------------------------
# 3.5 topGO graph visualization (optional - for top terms)
#-------------------------------------------------------------------------------

# Save GO graph for top BP terms
pdf("results/enrichment/plots/GO_BP_graph.pdf", width = 12, height = 10)
showSigOfNodes(GO_BP$GOdata, score(GO_BP$result_weight), firstSigNodes = 10, useInfo = "all")
dev.off()

cat("GO plots saved to results/enrichment/plots/\n")

#===============================================================================
# 4. KEGG PATHWAY ANALYSIS WITH clusterProfiler
#===============================================================================

cat("\n")
cat("================================================================\n")
cat("KEGG PATHWAY ENRICHMENT (clusterProfiler)\n")
cat("================================================================\n")

#-------------------------------------------------------------------------------
# 4.1 Convert gene IDs to Entrez IDs (required for KEGG)
#-------------------------------------------------------------------------------

# Map Ensembl IDs to Entrez IDs
ensembl_to_entrez <- AnnotationDbi::select(
 org.Hs.eg.db,
 keys = unique(results_genes$ensembl_clean),
 keytype = "ENSEMBL",
 columns = c("ENSEMBL", "ENTREZID", "SYMBOL")
)

# Remove NAs and duplicates
ensembl_to_entrez <- ensembl_to_entrez[!is.na(ensembl_to_entrez$ENTREZID), ]
ensembl_to_entrez <- ensembl_to_entrez[!duplicated(ensembl_to_entrez$ENSEMBL), ]

cat("Mapped", nrow(ensembl_to_entrez), "genes to Entrez IDs\n")

# Merge with results
results_with_entrez <- results_genes %>%
 left_join(ensembl_to_entrez, by = c("ensembl_clean" = "ENSEMBL"))

#-------------------------------------------------------------------------------
# 4.2 Prepare gene lists for KEGG (using Q-VALUE)
#-------------------------------------------------------------------------------

# All significant genes (Entrez IDs) - using q-value
sig_entrez <- results_with_entrez %>%
 filter(qval < qval_cutoff, !is.na(ENTREZID)) %>%
 pull(ENTREZID) %>%
 unique()

# Upregulated genes (q < 0.05)
up_entrez <- results_with_entrez %>%
 filter(qval < qval_cutoff, log2FC > 0, !is.na(ENTREZID)) %>%
 pull(ENTREZID) %>%
 unique()

# Downregulated genes (q < 0.05)
down_entrez <- results_with_entrez %>%
 filter(qval < qval_cutoff, log2FC < 0, !is.na(ENTREZID)) %>%
 pull(ENTREZID) %>%
 unique()

# Stringent upregulated (q < 0.05 & log2FC >= 1)
up_stringent_entrez <- results_with_entrez %>%
 filter(qval < qval_cutoff, log2FC >= log2fc_cutoff, !is.na(ENTREZID)) %>%
 pull(ENTREZID) %>%
 unique()

# Stringent downregulated (q < 0.05 & log2FC <= -1)
down_stringent_entrez <- results_with_entrez %>%
 filter(qval < qval_cutoff, log2FC <= -log2fc_cutoff, !is.na(ENTREZID)) %>%
 pull(ENTREZID) %>%
 unique()

# Background (all tested genes)
background_entrez <- results_with_entrez %>%
 filter(!is.na(ENTREZID)) %>%
 pull(ENTREZID) %>%
 unique()

cat("\nGenes for KEGG analysis (using q-value <", qval_cutoff, "):\n")
cat("- All significant:", length(sig_entrez), "\n")
cat("- Upregulated:", length(up_entrez), "\n")
cat("- Downregulated:", length(down_entrez), "\n")
cat("- Upregulated stringent (log2FC >=", log2fc_cutoff, "):", length(up_stringent_entrez), "\n")
cat("- Downregulated stringent (log2FC <= -", log2fc_cutoff, "):", length(down_stringent_entrez), "\n")
cat("- Background:", length(background_entrez), "\n")

#-------------------------------------------------------------------------------
# 4.3 Run KEGG Over-Representation Analysis (ORA)
#-------------------------------------------------------------------------------

cat("\nRunning KEGG ORA...\n")

# All significant genes
if (length(sig_entrez) >= 5) {
 kegg_ora_all <- enrichKEGG(
   gene = sig_entrez,
   universe = background_entrez,
   organism = "hsa",  # Human
   keyType = "ncbi-geneid",
   pAdjustMethod = "BH",
   pvalueCutoff = 0.05,
   qvalueCutoff = 0.2
 )
 cat("KEGG ORA (all):", nrow(as.data.frame(kegg_ora_all)), "pathways\n")
} else {
 kegg_ora_all <- NULL
 cat("KEGG ORA (all): Not enough genes\n")
}

# Downregulated genes only (main analysis since no significant upregulated)
if (length(down_entrez) >= 5) {
 kegg_ora_down <- enrichKEGG(
   gene = down_entrez,
   universe = background_entrez,
   organism = "hsa",
   keyType = "ncbi-geneid",
   pAdjustMethod = "BH",
   pvalueCutoff = 0.05,
   qvalueCutoff = 0.2
 )
 cat("KEGG ORA (downregulated):", nrow(as.data.frame(kegg_ora_down)), "pathways\n")
} else {
 kegg_ora_down <- NULL
 cat("KEGG ORA (downregulated): Not enough genes\n")
}

# Stringent downregulated
if (length(down_stringent_entrez) >= 5) {
 kegg_ora_down_stringent <- enrichKEGG(
   gene = down_stringent_entrez,
   universe = background_entrez,
   organism = "hsa",
   keyType = "ncbi-geneid",
   pAdjustMethod = "BH",
   pvalueCutoff = 0.05,
   qvalueCutoff = 0.2
 )
 cat("KEGG ORA (downregulated stringent):", nrow(as.data.frame(kegg_ora_down_stringent)), "pathways\n")
} else {
 kegg_ora_down_stringent <- NULL
 cat("KEGG ORA (downregulated stringent): Not enough genes (", length(down_stringent_entrez), ")\n")
}

# Stringent upregulated (q < 0.05 & log2FC >= 1)
if (length(up_stringent_entrez) >= 5) {
 kegg_ora_up_stringent <- enrichKEGG(
   gene = up_stringent_entrez,
   universe = background_entrez,
   organism = "hsa",
   keyType = "ncbi-geneid",
   pAdjustMethod = "BH",
   pvalueCutoff = 0.05,
   qvalueCutoff = 0.2
 )
 cat("KEGG ORA (upregulated stringent):", nrow(as.data.frame(kegg_ora_up_stringent)), "pathways\n")
} else {
 kegg_ora_up_stringent <- NULL
 cat("KEGG ORA (upregulated stringent): Not enough genes (", length(up_stringent_entrez), ")\n")
}

# Upregulated genes (likely empty or few)
if (length(up_entrez) >= 5) {
 kegg_ora_up <- enrichKEGG(
   gene = up_entrez,
   universe = background_entrez,
   organism = "hsa",
   keyType = "ncbi-geneid",
   pAdjustMethod = "BH",
   pvalueCutoff = 0.05,
   qvalueCutoff = 0.2
 )
 cat("KEGG ORA (upregulated):", nrow(as.data.frame(kegg_ora_up)), "pathways\n")
} else {
 kegg_ora_up <- NULL
 cat("KEGG ORA (upregulated): Not enough genes (", length(up_entrez), ")\n")
}

#-------------------------------------------------------------------------------
# 4.4 Run KEGG Gene Set Enrichment Analysis (GSEA)
#-------------------------------------------------------------------------------

cat("\nRunning KEGG GSEA...\n")

# GSEA requires a ranked gene list
# Create ranked list: sign(log2FC) * -log10(qval)
gene_list_gsea <- results_with_entrez %>%
 filter(!is.na(ENTREZID), !is.na(log2FC), !is.na(qval)) %>%
 mutate(rank_score = sign(log2FC) * -log10(qval + 1e-300)) %>%
 arrange(desc(rank_score)) %>%
 distinct(ENTREZID, .keep_all = TRUE)

# Create named vector
gsea_ranks <- gene_list_gsea$rank_score
names(gsea_ranks) <- gene_list_gsea$ENTREZID

# Sort descending
gsea_ranks <- sort(gsea_ranks, decreasing = TRUE)

cat("Genes in GSEA ranked list:", length(gsea_ranks), "\n")

# Run GSEA
kegg_gsea <- gseKEGG(
 geneList = gsea_ranks,
 organism = "hsa",
 keyType = "ncbi-geneid",
 minGSSize = 10,
 maxGSSize = 500,
 pvalueCutoff = 0.05,
 pAdjustMethod = "BH",
 verbose = FALSE
)

cat("KEGG GSEA:", nrow(as.data.frame(kegg_gsea)), "significant pathways\n")

#-------------------------------------------------------------------------------
# 4.5 Save KEGG results
#-------------------------------------------------------------------------------

# Save as CSV
if (!is.null(kegg_ora_all) && nrow(as.data.frame(kegg_ora_all)) > 0) {
 write.csv(as.data.frame(kegg_ora_all), "results/enrichment/KEGG/KEGG_ORA_all.csv", row.names = FALSE)
}
if (!is.null(kegg_ora_down) && nrow(as.data.frame(kegg_ora_down)) > 0) {
 write.csv(as.data.frame(kegg_ora_down), "results/enrichment/KEGG/KEGG_ORA_downregulated.csv", row.names = FALSE)
}
if (!is.null(kegg_ora_down_stringent) && nrow(as.data.frame(kegg_ora_down_stringent)) > 0) {
 write.csv(as.data.frame(kegg_ora_down_stringent), "results/enrichment/KEGG/KEGG_ORA_downregulated_stringent.csv", row.names = FALSE)
}
if (!is.null(kegg_ora_up_stringent) && nrow(as.data.frame(kegg_ora_up_stringent)) > 0) {
 write.csv(as.data.frame(kegg_ora_up_stringent), "results/enrichment/KEGG/KEGG_ORA_upregulated_stringent.csv", row.names = FALSE)
}
if (!is.null(kegg_ora_up) && nrow(as.data.frame(kegg_ora_up)) > 0) {
 write.csv(as.data.frame(kegg_ora_up), "results/enrichment/KEGG/KEGG_ORA_upregulated.csv", row.names = FALSE)
}
if (nrow(as.data.frame(kegg_gsea)) > 0) {
 write.csv(as.data.frame(kegg_gsea), "results/enrichment/KEGG/KEGG_GSEA.csv", row.names = FALSE)
}

cat("\nKEGG results saved to results/enrichment/KEGG/\n")

#-------------------------------------------------------------------------------
# 4.6 Visualize KEGG results
#-------------------------------------------------------------------------------

# Dot plot - ORA all
if (!is.null(kegg_ora_all) && nrow(as.data.frame(kegg_ora_all)) > 0) {
 p_kegg_dot <- dotplot(kegg_ora_all, showCategory = 20, 
                       title = paste("KEGG Pathway Enrichment (q <", qval_cutoff, ")"))
 ggsave("results/enrichment/plots/KEGG_ORA_dotplot.png", p_kegg_dot, width = 10, height = 10, dpi = 300)
 ggsave("results/enrichment/plots/KEGG_ORA_dotplot.pdf", p_kegg_dot, width = 10, height = 10)
}

# Bar plot
if (!is.null(kegg_ora_all) && nrow(as.data.frame(kegg_ora_all)) > 0) {
 p_kegg_bar <- barplot(kegg_ora_all, showCategory = 15, 
                       title = paste("KEGG Pathway Enrichment (q <", qval_cutoff, ")"))
 ggsave("results/enrichment/plots/KEGG_ORA_barplot.png", p_kegg_bar, width = 10, height = 8, dpi = 300)
}

# Downregulated pathways
if (!is.null(kegg_ora_down) && nrow(as.data.frame(kegg_ora_down)) > 0) {
 p_kegg_down <- dotplot(kegg_ora_down, showCategory = 20, 
                        title = "KEGG Pathways - Downregulated Genes (q < 0.05)")
 ggsave("results/enrichment/plots/KEGG_ORA_downregulated_dotplot.png", p_kegg_down, width = 10, height = 10, dpi = 300)
}

# Upregulated pathways
if (!is.null(kegg_ora_up) && nrow(as.data.frame(kegg_ora_up)) > 0) {
 p_kegg_up <- dotplot(kegg_ora_up, showCategory = 20, 
                      title = "KEGG Pathways - Upregulated Genes (q < 0.05)")
 ggsave("results/enrichment/plots/KEGG_ORA_upregulated_dotplot.png", p_kegg_up, width = 10, height = 10, dpi = 300)
}

# Stringent downregulated pathways
if (!is.null(kegg_ora_down_stringent) && nrow(as.data.frame(kegg_ora_down_stringent)) > 0) {
 p_kegg_down_str <- dotplot(kegg_ora_down_stringent, showCategory = 20, 
                            title = paste0("KEGG Pathways - Downregulated (q < ", qval_cutoff, ", log2FC <= -", log2fc_cutoff, ")"))
 ggsave("results/enrichment/plots/KEGG_ORA_downregulated_stringent_dotplot.png", p_kegg_down_str, width = 10, height = 10, dpi = 300)
}

# Stringent upregulated pathways
if (!is.null(kegg_ora_up_stringent) && nrow(as.data.frame(kegg_ora_up_stringent)) > 0) {
 p_kegg_up_str <- dotplot(kegg_ora_up_stringent, showCategory = 20, 
                          title = paste0("KEGG Pathways - Upregulated (q < ", qval_cutoff, ", log2FC >= ", log2fc_cutoff, ")"))
 ggsave("results/enrichment/plots/KEGG_ORA_upregulated_stringent_dotplot.png", p_kegg_up_str, width = 10, height = 10, dpi = 300)
}

# GSEA ridge plot
if (nrow(as.data.frame(kegg_gsea)) > 0) {
 p_gsea_ridge <- ridgeplot(kegg_gsea, showCategory = 15) +
   labs(title = "KEGG GSEA - Ridge Plot")
 ggsave("results/enrichment/plots/KEGG_GSEA_ridgeplot.png", p_gsea_ridge, width = 12, height = 10, dpi = 300)
 
 # GSEA dot plot
 p_gsea_dot <- dotplot(kegg_gsea, showCategory = 20, title = "KEGG GSEA")
 ggsave("results/enrichment/plots/KEGG_GSEA_dotplot.png", p_gsea_dot, width = 10, height = 10, dpi = 300)
 
 # GSEA running score plot for top pathway
 if (nrow(as.data.frame(kegg_gsea)) >= 1) {
   p_gsea_plot <- gseaplot2(kegg_gsea, geneSetID = 1:min(3, nrow(as.data.frame(kegg_gsea))), 
                            title = "KEGG GSEA - Top Pathways")
   ggsave("results/enrichment/plots/KEGG_GSEA_running_score.png", p_gsea_plot, width = 10, height = 8, dpi = 300)
 }
}

#-------------------------------------------------------------------------------
# 4.7 Pathway visualization with pathview
#-------------------------------------------------------------------------------

# Create gene data for pathview (log2FC values with Entrez IDs)
gene_data <- gene_list_gsea$log2FC
names(gene_data) <- gene_list_gsea$ENTREZID

# Visualize top pathways
if (!is.null(kegg_ora_all) && nrow(as.data.frame(kegg_ora_all)) > 0) {
 
 cat("\nGenerating pathway maps with pathview...\n")
 
 # Get top 5 pathways
 top_pathways <- head(as.data.frame(kegg_ora_all)$ID, 5)
 
 # Set output directory
 setwd("results/enrichment/KEGG")
 
 for (pathway in top_pathways) {
   tryCatch({
     pathview(
       gene.data = gene_data,
       pathway.id = pathway,
       species = "hsa",
       limit = list(gene = max(abs(gene_data)), cpd = 1),
       low = list(gene = "blue"),
       high = list(gene = "red"),
       mid = list(gene = "white")
     )
     cat("  Generated pathway map for:", pathway, "\n")
   }, error = function(e) {
     cat("  Could not generate pathway map for:", pathway, "\n")
   })
 }
 
 # Return to main directory
 setwd("../../..")
}

#===============================================================================
# 5. GO ENRICHMENT WITH clusterProfiler (Alternative to topGO)
#===============================================================================

cat("\n")
cat("================================================================\n")
cat("GO ENRICHMENT WITH clusterProfiler (using q-value)\n")
cat("================================================================\n")

# Run GO ORA with significant genes (q < 0.05)
if (length(sig_entrez) >= 5) {
 
 go_ora_bp <- enrichGO(
   gene = sig_entrez,
   universe = background_entrez,
   OrgDb = org.Hs.eg.db,
   keyType = "ENTREZID",
   ont = "BP",
   pAdjustMethod = "BH",
   pvalueCutoff = 0.05,
   qvalueCutoff = 0.2,
   readable = TRUE
 )
 
 go_ora_mf <- enrichGO(
   gene = sig_entrez,
   universe = background_entrez,
   OrgDb = org.Hs.eg.db,
   keyType = "ENTREZID",
   ont = "MF",
   pAdjustMethod = "BH",
   pvalueCutoff = 0.05,
   qvalueCutoff = 0.2,
   readable = TRUE
 )
 
 go_ora_cc <- enrichGO(
   gene = sig_entrez,
   universe = background_entrez,
   OrgDb = org.Hs.eg.db,
   keyType = "ENTREZID",
   ont = "CC",
   pAdjustMethod = "BH",
   pvalueCutoff = 0.05,
   qvalueCutoff = 0.2,
   readable = TRUE
 )
 
 cat("clusterProfiler GO ORA Results:\n")
 cat("- BP:", nrow(as.data.frame(go_ora_bp)), "terms\n")
 cat("- MF:", nrow(as.data.frame(go_ora_mf)), "terms\n")
 cat("- CC:", nrow(as.data.frame(go_ora_cc)), "terms\n")
 
 # Save results
 if (nrow(as.data.frame(go_ora_bp)) > 0) {
   write.csv(as.data.frame(go_ora_bp), "results/enrichment/GO/GO_BP_clusterProfiler.csv", row.names = FALSE)
 }
 if (nrow(as.data.frame(go_ora_mf)) > 0) {
   write.csv(as.data.frame(go_ora_mf), "results/enrichment/GO/GO_MF_clusterProfiler.csv", row.names = FALSE)
 }
 if (nrow(as.data.frame(go_ora_cc)) > 0) {
   write.csv(as.data.frame(go_ora_cc), "results/enrichment/GO/GO_CC_clusterProfiler.csv", row.names = FALSE)
 }
 
 #-------------------------------------------------------------------------------
 # 5.1 Advanced GO visualizations
 #-------------------------------------------------------------------------------
 
 # Enrichment map (shows relationships between terms)
 if (nrow(as.data.frame(go_ora_bp)) >= 2) {
   go_ora_bp_sim <- pairwise_termsim(go_ora_bp)
   
   p_emap <- emapplot(go_ora_bp_sim, showCategory = 30) +
     labs(title = "GO BP - Enrichment Map (q < 0.05)")
   ggsave("results/enrichment/plots/GO_BP_enrichment_map.png", p_emap, width = 14, height = 12, dpi = 300)
 }
 
 # Cnet plot (shows gene-concept relationships)
 if (nrow(as.data.frame(go_ora_bp)) > 0) {
   p_cnet <- cnetplot(go_ora_bp, 
                      showCategory = 10,
                      categorySize = "pvalue",
                      foldChange = gene_data) +
     labs(title = "GO BP - Gene-Concept Network")
   ggsave("results/enrichment/plots/GO_BP_cnetplot.png", p_cnet, width = 14, height = 12, dpi = 300)
 }
 
 # Upset plot (shows overlapping genes between terms)
 if (nrow(as.data.frame(go_ora_bp)) >= 5) {
   p_upset <- upsetplot(go_ora_bp, n = 10)
   ggsave("results/enrichment/plots/GO_BP_upset.png", p_upset, width = 12, height = 8, dpi = 300)
 }
 
 # Tree plot (hierarchical clustering of terms)
 if (nrow(as.data.frame(go_ora_bp)) >= 2) {
   go_ora_bp_sim <- pairwise_termsim(go_ora_bp)
   p_tree <- treeplot(go_ora_bp_sim, showCategory = 30)
   ggsave("results/enrichment/plots/GO_BP_treeplot.png", p_tree, width = 14, height = 10, dpi = 300)
 }
 
} else {
 cat("Not enough significant genes for clusterProfiler GO analysis\n")
 go_ora_bp <- NULL
 go_ora_mf <- NULL
 go_ora_cc <- NULL
}

#-------------------------------------------------------------------------------
# 5.2 GO GSEA
#-------------------------------------------------------------------------------

cat("\nRunning GO GSEA...\n")

go_gsea_bp <- gseGO(
 geneList = gsea_ranks,
 OrgDb = org.Hs.eg.db,
 keyType = "ENTREZID",
 ont = "BP",
 minGSSize = 10,
 maxGSSize = 500,
 pvalueCutoff = 0.05,
 pAdjustMethod = "BH",
 verbose = FALSE
)

cat("GO GSEA BP:", nrow(as.data.frame(go_gsea_bp)), "significant terms\n")

if (nrow(as.data.frame(go_gsea_bp)) > 0) {
 write.csv(as.data.frame(go_gsea_bp), "results/enrichment/GO/GO_GSEA_BP.csv", row.names = FALSE)
 
 p_go_gsea <- dotplot(go_gsea_bp, showCategory = 20, title = "GO GSEA - Biological Process")
 ggsave("results/enrichment/plots/GO_GSEA_BP_dotplot.png", p_go_gsea, width = 10, height = 10, dpi = 300)
}

#===============================================================================
# 6. SUMMARY AND SESSION INFO
#===============================================================================

cat("\n")
cat("================================================================\n")
cat("ENRICHMENT ANALYSIS COMPLETE\n")
cat("================================================================\n")

cat("\n--- Analysis Parameters ---\n")
cat("Significance threshold: q-value (FDR) <", qval_cutoff, "\n")
cat("Fold change threshold: |log2FC| >=", log2fc_cutoff, "\n")

cat("\n--- Input Summary ---\n")
cat("Total genes:", nrow(results_genes), "\n")
cat("Significant genes (q <", qval_cutoff, "):", nrow(sig_genes), "\n")
cat("  - Upregulated:", nrow(sig_up), "\n")
cat("  - Downregulated:", nrow(sig_down), "\n")

cat("\n--- GO Analysis (topGO) ---\n")
cat("  BP terms:", nrow(GO_BP$table), "\n")
cat("  MF terms:", nrow(GO_MF$table), "\n")
cat("  CC terms:", nrow(GO_CC$table), "\n")

cat("\n--- KEGG Analysis ---\n")
if (!is.null(kegg_ora_all)) cat("  ORA pathways (all):", nrow(as.data.frame(kegg_ora_all)), "\n")
if (!is.null(kegg_ora_down)) cat("  ORA pathways (down):", nrow(as.data.frame(kegg_ora_down)), "\n")
if (!is.null(kegg_ora_up)) cat("  ORA pathways (up):", nrow(as.data.frame(kegg_ora_up)), "\n")
if (!is.null(kegg_ora_down_stringent)) cat("  ORA pathways (down stringent):", nrow(as.data.frame(kegg_ora_down_stringent)), "\n")
if (!is.null(kegg_ora_up_stringent)) cat("  ORA pathways (up stringent):", nrow(as.data.frame(kegg_ora_up_stringent)), "\n")
cat("  GSEA pathways:", nrow(as.data.frame(kegg_gsea)), "\n")

cat("\n--- Output Directories ---\n")
cat("  results/enrichment/GO/\n")
cat("  results/enrichment/KEGG/\n")
cat("  results/enrichment/plots/\n")

cat("\n--- Key Output Files ---\n")
cat("  GO_all_results.csv - Combined GO results from topGO\n")
cat("  GO_BP/MF/CC_clusterProfiler.csv - GO results from clusterProfiler\n")
cat("  KEGG_ORA_all.csv - KEGG over-representation analysis\n")
cat("  KEGG_ORA_downregulated.csv - KEGG for downregulated genes\n")
cat("  KEGG_ORA_upregulated.csv - KEGG for upregulated genes\n")
cat("  KEGG_ORA_downregulated_stringent.csv - KEGG for stringent downregulated genes\n")
cat("  KEGG_ORA_upregulated_stringent.csv - KEGG for stringent upregulated genes\n")
cat("  KEGG_GSEA.csv - KEGG gene set enrichment analysis\n")

# Save session info
sink("results/enrichment/session_info.txt")
cat("Enrichment analysis completed:", as.character(Sys.time()), "\n\n")
cat("SIGNIFICANCE THRESHOLDS:\n")
cat("  q-value (FDR):", qval_cutoff, "\n")
cat("  |log2FC|:", log2fc_cutoff, "\n\n")
cat("INPUT SUMMARY:\n")
cat("  Total genes:", nrow(results_genes), "\n")
cat("  Significant (q <", qval_cutoff, "):", nrow(sig_genes), "\n")
cat("  Upregulated:", nrow(sig_up), "\n")
cat("  Downregulated:", nrow(sig_down), "\n\n")
sessionInfo()
sink()

# Save workspace
save.image("results/enrichment/enrichment_analysis.RData")

cat("\nDone!\n")
