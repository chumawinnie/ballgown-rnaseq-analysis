#===============================================================================
# INSTALL REQUIRED PACKAGES
# Run this script once before using the pipeline
#===============================================================================

cat("Installing required packages for RNA-seq DE Analysis Pipeline...\n\n")

#-------------------------------------------------------------------------------
# 1. Install BiocManager if not present
#-------------------------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

#-------------------------------------------------------------------------------
# 2. Bioconductor Packages
#-------------------------------------------------------------------------------

bioc_packages <- c(
  "ballgown",
  "genefilter",
  "biomaRt",
  "topGO",
  "clusterProfiler",
  "org.Hs.eg.db",
  "DOSE",
  "enrichplot",
  "pathview",
  "EnhancedVolcano",
  "AnnotationDbi"
)

cat("Installing Bioconductor packages...\n")
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("  Installing:", pkg, "\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat("  Already installed:", pkg, "\n")
  }
}

#-------------------------------------------------------------------------------
# 3. CRAN Packages
#-------------------------------------------------------------------------------

cran_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "pheatmap",
  "RColorBrewer",
  "ggrepel",
  "ggupset"
)

cat("\nInstalling CRAN packages...\n")
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("  Installing:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  } else {
    cat("  Already installed:", pkg, "\n")
  }
}

#-------------------------------------------------------------------------------
# 4. Verify Installation
#-------------------------------------------------------------------------------

cat("\n")
cat("================================================================\n")
cat("VERIFICATION\n")
cat("================================================================\n")

all_packages <- c(bioc_packages, cran_packages)
installed <- sapply(all_packages, function(pkg) {
  require(pkg, character.only = TRUE, quietly = TRUE)
})

if (all(installed)) {
  cat("\n✓ All packages installed successfully!\n")
  cat("\nYou can now run the pipeline:\n")
  cat("  1. source('ballgown_corrected.R')\n")
  cat("  2. source('functional_enrichment_qvalue.R')\n")
} else {
  cat("\n✗ Some packages failed to install:\n")
  print(names(installed)[!installed])
  cat("\nPlease install these manually.\n")
}

#-------------------------------------------------------------------------------
# 5. Session Info
#-------------------------------------------------------------------------------

cat("\n")
cat("================================================================\n")
cat("SESSION INFO\n")
cat("================================================================\n")
sessionInfo()
