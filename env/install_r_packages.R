#!/usr/bin/env Rscript
# ============================================================================
# Installation & Verification Script for Compositional FACS Analysis
# ============================================================================
# Purpose: Complete R package installation, verify compatibility, test loading
# Run: Rscript install_r_packages_FIXED.R
# FIXED: CRAN mirror configuration + correct R version expectations
# ============================================================================

cat("\n=== FACS Compositional Analysis - R Environment Setup ===\n\n")

# 0. CONFIGURE CRAN MIRROR (CRITICAL FIX)
cat("Step 0: Configuring CRAN repository...\n")
options(repos = c(CRAN = "https://cloud.r-project.org/"))
cat("  ✓ CRAN mirror set to: https://cloud.r-project.org/\n")

# 1. VERIFY R VERSION & BIOCONDUCTOR COMPATIBILITY
cat("\nStep 1: Checking R version and Bioconductor compatibility...\n")
r_version <- as.character(getRversion())
cat(sprintf("  R version: %s\n", r_version))

# FIXED: R 4.3.x is appropriate for Bioconductor 3.18/3.19
required_r_version <- "4.3.0"
expected_bioc <- "3.18"  # or 3.19 for R 4.3.x

if (utils::compareVersion(r_version, required_r_version) >= 0) {
  cat(sprintf("  ✓ R version %s meets requirement (>= %s)\n", r_version, required_r_version))
} else {
  stop(sprintf("ERROR: R version %s < required %s\n", r_version, required_r_version))
}

# 2. BIOCMANAGER SETUP
cat("\nStep 2: Setting up BiocManager...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("  Installing BiocManager from CRAN...\n")
  install.packages("BiocManager", quiet = TRUE)
}
library(BiocManager)

bioc_version <- as.character(BiocManager::version())
cat(sprintf("  Bioconductor version: %s\n", bioc_version))
cat(sprintf("  Expected for R %s: %s or 3.19\n", r_version, expected_bioc))

# Verify compatibility (accept both 3.18 and 3.19 for R 4.3.x)
if (bioc_version %in% c("3.18", "3.19")) {
  cat("  ✓ Bioconductor version compatible\n")
} else {
  cat(sprintf("  ⚠ WARNING: Bioconductor %s unexpected for R %s\n", bioc_version, r_version))
}

# 3. CORE PACKAGES - ALDEx2 & PROPELLER
cat("\nStep 3: Installing/Verifying CORE Differential Abundance packages...\n")

core_bioc_packages <- c(
  "ALDEx2",
  "speckle",
  "limma",
  "sva"
)

cat("  Installing Bioconductor packages:\n")
for (pkg in core_bioc_packages) {
  cat(sprintf("    - %s...", pkg))
  tryCatch({
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
    }
    cat(" ✓\n")
  }, error = function(e) {
    cat(sprintf(" ✗ ERROR: %s\n", conditionMessage(e)))
  })
}

# Note: zCompositions is CRAN, not Bioconductor
cat("\n  Installing CRAN compositional packages:\n")
cran_comp_packages <- c("zCompositions", "compositions", "Hmisc")
for (pkg in cran_comp_packages) {
  cat(sprintf("    - %s...", pkg))
  tryCatch({
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, quiet = TRUE)
    }
    cat(" ✓\n")
  }, error = function(e) {
    cat(sprintf(" ✗ ERROR: %s\n", conditionMessage(e)))
  })
}

# 4. NETWORK PACKAGES
cat("\nStep 4: Installing/Verifying NETWORK analysis packages...\n")

network_cran_packages <- c(
  "igraph",
  "tidygraph",
  "ggraph"
)

cat("  Installing CRAN packages:\n")
for (pkg in network_cran_packages) {
  cat(sprintf("    - %s...", pkg))
  tryCatch({
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, quiet = TRUE)
    }
    cat(" ✓\n")
  }, error = function(e) {
    cat(sprintf(" ✗ ERROR: %s\n", conditionMessage(e)))
  })
}

# 5. CLUSTERING & SURVIVAL PACKAGES
cat("\nStep 5: Installing/Verifying CLUSTERING & SURVIVAL packages...\n")

# ConsensusClusterPlus is Bioconductor
cat("  Installing Bioconductor clustering packages:\n")
clustering_bioc <- c("ConsensusClusterPlus")
for (pkg in clustering_bioc) {
  cat(sprintf("    - %s...", pkg))
  tryCatch({
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
    }
    cat(" ✓\n")
  }, error = function(e) {
    cat(sprintf(" ✗ ERROR: %s\n", conditionMessage(e)))
  })
}

# WGCNA is CRAN (not Bioconductor!)
cat("\n  Installing CRAN clustering packages:\n")
clustering_cran <- c("WGCNA")
for (pkg in clustering_cran) {
  cat(sprintf("    - %s...", pkg))
  tryCatch({
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, quiet = TRUE)
    }
    cat(" ✓\n")
  }, error = function(e) {
    cat(sprintf(" ✗ ERROR: %s\n", conditionMessage(e)))
  })
}

survival_cran_packages <- c(
  "survival",
  "survminer"
)

cat("\n  Installing CRAN survival packages:\n")
for (pkg in survival_cran_packages) {
  cat(sprintf("    - %s...", pkg))
  tryCatch({
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, quiet = TRUE)
    }
    cat(" ✓\n")
  }, error = function(e) {
    cat(sprintf(" ✗ ERROR: %s\n", conditionMessage(e)))
  })
}

# 6. VISUALIZATION & DATA WRANGLING
cat("\nStep 6: Installing/Verifying VISUALIZATION & DATA packages...\n")

viz_packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "tibble",
  "ggpubr",
  "cowplot",
  "viridis",
  "pheatmap",
  "reshape2",
  "data.table"
)

cat("  Installing CRAN packages:\n")
for (pkg in viz_packages) {
  cat(sprintf("    - %s...", pkg))
  tryCatch({
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, quiet = TRUE)
    }
    cat(" ✓\n")
  }, error = function(e) {
    cat(sprintf(" ✗ ERROR: %s\n", conditionMessage(e)))
  })
}

# 7. UTILITY PACKAGES
cat("\nStep 7: Installing/Verifying UTILITY packages...\n")

utility_packages <- c(
  "rmarkdown",
  "knitr",
  "pwr"  # For power analysis
)

cat("  Installing CRAN packages:\n")
for (pkg in utility_packages) {
  cat(sprintf("    - %s...", pkg))
  tryCatch({
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, quiet = TRUE)
    }
    cat(" ✓\n")
  }, error = function(e) {
    cat(sprintf(" ✗ ERROR: %s\n", conditionMessage(e)))
  })
}

# 8. COMPREHENSIVE COMPATIBILITY TEST
cat("\nStep 8: CRITICAL COMPATIBILITY TEST - Loading all packages...\n\n")

all_packages <- c(
  # Core DA
  "ALDEx2", "speckle", "limma", "compositions",
  # Batch correction
  "sva",
  # Network
  "igraph", "tidygraph", "ggraph", "Hmisc",
  # Clustering
  "ConsensusClusterPlus", "WGCNA",
  # Survival
  "survival", "survminer",
  # Visualization
  "ggplot2", "dplyr", "tidyr", "tibble", "ggpubr", "cowplot", "viridis",
  "pheatmap", "reshape2", "data.table"
)

failed_packages <- character()
loaded_packages <- character()

for (pkg in all_packages) {
  tryCatch({
    suppressPackageStartupMessages(library(pkg, character.only = TRUE, quietly = TRUE))
    cat(sprintf("  ✓ %s\n", pkg))
    loaded_packages <- c(loaded_packages, pkg)
  }, error = function(e) {
    cat(sprintf("  ✗ %s - ERROR: %s\n", pkg, conditionMessage(e)))
    failed_packages <<- c(failed_packages, pkg)
  })
}

# 9. SUMMARY & RECOMMENDATIONS
cat("\n\n=== SUMMARY ===\n")
cat(sprintf("Total packages tested: %d\n", length(all_packages)))
cat(sprintf("Successfully loaded: %d\n", length(loaded_packages)))
cat(sprintf("Failed packages: %d\n", length(failed_packages)))

if (length(failed_packages) > 0) {
  cat("\n⚠ FAILED PACKAGES:\n")
  for (pkg in failed_packages) {
    cat(sprintf("  - %s\n", pkg))
  }
  cat("\nRECOMMENDATIONS:\n")
  cat("  1. Check internet connection\n")
  cat("  2. Verify micromamba environment is activated\n")
  cat("  3. For specific package failures, run:\n")
  
  bioc_pkgs <- c("ALDEx2", "speckle", "limma", "sva", "ConsensusClusterPlus")
  for (pkg in failed_packages) {
    if (pkg %in% bioc_pkgs) {
      cat(sprintf("     BiocManager::install('%s', force = TRUE)\n", pkg))
    } else {
      cat(sprintf("     install.packages('%s')\n", pkg))
    }
  }
  
  cat("\n  4. If persistent errors, check conda environment:\n")
  cat("     micromamba list | grep 'bioconductor\\|r-'\n")
  
} else {
  cat("\n✓ ALL PACKAGES SUCCESSFULLY LOADED!\n")
  cat("\n🎯 Environment is ready for FACS compositional analysis:\n")
  cat("  ✓ ALDEx2 + Propeller DA testing (with scale models γ=0.5)\n")
  cat("  ✓ ComBat batch correction (sva package)\n")
  cat("  ✓ Spearman network analysis (igraph)\n")
  cat("  ✓ Consensus clustering (pan-cancer subtypes)\n")
  cat("  ✓ Survival analysis (Kaplan-Meier, Cox regression)\n")
  cat("  ✓ Full ggplot2/tidyverse visualization stack\n")
  cat("  ✓ WGCNA for module detection\n")
}

# 10. VERSION INFORMATION FOR REPRODUCIBILITY
cat("\n=== INSTALLED PACKAGE VERSIONS (for reproducibility) ===\n")
key_packages_versions <- c(
  "ALDEx2", "speckle", "sva", "igraph", "ConsensusClusterPlus",
  "WGCNA", "survival", "ggplot2", "dplyr", "compositions"
)

version_info <- data.frame(
  Package = character(),
  Version = character(),
  stringsAsFactors = FALSE
)

for (pkg in key_packages_versions) {
  tryCatch({
    ver <- as.character(packageVersion(pkg))
    version_info <- rbind(version_info, data.frame(Package = pkg, Version = ver))
  }, error = function(e) {
    version_info <<- rbind(version_info, data.frame(Package = pkg, Version = "NOT LOADED"))
  })
}

print(version_info, row.names = FALSE)

cat("\n=== SESSION INFO ===\n")
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("Bioconductor version: %s\n", BiocManager::version()))
cat(sprintf("Platform: %s\n", R.version$platform))
cat(sprintf("Date: %s\n", Sys.Date()))


# 11. WORKFLOW-SPECIFIC RECOMMENDATIONS
cat("\n=== WORKFLOW-SPECIFIC NOTES (from validation PDF) ===\n")
cat("📋 Critical reminders for analysis:\n\n")
cat("1. ALDEx2 Scale Models:\n")
cat("   - Use γ=0.5 (2024 update) instead of legacy γ=0:\n")
cat("     aldex_result <- aldex(data, conditions, mc.samples=128, gamma=0.5)\n\n")

cat("2. Batch Correction:\n")
cat("   - Apply ComBat to CLR-transformed data (NOT raw proportions):\n")
cat("     clr_data <- compositions::clr(proportions)\n")
cat("     corrected <- sva::ComBat(clr_data, batch=batch_var)\n\n")

cat("3. Network Construction:\n")
cat("   - Threshold: |ρ| > 0.6, p < 0.05 (validated for immune networks)\n")
cat("   - Bootstrap edges: B=1000, report only support > 70%\n\n")

cat("4. Sample Size (n=10-15 per group):\n")
cat("   - Use ALDEx2 + Propeller validation (intersection = high confidence)\n")
cat("   - Consider scCODA if rare cell types (<1%) are critical\n\n")

cat("5. Power Analysis:\n")
cat("   - Expected effect size > 1 (log scale) for reliable detection\n")
cat("   - Pre-study power analysis recommended (pwr package)\n\n")

cat("✓ Setup complete! Ready to analyze FACS data following the validated workflow.\n")
cat("📄 See workflow-validazione-composizionale-FACS.pdf for detailed methodology.\n\n")

# ============================================================================
# END INSTALLATION SCRIPT - FIXED VERSION
# ============================================================================