# ============================================================================
# FACS Compositional Analysis - Compact Pipeline
# Quality Control, Differential Abundance, PCA
# Input: list(abundance = df, metadata = df, conditions = factor)
# ============================================================================

library(tidyverse)
library(ALDEx2)
library(factoextra)

# ============================================================================
# QC_FACS: Quality Control for Compositional FACS Data
# ============================================================================
# Input: abundance dataframe (samples × populations, proportions 0-100)
# Output: list(clean_data, qc_report)

qc_facs <- function(abundance, threshold_missing = 0.4) {
  
  # Range check: [0, 100]
  if (min(abundance, na.rm = TRUE) < 0 | max(abundance, na.rm = TRUE) > 100) {
    stop("Abundance values outside [0, 100] range")
  }
  
  # Compositional sum check
  sums <- rowSums(abundance, na.rm = TRUE)
  if (any(sums < 80 | sums > 120)) {
    warning("Some samples have unusual compositional sums (check missing data)")
  }
  
  # Identify populations with excessive missingness
  missing_rate <- colSums(is.na(abundance)) / nrow(abundance)
  populations_drop <- names(missing_rate[missing_rate > threshold_missing])
  
  # Strict: remove populations and samples with any NA
  clean_data <- abundance %>%
    select(-all_of(populations_drop)) %>%
    filter(complete.cases(.))
  
  # Report
  qc_report <- list(
    n_samples_initial = nrow(abundance),
    n_samples_final = nrow(clean_data),
    n_populations_initial = ncol(abundance),
    n_populations_final = ncol(clean_data),
    populations_removed = populations_drop,
    samples_removed = rownames(abundance)[!rownames(abundance) %in% rownames(clean_data)]
  )
  
  cat("QC Report:\n")
  cat(sprintf("  Samples: %d → %d\n", qc_report$n_samples_initial, qc_report$n_samples_final))
  cat(sprintf("  Populations: %d → %d\n", qc_report$n_populations_initial, qc_report$n_populations_final))
  if (length(populations_drop) > 0) cat(sprintf("  Removed: %s\n", paste(populations_drop, collapse=", ")))
  
  return(list(data = clean_data, report = qc_report))
}

# ============================================================================
# CLR_TRANSFORM: Centered Log-Ratio Transformation
# ============================================================================
# For compositional data: CLR(x) = log(x / geometric_mean(x))
# Applied to abundance data (samples × populations)

clr_transform <- function(abundance, pseudocount = 0.1) {
  
  # Add pseudocount
  abundance_pc <- abundance + pseudocount
  
  # CLR: log(x / geometric mean)
  clr_data <- t(apply(abundance_pc, 1, function(row) {
    gm <- exp(mean(log(row)))
    return(log(row / gm))
  }))
  
  colnames(clr_data) <- colnames(abundance)
  rownames(clr_data) <- rownames(abundance)
  
  # Verify CLR centering: column means should be ≈0
  col_means <- colMeans(clr_data)
  if (max(abs(col_means)) > 1e-6) warning("CLR not properly centered")
  
  return(clr_data)
}

# ============================================================================
# PCA_FACS: PCA on CLR-Transformed Data
# ============================================================================
# Input: CLR-transformed abundance matrix
# Output: PCA object + variance explained

pca_facs <- function(clr_data) {
  
  # Handle NA: keep only complete cases
  clr_complete <- clr_data[complete.cases(clr_data), ]
  
  if (nrow(clr_complete) < nrow(clr_data)) {
    cat(sprintf("Dropped %d samples with incomplete data\n", nrow(clr_data) - nrow(clr_complete)))
  }
  
  # PCA with scaling
  pca <- prcomp(clr_complete, scale = TRUE, center = TRUE)
  
  # Variance explained
  var_exp <- summary(pca)$importance[2, ]
  cum_var <- cumsum(var_exp)
  
  cat("PCA Variance Explained:\n")
  cat(sprintf("  PC1-PC3: %.1f%%\n", cum_var[3] * 100))
  cat(sprintf("  PCs for 90%% variance: %d\n", min(which(cum_var >= 0.90))))
  
  return(pca)
}

# ============================================================================
# DA_FACS: ALDEx2 Differential Abundance Testing
# ============================================================================
# Input: abundance (proportions), conditions (factor with 2 levels)
# Output: ALDEx2 results with effect sizes, p-values

da_facs <- function(abundance, conditions, mc_samples = 128, fdr_threshold = 0.1) {
  
  # Convert proportions → pseudo-counts (standard for ALDEx2)
  # Scale: multiply by 10000 to avoid rounding issues
  counts <- round(abundance * 100)  # FACS proportions → scale 0-100
  counts <- round(counts * 100)  # Scale up for ALDEx2
  
  # Ensure integer matrix
  counts <- as.matrix(counts)
  mode(counts) <- "integer"
  
  # ALDEx2: robust DA testing
  # denom="iqlr": inter-quartile log-ratio (robust to outliers for small n)
  cat("Running ALDEx2 differential abundance...\n")
  
  aldex_obj <- aldex(counts,
                     conditions = conditions,
                     test = "t",
                     effect = TRUE,
                     denom = "iqlr",
                     mc.samples = mc_samples,
                     verbose = FALSE)
  
  # Extract and format results
  results <- aldex_obj %>%
    rownames_to_column("Population") %>%
    arrange(we.eBH) %>%
    mutate(
      Significant = we.eBH < fdr_threshold,
      EffectSize = effect
    )
  
  # Report
  n_sig <- sum(results$Significant)
  cat(sprintf("Significant populations (FDR < %.2f): %d\n", fdr_threshold, n_sig))
  
  if (n_sig > 0) {
    cat("\nTop significant populations:\n")
    print(results %>% filter(Significant) %>% select(Population, EffectSize, we.eBH))
  } else {
    cat("\nTop populations by effect size:\n")
    print(results %>% slice(1:5) %>% select(Population, EffectSize, we.eBH))
  }
  
  return(results)
}

# ============================================================================
# MAIN ANALYSIS WRAPPER
# ============================================================================

run_facs_analysis <- function(input_list) {
  
  # Input: list(abundance = df, metadata = df, conditions = factor)
  abundance <- input_list$abundance
  metadata <- input_list$metadata
  conditions <- input_list$conditions
  
  cat("\n===== FACS Compositional Analysis Pipeline =====\n\n")
  
  # Step 1: Quality Control
  cat("Step 1: Quality Control\n")
  qc_result <- qc_facs(abundance)
  abundance_clean <- qc_result$data
  
  # Step 2: CLR Transformation
  cat("\nStep 2: CLR Transformation\n")
  clr_data <- clr_transform(abundance_clean, pseudocount = 0.1)
  cat("✓ CLR transformation complete\n")
  
  # Step 3: PCA
  cat("\nStep 3: Principal Component Analysis\n")
  pca_result <- pca_facs(clr_data)
  cat("✓ PCA complete\n")
  
  # Step 4: Differential Abundance (ALDEx2)
  cat("\nStep 4: Differential Abundance Analysis (ALDEx2)\n")
  
  # Match conditions to clean samples
  cond_matched <- conditions[rownames(abundance_clean)]
  
  da_result <- da_facs(abundance_clean, cond_matched, mc_samples = 128, fdr_threshold = 0.1)
  
  cat("\n✓ Analysis complete\n")
  
  # Return results
  return(list(
    abundance_clean = abundance_clean,
    clr_data = clr_data,
    pca = pca_result,
    da_results = da_result,
    qc_report = qc_result$report
  ))
}

# ============================================================================
# USAGE EXAMPLE
# ============================================================================

# # Prepare input: list with abundance, metadata, conditions
# input_facs <- list(
#   abundance = abundance_df,      # samples × populations, proportions [0-100]
#   metadata = metadata_df,         # sample metadata
#   conditions = factor_conditions  # factor with 2 levels for DA comparison
# )
#
# # Run analysis
# results <- run_facs_analysis(input_facs)
#
# # Extract components
# clr_transformed <- results$clr_data
# pca_obj <- results$pca
# da_results <- results$da_results
