# ============================================================================
# FACS Compositional Analysis - Modular Step-by-Step Pipeline
# Quality Control, Exploration, CLR, PCA, Differential Abundance
# Designed for immunological FACS data with flexible NA handling
# ============================================================================

library(tidyverse)
library(ALDEx2)
library(pheatmap)
library(corrplot)
library(factoextra)

# ============================================================================
# STEP 1: DATA IMPORT AND INITIAL VALIDATION
# ============================================================================

step1_import_facs <- function(abundance, sample_ids = NULL) {
  
  cat("\n===== STEP 1: Import and Initial Validation =====\n\n")
  
  # Convert to dataframe if needed
  abundance <- as.data.frame(abundance)
  
  # Set rownames if sample_ids provided
  if (!is.null(sample_ids)) {
    rownames(abundance) <- sample_ids
  } else if (!"Patient_ID" %in% colnames(abundance) && is.null(rownames(abundance))) {
    rownames(abundance) <- paste0("Sample_", 1:nrow(abundance))
  }
  
  # Remove Patient_ID column if present (move to rownames)
  if ("Patient_ID" %in% colnames(abundance)) {
    rownames(abundance) <- abundance$Patient_ID
    abundance <- abundance %>% select(-Patient_ID)
  }
  
  # Convert all to numeric
  abundance <- abundance %>% mutate(across(everything(), as.numeric))
  
  # Basic statistics
  cat(sprintf("Samples: %d\n", nrow(abundance)))
  cat(sprintf("Populations: %d\n", ncol(abundance)))
  cat(sprintf("Total cells measured: %d\n", nrow(abundance) * ncol(abundance)))
  cat(sprintf("Missing values: %d (%.1f%%)\n", 
              sum(is.na(abundance)), 
              100 * sum(is.na(abundance)) / (nrow(abundance) * ncol(abundance))))
  
  # Range check
  range_min <- min(abundance, na.rm = TRUE)
  range_max <- max(abundance, na.rm = TRUE)
  cat(sprintf("Value range: [%.2f, %.2f]\n", range_min, range_max))
  
  if (range_min < 0) stop("ERROR: Negative values detected")
  if (range_max > 100) warning("WARNING: Values exceed 100% detected")
  
  # Compositional sums per sample
  sums <- rowSums(abundance, na.rm = TRUE)
  cat(sprintf("\nCompositional sums: [%.1f, %.1f]\n", min(sums), max(sums)))
  
  if (any(sums < 80 | sums > 120)) {
    cat("⚠ WARNING: Some samples have unusual compositional sums\n")
    cat("  This suggests incomplete gating or missing populations\n")
  }
  
  return(abundance)
}

# ============================================================================
# STEP 2: EXPLORATORY DATA ANALYSIS (EDA)
# ============================================================================

step2_explore_facs <- function(abundance, output_dir = NULL) {
  
  cat("\n===== STEP 2: Exploratory Data Analysis =====\n\n")
  
  # 2.1: Missing data pattern
  cat("2.1 Missing Data Pattern:\n")
  missing_by_pop <- colSums(is.na(abundance)) / nrow(abundance) * 100
  missing_by_sample <- rowSums(is.na(abundance)) / ncol(abundance) * 100
  
  cat(sprintf("  Populations with >50%% missing: %d\n", 
              sum(missing_by_pop > 50)))
  cat(sprintf("  Samples with >50%% missing: %d\n", 
              sum(missing_by_sample > 50)))
  
  if (sum(missing_by_pop > 50) > 0) {
    cat("\n  Populations with high missingness:\n")
    high_miss <- names(missing_by_pop[missing_by_pop > 50])
    for (pop in high_miss) {
      cat(sprintf("    %s: %.1f%%\n", pop, missing_by_pop[pop]))
    }
  }
  
  # 2.2: Population abundance distribution
  cat("\n2.2 Population Abundance Distribution:\n")
  medians <- apply(abundance, 2, median, na.rm = TRUE)
  rare_pops <- names(medians[medians < 0.5])
  abundant_pops <- names(medians[medians > 10])
  
  cat(sprintf("  Rare populations (<0.5%%): %d\n", length(rare_pops)))
  cat(sprintf("  Abundant populations (>10%%): %d\n", length(abundant_pops)))
  
  if (length(rare_pops) > 0 && length(rare_pops) <= 10) {
    cat("    Rare:", paste(rare_pops, collapse = ", "), "\n")
  }
  
  # 2.3: Identify potential outlier samples
  cat("\n2.3 Outlier Detection:\n")
  # Use compositional sum as simple outlier metric
  sums <- rowSums(abundance, na.rm = TRUE)
  outliers <- which(sums < quantile(sums, 0.25) - 1.5*IQR(sums) | 
                    sums > quantile(sums, 0.75) + 1.5*IQR(sums))
  
  if (length(outliers) > 0) {
    cat(sprintf("  Potential outlier samples: %s\n", 
                paste(rownames(abundance)[outliers], collapse = ", ")))
  } else {
    cat("  No clear outliers detected\n")
  }
  
  # 2.4: Correlation structure (only for populations with <50% missing)
  complete_pops <- names(missing_by_pop[missing_by_pop < 50])
  if (length(complete_pops) > 2) {
    cat("\n2.4 Computing correlations between populations...\n")
    
    # Use pairwise complete observations
    cor_mat <- cor(abundance[, complete_pops], 
                   use = "pairwise.complete.obs", 
                   method = "spearman")
    
    # Find highly correlated pairs
    high_cor <- which(abs(cor_mat) > 0.7 & cor_mat != 1, arr.ind = TRUE)
    if (nrow(high_cor) > 0) {
      cat("  High correlations detected (|r| > 0.7):\n")
      printed <- 0
      for (i in 1:nrow(high_cor)) {
        if (high_cor[i, 1] < high_cor[i, 2] && printed < 5) {  # Avoid duplicates
          cat(sprintf("    %s ↔ %s: r = %.2f\n", 
                      rownames(cor_mat)[high_cor[i, 1]],
                      colnames(cor_mat)[high_cor[i, 2]],
                      cor_mat[high_cor[i, 1], high_cor[i, 2]]))
          printed <- printed + 1
        }
      }
    }
  }
  
  # Store results
  eda_results <- list(
    missing_by_pop = missing_by_pop,
    missing_by_sample = missing_by_sample,
    medians = medians,
    rare_pops = rare_pops,
    abundant_pops = abundant_pops,
    outliers = outliers,
    correlation_matrix = if (exists("cor_mat")) cor_mat else NULL
  )
  
  cat("\n✓ Exploratory analysis complete\n")
  
  return(eda_results)
}

# ============================================================================
# STEP 3: QUALITY CONTROL AND FILTERING
# ============================================================================

step3_qc_facs <- function(abundance, 
                          pop_missing_threshold = 0.6,
                          sample_missing_threshold = 0.6,
                          min_samples = 5,
                          impute_method = "geometric_mean") {
  
  cat("\n===== STEP 3: Quality Control and Filtering =====\n\n")
  
  initial_samples <- nrow(abundance)
  initial_pops <- ncol(abundance)
  
  # 3.1: Remove populations with excessive missingness
  missing_by_pop <- colSums(is.na(abundance)) / nrow(abundance)
  pops_to_remove <- names(missing_by_pop[missing_by_pop > pop_missing_threshold])
  
  if (length(pops_to_remove) > 0) {
    cat(sprintf("Removing %d populations with >%.0f%% missing:\n", 
                length(pops_to_remove), pop_missing_threshold * 100))
    cat(sprintf("  %s\n", paste(pops_to_remove, collapse = ", ")))
    abundance <- abundance %>% select(-all_of(pops_to_remove))
  }
  
  # 3.2: Remove samples with excessive missingness
  missing_by_sample <- rowSums(is.na(abundance)) / ncol(abundance)
  samples_to_remove <- which(missing_by_sample > sample_missing_threshold)
  
  if (length(samples_to_remove) > 0) {
    cat(sprintf("Removing %d samples with >%.0f%% missing:\n", 
                length(samples_to_remove), sample_missing_threshold * 100))
    cat(sprintf("  %s\n", paste(rownames(abundance)[samples_to_remove], collapse = ", ")))
    abundance <- abundance[-samples_to_remove, ]
  }
  
  # 3.3: Check minimum sample requirement
  if (nrow(abundance) < min_samples) {
    stop(sprintf("ERROR: Only %d samples remain after filtering (minimum: %d)", 
                 nrow(abundance), min_samples))
  }
  
  # 3.4: Impute remaining NA values
  if (any(is.na(abundance))) {
    cat(sprintf("\nImputing %d remaining NA values using %s...\n", 
                sum(is.na(abundance)), impute_method))
    
    if (impute_method == "geometric_mean") {
      # Impute with geometric mean of non-NA values in that population
      for (col in names(abundance)) {
        na_idx <- is.na(abundance[[col]])
        if (any(na_idx)) {
          non_na_vals <- abundance[[col]][!na_idx]
          # Geometric mean (add small pseudocount for zeros)
          geom_mean <- exp(mean(log(non_na_vals + 0.001)))
          abundance[[col]][na_idx] <- geom_mean
        }
      }
    } else if (impute_method == "median") {
      for (col in names(abundance)) {
        na_idx <- is.na(abundance[[col]])
        if (any(na_idx)) {
          abundance[[col]][na_idx] <- median(abundance[[col]], na.rm = TRUE)
        }
      }
    }
  }
  
  # 3.5: Compositional closure (renormalize to 100%)
  cat("\nApplying compositional closure (renormalizing to 100%)...\n")
  row_sums <- rowSums(abundance)
  abundance_closed <- abundance / row_sums * 100
  
  # Verify closure
  new_sums <- rowSums(abundance_closed)
  if (max(abs(new_sums - 100)) > 0.01) {
    warning("Closure verification failed")
  }
  
  # QC Report
  cat("\n--- QC Summary ---\n")
  cat(sprintf("Samples: %d → %d (removed: %d)\n", 
              initial_samples, nrow(abundance_closed), 
              initial_samples - nrow(abundance_closed)))
  cat(sprintf("Populations: %d → %d (removed: %d)\n", 
              initial_pops, ncol(abundance_closed), 
              initial_pops - ncol(abundance_closed)))
  cat(sprintf("Missing values: 0 (%.1f%% imputed)\n", 
              100 * sum(is.na(abundance)) / (initial_samples * initial_pops)))
  
  qc_report <- list(
    n_samples_initial = initial_samples,
    n_samples_final = nrow(abundance_closed),
    n_populations_initial = initial_pops,
    n_populations_final = ncol(abundance_closed),
    populations_removed = pops_to_remove,
    samples_removed = if (length(samples_to_remove) > 0) 
      rownames(abundance)[samples_to_remove] else character(0)
  )
  
  cat("\n✓ Quality control complete\n")
  
  return(list(data = abundance_closed, report = qc_report))
}

# ============================================================================
# STEP 4: CLR TRANSFORMATION
# ============================================================================

step4_clr_transform <- function(abundance, pseudocount = 0.01) {
  
  cat("\n===== STEP 4: CLR Transformation =====\n\n")
  
  # Add pseudocount (smaller for already-closed data)
  abundance_pc <- abundance + pseudocount
  
  # CLR: log(x / geometric mean)
  clr_data <- t(apply(abundance_pc, 1, function(row) {
    gm <- exp(mean(log(row)))
    log(row / gm)
  }))
  
  colnames(clr_data) <- colnames(abundance)
  rownames(clr_data) <- rownames(abundance)
  
  # Verify centering
  row_means <- rowMeans(clr_data)
  if (max(abs(row_means)) > 1e-6) {
    warning("CLR not properly centered (row means != 0)")
  }
  
  cat(sprintf("CLR transformation applied with pseudocount = %.3f\n", pseudocount))
  cat(sprintf("CLR value range: [%.2f, %.2f]\n", 
              min(clr_data), max(clr_data)))
  
  cat("\n✓ CLR transformation complete\n")
  
  return(clr_data)
}

# ============================================================================
# STEP 5: PRINCIPAL COMPONENT ANALYSIS
# ============================================================================

step5_pca_facs <- function(clr_data, plot = TRUE) {
  
  cat("\n===== STEP 5: Principal Component Analysis =====\n\n")
  
  # PCA with scaling (standard for CLR data)
  pca <- prcomp(clr_data, scale. = TRUE, center = TRUE)
  
  # Variance explained
  var_exp <- summary(pca)$importance[2, ]
  cum_var <- cumsum(var_exp)
  
  cat("Variance Explained:\n")
  cat(sprintf("  PC1: %.1f%%\n", var_exp[1] * 100))
  cat(sprintf("  PC2: %.1f%%\n", var_exp[2] * 100))
  cat(sprintf("  PC3: %.1f%%\n", var_exp[3] * 100))
  cat(sprintf("  PC1-PC3 cumulative: %.1f%%\n", cum_var[3] * 100))
  cat(sprintf("  PCs needed for 90%% variance: %d\n", 
              min(which(cum_var >= 0.90))))
  
  # Top contributing populations
  cat("\nTop 5 populations contributing to PC1:\n")
  pc1_loadings <- abs(pca$rotation[, 1])
  top_pc1 <- names(sort(pc1_loadings, decreasing = TRUE)[1:5])
  for (pop in top_pc1) {
    cat(sprintf("  %s: %.3f\n", pop, pca$rotation[pop, 1]))
  }
  
  cat("\n✓ PCA complete\n")
  
  return(pca)
}

# ============================================================================
# STEP 6: DIFFERENTIAL ABUNDANCE (OPTIONAL - requires conditions)
# ============================================================================

step6_da_aldex2 <- function(abundance, 
                            conditions, 
                            mc_samples = 128, 
                            fdr_threshold = 0.1) {
  
  cat("\n===== STEP 6: Differential Abundance (ALDEx2) =====\n\n")
  
  # Validate conditions
  if (length(conditions) != nrow(abundance)) {
    stop("Length of conditions must match number of samples")
  }
  
  conditions <- as.factor(conditions)
  if (length(levels(conditions)) != 2) {
    stop("conditions must have exactly 2 levels for t-test")
  }
  
  cat(sprintf("Comparing: %s (n=%d) vs %s (n=%d)\n",
              levels(conditions)[1], sum(conditions == levels(conditions)[1]),
              levels(conditions)[2], sum(conditions == levels(conditions)[2])))
  
  # Convert to pseudo-counts for ALDEx2
  # Scale proportions to avoid integer rounding issues
  counts <- round(abundance * 100)
  counts <- as.matrix(counts)
  mode(counts) <- "integer"
  
  # ALDEx2 analysis
  cat("\nRunning ALDEx2 (this may take a minute)...\n")
  
  aldex_obj <- aldex(counts,
                     conditions = conditions,
                     test = "t",
                     effect = TRUE,
                     denom = "iqlr",  # Robust to outliers
                     mc.samples = mc_samples,
                     verbose = FALSE)
  
  # Format results
  results <- aldex_obj %>%
    rownames_to_column("Population") %>%
    arrange(we.eBH) %>%
    mutate(
      Significant = we.eBH < fdr_threshold,
      EffectSize = effect,
      Direction = ifelse(effect > 0, 
                        paste0("↑ in ", levels(conditions)[2]), 
                        paste0("↑ in ", levels(conditions)[1]))
    )
  
  # Summary
  n_sig <- sum(results$Significant)
  cat(sprintf("\n--- Differential Abundance Summary ---\n"))
  cat(sprintf("Significant populations (FDR < %.2f): %d / %d\n", 
              fdr_threshold, n_sig, nrow(results)))
  
  if (n_sig > 0) {
    cat("\nSignificant populations:\n")
    sig_results <- results %>% 
      filter(Significant) %>% 
      select(Population, EffectSize, we.eBH, Direction)
    print(sig_results, row.names = FALSE)
  } else {
    cat("\nNo significant differences detected.\n")
    cat("Top 5 populations by effect size:\n")
    top_results <- results %>% 
      slice(1:5) %>% 
      select(Population, EffectSize, we.eBH, Direction)
    print(top_results, row.names = FALSE)
  }
  
  cat("\n✓ Differential abundance analysis complete\n")
  
  return(results)
}

# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

plot_missing_pattern <- function(abundance) {
  missing_mat <- is.na(abundance) * 1
  pheatmap(missing_mat, 
           cluster_rows = TRUE, 
           cluster_cols = TRUE,
           color = c("lightblue", "darkred"),
           main = "Missing Data Pattern",
           legend_labels = c("Present", "Missing"))
}

plot_abundance_heatmap <- function(abundance, scale_method = "column") {
  pheatmap(t(abundance), 
           scale = scale_method,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "euclidean",
           main = "Population Abundance Heatmap",
           fontsize_row = 8)
}

plot_pca_biplot <- function(pca, conditions = NULL) {
  if (!is.null(conditions)) {
    fviz_pca_biplot(pca, 
                    geom.ind = "point",
                    col.ind = conditions,
                    addEllipses = TRUE,
                    repel = TRUE,
                    title = "PCA Biplot - Samples and Populations")
  } else {
    fviz_pca_biplot(pca, 
                    geom.ind = "point",
                    repel = TRUE,
                    title = "PCA Biplot - Samples and Populations")
  }
}

plot_da_volcano <- function(da_results, fdr_threshold = 0.1) {
  da_results %>%
    ggplot(aes(x = EffectSize, y = -log10(we.eBH))) +
    geom_point(aes(color = Significant), alpha = 0.7, size = 3) +
    geom_hline(yintercept = -log10(fdr_threshold), 
               linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-1, 1), 
               linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("gray50", "firebrick")) +
    labs(title = "Differential Abundance - Volcano Plot",
         x = "Effect Size",
         y = "-log10(FDR-adjusted p-value)") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

# ============================================================================
# COMPLETE PIPELINE WRAPPER (OPTIONAL - for bulk execution)
# ============================================================================

run_complete_pipeline <- function(abundance, 
                                  conditions = NULL,
                                  pop_missing_threshold = 0.6,
                                  sample_missing_threshold = 0.6,
                                  mc_samples = 128,
                                  fdr_threshold = 0.1) {
  
  cat("\n╔════════════════════════════════════════════════════╗\n")
  cat("║   FACS Compositional Analysis - Complete Pipeline ║\n")
  cat("╚════════════════════════════════════════════════════╝\n")
  
  # Step 1: Import
  abundance <- step1_import_facs(abundance)
  
  # Step 2: Explore
  eda_results <- step2_explore_facs(abundance)
  
  # Step 3: QC
  qc_result <- step3_qc_facs(abundance, 
                              pop_missing_threshold = pop_missing_threshold,
                              sample_missing_threshold = sample_missing_threshold)
  abundance_clean <- qc_result$data
  
  # Step 4: CLR
  clr_data <- step4_clr_transform(abundance_clean)
  
  # Step 5: PCA
  pca_result <- step5_pca_facs(clr_data)
  
  # Step 6: DA (optional)
  da_result <- NULL
  if (!is.null(conditions)) {
    # Match conditions to remaining samples
    cond_matched <- conditions[rownames(abundance_clean)]
    da_result <- step6_da_aldex2(abundance_clean, cond_matched, 
                                 mc_samples, fdr_threshold)
  } else {
    cat("\n⚠ Skipping differential abundance (no conditions provided)\n")
  }
  
  cat("\n╔════════════════════════════════════════════════════╗\n")
  cat("║              Pipeline Complete ✓                   ║\n")
  cat("╚════════════════════════════════════════════════════╝\n")
  
  return(list(
    abundance_clean = abundance_clean,
    clr_data = clr_data,
    pca = pca_result,
    da_results = da_result,
    qc_report = qc_result$report,
    eda = eda_results
  ))
}

# ============================================================================
# USAGE EXAMPLES
# ============================================================================

# # STEP-BY-STEP EXECUTION (RECOMMENDED)
# 
# # Load your data
# abundance_raw <- read.csv("your_data.csv", row.names = 1)
# 
# # Step 1: Import and validate
# abundance <- step1_import_facs(abundance_raw)
# 
# # Step 2: Explore (look at results before proceeding)
# eda <- step2_explore_facs(abundance)
# 
# # Step 3: QC and filter (adjust thresholds as needed)
# qc <- step3_qc_facs(abundance, 
#                     pop_missing_threshold = 0.6,
#                     sample_missing_threshold = 0.6)
# abundance_clean <- qc$data
# 
# # Step 4: CLR transform
# clr_data <- step4_clr_transform(abundance_clean)
# 
# # Step 5: PCA (works without conditions)
# pca <- step5_pca_facs(clr_data)
# plot_pca_biplot(pca)  # Visualize
# 
# # Step 6: DA (OPTIONAL - only if you have conditions)
# # conditions <- factor(c("Control", "Control", "Treatment", "Treatment", ...))
# # da <- step6_da_aldex2(abundance_clean, conditions)
# # plot_da_volcano(da)
# 
# # OR USE COMPLETE PIPELINE
# results <- run_complete_pipeline(abundance_raw, conditions = NULL)
