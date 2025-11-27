# ==============================================================================
# PROJECT: Pan-Cancer Differential Network Analysis
# PURPOSE: Compare immune cell composition networks between cancer and healthy
# METHOD: Network Comparison Test (NCT) with compositional data transformation
# ==============================================================================

# ------------------------------------------------------------------------------
# LIBRARIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(qgraph)
  library(NetworkComparisonTest)
  library(igraph)
  library(compositions)
  library(reshape2)
})

# ------------------------------------------------------------------------------
# CONFIGURATION
# ------------------------------------------------------------------------------
config <- list(
  # Paths
  data_path = "/mnt/c/Users/Davide/Desktop/Davide/bandi/0dottorato/projects/long survivors/DB_anonimo_standardizzato.xlsx",
  output_dir = paste0("Analysis_PanCancer_", format(Sys.time(), "%Y%m%d_%H%M")),
  
  # Group definitions
  sheets = list(
    nsclc = "NSCLC",
    hnscc = "HNSCC",
    healthy = "Healthy_Donors"
  ),
  
  # Data cleaning thresholds
  min_completeness = 0.40,      # Minimum 40% non-missing in combined tumor group
  clr_offset = 0.01,            # Offset for zero handling in CLR transform
  
  # Network parameters
  nct_iterations = 1000,
  nct_gamma = 0,                # BIC selection (gamma=0) for dense networks
  correlation_method = "spearman",
  
  # Filtering thresholds
  min_delta_r = 0.30,           # Minimum correlation difference to report
  alpha_level = 0.05,
  fdr_method = "fdr"            # Multiple testing correction method
)

# Create output directory
if (!dir.exists(config$output_dir)) {
  dir.create(config$output_dir, recursive = TRUE)
}

cat("=== Pan-Cancer Network Analysis Pipeline ===\n")
cat("Output directory:", config$output_dir, "\n\n")

# ------------------------------------------------------------------------------
# HELPER FUNCTIONS
# ------------------------------------------------------------------------------

#' Load Excel sheet and coerce to numeric matrix
#' @param path File path to Excel workbook
#' @param sheet Sheet name
#' @return Data frame with Patient_ID and numeric marker columns
load_excel_numeric <- function(path, sheet) {
  tryCatch({
    df <- read_excel(path, sheet = sheet)
    patient_ids <- df[[1]]
    
    # Force numeric conversion (suppress warnings for intended NAs)
    marker_data <- df[, -1] %>%
      lapply(function(x) suppressWarnings(as.numeric(as.character(x)))) %>%
      as.data.frame()
    
    marker_data$patient_id <- patient_ids
    return(marker_data)
  }, error = function(e) {
    stop(sprintf("Failed to load sheet '%s': %s", sheet, e$message))
  })
}

#' Calculate completeness rate for each marker
#' @param df Data frame
#' @return Named vector of completeness rates (0-1)
get_completeness <- function(df) {
  colMeans(!is.na(df))
}

#' Impute missing values using reference distribution or target median
#' @param target_df Target data frame with missing values
#' @param reference_df Reference data frame for distribution-based imputation
#' @param markers Character vector of marker names
#' @return Imputed data frame
impute_with_reference <- function(target_df, reference_df, markers) {
  imputed_df <- target_df
  
  for (marker in markers) {
    target_vec <- target_df[[marker]]
    ref_vec <- reference_df[[marker]]
    missing_idx <- is.na(target_vec)
    n_missing <- sum(missing_idx)
    
    if (n_missing == 0) next
    
    # Strategy A: Use target median if >30% data available
    if (mean(!missing_idx) > 0.30) {
      fill_value <- median(target_vec, na.rm = TRUE)
      imputed_df[missing_idx, marker] <- fill_value
      
    } else {
      # Strategy B: Use reference distribution with jitter
      ref_mean <- mean(ref_vec, na.rm = TRUE)
      ref_sd <- sd(ref_vec, na.rm = TRUE)
      
      # Safety: handle degenerate cases
      if (is.na(ref_sd) || ref_sd == 0) {
        ref_sd <- 0.001
        ref_mean <- ifelse(is.na(ref_mean), 0, ref_mean)
      }
      
      set.seed(123)
      synthetic_values <- rnorm(n_missing, mean = ref_mean, sd = ref_sd)
      synthetic_values[synthetic_values < 0] <- 0.001  # Enforce positivity
      
      imputed_df[missing_idx, marker] <- synthetic_values
    }
  }
  
  return(imputed_df)
}

#' Apply centered log-ratio transformation
#' @param df Data frame of compositional data
#' @param offset Small constant added before transformation
#' @return CLR-transformed data frame
apply_clr_transform <- function(df, offset = 0.01) {
  mat <- as.matrix(df) + offset
  clr_result <- compositions::clr(mat)
  return(as.data.frame(clr_result))
}

#' Export network edges to Cytoscape-compatible format
#' @param matrix_in Correlation/weight matrix
#' @param marker_names Node names
#' @param filename Output filename
#' @param output_dir Output directory
#' @param network_label Network type label
#' @param pvalues Optional p-value matrix
export_network_file <- function(matrix_in, marker_names, filename, output_dir, 
                                network_label, pvalues = NULL) {
  if (is.null(matrix_in) || sum(abs(matrix_in), na.rm = TRUE) == 0) {
    cat(sprintf("Skipped %s: empty network\n", filename))
    return(invisible(NULL))
  }
  
  rownames(matrix_in) <- marker_names
  colnames(matrix_in) <- marker_names
  
  # Convert to edge list (upper triangle only)
  edges <- reshape2::melt(as.matrix(matrix_in), value.name = "weight")
  colnames(edges)[1:2] <- c("source", "target")
  
  edges <- edges %>%
    filter(as.character(source) < as.character(target),  # Remove duplicates
           weight != 0) %>%                               # Remove zero weights
    mutate(
      network_type = network_label,
      interaction = ifelse(weight > 0, "positive", "negative"),
      abs_weight = abs(weight)
    )
  
  # Add p-values if available
  if (!is.null(pvalues)) {
    pval_edges <- reshape2::melt(as.matrix(pvalues), value.name = "p_value")
    colnames(pval_edges)[1:2] <- c("source", "target")
    edges <- left_join(edges, pval_edges, by = c("source", "target"))
  }
  
  if (nrow(edges) > 0) {
    write.table(edges, file.path(output_dir, filename), 
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat(sprintf("Exported %s: %d edges\n", filename, nrow(edges)))
  }
  
  return(edges)
}

# ------------------------------------------------------------------------------
# MAIN ANALYSIS PIPELINE
# ------------------------------------------------------------------------------

main_analysis <- function(config) {
  
  # ============================================================================
  # STEP 1: DATA LOADING
  # ============================================================================
  cat("\n[STEP 1] Loading data...\n")
  
  raw_nsclc <- load_excel_numeric(config$data_path, config$sheets$nsclc)
  raw_hnscc <- load_excel_numeric(config$data_path, config$sheets$hnscc)
  raw_healthy <- load_excel_numeric(config$data_path, config$sheets$healthy)
  
  # Identify common markers across all datasets
  common_markers <- Reduce(intersect, list(
    colnames(raw_nsclc),
    colnames(raw_hnscc),
    colnames(raw_healthy)
  ))
  marker_cols <- setdiff(common_markers, "patient_id")
  
  cat(sprintf("  Common markers identified: %d\n", length(marker_cols)))
  
  # ============================================================================
  # STEP 2: QUALITY CONTROL & MARKER FILTERING
  # ============================================================================
  cat("\n[STEP 2] Marker quality control...\n")
  
  # Combine tumor cohorts for unified filtering
  tumor_combined <- rbind(
    raw_nsclc[, marker_cols, drop = FALSE],
    raw_hnscc[, marker_cols, drop = FALSE]
  )
  
  completeness_rates <- get_completeness(tumor_combined)
  
  # Filter markers by completeness
  valid_markers <- names(completeness_rates)[completeness_rates >= config$min_completeness]
  dropped_markers <- names(completeness_rates)[completeness_rates < config$min_completeness]
  
  # Report dropped markers
  if (length(dropped_markers) > 0) {
    cat("\n  [DROPPED] High missingness (<%d%% complete):\n", 
        config$min_completeness * 100)
    for (marker in dropped_markers) {
      missing_pct <- (1 - completeness_rates[marker]) * 100
      cat(sprintf("    - %-20s: %.1f%% missing\n", marker, missing_pct))
    }
  }
  
  cat(sprintf("\n  Markers passing QC: %d\n", length(valid_markers)))
  
  # ============================================================================
  # STEP 3: IMPUTATION
  # ============================================================================
  cat("\n[STEP 3] Imputing missing values...\n")
  
  nsclc_imputed <- impute_with_reference(raw_nsclc, raw_nsclc, valid_markers)
  hnscc_imputed <- impute_with_reference(raw_hnscc, raw_nsclc, valid_markers)  # NSCLC as reference
  healthy_imputed <- impute_with_reference(raw_healthy, raw_healthy, valid_markers)
  
  # Merge tumor cohorts
  tumor_final <- rbind(
    nsclc_imputed[, valid_markers, drop = FALSE],
    hnscc_imputed[, valid_markers, drop = FALSE]
  )
  healthy_final <- healthy_imputed[, valid_markers, drop = FALSE]
  
  # ============================================================================
  # STEP 4: VARIANCE CHECK
  # ============================================================================
  cat("\n[STEP 4] Checking for zero-variance markers...\n")
  
  tumor_variance <- apply(tumor_final, 2, var)
  healthy_variance <- apply(healthy_final, 2, var)
  
  variance_valid <- (tumor_variance > 0) & (healthy_variance > 0)
  final_markers <- valid_markers[variance_valid]
  dropped_variance <- valid_markers[!variance_valid]
  
  if (length(dropped_variance) > 0) {
    cat("\n  [DROPPED] Zero variance:\n")
    for (marker in dropped_variance) {
      cat(sprintf("    - %-20s: Var(Tumor)=%.4f, Var(Healthy)=%.4f\n", 
                  marker, tumor_variance[marker], healthy_variance[marker]))
    }
  }
  
  # Apply final filter
  tumor_final <- tumor_final[, final_markers]
  healthy_final <- healthy_final[, final_markers]
  
  cat(sprintf("\n  Final marker set: %d markers\n", length(final_markers)))
  cat(sprintf("  Sample sizes: Tumor=%d, Healthy=%d\n", 
              nrow(tumor_final), nrow(healthy_final)))
  
  # ============================================================================
  # STEP 5: CLR TRANSFORMATION
  # ============================================================================
  cat("\n[STEP 5] Applying CLR transformation...\n")
  
  tumor_clr <- apply_clr_transform(tumor_final, config$clr_offset)
  healthy_clr <- apply_clr_transform(healthy_final, config$clr_offset)
  
  # Ensure column alignment
  tumor_clr <- tumor_clr[, final_markers]
  healthy_clr <- healthy_clr[, final_markers]
  
  # ============================================================================
  # STEP 6: NETWORK COMPARISON TEST
  # ============================================================================
  cat(sprintf("\n[STEP 6] Running Network Comparison Test (%d permutations)...\n", 
              config$nct_iterations))
  cat("  This may take 1-2 minutes...\n")
  
  set.seed(999)
  nct_result <- NCT(
    tumor_clr, 
    healthy_clr,
    it = config$nct_iterations,
    gamma = config$nct_gamma,
    test.edges = TRUE,
    test.centrality = TRUE,
    progressbar = TRUE,
    paired = FALSE,
    weighted = TRUE,
    AND = TRUE,
    abs = FALSE
  )
  
  cat("  NCT complete.\n")
  
  # ============================================================================
  # STEP 7: EXTRACT DIFFERENTIAL EDGES
  # ============================================================================
  cat("\n[STEP 7] Extracting differential edges...\n")
  
  # Compute correlation matrices
  cor_tumor <- cor(tumor_clr, method = config$correlation_method)
  cor_healthy <- cor(healthy_clr, method = config$correlation_method)
  cor_diff <- cor_tumor - cor_healthy
  
  # Extract upper triangle
  upper_idx <- upper.tri(cor_tumor)
  nodes <- colnames(cor_tumor)
  
  # Build edge data frame
  edge_df <- data.frame(
    node1 = nodes[row(cor_tumor)[upper_idx]],
    node2 = nodes[col(cor_tumor)[upper_idx]],
    r_tumor = cor_tumor[upper_idx],
    r_healthy = cor_healthy[upper_idx],
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      edge_id = paste(pmin(node1, node2), pmax(node1, node2), sep = "_"),
      delta_r = r_tumor - r_healthy
    )
  
  # Extract p-values from NCT
  if (!is.null(nct_result$einv.pvals)) {
    nct_pvals <- as.data.frame(nct_result$einv.pvals)
    
    if (ncol(nct_pvals) >= 3) {
      colnames(nct_pvals)[1:3] <- c("n1", "n2", "p_value")
      nct_pvals <- nct_pvals %>%
        mutate(
          edge_id = paste(pmin(as.character(n1), as.character(n2)),
                          pmax(as.character(n1), as.character(n2)), sep = "_"),
          p_value = as.numeric(as.character(p_value))
        ) %>%
        select(edge_id, p_value)
      
      # Merge with edge data
      edge_df <- left_join(edge_df, nct_pvals, by = "edge_id")
    }
  }
  
  # Apply FDR correction
  if ("p_value" %in% colnames(edge_df)) {
    edge_df <- edge_df %>%
      mutate(
        p_adj = p.adjust(p_value, method = config$fdr_method),
        significant = !is.na(p_adj) & p_adj < config$alpha_level
      )
  } else {
    edge_df$p_value <- NA
    edge_df$p_adj <- NA
    edge_df$significant <- FALSE
  }
  
  # Classify edges
  edge_df <- edge_df %>%
    mutate(
      category = case_when(
        delta_r > config$min_delta_r ~ "tumor_gain",
        delta_r < -config$min_delta_r ~ "healthy_gain",
        TRUE ~ "stable"
      )
    ) %>%
    arrange(p_adj, desc(abs(delta_r)))
  
  cat(sprintf("  Total edges analyzed: %d\n", nrow(edge_df)))
  cat(sprintf("  Significant edges (FDR < %.2f): %d\n", 
              config$alpha_level, sum(edge_df$significant, na.rm = TRUE)))
  
  # ============================================================================
  # STEP 8: EXPORT RESULTS
  # ============================================================================
  cat("\n[STEP 8] Exporting results...\n")
  
  # Significant edges
  sig_edges <- edge_df %>%
    filter(significant | abs(delta_r) > 0.4) %>%
    select(node1, node2, r_tumor, r_healthy, delta_r, p_value, p_adj, 
           significant, category)
  
  write.csv(sig_edges, 
            file.path(config$output_dir, "significant_differential_edges.csv"),
            row.names = FALSE)
  cat(sprintf("  Saved significant edges: %d\n", nrow(sig_edges)))
  
  # Cytoscape exports
  cat("\n  Exporting Cytoscape networks...\n")
  cytoscape_dir <- file.path(config$output_dir, "cytoscape")
  if (!dir.exists(cytoscape_dir)) dir.create(cytoscape_dir)
  
  export_network_file(nct_result$nw1, final_markers, 
                      "network_healthy.txt", cytoscape_dir, "healthy")
  
  export_network_file(nct_result$nw2, final_markers,
                      "network_tumor.txt", cytoscape_dir, "tumor")
  
  # Differential network (sorted by p-value and delta)
  diff_edges <- edge_df %>%
    filter(significant | abs(delta_r) > config$min_delta_r) %>%
    arrange(p_adj, desc(abs(delta_r))) %>%
    select(source = node1, target = node2, delta_r, p_value, p_adj, category)
  
  write.table(diff_edges, 
              file.path(cytoscape_dir, "network_differential.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  Saved differential network: %d edges\n", nrow(diff_edges)))
  
  # ============================================================================
  # STEP 9: VISUALIZATIONS
  # ============================================================================
  cat("\n[STEP 9] Generating visualizations...\n")
  
  pdf(file.path(config$output_dir, "network_plots.pdf"), width = 12, height = 8)
  
  # Compute shared layout
  layout_shared <- qgraph(cor_healthy, DoNotPlot = TRUE)$layout
  
  # Plot 1: Difference network
  qgraph(cor_diff, 
         layout = layout_shared,
         graph = "cor",
         title = "Differential Network (Tumor - Healthy)\nRed = Tumor-specific | Blue = Healthy-specific",
         cut = 0.2,
         minimum = 0.2,
         edge.labels = FALSE,
         posCol = "#D55E00",  # Colorblind-safe red
         negCol = "#0072B2",  # Colorblind-safe blue
         vsize = 6,
         borders = FALSE,
         labels = final_markers)
  
  # Plot 2 & 3: Side-by-side comparison
  par(mfrow = c(1, 2))
  qgraph(cor_tumor, 
         layout = layout_shared,
         title = sprintf("Tumor Network (N=%d)", nrow(tumor_clr)),
         graph = "cor",
         minimum = 0.35,
         posCol = "darkred",
         negCol = "darkblue",
         vsize = 5)
  
  qgraph(cor_healthy,
         layout = layout_shared,
         title = sprintf("Healthy Network (N=%d)", nrow(healthy_clr)),
         graph = "cor",
         minimum = 0.35,
         posCol = "darkred",
         negCol = "darkblue",
         vsize = 5)
  
  dev.off()
  
  # Scatter plot
  p_scatter <- ggplot(edge_df, aes(x = r_healthy, y = r_tumor)) +
    geom_point(aes(color = significant), alpha = 0.5, size = 2) +
    scale_color_manual(values = c("grey70", "red3")) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    theme_minimal(base_size = 12) +
    labs(
      title = "Correlation Stability Analysis",
      subtitle = "Points far from diagonal indicate altered immune co-expression",
      x = "Correlation in Healthy",
      y = "Correlation in Tumor",
      color = "Significant"
    )
  
  ggsave(file.path(config$output_dir, "scatter_correlation_shift.pdf"), 
         p_scatter, width = 7, height = 6)
  
  cat("  Visualizations saved.\n")
  
  # ============================================================================
  # STEP 10: SUMMARY REPORT
  # ============================================================================
  cat("\n[STEP 10] Writing summary report...\n")
  
  sink(file.path(config$output_dir, "analysis_report.txt"))
  cat("=================================================================\n")
  cat("        PAN-CANCER DIFFERENTIAL NETWORK ANALYSIS REPORT          \n")
  cat("=================================================================\n\n")
  
  cat("DATASET COMPOSITION\n")
  cat("-------------------\n")
  cat(sprintf("Tumor samples (NSCLC + HNSCC): %d\n", nrow(tumor_clr)))
  cat(sprintf("Healthy samples: %d\n", nrow(healthy_clr)))
  cat(sprintf("Markers analyzed: %d\n\n", length(final_markers)))
  
  cat("GLOBAL NETWORK STATISTICS\n")
  cat("-------------------------\n")
  cat(sprintf("Network structure invariance p-value: %.4f\n", nct_result$nwpval))
  cat(sprintf("Global strength difference p-value: %.4f\n", nct_result$glstrinv.pval))
  cat("(p < 0.05 indicates significant network reorganization)\n\n")
  
  cat("DIFFERENTIAL EDGE SUMMARY\n")
  cat("-------------------------\n")
  cat(sprintf("Total edges tested: %d\n", nrow(edge_df)))
  cat(sprintf("Significant edges (FDR < %.2f): %d\n", 
              config$alpha_level, sum(edge_df$significant, na.rm = TRUE)))
  cat(sprintf("High-magnitude shifts (|ΔR| > %.2f): %d\n\n",
              0.4, sum(abs(edge_df$delta_r) > 0.4, na.rm = TRUE)))
  
  cat("TOP 15 DIFFERENTIAL EDGES (by FDR-adjusted p-value)\n")
  cat("----------------------------------------------------\n")
  top_edges <- head(sig_edges, 15)
  print(top_edges[, c("node1", "node2", "delta_r", "p_adj", "category")], 
        row.names = FALSE)
  
  sink()
  
  cat("\n=== ANALYSIS COMPLETE ===\n")
  cat("All results saved to:", config$output_dir, "\n")
  
  return(list(
    edge_data = edge_df,
    nct_result = nct_result,
    markers = final_markers,
    config = config
  ))
}

# ==============================================================================
# EXECUTE PIPELINE
# ==============================================================================
results <- main_analysis(config)