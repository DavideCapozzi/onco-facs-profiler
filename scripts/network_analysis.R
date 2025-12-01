# ==============================================================================
# PROJECT: Robust Pan-Cancer Differential Network Analysis
# DATE: 2024-05-22
# CONTEXT: High-dimensional FACS data (Small N, Large P)
# METHOD: Spearman Correlation + Permutation Test + Smart Logging
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. LIBRARIES
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(qgraph)
  library(zCompositions)  # Best practice for zero handling
  library(VIM)            # For kNN imputation
  library(reshape2)
  library(foreach)        # For parallel processing
  library(doParallel)
  library(openxlsx)       # For Excel export
})

# ------------------------------------------------------------------------------
# 2. CONFIGURATION
# ------------------------------------------------------------------------------
config <- list(
  # Paths
  data_path = "/mnt/c/Users/Davide/Desktop/Davide/bandi/0dottorato/projects/long survivors/DB_anonimo_standardizzatov2.xlsx", 
  output_dir = "/home/davidec/projects/compositional_analysis/",
  
  # Analysis Parameters
  imputation_mode = "knn",      # "knn" (Smart) or "filter" (Strict)
  n_permutations = 1000,        # Number of permutations for p-values
  correlation_method = "spearman", 
  
  # Thresholds
  min_completeness_row = 0.50,  # Drop patient if >50% markers are missing
  min_completeness_col = 0.60,  # Drop marker if >40% patients are missing (completeness < 0.6)
  fdr_threshold = 0.05,
  min_delta_rho = 0.30,         # Minimum correlation difference to be relevant
  
  # Parallel Processing
  n_cores = parallel::detectCores() - 1
)

# Output directory setup
config$output_dir <- paste0(config$output_dir,"Results_", format(Sys.time(), "%Y%m%d_%H%M"))
if (!dir.exists(config$output_dir)) dir.create(config$output_dir)

# Initialize a log list to store report details
report_log <- list()

cat("=== PIPELINE START ===\n")
cat("Output Directory:", config$output_dir, "\n\n")

# ------------------------------------------------------------------------------
# 3. HELPER FUNCTIONS
# ------------------------------------------------------------------------------

#' Robust Data Loader
load_cohort <- function(path, sheet_name) {
  tryCatch({
    df <- read_excel(path, sheet = sheet_name)
    if (!"Patient_ID" %in% colnames(df)) stop("Patient_ID column missing!")
    
    ids <- df$Patient_ID
    markers <- df %>% select(-Patient_ID) %>%
      mutate(across(everything(), ~suppressWarnings(as.numeric(.)))) %>%
      as.data.frame()
    
    rownames(markers) <- ids
    cat(sprintf("  Loaded %s: %d samples, %d markers\n", sheet_name, nrow(markers), ncol(markers)))
    return(markers)
  }, error = function(e) {
    stop(paste("Error loading", sheet_name, ":", e$message))
  })
}

#' Verbose Filtering
#' Logs specific IDs and percentages of missingness for dropped items
clean_data_verbose <- function(df, max_na_row, max_na_col, group_name) {
  
  # 1. Filter Patients (Rows)
  row_missing_rate <- rowMeans(is.na(df))
  bad_rows_idx <- which(row_missing_rate > max_na_row)
  
  if (length(bad_rows_idx) > 0) {
    # Format: "ID (55%)"
    dropped_info <- paste0(
      rownames(df)[bad_rows_idx], 
      " (", round(row_missing_rate[bad_rows_idx] * 100, 1), "%)"
    )
    cat(sprintf("    [DROP] %s Patients (>%.0f%% missing): %s\n", 
                group_name, max_na_row*100, paste(dropped_info, collapse = ", ")))
    
    # Store for report
    report_log[[paste0("dropped_patients_", group_name)]] <<- dropped_info
    df <- df[-bad_rows_idx, ]
  }
  
  # 2. Filter Markers (Columns)
  # Max allowed missing is (1 - min_completeness)
  col_missing_rate <- colMeans(is.na(df))
  bad_cols_idx <- which(col_missing_rate > (1 - max_na_col))
  
  if (length(bad_cols_idx) > 0) {
    dropped_cols <- paste0(
      names(bad_cols_idx), 
      " (", round(col_missing_rate[bad_cols_idx] * 100, 1), "%)"
    )
    cat(sprintf("    [DROP] %s Markers (>%.0f%% missing): %s\n", 
                group_name, (1-max_na_col)*100, paste(dropped_cols, collapse = ", ")))
    
    report_log[[paste0("dropped_markers_", group_name)]] <<- dropped_cols
    df <- df[, -bad_cols_idx]
  }
  
  return(df)
}

#' Smart Imputation (kNN)
impute_knn_optimized <- function(df, k=3) {
  cat("    [IMPUTE] Running k-Nearest Neighbors (k=3)...\n")
  imputed_res <- VIM::kNN(df, k = k, imp_var = FALSE, trace = FALSE)
  return(imputed_res)
}

#' Compositional Transformation
transform_compositional <- function(df) {
  n_zeros <- sum(df == 0, na.rm = TRUE)
  if (n_zeros > 0) {
    cat(sprintf("    [TRANSFORM] Replacing %d zeros with CZM method...\n", n_zeros))
    df_clean <- zCompositions::cmultRepl(df, method = "CZM", output = "p-counts", suppress.print = TRUE)
  } else {
    df_clean <- df
  }
  cat("    [TRANSFORM] Applying CLR transformation...\n")
  clr_data <- compositions::clr(df_clean)
  return(as.data.frame(clr_data))
}

#' Permutation Test
run_permutation_test <- function(data_tumor, data_healthy, n_perms, n_cores) {
  
  cor_t <- cor(data_tumor, method = "spearman")
  cor_h <- cor(data_healthy, method = "spearman")
  diff_obs <- cor_t - cor_h
  
  combined <- rbind(data_tumor, data_healthy)
  n_t <- nrow(data_tumor)
  n_tot <- nrow(combined)
  
  cat(sprintf("    [TEST] Running %d permutations on %d cores...\n", n_perms, n_cores))
  
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  perm_results <- foreach(i = 1:n_perms, .combine = 'c', .packages = 'stats') %dopar% {
    perm_idx <- sample(1:n_tot)
    dt_perm <- combined[perm_idx[1:n_t], ]
    dh_perm <- combined[perm_idx[(n_t + 1):n_tot], ]
    c_t <- cor(dt_perm, method = "spearman")
    c_h <- cor(dh_perm, method = "spearman")
    list(c_t - c_h)
  }
  stopCluster(cl)
  
  cat("    [TEST] Calculating p-values...\n")
  p_mat <- matrix(0, nrow=nrow(diff_obs), ncol=ncol(diff_obs))
  for(r in 1:nrow(diff_obs)) {
    for(c in 1:ncol(diff_obs)) {
      obs_val <- abs(diff_obs[r,c])
      greater_counts <- sum(sapply(perm_results, function(x) abs(x[r,c]) >= obs_val))
      p_mat[r,c] <- (greater_counts + 1) / (n_perms + 1)
    }
  }
  
  rownames(p_mat) <- rownames(diff_obs)
  colnames(p_mat) <- colnames(diff_obs)
  
  return(list(diff_matrix = diff_obs, p_matrix = p_mat, cor_tumor = cor_t, cor_healthy = cor_h))
}

# ------------------------------------------------------------------------------
# 4. MAIN PIPELINE
# ------------------------------------------------------------------------------

# --- Step 1: Loading ---
cat("\n[STEP 1] Loading Data...\n")
raw_nsclc <- load_cohort(config$data_path, "NSCLC")
raw_hnscc <- load_cohort(config$data_path, "HNSCC")
raw_healthy <- load_cohort(config$data_path, "Healthy_Donors")

common_markers <- Reduce(intersect, list(colnames(raw_nsclc), colnames(raw_hnscc), colnames(raw_healthy)))
tumor_combined <- rbind(raw_nsclc[, common_markers], raw_hnscc[, common_markers])
healthy_clean <- raw_healthy[, common_markers]

report_log$n_tumor_raw <- nrow(tumor_combined)
report_log$n_healthy_raw <- nrow(healthy_clean)

# --- Step 2: Filtering & Imputation ---
cat("\n[STEP 2] Pre-processing...\n")

process_group <- function(df, name) {
  cat(sprintf("\n  Processing %s group:\n", name))
  df_filt <- clean_data_verbose(df, config$min_completeness_row, config$min_completeness_col, name)
  
  if (config$imputation_mode == "knn") {
    if(any(is.na(df_filt))) df_final <- impute_knn_optimized(df_filt, k = 3)
    else df_final <- df_filt
  } else {
    n_before <- nrow(df_filt)
    df_final <- na.omit(df_filt)
    if(nrow(df_final) < n_before) cat(sprintf("    [FILTER] Dropped %d incomplete cases\n", n_before - nrow(df_final)))
  }
  return(df_final)
}

tumor_proc <- process_group(tumor_combined, "Tumor")
healthy_proc <- process_group(healthy_clean, "Healthy")

report_log$n_tumor_final <- nrow(tumor_proc)
report_log$n_healthy_final <- nrow(healthy_proc)

final_markers <- intersect(colnames(tumor_proc), colnames(healthy_proc))
tumor_final <- tumor_proc[, final_markers]
healthy_final <- healthy_proc[, final_markers]
report_log$n_markers_final <- length(final_markers)
report_log$markers_list <- paste(final_markers, collapse=", ")

# --- Step 3: Transformation ---
cat("\n[STEP 3] Transformation...\n")
tumor_clr <- transform_compositional(tumor_final)
healthy_clr <- transform_compositional(healthy_final)

# --- Step 4: Analysis ---
cat("\n[STEP 4] Differential Correlation...\n")
results <- run_permutation_test(tumor_clr, healthy_clr, config$n_permutations, config$n_cores)

# FDR Correction
p_vec <- results$p_matrix[upper.tri(results$p_matrix)]
p_adj_vec <- p.adjust(p_vec, method = "fdr")
p_adj_mat <- matrix(NA, nrow=nrow(results$p_matrix), ncol=ncol(results$p_matrix))
p_adj_mat[upper.tri(p_adj_mat)] <- p_adj_vec
# Mirror for completeness if needed, but usually we just need upper tri
rownames(p_adj_mat) <- rownames(results$p_matrix)
colnames(p_adj_mat) <- colnames(results$p_matrix)

# --- Step 5: Export (Excel) ---
cat("\n[STEP 5] Exporting Results (XLSX)...\n")

edges <- reshape2::melt(results$diff_matrix)
colnames(edges) <- c("Node1", "Node2", "Delta_Rho")
edges$P_Value <- reshape2::melt(results$p_matrix)$value
edges$P_FDR <- reshape2::melt(p_adj_mat)$value

edges_clean <- edges %>% 
  filter(as.character(Node1) < as.character(Node2)) %>%
  mutate(
    Category = case_when(
      P_FDR < config$fdr_threshold & Delta_Rho > config$min_delta_rho ~ "Tumor_Gain",
      P_FDR < config$fdr_threshold & Delta_Rho < -config$min_delta_rho ~ "Healthy_Gain",
      TRUE ~ "Stable"
    ),
    Significant = Category != "Stable"
  ) %>%
  arrange(P_FDR, desc(abs(Delta_Rho)))

sig_edges <- edges_clean %>% filter(Significant)

# Create Excel Workbook
wb <- createWorkbook()

addWorksheet(wb, "Significant_Edges")
writeData(wb, "Significant_Edges", sig_edges)

addWorksheet(wb, "All_Edges")
writeData(wb, "All_Edges", edges_clean)

addWorksheet(wb, "Config_Info")
config_df <- data.frame(Parameter = names(config), Value = as.character(config))
writeData(wb, "Config_Info", config_df)

saveWorkbook(wb, file.path(config$output_dir, "Analysis_Results.xlsx"), overwrite = TRUE)
cat("  Saved: Analysis_Results.xlsx\n")

# --- Step 6: Visualization ---
cat("\n[STEP 6] Visualization...\n")

pdf(file.path(config$output_dir, "network_plots_corrected.pdf"), width = 14, height = 8)

# Layout: Calculate on Healthy (most stable)
L <- qgraph(results$cor_healthy, DoNotPlot=TRUE)$layout

# Plot 1: Differential Network (Corrected)
# We use graph="default" to avoid positive-definite checks, and diag=FALSE to ignore 0s
qgraph(results$diff_matrix,
       layout = L,
       graph = "default",      # FIX for "matrix not positive definite"
       diag = FALSE,           # FIX for "non-finite weights"
       title = "Differential Network (Tumor - Healthy)\nRed = Tumor | Blue = Healthy",
       cut = 0.2,              # Only show edges with |diff| > 0.2
       minimum = 0.2,
       posCol = "#D55E00", 
       negCol = "#0072B2",
       vsize = 4,
       labels = colnames(results$diff_matrix),
       label.scale = FALSE,
       label.cex = 0.8,
       edge.labels = FALSE)    # Too cluttered with labels

# Plot 2 & 3: Side by Side
par(mfrow=c(1,2))
qgraph(results$cor_tumor, layout=L, title="Tumor Network (Spearman)", 
       graph="cor", minimum=0.3, vsize=4, labels=colnames(results$diff_matrix),
       posCol = "darkred", negCol = "darkblue")
qgraph(results$cor_healthy, layout=L, title="Healthy Network (Spearman)", 
       graph="cor", minimum=0.3, vsize=4, labels=colnames(results$diff_matrix),
       posCol = "darkred", negCol = "darkblue")

dev.off()
cat("  Saved: network_plots_corrected.pdf\n")

# --- Step 7: Summary Report ---
cat("\n[STEP 7] Generating Analysis Report...\n")

report_file <- file.path(config$output_dir, "Analysis_Report.txt")
sink(report_file)

cat("===============================================================\n")
cat("       PAN-CANCER IMMUNE NETWORK ANALYSIS REPORT               \n")
cat("===============================================================\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

cat("1. PARAMETERS\n")
cat("-------------\n")
cat(sprintf("Imputation: %s\n", config$imputation_mode))
cat(sprintf("Permutations: %d\n", config$n_permutations))
cat(sprintf("Correlation: %s (CLR Transformed)\n", config$correlation_method))
cat(sprintf("FDR Threshold: %.2f\n", config$fdr_threshold))
cat(sprintf("Min Delta Rho: %.2f\n\n", config$min_delta_rho))

cat("2. DATA CLEANING SUMMARY\n")
cat("------------------------\n")
cat(sprintf("Tumor Samples: %d (Raw) -> %d (Final)\n", report_log$n_tumor_raw, report_log$n_tumor_final))
if(!is.null(report_log$dropped_patients_Tumor)) 
  cat("  DROPPED TUMOR PATIENTS:", report_log$dropped_patients_Tumor, "\n")

cat(sprintf("Healthy Samples: %d (Raw) -> %d (Final)\n", report_log$n_healthy_raw, report_log$n_healthy_final))
if(!is.null(report_log$dropped_patients_Healthy)) 
  cat("  DROPPED HEALTHY PATIENTS:", report_log$dropped_patients_Healthy, "\n")

cat(sprintf("\nMarkers Analyzed: %d\n", report_log$n_markers_final))
if(!is.null(report_log$dropped_markers_Tumor)) 
  cat("  DROPPED TUMOR MARKERS:", report_log$dropped_markers_Tumor, "\n")
if(!is.null(report_log$dropped_markers_Healthy)) 
  cat("  DROPPED HEALTHY MARKERS:", report_log$dropped_markers_Healthy, "\n")

cat("\n3. RESULTS SUMMARY\n")
cat("------------------\n")
cat(sprintf("Total Edges Tested: %d\n", nrow(edges_clean)))
cat(sprintf("Significant Edges (FDR < %.2f): %d\n", config$fdr_threshold, nrow(sig_edges)))
cat(sprintf("  - Gain in Tumor: %d\n", sum(sig_edges$Category == "Tumor_Gain")))
cat(sprintf("  - Gain in Healthy: %d\n", sum(sig_edges$Category == "Healthy_Gain")))

cat("\n4. TOP 10 DIFFERENTIAL INTERACTIONS\n")
cat("-----------------------------------\n")
if(nrow(sig_edges) > 0) {
  print(head(sig_edges[, c("Node1", "Node2", "Delta_Rho", "P_FDR", "Category")], 10), row.names=FALSE)
} else {
  cat("No significant edges found.\n")
}

cat("\n5. INTERPRETATION NOTES\n")
cat("-----------------------\n")
cat("- 'Delta_Rho' is Correlation(Tumor) - Correlation(Healthy).\n")
cat("- Positive Delta (Red edges in plot) means the interaction is stronger (or more positive) in Tumor.\n")
cat("- Negative Delta (Blue edges in plot) means the interaction is stronger in Healthy.\n")
cat("- Results are based on Spearman correlation on CLR-transformed data.\n")

sink()

cat("  Saved: Analysis_Report.txt\n")
cat("\n=== PIPELINE FINISHED SUCCESSFULLY ===\n")
