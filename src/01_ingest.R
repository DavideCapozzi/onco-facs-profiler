# src/01_ingest.R
# ==============================================================================
# STEP 01: DATA INGESTION & CODA PREPROCESSING
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
})

source("R/infrastructure.R")
source("R/modules_coda.R")

message("\n=== PIPELINE STEP 1: INGESTION & CODA ===")

# 1. Load Config & Data
config <- load_config("config/global_params.yml")
raw_data <- load_raw_data(config)

# 2. Setup Matrices
meta_cols <- c("Patient_ID", "Group")
marker_cols <- setdiff(names(raw_data), meta_cols)

mat_raw <- as.matrix(raw_data[, marker_cols])
rownames(mat_raw) <- raw_data$Patient_ID

message(sprintf("[Data] Initial Matrix: %d Samples x %d Markers", nrow(mat_raw), ncol(mat_raw)))

# 3. Quality Control (Filtering)
# ------------------------------------------------------------------------------
message("[QC] Running Quality Control...")

# Initialize QC Summary Object
qc_summary <- list(
  n_row_init = nrow(mat_raw),
  n_col_init = ncol(mat_raw),
  n_col_zerovar = 0,
  dropped_rows_detail = data.frame(Patient_ID=character(), NA_Percent=numeric()),
  dropped_cols_detail = data.frame(Marker=character(), NA_Percent=numeric())
)

# A. Remove Constant Columns (Zero Variance)
col_vars <- apply(mat_raw, 2, var, na.rm=TRUE)
const_cols <- which(col_vars == 0 | is.na(col_vars))

if (length(const_cols) > 0) {
  message(sprintf("   [QC] Dropping %d markers with zero variance.", length(const_cols)))
  qc_summary$n_col_zerovar <- length(const_cols)
  mat_raw <- mat_raw[, -const_cols]
}

# B. Filter by Missingness (Rows - Patients)
row_na_freq <- rowMeans(is.na(mat_raw))
drop_row_idx <- which(row_na_freq > config$qc$max_na_row_pct)

if (length(drop_row_idx) > 0) {
  # 1. Log to Console
  message(sprintf("   [QC] Dropping %d patients with >%.0f%% missingness:", 
                  length(drop_row_idx), config$qc$max_na_row_pct * 100))
  
  dropped_pats <- data.frame(
    Patient_ID = rownames(mat_raw)[drop_row_idx],
    NA_Percent = round(row_na_freq[drop_row_idx] * 100, 2)
  )
  
  # Print nicely formatted list
  print_df <- paste0("      - ", dropped_pats$Patient_ID, " (", dropped_pats$NA_Percent, "%)")
  cat(paste(print_df, collapse = "\n"), "\n")
  
  # 2. Store for Report
  qc_summary$dropped_rows_detail <- dropped_pats
  qc_summary$n_row_dropped <- length(drop_row_idx)
  
  # 3. Apply Filter
  mat_raw  <- mat_raw[-drop_row_idx, ]
  raw_data <- raw_data[-drop_row_idx, ] # Sync Metadata
} else {
  qc_summary$n_row_dropped <- 0
}

# C. Filter by Missingness (Cols - Markers)
col_na_freq <- colMeans(is.na(mat_raw))
drop_col_idx <- which(col_na_freq > config$qc$max_na_col_pct)

if (length(drop_col_idx) > 0) {
  # 1. Log to Console
  message(sprintf("   [QC] Dropping %d markers with >%.0f%% missingness:", 
                  length(drop_col_idx), config$qc$max_na_col_pct * 100))
  
  dropped_markers <- data.frame(
    Marker = colnames(mat_raw)[drop_col_idx],
    NA_Percent = round(col_na_freq[drop_col_idx] * 100, 2)
  )
  
  # Print nicely formatted list
  print_df <- paste0("      - ", dropped_markers$Marker, " (", dropped_markers$NA_Percent, "%)")
  cat(paste(print_df, collapse = "\n"), "\n")
  
  # 2. Store for Report
  qc_summary$dropped_cols_detail <- dropped_markers
  qc_summary$n_col_dropped <- length(drop_col_idx)
  
  # 3. Apply Filter
  mat_raw <- mat_raw[, -drop_col_idx]
} else {
  qc_summary$n_col_dropped <- 0
}

# Final Stats
marker_cols <- colnames(mat_raw)
raw_data <- raw_data[, c(meta_cols, marker_cols)]

qc_summary$n_row_final <- nrow(mat_raw)
qc_summary$n_col_final <- ncol(mat_raw)

message(sprintf("   [QC] Final Dimensions: %d Samples x %d Markers", nrow(mat_raw), ncol(mat_raw)))

# --- NEW: SAVE EXCEL REPORT ---
out_dir <- file.path(config$output_root, "01_processed")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

save_qc_report(qc_summary, file.path(out_dir, "QC_Filtering_Report.xlsx"))

if (nrow(mat_raw) < 10) warning("[QC] WARNING: Very low sample size. Analysis may be unstable.")


# 4. CoDa Pipeline Execution
# ------------------------------------------------------------------------------

# A. Zero Handling (Mask -> Fill -> Fix -> Restore)
mat_no_zeros <- replace_zeros_czm(mat_raw)

# B. NA Imputation (Optimization + Imputation)
final_k <- config$imputation$knn_k_default 

if (config$imputation$knn_optimize) {
  # Range 3 to 10 (or less if few samples)
  k_max <- min(10, floor(nrow(mat_raw) * 0.5)) # Safety cap for small N
  k_range <- 3:k_max
  
  if(length(k_range) > 1) {
    opt_results <- optimize_knn_k(mat_no_zeros, k_seq = k_range)
    final_k <- opt_results$best_k
    
    # Save Plot
    qc_dir <- file.path(config$output_root, "01_QC")
    if (!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE)
    ggsave(file.path(qc_dir, "imputation_optimization_k.pdf"), opt_results$plot, width=6, height=4)
  }
}

message(sprintf("[CoDa] Applying Imputation with k=%d...", final_k))
mat_imputed <- impute_knn_aitchison(mat_no_zeros, k = final_k)

# C. CLR Transformation
mat_clr <- transform_clr(mat_imputed)

# 5. Save Output
# ------------------------------------------------------------------------------
# Re-attach metadata
df_imputed <- cbind(raw_data[, meta_cols], as.data.frame(mat_imputed))
df_clr     <- cbind(raw_data[, meta_cols], as.data.frame(mat_clr))

# Create Master Object
processed_data <- list(
  metadata     = raw_data[, meta_cols],
  markers      = marker_cols,       # UPDATED LIST
  raw_matrix   = mat_raw,
  imputed_data = df_imputed,
  clr_data     = df_clr,
  config       = config
)

out_dir <- file.path(config$output_root, "01_QC")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

saveRDS(processed_data, file.path(out_dir, "data_processed.rds"))
message("=== STEP 1 COMPLETE ===\n")