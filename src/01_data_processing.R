# src/01_data_processing.R
# ==============================================================================
# STEP 01: DATA INGESTION + QC + HYBRID CODA PREPROCESSING 
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(pcaMethods)
})

source("R/utils_io.R")     
source("R/modules_coda.R")
source("R/modules_qc.R")   

message("\n=== PIPELINE STEP 1: INGESTION + QC + HYBRID TRANSFORM (BPCA) ===")

# 1. Load Config & Data
# ------------------------------------------------------------------------------
config <- load_config("config/global_params.yml")

if (is.null(config$hybrid_groups)) {
  stop("[FATAL] 'hybrid_groups' not found in config. Please update global_params.yml.")
}

raw_data <- load_raw_data(config)

# 2. Setup Matrices
# ------------------------------------------------------------------------------
# Pre-Split Metadata and Matrix
subgroup_col <- config$metadata$subgroup_col
meta_cols <- c("Patient_ID", "Group")
if (!is.null(subgroup_col)) meta_cols <- c(meta_cols, subgroup_col)

# Extract Metadata and Matrix aligned
mat_raw <- as.matrix(raw_data[, setdiff(names(raw_data), meta_cols)])
rownames(mat_raw) <- raw_data$Patient_ID
metadata_raw <- raw_data[, meta_cols]

message(sprintf("[Data] Initial Matrix: %d Samples x %d Markers", nrow(mat_raw), ncol(mat_raw)))

# --- Marker Filtering ---
target_markers <- colnames(mat_raw)
if (!is.null(config$marker_selection$whitelist) && length(config$marker_selection$whitelist) > 0) {
  target_markers <- intersect(target_markers, config$marker_selection$whitelist)
}
if (!is.null(config$marker_selection$blacklist) && length(config$marker_selection$blacklist) > 0) {
  target_markers <- setdiff(target_markers, config$marker_selection$blacklist)
}

mat_raw <- mat_raw[, target_markers, drop = FALSE]

# 3. Quality Control
# ------------------------------------------------------------------------------
# We pass metadata specifically for group-based outlier detection
qc_result <- run_qc_pipeline(mat_raw, metadata_raw, config$qc, dropped_apriori = data.frame())

mat_raw <- qc_result$data
raw_data <- cbind(qc_result$metadata, as.data.frame(mat_raw)) 

qc_summary <- qc_result$report
valid_patients <- rownames(mat_raw)
raw_data <- raw_data %>% filter(Patient_ID %in% valid_patients)

# Save QC Report
out_dir <- file.path(config$output_root, "01_data_processing")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
save_qc_report(qc_summary, file.path(out_dir, "QC_Filtering_Report.xlsx"))

# 4. Hybrid Logic: Split -> Transform -> Impute
# ------------------------------------------------------------------------------
# The transformation logic is now encapsulated in modules_coda.R
transform_results <- perform_hybrid_transformation(mat_raw, config)

# Extract results
mat_hybrid_raw <- transform_results$hybrid_data_raw
mat_hybrid_z   <- transform_results$hybrid_data_z
ilr_list       <- transform_results$ilr_balances
final_markers  <- transform_results$hybrid_markers

# 5. Save Output
# ------------------------------------------------------------------------------
# Bind metadata with transformed data
df_hybrid_raw <- cbind(raw_data[, meta_cols], as.data.frame(mat_hybrid_raw))
df_hybrid_z   <- cbind(raw_data[, meta_cols], as.data.frame(mat_hybrid_z))

processed_data <- list(
  metadata        = raw_data[, meta_cols],
  markers         = colnames(mat_raw),      
  raw_matrix      = mat_raw,
  # Primary outputs for downstream:
  hybrid_markers  = final_markers,
  hybrid_data_raw = df_hybrid_raw,  # Transformed (CLR/Logit), Imputed, NO Z-score
  hybrid_data_z   = df_hybrid_z,    # Transformed, Imputed, AND Z-scored
  ilr_balances    = ilr_list,       # For Permutation tests
  config          = config
)

saveRDS(processed_data, file.path(out_dir, "data_processed.rds"))
message("=== STEP 1 COMPLETE ===\n")