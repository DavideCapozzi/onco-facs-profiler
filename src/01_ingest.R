# src/01_ingest.R
# ==============================================================================
# STEP 01: DATA INGESTION & CODA PREPROCESSING
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
})

# Updated Source Imports
source("R/utils_io.R")     
source("R/modules_coda.R")
source("R/modules_qc.R")   

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

# --- Marker Whitelist/Blacklist ---
initial_markers <- colnames(mat_raw)
target_markers  <- initial_markers

if (!is.null(config$marker_selection$whitelist) && length(config$marker_selection$whitelist) > 0) {
  target_markers <- intersect(target_markers, config$marker_selection$whitelist)
  message(sprintf("   [Filter] Applied Whitelist: restricted to %d markers.", length(target_markers)))
}

if (!is.null(config$marker_selection$blacklist) && length(config$marker_selection$blacklist) > 0) {
  target_markers <- setdiff(target_markers, config$marker_selection$blacklist)
  message(sprintf("   [Filter] Applied Blacklist: removed %d markers.", length(intersect(initial_markers, config$marker_selection$blacklist))))
}

# Identify Excluded Markers (For Report)
dropped_apriori <- setdiff(initial_markers, target_markers)
df_dropped_apriori <- data.frame(
  Marker = dropped_apriori,
  Reason = rep("Excluded by Config (A Priori)", length(dropped_apriori)),
  stringsAsFactors = FALSE
)

if (length(target_markers) == 0) stop("[FATAL] All markers were filtered out! Check your config lists.")
mat_raw <- mat_raw[, target_markers, drop = FALSE]

# 3. Quality Control (Module Call)
# ------------------------------------------------------------------------------
# We pass the apriori dropped dataframe so it gets included in the final report
qc_result <- run_qc_pipeline(mat_raw, config$qc, dropped_apriori = df_dropped_apriori)

# Update local variables with clean data
mat_raw <- qc_result$data
qc_summary <- qc_result$report

# Update metadata to match filtered samples
valid_patients <- rownames(mat_raw)
raw_data <- raw_data %>% filter(Patient_ID %in% valid_patients)

# Save QC Report
out_dir <- file.path(config$output_root, "01_QC")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
save_qc_report(qc_summary, file.path(out_dir, "QC_Filtering_Report.xlsx"))

# 4. CoDa Pipeline Execution
# ------------------------------------------------------------------------------
# Zero Handling (Mask -> Fill -> CZM -> Restore)
mat_no_zeros <- coda_replace_zeros(mat_raw) # Assuming you might rename to standard convention

# NA Imputation
message("[CoDa] Handling remaining NAs (if any) via multiplicative replacement...")
mat_imputed <- impute_nas_simple(mat_no_zeros)

if (any(is.na(mat_imputed))) stop("[FATAL] Imputation failed: NAs remain.")
if (any(mat_imputed <= 0)) stop("[FATAL] Invalid values detected post-imputation.")

# Transformations
mat_clr <- coda_transform_clr(mat_imputed) # Rename if you applied standard naming
mat_ilr <- transform_ilr(mat_imputed)

# 5. Save Output
# ------------------------------------------------------------------------------
df_imputed <- cbind(raw_data[, meta_cols], as.data.frame(mat_imputed))
df_clr     <- cbind(raw_data[, meta_cols], as.data.frame(mat_clr))
df_ilr     <- cbind(raw_data[, meta_cols], as.data.frame(mat_ilr))

processed_data <- list(
  metadata     = raw_data[, meta_cols],
  markers      = colnames(mat_raw),      
  raw_matrix   = mat_raw,
  imputed_data = df_imputed,
  clr_data     = df_clr,
  ilr_data     = df_ilr,
  config       = config
)

saveRDS(processed_data, file.path(out_dir, "data_processed.rds"))
message("=== STEP 1 COMPLETE ===\n")