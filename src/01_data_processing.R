# src/01_data_processing.R
# ==============================================================================
# STEP 01: DATA PROCESSING
# Description: Data ingestion + QC + hybrid coda preprocessing
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(pcaMethods)
  library(zCompositions) 
})

source("R/utils_io.R")     
source("R/modules_coda.R") 
source("R/modules_qc.R")   

message("\n=== PIPELINE STEP 1: INGESTION + QC + HYBRID TRANSFORM ===")

# 1. Load Configuration & Data
# ------------------------------------------------------------------------------
config <- load_config("config/global_params.yml")

if (is.null(config$hybrid_groups)) {
  stop("[FATAL] 'hybrid_groups' not found in config. Please update global_params.yml.")
}

raw_data <- load_raw_data(config)

# 2. Setup Initial Matrices
# ------------------------------------------------------------------------------
# Split Metadata and Numeric Matrix
subgroup_col <- config$metadata$subgroup_col
meta_cols <- c("Patient_ID", "Group")
if (!is.null(subgroup_col)) meta_cols <- c(meta_cols, subgroup_col)

# Extract Matrix and Metadata aligned
mat_raw <- as.matrix(raw_data[, setdiff(names(raw_data), meta_cols)])
rownames(mat_raw) <- raw_data$Patient_ID
metadata_raw <- raw_data[, meta_cols]

message(sprintf("[Data] Initial Matrix: %d Samples x %d Markers", nrow(mat_raw), ncol(mat_raw)))

# --- Sample Filtering (Blacklist) ---
target_samples <- rownames(mat_raw)
dropped_samples_apriori <- data.frame()

if (!is.null(config$sample_selection$blacklist) && length(config$sample_selection$blacklist) > 0) {
  blacklist_samples <- config$sample_selection$blacklist
  
  # Identify samples to drop that actually exist in data
  samples_to_drop <- intersect(target_samples, blacklist_samples)
  
  if (length(samples_to_drop) > 0) {
    message(sprintf("[Data] Dropping %d samples defined in blacklist.", length(samples_to_drop)))
    
    # Store details for QC report
    dropped_indices <- which(rownames(mat_raw) %in% samples_to_drop)
    
    dropped_samples_apriori <- data.frame(
      Patient_ID = samples_to_drop,
      NA_Percent = NA, 
      Reason = "A Priori",
      Original_Source = if (subgroup_col %in% colnames(metadata_raw)) as.character(metadata_raw[dropped_indices, subgroup_col]) else "N/A",
      stringsAsFactors = FALSE
    )
    
    # Apply filter
    target_samples <- setdiff(target_samples, samples_to_drop)
    mat_raw <- mat_raw[target_samples, , drop = FALSE]
    metadata_raw <- metadata_raw[match(target_samples, metadata_raw$Patient_ID), , drop = FALSE]
    raw_data <- raw_data %>% filter(Patient_ID %in% target_samples)
  }
}

# --- Marker Filtering (Whitelist/Blacklist) ---
target_markers <- colnames(mat_raw)
if (!is.null(config$marker_selection$whitelist) && length(config$marker_selection$whitelist) > 0) {
  target_markers <- intersect(target_markers, config$marker_selection$whitelist)
}
if (!is.null(config$marker_selection$blacklist) && length(config$marker_selection$blacklist) > 0) {
  target_markers <- setdiff(target_markers, config$marker_selection$blacklist)
}

mat_raw <- mat_raw[, target_markers, drop = FALSE]

# 2.1 Data Assumption Validation
# ------------------------------------------------------------------------------
# We strictly assume input data are PERCENTAGES (0-100) 
# This check ensures we don't accidentally process Proportions (0-1) which would
# become tiny numbers when divided by 100.
max_val <- max(mat_raw, na.rm = TRUE)

if (max_val <= 1.05) { # Tolerance slightly above 1.0
  warning(sprintf(
    "\n[CAUTION] Data max value is %.2f. You stated input is PERCENTAGES (0-100).\nIf this is actually Proportions (0-1), the pipeline will divide by 100 again, causing errors.\nPlease check your raw input file.", 
    max_val
  ))
} else {
  message(sprintf("[Check] Data range confirmed compatible with Percentage mode (Max=%.2f).", max_val))
}

# 3. Geometric-Aware Quality Control
# ------------------------------------------------------------------------------
message("[QC] Generating geometric proxy for robust outlier detection...")

# A. Create Proxy Matrix (Fast Mode)
# Uses `perform_hybrid_transformation` which correctly handles percentages
# and uses robust epsilon replacement for QC purposes.
proxy_results <- perform_hybrid_transformation(mat_raw, config, mode = "fast")
mat_qc_proxy  <- proxy_results$hybrid_data_z

# B. Run QC on the Proxy Matrix
# This step identifies rows/cols to drop.
qc_result_proxy <- run_qc_pipeline(
  mat_raw = mat_qc_proxy, 
  metadata = metadata_raw, 
  qc_config = config$qc, 
  stratification_col = config$metadata$subgroup_col,
  dropped_markers_apriori = data.frame(), # If tracked previously, pass it here
  dropped_samples_apriori = dropped_samples_apriori
)

# C. Apply Filters to ORIGINAL Raw Data (Synchronization)
# this step ensures that `mat_raw` is clean before proceeding.
valid_patients <- rownames(qc_result_proxy$data)
valid_markers  <- colnames(qc_result_proxy$data)

dropped_pats <- setdiff(rownames(mat_raw), valid_patients)
dropped_mks  <- setdiff(colnames(mat_raw), valid_markers)

if(length(dropped_pats) > 0) message(sprintf("[QC] Removed %d patients (outlier/missingness).", length(dropped_pats)))
if(length(dropped_mks) > 0)  message(sprintf("[QC] Removed %d markers (low variance/missingness).", length(dropped_mks)))

# Update Raw Objects
mat_raw <- mat_raw[valid_patients, valid_markers, drop = FALSE]
raw_data <- raw_data %>% filter(Patient_ID %in% valid_patients)
metadata_raw <- metadata_raw[match(valid_patients, metadata_raw$Patient_ID), , drop = FALSE]

# Save QC Report
out_dir <- file.path(config$output_root, "01_data_processing")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
save_qc_report(qc_result_proxy$report, file.path(out_dir, "QC_Filtering_Report.xlsx"))

# 4. Final Hybrid Transformation (Complete Mode)
# ------------------------------------------------------------------------------
# Now running on strictly filtered data.
# Mode="complete" enables: CZM (zeros), BPCA (imputation), and CLR/ILR.

message("[CoDa] Running Final Hybrid Transformation (Complete Mode) on clean data...")
transform_results <- perform_hybrid_transformation(mat_raw, config, mode = "complete")

# Extract results
mat_hybrid_raw <- transform_results$hybrid_data_raw
mat_hybrid_z   <- transform_results$hybrid_data_z
ilr_list       <- transform_results$ilr_balances
final_markers  <- transform_results$hybrid_markers

# 5. Final Safety Polish
# ------------------------------------------------------------------------------
# Fallback for any residual NAs (e.g., if BPCA didn't converge for a specific point).
if (any(is.na(mat_hybrid_z))) {
  n_na_res <- sum(is.na(mat_hybrid_z))
  message(sprintf("[WARN] Detected %d residual NAs in final matrix. Applying median fallback...", n_na_res))
  
  mat_hybrid_z <- apply(mat_hybrid_z, 2, function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    return(x)
  })
}

# Check ILR balances
for (g_name in names(ilr_list)) {
  if (any(is.na(ilr_list[[g_name]]))) {
    ilr_list[[g_name]][is.na(ilr_list[[g_name]])] <- 0 
  }
}

# 6. Save Final Output
# ------------------------------------------------------------------------------
df_hybrid_raw <- cbind(raw_data[, meta_cols], as.data.frame(mat_hybrid_raw))
df_hybrid_z   <- cbind(raw_data[, meta_cols], as.data.frame(mat_hybrid_z))

processed_data <- list(
  metadata        = raw_data[, meta_cols],
  markers         = colnames(mat_raw),      
  raw_matrix      = mat_raw, 
  hybrid_markers  = final_markers,
  hybrid_data_raw = df_hybrid_raw,
  hybrid_data_z   = df_hybrid_z,
  ilr_balances    = ilr_list,
  config          = config
)

saveRDS(processed_data, file.path(out_dir, "data_processed.rds"))
message("=== STEP 1 COMPLETE ===\n")