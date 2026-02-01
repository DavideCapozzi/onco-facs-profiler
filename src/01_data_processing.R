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
    "\n[CAUTION] Data max value is %.2f. Expected input is PERCENTAGES (0-100).\nIf this is actually Proportions (0-1), the pipeline will divide by 100 again, causing errors.\nPlease check your raw input file.", 
    max_val
  ))
} else {
  message(sprintf("[Check] Data range confirmed compatible with Percentage mode (Max=%.2f).", max_val))
}

# 3. Geometric-Aware Quality Control (Split Strategy)
# ------------------------------------------------------------------------------
message("[QC] Running Split QC Strategy...")

# --- STEP A: Basic QC (Missingness & Zero Variance on RAW DATA) ---
# Temporarily disable outlier detection to focus strictly on cleaning NAs/Zeros first.
# This prevents epsilon-imputation from masking high-missingness samples.
qc_config_basic <- config$qc
qc_config_basic$remove_outliers <- FALSE

message("   [QC-A] Filtering Missingness and Zero Variance on RAW data...")
qc_result <- run_qc_pipeline(
  mat_raw = mat_raw, 
  metadata = metadata_raw, 
  qc_config = qc_config_basic, 
  stratification_col = config$metadata$subgroup_col,
  dropped_markers_apriori = data.frame(),
  dropped_samples_apriori = dropped_samples_apriori
)

# Extract cleaned datasets from Step A
mat_clean_basic  <- qc_result$data
meta_clean_basic <- qc_result$metadata

# --- STEP B: Geometric Outlier Detection (on PROXY DATA) ---
# Now that NAs are handled, we generate the proxy for multivariate checks.
if (config$qc$remove_outliers) {
  message("   [QC-B] Generating geometric proxy for Outlier Detection...")
  
  # 1. Generate Proxy (Fast Mode) on CLEAN data
  proxy_results <- perform_hybrid_transformation(mat_clean_basic, config, mode = "fast")
  mat_proxy <- proxy_results$hybrid_data_z
  
  # 2. Detect Outliers using the Proxy
  message("   [QC-B] Detecting multivariate outliers (PCA-based)...")
  
  # Determine correct grouping column
  group_col <- if(config$metadata$subgroup_col %in% colnames(meta_clean_basic)) config$metadata$subgroup_col else "Group"
  
  is_outlier <- detect_pca_outliers(
    mat = mat_proxy, 
    groups = meta_clean_basic[[group_col]], 
    conf_level = config$qc$outlier_conf_level
  )
  
  # 3. Apply Outlier Filter & Update Report
  if (any(is_outlier)) {
    out_pids <- rownames(mat_clean_basic)[is_outlier]
    message(sprintf("   [QC-B] Dropping %d multivariate outliers.", length(out_pids)))
    
    # Capture info for report
    group_info <- as.character(meta_clean_basic[[group_col]][is_outlier])
    
    # Calculate NA % for these outliers (on the basic clean matrix for transparency)
    na_pcts <- rowMeans(is.na(mat_clean_basic[is_outlier, , drop=FALSE])) * 100
    
    # Append to existing report details
    new_drops <- data.frame(
      Patient_ID = out_pids,
      NA_Percent = round(na_pcts, 2),
      Reason = "Outlier",
      Original_Source = group_info,
      stringsAsFactors = FALSE
    )
    
    qc_result$report$dropped_rows_detail <- rbind(qc_result$report$dropped_rows_detail, new_drops)
    qc_result$report$n_row_dropped <- nrow(qc_result$report$dropped_rows_detail)
    
    # Apply removal
    mat_clean_basic  <- mat_clean_basic[!is_outlier, , drop = FALSE]
    meta_clean_basic <- meta_clean_basic[!is_outlier, , drop = FALSE]
    
    # Update Final Dimensions in Report
    qc_result$report$n_row_final <- nrow(mat_clean_basic)
    
    # Recalculate Final Breakdown
    if (!is.null(qc_result$report$breakdown_final)) {
      qc_result$report$breakdown_final <- table(meta_clean_basic[[config$metadata$subgroup_col]])
    }
  }
}

# --- STEP C: Synchronization ---
# Sync main objects with the final result of the Split QC
mat_raw <- mat_clean_basic
metadata_raw <- meta_clean_basic
raw_data <- raw_data %>% filter(Patient_ID %in% rownames(mat_raw))

message(sprintf("[QC] Final Dimensions: %d Samples x %d Markers", nrow(mat_raw), ncol(mat_raw)))

# Save QC Report
out_dir <- file.path(config$output_root, "01_data_processing")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
save_qc_report(qc_result$report, file.path(out_dir, "QC_Filtering_Report.xlsx"))

# 4. Final Hybrid Transformation (Complete Mode)
# ------------------------------------------------------------------------------
# Now running on strictly filtered data.
# Mode="complete" enables: CZM (zeros), BPCA (imputation), and CLR/ILR.

na_rates <- rowMeans(is.na(mat_raw))
impute_indices <- which(na_rates > 0)

if (length(impute_indices) > 0) {
  # Sort indices by missingness (descending) to highlight heavily imputed samples first
  impute_indices <- impute_indices[order(na_rates[impute_indices], decreasing = TRUE)]
  
  # Format string: "PatientID (XX%)"
  log_entries <- sprintf("%s (%.0f%%)", 
                         names(impute_indices), 
                         na_rates[impute_indices] * 100)
  
  message(sprintf("\n[Impute] The following %d samples contain missing values and will undergo BPCA Imputation:", length(log_entries)))
  message(paste(log_entries, collapse = ", "))
  message("") # Add spacing for readability
} else {
  message("\n[Impute] No missing values detected in the filtered dataset. BPCA step will be skipped internally.")
}

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