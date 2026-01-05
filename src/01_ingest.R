# src/01_ingest.R
# ==============================================================================
# STEP 01: DATA INGESTION & HYBRID CODA PREPROCESSING (REFACTORED V2)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(pcaMethods)
})

source("R/utils_io.R")     
source("R/modules_coda.R")
source("R/modules_qc.R")   

message("\n=== PIPELINE STEP 1: INGESTION & HYBRID TRANSFORM (BPCA) ===")

# 1. Load Config & Data
# ------------------------------------------------------------------------------
config <- load_config("config/global_params.yml")

if (is.null(config$hybrid_groups)) {
  stop("[FATAL] 'hybrid_groups' not found in config. Please update global_params.yml.")
}

raw_data <- load_raw_data(config)

# 2. Setup Matrices
# ------------------------------------------------------------------------------
subgroup_col <- config$metadata$subgroup_col
meta_cols <- c("Patient_ID", "Group")

if (!is.null(subgroup_col) && subgroup_col %in% names(raw_data)) {
  meta_cols <- c(meta_cols, subgroup_col)
}

marker_cols <- setdiff(names(raw_data), meta_cols)
mat_raw <- as.matrix(raw_data[, marker_cols])
rownames(mat_raw) <- raw_data$Patient_ID

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
qc_result <- run_qc_pipeline(mat_raw, config$qc, dropped_apriori = data.frame())
mat_raw <- qc_result$data
qc_summary <- qc_result$report
valid_patients <- rownames(mat_raw)
raw_data <- raw_data %>% filter(Patient_ID %in% valid_patients)

# Save QC Report
out_dir <- file.path(config$output_root, "01_QC")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
save_qc_report(qc_summary, file.path(out_dir, "QC_Filtering_Report.xlsx"))

# 4. Hybrid Logic: Split -> Transform -> Impute
# ------------------------------------------------------------------------------
message("\n[Transform] Starting Hybrid Strategy (Compositional vs Functional)...")

all_current_markers <- colnames(mat_raw)
comp_markers_list <- unique(unlist(config$hybrid_groups))
comp_markers <- intersect(all_current_markers, comp_markers_list)
func_markers <- setdiff(all_current_markers, comp_markers)

message(sprintf("   [Strategy] Split: %d Compositional (CLR) vs %d Functional (Logit+BPCA)", 
                length(comp_markers), length(func_markers)))

# Storage for transformed parts
mat_comp_trans <- NULL
mat_func_trans <- NULL

# --- TRACK A: Compositional Markers (CZM -> CLR) ---
if (length(comp_markers) > 0) {
  message("   [Track A] Processing Compositional Groups...")
  mat_comp_raw <- mat_raw[, comp_markers, drop = FALSE]
  
  # 1. Zero Replacement (CZM) - Handles zeros in simplex
  mat_comp_clean <- coda_replace_zeros(mat_comp_raw)
  
  # 2. Handle remaining NAs (Rare, but must be done before CLR)
  # Using multRepl (multiplicative replacement) for compositional NAs
  if (any(is.na(mat_comp_clean))) {
    message("      -> Imputing remaining NAs in composition (multRepl)...")
    mat_comp_clean <- zCompositions::multRepl(mat_comp_clean, label = NA, imp.missing = TRUE)
  }
  
  # 3. CLR Transformation (Group-wise or Global? Strategy: Global for storage, Local for ILR stats)
  # Here we transform ALL compositional markers. Note: CLR on the whole set of comp markers 
  # preserves distances if they are analyzed together.
  # However, technically CLR should be per-simplex. 
  # Strategy: We apply CLR to the block. 
  mat_comp_trans <- coda_transform_clr(mat_comp_clean)
}

# --- TRACK B: Functional Markers (LOD -> Logit -> BPCA) ---
if (length(func_markers) > 0) {
  message("   [Track B] Processing Functional Markers...")
  mat_func_raw <- mat_raw[, func_markers, drop = FALSE]
  
  # 1. Zero Handling (LOD Proxy: Min/2)
  # Zeros here are LOD (Limit of Detection), not structural zeros.
  mat_func_lod <- mat_func_raw
  for (j in 1:ncol(mat_func_lod)) {
    col_vals <- mat_func_lod[, j]
    zero_idx <- which(col_vals == 0 & !is.na(col_vals))
    
    if (length(zero_idx) > 0) {
      pos_vals <- col_vals[col_vals > 0 & !is.na(col_vals)]
      lod_proxy <- if(length(pos_vals)>0) min(pos_vals)/2 else 1e-6
      mat_func_lod[zero_idx, j] <- lod_proxy
    }
  }
  
  # 2. Logit Transformation
  # We transform BEFORE imputing because BPCA assumes continuous/normal-ish data.
  # NAs are preserved here.
  mat_func_logit <- coda_transform_logit(mat_func_lod)
  
  # 3. BPCA Imputation
  # This reconstructs missing values based on covariance with other markers
  mat_func_trans <- impute_matrix_bpca(mat_func_logit, nPcs = 3)
}

# 5. Merge & Normalize
# ------------------------------------------------------------------------------
message("\n[Merge] Recombining and creating final matrices...")

# Merge columns (handling cases where one track might be empty)
list_parts <- list()
if (!is.null(mat_comp_trans)) list_parts[[1]] <- mat_comp_trans
if (!is.null(mat_func_trans)) list_parts[[2]] <- mat_func_trans

if (length(list_parts) == 0) stop("[FATAL] No data remaining after filtering.")

mat_hybrid_raw <- do.call(cbind, list_parts)

# Restore original column order for consistency
mat_hybrid_raw <- mat_hybrid_raw[, sort(colnames(mat_hybrid_raw)), drop = FALSE]

# Create Z-Scored Version (for PCA/Network)
message("   [Norm] Applying Z-Score Standardization...")
mat_hybrid_z <- scale(mat_hybrid_raw)
attr(mat_hybrid_z, "scaled:center") <- NULL
attr(mat_hybrid_z, "scaled:scale") <- NULL
mat_hybrid_z <- as.matrix(mat_hybrid_z)

# 6. ILR Balances (For Stats Only)
# ------------------------------------------------------------------------------
# We need an imputed count matrix for ILR. 
# Reconstructing "imputed counts" from BPCA logits is possible (inv.logit) but complex.
# For simplicity and robustness, we re-run simple imputation on the merged count matrix 
# specifically for the ILR step, OR we accept that ILR uses the CoDa-cleaned data.
# Best approach: Use the cleaned compositional block from Track A.

ilr_list <- list()
if (!is.null(mat_comp_trans)) {
  # We need the POSITIVE data (before CLR) for ILR. 
  # Recalculating from mat_comp_clean (Track A step 2)
  ilr_list <- coda_compute_local_ilr(mat_comp_clean, config$hybrid_groups)
}

# 7. Save Output
# ------------------------------------------------------------------------------
df_hybrid_raw <- cbind(raw_data[, meta_cols], as.data.frame(mat_hybrid_raw))
df_hybrid_z   <- cbind(raw_data[, meta_cols], as.data.frame(mat_hybrid_z))

# We construct a synthetic "imputed_data" for reference (approximate)
# This is less critical than the hybrid matrices
mat_imputed_counts <- mat_raw # Placeholder
# (In a production pipeline, we would invert-logit the functional parts to fill this)

processed_data <- list(
  metadata        = raw_data[, meta_cols],
  markers         = colnames(mat_raw),      
  raw_matrix      = mat_raw,
  # Primary outputs for downstream:
  hybrid_markers  = colnames(mat_hybrid_z),
  hybrid_data_raw = df_hybrid_raw,  # Transformed (CLR/Logit), Imputed, NO Z-score
  hybrid_data_z   = df_hybrid_z,    # Transformed, Imputed, AND Z-scored
  ilr_balances    = ilr_list,       # For Permutation tests
  config          = config
)

saveRDS(processed_data, file.path(out_dir, "data_processed.rds"))
message("=== STEP 1 COMPLETE ===\n")