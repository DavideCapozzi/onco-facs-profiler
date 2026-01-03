# src/01_ingest.R
# ==============================================================================
# STEP 01: DATA INGESTION & HYBRID CODA PREPROCESSING
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
})

source("R/utils_io.R")     
source("R/modules_coda.R")
source("R/modules_qc.R")   

message("\n=== PIPELINE STEP 1: INGESTION & HYBRID TRANSFORM ===")

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
  message(sprintf("[Data] Tracking extra metadata column: %s", subgroup_col))
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

# 4. Hybrid Transformation Logic & Zero Handling
# ------------------------------------------------------------------------------
message("\n[Transform] Starting Configuration-Driven Hybrid Strategy...")

# Identify Marker Types based on Config
all_current_markers <- colnames(mat_raw)
comp_markers_list <- unique(unlist(config$hybrid_groups))
comp_markers <- intersect(all_current_markers, comp_markers_list)
func_markers <- setdiff(all_current_markers, comp_markers)

message(sprintf("   [Strategy] Split Imputation: %d Compositional (CZM) vs %d Functional (LOD-proxy)", 
                length(comp_markers), length(func_markers)))

# --- A. Compositional Imputation (CZM) ---
# Applies only to markers defined in hybrid_groups to preserve sub-compositional geometry
if (length(comp_markers) > 0) {
  mat_comp_raw <- mat_raw[, comp_markers, drop = FALSE]
  # coda_replace_zeros uses zCompositions::cmultRepl internally
  mat_comp_imp <- coda_replace_zeros(mat_comp_raw)
} else {
  mat_comp_imp <- matrix(0, nrow = nrow(mat_raw), ncol = 0)
}

# --- B. Functional Imputation (Adaptive Min/2) ---
# Applies to independent functional markers. Replaces 0 with half the minimum positive value per marker.
# This mimics the Limit of Detection (LOD) without enforcing compositional constraints.
if (length(func_markers) > 0) {
  mat_func_raw <- mat_raw[, func_markers, drop = FALSE]
  mat_func_imp <- mat_func_raw
  
  for (j in 1:ncol(mat_func_imp)) {
    col_vals <- mat_func_imp[, j]
    # Identify zeros that are NOT NAs
    zero_idx <- which(col_vals == 0 & !is.na(col_vals))
    
    if (length(zero_idx) > 0) {
      # Calculate adaptive LOD: Half of the minimum positive value
      pos_vals <- col_vals[col_vals > 0 & !is.na(col_vals)]
      
      if (length(pos_vals) > 0) {
        lod_proxy <- min(pos_vals) / 2
        mat_func_imp[zero_idx, j] <- lod_proxy
      } else {
        # Fallback if column is all zeros (should be caught by QC, but purely defensive)
        mat_func_imp[zero_idx, j] <- 1e-6
      }
    }
  }
} else {
  mat_func_imp <- matrix(0, nrow = nrow(mat_raw), ncol = 0)
}

# --- C. Merge & Finalize Imputation ---
# Recombine matrices and handle any remaining true NAs (missing data)
mat_merged <- cbind(mat_comp_imp, mat_func_imp)
# Restore original column order for consistency
mat_merged <- mat_merged[, all_current_markers, drop = FALSE]

# Impute genuine NAs (not zeros) that might remain
mat_imputed <- impute_nas_simple(mat_merged)

# 5. Transformation (CLR & Logit)
# ------------------------------------------------------------------------------
all_markers <- colnames(mat_imputed)
transformed_list <- list()
assigned_markers <- c()

# A. Process Compositional Groups (CLR)
for (group_name in names(config$hybrid_groups)) {
  
  requested_mks <- config$hybrid_groups[[group_name]]
  present_mks <- intersect(requested_mks, all_markers)
  
  if (length(present_mks) == 0) {
    warning(sprintf("   [WARN] Group '%s' configured but no markers found. Skipping.", group_name))
    next
  }
  
  mat_sub <- mat_imputed[, present_mks, drop = FALSE]
  
  # Logic: CLR requires > 1 marker.
  if (ncol(mat_sub) > 1) {
    message(sprintf("   [Group: %s] %d markers found -> Applying CLR", group_name, length(present_mks)))
    mat_trans <- coda_transform_clr(mat_sub)
  } else {
    message(sprintf("   [Group: %s] 1 marker found (%s) -> Fallback to Log1p", group_name, present_mks))
    mat_trans <- coda_transform_logit(mat_sub)
  }
  
  transformed_list[[group_name]] <- mat_trans
  assigned_markers <- c(assigned_markers, present_mks)
}

# B. Process Functional/Remaining Markers (Logit)
# Note: These have been zero-imputed via Min/2 strategy above, suitable for Logit.
unassigned_markers <- setdiff(all_markers, assigned_markers)

if (length(unassigned_markers) > 0) {
  message(sprintf("   [Group: Functional/Other] %d markers -> Applying Logit", length(unassigned_markers)))
  mat_sub <- mat_imputed[, unassigned_markers, drop = FALSE]
  mat_trans <- coda_transform_logit(mat_sub)
  
  transformed_list[["Functional"]] <- mat_trans
}

# 6. Reconstruct & Normalization
# ------------------------------------------------------------------------------
# A. Merge pieces
mat_hybrid_tmp <- do.call(cbind, transformed_list)

# B. Clean column names (Remove "Group." prefix if added by list binding)
colnames(mat_hybrid_tmp) <- sub("^[^.]+\\.", "", colnames(mat_hybrid_tmp))
mat_hybrid_tmp <- mat_hybrid_tmp[, sort(colnames(mat_hybrid_tmp)), drop = FALSE]

# C. Save RAW Hybrid (Transformed but NOT Z-scored) - Useful for univariate tests
mat_hybrid_raw <- mat_hybrid_tmp

# D. Save Z-SCORED Hybrid - Useful for PCA, Heatmaps, Networks
message("   [Norm] Applying Z-Score Standardization to Hybrid Matrix...")
mat_hybrid_z <- scale(mat_hybrid_raw)
# ensure it's a plain matrix:
attr(mat_hybrid_z, "scaled:center") <- NULL
attr(mat_hybrid_z, "scaled:scale") <- NULL
mat_hybrid_z <- as.matrix(mat_hybrid_z)

# E. Compute Local ILR (For Statistical Inference only)
ilr_list <- coda_compute_local_ilr(mat_imputed, config$hybrid_groups)


# 7. Save Output
# ------------------------------------------------------------------------------
df_hybrid_raw <- cbind(raw_data[, meta_cols], as.data.frame(mat_hybrid_raw))
df_hybrid_z   <- cbind(raw_data[, meta_cols], as.data.frame(mat_hybrid_z))

processed_data <- list(
  metadata        = raw_data[, meta_cols],
  markers         = colnames(mat_raw),      
  raw_matrix      = mat_raw,
  imputed_data    = cbind(raw_data[, meta_cols], as.data.frame(mat_imputed)),
  hybrid_markers  = colnames(mat_hybrid_z),
  hybrid_data_raw = df_hybrid_raw,  # Transformed (CLR/Log), NO Z-score
  hybrid_data_z   = df_hybrid_z,    # Transformed AND Z-scored
  ilr_balances    = ilr_list,
  config          = config
)

saveRDS(processed_data, file.path(out_dir, "data_processed.rds"))
message("=== STEP 1 COMPLETE  ===\n")