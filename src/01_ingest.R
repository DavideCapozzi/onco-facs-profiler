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
config <- load_config("config/global_params.yml")

# CHECK DI SICUREZZA: Verifica che hybrid_groups esista
if (is.null(config$hybrid_groups)) {
  stop("[FATAL] 'hybrid_groups' not found in config. Please update global_params.yml.")
}

raw_data <- load_raw_data(config)

# 2. Setup Matrices
meta_cols <- c("Patient_ID", "Group")
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
qc_result <- run_qc_pipeline(mat_raw, config$qc, dropped_apriori = data.frame())
mat_raw <- qc_result$data
qc_summary <- qc_result$report
valid_patients <- rownames(mat_raw)
raw_data <- raw_data %>% filter(Patient_ID %in% valid_patients)

# Save QC Report
out_dir <- file.path(config$output_root, "01_QC")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
save_qc_report(qc_summary, file.path(out_dir, "QC_Filtering_Report.xlsx"))

# 4. Hybrid Transformation Logic
# ------------------------------------------------------------------------------
message("\n[Transform] Starting Configuration-Driven Hybrid Strategy...")

# A. Global Zero/NA Handling
mat_no_zeros <- coda_replace_zeros(mat_raw) 
mat_imputed  <- impute_nas_simple(mat_no_zeros)

# B. Split & Transform based on Config
all_markers <- colnames(mat_imputed)
transformed_list <- list()
assigned_markers <- c()

# Loop through Compositional Groups defined in Config
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
    message(sprintf("   [Group: %s] %d markers found -> Applying CLR (Compositional)", group_name, length(present_mks)))
    # Output visualizzato: NAIVE, CM, EFF, EM
    mat_trans <- coda_transform_clr(mat_sub)
  } else {
    message(sprintf("   [Group: %s] Only 1 marker found (%s) -> Fallback to Log1p", group_name, present_mks))
    mat_trans <- log1p(mat_sub)
  }
  
  transformed_list[[group_name]] <- mat_trans
  assigned_markers <- c(assigned_markers, present_mks)
}

# C. Process Functional/Remaining Markers (Log1p)
unassigned_markers <- setdiff(all_markers, assigned_markers)

if (length(unassigned_markers) > 0) {
  message(sprintf("   [Group: Functional/Other] %d markers -> Applying Log1p", length(unassigned_markers)))
  mat_sub <- mat_imputed[, unassigned_markers, drop = FALSE]
  mat_trans <- log1p(mat_sub)
  
  transformed_list[["Functional"]] <- mat_trans
}

# D. Reconstruct Hybrid Matrix
# do.call crea nomi tipo "Differentiation.NAIVE". Li dobbiamo pulire.
mat_hybrid <- do.call(cbind, transformed_list)

# FIX IMPORTANTE: Rimuovi i prefissi "Group." dai nomi delle colonne
# Questo garantisce che i nomi siano "NAIVE", "CD3", ecc., compatibili con gli script successivi.
colnames(mat_hybrid) <- sub("^[^.]+\\.", "", colnames(mat_hybrid))

# Riordina alfabeticamente
mat_hybrid <- mat_hybrid[, sort(colnames(mat_hybrid)), drop = FALSE]

# 5. Save Output
df_hybrid <- cbind(raw_data[, meta_cols], as.data.frame(mat_hybrid))

processed_data <- list(
  metadata     = raw_data[, meta_cols],
  markers      = colnames(mat_raw),      
  raw_matrix   = mat_raw,
  imputed_data = cbind(raw_data[, meta_cols], as.data.frame(mat_imputed)),
  clr_data     = df_hybrid,   # HYBRID MATRIX (Mixed CLR & Log1p)
  ilr_data     = df_hybrid,   # Placeholder
  config       = config
)

saveRDS(processed_data, file.path(out_dir, "data_processed.rds"))
message("=== STEP 1 COMPLETE (Hybrid) ===\n")