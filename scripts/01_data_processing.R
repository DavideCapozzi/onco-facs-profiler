# scripts/01_data_processing.R
# ==============================================================================
# STEP 01: DATA INGESTION, IMPUTATION & TRANSFORMATION
# ==============================================================================

# 1. Setup Environment
suppressPackageStartupMessages({
  library(yaml)
  library(readxl)
  library(dplyr)
  library(openxlsx)
  library(robCompositions)
})

# Load Configuration & Utilities
config <- read_yaml("config/global_params.yml")
source("R/utils_coda.R") # Loads CoDa functions (median fix, clr)

# Create Output Directories
processed_dir <- file.path("data", "processed") # As requested
qc_dir <- file.path(config$output_root, "01_QC")

if(!dir.exists(processed_dir)) dir.create(processed_dir, recursive = TRUE)
if(!dir.exists(qc_dir)) dir.create(qc_dir, recursive = TRUE)

message("=== PIPELINE STEP 1: DATA PREPARATION ===")

# ------------------------------------------------------------------------------
# 2. Load & Merge Data
# ------------------------------------------------------------------------------
message("[1] Loading Data from External Path...")
if(!file.exists(config$input_file)) stop("Input file not found check config.yml")

df_list <- list()

# Iterate over cohorts defined in config
for(cohort_name in names(config$cohorts)) {
  sheet_name <- config$cohorts[[cohort_name]]
  
  tryCatch({
    # Read Excel
    raw_tmp <- read_excel(config$input_file, sheet = sheet_name)
    
    # Standardize ID column (Assume 1st col is ID)
    colnames(raw_tmp)[1] <- "Patient_ID"
    
    # Force numeric cleaning
    clean_tmp <- raw_tmp %>%
      mutate(across(-1, ~suppressWarnings(as.numeric(as.character(.))))) %>%
      mutate(Group = cohort_name, .after = Patient_ID)
    
    df_list[[cohort_name]] <- clean_tmp
    message(sprintf("   -> Loaded %s: %d samples", cohort_name, nrow(clean_tmp)))
    
  }, error = function(e) {
    warning(sprintf("   [WARN] Could not load sheet '%s': %s", sheet_name, e$message))
  })
}

full_data <- bind_rows(df_list)
meta_cols <- c("Patient_ID", "Group")
marker_cols <- setdiff(names(full_data), meta_cols)

# ------------------------------------------------------------------------------
# 3. Quality Control (Filtering)
# ------------------------------------------------------------------------------
message("[2] Running QC & Filtering...")
df_vals <- full_data[, marker_cols]

# Calculate Missingness
row_na <- rowMeans(is.na(df_vals))
col_na <- colMeans(is.na(df_vals))

# Apply Thresholds from Config
bad_rows <- which(row_na > config$qc$max_na_row_pct)
bad_cols <- which(col_na > config$qc$max_na_col_pct)

# Save QC Report
qc_report <- data.frame(Patient_ID = full_data$Patient_ID, Missing_Pct = row_na)
write.xlsx(qc_report, file.path(qc_dir, "Missingness_Report.xlsx"))

# Filter
if(length(bad_rows) > 0) {
  message(sprintf("   [DROP] Removed %d patients (>%s%% NA)", length(bad_rows), config$qc$max_na_row_pct*100))
  full_data <- full_data[-bad_rows, ]
}
if(length(bad_cols) > 0) {
  message(sprintf("   [DROP] Removed %d markers (>%s%% NA)", length(bad_cols), config$qc$max_na_col_pct*100))
  full_data <- full_data[, -which(names(full_data) %in% names(bad_cols))]
}

# Update marker list after filtering
marker_cols <- setdiff(names(full_data), meta_cols)
message(sprintf("   -> Final Dimensions: %d Patients x %d Markers", nrow(full_data), length(marker_cols)))

# ------------------------------------------------------------------------------
# 4. Imputation & Transformation (The CoDa Engine)
# ------------------------------------------------------------------------------
message("[3] Running CoDa Imputation (Strategy: Median-Fix -> CZM -> kNN)...")

# We process the numeric matrix. 
# NOTE: 'run_complex_imputation' (from utils_coda.R) handles the median injection internally.
imputed_matrix <- run_complex_imputation(
  data_df = full_data[, marker_cols], 
  k = config$imputation$knn_k_default
)

message("[4] Applying CLR Transformation...")
# 'apply_clr' handles geometric mean and log-ratio
clr_matrix <- apply_clr(imputed_matrix)

# ------------------------------------------------------------------------------
# 5. Save Master Object
# ------------------------------------------------------------------------------
message("[5] Saving Processed Data...")

# Reconstruct DataFrames with Metadata
df_imputed <- cbind(full_data[, meta_cols], as.data.frame(imputed_matrix))
df_clr <- cbind(full_data[, meta_cols], as.data.frame(clr_matrix))

master_list <- list(
  metadata = full_data[, meta_cols],
  markers = marker_cols,
  raw_filtered = full_data,
  imputed = df_imputed,
  clr = df_clr,
  config_used = config
)

saveRDS(master_list, file.path(processed_dir, "clean_data.rds"))
message("=== DONE. Data saved to 'data/processed/clean_data.rds' ===")
