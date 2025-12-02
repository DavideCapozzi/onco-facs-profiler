# ==============================================================================
# SCRIPT 01: DATA PREPARATION & NORMALIZATION
# PURPOSE: Centralized cleaning, imputation, and transformation pipeline.
# INPUT: Raw Excel file
# OUTPUT: .rds object containing processed lists (Raw, Imputed, CLR)
# ==============================================================================

# 1. LIBRARIES -----------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(VIM)            # kNN Imputation
  library(zCompositions)  # Compositional Zero Handling
  library(openxlsx)
})

# 2. CONFIGURATION -------------------------------------------------------------
CONFIG <- list(
  input_file = "/mnt/c/Users/Davide/Desktop/Davide/bandi/0dottorato/projects/long survivors/DB_anonimo_standardizzatov2.xlsx",
  output_dir = "/home/davidec/projects/compositional_analysis/processed_data",
  qc_dir = "/home/davidec/projects/compositional_analysis/results/01_QC",
  
  # Cohorts to process (HNSCC included as requested)
  sheets = list(
    "NSCLC" = "NSCLC",
    "Healthy" = "Healthy_Donors",
    "HNSCC" = "HNSCC" 
  ),
  
  # QC Thresholds (RELAXED to visualize HNSCC)
  # Previous strict: 0.30. New relaxed: 0.55 to keep bad HNSCC samples for PCA check.
  max_na_row = 0.55,  
  max_na_col = 0.35   
)

# Create directories
if(!dir.exists(CONFIG$output_dir)) dir.create(CONFIG$output_dir, recursive = TRUE)
if(!dir.exists(CONFIG$qc_dir)) dir.create(CONFIG$qc_dir, recursive = TRUE)

cat("=== PIPELINE STEP 1: DATA PREPARATION ===\n")

# 3. HELPER FUNCTIONS ----------------------------------------------------------
clean_numeric <- function(x) {
  x <- as.character(x)
  x <- gsub(",", ".", x)
  x <- gsub("[<>]", "", x)
  x <- gsub("^(?i)(na|nd|n\\.a\\.|nan|missing|low)$", NA, x)
  suppressWarnings(as.numeric(x))
}

# 4. EXECUTION -----------------------------------------------------------------

# --- A. Load & Harmonize ---
cat("[1] Loading Data...\n")
df_list <- list()

for(grp_name in names(CONFIG$sheets)) {
  sheet <- CONFIG$sheets[[grp_name]]
  tryCatch({
    tmp <- read_excel(CONFIG$input_file, sheet = sheet)
    names(tmp)[1] <- "Patient_ID" # Normalize ID column
    
    # Clean numeric values
    tmp_clean <- tmp %>%
      mutate(across(-1, clean_numeric)) %>%
      mutate(Group = grp_name, .after = Patient_ID)
    
    df_list[[grp_name]] <- tmp_clean
    cat(sprintf("   -> Loaded %s: %d samples\n", grp_name, nrow(tmp_clean)))
  }, error = function(e) {
    # If a sheet is missing, warn but don't stop
    cat(sprintf("   [WARNING] Could not load sheet '%s'. Skipping.\n", sheet))
  })
}

full_data <- bind_rows(df_list)
meta_cols <- c("Patient_ID", "Group")
marker_cols <- setdiff(names(full_data), meta_cols)
# Remove columns that are 100% NA
full_data <- full_data[, colSums(is.na(full_data)) < nrow(full_data)]
marker_cols <- setdiff(names(full_data), meta_cols)

# --- B. Quality Control (Filtering) ---
cat("[2] Quality Control & Filtering...\n")
df_num <- full_data[, marker_cols]

# Calc Missingness
row_na <- rowMeans(is.na(df_num))
col_na <- colMeans(is.na(df_num))

# Save QC Report
qc_df <- data.frame(
  ID = full_data$Patient_ID,
  Group = full_data$Group,
  Missing_Pct = row_na
)
write.xlsx(qc_df, file.path(CONFIG$qc_dir, "Missingness_Report.xlsx"))

# Filter Rows (Patients)
bad_rows <- which(row_na > CONFIG$max_na_row)
if(length(bad_rows) > 0) {
  cat(sprintf("   [DROP] Removed %d patients (>%.0f%% NA): %s\n", 
              length(bad_rows), CONFIG$max_na_row*100, paste(full_data$Patient_ID[bad_rows], collapse=", ")))
  full_data <- full_data[-bad_rows, ]
}

# Filter Cols (Markers)
bad_cols <- which(col_na > CONFIG$max_na_col)
if(length(bad_cols) > 0) {
  cat(sprintf("   [DROP] Removed %d markers (>%.0f%% NA): %s\n", 
              length(bad_cols), CONFIG$max_na_col*100, paste(names(bad_cols), collapse=", ")))
  full_data <- full_data[, -which(names(full_data) %in% names(bad_cols))]
}

# Update marker list
marker_cols <- setdiff(names(full_data), meta_cols)
df_final_raw <- full_data 

# --- C. Imputation & Transformation ---
cat("[3] Imputation & Transformation...\n")
df_vals <- full_data[, marker_cols]

# 1. kNN Imputation (Crucial for HNSCC which has holes)
if(any(is.na(df_vals))) {
  cat("   -> Running kNN imputation (k=3)...\n")
  df_vals <- VIM::kNN(df_vals, k=3, imp_var=FALSE)
}

# 2. CZM (for Zeros)
if(any(df_vals == 0, na.rm=TRUE)) {
  cat("   -> Replacing zeros (CZM method)...\n")
  df_clean <- zCompositions::cmultRepl(df_vals, method="CZM", output="p-counts", suppress.print=TRUE)
} else {
  df_clean <- df_vals
}

# 3. CLR Transformation
cat("   -> Applying CLR transformation...\n")
gmean <- function(x) exp(mean(log(x)))
df_clr_vals <- t(apply(df_clean, 1, function(x) log(x / gmean(x))))
df_clr_vals <- as.data.frame(df_clr_vals)
colnames(df_clr_vals) <- colnames(df_clean)

# Combine back with metadata
df_clr_final <- cbind(full_data[, meta_cols], df_clr_vals)

# --- D. Save Master Object ---
cat("[4] Saving Master Data Object...\n")
master_data <- list(
  metadata = full_data[, meta_cols],
  raw_filtered = df_final_raw,      
  clr_transformed = df_clr_final,   
  markers = marker_cols,
  parameters = CONFIG
)

saveRDS(master_data, file.path(CONFIG$output_dir, "clean_data.rds"))
cat("   -> Done! Saved to 'processed_data/clean_data.rds'\n")
