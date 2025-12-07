# ==============================================================================
# SCRIPT 01: DATA PREPARATION & NORMALIZATION
# PURPOSE: Centralized cleaning, imputation (kNN optimization), and transformation.
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
  library(ggplot2)        # For optimization plots
  library(Metrics)        # For RMSE calculation
})

# 2. CONFIGURATION -------------------------------------------------------------
CONFIG <- list(
  input_file = "/mnt/c/Users/Davide/Desktop/Davide/bandi/0dottorato/projects/long survivors/DB_anonimo_standardizzatov3.xlsx",
  output_dir = "/home/davidec/projects/compositional_analysis/processed_data",
  qc_dir = "/home/davidec/projects/compositional_analysis/results_fixed/01_QC",
  
  # Cohorts to process
  sheets = list(
    "NSCLC" = "NSCLC",
    "Healthy" = "Healthy_Donors",
    "HNSCC" = "HNSCC" 
  ),
  
  # QC Thresholds
  max_na_row = 0.55,  
  max_na_col = 0.20,
  
  # kNN Optimization Parameters
  optimize_k = TRUE,
  k_default = 3,
  k_range = 3:10
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
    cat(sprintf("   [WARNING] Could not load sheet '%s'. Skipping.\n", sheet))
  })
}

full_data <- bind_rows(df_list)
meta_cols <- c("Patient_ID", "Group")
marker_cols <- setdiff(names(full_data), meta_cols)

# Remove columns that are 100% NA (Safety check)
full_data <- full_data[, colSums(is.na(full_data)) < nrow(full_data)]
marker_cols <- setdiff(names(full_data), meta_cols)

# --- B. Quality Control (Filtering) ---
cat("[2] Quality Control & Filtering...\n")
df_num <- full_data[, marker_cols]

# Calc Missingness
row_na <- rowMeans(is.na(df_num))
col_na <- colMeans(is.na(df_num))

# Save QC Report
qc_df <- data.frame(ID = full_data$Patient_ID, Group = full_data$Group, Missing_Pct = row_na)
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

# --- C. Parameter Optimization (kNN) ---
cat("[3] Optimizing Imputation Parameters (kNN)...\n")
k_final <- CONFIG$k_default # Initialize with default

if(CONFIG$optimize_k) {
  # Setup data for simulation (use numeric columns)
  test_data <- full_data[, marker_cols]
  complete_rows <- complete.cases(test_data)
  ground_truth <- test_data[complete_rows, ]
  
  if(nrow(ground_truth) < 20) {
    cat("   [WARNING] Not enough complete rows (<20) to optimize k. Using default k=3.\n")
  } else {
    cat(sprintf("   -> Running simulation on %d complete rows...\n", nrow(ground_truth)))
    
    results_k <- data.frame()
    n_reps <- 5
    prop_na <- 0.05
    
    # Simulation Loop
    for(k in CONFIG$k_range) {
      errors <- numeric(n_reps)
      for(i in 1:n_reps) {
        set.seed(123 + i)
        # Create artificial NAs
        na_mat <- as.matrix(ground_truth)
        n_vals <- prod(dim(na_mat))
        miss_idx <- sample(n_vals, floor(n_vals * prop_na))
        actual_vals <- na_mat[miss_idx]
        na_mat[miss_idx] <- NA
        
        # Impute
        imputed_df <- VIM::kNN(as.data.frame(na_mat), k = k, imp_var = FALSE)
        
        # Calculate NRMSE
        imp_vals <- as.matrix(imputed_df)[miss_idx]
        rmse <- sqrt(mean((imp_vals - actual_vals)^2))
        errors[i] <- rmse / (max(actual_vals) - min(actual_vals))
      }
      avg_err <- mean(errors)
      results_k <- rbind(results_k, data.frame(k = k, NRMSE = avg_err))
      cat(sprintf("     -> k=%2d | Avg NRMSE: %.4f\n", k, avg_err))
    }
    
    # Select Best k
    best_k <- results_k$k[which.min(results_k$NRMSE)]
    k_final <- best_k
    cat(sprintf("   -> [RESULT] Optimal k detected: %d\n", k_final))
    
    # Save Plot
    p_k <- ggplot(results_k, aes(x=k, y=NRMSE)) +
      geom_line(color="#2E8B57") + geom_point(size=3) +
      geom_vline(xintercept=best_k, linetype="dashed", color="red") +
      theme_bw() + labs(title="kNN Parameter Optimization", y="Imputation Error (NRMSE)")
    ggsave(file.path(CONFIG$qc_dir, "kNN_Optimization.png"), p_k, width=6, height=4)
  }
}

# --- D. Imputation & Transformation ---
cat(sprintf("[4] Imputation & Transformation (Using k=%d)...\n", k_final))
df_vals <- full_data[, marker_cols]

# 1. kNN Imputation
if(any(is.na(df_vals))) {
  cat(sprintf("   -> Running kNN imputation (k=%d)...\n", k_final))
  df_vals <- VIM::kNN(df_vals, k = k_final, imp_var=FALSE)
} else {
  cat("   -> No missing values detected. Skipping imputation.\n")
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

# --- E. Save Master Object ---
cat("[5] Saving Master Data Object...\n")
master_data <- list(
  metadata = full_data[, meta_cols],
  raw_filtered = df_final_raw,      
  clr_transformed = df_clr_final,   
  markers = marker_cols,
  parameters = CONFIG,
  optimal_k = k_final # Save the k used for reference
)

saveRDS(master_data, file.path(CONFIG$output_dir, "clean_data.rds"))
cat("   -> Done! Saved to 'processed_data/clean_data.rds'\n")
