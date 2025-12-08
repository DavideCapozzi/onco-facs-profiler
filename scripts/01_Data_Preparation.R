# ==============================================================================
# SCRIPT 01: DATA PREPARATION & NORMALIZATION
# PURPOSE: Centralized cleaning, CoDa imputation (robCompositions), and transformation.
# INPUT: Raw Excel file
# OUTPUT: .rds object containing processed lists (Raw, Imputed, CLR)
# ==============================================================================

# 1. LIBRARIES -----------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(robCompositions) # For Aitchison kNN Imputation (impKNNa)
  library(zCompositions)   # For Bayesian Zero Replacement (cmultRepl)
  library(openxlsx)
  library(ggplot2)         # For optimization plots
  library(Metrics)         # For RMSE calculation
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
  
  # kNN Optimization Parameters (robCompositions)
  optimize_k = TRUE,
  k_default = 3,
  k_range = 3:10
)

# Create directories
if(!dir.exists(CONFIG$output_dir)) dir.create(CONFIG$output_dir, recursive = TRUE)
if(!dir.exists(CONFIG$qc_dir)) dir.create(CONFIG$qc_dir, recursive = TRUE)

cat("=== PIPELINE STEP 1: DATA PREPARATION (CoDa Aware) ===\n")

# 3. HELPER FUNCTIONS ----------------------------------------------------------
clean_numeric <- function(x) {
  x <- as.character(x)
  x <- gsub(",", ".", x)
  x <- gsub("[<>]", "", x)
  x <- gsub("^(?i)(na|nd|n\\.a\\.|nan|missing|low)$", NA, x)
  suppressWarnings(as.numeric(x))
}

# Helper to format QC DROP messages with newlines and percentages
format_drop_msg <- function(indices, names_vec, pcts) {
  items <- paste0(names_vec[indices], " (", round(pcts[indices]*100, 1), "%)")
  paste(items, collapse = ",\n")
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
  msg_rows <- format_drop_msg(bad_rows, full_data$Patient_ID, row_na)
  cat(sprintf("   [DROP] Removed %d patients (>%.0f%% NA):\n%s\n", 
              length(bad_rows), CONFIG$max_na_row*100, msg_rows))
  full_data <- full_data[-bad_rows, ]
}

# Filter Cols (Markers)
bad_cols <- which(col_na > CONFIG$max_na_col)
if(length(bad_cols) > 0) {
  msg_cols <- format_drop_msg(bad_cols, names(col_na), col_na)
  cat(sprintf("   [DROP] Removed %d markers (>%.0f%% NA):\n%s\n", 
              length(bad_cols), CONFIG$max_na_col*100, msg_cols))
  full_data <- full_data[, -which(names(full_data) %in% names(bad_cols))]
}

# Update marker list
marker_cols <- setdiff(names(full_data), meta_cols)
df_final_raw <- full_data 

# --- C. Zero Handling (Strategy: Bridge Imputation -> CZM) ---
# We handle Zeros BEFORE imputation because impKNNa (Aitchison) fails on zeros.
# We handle Zeros BEFORE optimizing k, so the optimization runs on valid CoDa data.

cat("[3] Handling Zeros (CZM)...\n")
df_vals <- full_data[, marker_cols]

if(any(df_vals == 0, na.rm=TRUE)) {
  cat("   -> Detected zeros. Applying CZM with temporary NA handling...\n")
  
  # 1. Identify NAs positions
  na_mat <- is.na(df_vals)
  
  # 2. Temporary fill NAs (e.g., column median) to allow CZM to run
  # (CZM requires a complete matrix to estimate proportions correctly)
  df_temp <- df_vals
  for(j in 1:ncol(df_temp)) {
    # Use median of observed values as a neutral placeholder
    col_med <- median(df_temp[[j]], na.rm=TRUE)
    if(is.na(col_med)) col_med <- 0.001 # Fallback if col is empty (unlikely after QC)
    df_temp[is.na(df_temp[[j]]), j] <- col_med
  }
  
  # 3. Apply CZM (Bayesian-Multiplicative Replacement)
  # Now df_temp has NO NAs, so CZM runs optimally.
  df_no_zeros <- zCompositions::cmultRepl(df_temp, method="CZM", output="p-counts", suppress.print=TRUE)
  
  # 4. Restore NAs
  # We put the NAs back exactly where they were. 
  # Now we have a dataset with NO Zeros (handled by CZM) and Real NAs.
  df_ready_for_knn <- df_no_zeros
  df_ready_for_knn[na_mat] <- NA
  
} else {
  cat("   -> No zeros detected. Proceeding to imputation.\n")
  df_ready_for_knn <- df_vals
}

# --- D. Parameter Optimization (robCompositions::impKNNa) ---
cat("[4] Optimizing Imputation Parameters (Aitchison kNN)...\n")
k_final <- CONFIG$k_default 

if(CONFIG$optimize_k) {
  # Use the Zero-Free dataset for simulation
  test_data <- df_ready_for_knn
  complete_rows <- complete.cases(test_data)
  ground_truth <- test_data[complete_rows, ]
  
  if(nrow(ground_truth) < 20) {
    cat("   [WARNING] Not enough complete rows (<20) to optimize k. Using default k=3.\n")
  } else {
    cat(sprintf("   -> Running simulation on %d complete rows...\n", nrow(ground_truth)))
    
    results_k <- data.frame()
    n_reps <- 5
    prop_na <- 0.05
    
    for(k in CONFIG$k_range) {
      errors <- numeric(n_reps)
      for(i in 1:n_reps) {
        set.seed(123 + i)
        
        # Create artificial NAs
        na_mat_sim <- as.matrix(ground_truth)
        n_vals <- prod(dim(na_mat_sim))
        miss_idx <- sample(n_vals, floor(n_vals * prop_na))
        actual_vals <- na_mat_sim[miss_idx]
        na_mat_sim[miss_idx] <- NA
        
        tryCatch({
          # Impute using impKNNa (works because ground_truth has no zeros)
          # Removed 'print=FALSE' as it is not supported by robCompositions::impKNNa
          res_imp <- robCompositions::impKNNa(as.data.frame(na_mat_sim), k = k, method = "knn")
          imp_vals <- as.matrix(res_imp$xImp)[miss_idx]
          
          # Calculate NRMSE
          rmse <- sqrt(mean((imp_vals - actual_vals)^2))
          errors[i] <- rmse / (max(actual_vals) - min(actual_vals))
        }, error = function(e) {
          # Log error if optimization fails (to avoid silent 0.0000 bugs)
          cat(sprintf("      [Error in Sim] k=%d rep=%d: %s\n", k, i, e$message))
          errors[i] <- NA 
        })
      }
      
      avg_err <- mean(errors, na.rm=TRUE)
      results_k <- rbind(results_k, data.frame(k = k, NRMSE = avg_err))
      cat(sprintf("     -> k=%2d | Avg NRMSE: %.4f\n", k, avg_err))
    }
    
    if(!all(is.na(results_k$NRMSE)) && !all(is.nan(results_k$NRMSE))) {
      best_k <- results_k$k[which.min(results_k$NRMSE)]
      k_final <- best_k
      cat(sprintf("   -> [RESULT] Optimal k detected: %d\n", k_final))
      
      # Save Plot
      p_k <- ggplot(results_k, aes(x=k, y=NRMSE)) +
        geom_line(color="#2E8B57") + geom_point(size=3) +
        geom_vline(xintercept=best_k, linetype="dashed", color="red") +
        theme_bw() + labs(title="kNN (Aitchison) Parameter Optimization", y="Imputation Error (NRMSE)")
      ggsave(file.path(CONFIG$qc_dir, "kNN_Optimization.png"), p_k, width=6, height=4)
    } else {
      cat("   [WARNING] Optimization failed (all NRMSE are NaN). Reverting to default k=3.\n")
    }
  }
}

# --- E. Final Imputation & Transformation ---
cat(sprintf("[5] Final Imputation & Transformation (Using k=%d)...\n", k_final))

# 1. kNN Imputation (on the Zero-Free dataset from Step C)
if(any(is.na(df_ready_for_knn))) {
  cat(sprintf("   -> Running Aitchison kNN imputation (k=%d)...\n", k_final))
  # Data already free of zeros, safe for impKNNa. Removed 'print=FALSE'.
  res_imp <- robCompositions::impKNNa(df_ready_for_knn, k = k_final, method = "knn")
  df_clean <- res_imp$xImp
} else {
  df_clean <- df_ready_for_knn
}

# 2. CLR Transformation
cat("   -> Applying CLR transformation...\n")
gmean <- function(x) exp(mean(log(x)))
df_clr_vals <- t(apply(df_clean, 1, function(x) log(x / gmean(x))))
df_clr_vals <- as.data.frame(df_clr_vals)
colnames(df_clr_vals) <- colnames(df_clean)

# Combine back with metadata
df_clr_final <- cbind(full_data[, meta_cols], df_clr_vals)

# --- F. Save Master Object ---
cat("[6] Saving Master Data Object...\n")
master_data <- list(
  metadata = full_data[, meta_cols],
  raw_filtered = df_final_raw,      
  clr_transformed = df_clr_final,   
  markers = marker_cols,
  parameters = CONFIG,
  optimal_k = k_final 
)

saveRDS(master_data, file.path(CONFIG$output_dir, "clean_data.rds"))
cat("   -> Done! Saved to 'processed_data/clean_data.rds'\n")
