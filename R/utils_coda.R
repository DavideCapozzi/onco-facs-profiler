# R/utils_coda.R
# ==============================================================================
# COMPOSITIONAL DATA UTILITIES
# Description: Handling Zeros, NAs, and CLR Transformations.
# ==============================================================================

library(robCompositions)
library(zCompositions)
library(dplyr)

#' Geometric Mean
#' @param x Numeric vector
#' @return Geometric mean (scalar)
gmean <- function(x) {
  exp(mean(log(x), na.rm = TRUE))
}

#' Apply Centered Log-Ratio (CLR) Transformation
#' @param data_matrix Numeric matrix (n samples x p parts) - No NAs allowed
#' @return CLR transformed matrix
apply_clr <- function(data_matrix) {
  # Validate input
  if (any(is.na(data_matrix))) stop("CLR Error: Matrix contains NAs.")
  if (any(data_matrix <= 0)) stop("CLR Error: Matrix contains values <= 0.")
  
  # Apply CLR row-wise
  clr_mat <- t(apply(data_matrix, 1, function(x) {
    log(x / gmean(x))
  }))
  
  return(as.matrix(clr_mat))
}

#' Legacy Imputation Strategy (The "Median Hack")
#' @description 
#'   1. Detects NAs.
#'   2. Fills NAs with Column Medians (Temporary).
#'   3. Runs cmultRepl (CZM) to fix zeros (using the full matrix).
#'   4. Restores original NAs.
#'   5. Runs impKNNa (Aitchison kNN) to impute NAs properly.
#' 
#' @param data_df Numeric data.frame containing markers
#' @param k Integer, neighbors for kNN
#' @return Matrix of fully imputed data
run_complex_imputation <- function(data_df, k = 3) {
  
  # Convert to matrix for checks
  mat_vals <- as.matrix(data_df)
  na_map <- is.na(mat_vals)
  has_zeros <- any(mat_vals == 0, na.rm = TRUE)
  
  message(">>> [CoDa] Starting Imputation Pipeline...")
  
  # --- STEP A: Handle Zeros (if any) ---
  if (has_zeros) {
    message("    -> Zeros detected. Applying CZM with Median-Injection (Legacy Fix).")
    
    # 1. Temporary Median Fill
    temp_filled <- data_df
    for(j in 1:ncol(temp_filled)) {
      col_med <- median(temp_filled[[j]], na.rm=TRUE)
      # Fallback for completely empty cols (should be caught by QC, but safety first)
      if(is.na(col_med)) col_med <- 0.001 
      temp_filled[is.na(temp_filled[[j]]), j] <- col_med
    }
    
    # 2. Apply CZM
    # output="p-counts" returns the corrected counts/proportions
    clean_zeros <- zCompositions::cmultRepl(temp_filled, method="CZM", 
                                            output="p-counts", suppress.print=TRUE)
    
    # 3. Restore NAs
    # We now have a matrix with NO zeros, but median-filled NAs.
    # We must put the NAs back to let kNN handle them.
    clean_zeros[na_map] <- NA
    data_ready_for_knn <- clean_zeros
    
  } else {
    message("    -> No zeros detected. Skipping CZM.")
    data_ready_for_knn <- data_df
  }
  
  # --- STEP B: Handle NAs (kNN) ---
  # Check if NAs remain
  if (any(is.na(data_ready_for_knn))) {
    message(sprintf("    -> Running Aitchison kNN Imputation (k=%d)...", k))
    
    # impKNNa expects a data.frame
    res_imp <- robCompositions::impKNNa(as.data.frame(data_ready_for_knn), 
                                        k = k, method = "knn")
    final_mat <- as.matrix(res_imp$xImp)
    
  } else {
    final_mat <- as.matrix(data_ready_for_knn)
  }
  
  return(final_mat)
}
