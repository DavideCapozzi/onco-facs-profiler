# R/modules_qc.R
# ==============================================================================
# QUALITY CONTROL MODULE
# Description: pure functions for data filtering (Variance, Missingness).
# ==============================================================================

#' @title Run Quality Control Pipeline
#' @description 
#' Filters the data matrix based on variance and missingness (NA) thresholds.
#' Logs dropped items to the console and compiles a summary object for reporting.
#' 
#' @param mat_raw Raw numeric matrix (Samples x Markers).
#' @param qc_config List containing 'max_na_row_pct' and 'max_na_col_pct'.
#' @param dropped_apriori Dataframe of markers dropped by whitelist/blacklist (optional).
#' @return A list containing:
#'   - data: The filtered matrix.
#'   - report: A list compatible with save_qc_report().
run_qc_pipeline <- function(mat_raw, qc_config, dropped_apriori = NULL) {
  
  message("[QC] Running Quality Control...")
  
  # Initialize QC Summary Object
  qc_summary <- list(
    n_row_init = nrow(mat_raw),
    n_col_init = ncol(mat_raw),
    dropped_markers_apriori = dropped_apriori,
    n_col_zerovar = 0,
    dropped_rows_detail = data.frame(Patient_ID = character(), NA_Percent = numeric()),
    dropped_cols_detail = data.frame(Marker = character(), NA_Percent = numeric()),
    n_row_dropped = 0,
    n_col_dropped = 0
  )
  
  curr_mat <- mat_raw
  
  # 1. Remove Constant Columns (Zero Variance)
  col_vars <- apply(curr_mat, 2, var, na.rm = TRUE)
  # Check for 0 or NA (if all NAs, var is NA)
  const_cols <- which(col_vars == 0 | is.na(col_vars))
  
  if (length(const_cols) > 0) {
    message(sprintf("   [QC] Dropping %d markers with zero variance.", length(const_cols)))
    qc_summary$n_col_zerovar <- length(const_cols)
    curr_mat <- curr_mat[, -const_cols, drop = FALSE]
  }
  
  # 2. Filter by Missingness (Rows - Patients)
  row_na_freq <- rowMeans(is.na(curr_mat))
  drop_row_idx <- which(row_na_freq > qc_config$max_na_row_pct)
  
  if (length(drop_row_idx) > 0) {
    message(sprintf("   [QC] Dropping %d patients with >%.0f%% missingness:", 
                    length(drop_row_idx), qc_config$max_na_row_pct * 100))
    
    dropped_pats <- data.frame(
      Patient_ID = rownames(curr_mat)[drop_row_idx],
      NA_Percent = round(row_na_freq[drop_row_idx] * 100, 2)
    )
    
    # Log to console
    print_df <- paste0("      - ", dropped_pats$Patient_ID, " (", dropped_pats$NA_Percent, "%)")
    cat(paste(print_df, collapse = "\n"), "\n")
    
    qc_summary$dropped_rows_detail <- dropped_pats
    qc_summary$n_row_dropped <- length(drop_row_idx)
    
    curr_mat <- curr_mat[-drop_row_idx, , drop = FALSE]
  }
  
  # 3. Filter by Missingness (Cols - Markers)
  col_na_freq <- colMeans(is.na(curr_mat))
  drop_col_idx <- which(col_na_freq > qc_config$max_na_col_pct)
  
  if (length(drop_col_idx) > 0) {
    message(sprintf("   [QC] Dropping %d markers with >%.0f%% missingness:", 
                    length(drop_col_idx), qc_config$max_na_col_pct * 100))
    
    dropped_markers <- data.frame(
      Marker = colnames(curr_mat)[drop_col_idx],
      NA_Percent = round(col_na_freq[drop_col_idx] * 100, 2)
    )
    
    # Log to console
    print_df <- paste0("      - ", dropped_markers$Marker, " (", dropped_markers$NA_Percent, "%)")
    cat(paste(print_df, collapse = "\n"), "\n")
    
    qc_summary$dropped_cols_detail <- dropped_markers
    qc_summary$n_col_dropped <- length(drop_col_idx)
    
    curr_mat <- curr_mat[, -drop_col_idx, drop = FALSE]
  }
  
  # Final Stats
  qc_summary$n_row_final <- nrow(curr_mat)
  qc_summary$n_col_final <- ncol(curr_mat)
  
  message(sprintf("   [QC] Final Dimensions: %d Samples x %d Markers", nrow(curr_mat), ncol(curr_mat)))
  
  if (nrow(curr_mat) < 10) warning("[QC] WARNING: Very low sample size. Analysis may be unstable.")
  
  return(list(
    data = curr_mat,
    report = qc_summary
  ))
}