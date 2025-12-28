# R/modules_hypothesis.R
# ==============================================================================
# HYPOTHESIS TESTING MODULE
# Description: MANOVA, PERMANOVA, and Assumption Checks for CoDa.
# Dependencies: vegan, MASS, dplyr
# ==============================================================================

library(dplyr)

#' @title Run PERMANOVA (Robust)
#' @description Performs Permutational Multivariate Analysis of Variance (adonis2).
#' @param data_input A dataframe containing markers/coordinates AND metadata.
#' @param group_col Name of the grouping column.
#' @param metadata_cols Vector of columns to exclude from the distance calculation.
#' @param n_perm Number of permutations.
test_coda_permanova <- function(data_input, group_col = "Group", metadata_cols = c("Patient_ID"), n_perm = 999) {
  requireNamespace("vegan", quietly = TRUE)
  
  # Ensure Group exists
  if (!group_col %in% names(data_input)) stop(sprintf("Column '%s' not found in input.", group_col))
  
  # Prepare Numeric Matrix for Distance
  # Exclude Group and other metadata
  cols_to_exclude <- c(metadata_cols, group_col)
  numeric_data <- data_input[, !names(data_input) %in% cols_to_exclude, drop = FALSE]
  
  # Check if numeric
  if (!all(sapply(numeric_data, is.numeric))) {
    warning("Non-numeric columns found in data matrix. Attempting conversion.")
  }
  mat <- as.matrix(numeric_data)
  
  groups <- data_input[[group_col]]
  
  # Clean NA if present (adonis2 fails with NA)
  if (any(is.na(mat))) {
    stop("Input data for PERMANOVA contains NAs. Please check imputation.")
  }
  
  message(sprintf("   [Stats] Running PERMANOVA on %d vars x %d samples (%d perms)...", 
                  ncol(mat), nrow(mat), n_perm))
  
  # Run adonis2 (Euclidean distance is appropriate for CLR/ILR/Z-score)
  res_adonis <- vegan::adonis2(mat ~ groups, method = "euclidean", permutations = n_perm)
  
  return(res_adonis)
}
