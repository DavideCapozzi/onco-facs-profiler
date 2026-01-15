# R/modules_hypothesis.R
# ==============================================================================
# HYPOTHESIS TESTING MODULE
# Description: MANOVA, PERMANOVA, and Assumption Checks for CoDa.
# Dependencies: vegan, MASS, dplyr
# ==============================================================================

library(dplyr)

#' @title Test Multivariate Homogeneity of Group Dispersions (Betadisper)
#' @description 
#' Implements Anderson's PERMDISP2 procedure for the analysis of multivariate homogeneity of group dispersions.
#' Ideally, this should be non-significant to validate PERMANOVA results.
#' 
#' @param data_input A dataframe containing markers/coordinates AND metadata.
#' @param group_col Name of the grouping column.
#' @param metadata_cols Vector of columns to exclude from the distance calculation.
#' @param n_perm Number of permutations for the significance test.
#' @return A list containing the test summary table and group average distances.
test_coda_dispersion <- function(data_input, group_col = "Group", metadata_cols = c("Patient_ID"), n_perm = 999) {
  requireNamespace("vegan", quietly = TRUE)
  
  # Ensure Group exists
  if (!group_col %in% names(data_input)) stop(sprintf("Column '%s' not found in input.", group_col))
  
  # Prepare Numeric Matrix
  cols_to_exclude <- c(metadata_cols, group_col)
  numeric_data <- data_input[, !names(data_input) %in% cols_to_exclude, drop = FALSE]
  
  # Check numeric integrity
  if (!all(sapply(numeric_data, is.numeric))) {
    warning("Non-numeric columns found in data matrix. Attempting conversion.")
  }
  mat <- as.matrix(numeric_data)
  
  groups <- as.factor(data_input[[group_col]])
  
  message(sprintf("   [Stats] Running Beta-Dispersion Check on %d vars x %d samples...", 
                  ncol(mat), nrow(mat)))
  
  # 1. Calculate Euclidean Distance (Equivalent to Aitchison if CLR input)
  dist_mat <- vegan::vegdist(mat, method = "euclidean")
  
  # 2. Calculate Multivariate Dispersions
  mod_disp <- vegan::betadisper(dist_mat, groups, type = "median")
  
  # 3. Permutation Test for Significance
  res_perm <- vegan::permutest(mod_disp, permutations = n_perm, pairwise = FALSE)
  
  # 4. Format Output for Export
  # Extract generic ANOVA table from the permutest result
  res_table <- as.data.frame(res_perm$tab)
  
  # Extract average distance to median per group (Interpretation: Higher = More heterogeneous)
  group_distances <- data.frame(
    Group = names(mod_disp$group.distances),
    Avg_Distance_to_Median = as.numeric(mod_disp$group.distances)
  )
  
  return(list(
    anova_table = res_table,
    group_distances = group_distances
  ))
}

#' @title Run PERMANOVA 
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
