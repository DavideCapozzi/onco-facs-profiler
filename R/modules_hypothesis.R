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

#' @title Run Pairwise Beta-Dispersion Test
#' @description 
#' Performs pairwise comparisons of multivariate dispersion between groups.
#' Uses vegan::betadisper followed by permutest.
#' 
#' @param data_input Dataframe with metadata and numeric data.
#' @param group_col String. Name of the grouping column.
#' @param metadata_cols Vector of columns to exclude.
#' @param n_perm Integer. Number of permutations.
#' @param min_n Integer. Minimum sample size per group required to run the test.
#' @return A dataframe of pairwise results.
run_pairwise_betadisper <- function(data_input, group_col = "Group", 
                                    metadata_cols = c("Patient_ID"), 
                                    n_perm = 999, min_n = 5) {
  requireNamespace("vegan", quietly = TRUE)
  requireNamespace("utils", quietly = TRUE)
  
  # Setup Data
  cols_to_exclude <- c(metadata_cols, group_col)
  numeric_data <- data_input[, !names(data_input) %in% cols_to_exclude, drop = FALSE]
  
  # Check numeric integrity
  if (!all(sapply(numeric_data, is.numeric))) {
    warning("Non-numeric columns found in data matrix. Attempting conversion.")
  }
  mat <- as.matrix(numeric_data)
  
  groups <- as.factor(data_input[[group_col]])
  levels_vec <- levels(groups)
  
  if (length(levels_vec) < 2) return(NULL)
  
  # Generate all unique pairs
  pairs <- utils::combn(levels_vec, 2, simplify = FALSE)
  results <- data.frame()
  
  message(sprintf("   [Stats] Running Pairwise Beta-Dispersion on %d pairs (Min N=%d)...", 
                  length(pairs), min_n))
  
  for (pair in pairs) {
    g1 <- pair[1]
    g2 <- pair[2]
    
    # Subset data for the current pair
    idx <- groups %in% c(g1, g2)
    sub_mat <- mat[idx, , drop = FALSE]
    sub_grps <- droplevels(groups[idx])
    
    # Check sample sizes
    n_g1 <- sum(sub_grps == g1)
    n_g2 <- sum(sub_grps == g2)
    
    # Skip if insufficient data
    if (n_g1 < min_n || n_g2 < min_n) {
      next
    }
    
    tryCatch({
      # 1. Calculate Distance Matrix for the pair
      # Note: Recalculating distance on subset ensures the geometric space is specific to the pair
      dist_mat <- vegan::vegdist(sub_mat, method = "euclidean")
      
      # 2. Run Betadisper
      mod <- vegan::betadisper(dist_mat, sub_grps, type = "median")
      
      # 3. Permutation Test
      perm_res <- vegan::permutest(mod, permutations = n_perm, pairwise = FALSE)
      
      # Extract Stats
      p_val <- perm_res$tab["Groups", "Pr(>F)"]
      f_val <- perm_res$tab["Groups", "F"]
      
      # Get Average Distances (to understand who has more dispersion)
      dist_g1 <- mod$group.distances[g1]
      dist_g2 <- mod$group.distances[g2]
      
      results <- rbind(results, data.frame(
        Group1 = g1, 
        Group2 = g2,
        Avg_Dist_G1 = round(dist_g1, 3),
        Avg_Dist_G2 = round(dist_g2, 3),
        F_Score = round(f_val, 3),
        P_Value = p_val,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      warning(paste("      [Fail] Pairwise Disp-Test failed for", g1, "vs", g2))
    })
  }
  
  # Apply FDR Correction
  if (nrow(results) > 0) {
    results$FDR <- p.adjust(results$P_Value, method = "BH")
    results <- results[order(results$P_Value), ]
  }
  
  return(results)
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

#' @title Run Pairwise PERMANOVA
#' @description 
#' Performs pairwise comparisons between all levels of a grouping factor using adonis2.
#' Automatically filters comparisons where sample size is insufficient.
#' Includes Benjamini-Hochberg FDR correction.
#' 
#' @param data_input Dataframe with metadata and numeric data.
#' @param group_col String. Name of the grouping column.
#' @param metadata_cols Vector of columns to exclude.
#' @param n_perm Integer. Number of permutations.
#' @param min_n Integer. Minimum sample size per group required to run the test.
#' @return A dataframe of pairwise results.
run_pairwise_permanova <- function(data_input, group_col = "Group", 
                                   metadata_cols = c("Patient_ID"), 
                                   n_perm = 999, min_n = 5) {
  requireNamespace("vegan", quietly = TRUE)
  requireNamespace("utils", quietly = TRUE)
  
  # Setup Data
  cols_to_exclude <- c(metadata_cols, group_col)
  numeric_data <- data_input[, !names(data_input) %in% cols_to_exclude, drop = FALSE]
  mat <- as.matrix(numeric_data)
  groups <- as.factor(data_input[[group_col]])
  levels_vec <- levels(groups)
  
  if (length(levels_vec) < 2) return(NULL)
  
  # Generate all unique pairs
  pairs <- utils::combn(levels_vec, 2, simplify = FALSE)
  results <- data.frame()
  
  message(sprintf("   [Stats] Running Pairwise PERMANOVA on %d pairs (Min N=%d)...", 
                  length(pairs), min_n))
  
  for (pair in pairs) {
    g1 <- pair[1]
    g2 <- pair[2]
    
    # Subset data for the current pair
    idx <- groups %in% c(g1, g2)
    sub_mat <- mat[idx, , drop = FALSE]
    sub_grps <- droplevels(groups[idx])
    
    # Check sample sizes
    n_g1 <- sum(sub_grps == g1)
    n_g2 <- sum(sub_grps == g2)
    
    # Skip if insufficient data (Robustness check for small groups like HNSCC_LS)
    if (n_g1 < min_n || n_g2 < min_n) {
      warning(sprintf("      [Skip] %s vs %s (N=%d/%d < %d)", g1, g2, n_g1, n_g2, min_n))
      next
    }
    
    tryCatch({
      # Run adonis2 on the subset
      res <- vegan::adonis2(sub_mat ~ sub_grps, method = "euclidean", permutations = n_perm)
      
      p_val <- res$`Pr(>F)`[1]
      r2    <- res$R2[1]
      f_val <- res$F[1]
      
      results <- rbind(results, data.frame(
        Group1 = g1, 
        Group2 = g2,
        N1 = n_g1,
        N2 = n_g2,
        R2_Percent = round(r2 * 100, 2),
        F_Model = round(f_val, 2),
        P_Value = p_val,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      warning(paste("      [Fail] Pairwise PERMANOVA failed for", g1, "vs", g2))
    })
  }
  
  # Apply FDR Correction (Benjamini-Hochberg) across all tests performed
  if (nrow(results) > 0) {
    results$FDR <- p.adjust(results$P_Value, method = "BH")
    # Reorder by significance
    results <- results[order(results$P_Value), ]
  }
  
  return(results)
}
