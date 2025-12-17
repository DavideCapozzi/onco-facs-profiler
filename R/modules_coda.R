# R/modules_coda.R
# ==============================================================================
# COMPOSITIONAL DATA ANALYSIS MODULE
# Description: Pure functions for Zero Replacement (CZM), Imputation, and CLR/ILR.
# Dependencies: zCompositions, dplyr, compositions
# ==============================================================================

library(zCompositions)
library(dplyr)
library(compositions)

#' @title Zero Replacement via Count Zero Multiplicative (CZM) method
#' @description 
#' Replaces zeros in a compositional matrix using Bayesian-Multiplicative replacement.
#' Handles the "NA conflict" by temporarily masking NAs.
#' 
#' @param mat A numeric matrix (samples x markers).
#' @return A matrix with zeros replaced by small estimated values.
coda_replace_zeros <- function(mat) {
  
  # Validation
  if (!is.matrix(mat)) stop("Input must be a matrix.")
  
  if (!any(mat == 0, na.rm = TRUE)) {
    message("   [CoDa] No zeros detected. Skipping zero replacement.")
    return(mat)
  }
  
  message("   [CoDa] Zeros detected. Applying zCompositions::cmultRepl (CZM)...")
  
  # --- NA HANDLING STRATEGY ---
  # cmultRepl fails with NAs. We perform a "Mask & Restore" strategy.
  
  na_locs <- is.na(mat)
  has_na <- any(na_locs)
  mat_temp <- mat
  
  # 1. Temporary Fill (if NAs exist) to allow cmultRepl to run on Zeros
  if (has_na) {
    gmean_safe <- function(x) {
      vals <- x[!is.na(x) & x > 0]
      if(length(vals) == 0) return(1e-6)
      return(exp(mean(log(vals))))
    }
    
    for(j in 1:ncol(mat_temp)) {
      if(any(na_locs[,j])) {
        fill_val <- gmean_safe(mat_temp[,j])
        mat_temp[na_locs[,j], j] <- fill_val
      }
    }
  }
  
  # 2. Apply cmultRepl
  tryCatch({
    clean_mat <- zCompositions::cmultRepl(
      X = mat_temp,
      label = 0,
      method = "CZM", 
      output = "p-counts",
      suppress.print = TRUE 
    )
    
    # 3. Restore NAs (Actual imputation happens in next step)
    if (has_na) {
      clean_mat[na_locs] <- NA
    }
    
    return(clean_mat)
    
  }, error = function(e) {
    stop(paste("Zero replacement via CZM failed:", e$message))
  })
}

#' @title Imputation for Remaining NAs
#' @description 
#' Handles remaining NAs with Multiplicative Replacement (preserving ratios).
impute_nas_simple <- function(mat) {
  
  if (!any(is.na(mat))) return(mat)
  
  message("   [CoDa] Imputing remaining NAs using Multiplicative Replacement...")
  
  tryCatch({
    # imp.missing = TRUE imputes NAs specifically
    result <- zCompositions::multRepl(X = mat, label = NA, imp.missing = TRUE)
    return(as.matrix(result))
    
  }, error = function(e) {
    warning(paste("multRepl failed, using geometric mean fallback. Error:", e$message))
    
    mat_fixed <- mat
    for (j in 1:ncol(mat_fixed)) {
      na_idx <- is.na(mat_fixed[,j])
      if (any(na_idx)) {
        valid <- mat_fixed[!na_idx, j]
        if (length(valid) > 0) {
          gmean <- exp(mean(log(valid[valid > 0])))
          mat_fixed[na_idx, j] <- gmean
        }
      }
    }
    return(mat_fixed)
  })
}

#' @title Centered Log-Ratio Transformation (CLR)
#' @description Projects compositional data from simplex to real space (Euclidean).
coda_transform_clr <- function(mat) {
  
  if (any(mat <= 0, na.rm = TRUE)) stop("CLR input must be strictly positive (no zeros).")
  
  message("   [CoDa] Applying CLR transformation...")
  
  gmean <- function(x) exp(mean(log(x)))
  clr_mat <- t(apply(mat, 1, function(x) log(x / gmean(x))))
  
  return(clr_mat)
}

#' @title Isometric Log-Ratio Transformation (ILR)
#' @description Transforms data into D-1 orthogonal coordinates for MANOVA.
transform_ilr <- function(mat) {
  
  if (any(mat <= 0, na.rm = TRUE)) stop("ILR input must be strictly positive.")
  
  message("   [CoDa] Applying ILR transformation (D-1 coordinates)...")
  
  ilr_obj <- compositions::ilr(compositions::acomp(mat))
  ilr_mat <- as.matrix(ilr_obj)
  rownames(ilr_mat) <- rownames(mat)
  
  return(ilr_mat)
}