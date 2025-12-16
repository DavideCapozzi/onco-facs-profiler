# R/modules_coda.R
# ==============================================================================
# COMPOSITIONAL DATA ANALYSIS MODULE
# Description: Pure functions for Zero Replacement (CZM), Imputation, and CLR.
# Dependencies: zCompositions, dplyr
# ==============================================================================

library(zCompositions)
library(dplyr)
library(compositions)

#' @title Zero Replacement via Count Zero Multiplicative (CZM) method
#' @description 
#' Replaces zeros in a compositional matrix using Bayesian-Multiplicative replacement (CZM).
#' This is the preferred method for rounded zeros in flow cytometry data.
#' Handles the "NA conflict" by temporarily masking NAs to allow cmultRepl to run.
#' 
#' @param mat A numeric matrix (samples x markers).
#' @return A matrix with zeros replaced by small estimated values (pseudo-counts), and NAs preserved.
replace_zeros_czm <- function(mat) {
  
  # Validation: Ensure input is a matrix
  if (!is.matrix(mat)) stop("Input must be a matrix.")
  
  # Check if zeros exist
  if (!any(mat == 0, na.rm = TRUE)) {
    message("   [CoDa] No zeros detected. Skipping zero replacement.")
    return(mat)
  }
  
  message("   [CoDa] Zeros detected. Applying zCompositions::cmultRepl (CZM)...")
  
  # --- NA HANDLING STRATEGY ---
  # cmultRepl fails if it encounters NAs. We perform a "Mask & Restore" strategy.
  
  # 1. Identify NA locations
  na_locs <- is.na(mat)
  has_na <- any(na_locs)
  
  # 2. Temporary Fill (if NAs exist)
  # We fill NAs with the column geometric mean (neutral value in CoDa) 
  # just to let cmultRepl focus on the Zeros without crashing.
  mat_temp <- mat
  
  if (has_na) {
    # Helper for geometric mean ignoring NA and <= 0
    gmean_safe <- function(x) {
      vals <- x[!is.na(x) & x > 0]
      if(length(vals) == 0) return(1e-6) # Fallback if col is empty/zeros
      return(exp(mean(log(vals))))
    }
    
    # Fill NAs column-wise
    for(j in 1:ncol(mat_temp)) {
      if(any(na_locs[,j])) {
        fill_val <- gmean_safe(mat_temp[,j])
        mat_temp[na_locs[,j], j] <- fill_val
      }
    }
  }
  
  # 3. Apply cmultRepl on the "complete" (masked) matrix
  # output = "p-counts" preserves the scale of the data (pseudo-counts)
  tryCatch({
    clean_mat <- zCompositions::cmultRepl(
      X = mat_temp,
      label = 0,
      method = "CZM", 
      output = "p-counts",
      suppress.print = TRUE 
    )
    
    # 4. Restore NAs (The actual Imputation will be done in the next step)
    if (has_na) {
      clean_mat[na_locs] <- NA
    }
    
    return(clean_mat)
    
  }, error = function(e) {
    stop(paste("Zero replacement via CZM failed:", e$message))
  })
}

#' @title Simple Imputation for Remaining NAs
#' @description 
#' After zero replacement, handles any remaining NAs with a stable, 
#' non-stochastic Multiplicative Replacement method.
#' 
#' @param mat A numeric matrix (no zeros, but may contain NAs).
#' @return A complete matrix with no NAs and strictly positive values.
impute_nas_simple <- function(mat) {
  
  # If no NAs, return immediately
  if (!any(is.na(mat))) return(mat)
  
  message("   [CoDa] Imputing remaining NAs using Multiplicative Replacement...")
  
  # Option 1: Use multRepl with imp.missing=TRUE
  # This imputes missing values preserving the log-ratio structure
  tryCatch({
    # FIX: removed 'suppress.print' which caused the error
    result <- zCompositions::multRepl(
      X = mat, 
      label = NA, 
      imp.missing = TRUE
    )
    return(as.matrix(result))
    
  }, error = function(e) {
    # Fallback: column geometric mean imputation
    warning(paste("multRepl failed for NAs, using geometric mean fallback. Error:", e$message))
    
    mat_fixed <- mat
    for (j in 1:ncol(mat_fixed)) {
      na_idx <- is.na(mat_fixed[,j])
      if (any(na_idx)) {
        valid <- mat_fixed[!na_idx, j]
        # Calculate geometric mean of valid parts
        if (length(valid) > 0) {
          gmean <- exp(mean(log(valid[valid > 0])))
          mat_fixed[na_idx, j] <- gmean
        } else {
          stop(sprintf("Column %d contains only NAs or Zeros. Cannot impute.", j))
        }
      }
    }
    return(mat_fixed)
  })
}

#' @title Centered Log-Ratio Transformation (CLR)
#' @description 
#' Projects the compositional data from the simplex to real space (Euclidean).
#' 
#' @param mat A numeric matrix (positive values only).
#' @return A CLR-transformed matrix.
transform_clr <- function(mat) {
  
  if (any(mat <= 0, na.rm = TRUE)) stop("CLR input must be strictly positive (no zeros).")
  if (any(is.na(mat))) stop("CLR input must not contain NAs.")
  
  message("   [CoDa] Applying CLR transformation...")
  
  # Define Geometric Mean function
  gmean <- function(x) exp(mean(log(x)))
  
  # Apply CLR row-wise
  # formula: clr(x) = [ln(x1/g(x)), ..., ln(xp/g(x))]
  clr_mat <- t(apply(mat, 1, function(x) log(x / gmean(x))))
  
  return(clr_mat)
}

#' @title Isometric Log-Ratio Transformation (ILR)
#' @description 
#' Transforms data into D-1 orthogonal coordinates suitable for standard 
#' multivariate statistics (MANOVA, Regression).
#' Note: ILR variables do not map 1:1 to markers; they represent balances.
#' 
#' @param mat A numeric matrix (positive values only).
#' @return A matrix of ILR coordinates (n x D-1).
transform_ilr <- function(mat) {
  
  if (any(mat <= 0, na.rm = TRUE)) stop("ILR input must be strictly positive.")
  
  message("   [CoDa] Applying ILR transformation (D-1 coordinates)...")
  
  # Use the compositions package for robust basis construction
  # We use the default basis (balanced) as MANOVA is invariant to rotation
  ilr_obj <- compositions::ilr(compositions::acomp(mat))
  
  # Convert back to standard matrix
  ilr_mat <- as.matrix(ilr_obj)
  rownames(ilr_mat) <- rownames(mat)
  
  return(ilr_mat)
}