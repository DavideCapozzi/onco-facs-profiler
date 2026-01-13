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

#' @title Compute Local ILR for Compositional Subgroups
#' @description 
#' Iterates through defined groups in config. If a group has >1 marker, 
#' it treats it as a sub-composition and calculates Isometric Log-Ratio coordinates.
#' 
#' @param mat_imputed A matrix of imputed counts/concentrations (strictly positive).
#' @param hybrid_groups A list from config defining marker groups.
#' @return A named list of matrices, where each matrix contains ILR coordinates for a group.
coda_compute_local_ilr <- function(mat_imputed, hybrid_groups) {
  
  if (any(mat_imputed <= 0, na.rm = TRUE)) stop("ILR input must be strictly positive.")
  
  ilr_results <- list()
  
  message("   [CoDa] Computing Local ILR balances for sub-compositions...")
  
  for (group_name in names(hybrid_groups)) {
    target_mks <- hybrid_groups[[group_name]]
    present_mks <- intersect(target_mks, colnames(mat_imputed))
    
    # Check 1: Markers must exist
    if (length(present_mks) == 0) next
    
    # Check 2: Must have at least 2 markers for ILR (D parts -> D-1 balances)
    if (length(present_mks) < 2) {
      # message(sprintf("      -> Skip '%s': < 2 markers (cannot calculate balances).", group_name))
      next
    }
    
    # Subset and Closure (acomp does closure automatically)
    mat_sub <- mat_imputed[, present_mks, drop = FALSE]
    
    # Calculate ILR (default basis)
    # Note: Column names will be V1, V2... Interpretability requires looking at the basis.
    ilr_obj <- compositions::ilr(compositions::acomp(mat_sub))
    ilr_mat <- as.matrix(ilr_obj)
    
    # Preserve Rownames
    rownames(ilr_mat) <- rownames(mat_imputed)
    
    ilr_results[[group_name]] <- ilr_mat
    message(sprintf("      -> '%s': Computed %d ILR balances.", group_name, ncol(ilr_mat)))
  }
  
  return(ilr_results)
}

#' @title Robust Logit Transformation
#' @description 
#' Transforms percentages (0-100) or proportions (0-1) into Logit space (Log-Odds).
#' Crucially, it handles boundary conditions (0% and 100%) by clamping values 
#' to a safety range [epsilon, 1-epsilon] to avoid infinite results.
#' 
#' @param mat A numeric matrix or vector.
#' @param epsilon A small numeric constant to prevent division by zero or log(0). 
#'        Default is 1e-6 (0.0001%).
#' @return A matrix of the same dimensions, transformed to real space (-Inf, +Inf).
coda_transform_logit <- function(mat, epsilon = 1e-6) {
  
  # 1. Validation and Auto-Scaling
  # Check if data looks like percentages (0-100) or proportions (0-1).
  # If max value > 1, assume percentages and divide by 100.
  max_val <- max(mat, na.rm = TRUE)
  
  p <- if (max_val > 1.0) {
    message("   [CoDa] Detected percentages (0-100). Scaling to proportions (0-1)...")
    mat / 100
  } else {
    mat
  }
  
  # 2. Boundary Protection (Clamping)
  # Logit is undefined for p=0 (gives -Inf) and p=1 (gives +Inf).
  # We clamp values within [epsilon, 1-epsilon].
  # This preserves the distribution shape without breaking downstream PCA/Networks.
  
  # Count boundary hits for logging
  n_low  <- sum(p < epsilon, na.rm = TRUE)
  n_high <- sum(p > (1 - epsilon), na.rm = TRUE)
  
  if (n_low > 0 || n_high > 0) {
    message(sprintf("   [CoDa] Clamped %d boundary values (0 or 1) to epsilon range.", n_low + n_high))
  }
  
  p[p < epsilon] <- epsilon
  p[p > (1 - epsilon)] <- 1 - epsilon
  
  # 3. Logit Transformation
  # Formula: log( p / (1 - p) )
  logit_mat <- log(p / (1 - p))
  
  return(logit_mat)
}

#' @title Perform Hybrid Transformation Pipeline
#' @description 
#' Orchestrates the complete transformation logic:
#' 1. Splits markers into Compositional (Track A) and Functional (Track B) based on config.
#' 2. Track A: Zero replacement (CZM) -> CLR transformation.
#' 3. Track B: LOD Zero handling -> Logit transformation -> BPCA imputation.
#' 4. Merges tracks and performs Z-score standardization.
#' 5. Computes ILR balances for compositional groups.
#' 
#' @param mat_raw Raw numeric matrix (Samples x Markers). Post-QC.
#' @param config The configuration list containing 'hybrid_groups'.
#' @return A list containing:
#'   - hybrid_data_raw: Merged matrix (Transformed/Imputed but NOT scaled).
#'   - hybrid_data_z: Merged matrix (Transformed/Imputed AND Z-scored).
#'   - ilr_balances: List of ILR coordinates for defined groups.
#'   - hybrid_markers: Vector of final marker names.
#' @export
perform_hybrid_transformation <- function(mat_raw, config) {
  
  message("\n[CoDa] Starting Hybrid Strategy (Compositional vs Functional)...")
  
  # 1. Split Markers
  all_current_markers <- colnames(mat_raw)
  comp_markers_list <- unique(unlist(config$hybrid_groups))
  comp_markers <- intersect(all_current_markers, comp_markers_list)
  func_markers <- setdiff(all_current_markers, comp_markers)
  
  message(sprintf("   [Strategy] Split: %d Compositional (CLR) vs %d Functional (Logit+BPCA)", 
                  length(comp_markers), length(func_markers)))
  
  mat_comp_trans <- NULL
  mat_func_trans <- NULL
  ilr_results <- list()
  
  # --- TRACK A: Compositional Markers (CZM -> CLR) ---
  if (length(comp_markers) > 0) {
    message("   [Track A] Processing Compositional Groups...")
    mat_comp_raw <- mat_raw[, comp_markers, drop = FALSE]
    
    # A1. Zero Replacement (CZM)
    mat_comp_clean <- coda_replace_zeros(mat_comp_raw)
    
    # A2. Handle remaining NAs (Multiplicative Replacement)
    if (any(is.na(mat_comp_clean))) {
      message("      -> Imputing remaining NAs in composition (multRepl)...")
      mat_comp_clean <- zCompositions::multRepl(mat_comp_clean, label = NA, imp.missing = TRUE)
    }
    
    # A3. CLR Transformation (Global)
    mat_comp_trans <- coda_transform_clr(mat_comp_clean)
    
    # A4. Compute ILR Balances (for Stats)
    # We use the cleaned counts (mat_comp_clean) before CLR
    ilr_results <- coda_compute_local_ilr(mat_comp_clean, config$hybrid_groups)
  }
  
  # --- TRACK B: Functional Markers (LOD -> Logit -> BPCA) ---
  if (length(func_markers) > 0) {
    message("   [Track B] Processing Functional Markers...")
    mat_func_raw <- mat_raw[, func_markers, drop = FALSE]
    
    # B1. Zero Handling (LOD Proxy: Min/2)
    mat_func_lod <- mat_func_raw
    for (j in 1:ncol(mat_func_lod)) {
      col_vals <- mat_func_lod[, j]
      zero_idx <- which(col_vals == 0 & !is.na(col_vals))
      
      if (length(zero_idx) > 0) {
        pos_vals <- col_vals[col_vals > 0 & !is.na(col_vals)]
        lod_proxy <- if(length(pos_vals) > 0) min(pos_vals)/2 else 1e-6
        mat_func_lod[zero_idx, j] <- lod_proxy
      }
    }
    
    # B2. Logit Transformation
    mat_func_logit <- coda_transform_logit(mat_func_lod)
    
    # B3. BPCA Imputation
    mat_func_trans <- impute_matrix_bpca(mat_func_logit, nPcs = 8)
  }
  
  # --- MERGE & NORMALIZE ---
  message("[Merge] Recombining and creating final matrices...")
  
  list_parts <- list()
  if (!is.null(mat_comp_trans)) list_parts[[1]] <- mat_comp_trans
  if (!is.null(mat_func_trans)) list_parts[[2]] <- mat_func_trans
  
  if (length(list_parts) == 0) stop("[FATAL] No data remaining after filtering.")
  
  mat_hybrid_raw <- do.call(cbind, list_parts)
  
  # Restore alphabetic column order
  mat_hybrid_raw <- mat_hybrid_raw[, sort(colnames(mat_hybrid_raw)), drop = FALSE]
  
  # Z-Score Standardization
  message("   [Norm] Applying Z-Score Standardization...")
  mat_hybrid_z <- scale(mat_hybrid_raw)
  attr(mat_hybrid_z, "scaled:center") <- NULL
  attr(mat_hybrid_z, "scaled:scale") <- NULL
  mat_hybrid_z <- as.matrix(mat_hybrid_z)
  
  return(list(
    hybrid_data_raw = mat_hybrid_raw,
    hybrid_data_z = mat_hybrid_z,
    ilr_balances = ilr_results,
    hybrid_markers = colnames(mat_hybrid_z)
  ))
}

#' @title Bayesian PCA Imputation (Wrapper)
#' @description 
#' Imputes missing values using Bayesian PCA (pcaMethods::bpca).
#' Ideal for datasets with high collinearity and N approx P.
#' Works best on continuous/transformed data (e.g., after Logit/CLR).
#' 
#' @title Bayesian PCA Imputation (Wrapper)
#' @description 
#' Imputes missing values using Bayesian PCA (pcaMethods::bpca).
#' Includes automatic estimation of optimal components (kEstimate) if nPcs is not fixed.
#' Ideal for datasets with high collinearity.
#' 
#' @param mat Numeric matrix (samples x markers) containing NAs.
#' @param nPcs Number of components to estimate. Can be an integer or "auto". 
#'        If "auto", it uses Cross-Validation to find the best K.
#' @return A numeric matrix with NAs filled.
#' @export
impute_matrix_bpca <- function(mat, nPcs = "auto") {
  
  if (!any(is.na(mat))) {
    message("   [Impute] No NAs found. Skipping BPCA.")
    return(mat)
  }
  
  # Check dimensions
  if (nrow(mat) < 5) stop("Too few samples for BPCA imputation (min 5).")
  
  # Requirement check
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
    stop("[FATAL] Package 'pcaMethods' (Bioconductor) is required for BPCA imputation.\n",
         "Please install it via: BiocManager::install('pcaMethods')")
  }
  
  # Logic to determine nPcs
  k_used <- 3 # Default fallback
  
  if (is.null(nPcs) || nPcs == "auto") {
    message("   [Impute] Estimating optimal K using Cross-Validation (kEstimate)...")
    tryCatch({
      # Test range from 2 up to min(15, dim-1)
      max_k <- min(15, ncol(mat) - 1, nrow(mat) - 1)
      if (max_k < 2) max_k <- 2
      
      k_est <- pcaMethods::kEstimate(mat, method = "bpca", evalPcs = 2:max_k, verbose = FALSE)
      k_used <- k_est$bestNPcs
      message(sprintf("      -> Estimated Optimal K: %d components", k_used))
      
    }, error = function(e) {
      warning(paste("   [Impute] kEstimate failed, falling back to default (3). Error:", e$message))
      k_used <- 3
    })
  } else {
    k_used <- as.numeric(nPcs)
  }
  
  message(sprintf("   [Impute] Running Bayesian PCA (nPcs=%d) on %d x %d matrix...", 
                  k_used, nrow(mat), ncol(mat)))
  
  tryCatch({
    # BPCA requires samples in rows, variables in columns (standard R)
    res_bpca <- pcaMethods::pca(mat, method = "bpca", nPcs = k_used, verbose = FALSE)
    mat_imputed <- pcaMethods::completeObs(res_bpca)
    
    # Validation: Check for extreme outliers generated by imputation
    range_orig <- range(mat, na.rm = TRUE)
    range_imp  <- range(mat_imputed)
    
    if (any(abs(range_imp) > 10 * max(abs(range_orig)))) {
      warning("[Impute] BPCA produced values far outside original range. Check data distribution.")
    }
    
    return(as.matrix(mat_imputed))
    
  }, error = function(e) {
    stop(paste("BPCA Imputation failed:", e$message))
  })
}