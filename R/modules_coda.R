# R/modules_coda.R
# ==============================================================================
# COMPOSITIONAL DATA ANALYSIS MODULE
# Description: Pure functions for Zero Replacement (CZM), Imputation, and CLR/ILR.
# Dependencies: zCompositions, dplyr, compositions, pcaMethods
# ==============================================================================

library(zCompositions)
library(dplyr)
library(compositions)

# --- HELPER FUNCTIONS ---

#' @title Calculate Geometric Mean (Safe)
#' @description Computes geometric mean handling positive values only.
#' @param x Numeric vector.
#' @return Numeric scalar.
calc_gmean_safe <- function(x) {
  vals <- x[!is.na(x) & x > 0]
  if(length(vals) == 0) return(1e-6) # Fallback for empty/zero vectors
  return(exp(mean(log(vals))))
}

#' @title Internal Median Fallback Helper
#' @description Fills NAs with column medians. Used when BPCA fails.
#' @param mat Numeric matrix.
#' @return Numeric matrix with NAs filled.
.impute_fallback_median <- function(mat) {
  apply(mat, 2, function(x) {
    if(all(is.na(x))) return(rep(0, length(x))) # Handle 100% NA
    x[is.na(x)] <- median(x, na.rm = TRUE)
    return(x)
  })
}

# --- MAIN FUNCTIONS ---

#' @title Zero Replacement via Count Zero Multiplicative (CZM) method
#' @description 
#' Replaces zeros in a compositional matrix using Bayesian-Multiplicative replacement.
#' Handles the "NA conflict" by temporarily masking NAs.
coda_replace_zeros <- function(mat) {
  
  if (!is.matrix(mat)) stop("Input must be a matrix.")
  
  if (!any(mat == 0, na.rm = TRUE)) {
    message("   [CoDa] No zeros detected. Skipping zero replacement.")
    return(mat)
  }
  
  message("   [CoDa] Zeros detected. Applying zCompositions::cmultRepl (CZM)...")
  
  na_locs <- is.na(mat)
  has_na <- any(na_locs)
  mat_temp <- mat
  
  # 1. Temporary Fill (if NAs exist) to allow cmultRepl to run on Zeros
  if (has_na) {
    for(j in 1:ncol(mat_temp)) {
      if(any(na_locs[,j])) {
        fill_val <- calc_gmean_safe(mat_temp[,j])
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
        # Refactored: Use shared helper logic
        mat_fixed[na_idx, j] <- calc_gmean_safe(mat_fixed[,j])
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
  
  # Optimization: Use vectorized math to prevent 'apply' dimension dropping on 1-col matrices
  # Algebra: log(x / gmean(x)) == log(x) - mean(log(x))
  log_mat <- log(mat)
  
  # Calculate row means (equivalent to log of geometric means), ignoring NAs
  row_means <- rowMeans(log_mat, na.rm = TRUE)
  
  # R handles row-wise subtraction via recycling
  clr_mat <- log_mat - row_means
  
  return(clr_mat)
}

#' @title Compute Local ILR for Compositional Subgroups
#' @description 
#' Iterates through defined groups. If a group has >1 marker, 
#' calculates Isometric Log-Ratio coordinates.
coda_compute_local_ilr <- function(mat_imputed, hybrid_groups) {
  
  if (any(mat_imputed <= 0, na.rm = TRUE)) stop("ILR input must be strictly positive.")
  
  ilr_results <- list()
  message("   [CoDa] Computing Local ILR balances for sub-compositions...")
  
  for (group_name in names(hybrid_groups)) {
    target_mks <- hybrid_groups[[group_name]]
    present_mks <- intersect(target_mks, colnames(mat_imputed))
    
    if (length(present_mks) < 2) next
    
    mat_sub <- mat_imputed[, present_mks, drop = FALSE]
    
    ilr_obj <- compositions::ilr(compositions::acomp(mat_sub))
    ilr_mat <- as.matrix(ilr_obj)
    rownames(ilr_mat) <- rownames(mat_imputed)
    
    ilr_results[[group_name]] <- ilr_mat
    message(sprintf("      -> '%s': Computed %d ILR balances.", group_name, ncol(ilr_mat)))
  }
  
  return(ilr_results)
}

#' @title Robust Logit Transformation (Strict Percentage Mode)
#' @description Transforms percentages (0-100) into Logit space (Log-Odds).
coda_transform_logit <- function(mat, epsilon = 1e-6, input_type = "percentage") {
  
  p <- mat
  
  if (input_type == "percentage") {
    p <- mat / 100
  }
  
  # Boundary Protection (Clamping)
  n_low  <- sum(p < epsilon, na.rm = TRUE)
  n_high <- sum(p > (1 - epsilon), na.rm = TRUE)
  
  if (n_low > 0 || n_high > 0) {
    message(sprintf("   [CoDa] Clamped %d boundary values (0 or 1) to epsilon range.", n_low + n_high))
  }
  
  p[p < epsilon] <- epsilon
  p[p > (1 - epsilon)] <- 1 - epsilon
  
  logit_mat <- log(p / (1 - p))
  return(logit_mat)
}

#' @title Bayesian PCA Imputation (Robust & Fail-Safe)
#' @description 
#' Imputes missing values using Bayesian PCA.
#' Implements a multi-stage safety mechanism:
#' 1. Ignores strictly constant columns.
#' 2. Catches errors/ambiguities in kEstimate (CV).
#' 3. Catches errors in PCA execution.
#' 
#' @param mat Numeric matrix (samples x markers).
#' @param nPcs Number of components. "auto" uses kEstimate (CV).
#' @param seed Integer. Seed for reproducibility.
#' @return A numeric matrix with NAs filled.
#' @export
impute_matrix_bpca <- function(mat, nPcs = "auto", seed = 123) {
  
  # 0. Basic Validation
  if (!any(is.na(mat)) && !any(is.infinite(mat))) return(mat)
  
  if (nrow(mat) < 5 || ncol(mat) < 2) {
    warning("[Impute] Matrix too small (<5 samples or <2 cols). Falling back to Median.")
    return(.impute_fallback_median(mat)) 
  }
  
  mat[is.infinite(mat)] <- NA
  
  if (!requireNamespace("pcaMethods", quietly = TRUE)) {
    stop("Package 'pcaMethods' required.")
  }
  
  # 1. Pre-Processing: Handle Constant Columns
  vars <- apply(mat, 2, var, na.rm = TRUE)
  is_const <- !is.na(vars) & (vars < 1e-12)
  const_cols <- names(which(is_const))
  var_cols   <- names(which(!is_const))
  
  mat_bpca_in <- mat[, var_cols, drop = FALSE]
  
  if (ncol(mat_bpca_in) < 2) {
    warning("[Impute] Too few variable columns after filtering. Using Median fallback.")
    return(.impute_fallback_median(mat))
  }
  
  # 2. Estimate K (Cross-Validation) - The "Fragile" Step
  # Define the mathematical hard limit: min(n, p) - 1
  hard_limit_k <- min(nrow(mat_bpca_in), ncol(mat_bpca_in)) - 1
  if (hard_limit_k < 1) hard_limit_k <- 1
  
  k_used <- 2 # Default safe value
  
  if (is.null(nPcs) || nPcs == "auto") {
    tryCatch({
      set.seed(seed)
      max_k <- min(5, hard_limit_k)
      
      utils::capture.output({
        k_est <- pcaMethods::kEstimate(mat_bpca_in, method = "bpca", evalPcs = 1:max_k, verbose = FALSE, scale = "none") 
      })
      suggested_k <- k_est$bestNPcs
      
      # --- SMART FIX: Handling Ambiguity ---
      # Case 1: kEstimate returns multiple values (e.g., 4, 5) -> Take min
      if (length(suggested_k) > 1) {
        k_used <- min(suggested_k)
        warning(sprintf("   [Impute] kEstimate ambiguous (returned %s). Using conservative K=%d.", 
                        paste(suggested_k, collapse=","), k_used))
        
        # Case 2: kEstimate returns invalid/null -> Fallback to 2
      } else if (length(suggested_k) == 0 || is.na(suggested_k)) {
        warning("[Impute] kEstimate returned NA/NULL. Defaulting to K=2.")
        k_used <- 2
        
        # Case 3: kEstimate > hard limit -> Clamp
      } else if (suggested_k > hard_limit_k) {
        warning(sprintf("   [Impute] kEstimate returned invalid K=%d (Max allowed=%d). Clamping.", 
                        suggested_k, hard_limit_k))
        k_used <- hard_limit_k
        
        # Case 4: Valid single K
      } else {
        k_used <- suggested_k
      }
      
    }, error = function(e) {
      warning(sprintf("   [Impute] kEstimate (CV) failed: '%s'. Defaulting to K=2.", e$message))
      k_used <<- 2 
    })
  } else {
    k_used <- as.numeric(nPcs)
  }
  
  # --- FINAL SAFETY CLAMP ---
  if (k_used > hard_limit_k) {
    message(sprintf("   [Impute] Clamping K from %d to %d (max rank).", k_used, hard_limit_k))
    k_used <- hard_limit_k
  }
  
  # 3. Run BPCA (The Execution)
  mat_bpca_out <- mat_bpca_in 
  
  tryCatch({
    set.seed(seed)
    res_bpca <- pcaMethods::pca(mat_bpca_in, method = "bpca", nPcs = k_used, verbose = FALSE)
    mat_bpca_out <- as.matrix(pcaMethods::completeObs(res_bpca))
    
  }, error = function(e) {
    warning(sprintf("   [Impute] BPCA execution failed (K=%d): '%s'. Applying Median Fallback.", k_used, e$message))
    mat_bpca_out <<- .impute_fallback_median(mat_bpca_in)
  })
  
  # 4. Reconstruct Full Matrix
  mat_final <- mat 
  mat_final[, var_cols] <- mat_bpca_out
  
  if (length(const_cols) > 0) {
    for (col in const_cols) {
      val <- mean(mat[, col], na.rm = TRUE)
      if (is.na(val)) val <- 0 
      mat_final[is.na(mat_final[, col]), col] <- val
    }
  }
  
  return(mat_final)
}

#' @title Perform Hybrid Transformation Pipeline
#' @export
perform_hybrid_transformation <- function(mat_raw, config, mode = "complete") {
  
  if (!mode %in% c("complete", "fast")) stop("Mode must be 'complete' or 'fast'.")
  
  message(sprintf("\n[CoDa] Starting Hybrid Strategy (%s mode)...", mode))
  
  rn_safe <- rownames(mat_raw)
  if (is.null(rn_safe)) stop("[FATAL] Input matrix 'mat_raw' lacks rownames (Patient IDs).")
  
  all_current_markers <- colnames(mat_raw)
  comp_markers_list <- unique(unlist(config$hybrid_groups))
  comp_markers <- intersect(all_current_markers, comp_markers_list)
  func_markers <- setdiff(all_current_markers, comp_markers)
  if (length(func_markers) > 0) {
    message(paste("[CoDa] The following markers will be considered as Functional (Logit):", 
                  paste(func_markers, collapse=", ")))
  }
  
  mat_comp_trans <- NULL
  mat_func_trans <- NULL
  ilr_results <- list()
  
  # --- TRACK A: Compositional Markers (CLR) ---
  if (length(comp_markers) > 0) {
    if (mode == "complete") message("   [Track A] Processing Compositional Groups (CZM -> CLR)...")
    mat_comp_raw <- mat_raw[, comp_markers, drop = FALSE]
    
    # Drop 100% NA cols
    na_counts <- colSums(is.na(mat_comp_raw))
    empty_cols <- which(na_counts == nrow(mat_comp_raw))
    if (length(empty_cols) > 0) {
      mat_comp_raw <- mat_comp_raw[, -empty_cols, drop = FALSE]
    }
    
    if (ncol(mat_comp_raw) > 0) {
      if (mode == "complete") {
        mat_comp_clean <- coda_replace_zeros(mat_comp_raw)
        if (any(is.na(mat_comp_clean))) {
          mat_comp_clean <- impute_nas_simple(mat_comp_clean)
        }
      } else {
        mat_comp_clean <- mat_comp_raw
        mat_comp_clean[mat_comp_clean <= 0 | is.na(mat_comp_clean)] <- 1e-6
      }
      rownames(mat_comp_clean) <- rn_safe 
      
      mat_comp_trans <- coda_transform_clr(mat_comp_clean)
      
      if (mode == "complete") {
        ilr_results <- coda_compute_local_ilr(mat_comp_clean, config$hybrid_groups)
      }
    }
  }
  
  # --- TRACK B: Functional Markers (Logit) ---
  if (length(func_markers) > 0) {
    if (mode == "complete") message("   [Track B] Processing Functional Markers (LOD -> Logit -> BPCA)...")
    mat_func_raw <- mat_raw[, func_markers, drop = FALSE]
    
    na_counts <- colSums(is.na(mat_func_raw))
    empty_cols <- which(na_counts == nrow(mat_func_raw))
    if (length(empty_cols) > 0) {
      mat_func_raw <- mat_func_raw[, -empty_cols, drop = FALSE]
    }
    
    if (ncol(mat_func_raw) > 0) {
      mat_func_lod <- mat_func_raw
      
      if (mode == "complete") {
        mins <- apply(mat_func_lod, 2, function(x) {
          pos <- x[x > 0 & !is.na(x)]
          if(length(pos) > 0) min(pos)/2 else 1e-6
        })
        for (j in 1:ncol(mat_func_lod)) {
          zeros <- which(mat_func_lod[, j] == 0 & !is.na(mat_func_lod[, j]))
          if (length(zeros) > 0) mat_func_lod[zeros, j] <- mins[j]
        }
      } else {
        mat_func_lod[mat_func_lod <= 0] <- 1e-6
      }
      
      mat_func_logit <- coda_transform_logit(mat_func_lod, input_type = "percentage")
      
      if (mode == "complete") {
        if (any(colSums(!is.na(mat_func_logit)) == 0)) {
          stop("[FATAL] Functional matrix has empty columns after Logit.")
        }
        
        seed_val <- if(!is.null(config$stats$seed)) config$stats$seed else 123
        mat_func_trans <- impute_matrix_bpca(mat_func_logit, nPcs = "auto", seed = seed_val)
      } else {
        mat_func_trans <- mat_func_logit
      }
      rownames(mat_func_trans) <- rn_safe
    }
  }
  
  # --- MERGE & FINALIZE ---
  list_parts <- list()
  if (!is.null(mat_comp_trans)) list_parts[[1]] <- mat_comp_trans
  if (!is.null(mat_func_trans)) list_parts[[2]] <- mat_func_trans
  
  if (length(list_parts) == 0) stop("[FATAL] No data remaining after filtering.")
  
  mat_hybrid_raw <- do.call(cbind, list_parts)
  mat_hybrid_raw <- mat_hybrid_raw[, sort(colnames(mat_hybrid_raw)), drop = FALSE]
  
  mat_hybrid_z <- NULL
  if (mode == "complete") {
    message("   [Norm] Applying Global Z-Score Standardization...")
    mat_hybrid_z <- scale(mat_hybrid_raw)
    attr(mat_hybrid_z, "scaled:center") <- NULL
    attr(mat_hybrid_z, "scaled:scale") <- NULL
    mat_hybrid_z <- as.matrix(mat_hybrid_z)
    rownames(mat_hybrid_z) <- rn_safe
  } else {
    mat_hybrid_z <- mat_hybrid_raw
  }
  
  return(list(
    hybrid_data_raw = mat_hybrid_raw,
    hybrid_data_z = mat_hybrid_z,
    ilr_balances = ilr_results,
    hybrid_markers = colnames(mat_hybrid_z)
  ))
}