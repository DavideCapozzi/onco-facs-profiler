# R/modules_coda.R
# ==============================================================================
# DATA TRANSFORMATION MODULE
# Description: Handling imputation (BPCA) and transformation (Logit).
# Dependencies: dplyr, pcaMethods
# ==============================================================================

library(dplyr)

# --- HELPER FUNCTIONS ---

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

#' @title Robust Logit Transformation
#' @description 
#' Transforms values into Logit space (Log-Odds).
#' Handles scale normalization based on input type (Percentage vs Proportion).
#' 
#' @param mat Numeric matrix.
#' @param epsilon Boundary protection value.
#' @param input_type String. "percentage" (divides by 100) or "proportion" (keeps as is).
#' @return Numeric matrix in logit space.
coda_transform_logit <- function(mat, epsilon = 1e-6, input_type = "percentage") {
  
  p <- mat
  
  # 1. Normalize Scale
  if (input_type == "percentage") {
    if (any(mat > 100, na.rm = TRUE)) {
      stop("[Transform] Critical Error: Values > 100 detected in percentage mode.")
    }
    p <- mat / 100
  } else {
    if (any(mat > 1.0 + epsilon, na.rm = TRUE)) {
      stop("[Transform] Critical Error: Values > 1.0 detected in proportion mode.")
    }
  }
  
  # 2. Boundary Protection (Clamping)
  n_low  <- sum(p < epsilon, na.rm = TRUE)
  n_high <- sum(p > (1 - epsilon), na.rm = TRUE)
  
  if (n_low > 0 || n_high > 0) {
    message(sprintf("   [Transform] Clamped %d boundary values.", n_low + n_high))
  }
  
  # Explicitly track clamped values
  is_clamped <- !is.na(p) & (p < epsilon | p > (1 - epsilon))
  
  p[p < epsilon] <- epsilon
  p[p > (1 - epsilon)] <- 1 - epsilon
  
  # 3. Transform
  logit_mat <- log(p / (1 - p))
  
  return(list(
    data = logit_mat,
    clamped_mask = is_clamped
  ))
}

#' @title Bayesian PCA Imputation (Robust & Fail-Safe)
#' @description Imputes missing values using Bayesian PCA.
#' @param mat Numeric matrix (samples x markers).
#' @param nPcs Number of components. "auto" uses kEstimate (CV).
#' @param seed Integer. Seed for reproducibility.
#' @return A numeric matrix with NAs filled.
#' @export
impute_matrix_bpca <- function(mat, nPcs = "auto", seed = 123) {
  
  if (!any(is.na(mat)) && !any(is.infinite(mat))) return(mat)
  
  if (nrow(mat) < 5 || ncol(mat) < 2) {
    warning("[Impute] Matrix too small (<5 samples or <2 cols). Falling back to Median.")
    return(.impute_fallback_median(mat)) 
  }
  
  mat[is.infinite(mat)] <- NA
  
  if (!requireNamespace("pcaMethods", quietly = TRUE)) stop("Package 'pcaMethods' required.")
  
  vars <- apply(mat, 2, var, na.rm = TRUE)
  is_const <- !is.na(vars) & (vars < 1e-12)
  const_cols <- names(which(is_const))
  var_cols   <- names(which(!is_const))
  
  mat_bpca_in <- mat[, var_cols, drop = FALSE]
  
  if (ncol(mat_bpca_in) < 2) {
    warning("[Impute] Too few variable columns after filtering. Using Median fallback.")
    return(.impute_fallback_median(mat))
  }
  
  hard_limit_k <- min(nrow(mat_bpca_in), ncol(mat_bpca_in)) - 1
  if (hard_limit_k < 1) hard_limit_k <- 1
  
  k_used <- 2 
  
  if (is.null(nPcs) || nPcs == "auto") {
    k_used <- tryCatch({
      set.seed(seed)
      max_k <- min(5, hard_limit_k)
      
      utils::capture.output({
        k_est <- pcaMethods::kEstimate(mat_bpca_in, method = "bpca", evalPcs = 1:max_k, verbose = FALSE, scale = "none") 
      })
      suggested_k <- k_est$bestNPcs
      
      # Return suggested value
      if (length(suggested_k) > 1) {
        warning(sprintf("   [Impute] kEstimate ambiguous. Using conservative K=%d.", min(suggested_k)))
        min(suggested_k)
      } else if (length(suggested_k) == 0 || is.na(suggested_k)) {
        2
      } else if (suggested_k > hard_limit_k) {
        hard_limit_k
      } else {
        suggested_k
      }
    }, error = function(e) {
      warning("   [Impute] kEstimate failed. Falling back to K=2.")
      2 
    })
  } else {
    k_used <- as.numeric(nPcs)
  }
  
  if (k_used > hard_limit_k) k_used <- hard_limit_k
  
  mat_bpca_out <- mat_bpca_in 
  
  tryCatch({
    set.seed(seed)
    res_bpca <- pcaMethods::pca(mat_bpca_in, method = "bpca", nPcs = k_used, verbose = FALSE)
    mat_bpca_out <- as.matrix(pcaMethods::completeObs(res_bpca))
    
  }, error = function(e) {
    warning(sprintf("   [Impute] BPCA execution failed. Applying Median Fallback."))
    mat_bpca_out <<- .impute_fallback_median(mat_bpca_in)
  })
  
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

#' @title Perform Data Transformation Pipeline
#' @description 
#' Orchestrates the transformation:
#' - Unified scale management (Raw vs Proportion)
#' - Robust Auto-Epsilon calculation (Median of bottom 1% of strictly positive values)
#' - Handles LOD (Limit of Detection) Zeros dynamically per column
#' - Applies Logit transformation
#' - Performs BPCA Imputation
#' - Z-Score Standardization
#' @export
perform_data_transformation <- function(mat_raw, config, mode = "complete") {
  
  if (!mode %in% c("complete", "fast")) stop("Mode must be 'complete' or 'fast'.")
  
  rn_safe <- rownames(mat_raw)
  input_fmt <- if(!is.null(config$input_format)) config$input_format else "percentage"
  trans_method <- if(!is.null(config$transformation$method)) config$transformation$method else "logit"
  
  # Remove 100% NA columns before any operation
  na_counts <- colSums(is.na(mat_raw))
  empty_cols <- which(na_counts == nrow(mat_raw))
  if (length(empty_cols) > 0) {
    mat_raw <- mat_raw[, -empty_cols, drop = FALSE]
  }
  
  if (ncol(mat_raw) == 0) stop("[FATAL] All columns are empty.")
  
  # --------------------------------------------------------------------------
  # 1. ROBUST AUTO-EPSILON CALCULATION
  # --------------------------------------------------------------------------
  raw_eps_conf <- if(!is.null(config$imputation$epsilon)) config$imputation$epsilon else 1e-6
  
  # Temporary conversion to proportion strictly for robust boundary estimation
  mat_prop_temp <- mat_raw
  if (input_fmt == "percentage") {
    mat_prop_temp <- mat_raw / 100
  }
  
  if (is.character(raw_eps_conf) && tolower(raw_eps_conf) == "auto") {
    pos_vals <- mat_prop_temp[mat_prop_temp > 0 & !is.na(mat_prop_temp)]
    
    if (length(pos_vals) > 0) {
      threshold_1pct <- quantile(pos_vals, 0.01, na.rm = TRUE)
      bottom_1pct_vals <- pos_vals[pos_vals <= threshold_1pct]
      
      # Divide by 2 (standard for LOD imputation) and enforce a strict mathematical floor
      calculated_eps <- as.numeric(median(bottom_1pct_vals) / 2)
      eps_prop <- max(calculated_eps, 1e-8) 
    } else {
      eps_prop <- 1e-6
    }
    
    eps_raw_display <- eps_prop * (if(input_fmt == "percentage") 100 else 1)
    message(sprintf("   [Transform] Epsilon 'auto' selected. Robust LOD boundary set to %s (Proportion space: %s)", 
                    format(eps_raw_display, scientific = TRUE), format(eps_prop, scientific = TRUE)))
  } else {
    # If user provided a numeric value, it is treated as the target proportion epsilon
    eps_prop <- as.numeric(raw_eps_conf)
  }
  
  # Scale epsilon back to raw space for LOD matrix adjustments (preventing scale-mismatch in "fast" mode)
  eps_raw <- eps_prop * (if(input_fmt == "percentage") 100 else 1)
  
  message(sprintf("\n[Transform] Starting Strategy (%s mode). Method: %s. Input Format: %s...", mode, trans_method, input_fmt))
  
  # Initialize clamping mask for tracking
  global_clamped_mask <- matrix(FALSE, nrow = nrow(mat_raw), ncol = ncol(mat_raw))
  rownames(global_clamped_mask) <- rownames(mat_raw)
  colnames(global_clamped_mask) <- colnames(mat_raw)
  
  mat_trans <- NULL
  
  # --------------------------------------------------------------------------
  # 2. LOD HANDLING & TRANSFORMATION
  # --------------------------------------------------------------------------
  if (trans_method == "logit") {
    
    mat_lod <- mat_raw
    
    # Handle Zeros / LOD specifically for logit, in RAW space
    if (mode == "complete") {
      mins_raw <- apply(mat_lod, 2, function(x) {
        pos <- x[x > 0 & !is.na(x)]
        if(length(pos) == 0) return(eps_raw)
        return(min(pos) / 2)
      })
      
      for (j in 1:ncol(mat_lod)) {
        zeros <- which(mat_lod[, j] == 0 & !is.na(mat_lod[, j]))
        if (length(zeros) > 0) {
          mat_lod[zeros, j] <- mins_raw[j]
          global_clamped_mask[zeros, j] <- TRUE
        }
      }
    } else {
      # Fast mode uses the unified eps_raw to avoid scale-crushing bugs
      zeros <- which(mat_lod <= 0 & !is.na(mat_lod))
      if (length(zeros) > 0) global_clamped_mask[zeros] <- TRUE
      mat_lod[mat_lod <= 0] <- eps_raw 
    }
    
    # Apply Logit (API naturally expects mat_lod in raw space and epsilon in proportion space)
    logit_res <- coda_transform_logit(
      mat = mat_lod, 
      epsilon = eps_prop, 
      input_type = input_fmt 
    )
    mat_transformed <- logit_res$data
    
    if (!is.null(logit_res$clamped_mask)) {
      global_clamped_mask <- global_clamped_mask | logit_res$clamped_mask
    }
    
  } else if (trans_method == "none") {
    mat_transformed <- mat_raw
  } else {
    stop(sprintf("[Transform] Unknown transformation method: %s", trans_method))
  }
  
  # --------------------------------------------------------------------------
  # 3. IMPUTATION & SCALING
  # --------------------------------------------------------------------------
  if (mode == "complete") {
    seed_val <- if(!is.null(config$stats$seed)) config$stats$seed else 123
    mat_trans <- impute_matrix_bpca(mat_transformed, nPcs = "auto", seed = seed_val)
  } else {
    mat_trans <- mat_transformed
  }
  rownames(mat_trans) <- rn_safe
  
  # The output raw matrix must reflect the LOD adjustments but remain in raw scale
  mat_final_raw <- mat_lod[, sort(colnames(mat_lod)), drop = FALSE] 
  mat_trans <- mat_trans[, sort(colnames(mat_trans)), drop = FALSE]
  
  mat_final_z <- NULL
  
  if (mode == "complete") {
    # Protect scaling from zero-variance columns post-imputation
    col_vars_post <- apply(mat_trans, 2, var, na.rm = TRUE)
    valid_cols_post <- !is.na(col_vars_post) & (col_vars_post > 1e-12)
    
    if (sum(!valid_cols_post) > 0) {
      n_dropped <- sum(!valid_cols_post)
      message(sprintf("   [Transform] Dropping %d columns with zero variance post-imputation to protect scaling.", n_dropped))
      
      mat_final_raw <- mat_final_raw[, valid_cols_post, drop = FALSE]
      mat_trans <- mat_trans[, valid_cols_post, drop = FALSE]
      global_clamped_mask <- global_clamped_mask[, valid_cols_post, drop = FALSE]
    }
    
    if (ncol(mat_trans) == 0) stop("[FATAL] All columns dropped due to zero variance post-imputation.")
    
    # Safe scaling
    mat_final_z <- scale(mat_trans)
    attr(mat_final_z, "scaled:center") <- NULL
    attr(mat_final_z, "scaled:scale") <- NULL
    mat_final_z <- as.matrix(mat_final_z)
    rownames(mat_final_z) <- rn_safe
  } else {
    # Fast mode strictly serves as a geometric proxy; scaling is bypassed intentionally.
    mat_final_z <- mat_trans
    message("   [Transform] Fast mode active: Outputting proxy matrix without Z-scoring.")
  }
  
  return(list(
    hybrid_data_raw = mat_final_raw,
    hybrid_data_z = mat_final_z,
    hybrid_markers = colnames(mat_final_z),
    clamped_mask = global_clamped_mask
  ))
}