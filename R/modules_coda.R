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
  
  p[p < epsilon] <- epsilon
  p[p > (1 - epsilon)] <- 1 - epsilon
  
  # 3. Transform
  logit_mat <- log(p / (1 - p))
  return(logit_mat)
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
    tryCatch({
      set.seed(seed)
      max_k <- min(5, hard_limit_k)
      
      utils::capture.output({
        k_est <- pcaMethods::kEstimate(mat_bpca_in, method = "bpca", evalPcs = 1:max_k, verbose = FALSE, scale = "none") 
      })
      suggested_k <- k_est$bestNPcs
      
      if (length(suggested_k) > 1) {
        k_used <- min(suggested_k)
        warning(sprintf("   [Impute] kEstimate ambiguous. Using conservative K=%d.", k_used))
      } else if (length(suggested_k) == 0 || is.na(suggested_k)) {
        k_used <- 2
      } else if (suggested_k > hard_limit_k) {
        k_used <- hard_limit_k
      } else {
        k_used <- suggested_k
      }
      
    }, error = function(e) {
      k_used <<- 2 
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
#' - Handles LOD (Limit of Detection) Zeros
#' - Applies Logit transformation (if enabled via config)
#' - Performs BPCA Imputation
#' - Z-Score Standardization
#' @export
perform_data_transformation <- function(mat_raw, config, mode = "complete") {
  
  if (!mode %in% c("complete", "fast")) stop("Mode must be 'complete' or 'fast'.")
  
  eps_val <- if(!is.null(config$imputation$epsilon)) config$imputation$epsilon else 1e-6
  input_fmt <- if(!is.null(config$input_format)) config$input_format else "percentage"
  
  # Future-proof architecture: check if transformation method is defined, default to logit
  trans_method <- if(!is.null(config$transformation$method)) config$transformation$method else "logit"
  
  message(sprintf("\n[Transform] Starting Strategy (%s mode). Method: %s. Input Format: %s...", mode, trans_method, input_fmt))
  
  rn_safe <- rownames(mat_raw)
  mat_trans <- NULL
  
  # Remove 100% NA columns
  na_counts <- colSums(is.na(mat_raw))
  empty_cols <- which(na_counts == nrow(mat_raw))
  if (length(empty_cols) > 0) {
    mat_raw <- mat_raw[, -empty_cols, drop = FALSE]
  }
  
  if (ncol(mat_raw) > 0) {
    mat_lod <- mat_raw
    
    if (trans_method == "logit") {
      # Handle Zeros / LOD specifically for logit
      if (mode == "complete") {
        mins <- apply(mat_lod, 2, function(x) {
          pos <- x[x > 0 & !is.na(x)]
          if(length(pos) == 0) return(eps_val)
          return(min(pos) / 2)
        })
        for (j in 1:ncol(mat_lod)) {
          zeros <- which(mat_lod[, j] == 0 & !is.na(mat_lod[, j]))
          if (length(zeros) > 0) mat_lod[zeros, j] <- mins[j]
        }
      } else {
        mat_lod[mat_lod <= 0] <- eps_val 
      }
      
      # Apply Logit
      mat_transformed <- coda_transform_logit(
        mat = mat_lod, 
        epsilon = eps_val, 
        input_type = input_fmt 
      )
    } else if (trans_method == "none") {
      # Pass-through for unbounded data (e.g. pg/mL Cytokines)
      mat_transformed <- mat_lod
    } else {
      stop(sprintf("[Transform] Unknown transformation method: %s", trans_method))
    }
    
    # Imputation
    if (mode == "complete") {
      seed_val <- if(!is.null(config$stats$seed)) config$stats$seed else 123
      mat_trans <- impute_matrix_bpca(mat_transformed, nPcs = "auto", seed = seed_val)
    } else {
      mat_trans <- mat_transformed
    }
    rownames(mat_trans) <- rn_safe
  }
  
  if (is.null(mat_trans)) stop("[FATAL] No data remaining after filtering.")
  
  mat_final_raw <- mat_trans[, sort(colnames(mat_trans)), drop = FALSE]
  mat_final_z <- NULL
  
  if (mode == "complete") {
    mat_final_z <- scale(mat_final_raw)
    attr(mat_final_z, "scaled:center") <- NULL
    attr(mat_final_z, "scaled:scale") <- NULL
    mat_final_z <- as.matrix(mat_final_z)
    rownames(mat_final_z) <- rn_safe
  } else {
    mat_final_z <- mat_final_raw
  }
  
  return(list(
    hybrid_data_raw = mat_final_raw,
    hybrid_data_z = mat_final_z,
    hybrid_markers = colnames(mat_final_z)
  ))
}