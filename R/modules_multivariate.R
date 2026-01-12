# R/modules_multivariate.R
# ==============================================================================
# MULTIVARIATE ANALYSIS MODULE
# Description: Wrappers for PLS-DA and Supervised Feature Selection.
# Dependencies: mixOmics, dplyr
# ==============================================================================

library(dplyr)

#' @title Run PLS-DA (Partial Least Squares Discriminant Analysis)
#' @description 
#' Fits a supervised PLS-DA model. Includes SAFE cross-validation.
#' If CV fails, it returns the model anyway so we don't lose the drivers.
run_plsda_model <- function(data_z, metadata, group_col = "Group", n_comp = 2) {
  
  # Assumption: mixOmics is already loaded by the main script
  
  Y <- as.factor(metadata[[group_col]])
  X <- as.matrix(data_z)
  
  if (nrow(X) != length(Y)) stop("Dimension mismatch between Data and Group labels.")
  
  message(sprintf("   [PLS-DA] Fitting model with %d components...", n_comp))
  
  # 1. Run Model
  model <- mixOmics::plsda(X, Y, ncomp = n_comp)
  
  # 2. Cross-Validation
  # Wrap in tryCatch to prevent crashing if CV fails
  perf_pls <- tryCatch({
    val_method <- if(nrow(X) > 20) "Mfold" else "loo"
    n_folds <- if(val_method == "Mfold") 5 else NULL
    
    # Check for too few samples for Mfold
    min_samples <- min(table(Y))
    if (val_method == "Mfold" && min_samples < n_folds) {
      n_folds <- min_samples
    }
    
    mixOmics::perf(model, validation = val_method, folds = n_folds, 
                   progressBar = FALSE, nrepeat = 10, auc = TRUE)
  }, error = function(e) {
    message(paste("      [WARN] Cross-validation (BER) failed:", e$message))
    return(NULL) 
  })
  
  return(list(model = model, performance = perf_pls))
}

#' @title Extract PLS-DA Performance Metrics
#' @description Returns BER if available, otherwise a placeholder.
extract_plsda_performance <- function(pls_result) {
  perf <- pls_result$performance
  n_comp <- pls_result$model$ncomp
  
  # Handle case where perf failed (NULL)
  if (is.null(perf)) {
    return(data.frame(
      Component = 1:n_comp,
      Overall_BER = NA,
      Validation_Method = "FAILED"
    ))
  }
  
  # Extract BER safely
  ber_vals <- rep(NA, n_comp)
  
  # Robust extraction of BER
  tryCatch({
    if (!is.null(perf$error.rate$BER)) {
      if (is.matrix(perf$error.rate$BER)) {
        raw_ber <- perf$error.rate$BER[, "max.dist"]
        if(length(raw_ber) >= n_comp) ber_vals <- raw_ber[1:n_comp]
      } else if (is.vector(perf$error.rate$BER)) {
        ber_vals[1] <- perf$error.rate$BER["max.dist"]
      }
    }
  }, error = function(e) {
    message("      [WARN] BER extraction encountered an issue. Using NAs.")
  })
  
  # [FIXED] Ensure Validation Method handles empty vectors (length 0)
  val_method <- "Unknown"
  if (!is.null(perf$validation)) {
    if (length(perf$validation) > 0) {
      val_method <- as.character(perf$validation)[1] # Force scalar
    }
  }
  
  df_perf <- data.frame(
    Component = 1:n_comp,
    Overall_BER = ber_vals,
    Validation_Method = val_method 
  )
  
  return(df_perf)
}

#' @title Extract Top Discriminant Features (Loadings)
extract_plsda_loadings <- function(pls_result, top_n = 15) {
  model <- pls_result$model
  if(is.null(model$loadings$X)) return(data.frame())
  
  loadings_df <- data.frame(
    Marker = rownames(model$loadings$X),
    Comp1_Weight = model$loadings$X[, 1]
  )
  
  if(model$ncomp > 1) {
    loadings_df$Comp2_Weight <- model$loadings$X[, 2]
  }
  
  loadings_df$Importance <- abs(loadings_df$Comp1_Weight)
  
  top_drivers <- loadings_df %>%
    arrange(desc(Importance)) %>%
    head(top_n) %>%
    mutate(
      Direction = ifelse(Comp1_Weight > 0, "Positive_Assoc", "Negative_Assoc")
    )
  
  return(top_drivers)
}