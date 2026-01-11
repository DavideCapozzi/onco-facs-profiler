# R/modules_multivariate.R
# ==============================================================================
# MULTIVARIATE ANALYSIS MODULE
# Description: Wrappers for PLS-DA and Supervised Feature Selection.
# Dependencies: mixOmics, dplyr
# ==============================================================================

library(dplyr)

#' @title Run PLS-DA
run_plsda_model <- function(data_z, metadata, group_col = "Group", n_comp = 2) {
  
  # Y must be a factor
  Y <- as.factor(metadata[[group_col]])
  X <- data_z
  
  # Basic safety check for N
  if (nrow(X) < 6) {
    warning("[PLS-DA] Very low sample size. Cross-validation may fail.")
  }
  
  message(sprintf("   [PLS-DA] Fitting model with %d components...", n_comp))
  
  model <- mixOmics::plsda(X, Y, ncomp = n_comp)
  
  # Cross Validation (Critical for validity)
  # Dynamic folds: 5-fold if N>20, else Leave-One-Out (LOO)
  val_method <- if(nrow(X) > 20) "Mfold" else "loo"
  n_folds <- if(val_method == "Mfold") 5 else NULL
  
  perf_pls <- mixOmics::perf(model, validation = val_method, folds = n_folds, 
                             progressBar = FALSE, nrepeat = 10, auc = TRUE)
  
  return(list(model = model, performance = perf_pls))
}

#' @title Extract PLS-DA Performance Metrics
#' @description Returns the Balanced Error Rate (BER) and AUC to check model validity.
extract_plsda_performance <- function(pls_result) {
  perf <- pls_result$performance
  
  # Extract BER for the chosen component number
  n_comp <- pls_result$model$ncomp
  ber_vals <- perf$error.rate$BER[, "max.dist"] # max.dist is standard for PLSDA
  
  df_perf <- data.frame(
    Component = 1:n_comp,
    Overall_BER = ber_vals[1:n_comp],
    Validation_Method = perf$validation
  )
  
  return(df_perf)
}

#' @title Extract Top Discriminant Features (Loadings)
extract_plsda_loadings <- function(pls_result, top_n = 15) {
  
  model <- pls_result$model
  
  # Safe extraction (mixOmics structure can vary)
  if(is.null(model$loadings$X)) return(data.frame())
  
  # Extract Loadings for Component 1
  loadings_df <- data.frame(
    Marker = rownames(model$loadings$X),
    Comp1_Weight = model$loadings$X[, 1]
  )
  
  # Add Comp2 if exists
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