# R/modules_interpretation.R
# ==============================================================================
# INTERPRETATION & DECODING MODULE
# Description: Decodes ILR balances by correlating them with Local CLR parts.
# Dependencies: dplyr, stats
# ==============================================================================

library(dplyr)
library(stats)

#' @title Decode ILR Balances (The Rosetta Stone)
#' @description 
#' Correlates abstract ILR coordinates with the CLR-transformed original markers.
#' CRITICAL: 'raw_parts_mat' must be RAW percentages/counts, not already transformed data.
decode_ilr_to_clr <- function(ilr_mat, raw_parts_mat, p_threshold = 0.05) {
  
  # Validation
  if(nrow(ilr_mat) != nrow(raw_parts_mat)) {
    warning("[Decoder] Row mismatch between ILR and Raw data. Skipping.")
    return(data.frame())
  }
  
  # 1. Compute Local CLR for the raw parts
  # gmean function robust to small zeros (add epsilon if strict zero)
  gmean_safe <- function(x) {
    x_safe <- x
    x_safe[x <= 0 | is.na(x)] <- 1e-6 # Epsilon for zeros/NAs
    exp(mean(log(x_safe)))
  }
  
  # Apply CLR: log( x / gmean(x) )
  # We do this sample-by-sample (row-wise)
  clr_mat <- t(apply(raw_parts_mat, 1, function(x) {
    if(all(is.na(x))) return(rep(NA, length(x)))
    gm <- gmean_safe(x)
    return(log((x + 1e-6) / gm)) # Epsilon added to numerator to prevent log(0)
  }))
  colnames(clr_mat) <- colnames(raw_parts_mat)
  
  # 2. Correlate ILR Cols with CLR Cols
  results <- data.frame()
  balances <- colnames(ilr_mat)
  markers  <- colnames(clr_mat)
  
  for (bal in balances) {
    for (mk in markers) {
      
      # Pairwise complete obs handles NAs if present in raw data
      test <- tryCatch({
        cor.test(ilr_mat[, bal], clr_mat[, mk], 
                 method = "spearman", 
                 use = "pairwise.complete.obs")
      }, error = function(e) NULL)
      
      if (!is.null(test) && !is.na(test$p.value) && test$p.value < p_threshold) {
        results <- rbind(results, data.frame(
          Balance_ID = bal,
          Marker = mk,
          Correlation = test$estimate,
          P_Value = test$p.value,
          Interpretation = ifelse(test$estimate > 0, 
                                  "Positively correlates (Increases Balance)", 
                                  "Negatively correlates (Decreases Balance)")
        ))
      }
    }
  }
  
  if (nrow(results) > 0) {
    results <- results %>% 
      arrange(Balance_ID, desc(abs(Correlation)))
  }
  
  return(results)
}