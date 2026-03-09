# R/workflows.R
# ==============================================================================
# WORKFLOWS & MACROS MODULE
# Description: Orchestration functions that combine multiple analytical steps 
#              (Filtering, Hypothesis Testing, Multivariate, Viz).
# Dependencies: dplyr, modules_hypothesis, modules_multivariate, modules_viz
# ==============================================================================

library(dplyr)

#' @title Execute Complete Statistical Pipeline for a Single Scenario
#' @description 
#' Filters data, handles zero-variance, checks dispersion assumptions, 
#' runs PERMANOVA, fits sPLS-DA, and generates plots.
#' 
#' @param scenario List containing scenario definitions from config.
#' @param full_mat Global numeric matrix (Z-scored).
#' @param full_meta Global metadata dataframe.
#' @param config The complete configuration object.
#' @param out_dir Directory path to save scenario-specific PDF plots.
#' @return A list of dataframes: 'dispersion', 'permanova', and 'drivers'.
run_scenario_pipeline <- function(scenario, full_mat, full_meta, config, out_dir) {
  
  # 1. Subsetting Data
  target_groups <- unlist(c(scenario$case_groups, scenario$control_groups))
  
  sub_meta <- full_meta %>% 
    dplyr::filter(Group %in% target_groups) %>% 
    dplyr::mutate(Analysis_Group = ifelse(Group %in% unlist(scenario$case_groups), 
                                          scenario$case_label, 
                                          scenario$control_label)) %>% 
    dplyr::mutate(Analysis_Group = factor(Analysis_Group, 
                                          levels = c(scenario$control_label, scenario$case_label)))
  
  sub_mat <- full_mat[sub_meta$Patient_ID, , drop = FALSE]
  
  if (nrow(sub_mat) == 0) {
    message("      [Skip] No samples found after filtering.")
    return(NULL)
  }
  
  counts <- table(sub_meta$Analysis_Group)
  message(sprintf("      [Info] Samples: %s=%d vs %s=%d", 
                  names(counts)[1], counts[1], names(counts)[2], counts[2]))
  
  # 2. Zero-Variance Filtering (Universal Space Protection)
  var_thresh <- if(!is.null(config$stats$min_variance)) config$stats$min_variance else 1e-6
  vars <- apply(sub_mat, 2, var, na.rm = TRUE)
  low_vars <- vars < var_thresh | is.na(vars)
  
  if (sum(low_vars) > 0) {
    dropped_names <- names(vars)[low_vars]
    message(sprintf("      [Prep] Dropping %d flat features with zero pooled variance: %s", 
                    sum(low_vars), paste(dropped_names, collapse = ", ")))
    sub_mat <- sub_mat[, !low_vars, drop = FALSE]
  }
  
  if (ncol(sub_mat) < 2) {
    message("      [Skip] Less than 2 features remaining after variance filter.")
    return(NULL)
  }
  
  sub_data_input <- cbind(sub_meta, as.data.frame(sub_mat))
  res_list <- list(dispersion = NULL, permanova = NULL, drivers = NULL)
  
  # 3. Beta-Dispersion Check (Strictly Scenario-Specific)
  tryCatch({
    disp_obj <- test_coda_dispersion(
      data_input = sub_data_input, 
      group_col = "Analysis_Group", 
      metadata_cols = colnames(sub_meta),
      n_perm = if(!is.null(config$stats$n_perm)) config$stats$n_perm else 999
    )
    res_list$dispersion <- as.data.frame(disp_obj$anova_table)
  }, error = function(e) message(sprintf("      [Fail] Dispersion test failed: %s", e$message)))
  
  # 4. PERMANOVA
  tryCatch({
    perm_obj <- test_coda_permanova(
      data_input = sub_data_input, 
      group_col = "Analysis_Group", 
      metadata_cols = colnames(sub_meta),
      n_perm = if(!is.null(config$stats$n_perm)) config$stats$n_perm else 999
    )
    res_list$permanova <- as.data.frame(perm_obj)
  }, error = function(e) message(sprintf("      [Fail] PERMANOVA failed: %s", e$message)))
  
  # 5. sPLS-DA & Visualization
  valid_levels <- names(counts)[counts > 0]
  if (length(valid_levels) >= 2 && min(counts[valid_levels]) >= 2) {
    tryCatch({
      n_comp <- if(!is.null(config$multivariate$n_comp)) config$multivariate$n_comp else 2
      cv_folds <- if(!is.null(config$multivariate$validation_folds)) config$multivariate$validation_folds else 5
      actual_folds <- min(cv_folds, min(counts[valid_levels]))
      if (actual_folds < 2) actual_folds <- 2
      
      spls_res <- run_splsda_model(
        data_z = sub_mat,
        metadata = sub_meta, 
        group_col = "Analysis_Group", 
        n_comp = n_comp, 
        folds = actual_folds, 
        n_repeat = if(!is.null(config$multivariate$n_repeat_cv)) config$multivariate$n_repeat_cv else 10
      )
      
      drivers_df <- extract_plsda_loadings(spls_res)
      
      if (nrow(drivers_df) > 0) {
        # Setup specific palette for plotting
        scen_palette <- c()
        scen_palette[scenario$control_label] <- config$colors$control
        
        # Safe extraction for case color to avoid out of bounds
        case_col <- "firebrick"
        if(!is.null(config$colors$groups[[scenario$case_label]])) {
          case_col <- config$colors$groups[[scenario$case_label]]
        } else if (length(config$colors$cases) > 0) {
          case_col <- config$colors$cases[[1]]
        }
        scen_palette[scenario$case_label] <- case_col
        
        viz_report_plsda(
          pls_res = spls_res, 
          drivers_df = drivers_df, 
          metadata_viz = sub_meta, 
          colors_viz = scen_palette, 
          out_path = file.path(out_dir, paste0(scenario$id, "_Plot_sPLSDA.pdf")),
          group_col = "Analysis_Group"
        )
        
        # Dynamic renaming for Excel readability
        contrast_lbl <- paste0(scenario$case_label, "_VS_", scenario$control_label)
        drivers_df <- drivers_df %>% 
          dplyr::rename_with(.fn = ~ paste0("Weight_", contrast_lbl, "_PC1"), .cols = dplyr::matches("^Comp1_Weight$")) %>%
          dplyr::rename_with(.fn = ~ paste0("Weight_", contrast_lbl, "_PC2"), .cols = dplyr::matches("^Comp2_Weight$"))
        
        res_list$drivers <- drivers_df
      }
    }, error = function(e) message(sprintf("      [Fail] sPLS-DA failed: %s", e$message)))
  } else {
    message("      [Skip] sPLS-DA skipped: Class size too small for CV (<2).")
  }
  
  return(res_list)
}