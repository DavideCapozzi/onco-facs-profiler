# R/workflows.R
# ==============================================================================
# WORKFLOW ORCHESTRATION
# ==============================================================================

#' @title Run Comparative Analysis Workflow
#' @description 
#' Executes the full statistical pipeline for a specific contrast.
#' 
#' @param data_list The processed data object (output of step 01).
#' @param scenario A list containing the specific scenario config (cases/controls).
#' @param config The global configuration object.
#' @param output_root Directory for saving plot-based outputs (PDFs).
#' @return A list containing results: $id, $permanova, $drivers, $network.
run_comparative_workflow <- function(data_list, scenario, config, output_root) {
  
  # ============================================================================
  # 1.0 INITIALIZATION
  # ============================================================================
  
  # 1.1 Setup Output Directory
  out_dir <- file.path(output_root, scenario$id)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  message(sprintf("\n[Workflow] Starting Scenario: %s", scenario$id))
  
  # ============================================================================
  # 2.0 DATA PREPARATION
  # ============================================================================
  
  # 2.1 Extract Metadata and Matrix
  full_meta <- data_list$metadata
  full_mat  <- data_list$hybrid_data_z %>% 
    dplyr::select(dplyr::all_of(data_list$hybrid_markers)) %>% 
    as.matrix()
  rownames(full_mat) <- data_list$metadata$Patient_ID
  
  # 2.2 Apply Stratification Logic derived from Config
  if (!is.null(config$stratification$mode) && config$stratification$mode == "selection") {
    strat_col <- config$stratification$column
    if (strat_col %in% colnames(full_meta)) {
      full_meta$Group <- full_meta[[strat_col]]
    } else {
      warning(sprintf("   [Warn] Stratification column '%s' not found. Using existing 'Group'.", strat_col))
    }
  }
  
  target_groups <- unlist(c(scenario$case_groups, scenario$control_groups))
  
  # 2.3 Filter Metadata
  sub_meta <- full_meta %>% 
    dplyr::filter(Group %in% target_groups) %>% 
    dplyr::mutate(Analysis_Group = ifelse(Group %in% unlist(scenario$case_groups), 
                                          scenario$case_label, 
                                          scenario$control_label)) %>% 
    dplyr::mutate(Analysis_Group = factor(Analysis_Group, 
                                          levels = c(scenario$control_label, scenario$case_label)))
  
  # 2.4 Filter Matrix
  sub_mat <- full_mat[sub_meta$Patient_ID, , drop = FALSE]
  
  # 2.5 Validation: Check if samples exist
  if (nrow(sub_mat) == 0) {
    message("   [Stop] No samples found after filtering. Check Scenario groups.")
    return(NULL)
  }
  
  counts <- table(sub_meta$Analysis_Group)
  message(sprintf("   [Info] Samples found: %s=%d vs %s=%d", 
                  names(counts)[1], counts[1], names(counts)[2], counts[2]))
  
  # ============================================================================
  # 3.0 DATA SANITY CHECKS & UNIVERSAL SPACE PRESERVATION
  # ============================================================================
  
  # 3.1 Variance Sanity Check (Micro-Jitter Approach)
  # Instead of dropping low-variance markers and breaking the Universal Feature Space,
  # we inject a microscopic noise to prevent singularity errors in correlation matrices.
  var_thresh <- if(!is.null(config$stats$min_variance)) config$stats$min_variance else 1e-6
  
  vars <- apply(sub_mat, 2, var, na.rm = TRUE)
  low_vars <- vars < var_thresh | is.na(vars)
  
  if (sum(low_vars) > 0) {
    dropped_names <- names(vars)[low_vars]
    message(sprintf("   [Prep] Warning: %d features have variance < %s.", sum(low_vars), as.character(var_thresh)))
    message("   [Prep] Injecting micro-jitter to preserve Universal Feature Space and prevent mathematical singularities.")
    
    # Inject microscopic noise (sd = 1e-8) to prevent division by zero in cor/pcor
    # This ensures the correlation calculation completes and returns ~0 for these flat features.
    set.seed(if(!is.null(config$stats$seed)) config$stats$seed else 123)
    for (col in dropped_names) {
      sub_mat[, col] <- sub_mat[, col] + rnorm(nrow(sub_mat), mean = 0, sd = 1e-8)
    }
  }
  
  # ============================================================================
  # 4.0 STATISTICAL TESTS
  # ============================================================================
  
  # 4.1 PERMANOVA Inference
  message("   [Stats] Running PERMANOVA...")
  perm_res <- NULL
  tryCatch({
    sub_data_input <- cbind(sub_meta, as.data.frame(sub_mat))
    meta_cols_to_exclude <- colnames(sub_meta)
    
    perm_obj <- test_coda_permanova(
      sub_data_input, 
      group_col = "Analysis_Group", 
      metadata_cols = meta_cols_to_exclude,
      n_perm = if(!is.null(config$stats$n_perm)) config$stats$n_perm else 999
    )
    perm_res <- as.data.frame(perm_obj)
  }, error = function(e) {
    message(sprintf("      [Fail] PERMANOVA failed: %s", e$message))
  })
  
  # ============================================================================
  # 5.0 MULTIVARIATE ANALYSIS
  # ============================================================================
  
  # 5.1 sPLS-DA Modeling
  message("   [Stats] Running sPLS-DA...")
  spls_drivers <- NULL
  
  valid_levels <- names(counts)[counts > 0]
  
  if (length(valid_levels) < 2) {
    message(sprintf("      [Skip] sPLS-DA skipped: Only %d group(s) present.", length(valid_levels)))
  } else {
    tryCatch({
      n_comp <- if(!is.null(config$multivariate$n_comp)) config$multivariate$n_comp else 2
      cv_folds <- if(!is.null(config$multivariate$validation_folds)) config$multivariate$validation_folds else 5
      n_rep <- if(!is.null(config$multivariate$n_repeat_cv)) config$multivariate$n_repeat_cv else 10
      
      min_class_n <- min(counts[valid_levels])
      actual_folds <- min(cv_folds, min_class_n)
      
      if (min_class_n < 2) {
        message("      [Skip] sPLS-DA skipped: Class size too small for CV (<2).")
      } else {
        if (actual_folds < 2) actual_folds <- 2
        
        spls_res <- run_splsda_model(
          data_z = sub_mat,
          metadata = sub_meta, 
          group_col = "Analysis_Group", 
          n_comp = n_comp, 
          folds = actual_folds, 
          n_repeat = n_rep
        )
        
        spls_drivers <- extract_plsda_loadings(spls_res)
        
        if (nrow(spls_drivers) > 0) {
          
          # 5.2 Dynamic Color Mapping
          scen_palette <- c()
          scen_palette[scenario$control_label] <- config$colors$control
          scen_palette[scenario$case_label] <- if(length(config$colors$cases) > 0) config$colors$cases[[1]] else "firebrick"
          
          viz_report_plsda(
            pls_res = spls_res, 
            drivers_df = spls_drivers, 
            metadata_viz = sub_meta, 
            colors_viz = scen_palette, 
            out_path = file.path(out_dir, paste0(scenario$id, "_Plot_sPLSDA.pdf")),
            group_col = "Analysis_Group"
          )
          
          # 5.3 Dynamic Renaming of Extractable Drivers
          contrast_lbl <- paste0(scenario$case_label, "_VS_", scenario$control_label)
          
          spls_drivers <- spls_drivers %>% 
            dplyr::rename_with(
              .fn = ~ paste0("Weight_", contrast_lbl, "_PC1"), 
              .cols = dplyr::matches("^Comp1_Weight$")
            ) %>%
            dplyr::rename_with(
              .fn = ~ paste0("Weight_", contrast_lbl, "_PC2"), 
              .cols = dplyr::matches("^Comp2_Weight$")
            )
        }
      }
    }, error = function(e) {
      message(sprintf("      [Fail] sPLS-DA failed: %s", e$message))
    })
  }
  
  return(list(
    id = scenario$id,
    permanova = perm_res,
    drivers = spls_drivers
  ))
}