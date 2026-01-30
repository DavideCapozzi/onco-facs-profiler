# R/workflows.R
# ==============================================================================
# WORKFLOW ORCHESTRATION
# Description: Orchestrates analysis for specific scenarios defined in config.
#              Connects data preparation with statistical modules.
# ==============================================================================

source("R/modules_hypothesis.R")
source("R/modules_multivariate.R")
source("R/modules_network.R")
source("R/modules_viz.R")

#' @title Run Comparative Analysis Workflow
#' @description 
#' Executes the full statistical pipeline for a specific contrast.
#' 1. Subsets data based on scenario groups.
#' 2. Runs PERMANOVA (Global multivariate difference).
#' 3. Runs sPLS-DA (Driver identification).
#' 4. Runs Robust Differential Network Inference (Bootstrap + Permutation).
#' 
#' @param data_list The processed data object (output of step 01).
#' @param scenario A list containing the specific scenario config (cases/controls).
#' @param config The global configuration object.
#' @param output_root Directory for saving plot-based outputs (PDFs).
#' @return A list containing results: $id, $permanova, $drivers, $network.
run_comparative_workflow <- function(data_list, scenario, config, output_root) {
  
  # 1. Setup Output Directory for Plots
  out_dir <- file.path(output_root, scenario$id)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  message(sprintf("\n[Workflow] Starting Scenario: %s", scenario$id))
  
  # 2. Data Preparation & Subsetting
  full_meta <- data_list$metadata
  # Use hybrid_data_z for statistics (Z-scored data)
  full_mat  <- data_list$hybrid_data_z %>% 
    dplyr::select(dplyr::all_of(data_list$hybrid_markers)) %>% 
    as.matrix()
  rownames(full_mat) <- data_list$metadata$Patient_ID
  
  target_groups <- c(scenario$case_groups, scenario$control_groups)
  
  # Filter Metadata and Create Analysis_Group Factor
  # Levels are set so Control is reference (first level)
  sub_meta <- full_meta %>% 
    dplyr::filter(Group %in% target_groups) %>%
    dplyr::mutate(Analysis_Group = ifelse(Group %in% scenario$case_groups, 
                                          scenario$case_label, 
                                          scenario$control_label)) %>%
    dplyr::mutate(Analysis_Group = factor(Analysis_Group, 
                                          levels = c(scenario$control_label, scenario$case_label)))
  
  # Filter Matrix to match Metadata
  sub_mat <- full_mat[sub_meta$Patient_ID, , drop = FALSE]
  
  # Check if we have data after filtering
  if (nrow(sub_mat) == 0) {
    message("   [Stop] No samples found for this scenario after filtering.")
    return(NULL)
  }
  
  # 3. Variance Sanity Check (Dynamic Config)
  # Uses config parameter or defaults to 1e-6 if missing
  var_thresh <- if(!is.null(config$stats$min_variance)) config$stats$min_variance else 1e-6
  
  vars <- apply(sub_mat, 2, var, na.rm = TRUE)
  keep_vars <- vars > var_thresh & !is.na(vars)
  
  if (sum(!keep_vars) > 0) {
    dropped <- names(vars)[!keep_vars]
    message(sprintf("   [Prep] Dropping %d features with variance < %s in this subset.", length(dropped), as.character(var_thresh)))
    sub_mat <- sub_mat[, keep_vars, drop = FALSE]
  }
  
  # 4. Statistical Tests (PERMANOVA)
  message("   [Stats] Running PERMANOVA...")
  perm_res <- NULL
  tryCatch({
    # Prepare input dataframe
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
  
  # 5. Multivariate Analysis (sPLS-DA)
  message("   [Stats] Running sPLS-DA...")
  spls_drivers <- NULL
  
  # Pre-check for group levels to avoid "Y should be a factor" error
  n_levels_y <- length(unique(as.character(sub_meta$Analysis_Group)))
  
  if (n_levels_y < 2) {
    message(sprintf("      [Skip] sPLS-DA skipped: Only %d group(s) present (needs 2).", n_levels_y))
  } else {
    tryCatch({
      # Config parameters
      n_comp <- if(!is.null(config$multivariate$n_comp)) config$multivariate$n_comp else 2
      cv_folds <- if(!is.null(config$multivariate$validation_folds)) config$multivariate$validation_folds else 5
      n_rep <- if(!is.null(config$multivariate$n_repeat_cv)) config$multivariate$n_repeat_cv else 10
      
      # Adjust folds if sample size is too small
      min_class_n <- min(table(sub_meta$Analysis_Group))
      # Ensure folds is at least 2 and not larger than N
      if (cv_folds > min_class_n) cv_folds <- max(2, min_class_n)
      
      # Double check if CV is possible (need at least 2 samples per class for LOO/Mfold)
      if (min_class_n < 2) {
        message("      [Skip] sPLS-DA skipped: Class size too small for CV (<2).")
      } else {
        spls_res <- run_splsda_model(
          data_z = sub_mat,          
          metadata = sub_meta,       
          group_col = "Analysis_Group",
          n_comp = n_comp, 
          folds = cv_folds,
          n_repeat = n_rep 
        )
        
        spls_drivers <- extract_plsda_loadings(spls_res)
        
        # Plotting is handled here as it generates a PDF file
        if (nrow(spls_drivers) > 0) {
          viz_report_plsda(
            pls_res = spls_res,
            drivers_df = spls_drivers,
            metadata_viz = sub_meta,
            colors_viz = get_palette(config),
            out_path = file.path(out_dir, "Plot_sPLSDA.pdf"),
            group_col = "Analysis_Group"
          )
        }
      }
    }, error = function(e) {
      message(sprintf("      [Fail] sPLS-DA failed: %s", e$message))
    })
  }
  
  # 6. Network Inference (Robust Differential)
  message("   [Network] Inferring Differential Networks (Robust PCOR)...")
  net_res <- NULL
  
  tryCatch({
    mat_case <- sub_mat[sub_meta$Analysis_Group == scenario$case_label, , drop=FALSE]
    mat_ctrl <- sub_mat[sub_meta$Analysis_Group == scenario$control_label, , drop=FALSE]
    
    # --- DYNAMIC CONFIG CHECK ---
    # Retrieve minimum N from config or default to 5
    net_min_n <- if(!is.null(config$stats$min_network_n)) config$stats$min_network_n else 5
    
    # Check minimum sample size for network inference using dynamic threshold
    if(nrow(mat_case) < net_min_n || nrow(mat_ctrl) < net_min_n) {
      message(sprintf("      [Skip] Network inference skipped (N < %d in one group: Case=%d, Ctrl=%d).", 
                      net_min_n, nrow(mat_case), nrow(mat_ctrl)))
    } else {
      
      # Configuration parameters extraction
      n_boot_cfg <- if(!is.null(config$stats$n_boot)) config$stats$n_boot else 100
      n_perm_cfg <- if(!is.null(config$stats$n_perm)) config$stats$n_perm else 1000
      fdr_cfg <- if(!is.null(config$stats$fdr_threshold)) config$stats$fdr_threshold else 0.1
      seed_cfg <- if(!is.null(config$stats$seed)) config$stats$seed else 123
      n_cores_cfg <- if(!is.null(config$stats$n_cores) && config$stats$n_cores != "auto") config$stats$n_cores else 1
      
      # Calling the Robust Module (Bootstrap + Permutation)
      net_res <- run_differential_network(
        mat_ctrl = mat_ctrl,
        mat_case = mat_case,
        n_boot = n_boot_cfg,
        n_perm = n_perm_cfg,
        seed = seed_cfg,
        n_cores = n_cores_cfg,
        fdr_thresh = fdr_cfg,
        stability_thresh = 0.8 
      )
    }
  }, error = function(e) {
    message(sprintf("      [Fail] Network Inference failed: %s", e$message))
  })
  
  # Return structured list for Main Orchestrator
  return(list(
    id = scenario$id,
    permanova = perm_res,
    drivers = spls_drivers,
    network = net_res # Contains $edges_table and $networks
  ))
}