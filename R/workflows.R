# R/workflows.R
# ==============================================================================
# WORKFLOW ORCHESTRATION
# Description: High-level functions to execute analysis pipelines based on config.
# ==============================================================================

source("R/modules_hypothesis.R")
source("R/modules_multivariate.R")
source("R/modules_network.R")
source("R/modules_viz.R")

#' @title Run Comparative Analysis Workflow
#' @description 
#' Executes the full statistical pipeline for a specific contrast defined in the config.
#' Handles subsetting, zero-variance filtering, PERMANOVA, sPLS-DA, and Network Inference.
#' 
#' @param data_list The processed data object (output of step 01).
#' @param scenario A list containing the specific scenario config.
#' @param config The global configuration object.
#' @param output_root The root directory where scenario-specific subfolders will be created.
#' @return A list containing the results (drivers, stats, networks).
run_comparative_workflow <- function(data_list, scenario, config, output_root) {
  
  # 1. Setup Output Directory
  out_dir <- file.path(output_root, scenario$id)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  message(sprintf("\n[Workflow] Starting Scenario: %s", scenario$id))
  message(sprintf("           Comparison: %s (%s) vs %s (%s)", 
                  scenario$case_label, paste(scenario$case_groups, collapse="+"),
                  scenario$control_label, paste(scenario$control_groups, collapse="+")))
  
  # 2. Data Preparation & Subsetting
  full_meta <- data_list$metadata
  full_mat  <- data_list$hybrid_data_z %>% 
    dplyr::select(dplyr::all_of(data_list$hybrid_markers)) %>% 
    as.matrix()
  rownames(full_mat) <- data_list$metadata$Patient_ID
  
  target_groups <- c(scenario$case_groups, scenario$control_groups)
  
  # Filter Metadata and Create Analysis_Group
  sub_meta <- full_meta %>% 
    dplyr::filter(Group %in% target_groups) %>%
    dplyr::mutate(Analysis_Group = ifelse(Group %in% scenario$case_groups, 
                                          scenario$case_label, 
                                          scenario$control_label)) %>%
    dplyr::mutate(Analysis_Group = factor(Analysis_Group, levels = c(scenario$control_label, scenario$case_label)))
  
  # Filter Matrix to match Metadata
  sub_mat <- full_mat[sub_meta$Patient_ID, , drop = FALSE]
  
  # 3. Variance Sanity Check
  vars <- apply(sub_mat, 2, var, na.rm = TRUE)
  keep_vars <- vars > 1e-6 & !is.na(vars)
  
  if (sum(!keep_vars) > 0) {
    dropped <- names(vars)[!keep_vars]
    message(sprintf("   [Prep] Dropping %d features with zero variance in this subset: %s", 
                    length(dropped), paste(dropped, collapse=", ")))
    sub_mat <- sub_mat[, keep_vars, drop = FALSE]
  }
  
  # 4. Construct Clean Input for Modules (Dataframe version for PERMANOVA)
  sub_data_input <- cbind(sub_meta, as.data.frame(sub_mat))
  
  # Robustly identify metadata columns to exclude
  meta_cols_to_exclude <- colnames(sub_meta)
  
  # 5. Statistical Tests (PERMANOVA)
  message("   [Stats] Running PERMANOVA...")
  tryCatch({
    res_perm <- test_coda_permanova(
      sub_data_input, 
      group_col = "Analysis_Group", 
      metadata_cols = meta_cols_to_exclude, 
      n_perm = 999
    )
    write.csv(as.data.frame(res_perm), file.path(out_dir, "PERMANOVA_Result.csv"))
  }, error = function(e) {
    message(sprintf("      [Fail] PERMANOVA failed: %s", e$message))
  })
  
  # 6. Multivariate Analysis (sPLS-DA)
  message("   [Stats] Running sPLS-DA...")
  
  local_colors <- get_palette(config)
  spls_drivers <- NULL
  
  tryCatch({
    # CORRECTED CALL: Pass data_z (matrix) and metadata separately
    spls_res <- run_splsda_model(
      data_z = sub_mat,          # Numeric Matrix ONLY
      metadata = sub_meta,       # Metadata DataFrame
      group_col = "Analysis_Group",
      n_comp = 2, 
      folds = 5,
      n_repeat = 10 # Reduced for speed, increase for final run
    )
    
    # CORRECTED EXTRACTION: Must use extract_plsda_loadings
    spls_drivers <- extract_plsda_loadings(spls_res)
    
    if (nrow(spls_drivers) > 0) {
      write.csv(spls_drivers, file.path(out_dir, "sPLSDA_Drivers.csv"), row.names = FALSE)
      
      # Visualization
      viz_report_plsda(
        pls_res = spls_res,
        drivers_df = spls_drivers,
        metadata_viz = sub_meta,
        colors_viz = local_colors,
        out_path = file.path(out_dir, "Plot_sPLSDA.pdf"),
        group_col = "Analysis_Group"
      )
    } else {
      message("      [Info] No significant drivers found.")
    }
  }, error = function(e) {
    message(sprintf("      [Fail] sPLS-DA failed: %s", e$message))
  })
  
  # 7. Network Inference (Differential)
  message("   [Network] Inferring Differential Networks...")
  
  net_res <- NULL
  tryCatch({
    # Use clean numeric matrices
    mat_case <- sub_mat[sub_meta$Analysis_Group == scenario$case_label, , drop=FALSE]
    mat_ctrl <- sub_mat[sub_meta$Analysis_Group == scenario$control_label, , drop=FALSE]
    
    n_boot <- if(!is.null(config$stats$n_boot)) config$stats$n_boot else 100
    
    # Now calls the function we added to modules_network.R
    net_res <- infer_differential_network(
      mat_ctrl = mat_ctrl,
      mat_case = mat_case,
      method = "spearman",
      n_bootstrap = n_boot
    )
    
    if (!is.null(net_res$diff_edges)) {
      write.csv(net_res$diff_edges, file.path(out_dir, "Differential_Network_Edges.csv"), row.names = FALSE)
    }
  }, error = function(e) {
    message(sprintf("      [Fail] Network Inference failed: %s", e$message))
  })
  
  return(list(
    id = scenario$id,
    drivers = spls_drivers,
    network = net_res
  ))
}