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
  
  # 1. Setup Output Directory
  out_dir <- file.path(output_root, scenario$id)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  message(sprintf("\n[Workflow] Starting Scenario: %s", scenario$id))
  
  # 2. Data Preparation & Stratification Logic
  # ----------------------------------------------------------------------------
  full_meta <- data_list$metadata
  full_mat  <- data_list$hybrid_data_z %>% 
    dplyr::select(dplyr::all_of(data_list$hybrid_markers)) %>% 
    as.matrix()
  rownames(full_mat) <- data_list$metadata$Patient_ID
  
  # Apply Stratification Logic derived from Config
  if (!is.null(config$stratification$mode) && config$stratification$mode == "selection") {
    strat_col <- config$stratification$column
    if (strat_col %in% colnames(full_meta)) {
      full_meta$Group <- full_meta[[strat_col]]
    } else {
      warning(sprintf("   [Warn] Stratification column '%s' not found. Using existing 'Group'.", strat_col))
    }
  }
  
  # CRITICAL FIX: unlist() ensures we get a vector of strings, not a list
  target_groups <- unlist(c(scenario$case_groups, scenario$control_groups))
  
  # Filter Metadata
  sub_meta <- full_meta %>% 
    dplyr::filter(Group %in% target_groups) %>% 
    dplyr::mutate(Analysis_Group = ifelse(Group %in% unlist(scenario$case_groups), 
                                          scenario$case_label, 
                                          scenario$control_label)) %>% 
    dplyr::mutate(Analysis_Group = factor(Analysis_Group, 
                                          levels = c(scenario$control_label, scenario$case_label)))
  
  # Filter Matrix
  sub_mat <- full_mat[sub_meta$Patient_ID, , drop = FALSE]
  
  # Check if we have data after filtering
  if (nrow(sub_mat) == 0) {
    message("   [Stop] No samples found after filtering. Check Scenario groups.")
    return(NULL)
  }
  
  # Log found samples
  counts <- table(sub_meta$Analysis_Group)
  message(sprintf("   [Info] Samples found: %s=%d vs %s=%d", 
                  names(counts)[1], counts[1], names(counts)[2], counts[2]))
  
  # 3. Variance Sanity Check
  # ----------------------------------------------------------------------------
  var_thresh <- if(!is.null(config$stats$min_variance)) config$stats$min_variance else 1e-6
  
  vars <- apply(sub_mat, 2, var, na.rm = TRUE)
  keep_vars <- vars > var_thresh & !is.na(vars)
  
  if (sum(!keep_vars) > 0) {
    dropped <- names(vars)[!keep_vars]
    message(sprintf("   [Prep] Dropping %d features with variance < %s.", length(dropped), as.character(var_thresh)))
    sub_mat <- sub_mat[, keep_vars, drop = FALSE]
  }
  
  # 4. Statistical Tests (PERMANOVA)
  # ----------------------------------------------------------------------------
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
  
  # 5. Multivariate Analysis (sPLS-DA)
  # ----------------------------------------------------------------------------
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
          
          # --- DYNAMIC COLOR MAPPING ---
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
          
          # --- DYNAMIC RENAMING (Drivers) ---
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
  
  # 6. Network Inference (Robust Differential)
  # ----------------------------------------------------------------------------
  message("   [Network] Inferring Differential Networks...")
  net_res <- NULL
  
  tryCatch({
    mat_case <- sub_mat[sub_meta$Analysis_Group == scenario$case_label, , drop=FALSE]
    mat_ctrl <- sub_mat[sub_meta$Analysis_Group == scenario$control_label, , drop=FALSE]
    
    net_min_n <- if(!is.null(config$stats$min_network_n)) config$stats$min_network_n else 5
    
    if(nrow(mat_case) < net_min_n || nrow(mat_ctrl) < net_min_n) {
      message(sprintf("      [Skip] Network skipped (N < %d in one group).", net_min_n))
    } else {
      
      n_boot_cfg <- if(!is.null(config$stats$n_boot)) config$stats$n_boot else 100
      n_perm_cfg <- if(!is.null(config$stats$n_perm)) config$stats$n_perm else 1000
      
      # Correct parameter reading
      pval_threshold_cfg <- if(!is.null(config$stats$pvalue_threshold)) config$stats$pvalue_threshold else 0.05
      edge_threshold_cfg <- if(!is.null(config$stats$network_edge_threshold)) config$stats$network_edge_threshold else 0.15
      
      seed_cfg <- if(!is.null(config$stats$seed)) config$stats$seed else 123
      n_cores_cfg <- if(!is.null(config$stats$n_cores) && config$stats$n_cores != "auto") config$stats$n_cores else 1
      
      net_res <- run_differential_network(
        mat_ctrl = mat_ctrl, 
        mat_case = mat_case, 
        n_boot = n_boot_cfg, 
        n_perm = n_perm_cfg, 
        seed = seed_cfg, 
        n_cores = n_cores_cfg, 
        pvalue_thresh = pval_threshold_cfg,
        stability_thresh = 0.8,
        label_ctrl = scenario$control_label, 
        label_case = scenario$case_label,
        magnitude_threshold = edge_threshold_cfg
      )
      
      if (!is.null(net_res)) {
        
        # --- 1. Save RDS (Raw Object) ---
        rds_path <- file.path(out_dir, paste0(scenario$id, "_Differential_Network.rds"))
        saveRDS(net_res, rds_path)
        message(sprintf("      [Output] Network RDS saved: %s", basename(rds_path)))
        
        # --- 2. Advanced Topology Integration ---
        topo_xlsx_path <- file.path(out_dir, paste0(scenario$id, "_Topology_Metrics.xlsx"))
        message(sprintf("      [Export] Generating Advanced Topology Report: %s", basename(topo_xlsx_path)))
        
        wb_topo <- createWorkbook()
        has_topo_data <- FALSE
        
        # A. Calculate Topologies
        topo_ctrl <- NULL
        if (sum(net_res$stability$ctrl) > 0) {
          topo_ctrl <- calculate_node_topology(net_res$stability$ctrl)
          if (!is.null(topo_ctrl)) {
            sh_name <- substr(paste0("Topology_", scenario$control_label), 1, 31)
            addWorksheet(wb_topo, sh_name)
            writeData(wb_topo, sh_name, topo_ctrl)
            has_topo_data <- TRUE
          }
        }
        
        topo_case <- NULL
        if (sum(net_res$stability$case) > 0) {
          topo_case <- calculate_node_topology(net_res$stability$case)
          if (!is.null(topo_case)) {
            sh_name <- substr(paste0("Topology_", scenario$case_label), 1, 31)
            addWorksheet(wb_topo, sh_name)
            writeData(wb_topo, sh_name, topo_case)
            has_topo_data <- TRUE
          }
        }
        
        # B. Node Rewiring (Jaccard)
        if (sum(net_res$stability$ctrl) > 0 && sum(net_res$stability$case) > 0) {
          rewiring_df <- calculate_jaccard_rewiring(net_res$stability$ctrl, net_res$stability$case)
          if (!is.null(rewiring_df)) {
            col_ctrl <- paste0("Degree_", scenario$control_label)
            col_case <- paste0("Degree_", scenario$case_label)
            
            rewiring_df <- rewiring_df %>% 
              dplyr::rename(!!col_ctrl := Degree_Ctrl, 
                            !!col_case := Degree_Case)
            
            addWorksheet(wb_topo, "Rewiring_Analysis")
            writeData(wb_topo, "Rewiring_Analysis", rewiring_df)
            has_topo_data <- TRUE
          }
        }
        
        # C. Hub-Driver Integration
        pdf_path <- file.path(out_dir, paste0(scenario$id, "_Hub_Driver_Report_4Page.pdf"))
        if (!is.null(spls_drivers) && nrow(spls_drivers) > 0) {
          pdf(pdf_path, width = 11, height = 8)
          
          # Helper to plot and save
          process_hub_driver <- function(topo, label, metric, suffix) {
            if (!is.null(topo)) {
              hd_res <- integrate_hub_drivers(spls_drivers, topo, topo_col = metric)
              if (!is.null(hd_res)) {
                sh_name <- substr(paste0("Hub_Driver_", label, "_", substr(metric,1,3)), 1, 31)
                if(!sh_name %in% names(wb_topo)) {
                  addWorksheet(wb_topo, sh_name)
                  writeData(wb_topo, sh_name, hd_res)
                }
                p <- plot_hub_driver_quadrant(hd_res, y_label = metric, title_suffix = suffix)
                print(p)
              }
            }
          }
          
          process_hub_driver(topo_ctrl, scenario$control_label, "Degree", paste0("\nNetwork: ", scenario$control_label))
          process_hub_driver(topo_case, scenario$case_label, "Degree", paste0("\nNetwork: ", scenario$case_label))
          process_hub_driver(topo_ctrl, scenario$control_label, "Betweenness", paste0("\nNetwork: ", scenario$control_label))
          process_hub_driver(topo_case, scenario$case_label, "Betweenness", paste0("\nNetwork: ", scenario$case_label))
          
          dev.off()
          message(sprintf("      [Viz] Saved 4-Page Hub-Driver Report: %s", basename(pdf_path)))
          has_topo_data <- TRUE # If we wrote sheets
        }
        
        # D. Differential Edges
        if (!is.null(net_res$edges_table) && nrow(net_res$edges_table) > 0) {
          addWorksheet(wb_topo, "Differential_Edges")
          writeData(wb_topo, "Differential_Edges", net_res$edges_table)
          has_topo_data <- TRUE
        }
        
        if (has_topo_data) {
          saveWorkbook(wb_topo, topo_xlsx_path, overwrite = TRUE)
        }
        
        # --- 3. Visualization: Edge Density (Distribution) ---
        if (!is.null(net_res$networks)) {
          dist_pdf_path <- file.path(out_dir, paste0(scenario$id, "_Edge_Distribution.pdf"))
          pdf(dist_pdf_path, width = 8, height = 6)
          
          if(sum(abs(net_res$networks$ctrl) > 0) > 0) {
            print(viz_plot_edge_density(net_res$networks$ctrl, 
                                        threshold = edge_threshold_cfg, 
                                        group_label = scenario$control_label))
          }
          if(sum(abs(net_res$networks$case) > 0) > 0) {
            print(viz_plot_edge_density(net_res$networks$case, 
                                        threshold = edge_threshold_cfg, 
                                        group_label = scenario$case_label))
          }
          dev.off()
          message(sprintf("      [Viz] Saved Edge Distribution: %s", basename(dist_pdf_path)))
        }
        
        # --- 4. Network Graph Visualization ---
        if ((sum(net_res$stability$ctrl) > 0 || sum(net_res$stability$case) > 0)) {
          viz_path <- file.path(out_dir, paste0(scenario$id, "_Plot_Networks.pdf"))
          pdf(viz_path, width = 12, height = 6)
          try({
            if(sum(net_res$stability$ctrl) > 0) {
              print(plot_network_structure(net_res$stability$ctrl, net_res$networks$ctrl, 
                                           title = paste("Control:", scenario$control_label), min_cor = 0))
            }
            if(sum(net_res$stability$case) > 0) {
              print(plot_network_structure(net_res$stability$case, net_res$networks$case, 
                                           title = paste("Case:", scenario$case_label), min_cor = 0))
            }
          })
          dev.off()
          message("      [Viz] Network plots saved.")
        }
        
        # --- 5. Cytoscape Export ---
        message("      [Export] Generating Rich Cytoscape Tables...")
        cyto_dir <- file.path(out_dir, "cytoscape_export")
        if (!dir.exists(cyto_dir)) dir.create(cyto_dir)
        
        if (sum(net_res$stability$ctrl) > 0) {
          tbl_ctrl <- get_rich_network_table(net_res$stability$ctrl, net_res$networks$ctrl, group_label = scenario$control_label)
          if (!is.null(tbl_ctrl)) write.csv(tbl_ctrl, file.path(cyto_dir, paste0(scenario$control_label, "_RichTable.csv")), row.names = FALSE)
        }
        if (sum(net_res$stability$case) > 0) {
          tbl_case <- get_rich_network_table(net_res$stability$case, net_res$networks$case, group_label = scenario$case_label)
          if (!is.null(tbl_case)) write.csv(tbl_case, file.path(cyto_dir, paste0(scenario$case_label, "_RichTable.csv")), row.names = FALSE)
        }
      }
    }
  }, error = function(e) {
    message(sprintf("      [Fail] Network Inference failed: %s", e$message))
  })
  
  return(list(
    id = scenario$id,
    permanova = perm_res,
    drivers = spls_drivers,
    network = net_res
  ))
}