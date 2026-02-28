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
  # 3.0 QUALITY CONTROL
  # ============================================================================
  
  # 3.1 Variance Sanity Check
  var_thresh <- if(!is.null(config$stats$min_variance)) config$stats$min_variance else 1e-6
  
  vars <- apply(sub_mat, 2, var, na.rm = TRUE)
  keep_vars <- vars > var_thresh & !is.na(vars)
  
  if (sum(!keep_vars) > 0) {
    dropped <- names(vars)[!keep_vars]
    message(sprintf("   [Prep] Dropping %d features with variance < %s.", length(dropped), as.character(var_thresh)))
    sub_mat <- sub_mat[, keep_vars, drop = FALSE]
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
  
  # ============================================================================
  # 6.0 NETWORK INFERENCE & CYTOSCAPE EXPORT
  # ============================================================================
  
  # 6.1 Infer Differential Networks
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
        label_ctrl = scenario$control_label, 
        label_case = scenario$case_label,
        threshold_type = "percentile", 
        threshold_value = 0.85
      )
      
      if (!is.null(net_res)) {
        
        # 6.2 Save Raw RDS Object
        rds_path <- file.path(out_dir, paste0(scenario$id, "_Differential_Network.rds"))
        saveRDS(net_res, rds_path)
        message(sprintf("      [Output] Network RDS saved: %s", basename(rds_path)))
        
        # 6.3 Topology Generation and Integration
        topo_xlsx_path <- file.path(out_dir, paste0(scenario$id, "_Topology_Metrics.xlsx"))
        wb_topo <- createWorkbook()
        has_topo_data <- FALSE
        
        # 6.3.1 Topologies
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
        
        # 6.3.2 Rewiring Architecture (Jaccard)
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
        
        # 6.3.3 Hub-Driver Integration Visuals
        pdf_path <- file.path(out_dir, paste0(scenario$id, "_Hub_Driver_Report_4Page.pdf"))
        if (!is.null(spls_drivers) && nrow(spls_drivers) > 0) {
          pdf(pdf_path, width = 11, height = 8)
          
          process_hub_driver <- function(topo, label, metric, suffix) {
            if (!is.null(topo)) {
              hd_res <- integrate_hub_drivers(spls_drivers, topo, topo_col = metric)
              if (!is.null(hd_res)) {
                sh_name <- substr(paste0("Hub_Driver_", label, "_", substr(metric,1,3)), 1, 31)
                if(metric == "Degree" && !sh_name %in% names(wb_topo)) {
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
          has_topo_data <- TRUE
        }
        
        # 6.3.4 Edge Table Extractor
        if (!is.null(net_res$edges_table) && nrow(net_res$edges_table) > 0) {
          addWorksheet(wb_topo, "Differential_Edges")
          writeData(wb_topo, "Differential_Edges", net_res$edges_table)
          
          sig_edges_strict <- net_res$edges_table %>%
            dplyr::filter(Edge_Category != "Weak" & P_Value < 0.05) %>%
            dplyr::arrange(P_Value)
          
          if (nrow(sig_edges_strict) > 0) {
            addWorksheet(wb_topo, "Significant_Edges")
            writeData(wb_topo, "Significant_Edges", sig_edges_strict)
          }
          has_topo_data <- TRUE
        }
        
        if (has_topo_data) {
          saveWorkbook(wb_topo, topo_xlsx_path, overwrite = TRUE)
        }
        
        # 6.4 Distribution Viz
        if (!is.null(net_res$networks)) {
          dist_pdf_path <- file.path(out_dir, paste0(scenario$id, "_Edge_Distribution.pdf"))
          pdf(dist_pdf_path, width = 8, height = 6)
          
          thresh_to_plot <- if(!is.null(net_res$applied_threshold)) net_res$applied_threshold else edge_threshold_cfg
          
          if(sum(abs(net_res$networks$ctrl) > 0) > 0) {
            print(viz_plot_edge_density(net_res$networks$ctrl, 
                                        adj_mat = net_res$stability$ctrl,
                                        threshold = thresh_to_plot, 
                                        group_label = scenario$control_label))
            
            if (!is.null(net_res$raw_cor$ctrl)) {
              print(viz_plot_edge_density_overlay(pcor_mat = net_res$networks$ctrl, 
                                                  cor_mat = net_res$raw_cor$ctrl,
                                                  adj_mat = net_res$stability$ctrl,
                                                  threshold = thresh_to_plot, 
                                                  group_label = scenario$control_label))
            }
          }
          
          if(sum(abs(net_res$networks$case) > 0) > 0) {
            print(viz_plot_edge_density(net_res$networks$case, 
                                        adj_mat = net_res$stability$case,
                                        threshold = thresh_to_plot, 
                                        group_label = scenario$case_label))
            
            if (!is.null(net_res$raw_cor$case)) {
              print(viz_plot_edge_density_overlay(pcor_mat = net_res$networks$case, 
                                                  cor_mat = net_res$raw_cor$case,
                                                  adj_mat = net_res$stability$case,
                                                  threshold = thresh_to_plot, 
                                                  group_label = scenario$case_label))
            }
          }
          dev.off()
        }
        
        # 6.5 Network Core Structure Viz
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
        }
        
        # 6.6 Explicit Cytoscape Formatter
        message("      [Export] Generating Explicit Cytoscape Networks...")
        cyto_dir <- file.path(out_dir, "cytoscape_export")
        if (!dir.exists(cyto_dir)) dir.create(cyto_dir)
        
        if (!is.null(net_res$edges_table) && nrow(net_res$edges_table) > 0) {
          
          base_edges <- net_res$edges_table %>%
            dplyr::filter(Edge_Category != "Weak")
          
          if (nrow(base_edges) > 0) {
            
            col_pcor_ctrl <- paste0("Pcor_", scenario$control_label)
            col_pcor_case <- paste0("Pcor_", scenario$case_label)
            col_stab_ctrl <- paste0("Is_Stable_", scenario$control_label)
            col_stab_case <- paste0("Is_Stable_", scenario$case_label)
            
            # 6.6.1 Control Target Matrix
            ctrl_net <- base_edges %>%
              dplyr::select(
                Source = Node1,
                Target = Node2,
                Pcor = .data[[col_pcor_ctrl]],
                Is_Stable = .data[[col_stab_ctrl]],
                Edge_Category,
                Significant,
                P_Value
              ) %>%
              dplyr::mutate(Interaction = Edge_Category)
            
            readr::write_csv(ctrl_net, file.path(cyto_dir, paste0(scenario$id, "_", scenario$control_label, "_network.csv")))
            
            # 6.6.2 Case Target Matrix
            case_net <- base_edges %>%
              dplyr::select(
                Source = Node1,
                Target = Node2,
                Pcor = .data[[col_pcor_case]],
                Is_Stable = .data[[col_stab_case]],
                Edge_Category,
                Significant,
                P_Value
              ) %>%
              dplyr::mutate(Interaction = Edge_Category)
            
            readr::write_csv(case_net, file.path(cyto_dir, paste0(scenario$id, "_", scenario$case_label, "_network.csv")))
            
            # 6.6.3 Differential Matrix
            diff_net <- base_edges %>%
              dplyr::select(
                Source = Node1,
                Target = Node2,
                Weight = Diff_Score,
                Pcor_Control = .data[[col_pcor_ctrl]],
                Pcor_Case = .data[[col_pcor_case]],
                Edge_Category,
                Significant,
                P_Value
              ) %>%
              dplyr::mutate(Interaction = Edge_Category)
            
            readr::write_csv(diff_net, file.path(cyto_dir, paste0(scenario$id, "_diff_network.csv")))
            
            # 6.6.4 Global Node Attributes
            all_nodes <- unique(c(base_edges$Node1, base_edges$Node2))
            cyto_nodes <- data.frame(Node = all_nodes, stringsAsFactors = FALSE)
            
            if (!is.null(topo_ctrl)) {
              t_ctrl <- topo_ctrl %>% 
                dplyr::select(Node, Degree, Betweenness) %>%
                dplyr::rename(!!paste0("Degree_", scenario$control_label) := Degree,
                              !!paste0("Betweenness_", scenario$control_label) := Betweenness)
              cyto_nodes <- dplyr::left_join(cyto_nodes, t_ctrl, by = "Node")
            }
            
            if (!is.null(topo_case)) {
              t_case <- topo_case %>% 
                dplyr::select(Node, Degree, Betweenness) %>%
                dplyr::rename(!!paste0("Degree_", scenario$case_label) := Degree,
                              !!paste0("Betweenness_", scenario$case_label) := Betweenness)
              cyto_nodes <- dplyr::left_join(cyto_nodes, t_case, by = "Node")
            }
            
            if (!is.null(spls_drivers) && nrow(spls_drivers) > 0) {
              drv_summ <- spls_drivers %>%
                dplyr::group_by(Marker) %>%
                dplyr::slice_max(order_by = Importance, n = 1, with_ties = FALSE) %>%
                dplyr::ungroup() %>%
                dplyr::select(Marker, Importance, Direction) %>%
                dplyr::rename(PLSDA_Importance = Importance, PLSDA_Direction = Direction)
              
              cyto_nodes <- dplyr::left_join(cyto_nodes, drv_summ, by = c("Node" = "Marker"))
            }
            
            numeric_cols <- names(cyto_nodes)[sapply(cyto_nodes, is.numeric)]
            cyto_nodes[numeric_cols] <- lapply(cyto_nodes[numeric_cols], function(x) replace(x, is.na(x), 0))
            
            readr::write_csv(cyto_nodes, file.path(cyto_dir, paste0(scenario$id, "_node_attributes.csv")))
          }
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