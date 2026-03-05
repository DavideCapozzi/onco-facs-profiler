# src/04_network_analysis.R
# ==============================================================================
# STEP 04: NETWORK ANALYSIS & META-ANALYSIS
# Description: Orchestrates Universal Baselines, computes Differential Overlays
#              for all scenarios, saves topologies, appends to Master Report, 
#              and generates Venn/UpSet plots.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(ComplexHeatmap)
})

source("R/utils_io.R")
source("R/modules_viz.R")
source("R/modules_network.R")

message("\n=== PIPELINE STEP 4: NETWORK ANALYSIS & DIFFERENTIAL META-ANALYSIS ===")

config <- load_config("config/global_params.yml")
data_file <- file.path(config$output_root, "01_data_processing", "data_processed.rds")
if (!file.exists(data_file)) stop("[Fatal] Step 01 output not found.")
DATA <- readRDS(data_file)

results_dir <- file.path(config$output_root, "results_analysis")
master_report_path <- file.path(results_dir, "Multi_Scenario_Analysis_Report.xlsx")
if (!file.exists(master_report_path)) stop("[Fatal] Master Report from Step 03 not found.")

# ==============================================================================
# PHASE 1: COMPUTE UNIVERSAL BASELINES & CENTRALIZED EXPORT
# ==============================================================================
message("\n--- PHASE 1: GENERATING UNIVERSAL BASELINES ---")

# Setup centralized Cytoscape folder for absolute baseline networks
cyto_base_dir <- file.path(results_dir, "cytoscape_baselines")
if (!dir.exists(cyto_base_dir)) dir.create(cyto_base_dir, recursive = TRUE)

macro_groups <- list()
for (scen in config$analysis_scenarios) {
  macro_groups[[scen$control_label]] <- unique(c(macro_groups[[scen$control_label]], scen$control_groups))
  macro_groups[[scen$case_label]]    <- unique(c(macro_groups[[scen$case_label]], scen$case_groups))
}

baselines <- list()
strat_col <- config$stratification$column

for (label in names(macro_groups)) {
  groups_in_macro <- macro_groups[[label]]
  
  sub_meta <- DATA$metadata %>% dplyr::filter(.data[[strat_col]] %in% groups_in_macro)
  sub_mat <- DATA$hybrid_data_z[sub_meta$Patient_ID, DATA$hybrid_markers, drop=FALSE]
  
  if (nrow(sub_mat) < if(!is.null(config$stats$min_network_n)) config$stats$min_network_n else 5) {
    message(sprintf("   [Skip] Macro-Group '%s' skipped (insufficient samples).", label))
    next
  }
  
  # Zero-Variance Detection & Subset
  var_thresh <- if(!is.null(config$stats$min_variance)) config$stats$min_variance else 1e-6
  vars <- apply(sub_mat, 2, var, na.rm = TRUE)
  low_vars <- vars < var_thresh | is.na(vars)
  
  valid_cols <- names(vars)[!low_vars]
  dropped_cols <- names(vars)[low_vars]
  
  if (length(dropped_cols) > 0) {
    message(sprintf("   [%s] Detected %d flat features. Isolating topographically...", label, length(dropped_cols)))
  }
  
  sub_mat_valid <- sub_mat[, valid_cols, drop = FALSE]
  
  if (ncol(sub_mat_valid) < 2) {
    message(sprintf("   [Skip] Macro-Group '%s' skipped (insufficient varying features).", label))
    next
  }
  
  # Compute Baseline ONLY on valid varying features
  base_obj_raw <- compute_universal_baseline(
    mat = sub_mat_valid,
    label = label,
    n_boot = if(!is.null(config$stats$n_boot)) config$stats$n_boot else 100,
    seed = if(!is.null(config$stats$seed)) config$stats$seed else 123,
    n_cores = if(!is.null(config$stats$n_cores) && config$stats$n_cores != "auto") config$stats$n_cores else 1,
    threshold_type = if(!is.null(config$stats$network_threshold_type)) config$stats$network_threshold_type else "percentile",
    threshold_value = if(!is.null(config$stats$network_edge_threshold)) config$stats$network_edge_threshold else 0.15
  )
  
  # Topological Padding: Rebuild full matrices with zeros for dropped features
  all_features <- colnames(sub_mat)
  p <- length(all_features)
  
  # Helper function to pad matrices
  pad_matrix <- function(small_mat, full_names, default_val = 0) {
    big_mat <- matrix(default_val, nrow = p, ncol = p, dimnames = list(full_names, full_names))
    if (!is.null(small_mat) && nrow(small_mat) > 0) {
      valid_n <- rownames(small_mat)
      big_mat[valid_n, valid_n] <- small_mat
    }
    return(big_mat)
  }
  
  base_obj <- list(
    label = base_obj_raw$label,
    mat = sub_mat, # Keep original matrix for reference
    pcor = pad_matrix(base_obj_raw$pcor, all_features, 0),
    raw_cor = pad_matrix(base_obj_raw$raw_cor, all_features, 0),
    stability = pad_matrix(base_obj_raw$stability, all_features, 0),
    adj_final = pad_matrix(base_obj_raw$adj_final, all_features, FALSE),
    applied_threshold = base_obj_raw$applied_threshold
  )
  # End of Topological Padding substitution
  baselines[[label]] <- base_obj
  
  # Centralized Cytoscape Baseline Export (Only valid edges, inherently drops "weak" concepts)
  adj <- base_obj$adj_final
  edges_idx <- which(upper.tri(adj) & adj, arr.ind = TRUE)
  
  if (nrow(edges_idx) > 0) {
    nodes <- colnames(adj)
    base_edges <- data.frame(
      Source = nodes[edges_idx[,1]],
      Target = nodes[edges_idx[,2]],
      Pcor = base_obj$pcor[edges_idx],
      Is_Stable = TRUE,
      stringsAsFactors = FALSE
    )
    readr::write_csv(base_edges, file.path(cyto_base_dir, paste0(label, "_baseline_network.csv")))
    
    topo_base <- calculate_node_topology(adj)
    if (!is.null(topo_base)) {
      readr::write_csv(topo_base, file.path(cyto_base_dir, paste0(label, "_node_attributes.csv")))
    }
  }
}

# ==============================================================================
# PHASE 2: DIFFERENTIAL OVERLAYS & EXPORTS
# ==============================================================================
message("\n--- PHASE 2: DIFFERENTIAL OVERLAYS & CYTOSCAPE EXPORT ---")

wb_master <- loadWorkbook(master_report_path)
diff_edges_list <- list()
scenario_colors <- c()

for (scen in config$analysis_scenarios) {
  message(sprintf("\n   [Scenario] Processing %s...", scen$id))
  
  base_ctrl <- baselines[[scen$control_label]]
  base_case <- baselines[[scen$case_label]]
  
  if (is.null(base_ctrl) || is.null(base_case)) {
    message(sprintf("      [Skip] Missing baseline objects for %s.", scen$id))
    next
  }
  
  scen_dir <- file.path(results_dir, scen$id)
  if (!dir.exists(scen_dir)) dir.create(scen_dir, recursive = TRUE)
  
  net_res <- compute_differential_overlay(
    base_ctrl = base_ctrl,
    base_case = base_case,
    n_perm = if(!is.null(config$stats$n_perm)) config$stats$n_perm else 1000,
    seed = if(!is.null(config$stats$seed)) config$stats$seed else 123,
    n_cores = if(!is.null(config$stats$n_cores) && config$stats$n_cores != "auto") config$stats$n_cores else 1,
    pvalue_thresh = if(!is.null(config$stats$pvalue_threshold)) config$stats$pvalue_threshold else 0.05
  )
  
  if (!is.null(net_res) && nrow(net_res$edges_table) > 0) {
    
    # Restored: Save Raw RDS Object
    rds_path <- file.path(scen_dir, paste0(scen$id, "_Differential_Network.rds"))
    saveRDS(net_res, rds_path)
    message(sprintf("      [Output] Network RDS saved: %s", basename(rds_path)))
    
    # Extract sPLS-DA drivers from Master Report
    drv_sheet <- substr(paste0(scen$id, "_Drv"), 1, 31)
    spls_drivers <- NULL
    if (drv_sheet %in% names(wb_master)) {
      spls_drivers <- read.xlsx(master_report_path, sheet = drv_sheet)
    }
    
    # Topologies & Rewiring Export
    topo_xlsx_path <- file.path(scen_dir, paste0(scen$id, "_Topology_Metrics.xlsx"))
    wb_topo <- createWorkbook()
    
    topo_ctrl <- if (sum(net_res$adj_final$ctrl)>0) calculate_node_topology(net_res$adj_final$ctrl) else NULL
    topo_case <- if (sum(net_res$adj_final$case)>0) calculate_node_topology(net_res$adj_final$case) else NULL
    
    if (!is.null(topo_ctrl)) {
      sh_name <- substr(paste0("Topology_", scen$control_label), 1, 31)
      addWorksheet(wb_topo, sh_name); writeData(wb_topo, sh_name, topo_ctrl)
    }
    if (!is.null(topo_case)) {
      sh_name <- substr(paste0("Topology_", scen$case_label), 1, 31)
      addWorksheet(wb_topo, sh_name); writeData(wb_topo, sh_name, topo_case)
    }
    
    if (!is.null(topo_ctrl) && !is.null(topo_case)) {
      rewiring_df <- calculate_jaccard_rewiring(net_res$adj_final$ctrl, net_res$adj_final$case)
      addWorksheet(wb_topo, "Rewiring_Analysis"); writeData(wb_topo, "Rewiring_Analysis", rewiring_df)
    }
    
    addWorksheet(wb_topo, "Differential_Edges"); writeData(wb_topo, "Differential_Edges", net_res$edges_table)
    saveWorkbook(wb_topo, topo_xlsx_path, overwrite = TRUE)
    
    # Restored: Save Hub-Driver PDF with exact original name
    pdf_hd_path <- file.path(scen_dir, paste0(scen$id, "_Hub_Driver_Report_4Page.pdf"))
    if (!is.null(spls_drivers) && nrow(spls_drivers) > 0 && (!is.null(topo_ctrl) || !is.null(topo_case))) {
      pdf(pdf_hd_path, width = 11, height = 8)
      if(!is.null(topo_ctrl)) print(plot_hub_driver_quadrant(integrate_hub_drivers(spls_drivers, topo_ctrl, "Degree"), "Degree", paste("\nNetwork:", scen$control_label)))
      if(!is.null(topo_case)) print(plot_hub_driver_quadrant(integrate_hub_drivers(spls_drivers, topo_case, "Degree"), "Degree", paste("\nNetwork:", scen$case_label)))
      if(!is.null(topo_ctrl)) print(plot_hub_driver_quadrant(integrate_hub_drivers(spls_drivers, topo_ctrl, "Betweenness"), "Betweenness", paste("\nNetwork:", scen$control_label)))
      if(!is.null(topo_case)) print(plot_hub_driver_quadrant(integrate_hub_drivers(spls_drivers, topo_case, "Betweenness"), "Betweenness", paste("\nNetwork:", scen$case_label)))
      dev.off()
    }
    
    # Restored: Edge Distribution Visualization
    dist_pdf_path <- file.path(scen_dir, paste0(scen$id, "_Edge_Distribution.pdf"))
    pdf(dist_pdf_path, width = 8, height = 6)
    thresh_to_plot <- net_res$applied_threshold
    
    if(sum(abs(net_res$networks$ctrl) > 0) > 0) {
      print(viz_plot_edge_density(net_res$networks$ctrl, adj_mat = net_res$adj_final$ctrl, threshold = thresh_to_plot, group_label = scen$control_label))
      if (!is.null(net_res$raw_cor$ctrl)) {
        print(viz_plot_edge_density_overlay(pcor_mat = net_res$networks$ctrl, cor_mat = net_res$raw_cor$ctrl, adj_mat = net_res$adj_final$ctrl, threshold = thresh_to_plot, group_label = scen$control_label))
      }
    }
    if(sum(abs(net_res$networks$case) > 0) > 0) {
      print(viz_plot_edge_density(net_res$networks$case, adj_mat = net_res$adj_final$case, threshold = thresh_to_plot, group_label = scen$case_label))
      if (!is.null(net_res$raw_cor$case)) {
        print(viz_plot_edge_density_overlay(pcor_mat = net_res$networks$case, cor_mat = net_res$raw_cor$case, adj_mat = net_res$adj_final$case, threshold = thresh_to_plot, group_label = scen$case_label))
      }
    }
    dev.off()
    
    pdf_net_path <- file.path(scen_dir, paste0(scen$id, "_Plot_Networks.pdf"))
    pdf(pdf_net_path, width = 12, height = 6)
    if(!is.null(topo_ctrl)) print(plot_network_structure(net_res$adj_final$ctrl, net_res$networks$ctrl, title = paste("Control:", scen$control_label), min_cor = 0))
    if(!is.null(topo_case)) print(plot_network_structure(net_res$adj_final$case, net_res$networks$case, title = paste("Case:", scen$case_label), min_cor = 0))
    dev.off()
    
    # Cytoscape EXPORT (Differential Only)
    cyto_dir <- file.path(scen_dir, "cytoscape_export")
    if (!dir.exists(cyto_dir)) dir.create(cyto_dir)
    
    diff_edges <- net_res$edges_table %>% filter(Significant == TRUE & Edge_Category != "Weak")
    if (nrow(diff_edges) > 0) {
      
      # 1. Export Edges (Clean format: only Differential Weight and exact Interaction)
      diff_net <- diff_edges %>%
        dplyr::select(Source = Node1, Target = Node2, Weight = Diff_Score,
                      Significant, P_Value, Edge_Category) %>%
        dplyr::rename(Interaction = Edge_Category)
      readr::write_csv(diff_net, file.path(cyto_dir, paste0(scen$id, "_diff_network.csv")))
      
      # 2. Export Nodes (Calculate topology strictly on the differential network structure)
      diff_nodes <- unique(c(diff_edges$Node1, diff_edges$Node2))
      
      # Build unweighted adjacency matrix for the differential sub-network
      adj_diff <- matrix(0, nrow = length(diff_nodes), ncol = length(diff_nodes), 
                         dimnames = list(diff_nodes, diff_nodes))
      
      for(k in 1:nrow(diff_edges)) {
        n1 <- diff_edges$Node1[k]
        n2 <- diff_edges$Node2[k]
        adj_diff[n1, n2] <- 1
        adj_diff[n2, n1] <- 1
      }
      
      # Calculate topology features for differential nodes
      topo_diff <- calculate_node_topology(adj_diff)
      
      if (!is.null(topo_diff)) {
        all_nodes <- topo_diff %>% 
          dplyr::select(Node, Degree, Betweenness, Closeness, Eigen_Centrality)
        
        # Merge with sPLS-DA drivers if available
        if (!is.null(spls_drivers) && nrow(spls_drivers) > 0) {
          drv_summ <- spls_drivers %>% 
            dplyr::group_by(Marker) %>% 
            dplyr::slice_max(order_by = Importance, n = 1, with_ties = FALSE) %>%
            dplyr::ungroup() %>% 
            dplyr::select(Marker, Importance, Direction) %>% 
            dplyr::rename(PLSDA_Importance = Importance, PLSDA_Direction = Direction)
          
          all_nodes <- dplyr::left_join(all_nodes, drv_summ, by = c("Node" = "Marker"))
        }
        
        # Safely handle NAs for nodes that are in the network but are not PLS-DA drivers
        numeric_cols <- names(all_nodes)[sapply(all_nodes, is.numeric)]
        all_nodes[numeric_cols] <- lapply(all_nodes[numeric_cols], function(x) replace(x, is.na(x), 0))
        
        readr::write_csv(all_nodes, file.path(cyto_dir, paste0(scen$id, "_node_attributes.csv")))
      }
    }
    
    # Append to Master Report
    net_sheet <- substr(paste0(scen$id, "_Net"), 1, 31)
    if(!net_sheet %in% names(wb_master)) addWorksheet(wb_master, net_sheet)
    writeData(wb_master, net_sheet, net_res$edges_table)
    
    sig_edges <- net_res$edges_table %>% filter(Edge_Category != "Weak" & Significant == TRUE) %>% arrange(P_Value)
    if(nrow(sig_edges) > 0) {
      sig_sheet <- substr(paste0(scen$id, "_SigNet"), 1, 31)
      if(!sig_sheet %in% names(wb_master)) addWorksheet(wb_master, sig_sheet)
      writeData(wb_master, sig_sheet, sig_edges)
      
      edge_ids <- apply(sig_edges, 1, function(r) paste(sort(c(r["Node1"], r["Node2"])), collapse="~"))
      diff_edges_list[[scen$id]] <- unique(edge_ids)
      
      tgt_grp <- scen$case_label
      scenario_colors[scen$id] <- if(!is.null(config$colors$groups[[tgt_grp]])) config$colors$groups[[tgt_grp]] else "grey50"
    }
  }
}

saveWorkbook(wb_master, master_report_path, overwrite = TRUE)
message("   [Output] Network sheets appended to Master Report.")

# ==============================================================================
# PHASE 3: META-ANALYSIS (Overlap)
# ==============================================================================
message("\n--- PHASE 3: DIFFERENTIAL OVERLAP META-ANALYSIS ---")

meta_dir <- file.path(results_dir, "Meta_Analysis")
if (!dir.exists(meta_dir)) dir.create(meta_dir, recursive = TRUE)

if (length(diff_edges_list) >= 2) {
  pdf_path <- file.path(meta_dir, "Differential_Networks_Overlap.pdf")
  pdf(pdf_path, width = 8, height = 8)
  tryCatch({
    viz_plot_differential_overlap(edge_list = diff_edges_list, fill_colors = scenario_colors, title = "Differential Edges Overlap")
  }, error = function(e) { plot.new(); text(0.5,0.5, "Plot Error") })
  dev.off()
  
  wb_meta <- createWorkbook()
  addWorksheet(wb_meta, "Summary")
  m_comb <- ComplexHeatmap::make_comb_mat(diff_edges_list)
  comb_names <- ComplexHeatmap::comb_name(m_comb)
  summary_df <- data.frame()
  
  for (nm in comb_names) {
    edges_in_comb <- ComplexHeatmap::extract_comb(m_comb, nm)
    if (length(edges_in_comb) > 0) {
      set_indices <- as.numeric(strsplit(nm, "")[[1]])
      involved_sets <- names(diff_edges_list)[which(set_indices == 1)]
      set_label <- paste(involved_sets, collapse = " & ")
      sheet_id <- paste0("Int_", nrow(summary_df) + 1)
      
      summary_df <- rbind(summary_df, data.frame(Intersection = set_label, Count = length(edges_in_comb), Sheet_Link = sheet_id))
      addWorksheet(wb_meta, sheet_id)
      details <- data.frame(Edge_ID = edges_in_comb) %>% separate(Edge_ID, into = c("Node_A", "Node_B"), sep = "~")
      writeData(wb_meta, sheet_id, details)
    }
  }
  writeData(wb_meta, "Summary", summary_df)
  saveWorkbook(wb_meta, file.path(meta_dir, "Differential_Intersections_List.xlsx"), overwrite = TRUE)
  message(sprintf("   [Output] Meta-Analysis exported to %s", basename(meta_dir)))
} else {
  message("   [Skip] Less than 2 scenarios with significant edges. Meta-Analysis skipped.")
}

message("\n=== STEP 4 COMPLETE ===\n")