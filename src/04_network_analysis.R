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
# PHASE 1: COMPUTE UNIVERSAL BASELINES
# ==============================================================================
message("\n--- PHASE 1: GENERATING UNIVERSAL BASELINES ---")

# Collect unique macro-groups required across all scenarios
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
  
  # Inject Micro-Jitter
  var_thresh <- if(!is.null(config$stats$min_variance)) config$stats$min_variance else 1e-6
  vars <- apply(sub_mat, 2, var, na.rm = TRUE)
  low_vars <- vars < var_thresh | is.na(vars)
  
  if (sum(low_vars) > 0) {
    dropped_names <- names(vars)[low_vars]
    message(sprintf("   [%s] Injecting micro-jitter for %d flat features.", label, sum(low_vars)))
    set.seed(if(!is.null(config$stats$seed)) config$stats$seed else 123)
    for (col in dropped_names) {
      sub_mat[, col] <- sub_mat[, col] + rnorm(nrow(sub_mat), mean = 0, sd = 1e-8)
    }
  }
  
  baselines[[label]] <- compute_universal_baseline(
    mat = sub_mat,
    label = label,
    n_boot = if(!is.null(config$stats$n_boot)) config$stats$n_boot else 100,
    seed = if(!is.null(config$stats$seed)) config$stats$seed else 123,
    n_cores = if(!is.null(config$stats$n_cores) && config$stats$n_cores != "auto") config$stats$n_cores else 1,
    threshold_type = if(!is.null(config$stats$network_threshold_type)) config$stats$network_threshold_type else "percentile",
    threshold_value = if(!is.null(config$stats$network_edge_threshold)) config$stats$network_edge_threshold else 0.15
  )
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
  
  # Run Overlay
  net_res <- compute_differential_overlay(
    base_ctrl = base_ctrl,
    base_case = base_case,
    n_perm = if(!is.null(config$stats$n_perm)) config$stats$n_perm else 1000,
    seed = if(!is.null(config$stats$seed)) config$stats$seed else 123,
    n_cores = if(!is.null(config$stats$n_cores) && config$stats$n_cores != "auto") config$stats$n_cores else 1,
    pvalue_thresh = if(!is.null(config$stats$pvalue_threshold)) config$stats$pvalue_threshold else 0.05
  )
  
  if (!is.null(net_res) && nrow(net_res$edges_table) > 0) {
    
    # Extract sPLS-DA drivers from Master Report for Hub-Driver mapping
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
    
    # Save Plot PDFs
    pdf_hd_path <- file.path(scen_dir, paste0(scen$id, "_Hub_Driver_Report.pdf"))
    if (!is.null(spls_drivers) && nrow(spls_drivers) > 0 && (!is.null(topo_ctrl) || !is.null(topo_case))) {
      pdf(pdf_hd_path, width = 11, height = 8)
      if(!is.null(topo_ctrl)) print(plot_hub_driver_quadrant(integrate_hub_drivers(spls_drivers, topo_ctrl, "Degree"), "Degree", paste("Network:", scen$control_label)))
      if(!is.null(topo_case)) print(plot_hub_driver_quadrant(integrate_hub_drivers(spls_drivers, topo_case, "Degree"), "Degree", paste("Network:", scen$case_label)))
      dev.off()
    }
    
    pdf_net_path <- file.path(scen_dir, paste0(scen$id, "_Plot_Networks.pdf"))
    pdf(pdf_net_path, width = 12, height = 6)
    if(!is.null(topo_ctrl)) print(plot_network_structure(net_res$adj_final$ctrl, net_res$networks$ctrl, title = paste("Control:", scen$control_label), min_cor = 0))
    if(!is.null(topo_case)) print(plot_network_structure(net_res$adj_final$case, net_res$networks$case, title = paste("Case:", scen$case_label), min_cor = 0))
    dev.off()
    
    # Cytoscape EXPORT
    cyto_dir <- file.path(scen_dir, "cytoscape_export")
    if (!dir.exists(cyto_dir)) dir.create(cyto_dir)
    readr::write_csv(net_res$edges_table %>% filter(Significant == TRUE), file.path(cyto_dir, paste0(scen$id, "_diff_network.csv")))
    
    all_nodes <- data.frame(Node = unique(c(net_res$edges_table$Node1, net_res$edges_table$Node2)), stringsAsFactors=F)
    if(!is.null(topo_ctrl)) all_nodes <- left_join(all_nodes, topo_ctrl %>% select(Node, Degree) %>% rename(!!paste0("Degree_", scen$control_label) := Degree), by="Node")
    if(!is.null(topo_case)) all_nodes <- left_join(all_nodes, topo_case %>% select(Node, Degree) %>% rename(!!paste0("Degree_", scen$case_label) := Degree), by="Node")
    readr::write_csv(all_nodes, file.path(cyto_dir, paste0(scen$id, "_node_attributes.csv")))
    
    # Append to Master Report
    net_sheet <- substr(paste0(scen$id, "_Net"), 1, 31)
    if(!net_sheet %in% names(wb_master)) addWorksheet(wb_master, net_sheet)
    writeData(wb_master, net_sheet, net_res$edges_table)
    
    sig_edges <- net_res$edges_table %>% filter(Edge_Category != "Weak" & Significant == TRUE) %>% arrange(P_Value)
    if(nrow(sig_edges) > 0) {
      sig_sheet <- substr(paste0(scen$id, "_SigNet"), 1, 31)
      if(!sig_sheet %in% names(wb_master)) addWorksheet(wb_master, sig_sheet)
      writeData(wb_master, sig_sheet, sig_edges)
      
      # Prepare for Meta-Analysis
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