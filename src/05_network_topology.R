# src/05_network_topology.R
# ==============================================================================
# STEP 05: NETWORK TOPOLOGY & CYTOSCAPE EXPORT
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(igraph)
  library(ggraph)
})

source("R/utils_io.R")          
source("R/modules_network.R")   
source("R/modules_viz.R")       

message("\n=== PIPELINE STEP 5: NETWORK ANALYSIS & EXPORT ===")

# 1. Load Config & Step 4 Results
config <- load_config("config/global_params.yml")
input_file <- file.path(config$output_root, "04_network_inference", "inference_results.rds")

if (!file.exists(input_file)) stop("Step 04 output not found. Run src/04_network_inference.R first.")

INFERENCE <- readRDS(input_file)
res_ctrl <- INFERENCE$ctrl_network
res_case <- INFERENCE$case_network

# 2. Prepare Output Directory
out_dir <- file.path(config$output_root, "05_network_topology")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 3. Network Topology Analysis (Excel)
# ------------------------------------------------------------------------------
message("[Analysis] Calculating topology metrics...")

# Function name get_topology_metrics is unchanged in modules_network.R
topo_ctrl <- get_topology_metrics(res_ctrl$adj, res_ctrl$weights)
topo_case <- get_topology_metrics(res_case$adj, res_case$weights)

# Add group label for potential merging later
if(!is.null(topo_ctrl)) topo_ctrl$Group <- "Control"
if(!is.null(topo_case)) topo_case$Group <- "Case"

# Setup Excel Workbook
wb <- createWorkbook()
addWorksheet(wb, "Control_Topology")
addWorksheet(wb, "Case_Topology")

if(!is.null(topo_ctrl)) writeData(wb, "Control_Topology", topo_ctrl)
if(!is.null(topo_case)) writeData(wb, "Case_Topology", topo_case)

topo_file <- file.path(out_dir, "Network_Topology_Metrics.xlsx")
saveWorkbook(wb, topo_file, overwrite = TRUE)
message(sprintf("   -> Topology metrics saved to: %s", topo_file))

# 4. Visualization (PDF)
# ------------------------------------------------------------------------------
message("[Viz] Generating network plots...")

viz_thresh <- if (!is.null(config$viz$min_edge_weight)) config$viz$min_edge_weight else 0

pdf_path <- file.path(out_dir, "Network_Structures.pdf")
pdf(pdf_path, width = 10, height = 8)

# Plot Control
if (sum(res_ctrl$adj) > 0) {
  p_ctrl <- plot_network_structure(
    res_ctrl$adj, 
    res_ctrl$weights, 
    title = paste("Control Network:", config$control_group),
    min_cor = viz_thresh
  )
  print(p_ctrl)
} else {
  plot.new()
  text(0.5, 0.5, "No stable edges found in Control group.")
}

# Plot Case
if (sum(res_case$adj) > 0) {
  p_case <- plot_network_structure(
    res_case$adj, 
    res_case$weights, 
    title = paste("Case Network:", paste(config$case_groups, collapse="+")),
    min_cor = viz_thresh
  )
  print(p_case)
} else {
  plot.new()
  text(0.5, 0.5, "No stable edges found in Case group.")
}

dev.off()
message(sprintf("   -> Network plots saved to: %s", pdf_path))

# 5. Cytoscape Export (CSV/TSV)
# ------------------------------------------------------------------------------
message("[Export] Generating Cytoscape input files...")

cyto_dir <- file.path(out_dir, "cytoscape_inputs")
if (!dir.exists(cyto_dir)) dir.create(cyto_dir, recursive = TRUE)

# Format Edges - UPDATED FUNCTION NAME
edges_ctrl <- export_cytoscape_edges(res_ctrl$adj, res_ctrl$weights)
edges_case <- export_cytoscape_edges(res_case$adj, res_case$weights)

# Write to CSV
write.csv(edges_ctrl, file.path(cyto_dir, "edges_control.csv"), row.names = FALSE)
write.csv(edges_case, file.path(cyto_dir, "edges_case.csv"), row.names = FALSE)

# Also create a Node Attributes file
if(!is.null(topo_ctrl)) {
  write.csv(topo_ctrl, file.path(cyto_dir, "nodes_attributes_control.csv"), row.names = FALSE)
}
if(!is.null(topo_case)) {
  write.csv(topo_case, file.path(cyto_dir, "nodes_attributes_case.csv"), row.names = FALSE)
}

message(sprintf("   -> Cytoscape files saved to: %s", cyto_dir))
message("=== STEP 5 COMPLETE ===\n")