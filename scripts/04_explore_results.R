# scripts/04_explore_results.R
# ==============================================================================
# STEP 04: RESULT EXPLORATION & VISUALIZATION
# ==============================================================================

library(tidyverse)
library(igraph)
library(yaml)

# 1. Load Configuration & Data
# ------------------------------------------------------------------------------
config <- read_yaml("config/global_params.yml")
output_dir <- file.path(config$output_root, "03_Networks")

rds_path <- file.path(output_dir, "inference_results.rds")
csv_path <- file.path(output_dir, "Differential_Network_Stats.csv")

if(!file.exists(rds_path)) stop("Results not found. Run Step 03 first.")

# Load the full R object and the stats table
final_res <- readRDS(rds_path)
edge_stats <- read.csv(csv_path)

message("=== RESULTS SUMMARY ===")

# 2. Statistical Overview
# ------------------------------------------------------------------------------
# Filter for significant edges (FDR < threshold)
sig_edges <- edge_stats %>% 
  filter(Significant == TRUE)

cat(sprintf("Total Edges Tested: %d\n", nrow(edge_stats)))
cat(sprintf("Significant Differential Edges (FDR < %.2f): %d\n", 
            config$stats$fdr_threshold, nrow(sig_edges)))

if(nrow(sig_edges) > 0) {
  # Top 5 most significant changes
  message("\nTop 5 Most Significant Differences:")
  print(head(sig_edges %>% arrange(Q_Value), 5))
}

# 3. Network Visualization (Differential Network)
# ------------------------------------------------------------------------------
# We will plot ONLY the edges that are significantly different.
# Edge Color: Red (Higher in Case), Blue (Higher in Control)
# Edge Width: Proportional to the magnitude of the difference

if(nrow(sig_edges) > 0) {
  
  # Create graph from data frame
  g <- graph_from_data_frame(sig_edges, directed = FALSE)
  
  # Define plotting attributes based on data
  # Calculate raw difference for coloring
  # (Using the observed weights we saved in the CSV)
  diff_vals <- sig_edges$Weight_Case_Obs - sig_edges$Weight_Ctrl_Obs
  
  # Color: Red if Case > Control, Blue otherwise
  E(g)$color <- ifelse(diff_vals > 0, "#D62728", "#1F77B4") # Red / Blue
  
  # Width: based on absolute difference
  E(g)$width <- abs(diff_vals) * 5
  
  # Node styling
  V(g)$color <- "white"
  V(g)$frame.color <- "black"
  V(g)$label.color <- "black"
  V(g)$size <- 15
  
  # Plot
  pdf(file.path(output_dir, "Differential_Network_Plot.pdf"), width = 10, height = 10)
  plot(g, layout = layout_with_fr,
       main = paste("Differential Network (FDR <", config$stats$fdr_threshold, ")"),
       sub = "Red: Stronger in Case | Blue: Stronger in Control")
  legend("bottomright", legend=c("Up in Case", "Up in Control"), 
         col=c("#D62728", "#1F77B4"), lwd=3, bty="n")
  dev.off()
  
  message(sprintf("\n[Output] Plot saved to: %s", file.path(output_dir, "Differential_Network_Plot.pdf")))
  
} else {
  message("\n[Info] No significant edges to plot.")
}

# 4. Accessing Specific Matrices (Examples)
# ------------------------------------------------------------------------------
# If you need the raw partial correlation matrices for other heatmaps:
mat_ctrl <- final_res$obs_matrices$ctrl
mat_case <- final_res$obs_matrices$case

# Example: Check connection between two specific genes/markers
# nodeA <- "Marker1" # Replace with real names
# nodeB <- "Marker2"
# print(mat_ctrl[nodeA, nodeB])
