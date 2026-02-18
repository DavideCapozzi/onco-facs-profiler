# src/04_network_meta_analysis.R
# ==============================================================================
# STEP 04: NETWORK META-ANALYSIS & DIFFERENTIAL OVERLAP
# Description: Compares Differential Networks across analysis scenarios using 
#              the Master Report. Generates Venn/UpSet plots.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(ComplexHeatmap) 
})

# Load Modules (Paths relative to project root)
source("R/utils_io.R")
source("R/modules_viz.R")

message("\n=== PIPELINE STEP 4: DIFFERENTIAL META-ANALYSIS ===")

# 1. Setup & Configuration
config <- load_config("config/global_params.yml")

# Define Paths - POINTING TO RESULTS_ANALYSIS WHERE MAIN.R SAVES THE REPORT
results_dir <- file.path(config$output_root, "results_analysis") 
master_report_path <- file.path(results_dir, "Multi_Scenario_Analysis_Report.xlsx")
out_dir <- file.path(results_dir, "Meta_Analysis")

# Ensure output directory exists
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Validate Input
if (!file.exists(master_report_path)) {
  stop(sprintf("[Fatal] Master Report not found at: %s", master_report_path))
}

message(sprintf("[Meta] Reading Master Report: %s", basename(master_report_path)))

# 2. Harvest Data
scenarios <- config$analysis_scenarios
diff_edges_list <- list()
scenario_colors <- c()
available_sheets <- getSheetNames(master_report_path)

for (scen in scenarios) {
  sheet_name <- substr(paste0(scen$id, "_Net"), 1, 31)
  label <- scen$id
  
  # Color Logic
  target_group <- scen$case_label
  assigned_color <- "grey50"
  if (!is.null(config$colors$groups[[target_group]])) {
    assigned_color <- config$colors$groups[[target_group]]
  } else if (!is.null(config$colors$cases)) {
    assigned_color <- config$colors$cases[[1]]
  }
  
  if (sheet_name %in% available_sheets) {
    df_net <- read.xlsx(master_report_path, sheet = sheet_name)
    
    if (nrow(df_net) > 0) {
      sig_edges <- data.frame()
      
      # --- ROBUST FILTERING LOGIC ---
      # 1. Edge Category (Best)
      if ("Edge_Category" %in% names(df_net)) {
        sig_edges <- df_net %>% filter(Edge_Category != "Weak")
        message(sprintf("   -> Scenario '%s': Filtering via Edge_Category.", label))
        
        # 2. Significant Boolean (Legacy)
      } else if ("Significant" %in% names(df_net)) {
        sig_edges <- df_net %>% filter(Significant == TRUE)
        message(sprintf("   -> Scenario '%s': Filtering via 'Significant' column.", label))
        
        # 3. P-Value Fallback
      } else {
        sig_edges <- df_net %>% filter(P_Value < 0.05)
        message(sprintf("   -> Scenario '%s': Filtering via P-Value < 0.05.", label))
      }
      
      if (nrow(sig_edges) > 0) {
        edge_ids <- apply(sig_edges, 1, function(r) {
          nodes <- sort(c(r["Node1"], r["Node2"]))
          paste(nodes, collapse = "~")
        })
        
        diff_edges_list[[label]] <- unique(edge_ids)
        scenario_colors[label] <- assigned_color
      }
    }
  }
}

# 3. Visualization & Export
if (length(diff_edges_list) >= 2) {
  
  # A. Plot
  pdf_path <- file.path(out_dir, "Differential_Networks_Overlap.pdf")
  pdf(pdf_path, width = 8, height = 8)
  tryCatch({
    viz_plot_differential_overlap(
      edge_list = diff_edges_list, 
      fill_colors = scenario_colors,
      title = "Differential Edges Overlap"
    )
  }, error = function(e) {
    plot.new(); text(0.5,0.5, "Plot Error")
  })
  dev.off()
  message(sprintf("   [Output] Plot saved: %s", basename(pdf_path)))
  
  # B. Intersection Table
  wb <- createWorkbook()
  addWorksheet(wb, "Summary")
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
      
      summary_df <- rbind(summary_df, data.frame(
        Intersection = set_label, Count = length(edges_in_comb), Sheet_Link = sheet_id
      ))
      
      addWorksheet(wb, sheet_id)
      details <- data.frame(Edge_ID = edges_in_comb) %>% 
        separate(Edge_ID, into = c("Node_A", "Node_B"), sep = "~")
      writeData(wb, sheet_id, details)
    }
  }
  writeData(wb, "Summary", summary_df)
  saveWorkbook(wb, file.path(out_dir, "Differential_Intersections_List.xlsx"), overwrite = TRUE)
  
} else {
  message("[Meta] Less than 2 scenarios with significant edges. Venn skipped.")
}

message("=== STEP 4 COMPLETE ===\n")