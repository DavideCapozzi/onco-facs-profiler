# main.R
# ==============================================================================
# MAIN PIPELINE ORCHESTRATOR
# Description: End-to-end execution of the Immunological Analysis Pipeline.
#              Orchestrates ETL, Visualization, and Multi-Scenario Analysis.
# ==============================================================================

# 1. Environment Setup
suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(here)
  library(openxlsx)
})

# Load Modules
list.files("R", pattern = "\\.R$", full.names = TRUE) %>% walk(source)

# Load Configuration
config <- load_config("config/global_params.yml")
message(sprintf("[Main] Project: %s | Output: %s", 
                config$project_name, config$output_root))

# Create Output Root
if (!dir.exists(config$output_root)) dir.create(config$output_root, recursive = TRUE)


# ==============================================================================
# STEP 01 & 02: ETL & VISUALIZATION (Core Data Processing)
# ==============================================================================
message("\n--- PHASE 1: DATA PROCESSING & VIZ ---")
source("src/01_data_processing.R")
source("src/02_visualization.R")


# ==============================================================================
# STEP 03: ANALYTICAL WORKFLOWS (Multi-Scenario)
# ==============================================================================
message("\n--- PHASE 2: COMPARATIVE ANALYSIS (SCENARIOS) ---")

# Load processed data
data_file <- file.path(config$output_root, "01_data_processing", "data_processed.rds")
if (!file.exists(data_file)) stop("CRITICAL: Processed data not found.")
processed_data <- readRDS(data_file)

# Define results directory
results_root <- file.path(config$output_root, "results_analysis")
if (!dir.exists(results_root)) dir.create(results_root, recursive = TRUE)

if (is.null(config$analysis_scenarios) || length(config$analysis_scenarios) == 0) {
  warning("[Main] No 'analysis_scenarios' found in config. Skipping Phase 2.")
} else {
  
  # 1. Execute Workflows
  # ----------------------------------------------------------------------------
  all_results <- list()
  
  for (scenario in config$analysis_scenarios) {
    tryCatch({
      # run_comparative_workflow returns a list object, does not write tables
      res <- run_comparative_workflow(
        data_list = processed_data, 
        scenario = scenario, 
        config = config, 
        output_root = results_root
      )
      all_results[[scenario$id]] <- res
      
    }, error = function(e) {
      message(sprintf("!!! [Error] Scenario '%s' failed: %s", scenario$id, e$message))
    })
  }
  
  # 2. Consolidated Reporting (Excel)
  # ----------------------------------------------------------------------------
  if (length(all_results) > 0) {
    message("\n[Main] Generating Consolidated Multi-Scenario Report...")
    
    wb_master <- createWorkbook()
    
    # A. Global Drivers Summary (Concatenated)
    combined_drivers <- map_dfr(all_results, ~ .x$drivers, .id = "Scenario")
    if (nrow(combined_drivers) > 0) {
      addWorksheet(wb_master, "All_Drivers_Summary")
      writeData(wb_master, "All_Drivers_Summary", combined_drivers)
    }
    
    # B. Per-Scenario Detailed Sheets
    for (scen_id in names(all_results)) {
      res <- all_results[[scen_id]]
      
      # 1. PERMANOVA
      if (!is.null(res$permanova)) {
        sheet_name <- substr(paste0(scen_id, "_Perm"), 1, 31)
        addWorksheet(wb_master, sheet_name)
        writeData(wb_master, sheet_name, res$permanova, rowNames = TRUE)
      }
      
      # 2. sPLS-DA Drivers
      if (!is.null(res$drivers) && nrow(res$drivers) > 0) {
        sheet_name <- substr(paste0(scen_id, "_Drivers"), 1, 31)
        addWorksheet(wb_master, sheet_name)
        writeData(wb_master, sheet_name, res$drivers)
      }
      
      # 3. Differential Network Edges
      if (!is.null(res$network$edges_table) && nrow(res$network$edges_table) > 0) {
        sheet_name <- substr(paste0(scen_id, "_NetDiff"), 1, 31)
        addWorksheet(wb_master, sheet_name)
        writeData(wb_master, sheet_name, res$network$edges_table)
      }
    }
    
    # Save Master File
    out_file <- file.path(results_root, "Multi_Scenario_Analysis_Report.xlsx")
    saveWorkbook(wb_master, out_file, overwrite = TRUE)
    message(sprintf("       -> Saved Master Report: %s", out_file))
  }
}

message("\n=== PIPELINE COMPLETED SUCCESSFULLY ===")