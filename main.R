# main.R
# ==============================================================================
# MAIN PIPELINE ORCHESTRATOR
# Description: End-to-end execution with Real-Time TEE logging & Scenarios.
# ==============================================================================

# Clean environment
rm(list = ls())
graphics.off()

# 1. Setup Environment & Configuration
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(here)
  library(openxlsx)
})

# Disable ANSI colors for clean log files
options(crayon.enabled = FALSE)

# Load Config
config_path <- "config/global_params.yml"
if (!file.exists(config_path)) stop("Config file not found!")
config <- yaml::read_yaml(config_path)

# 2. Logging Setup
# ------------------------------------------------------------------------------
# Define Output Root
out_root <- config$output_root
if (!is.null(config$project_name) && config$project_name != "") {
  out_root <- paste0(out_root, "_", config$project_name)
}
if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)

# Define Log File
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(out_root, paste0("pipeline_log_", timestamp, ".txt"))

# Initialize Log
cat(sprintf("=== PIPELINE STARTED: %s ===\n", Sys.time()), file = log_file)
message(sprintf("[System] Saving full log to: %s", log_file))

# 3. Execution Wrapper
# ------------------------------------------------------------------------------
# Redirects output to both Console and File
sink(file = log_file, append = TRUE, split = TRUE)

tryCatch({
  
  # Load Modules inside the tryCatch to log loading errors
  cat("\n>>> LOADING MODULES <<<\n")
  list.files("R", pattern = "\\.R$", full.names = TRUE) %>% walk(source)
  message("[System] Modules loaded successfully.")
  
  # ==============================================================================
  # STEP 01 & 02: ETL & VISUALIZATION
  # ==============================================================================
  cat("\n>>> RUNNING PHASE 1: DATA PROCESSING & VIZ <<<\n")
  
  # Execute scripts using source() to keep variables in environment
  source("src/01_data_processing.R", echo = FALSE)
  source("src/02_visualization.R", echo = FALSE)
  
  # ==============================================================================
  # STEP 03: ANALYTICAL WORKFLOWS (Multi-Scenario)
  # ==============================================================================
  cat("\n>>> RUNNING PHASE 2: COMPARATIVE ANALYSIS (SCENARIOS) <<<\n")
  
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
    
    # Execute Scenarios
    all_results <- list()
    
    for (scenario in config$analysis_scenarios) {
      cat(sprintf("\n--- Executing Scenario: %s ---\n", scenario$id))
      
      tryCatch({
        res <- run_comparative_workflow(
          data_list = processed_data, 
          scenario = scenario, 
          config = config, 
          output_root = results_root
        )
        all_results[[scenario$id]] <- res
        
      }, error = function(e) {
        # Log specific scenario error but allow loop to continue
        message(sprintf("!!! [Error] Scenario '%s' failed: %s", scenario$id, e$message))
      })
    }
    
    # Consolidated Reporting
    if (length(all_results) > 0) {
      message("\n[Main] Generating Consolidated Multi-Scenario Report...")
      
      wb_master <- createWorkbook()
      
      # A. Global Drivers
      combined_drivers <- map_dfr(all_results, ~ .x$drivers, .id = "Scenario")
      if (nrow(combined_drivers) > 0) {
        addWorksheet(wb_master, "All_Drivers_Summary")
        writeData(wb_master, "All_Drivers_Summary", combined_drivers)
      }
      
      # B. Scenario Sheets
      for (scen_id in names(all_results)) {
        res <- all_results[[scen_id]]
        if (is.null(res)) next
        
        # PERMANOVA
        if (!is.null(res$permanova)) {
          sheet_name <- substr(paste0(scen_id, "_Perm"), 1, 31)
          if(!sheet_name %in% names(wb_master)) addWorksheet(wb_master, sheet_name)
          writeData(wb_master, sheet_name, res$permanova, rowNames = TRUE)
        }
        
        # Drivers
        if (!is.null(res$drivers) && nrow(res$drivers) > 0) {
          sheet_name <- substr(paste0(scen_id, "_Drv"), 1, 31)
          if(!sheet_name %in% names(wb_master)) addWorksheet(wb_master, sheet_name)
          writeData(wb_master, sheet_name, res$drivers)
        }
        
        # Network Diff
        if (!is.null(res$network$edges_table) && nrow(res$network$edges_table) > 0) {
          sheet_name <- substr(paste0(scen_id, "_Net"), 1, 31)
          if(!sheet_name %in% names(wb_master)) addWorksheet(wb_master, sheet_name)
          writeData(wb_master, sheet_name, res$network$edges_table)
        }
      }
      
      out_file <- file.path(results_root, "Multi_Scenario_Analysis_Report.xlsx")
      saveWorkbook(wb_master, out_file, overwrite = TRUE)
      message(sprintf("       -> Saved Master Report: %s", out_file))
    }
  }
  
  final_msg <- sprintf("\n=== PIPELINE FINISHED SUCCESSFULLY: %s ===", Sys.time())
  message(final_msg)
  
}, error = function(e) {
  # Fatal Error Handler
  cat(paste0("\n[FATAL ERROR] Pipeline stopped: ", e$message, "\n"), file = log_file, append = TRUE)
  stop(e$message)
  
}, finally = {
  # Cleanup sink
  sink()
  options(crayon.enabled = TRUE)
})