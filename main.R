# main.R
# ==============================================================================
# MAIN PIPELINE ORCHESTRATOR
# Description: End-to-end execution of the Immunological Analysis Pipeline.
#              Orchestrates ETL, Visualization, and Multi-Scenario Analysis.
# ==============================================================================

# 1. Environment Setup
rm(list = ls())
graphics.off()

suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(here)
  library(openxlsx)
})

# Disable ANSI colors for clean log files
options(crayon.enabled = FALSE)

# Load Configuration
config_path <- here("config/global_params.yml")
if (!file.exists(config_path)) stop("Config file not found!")
config <- yaml::read_yaml(config_path)

# 2. Logging Setup
# ------------------------------------------------------------------------------
out_root <- here(config$output_root)
if (!is.null(config$project_name) && config$project_name != "") {
  out_root <- paste0(out_root, "_", config$project_name)
}
if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(out_root, paste0("pipeline_log_", timestamp, ".txt"))

# Initialize Log
cat(sprintf("=== PIPELINE STARTED: %s ===\n", Sys.time()), file = log_file)
message(sprintf("[System] Saving full log to: %s", log_file))

# Helper to log messages to file AND console
log_handler <- function(m) {
  cat(conditionMessage(m), file = log_file, append = TRUE, sep = "\n")
  # Message bubbles up to console automatically
}

warn_handler <- function(w) {
  cat(paste0("WARNING: ", conditionMessage(w)), file = log_file, append = TRUE, sep = "\n")
  # Warning bubbles up to console automatically
}

# 3. Pipeline Execution
# ------------------------------------------------------------------------------
# We wrap the entire execution in withCallingHandlers to capture output reliably
tryCatch({
  withCallingHandlers({
    
    # --- Module Loading ---
    message("\n>>> LOADING MODULES <<<")
    list.files(here("R"), pattern = "\\.R$", full.names = TRUE) %>% walk(source)
    message("[System] Modules loaded successfully.")
    
    # --- PHASE 1: DATA PROCESSING & VIZ ---
    message("\n>>> RUNNING PHASE 1: DATA PROCESSING & VIZ <<<")
    
    # Sourcing these files executes the steps immediately
    source(here("src/01_data_processing.R"), echo = FALSE)
    source(here("src/02_visualization.R"), echo = FALSE)
    
    # --- PHASE 2: COMPARATIVE ANALYSIS (SCENARIOS) ---
    message("\n>>> RUNNING PHASE 2: COMPARATIVE ANALYSIS (SCENARIOS) <<<")
    
    data_file <- file.path(here(config$output_root), "01_data_processing", "data_processed.rds")
    if (!file.exists(data_file)) stop("CRITICAL: Processed data not found.")
    processed_data <- readRDS(data_file)
    
    results_root <- file.path(here(config$output_root), "results_analysis")
    if (!dir.exists(results_root)) dir.create(results_root, recursive = TRUE)
    
    if (is.null(config$analysis_scenarios) || length(config$analysis_scenarios) == 0) {
      warning("[Main] No 'analysis_scenarios' found in config. Skipping Phase 2.")
    } else {
      
      all_results <- list()
      
      for (scenario in config$analysis_scenarios) {
        message(sprintf("\n--- Executing Scenario: %s ---", scenario$id))
        
        tryCatch({
          res <- run_comparative_workflow(
            data_list = processed_data, 
            scenario = scenario, 
            config = config, 
            output_root = results_root
          )
          all_results[[scenario$id]] <- res
          
        }, error = function(e) {
          # Log error but continue pipeline
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
    
  }, message = log_handler, warning = warn_handler)
  
}, error = function(e) {
  # Fatal Error Handler
  err_msg <- paste0("\n[FATAL ERROR] Pipeline stopped: ", e$message)
  cat(err_msg, file = log_file, append = TRUE, sep = "\n")
  stop(e$message)
  
}, finally = {
  options(crayon.enabled = TRUE)
})