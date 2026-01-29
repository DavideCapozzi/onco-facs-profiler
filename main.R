# main.R
# ==============================================================================
# MAIN PIPELINE ORCHESTRATOR
# Description: End-to-end execution of the Immunological Analysis Pipeline.
#              Supports single-run and multi-scenario campaigns via config.
# ==============================================================================

# 1. Environment Setup
suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(here)
})

here()
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
# These steps generate the "Master Dataset" used by all downstream analyses.
# We execute them as scripts to maintain the file-based checkpointing.

message("\n--- PHASE 1: DATA PROCESSING & VIZ ---")
source("src/01_data_processing.R")
source("src/02_visualization.R")


# ==============================================================================
# STEP 03: ANALYTICAL WORKFLOWS (Config-Driven)
# ==============================================================================
# Instead of static scripts, we now iterate over the scenarios defined in config.
# This handles the "Standard" analysis and the "LS Characterization" uniformly.

message("\n--- PHASE 2: STATISTICAL & NETWORK ANALYSIS ---")

# Load the processed data once
data_file <- file.path(config$output_root, "01_data_processing", "data_processed.rds")
if (!file.exists(data_file)) stop("CRITICAL: Processed data not found.")
processed_data <- readRDS(data_file)

# Define where results go
results_root <- file.path(config$output_root, "results_analysis")

# Check for scenarios in config
if (is.null(config$analysis_scenarios) || length(config$analysis_scenarios) == 0) {
  warning("[Main] No 'analysis_scenarios' found in config. Skipping Phase 2.")
} else {
  
  # Iterate and Execute
  all_results <- list()
  
  for (scenario in config$analysis_scenarios) {
    tryCatch({
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
  
  # Optional: Global Synthesis (e.g. Venn Diagram of Drivers)
  # If we have multiple scenarios, we can compare them here.
  if (length(all_results) > 1) {
    message("\n[Main] Generating Global Synthesis across scenarios...")
    # Example: Extract drivers from all scenarios and save a summary table
    combined_drivers <- map_dfr(all_results, ~ .x$drivers, .id = "Scenario")
    if (nrow(combined_drivers) > 0) {
      write.csv(combined_drivers, file.path(results_root, "Global_Drivers_Summary.csv"), row.names = FALSE)
      message("       -> Saved 'Global_Drivers_Summary.csv'")
    }
  }
}

message("\n=== PIPELINE COMPLETED SUCCESSFULLY ===")