# main.R
# ==============================================================================
# MAIN PIPELINE ORCHESTRATOR
# Description: Executes steps 01-05 sequentially with Real-Time TEE logging.
# ==============================================================================

# Clean environment
rm(list = ls())
graphics.off()

# 1. Setup Logging & Clean Output
# ------------------------------------------------------------------------------
library(yaml)

# Disable ANSI colors
options(crayon.enabled = FALSE)

config_path <- "config/global_params.yml"
if (!file.exists(config_path)) stop("Config file not found!")

config <- yaml::read_yaml(config_path)

# Define Output Structure
out_root <- config$output_root
if (!is.null(config$project_name) && config$project_name != "") {
  out_root <- paste0(out_root, "_", config$project_name)
}

if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE)

# Define Log File
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(out_root, paste0("pipeline_log_", timestamp, ".txt"))

# Initialize Log File
cat(sprintf("=== PIPELINE STARTED: %s ===\n", Sys.time()), file = log_file)
message(sprintf("[System] Saving full log to: %s", log_file))

# 2. Advanced Logger Function
# ------------------------------------------------------------------------------
# Captures messages/warnings and writes them to file
run_pipeline_step <- function(script_path) {
  
  if (!file.exists(script_path)) stop(paste("Script missing:", script_path))
  
  # Append separator to log
  cat(paste0("\n>>> RUNNING: ", basename(script_path), " <<<\n"), 
      file = log_file, append = TRUE)
  
  # Execute with Handlers
  withCallingHandlers({
    source(script_path, echo = FALSE)
  }, 
  # Handler for Messages (Green text in console)
  message = function(m) {
    cat(conditionMessage(m), file = log_file, append = TRUE, sep = "\n")
  },
  # Handler for Warnings
  warning = function(w) {
    cat(paste0("WARNING: ", conditionMessage(w)), file = log_file, append = TRUE, sep = "\n")
  })
}

# 3. Execution Block 
# ------------------------------------------------------------------------------
# sink(..., split=TRUE) handles standard prints/cats, sending them to both console and file
sink(file = log_file, append = TRUE, split = TRUE)

tryCatch({
  
  # Step 1: Data Processing & QC
  run_pipeline_step("src/01_data_processing.R")
  
  # Step 2: Visualization
  run_pipeline_step("src/02_visualization.R")
  
  # Step 3: Statistical Analysis
  run_pipeline_step("src/03_statistical_analysis.R")
  
  # Step 4: Network Inference
  run_pipeline_step("src/04_network_inference.R")
  
  # Step 5: Network Topology
  run_pipeline_step("src/05_network_topology.R")
  
  final_msg <- sprintf("\n=== PIPELINE FINISHED SUCCESSFULLY: %s ===", Sys.time())
  message(final_msg)
  
}, error = function(e) {
  # Error Handling
  err_msg <- paste("\n[FATAL ERROR] Pipeline stopped unexpectedly.",
                   "Error Message:", e$message, sep = "\n")
  
  # Write to file explicitly because sink might be unstable during error
  cat(err_msg, file = log_file, append = TRUE)
  # Print to console (red)
  stop(err_msg)
  
}, finally = {
  # 4. Cleanup
  # ------------------------------------------------------------------------------
  sink() # Turn off stdout diversion
  options(crayon.enabled = TRUE) 
})