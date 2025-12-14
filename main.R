# main.R
# ==============================================================================
# MASTER PIPELINE RUNNER
# ==============================================================================
# 0. Set working directory
setwd("/home/davidec/projects/compositional_analysis/")

# 1. Environment Check
if(!file.exists("config/global_params.yml")) stop("Config file missing!")
if(!dir.exists("data")) dir.create("data")

# 2. Run Modules sequentially
run_step <- function(script_path) {
  message(paste0("\n\n>>> RUNNING: ", script_path))
  source(script_path, local = new.env()) # Run in isolated env to prevent leaks
}

# Step 1: Data Processing (Ingest -> Impute -> CLR)
run_step("scripts/01_data_processing.R")

# Step 2: Exploratory (PCA -> Plots)
run_step("scripts/02_exploratory.R")

# Step 3: Network Inference (Bootstrap -> Permutation)
# WARNING: This step takes time!
run_step("scripts/03_inference.R")

message("\n\n>>> PIPELINE FINISHED SUCCESSFULLY.")
