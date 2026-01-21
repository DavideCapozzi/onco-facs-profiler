# main.R
# ==============================================================================
# PROJECT: ROBUST COMPOSITIONAL DATA ANALYSIS PIPELINE (FACS/CYTOF)
# AUTHOR: Davide Capozzi
# DATE: 21/01/2026
# ==============================================================================
setwd("/home/davidec/projects/compositional_analysis")

# A. Global Environment Setup
# Clean workspace to ensure reproducibility
rm(list = ls())
graphics.off()

# Load Infrastructure Module (for logging and config)
source("R/utils_io.R")

# Initialize Logger (Optional but recommended, otherwise use base message)
# If 'logger' package is not installed, simple messages will work via the scripts
message("[MAIN] Initializing Project Environment...")

# B. Configuration Check
# Defines the single source of truth for the entire run
config_path <- "config/global_params.yml"
if (!file.exists(config_path)) stop("Configuration file missing!")

# C. Pipeline Orchestration
# --- STEP 1: Data Ingestion & CoDa Preprocessing ---
# Handles: Excel loading -> Filter -> Zero Repl (CZM) ->  Imputation -> CLR/logit
message("\n========================================================\nSTEP 1: STARTING DATA PROCESSING\n========================================================\n")
source("src/01_data_processing.R", echo = FALSE)

# --- STEP 2: Visualization ---
# Handles: PCA, Scree Plots, Outlier visualization
message("\n========================================================\nSTEP 2: STARTING VISUALIZATION\n========================================================\n")
source("src/02_visualization.R", echo = FALSE)

# --- STEP 3: Statistical Analysis ---
# Handles: PERMANOVA, sPLS-DA, ILR Decoding
message("\n========================================================\nSTEP 3: STARTING STATISTICAL ANALYSIS\n========================================================\n")
source("src/03_statistical_analysis.R", echo = FALSE)

# --- STEP 4: Network Inference ---
# Handles: Bootstrap Networks, Permutation Tests, FDR
message("\n========================================================\nSTEP 4: STARTING NETWORK INFERENCE\n========================================================\n")
source("src/04_network_inference.R", echo = FALSE)

# --- STEP 5: Statistical Inference ---
# Handles: Network Topology Analysis
message("\n========================================================\nSTEP 5: STARTING NETWORK TOPOLOGY\n========================================================\n")
source("src/05_network_topology.R", echo = FALSE)

# D. Completion
message("\n========================================================\n[SUCCESS] PIPELINE FINISHED.\n========================================================\n")