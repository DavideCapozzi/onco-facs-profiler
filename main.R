# main.R
# ==============================================================================
# PROJECT: ROBUST COMPOSITIONAL DATA ANALYSIS PIPELINE (FACS/CYTOF)
# AUTHOR: [Your Name/Lab]
# DATE: 2025-12-07
# ==============================================================================
setwd("/home/davidec/projects/compositional_analysis")

# 1. Global Environment Setup
# Clean workspace to ensure reproducibility
rm(list = ls())
graphics.off()

# Load Infrastructure Module (for logging and config)
source("R/infrastructure.R")

# Initialize Logger (Optional but recommended, otherwise use base message)
# If 'logger' package is not installed, simple messages will work via the scripts
message("[MAIN] Initializing Project Environment...")

# 2. Configuration Check
# Defines the single source of truth for the entire run
config_path <- "config/global_params.yml"
if (!file.exists(config_path)) stop("Configuration file missing!")

# 3. Pipeline Orchestration

# --- STEP 1: Data Ingestion & CoDa Preprocessing ---
# Handles: Excel loading -> Filter -> Zero Repl (CZM) -> kNN Opt -> Imputation -> CLR
message("\n========================================================")
message(" STEP 1: STARTING INGESTION & PREPROCESSING")
message("========================================================")
source("src/01_ingest.R", echo = FALSE)


# --- STEP 2: Exploratory Data Analysis (EDA) ---
# Handles: PCA, Scree Plots, Outlier visualization
message("\n========================================================")
message(" STEP 2: STARTING EXPLORATORY ANALYSIS (PCA)")
message("========================================================")
source("src/02_exploratory.R", echo = FALSE)


# --- STEP 3: Statistical Inference ---
# Handles: Bootstrap Networks, Permutation Tests, FDR
message("\n========================================================")
message(" STEP 3: STARTING INFERENCE ENGINE")
message("========================================================")
source("src/03_inference.R", echo = FALSE)

# --- STEP 3: Statistical Inference ---
# Handles: Bootstrap Networks, Permutation Tests, FDR
message("\n========================================================")
message(" STEP 4: STARTING NETWORK ANALYSIS")
message("========================================================")
source("src/04_network_analysis.R", echo = FALSE)

# 5. Completion
message("\n========================================================")
message(" [SUCCESS] PIPELINE FINISHED.")
message(" Results stored in: results/")
message("========================================================")