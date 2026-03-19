# test/test_epsilon_robustness.R
# ==============================================================================
# EPSILON DIAGNOSTICS & ROBUSTNESS TESTING
# Description: Standalone diagnostic tool to evaluate the impact of different
#              Logit epsilon values (pseudo-counts) on the final multivariate
#              feature selection. Run this before finalizing global_params.yml.
# ==============================================================================

# 1. Environment Setup
rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# Source required pipeline modules
source("R/utils_io.R")
source("R/modules_coda.R")
source("R/modules_multivariate.R")
source("R/modules_viz.R")

message("\n=== DIAGNOSTICS: EPSILON SENSITIVITY ANALYSIS ===")

# 2. Configuration & Data Loading
config <- load_config("config/global_params.yml")
input_file <- file.path(config$output_root, "01_data_processing", "data_processed.rds")

if (!file.exists(input_file)) {
  stop("[Fatal] Step 01 output not found. Run src/01_data_processing.R first to obtain clean QC data.")
}
DATA <- readRDS(input_file)

# Dedicated diagnostic output directory
diag_dir <- file.path(config$output_root, "results_diagnostics", "Epsilon_Robustness")
if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)

# 3. Define Grid of Values to Test
# Includes standard magnitudes and the data-driven "auto" mode
eps_grid <- list(
  "Eps_1e-3" = 1e-3,
  "Eps_1e-4" = 1e-4,
  "Eps_1e-5" = 1e-5,
  "Eps_1e-6" = 1e-6,
  "Eps_Auto" = "auto"
)

# Use raw matrix (post-QC, pre-transformation)
mat_raw <- DATA$raw_matrix

# Determine global contrast groups for sPLS-DA testing
target_all <- unique(c(config$control_group, config$case_groups))
strat_col <- config$stratification$column

meta_test <- DATA$metadata %>%
  dplyr::filter(!!rlang::sym(strat_col) %in% target_all) %>%
  dplyr::mutate(Group = !!rlang::sym(strat_col))

# Align matrix
valid_samples <- intersect(rownames(mat_raw), meta_test$Patient_ID)
mat_raw <- mat_raw[valid_samples, , drop = FALSE]
meta_test <- meta_test[match(valid_samples, meta_test$Patient_ID), , drop = FALSE]

# 4. Iterative Testing
top_features_list <- list()
message(sprintf("   [Diagnostics] Testing %d epsilon configurations...", length(eps_grid)))

for (eps_name in names(eps_grid)) {
  current_eps <- eps_grid[[eps_name]]
  message(sprintf("\n      -> Running configuration: %s", eps_name))
  
  # Inject temporary parameter
  temp_config <- config
  temp_config$imputation$epsilon <- current_eps
  
  # Step A: Transform Data
  # We suppress messages to keep console clean during diagnostics
  transform_res <- suppressMessages(perform_data_transformation(mat_raw, temp_config, mode = "complete"))
  mat_z <- transform_res$hybrid_data_z
  
  # Step B: Fit Global sPLS-DA Model
  set.seed(temp_config$stats$seed)
  tryCatch({
    pls_res <- suppressMessages(run_splsda_model(
      data_z = mat_z, 
      metadata = meta_test, 
      group_col = "Group", 
      n_comp = 1, # Fast diagnostic check on PC1
      folds = min(5, min(table(meta_test$Group))), 
      n_repeat = 5 # Reduced repeats for diagnostic speed
    ))
    
    # Step C: Extract Top Drivers
    drivers <- extract_plsda_loadings(pls_res)
    
    if (nrow(drivers) > 0) {
      # Extract Top 15 robust markers for overlap
      top_n <- min(15, nrow(drivers))
      top_features_list[[eps_name]] <- head(drivers$Marker, top_n)
    }
  }, error = function(e) {
    message(sprintf("         [Error] Model failed for %s: %s", eps_name, e$message))
  })
}

# 5. Visual Reporting (Overlap Analysis)
message("\n   [Diagnostics] Generating Overlap Visualizations...")

if (length(top_features_list) >= 2) {
  
  pdf(file.path(diag_dir, "Top_Drivers_Overlap.pdf"), width = 9, height = 7)
  
  # Leverage the existing pipeline module for dynamic Venn/UpSet generation
  tryCatch({
    p_overlap <- viz_plot_differential_overlap(
      edge_list = top_features_list, 
      title = "Robustness Analysis: Top 15 sPLS-DA Drivers Overlap"
    )
    # If the function returns a plot object (Venn), print it. 
    # If UpSet, it draws internally via ComplexHeatmap.
    if (inherits(p_overlap, "gg")) print(p_overlap)
    
  }, error = function(e) message("      [Warn] Overlap plotting failed: ", e$message))
  
  dev.off()
  
  message(sprintf("   [Output] Diagnostic PDF saved to: %s", diag_dir))
  message("   [Interpretation] High intersection indicates stable biological signals independent of Logit boundaries.")
  
} else {
  message("   [Warn] Insufficient successful runs to perform overlap diagnostics.")
}

message("\n=== DIAGNOSTICS COMPLETE ===\n")
