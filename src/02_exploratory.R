# src/02_exploratory.R
# ==============================================================================
# STEP 02: EXPLORATORY DATA ANALYSIS (PCA)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(FactoMineR)
  library(factoextra)
})

# Source modules
source("R/infrastructure.R")
source("R/modules_viz.R")

message("\n=== PIPELINE STEP 2: EXPLORATORY ANALYSIS ===")

# 1. Load Config & Data
config <- load_config("config/global_params.yml")
input_file <- file.path(config$output_root, "01_processed", "data_processed.rds")

if (!file.exists(input_file)) stop("Step 01 output not found. Run src/01_ingest.R first.")

DATA <- readRDS(input_file)
df_clr <- DATA$clr_data
markers <- DATA$markers

# 2. Prepare Output Directory
out_dir <- file.path(config$output_root, "02_exploratory")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(sprintf("[Data] Analyzable Set: %d Samples x %d Markers", nrow(df_clr), length(markers)))

# 3. PCA Analysis
message("[PCA] Running PCA on CLR transformed data...")

# Isolate numeric matrix for PCA
# Note: df_clr contains Metadata columns + Marker columns. We select only markers.
pca_input <- df_clr %>% 
  select(all_of(markers)) %>% 
  as.matrix()

res_pca <- run_coda_pca(pca_input)

# 4. Visualization
message("[Viz] Generating plots...")
my_colors <- get_palette(config)

# A. Main PCA Plot (Biplot-like with ellipses)
# CHANGED: Added show_labels = TRUE
p_pca <- plot_pca_custom(res_pca, df_clr[, c("Patient_ID", "Group")], my_colors, show_labels = TRUE)
ggsave(file.path(out_dir, "PCA_Score_Plot.pdf"), p_pca, width = 8, height = 6)

# B. Scree Plot (Variance Explained)
p_scree <- fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 50)) + 
  theme_coda() +
  labs(title = "Scree Plot: Variance Explained")
ggsave(file.path(out_dir, "PCA_Scree_Plot.pdf"), p_scree, width = 6, height = 4)

# C. Variable Contribution (Loadings)
p_vars <- plot_pca_variables(res_pca)
ggsave(file.path(out_dir, "PCA_Variable_Contrib.pdf"), p_vars, width = 8, height = 6)

message(sprintf("[Output] Plots saved to: %s", out_dir))
message("=== STEP 2 COMPLETE ===\n")