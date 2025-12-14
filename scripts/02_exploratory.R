# scripts/02_exploratory.R
# ==============================================================================
# STEP 02: EXPLORATORY DATA ANALYSIS (PCA & Stats)
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
  library(FactoMineR)
  library(factoextra)
})

# Load Config & Utils
config <- read_yaml("config/global_params.yml")
source("R/utils_viz.R")

# Paths
input_rds <- file.path("data", "processed", "clean_data.rds")
output_dir <- file.path(config$output_root, "02_Exploratory")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message("=== PIPELINE STEP 2: EXPLORATORY ANALYSIS ===")

# 1. Load Data
if(!file.exists(input_rds)) stop("Run Step 01 first!")
DATA <- readRDS(input_rds)
df_clr <- DATA$clr
markers <- DATA$markers

message(sprintf("   -> Loaded Data: %d Samples x %d Markers", nrow(df_clr), length(markers)))

# 2. PCA Analysis
message("[1] Running PCA on CLR data...")
# PCA on markers only, scale.unit = FALSE because CLR is already akin to scaling 
# (though variance scaling is debated in CoDa, usually unscaled is safer for geometry preservation)
res_pca <- PCA(df_clr[, markers], scale.unit = FALSE, graph = FALSE)

# 3. Generate & Save Plots
message("[2] Generating Plots...")
my_colors <- get_group_colors(config)

# A. PCA Plot
p_pca <- plot_pca_custom(res_pca, df_clr, my_colors)
ggsave(file.path(output_dir, "PCA_Plot.pdf"), p_pca, width = 8, height = 6)

# B. Scree Plot (Variance Explained)
p_scree <- fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 50)) + theme_coda()
ggsave(file.path(output_dir, "PCA_ScreePlot.pdf"), p_scree, width = 6, height = 4)

# C. Variable Contribution (Biplot arrows style)
p_var <- fviz_pca_var(res_pca, col.var = "contrib", 
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                      repel = TRUE) + 
         theme_coda() + labs(title = "Marker Contribution to Variance")
ggsave(file.path(output_dir, "PCA_Variables.pdf"), p_var, width = 8, height = 6)

message("=== PIPELINE STEP 2 COMPLETE ===")
