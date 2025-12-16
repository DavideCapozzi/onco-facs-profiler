# src/02_exploratory.R
# ==============================================================================
# STEP 02: EXPLORATORY DATA ANALYSIS (PCA)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(FactoMineR)
  library(factoextra)
  library(openxlsx)
})

# Source modules
source("R/infrastructure.R")
source("R/modules_stats.R")
source("R/modules_viz.R")


message("\n=== PIPELINE STEP 2: EXPLORATORY ANALYSIS ===")

# 1. Load Config & Data
config <- load_config("config/global_params.yml")
input_file <- file.path(config$output_root, "01_QC", "data_processed.rds")

if (!file.exists(input_file)) stop("Step 01 output not found. Run src/01_ingest.R first.")

DATA <- readRDS(input_file)
df_clr <- DATA$clr_data
df_ilr <- DATA$ilr_data
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

message("[Stats] Running MANOVA on ILR transformed data...")

# 5. Execute MANOVA
# Note: ILR is used because it has a full-rank covariance matrix (unlike CLR)
manova_results <- run_coda_manova(df_ilr, df_clr, metadata_cols = c("Patient_ID", "Group"))

message(sprintf("   -> MANOVA Result: Pillai = %.3f, p-value = %.5f", 
                manova_results$pillai_stat, manova_results$p_value))

if (manova_results$p_value < 0.05) {
  message("   -> SIGNIFCANT global difference detected between groups.")
} else {
  message("   -> NO significant global difference detected.")
}

# 6. Visualization (Boxplot of Canonical Variate)
p_manova <- plot_manova_results(manova_results, my_colors)
ggsave(file.path(out_dir, "MANOVA_Separation_Boxplot.pdf"), p_manova, width = 6, height = 5)

# 6b. Visualization of Drivers (Loadings)
p_loadings <- plot_manova_loadings(manova_results, top_n = 20)
ggsave(file.path(out_dir, "MANOVA_Drivers_Lollipop.pdf"), p_loadings, width = 7, height = 6)

# 7. Save Detailed Excel Report
wb <- createWorkbook()

addWorksheet(wb, "MANOVA_Results")
writeData(wb, "MANOVA_Results", as.data.frame(manova_results$summary), rowNames = TRUE)

addWorksheet(wb, "Canonical_Scores")
writeData(wb, "Canonical_Scores", manova_results$plot_data)

addWorksheet(wb, "Marker_Loadings")
writeData(wb, "Marker_Loadings", manova_results$loadings)

excel_path <- file.path(out_dir, "Statistical_Test_Results.xlsx")
saveWorkbook(wb, excel_path, overwrite = TRUE)

message(sprintf("[Output] Plots and Excel saved to: %s", out_dir))
message("=== STEP 2 COMPLETE ===\n")