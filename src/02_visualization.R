# src/02_visualization.R
# ==============================================================================
# STEP 02: VISUALIZATION & DATA EXPLORATION
# Description: Descriptive analysis (Raw Distributions, PCA, Heatmaps).
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(FactoMineR)
  library(factoextra)
  library(ComplexHeatmap)
  library(circlize)
})

source("R/utils_io.R")          
source("R/modules_viz.R")

message("\n=== PIPELINE STEP 2: VISUALIZATION ===")

# 1. Load Config & Data
# ------------------------------------------------------------------------------
config <- load_config("config/global_params.yml")
input_file <- file.path(config$output_root, "01_data_processing", "data_processed.rds")

if (!file.exists(input_file)) stop("Step 01 output not found. Run src/01_data_processing.R first.")

DATA <- readRDS(input_file)

# Prepare Data
meta_viz <- DATA$metadata
safe_markers <- DATA$hybrid_markers
raw_matrix   <- DATA$raw_matrix
# For visualization (PCA/Heatmap), we use the Z-scored hybrid data
mat_z_global <- as.matrix(DATA$hybrid_data_z[, safe_markers]) 

# Colors
colors_viz <- get_palette(config)
hl_pattern <- if(!is.null(config$viz$highlight_pattern)) config$viz$highlight_pattern else ""

# Prepare Output Directory
out_dir <- file.path(config$output_root, "02_visualization")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 2. Raw Distributions (PDF)
# ------------------------------------------------------------------------------
message("[Viz] Generating unified distribution plots for RAW data...")

# Combine metadata with raw matrix for plotting
meta_ordered <- meta_viz[match(rownames(raw_matrix), meta_viz$Patient_ID), ]
df_raw_viz <- cbind(meta_ordered, as.data.frame(raw_matrix))

raw_pdf_path <- file.path(out_dir, "Distributions_Raw_Distinct.pdf")

viz_save_distribution_report(
  data_df = df_raw_viz, 
  markers = DATA$markers, 
  file_path = raw_pdf_path, 
  colors = colors_viz, 
  hl_pattern = hl_pattern
)

# 3. PCA (Exploratory) - Multi-Page
# ------------------------------------------------------------------------------
message("[PCA] Running PCA on Global Z-Scored Matrix...")

res_pca <- PCA(mat_z_global, scale.unit = FALSE, graph = FALSE)

# Setup Highlight Pattern Map
highlight_map <- c()
if (!is.null(hl_pattern) && hl_pattern != "") {
  highlight_map[[hl_pattern]] <- "#FFD700"
}

# A. Individuals PCA plots (Multi-Page)
message("[PCA] Saving Individuals plots...")
pdf_ind_path <- file.path(out_dir, "PCA_Global_Individuals.pdf")
pdf(pdf_ind_path, width = 9, height = 7)

# Loop common component pairs
for (dims in list(c(1,2), c(1,3), c(2,3))) {
  print(plot_pca_custom(res_pca, meta_viz, colors_viz, dims = dims, 
                        show_labels = TRUE, highlight_patterns = highlight_map))
}
dev.off() 

# B. Variables/Loadings plots
message("[PCA] Saving Variables/Loadings plots...")
pdf_var_path <- file.path(out_dir, "PCA_Global_Variables.pdf")
pdf(pdf_var_path, width = 8, height = 6)

plot_vars_dims <- function(pca_res, axes_vec) {
  fviz_pca_var(pca_res, axes = axes_vec, col.var = "contrib", 
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) + 
    labs(title = sprintf("Variables - PCA (PC%d vs PC%d)", axes_vec[1], axes_vec[2])) +
    theme_bw()
}

for (dims in list(c(1,2), c(1,3), c(2,3))) {
  print(plot_vars_dims(res_pca, dims))
}
dev.off()

# C. Scree Plot
pdf(file.path(out_dir, "PCA_Global_Scree.pdf"), width = 6, height = 4)
print(fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 50)) + 
        theme_minimal() + ggtitle("Scree Plot (Global)"))
dev.off()

# 4. Patient Stratification Heatmap
# ------------------------------------------------------------------------------
message("[Viz] Generating Patient Stratification Heatmap...")

# A. Dynamic Annotation Colors
target_cols <- config$viz$heatmap_metadata
if (is.null(target_cols)) target_cols <- c("Group")
valid_cols <- intersect(target_cols, names(meta_viz))
meta_heatmap <- meta_viz[, valid_cols, drop = FALSE]

annotation_colors <- viz_generate_complex_heatmap_colors(
  metadata = meta_heatmap,
  group_col = "Group",
  base_colors = colors_viz,
  hl_pattern = hl_pattern
)

# B. Output PDF
pdf_heat_path <- file.path(out_dir, "Heatmap_Stratification.pdf")
pdf(pdf_heat_path, width = 10, height = 8)

tryCatch({
  ht_obj <- plot_stratification_heatmap(
    mat_z = mat_z_global,
    metadata = meta_heatmap,
    annotation_colors_list = annotation_colors,
    title = "Clustering (Global)"
  )
  draw(ht_obj, merge_legend = TRUE)
}, error = function(e) warning(paste("Heatmap failed:", e$message)))

dev.off()

message("=== STEP 2 COMPLETE ===\n")