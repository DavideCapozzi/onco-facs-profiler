# src/02_visualization.R
# ==============================================================================
# STEP 02: VISUALIZATION
# Description: Generates overview plots (Distributions, PCA, Heatmap)
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
config <- load_config("config/global_params.yml")
input_file <- file.path(config$output_root, "01_data_processing", "data_processed.rds")
if (!file.exists(input_file)) stop("Step 01 output not found.")

DATA <- readRDS(input_file)

# 2. Data Setup & Filtering
# ------------------------------------------------------------------------------
strat_col <- config$stratification$column
if (!strat_col %in% names(DATA$metadata)) {
  stop(sprintf("Column '%s' not found in metadata.", strat_col))
}

# Define target subgroups 
target_subgroups <- unique(c(config$control_group, config$case_groups))
message(sprintf("[Viz] Filtering dataset for targets: %s", paste(target_subgroups, collapse=", ")))

# Filter metadata 
meta_viz <- DATA$metadata %>% 
  filter(.data[[strat_col]] %in% target_subgroups)

# Align matrices
mat_z_global <- DATA$hybrid_data_z %>%
  filter(Patient_ID %in% meta_viz$Patient_ID) %>%
  select(all_of(DATA$hybrid_markers)) %>%
  as.matrix()

raw_matrix <- DATA$raw_matrix[meta_viz$Patient_ID, , drop = FALSE]

# Validation
if (nrow(meta_viz) < 3) stop("Insufficient samples after filtering (<3). Check config.")
message(sprintf("[Viz] Analyzed Set: %d Samples across %d Macro-Groups", 
                nrow(meta_viz), length(unique(meta_viz$Group))))

# 3. Colors Setup 
# ------------------------------------------------------------------------------
# Ensure no NAs are passed to the palette generator
unique_groups <- unique(as.character(na.omit(meta_viz$Group)))

# Pass match_groups to activate root matching
full_palette <- get_palette(config, match_groups = unique_groups)

# Subset safely
colors_viz <- full_palette[unique_groups]

names(colors_viz) <- unique_groups 

# Fallback: Handle groups present in data but undefined in config
if (any(is.na(colors_viz))) {
  missing_grps <- unique_groups[is.na(colors_viz)]
  warning(sprintf("[Viz] Groups found in data but missing in config/palette: %s", 
                  paste(missing_grps, collapse=", ")))
  colors_viz[missing_grps] <- "grey50"
}

hl_pattern <- if(!is.null(config$viz$highlight_pattern)) config$viz$highlight_pattern else ""

hl_pattern <- if(!is.null(config$viz$highlight_pattern)) config$viz$highlight_pattern else ""

out_dir <- file.path(config$output_root, "02_visualization")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 4. Raw Distributions
# ------------------------------------------------------------------------------
message("[Viz] Generating distribution plots...")

df_raw_viz <- cbind(meta_viz, as.data.frame(raw_matrix))

viz_save_distribution_report(
  data_df = df_raw_viz, 
  markers = DATA$markers, 
  file_path = file.path(out_dir, "Distributions_Raw_Distinct.pdf"), 
  colors = colors_viz, 
  hl_pattern = hl_pattern 
)

# 5. PCA Analysis
# ------------------------------------------------------------------------------
message("[PCA] Running PCA...")
res_pca <- PCA(mat_z_global, scale.unit = FALSE, graph = FALSE)

highlight_map <- c()
if (hl_pattern != "") highlight_map[[hl_pattern]] <- config$colors$highlight

pdf(file.path(out_dir, "PCA_Global_Individuals.pdf"), width = 9, height = 7)

# Define dimensions to explore
target_dims_list <- list(c(1,2), c(1,3), c(2,3))

# Loop A: Individuals Plots
for (dims in target_dims_list) {
  print(plot_pca_custom(res_pca, meta_viz, colors_viz, dims = dims, 
                        show_labels = TRUE, highlight_patterns = highlight_map))
}

# Loop B: Variables Plots (Correlation Circle) [NEW ADDITION]
for (dims in target_dims_list) {
  # We use the standard function from modules_viz but explicitly setting axes
  p_var <- plot_pca_variables(res_pca) +
    labs(title = sprintf("Variables PCA (PC%d vs PC%d)", dims[1], dims[2]))
  
  p_var_custom <- fviz_pca_var(res_pca, 
                               axes = dims,
                               col.var = "contrib", 
                               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                               repel = TRUE) +
    labs(title = sprintf("Marker Contribution (PC%d vs PC%d)", dims[1], dims[2])) +
    theme_coda()
  
  print(p_var_custom)
}

dev.off()

# 6. Heatmap Analysis
# ------------------------------------------------------------------------------
message("[Viz] Generating Heatmap...")

target_cols <- config$viz$heatmap_metadata
if (is.null(target_cols)) target_cols <- c("Group")
valid_cols <- intersect(target_cols, names(meta_viz))
meta_heatmap <- meta_viz[, valid_cols, drop = FALSE]

annotation_colors <- list(Group = colors_viz)

if ("Original_Source" %in% names(meta_heatmap) && !identical(meta_heatmap$Group, meta_heatmap$Original_Source)) {
  annotation_colors[["Original_Source"]] <- viz_generate_complex_heatmap_colors(
    metadata = meta_heatmap,
    group_col = "Group",
    base_colors = colors_viz,
    hl_pattern = hl_pattern,
    hl_color = config$colors$highlight
  )[["Original_Source"]]
}

pdf(file.path(out_dir, "Heatmap_Stratification.pdf"), width = 10, height = 8)
tryCatch({
  ht_obj <- plot_stratification_heatmap(
    mat_z = mat_z_global,
    metadata = meta_heatmap,
    annotation_colors_list = annotation_colors,
    title = "Clustering Analysis"
  )
  draw(ht_obj, merge_legend = TRUE)
}, error = function(e) warning(paste("Heatmap failed:", e$message)))
dev.off()

message("=== STEP 2 COMPLETE ===\n")