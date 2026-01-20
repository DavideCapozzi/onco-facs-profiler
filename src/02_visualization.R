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
# Update the main grouping variable based on the stratification column
strat_col <- config$stratification$column
if (!strat_col %in% names(DATA$metadata)) {
  stop(sprintf("Column '%s' not found in metadata.", strat_col))
}

# Overwrite Group with the specific analysis column
DATA$metadata$Group <- DATA$metadata[[strat_col]]
DATA$hybrid_data_z$Group <- DATA$hybrid_data_z[[strat_col]]

# Define target groups from configuration
target_groups <- c(config$control_group, config$case_groups)
message(sprintf("[Viz] Filtering dataset for targets: %s", paste(target_groups, collapse=", ")))

# Filter Metadata
meta_viz <- DATA$metadata %>% 
  filter(Group %in% target_groups) %>%
  mutate(Group = factor(Group, levels = target_groups))

# Filter Z-Score Matrix
mat_z_global <- DATA$hybrid_data_z %>%
  filter(Patient_ID %in% meta_viz$Patient_ID) %>%
  select(all_of(DATA$hybrid_markers)) %>%
  as.matrix()

# Filter Raw Matrix
raw_matrix <- DATA$raw_matrix[meta_viz$Patient_ID, , drop = FALSE]

# Validation
if (nrow(meta_viz) < 3) stop("Insufficient samples after filtering (<3). Check config group names.")
message(sprintf("[Viz] Analyzed Set: %d Samples across %d Groups", nrow(meta_viz), length(unique(meta_viz$Group))))

# 3. Colors Setup
# ------------------------------------------------------------------------------
assign_color <- function(grp_name, cfg) {
  # 1. Check specific definition in config
  if (!is.null(cfg$colors$groups[[grp_name]])) return(cfg$colors$groups[[grp_name]])
  # 2. Check control definition
  if (grp_name == cfg$control_group) return(cfg$colors$control)
  # 3. Fallback
  return("grey50")
}

unique_groups <- levels(meta_viz$Group)
colors_viz <- setNames(sapply(unique_groups, function(g) assign_color(g, config)), unique_groups)

# Highlight Pattern (Optional)
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

# Setup Highlight Map for PCA borders
highlight_map <- c()
if (hl_pattern != "") highlight_map[[hl_pattern]] <- config$colors$highlight

pdf(file.path(out_dir, "PCA_Global_Individuals.pdf"), width = 9, height = 7)
for (dims in list(c(1,2), c(1,3), c(2,3))) {
  print(plot_pca_custom(res_pca, meta_viz, colors_viz, dims = dims, 
                        show_labels = TRUE, highlight_patterns = highlight_map))
}
dev.off() 

# 6. Heatmap Analysis
# ------------------------------------------------------------------------------
message("[Viz] Generating Heatmap...")

# Select metadata for annotation
target_cols <- config$viz$heatmap_metadata
if (is.null(target_cols)) target_cols <- c("Group")

valid_cols <- intersect(target_cols, names(meta_viz))
meta_heatmap <- meta_viz[, valid_cols, drop = FALSE]

# Define Annotation Colors
annotation_colors <- list(Group = colors_viz)

# If 'Original_Source' is present and distinct from Group, map it using the same palette if applicable
if ("Original_Source" %in% names(meta_heatmap) && !identical(meta_heatmap$Group, meta_heatmap$Original_Source)) {
  annotation_colors[["Original_Source"]] <- colors_viz 
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