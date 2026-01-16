# src/02_visualization.R
# ==============================================================================
# STEP 02: VISUALIZATION (GLOBAL OVERVIEW)
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

message("\n=== PIPELINE STEP 2: VISUALIZATION (GLOBAL) ===")

# 1. Load Config & Data
config <- load_config("config/global_params.yml")
input_file <- file.path(config$output_root, "01_data_processing", "data_processed.rds")
if (!file.exists(input_file)) stop("Step 01 output not found.")

DATA <- readRDS(input_file)

# 2. Setup Data (NO STRATIFICATION HERE - WE WANT RAW GROUPS)
# ------------------------------------------------------------------------------
# We use the ORIGINAL metadata names (NSCLC, HNSCC) for the overview plots.
meta_viz <- DATA$metadata
safe_markers <- DATA$hybrid_markers
raw_matrix   <- DATA$raw_matrix
mat_z_global <- as.matrix(DATA$hybrid_data_z[, safe_markers]) 

# 3. Colors Setup (Robust)
# ------------------------------------------------------------------------------
# Function to safely find a color for a group name
get_color_safe <- function(grp_name, cfg) {
  # 1. Check specific assignment in config
  if (!is.null(cfg$colors$groups[[grp_name]])) return(cfg$colors$groups[[grp_name]])
  # 2. Check control
  if (grp_name == cfg$control_group) return(cfg$colors$control)
  # 3. Fallback
  return("grey50")
}

# 3a. Palette for Main Groups (Healthy, NSCLC, HNSCC)
unique_groups <- unique(meta_viz$Group)
colors_viz <- setNames(sapply(unique_groups, function(g) get_color_safe(g, config)), unique_groups)

# 3b. [NEW] Palette for Original_Source (Subgroups like NSCLC_LS)
# We map these sub-labels to the colors defined in config$colors$groups
if ("Original_Source" %in% colnames(meta_viz)) {
  unique_sources <- unique(meta_viz$Original_Source)
  colors_source <- setNames(sapply(unique_sources, function(s) {
    # Try to find direct match in config (e.g., "NSCLC_LS")
    # If not found, try to match by prefix or fallback to gray
    if (!is.null(config$colors$groups[[s]])) {
      return(config$colors$groups[[s]])
    } else if (grepl(config$control_group, s)) {
      return(config$colors$control)
    } else {
      # Heuristic: if 'NSCLC' is in the source name, use NSCLC color, etc.
      parent <- unique_groups[sapply(unique_groups, function(g) grepl(g, s))]
      if(length(parent) > 0) return(colors_viz[[parent[1]]])
      return("grey70")
    }
  }), unique_sources)
} else {
  colors_source <- NULL
}

# Highlight Pattern (The Golden Dots)
hl_pattern <- if(!is.null(config$viz$highlight_pattern)) config$viz$highlight_pattern else ""

out_dir <- file.path(config$output_root, "02_visualization")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 4. Raw Distributions (Unified by Cohort)
# ------------------------------------------------------------------------------
message("[Viz] Generating unified distribution plots...")

meta_ordered <- meta_viz[match(rownames(raw_matrix), meta_viz$Patient_ID), ]
df_raw_viz <- cbind(meta_ordered, as.data.frame(raw_matrix))

viz_save_distribution_report(
  data_df = df_raw_viz, 
  markers = DATA$markers, 
  file_path = file.path(out_dir, "Distributions_Raw_Distinct.pdf"), 
  colors = colors_viz, 
  hl_pattern = hl_pattern 
)

# 5. PCA (Global)
# ------------------------------------------------------------------------------
message("[PCA] Running Global PCA...")
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

# 6. Heatmap (Global)
# ------------------------------------------------------------------------------
message("[Viz] Generating Global Heatmap...")

target_cols <- config$viz$heatmap_metadata
if (is.null(target_cols)) target_cols <- c("Group")
valid_cols <- intersect(target_cols, names(meta_viz))
meta_heatmap <- meta_viz[, valid_cols, drop = FALSE]

# [FIX] Add colors_source to the annotation list
annotation_colors <- list(Group = colors_viz)
if (!is.null(colors_source) && "Original_Source" %in% names(meta_heatmap)) {
  annotation_colors[["Original_Source"]] = colors_source
}

pdf(file.path(out_dir, "Heatmap_Stratification.pdf"), width = 10, height = 8)
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