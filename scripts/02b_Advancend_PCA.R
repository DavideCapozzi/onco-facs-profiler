# ==============================================================================
# SCRIPT 02 (FINAL): GLOBAL, SENSITIVITY & FOCUSED PCA ANALYSIS
# PURPOSE: PCA on All markers, Leave-One-Out, and Single-Group-Only
# INPUT: clean_data.rds
# OUTPUT: PDFs for Global, Without-Group, and Only-Group strategies
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(FactoMineR)
  library(factoextra)
  library(vegan)        
  library(ggpubr)
  library(openxlsx)
  library(ggrepel)
})

# 1. CONFIGURATION -------------------------------------------------------------
CONFIG <- list(
  input_file = "/home/davidec/projects/compositional_analysis/processed_data/clean_data.rds",
  base_output_dir = "/home/davidec/projects/compositional_analysis/results/02_Stats",
  # Sub-directories for specific analyses
  dir_sensitivity = "pca_without_single_groups", # All minus One
  dir_focused     = "pca_single_groups_only",    # Only One
  
  colors = c("Healthy" = "#2E8B57", "NSCLC" = "#4682B4", "HNSCC" = "#CD5C5C") 
)

# Setup directories
path_sensitivity <- file.path(CONFIG$base_output_dir, CONFIG$dir_sensitivity)
path_focused     <- file.path(CONFIG$base_output_dir, CONFIG$dir_focused)

if(!dir.exists(CONFIG$base_output_dir)) dir.create(CONFIG$base_output_dir, recursive = TRUE)
if(!dir.exists(path_sensitivity)) dir.create(path_sensitivity, recursive = TRUE)
if(!dir.exists(path_focused)) dir.create(path_focused, recursive = TRUE)

cat("=== PIPELINE STEP 2: ADVANCED PCA ANALYSIS ===\n")
cat("Output Root:", CONFIG$base_output_dir, "\n")

# 2. LOAD DATA -----------------------------------------------------------------
if(!file.exists(CONFIG$input_file)) stop("Step 1 output not found! Check input path.")
DATA <- readRDS(CONFIG$input_file)
df_clr <- DATA$clr_transformed

# Validate columns
cat("-> Data loaded. N Patients:", nrow(DATA$metadata), "\n")

# 3. DEFINE MARKER GROUPS ------------------------------------------------------
# Definition based on immunological function
marker_groups <- list(
  tcells = c("CD3", "CD4", "CD8"),
  tregs  = c("TREG", "RESTING", "ACTIVE", "NON-SUPPRESSIVE"),
  memory = c("CM", "EFF", "EM", "NAIVE"),
  cd137  = c("CD137", "CD137CD4", "CD137CD8", "CD137KI67", "CD137PD1", 
             "CD137CM", "CD137EFF", "CD137EMRA", "CD137NAIVE"),
  ki67   = c("KI67", "KI67CD4", "KI67CD8"),
  pd1    = c("PD1", "PD1CD4", "PD1CD8", "PD1CM", "PD1EFF", "PD1EMRA", "PD1NAIVE")
)

# Collect ALL marker columns to be used in the Global analysis
global_markers <- DATA$markers 

cat("-> Marker Groups defined.\n")

# 4. HELPER FUNCTIONS ----------------------------------------------------------

plot_pca_dims <- function(pca_data, x_dim, y_dim, eigen_data, group_colors, title_suffix) {
  
  x_col <- paste0("Dim.", x_dim)
  y_col <- paste0("Dim.", y_dim)
  
  # Check if dimensions exist (crucial for small groups like Ki67 which has only 3 markers)
  if(!x_col %in% colnames(pca_data) || !y_col %in% colnames(pca_data)) return(NULL)
  
  var_x <- eigen_data[x_dim, 2]
  var_y <- eigen_data[y_dim, 2]
  
  ggplot(pca_data, aes(x = .data[[x_col]], y = .data[[y_col]], color = Group, fill = Group)) +
    geom_point(size = 4, alpha = 0.8) +
    stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +
    geom_text_repel(aes(label = ID), size = 3, show.legend = FALSE, max.overlaps = 10) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    labs(title = paste0("PCA: ", title_suffix),
         subtitle = sprintf("Dim%d: %.1f%% | Dim%d: %.1f%%", x_dim, var_x, y_dim, var_y),
         caption = "Method: CLR Transform -> PCA",
         x = paste0("PC ", x_dim),
         y = paste0("PC ", y_dim)) +
    theme_bw()
}

run_pca_analysis <- function(data_df, markers_subset, file_name, plot_title_suffix, output_path) {
  
  # Check if markers actually exist in data
  valid_markers <- markers_subset[markers_subset %in% colnames(data_df)]
  
  # We need at least 3 markers to generate up to PC3 effectively for the plots
  if(length(valid_markers) < 2) {
    cat("   [!] Skipping", plot_title_suffix, "- Not enough markers (Need 2+).\n")
    return(NULL)
  }
  
  # 1. Compute PCA
  # We set ncp=5, but FactoMineR caps it automatically at min(nrow, ncol)
  pca_res <- PCA(data_df[, valid_markers], graph = FALSE, scale.unit = FALSE, ncp = 5)
  
  # 2. Prepare Data for Plotting
  pca_ind <- data.frame(pca_res$ind$coord, Group = data_df$Group, ID = data_df$Patient_ID)
  
  # 3. Generate Plots (Handling potential missing dims for small groups)
  plot_list <- list()
  
  p1 <- plot_pca_dims(pca_ind, 1, 2, pca_res$eig, CONFIG$colors, plot_title_suffix)
  if(!is.null(p1)) plot_list[[1]] <- p1
  
  p2 <- plot_pca_dims(pca_ind, 1, 3, pca_res$eig, CONFIG$colors, plot_title_suffix)
  if(!is.null(p2)) plot_list[[2]] <- p2
  
  p3 <- plot_pca_dims(pca_ind, 2, 3, pca_res$eig, CONFIG$colors, plot_title_suffix)
  if(!is.null(p3)) plot_list[[3]] <- p3
  
  # 4. Save to PDF
  full_path <- file.path(output_path, paste0(file_name, ".pdf"))
  pdf(full_path, width = 11, height = 7)
  for(p in plot_list) print(p)
  dev.off()
  
  return(pca_res)
}

# 5. EXECUTION -----------------------------------------------------------------

# --- A. GLOBAL ANALYSIS (All Markers) ---
cat("\n[1] Running GLOBAL PCA (All Markers)...\n")
res_global <- run_pca_analysis(df_clr, global_markers, "PCA_Global", "Global Landscape", CONFIG$base_output_dir)

# Export Detailed Stats (Only for Global)
wb <- createWorkbook()
addWorksheet(wb, "Eigenvalues"); writeData(wb, "Eigenvalues", as.data.frame(res_global$eig), rowNames = T)
addWorksheet(wb, "Var_Coords"); writeData(wb, "Var_Coords", as.data.frame(res_global$var$coord), rowNames = T)
if(length(unique(df_clr$Group)) > 1) {
  dist_mat <- dist(df_clr[, global_markers], method = "euclidean")
  adonis_res <- adonis2(dist_mat ~ Group, data = df_clr, permutations = 999)
  addWorksheet(wb, "PERMANOVA"); writeData(wb, "PERMANOVA", as.data.frame(adonis_res), rowNames = T)
}
saveWorkbook(wb, file.path(CONFIG$base_output_dir, "Detailed_Stats_Global.xlsx"), overwrite = TRUE)


# --- B. SENSITIVITY ANALYSIS (Leave-One-Out) ---
cat("\n[2] Running SENSITIVITY ANALYSIS (Exclude One Group)...\n")

for(group_name in names(marker_groups)) {
  
  vars_to_exclude <- marker_groups[[group_name]]
  # Global minus Group
  vars_to_keep <- setdiff(global_markers, vars_to_exclude)
  
  cat(sprintf("   -> Processing: Without '%s' (%d markers)\n", group_name, length(vars_to_keep)))
  
  file_name <- paste0("PCA_without_", group_name)
  title_suffix <- paste0("Excluding ", group_name)
  
  run_pca_analysis(df_clr, vars_to_keep, file_name, title_suffix, path_sensitivity)
}


# --- C. FOCUSED ANALYSIS (Single Group Only) ---
cat("\n[3] Running FOCUSED ANALYSIS (Single Group Only)...\n")

for(group_name in names(marker_groups)) {
  
  # Only Group (No spurious cols, no other groups)
  vars_to_keep <- marker_groups[[group_name]]
  
  cat(sprintf("   -> Processing: Only '%s' (%d markers)\n", group_name, length(vars_to_keep)))
  
  file_name <- paste0("PCA_only_", group_name)
  title_suffix <- paste0("Only ", group_name)
  
  run_pca_analysis(df_clr, vars_to_keep, file_name, title_suffix, path_focused)
}

# 6. MARKER BOXPLOTS (Global) --------------------------------------------------
cat("\n[4] Generating Marker Boxplots (Global)...\n")
df_raw <- DATA$raw_imputed
long_df <- df_raw %>%
  dplyr::select(Patient_ID, Group, all_of(global_markers)) %>%
  pivot_longer(cols = -c(Patient_ID, Group), names_to = "Marker", values_to = "Value")

p_box <- ggplot(long_df, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~Marker, scales = "free_y") +
  stat_compare_means(label = "p.signif", label.y.npc = "top") + 
  scale_fill_manual(values = CONFIG$colors) +
  theme_bw() +
  labs(y = "% Frequency (Raw)", title = "Marker Expression Levels")

ggsave(file.path(CONFIG$base_output_dir, "Marker_Boxplots.pdf"), p_box, width = 16, height = 12)

cat("\n=== ANALYSIS COMPLETE ===\n")