# ==============================================================================
# SCRIPT 02: EXPLORATORY & STATISTICAL ANALYSIS
# PURPOSE: PCA, Batch Analysis, and Marker Boxplots
# INPUT: clean_data.rds (from Step 1)
# OUTPUT: PCA plots, Excel Statistics
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
  output_dir = "/home/davidec/projects/compositional_analysis/results_fixed/02_Stats",
  # Added color for HNSCC
  colors = c("Healthy" = "#2E8B57", "NSCLC" = "#4682B4", "HNSCC" = "#CD5C5C") 
)

if(!dir.exists(CONFIG$output_dir)) dir.create(CONFIG$output_dir, recursive = TRUE)
cat("=== PIPELINE STEP 2: STATISTICAL ANALYSIS ===\n")

# 2. LOAD DATA -----------------------------------------------------------------
if(!file.exists(CONFIG$input_file)) stop("Step 1 output not found! Run 01_Data_Preparation.R first.")
DATA <- readRDS(CONFIG$input_file)
cat("   -> Data loaded. N Patients:", nrow(DATA$metadata), "\n")
cat("   -> Groups present:", paste(unique(DATA$metadata$Group), collapse=", "), "\n")

# 3. PCA & BATCH ANALYSIS ------------------------------------------------------
cat("[1] Running PCA & PERMANOVA...\n")

df_clr <- DATA$clr_transformed
markers <- DATA$markers

# A. PCA
pca_res <- PCA(df_clr[, markers], graph = FALSE, scale.unit = FALSE) 
pca_ind <- data.frame(pca_res$ind$coord, Group = df_clr$Group, ID = df_clr$Patient_ID)

# Plot (Dynamic colors based on groups present)
p_pca <- ggplot(pca_ind, aes(x = Dim.1, y = Dim.2, color = Group, fill = Group)) +
  geom_point(size = 4, alpha = 0.8) +
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +
  geom_text_repel(aes(label = ID), size = 3, show.legend = FALSE, max.overlaps = 20) +
  scale_color_manual(values = CONFIG$colors) +
  scale_fill_manual(values = CONFIG$colors) +
  labs(title = "PCA: Immunological Landscape",
       subtitle = sprintf("Dim1: %.1f%% | Dim2: %.1f%%", pca_res$eig[1,2], pca_res$eig[2,2]),
       caption = "Method: CLR Transform -> PCA") +
  theme_bw()

ggsave(file.path(CONFIG$output_dir, "PCA_Plot.pdf"), p_pca, width = 11, height = 7)

# B. PERMANOVA (Adonis)
# Only run if we have more than 1 group
if(length(unique(df_clr$Group)) > 1) {
  dist_mat <- dist(df_clr[, markers], method = "euclidean")
  adonis_res <- adonis2(dist_mat ~ Group, data = df_clr, permutations = 999)
  
  # C. BETADISPER (Homogeneity)
  beta_mod <- betadisper(dist_mat, df_clr$Group)
  
  # Export Stats
  cat("[2] Exporting Statistical Tables...\n")
  wb <- createWorkbook()
  addWorksheet(wb, "PERMANOVA")
  writeData(wb, "PERMANOVA", as.data.frame(adonis_res))
  addWorksheet(wb, "Dispersion")
  disp_res <- data.frame(Group = df_clr$Group, Dist_to_Centroid = beta_mod$distances)
  writeData(wb, "Dispersion", disp_res)
  saveWorkbook(wb, file.path(CONFIG$output_dir, "Detailed_Stats.xlsx"), overwrite = TRUE)
}

# 4. MARKER BOXPLOTS -----------------------------------------------------------
cat("[3] Generating Marker Boxplots (Raw Data)...\n")

df_raw <- DATA$raw_filtered
long_df <- df_raw %>%
  dplyr::select(Patient_ID, Group, all_of(markers)) %>%
  pivot_longer(cols = -c(Patient_ID, Group), names_to = "Marker", values_to = "Value")

p_box <- ggplot(long_df, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~Marker, scales = "free_y") +
  # Add global p-value (Kruskal-Wallis)
  stat_compare_means(label = "p.signif", label.y.npc = "top") + 
  scale_fill_manual(values = CONFIG$colors) +
  theme_bw() +
  labs(y = "% Frequency (Raw)", title = "Marker Expression Levels")

ggsave(file.path(CONFIG$output_dir, "Marker_Boxplots.pdf"), p_box, width = 16, height = 12)
cat("   -> Done! Results in 'results/02_Stats/'\n")
