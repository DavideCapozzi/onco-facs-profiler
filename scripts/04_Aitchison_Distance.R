# ==============================================================================
# SCRIPT 04: AITCHISON DISTANCE ANALYSIS
# PURPOSE: Quantify immunological similarity to Healthy state using 
#          multivariate Aitchison distance (Euclidean on CLR data)
# INPUT: clean_data.rds (from Step 1)
# OUTPUT: Distance metrics, statistical comparisons, visualizations
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(openxlsx)
  library(reshape2)
})

# 1. CONFIGURATION -------------------------------------------------------------
CONFIG <- list(
  input_file = "/home/davidec/projects/compositional_analysis/processed_data/clean_data.rds",
  output_dir = "/home/davidec/projects/compositional_analysis/results/04_Aitchison",
  
  # Reference group for centroid calculation
  reference_group = "Healthy",
  
  # Visualization
  colors = c("Healthy" = "#2E8B57", "NSCLC" = "#4682B4", "HNSCC" = "#CD5C5C")
)

if(!dir.exists(CONFIG$output_dir)) dir.create(CONFIG$output_dir, recursive = TRUE)
cat("=== PIPELINE STEP 4: AITCHISON DISTANCE ANALYSIS ===\n")

# 2. LOAD DATA -----------------------------------------------------------------
if(!file.exists(CONFIG$input_file)) stop("Step 1 output not found! Run 01_Data_Preparation.R first.")
DATA <- readRDS(CONFIG$input_file)
cat("   -> Data loaded. N Patients:", nrow(DATA$metadata), "\n")

df_clr <- DATA$clr_transformed
markers <- DATA$markers

# 3. CALCULATE AITCHISON DISTANCES ---------------------------------------------
cat("[1] Calculating Aitchison Distances from Healthy Centroid...\n")

# A. Compute Healthy Centroid (multivariate mean in CLR space)
healthy_mask <- df_clr$Group == CONFIG$reference_group
healthy_clr <- df_clr[healthy_mask, markers]
healthy_centroid <- colMeans(healthy_clr)

cat(sprintf("   -> Healthy centroid calculated from %d samples\n", sum(healthy_mask)))

# B. Calculate Euclidean distance (= Aitchison distance on CLR data)
calculate_aitchison_dist <- function(patient_profile, centroid) {
  sqrt(sum((patient_profile - centroid)^2))
}

distances <- apply(df_clr[, markers], 1, function(row) {
  calculate_aitchison_dist(row, healthy_centroid)
})

# C. Create results table
results_df <- data.frame(
  Patient_ID = df_clr$Patient_ID,
  Group = df_clr$Group,
  Aitchison_Distance = distances,
  stringsAsFactors = FALSE
)

# 4. STATISTICAL COMPARISON ----------------------------------------------------
cat("[2] Performing Statistical Comparisons...\n")

# Summary statistics by group
summary_stats <- results_df %>%
  group_by(Group) %>%
  summarise(
    N = n(),
    Mean_Distance = mean(Aitchison_Distance),
    Median_Distance = median(Aitchison_Distance),
    SD_Distance = sd(Aitchison_Distance),
    Min_Distance = min(Aitchison_Distance),
    Max_Distance = max(Aitchison_Distance),
    .groups = 'drop'
  )

print(summary_stats)

# Pairwise Wilcoxon tests vs Healthy (if multiple groups exist)
if(length(unique(results_df$Group)) > 1) {
  pairwise_tests <- results_df %>%
    filter(Group != CONFIG$reference_group) %>%
    group_by(Group) %>%
    summarise(
      vs_Healthy_P = wilcox.test(
        Aitchison_Distance, 
        results_df$Aitchison_Distance[results_df$Group == CONFIG$reference_group],
        exact = FALSE
      )$p.value,
      .groups = 'drop'
    )
  
  print(pairwise_tests)
}

# 5. VISUALIZATIONS ------------------------------------------------------------
cat("[3] Generating Visualizations...\n")

# A. Boxplot with statistical comparisons
p_box <- ggplot(results_df, aes(x = Group, y = Aitchison_Distance, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.7) +
  stat_compare_means(
    method = "wilcox.test", 
    ref.group = CONFIG$reference_group,
    label = "p.format",
    label.y.npc = "top"
  ) +
  scale_fill_manual(values = CONFIG$colors) +
  theme_bw(base_size = 12) +
  labs(
    title = "Immunological Divergence from Healthy State",
    subtitle = "Metric: Aitchison Distance (Euclidean on CLR-transformed data)",
    y = "Aitchison Distance",
    x = "",
    caption = sprintf("Reference: %s Centroid (N=%d)", 
                      CONFIG$reference_group, sum(healthy_mask))
  ) +
  theme(legend.position = "none")

ggsave(file.path(CONFIG$output_dir, "Aitchison_Distance_Boxplot.pdf"), 
       p_box, width = 7, height = 6)

# B. Density plot (distribution comparison)
p_density <- ggplot(results_df, aes(x = Aitchison_Distance, fill = Group)) +
  geom_density(alpha = 0.5, color = NA) +
  geom_vline(data = summary_stats, aes(xintercept = Median_Distance, color = Group),
             linetype = "dashed", size = 1) +
  scale_fill_manual(values = CONFIG$colors) +
  scale_color_manual(values = CONFIG$colors) +
  theme_bw(base_size = 12) +
  labs(
    title = "Distribution of Aitchison Distances",
    subtitle = "Dashed lines indicate group medians",
    x = "Aitchison Distance from Healthy Centroid",
    y = "Density"
  )

ggsave(file.path(CONFIG$output_dir, "Aitchison_Distance_Density.pdf"), 
       p_density, width = 8, height = 5)

# C. Individual patient plot (ordered by distance)
results_df_sorted <- results_df %>%
  arrange(Aitchison_Distance) %>%
  mutate(Rank = row_number())

p_individual <- ggplot(results_df_sorted, 
                       aes(x = Rank, y = Aitchison_Distance, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = median(results_df$Aitchison_Distance[healthy_mask]),
             linetype = "dashed", color = CONFIG$colors["Healthy"], size = 1) +
  scale_color_manual(values = CONFIG$colors) +
  theme_bw(base_size = 12) +
  labs(
    title = "Patient Proximity to Healthy Immunological State",
    subtitle = "Each point represents one patient (ranked by distance)",
    x = "Patient Rank (ordered by distance)",
    y = "Aitchison Distance",
    caption = "Dashed line: Median Healthy distance (intra-group dispersion)"
  )

ggsave(file.path(CONFIG$output_dir, "Aitchison_Distance_Individual.pdf"), 
       p_individual, width = 10, height = 6)

# 6. EXPORT RESULTS ------------------------------------------------------------
cat("[4] Exporting Results...\n")

# Create Excel workbook with multiple sheets
wb <- createWorkbook()

# Sheet 1: Individual distances
addWorksheet(wb, "Patient_Distances")
writeData(wb, "Patient_Distances", results_df)

# Sheet 2: Summary statistics
addWorksheet(wb, "Summary_Statistics")
writeData(wb, "Summary_Statistics", summary_stats)

# Sheet 3: Pairwise tests (if applicable)
if(exists("pairwise_tests")) {
  addWorksheet(wb, "Pairwise_Tests")
  writeData(wb, "Pairwise_Tests", pairwise_tests)
}

saveWorkbook(wb, file.path(CONFIG$output_dir, "Aitchison_Distance_Results.xlsx"), 
             overwrite = TRUE)

# Also save CSV for easy import
write.csv(results_df, 
          file.path(CONFIG$output_dir, "Patient_Distances.csv"), 
          row.names = FALSE)

# 7. INTERPRETATION REPORT -----------------------------------------------------
cat("\n=== DISTANCE ANALYSIS SUMMARY ===\n")
cat(sprintf("Reference Group: %s (N=%d)\n", CONFIG$reference_group, sum(healthy_mask)))
cat("\nGroup Medians:\n")
for(i in 1:nrow(summary_stats)) {
  cat(sprintf("  %s: %.3f (SD=%.3f)\n", 
              summary_stats$Group[i], 
              summary_stats$Median_Distance[i],
              summary_stats$SD_Distance[i]))
}

# Biological interpretation helper
cat("\n[INTERPRETATION GUIDE]\n")
cat("- Lower distance = Greater similarity to Healthy immune profile\n")
cat("- Healthy group distance reflects intra-group variability (baseline)\n")
cat("- Long Survivors with distances near Healthy median suggest 'Healthy-like' state\n")
cat("- Statistical significance tests compare group distributions\n")

cat(sprintf("\n[SUCCESS] All results saved to: %s\n", CONFIG$output_dir))
cat("=== PIPELINE STEP 4 COMPLETED ===\n")