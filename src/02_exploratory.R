# src/02_exploratory.R
# ==============================================================================
# STEP 02: EXPLORATORY DATA ANALYSIS & STATISTICAL ASSUMPTION CHECK
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(FactoMineR)
  library(factoextra)
  library(openxlsx)
  library(vegan) # Essential for robust testing
})

# Source modules
source("R/infrastructure.R")
source("R/modules_stats.R")
source("R/modules_viz.R")


message("\n=== PIPELINE STEP 2: EXPLORATORY ANALYSIS & ASSUMPTION CHECKS ===")

# 1. Load Config & Data
config <- load_config("config/global_params.yml")
input_file <- file.path(config$output_root, "01_QC", "data_processed.rds")

if (!file.exists(input_file)) stop("Step 01 output not found. Run src/01_ingest.R first.")

DATA <- readRDS(input_file)
df_clr <- DATA$clr_data
df_ilr <- DATA$ilr_data
markers <- DATA$markers
my_colors <- get_palette(config)

# 2. Prepare Output Directory
out_dir <- file.path(config$output_root, "02_exploratory")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(sprintf("[Data] Analyzable Set: %d Samples x %d Markers", nrow(df_clr), length(markers)))

# ==============================================================================
# PART A: PCA (Visualization)
# ==============================================================================
message("[PCA] Running PCA on CLR transformed data...")

pca_input <- df_clr %>% select(all_of(markers)) %>% as.matrix()
res_pca <- run_coda_pca(pca_input)

# Plotting
p_pca <- plot_pca_custom(res_pca, df_clr[, c("Patient_ID", "Group")], my_colors, show_labels = TRUE)
ggsave(file.path(out_dir, "PCA_Score_Plot.pdf"), p_pca, width = 8, height = 6)

p_scree <- fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 50)) + theme_coda()
ggsave(file.path(out_dir, "PCA_Scree_Plot.pdf"), p_scree, width = 6, height = 4)

p_vars <- plot_pca_variables(res_pca)
ggsave(file.path(out_dir, "PCA_Variable_Contrib.pdf"), p_vars, width = 8, height = 6)


# ==============================================================================
# PART B: STATISTICAL ASSUMPTION CHECKING
# ==============================================================================
message("[Stats] Checking MANOVA Assumptions on ILR Data...")
# We use ILR because it is Euclidean and mathematically proper for these tests.

check_res <- check_manova_assumptions(df_ilr, group_col = "Group", metadata_cols = c("Patient_ID", "Group"))

# 1. Viz: Homogeneity
p_homog <- plot_homogeneity(check_res, my_colors)
ggsave(file.path(out_dir, "Assumption_Homogeneity.pdf"), p_homog, width = 6, height = 5)

# 2. Viz: Normality
p_mvn <- plot_mvn_check(check_res)
ggsave(file.path(out_dir, "Assumption_MVN_QQ.pdf"), p_mvn, width = 6, height = 5)

# 3. Interpretation
is_homog_ok <- check_res$homogeneity_pval > 0.05
message(sprintf("   -> Homogeneity of Dispersion: p = %.4f [%s]", 
                check_res$homogeneity_pval, 
                ifelse(is_homog_ok, "PASS", "FAIL - Groups have different spread")))

if (!is_homog_ok) {
  warning("[Stats] Homogeneity assumption violated. Standard MANOVA results may be unreliable (high Type I error).")
  warning("[Stats] Proceeding with Robust PERMANOVA.")
}

# ==============================================================================
# PART C: HYPOTHESIS TESTING (MANOVA vs PERMANOVA)
# ==============================================================================

# 1. Standard MANOVA (Pillai's Trace - most robust parametric option)
message("[Stats] Running Parametric MANOVA (Pillai)...")
manova_results <- run_coda_manova(df_ilr, df_clr, metadata_cols = c("Patient_ID", "Group"))

# 2. Robust PERMANOVA (Distance-based, Non-parametric)
message("[Stats] Running Robust PERMANOVA (adonis2)...")
permanova_res <- run_coda_permanova(df_ilr, group_col = "Group", n_perm = config$stats$n_perm)

# Extract PERMANOVA stats
perm_pval <- permanova_res$`Pr(>F)`[1]
perm_r2   <- permanova_res$R2[1]

message(sprintf("   -> MANOVA (Parametric): p = %.5f", manova_results$p_value))
message(sprintf("   -> PERMANOVA (Robust):  p = %.5f (R2 = %.2f%%)", perm_pval, perm_r2 * 100))

# 3. Visualization of Separation
p_manova <- plot_manova_results(manova_results, my_colors) +
  labs(subtitle = sprintf("MANOVA p=%.4f | PERMANOVA p=%.4f (R2=%.1f%%)", 
                          manova_results$p_value, perm_pval, perm_r2*100))
ggsave(file.path(out_dir, "Global_Separation_Boxplot.pdf"), p_manova, width = 6, height = 5)

# Loadings
p_loadings <- plot_manova_loadings(manova_results, top_n = 20)
ggsave(file.path(out_dir, "Separation_Drivers_Lollipop.pdf"), p_loadings, width = 7, height = 6)


# ==============================================================================
# PART D: REPORTING
# ==============================================================================

wb <- createWorkbook()

# Sheet 1: Assumptions
addWorksheet(wb, "Assumptions")
assumption_df <- data.frame(
  Test = c("Homogeneity of Dispersion (Betadisper)"),
  P_Value = c(check_res$homogeneity_pval),
  Result = c(ifelse(is_homog_ok, "Pass", "Fail (Unequal Variance)"))
)
writeData(wb, "Assumptions", assumption_df)

# Sheet 2: MANOVA
addWorksheet(wb, "MANOVA_Parametric")
writeData(wb, "MANOVA_Parametric", as.data.frame(manova_results$summary), rowNames = TRUE)

# Sheet 3: PERMANOVA
addWorksheet(wb, "PERMANOVA_Robust")
writeData(wb, "PERMANOVA_Robust", as.data.frame(permanova_res), rowNames = TRUE)

# Sheet 4: Loadings
addWorksheet(wb, "Marker_Loadings")
writeData(wb, "Marker_Loadings", manova_results$loadings)

excel_path <- file.path(out_dir, "Statistical_Test_Results.xlsx")
saveWorkbook(wb, excel_path, overwrite = TRUE)

message(sprintf("[Output] All exploratory and statistical results saved to: %s", out_dir))
message("=== STEP 2 COMPLETE ===\n")