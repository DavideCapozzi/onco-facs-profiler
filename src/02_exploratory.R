# src/02_exploratory.R
# ==============================================================================
# STEP 02: EXPLORATORY DATA ANALYSIS & STATISTICAL ASSUMPTION CHECK
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(FactoMineR)
  library(factoextra)
  library(openxlsx)
  library(vegan)
})

# Source modules
source("R/utils_io.R")          
source("R/modules_hypothesis.R") 
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
# PART B: HYPOTHESIS TESTING (PERMANOVA)
# ==============================================================================

# 1. Robust PERMANOVA
message("[Stats] Running Robust PERMANOVA (adonis2)...")
# UPDATED FUNCTION NAME: test_coda_permanova
permanova_res <- test_coda_permanova(df_clr, group_col = "Group", n_perm = config$stats$n_perm)

perm_pval <- permanova_res$`Pr(>F)`[1]
perm_r2   <- permanova_res$R2[1]

message(sprintf("   -> PERMANOVA (Robust):  p = %.5f (R2 = %.2f%%)", perm_pval, perm_r2 * 100))

# ==============================================================================
# PART C: REPORTING
# ==============================================================================

wb <- createWorkbook()

# Salviamo solo PERMANOVA
addWorksheet(wb, "PERMANOVA_Robust")
writeData(wb, "PERMANOVA_Robust", as.data.frame(permanova_res), rowNames = TRUE)

excel_path <- file.path(out_dir, "Statistical_Test_Results.xlsx")
saveWorkbook(wb, excel_path, overwrite = TRUE)

message(sprintf("[Output] All exploratory and statistical results saved to: %s", out_dir))
message("=== STEP 2 COMPLETE ===\n")