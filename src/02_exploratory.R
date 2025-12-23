# src/02_exploratory.R
# ==============================================================================
# STEP 02: EXPLORATORY DATA ANALYSIS & ROBUST HYPOTHESIS TESTING
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

message("\n=== PIPELINE STEP 2: EXPLORATORY & STATS (Hybrid/Z-Score) ===")

# 1. Load Config & Data
config <- load_config("config/global_params.yml")
input_file <- file.path(config$output_root, "01_QC", "data_processed.rds")

if (!file.exists(input_file)) stop("Step 01 output not found. Run src/01_ingest.R first.")

DATA <- readRDS(input_file)

# MAPPING NEW DATA STRUCTURE (from Step 1)
# ------------------------------------------------------------------------------
df_global <- DATA$hybrid_data_z 
ilr_list <- DATA$ilr_balances 
metadata  <- DATA$metadata
my_colors <- get_palette(config)
safe_markers <- DATA$hybrid_markers

out_dir <- file.path(config$output_root, "02_exploratory")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(sprintf("[Data] Global Matrix: %d Samples x %d Markers (Z-Scored)", 
                nrow(df_global), ncol(df_global) - length(meta_cols)))

# ==============================================================================
# PART A: PCA (Global Hybrid View)
# ==============================================================================
message("[PCA] Running PCA on Global Z-Scored Matrix...")

# Prepare Matrix (Exclude metadata)
pca_input <- df_global[, safe_markers]

# Run PCA (FactoMineR)
res_pca <- PCA(pca_input, scale.unit = FALSE, graph = FALSE)

highlight_map <- c("_LS" = "#FFD700")

# --- 1. INDIVIDUALS PCA PLOT ---
message("[PCA] Saving Individuals plots to single PDF...")

pdf_ind_path <- file.path(out_dir, "PCA_Global_Individuals_MultiPage.pdf")
pdf(pdf_ind_path, width = 8, height = 6)

# Configurazione Label
SHOW_LABELS_PCA <- TRUE 

# Page 1: PC1 vs PC2
print(plot_pca_custom(res_pca, df_global, my_colors, dims = c(1, 2), show_labels = SHOW_LABELS_PCA, highlight_patterns = highlight_map))

# Page 2: PC1 vs PC3
print(plot_pca_custom(res_pca, df_global, my_colors, dims = c(1, 3), show_labels = SHOW_LABELS_PCA, highlight_patterns = highlight_map))

# Page 3: PC2 vs PC3
print(plot_pca_custom(res_pca, df_global, my_colors, dims = c(2, 3), show_labels = SHOW_LABELS_PCA, highlight_patterns = highlight_map))

dev.off() 


# --- 2. VARIABLES/LOADINGS PLOT  ---
message("[PCA] Saving Variables/Loadings plots to single PDF...")

pdf_var_path <- file.path(out_dir, "PCA_Global_Variables_MultiPage.pdf")
pdf(pdf_var_path, width = 8, height = 6)

plot_vars_dims <- function(pca_res, axes_vec) {
  fviz_pca_var(pca_res, 
               axes = axes_vec,
               col.var = "contrib", 
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE) + 
    labs(title = sprintf("Variables - PCA (PC%d vs PC%d)", axes_vec[1], axes_vec[2])) +
    theme_coda()
}

# Page 1: PC1 vs PC2
print(plot_vars_dims(res_pca, c(1, 2)))

# Page 2: PC1 vs PC3
print(plot_vars_dims(res_pca, c(1, 3)))

# Page 3: PC2 vs PC3
print(plot_vars_dims(res_pca, c(2, 3)))

dev.off()


# --- 3. SCREE PLOT ---
p_scree <- fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 50)) + 
  theme_minimal() + ggtitle("Scree Plot (Global)")
ggsave(file.path(out_dir, "PCA_Global_Scree.pdf"), p_scree, width = 6, height = 4)

# ==============================================================================
# PART B: STATISTICAL TESTING (PERMANOVA)
# ==============================================================================
wb <- createWorkbook()

# 1. GLOBAL PERMANOVA
# ------------------------------------------------------------------------------
message("\n[Stats] Running Global PERMANOVA (All Markers)...")

# [ROBUST FIX] Construct input data using the whitelist.
# We explicitly select ONLY the Group column and the validated numeric markers.
# This automatically drops 'Patient_ID', 'Original_Source', or any other metadata.
df_stats_global <- df_global[, c("Group", safe_markers)]

# Now we run the test. 
# metadata_cols is set to NULL because we already cleaned the data above.
perm_global <- test_coda_permanova(
  data_input = df_stats_global, 
  group_col = "Group", 
  metadata_cols = NULL, 
  n_perm = config$stats$n_perm
)

message(sprintf("   -> Global p-value: %.5f (R2: %.2f%%)", 
                perm_global$`Pr(>F)`[1], perm_global$R2[1] * 100))

addWorksheet(wb, "Global_PERMANOVA")
writeData(wb, "Global_PERMANOVA", as.data.frame(perm_global), rowNames = TRUE)


# 2. LOCAL PERMANOVA (Sub-groups ILR Balances)
# ------------------------------------------------------------------------------
if (length(ilr_list) > 0) {
  message("\n[Stats] Running Local PERMANOVA on Compositional Groups (ILR)...")
  
  summary_local <- data.frame()
  
  for (grp_name in names(ilr_list)) {
    
    mat_ilr <- ilr_list[[grp_name]]
    
    # [ROBUST FIX] Construct input data specifically for this group.
    # We bind the 'Group' column from metadata directly with the numeric ILR matrix.
    # This ensures no other metadata columns (like Original_Source) are included.
    df_stats_local <- cbind(
      metadata[, "Group", drop = FALSE], 
      as.data.frame(mat_ilr)
    )
    
    # Run Test
    res_perm <- test_coda_permanova(
      data_input = df_stats_local, 
      group_col = "Group", 
      metadata_cols = NULL, 
      n_perm = config$stats$n_perm
    )
    
    # Store Summary
    pval <- res_perm$`Pr(>F)`[1]
    r2   <- res_perm$R2[1]
    
    message(sprintf("   -> Group '%s': p = %.5f", grp_name, pval))
    
    summary_local <- rbind(summary_local, data.frame(
      SubGroup = grp_name,
      P_Value = pval,
      R2_Percent = r2 * 100,
      F_Model = res_perm$F[1]
    ))
    
    sheet_name <- substr(paste0("ILR_", grp_name), 1, 31) 
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, as.data.frame(res_perm), rowNames = TRUE)
  }
  
  addWorksheet(wb, "Summary_Local_Tests")
  writeData(wb, "Summary_Local_Tests", summary_local)
}

# ==============================================================================
# PART C: SAVE OUTPUT
# ==============================================================================

excel_path <- file.path(out_dir, "Statistical_Test_Results.xlsx")
saveWorkbook(wb, excel_path, overwrite = TRUE)

message(sprintf("\n[Output] Results saved to: %s", out_dir))
message("=== STEP 2 COMPLETE ===\n")