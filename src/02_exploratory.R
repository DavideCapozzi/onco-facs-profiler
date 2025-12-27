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

# MAPPING DATA STRUCTURE (from Step 1)
# ------------------------------------------------------------------------------
df_global <- DATA$hybrid_data_z 
ilr_list <- DATA$ilr_balances 
metadata  <- DATA$metadata
my_colors <- get_palette(config)
safe_markers <- DATA$hybrid_markers
raw_matrix <- DATA$raw_matrix

out_dir <- file.path(config$output_root, "02_exploratory")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(sprintf("[Data] Global Matrix: %d Samples x %d Markers (Z-Scored)", 
                nrow(df_global), ncol(df_global) - length(meta_cols)))

# ==============================================================================
# PART A: RAW DATA DIAGNOSTICS (Unified Distribution Checks)
# ==============================================================================
message("[Diagnostics] Generating unified distribution plots for RAW data...")

# 1. Prepare Raw Data with Metadata
if (!all(rownames(raw_matrix) %in% metadata$Patient_ID)) {
  warning("[Diagnostics] Mismatch between Raw Matrix IDs and Metadata IDs.")
}

meta_ordered <- metadata[match(rownames(raw_matrix), metadata$Patient_ID), ]
df_raw_viz <- cbind(meta_ordered, as.data.frame(raw_matrix))

# --- LOGIC CHANGE: Merge Tumor Groups for Visualization ---
# We create a temporary dataframe for plotting to preserve original data for stats.
df_viz_merged <- df_raw_viz
control_grp <- config$control_group  # e.g., "Healthy"

# If Group is NOT Control, rename to "Tumor"
# (Assumes config$control_group matches exactly the string in metadata)
df_viz_merged$Group <- ifelse(
  df_viz_merged$Group == control_grp, 
  control_grp, 
  "Tumor"
)

# Define a visualization palette specifically for this 2-group view
# We use the Healthy color from config, and a generic color (e.g., SteelBlue) for Tumor
viz_colors <- c()
viz_colors[[control_grp]] <- my_colors[[control_grp]] 
viz_colors[["Tumor"]] <- "#4682B4" # Manual assignment or pick from config

# Handle case where control_group might be named differently in config vs data
if (is.null(viz_colors[[control_grp]])) viz_colors[[control_grp]] <- "green4"

message(sprintf("   -> Groups merged into: %s", paste(names(viz_colors), collapse = ", ")))

# 2. Setup Output - Single PDF
raw_pdf_path <- file.path(out_dir, "Distributions_Raw_Merged_Stats.pdf")

# 3. Generate Plots
message(sprintf("   -> Saving plots to: %s", raw_pdf_path))

pdf(raw_pdf_path, width = 8, height = 6) 

# Loop over ALL markers
for (marker in DATA$markers) {
  
  if (marker %in% names(df_viz_merged)) {
    
    p <- plot_raw_distribution_merged(
      data_df = df_viz_merged, # Passing the merged dataframe
      marker_name = marker,
      colors = viz_colors,     # Passing the simplified palette
      highlight_pattern = "_LS",   
      highlight_color = "#FFD700"  
    )
    
    if (!is.null(p)) print(p)
  }
}

dev.off()
message("[Diagnostics] Raw distribution analysis complete.\n")

# ==============================================================================
# PART B1: PCA (Global Hybrid View)
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
# PART B2: STRATIFICATION HEATMAP (ComplexHeatmap)
# ==============================================================================
message("[Viz] Generating Stratification Heatmap...")

# 1. Prepare Data
# We rely on 'df_global' (Z-scored) and 'safe_markers' defined earlier
mat_heatmap_z <- as.matrix(df_global[, safe_markers])
rownames(mat_heatmap_z) <- df_global$Patient_ID

# Metadata aligned
meta_heatmap <- df_global[, meta_cols]

# 2. Output PDF
pdf_heat_path <- file.path(out_dir, "Heatmap_Stratification_Clustered.pdf")

# Note: ComplexHeatmap draws directly to the device
pdf(pdf_heat_path, width = 10, height = 8)

tryCatch({
  ht_obj <- plot_stratification_heatmap(
    mat_z = mat_heatmap_z,
    metadata = meta_heatmap,
    group_colors = my_colors,
    title = "Global Clustering (Hybrid Z-Score)"
  )
  draw(ht_obj, merge_legend = TRUE)
  
}, error = function(e) {
  warning(paste("Heatmap generation failed:", e$message))
  plot.new()
  text(0.5, 0.5, "Heatmap Failed (Check ComplexHeatmap installation)")
})

dev.off()
message(sprintf("   -> Heatmap saved to: %s", pdf_heat_path))

# ==============================================================================
# PART C: STATISTICAL TESTING (PERMANOVA)
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
# PART D: SAVE OUTPUT
# ==============================================================================

excel_path <- file.path(out_dir, "Statistical_Test_Results.xlsx")
saveWorkbook(wb, excel_path, overwrite = TRUE)

message(sprintf("\n[Output] Results saved to: %s", out_dir))
message("=== STEP 2 COMPLETE ===\n")