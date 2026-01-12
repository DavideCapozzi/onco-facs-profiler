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
  library(ComplexHeatmap)
  library(mixOmics)
})

source("R/utils_io.R")          
source("R/modules_hypothesis.R") 
source("R/modules_viz.R")
source("R/modules_multivariate.R")   
source("R/modules_interpretation.R") 

message("\n=== PIPELINE STEP 2: EXPLORATORY & STATS (Hybrid/Z-Score) ===")

# 1. Load Config & Data
# ------------------------------------------------------------------------------
config <- load_config("config/global_params.yml")
input_file <- file.path(config$output_root, "01_QC", "data_processed.rds")

if (!file.exists(input_file)) stop("Step 01 output not found. Run src/01_ingest.R first.")

DATA <- readRDS(input_file)

# Mapping data structure
df_global <- DATA$hybrid_data_z 
ilr_list <- DATA$ilr_balances 
metadata  <- DATA$metadata
meta_cols <- colnames(metadata)
my_colors <- get_palette(config)
safe_markers <- DATA$hybrid_markers
raw_matrix <- DATA$raw_matrix

out_dir <- file.path(config$output_root, "02_exploratory")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(sprintf("[Data] Global Matrix: %d Samples x %d Markers (Z-Scored)", 
                nrow(df_global), ncol(df_global) - length(meta_cols)))

# 2. Data diagnostics (Refactored)
# ------------------------------------------------------------------------------
message("[Diagnostics] Generating unified distribution plots for RAW data...")

meta_ordered <- metadata[match(rownames(raw_matrix), metadata$Patient_ID), ]
df_raw_viz <- cbind(meta_ordered, as.data.frame(raw_matrix))

# Visualization merging logic
df_viz_merged <- df_raw_viz
viz_colors <- my_colors
ctrl_grp <- config$control_group 
should_merge <- if (!is.null(config$viz$merge_case_groups)) config$viz$merge_case_groups else FALSE

if (should_merge) {
  case_label_str <- paste(config$case_groups, collapse = " + ")
  df_viz_merged$Plot_Label <- ifelse(df_viz_merged$Group == ctrl_grp, ctrl_grp, case_label_str)
  df_viz_merged$Group <- df_viz_merged$Plot_Label
  
  c_case <- if (!is.null(config$colors$Case)) config$colors$Case else config$colors$cases[1]
  viz_colors <- character()
  viz_colors[[ctrl_grp]] <- config$colors$control
  viz_colors[[case_label_str]] <- c_case
}

# Prepare Distribution Colors from Config
dist_cols <- c(
  mean = if(!is.null(config$colors$distribution$mean)) config$colors$distribution$mean else "darkred",
  median = if(!is.null(config$colors$distribution$median)) config$colors$distribution$median else "black"
)

# Output Distribution PDF
raw_pdf_path <- file.path(out_dir, "Distributions_Raw_Merged_Stats.pdf")
hl_pattern <- if(!is.null(config$viz$highlight_pattern)) config$viz$highlight_pattern else ""

viz_save_distribution_report(
  data_df = df_viz_merged,
  markers = DATA$markers,
  file_path = raw_pdf_path,
  colors = viz_colors,
  hl_pattern = hl_pattern,
  stat_colors = dist_cols
)
message("[Diagnostics] Raw distribution analysis complete.\n")

# 3. PCA
# ------------------------------------------------------------------------------
message("[PCA] Running PCA on Global Z-Scored Matrix...")

pca_input <- df_global[, safe_markers]
res_pca <- PCA(pca_input, scale.unit = FALSE, graph = FALSE)
highlight_map <- c("_LS" = "#FFD700") # Could also be moved to config if needed
SHOW_LABELS_PCA <- TRUE

# A. Individuals PCA plots
pdf_ind_path <- file.path(out_dir, "PCA_Global_Individuals_MultiPage.pdf")
pdf(pdf_ind_path, width = 8, height = 6)
print(plot_pca_custom(res_pca, df_global, my_colors, dims = c(1, 2), show_labels = SHOW_LABELS_PCA, highlight_patterns = highlight_map))
print(plot_pca_custom(res_pca, df_global, my_colors, dims = c(1, 3), show_labels = SHOW_LABELS_PCA, highlight_patterns = highlight_map))
print(plot_pca_custom(res_pca, df_global, my_colors, dims = c(2, 3), show_labels = SHOW_LABELS_PCA, highlight_patterns = highlight_map))
dev.off() 

# B. Variables/Loadings plots
pdf_var_path <- file.path(out_dir, "PCA_Global_Variables_MultiPage.pdf")
pdf(pdf_var_path, width = 8, height = 6)
plot_vars_dims <- function(pca_res, axes_vec) {
  fviz_pca_var(pca_res, axes = axes_vec, col.var = "contrib", 
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE) + 
    labs(title = sprintf("Variables - PCA (PC%d vs PC%d)", axes_vec[1], axes_vec[2])) +
    theme_coda()
}
print(plot_vars_dims(res_pca, c(1, 2)))
print(plot_vars_dims(res_pca, c(1, 3)))
print(plot_vars_dims(res_pca, c(2, 3)))
dev.off()

# C. Scree Plot
pdf(file.path(out_dir, "PCA_Global_Scree.pdf"), width = 6, height = 4)
print(fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 50)) + theme_minimal() + ggtitle("Scree Plot (Global)"))
dev.off()

# 4. Patient Stratification Heatmap
# ------------------------------------------------------------------------------
message("[Viz] Generating Patient Stratification Heatmap...")

mat_heatmap_z <- as.matrix(df_global[, safe_markers])
rownames(mat_heatmap_z) <- df_global$Patient_ID

target_cols <- intersect(if(!is.null(config$viz$heatmap_metadata)) config$viz$heatmap_metadata else "Group", names(df_global))
meta_heatmap <- df_global[, target_cols, drop = FALSE]

annotation_colors <- viz_generate_complex_heatmap_colors(
  metadata = meta_heatmap,
  group_col = "Group",
  base_colors = my_colors,
  hl_pattern = config$viz$highlight_pattern,
  hl_color = config$colors$highlight
)

# Prepare Heatmap Gradient Colors from Config
hm_cols <- c(
  low = if(!is.null(config$colors$heatmap$low)) config$colors$heatmap$low else "#2166AC",
  mid = if(!is.null(config$colors$heatmap$mid)) config$colors$heatmap$mid else "#F7F7F7",
  high = if(!is.null(config$colors$heatmap$high)) config$colors$heatmap$high else "#B2182B"
)

pdf_heat_path <- file.path(out_dir, "Heatmap_Stratification_Clustered.pdf")
pdf(pdf_heat_path, width = 10, height = 8)

tryCatch({
  ht_obj <- plot_stratification_heatmap(
    mat_z = mat_heatmap_z,
    metadata = meta_heatmap,
    annotation_colors_list = annotation_colors,
    title = "Global Clustering (Hybrid Z-Score)",
    gradient_colors = hm_cols
  )
  draw(ht_obj, merge_legend = TRUE)
}, error = function(e) {
  warning(paste("Heatmap generation failed:", e$message))
  plot.new(); text(0.5, 0.5, paste("Heatmap Failed:", e$message), cex = 0.8)
})
dev.off()
message(sprintf("   -> Heatmap saved to: %s", pdf_heat_path))

# 5. Statistical Testing & Advanced Interpretation
# ------------------------------------------------------------------------------
wb <- createWorkbook()

message("\n[Stats] Running Global PERMANOVA & Driver Analysis...")

# 1. PERMANOVA
perm_global <- test_coda_permanova(
  data_input = df_global[, c("Group", safe_markers)], 
  group_col = "Group", 
  n_perm = config$stats$n_perm
)
addWorksheet(wb, "Global_PERMANOVA")
writeData(wb, "Global_PERMANOVA", as.data.frame(perm_global), rowNames = TRUE)

# 2. sPLS-DA Driver Analysis (Sparse)
run_pls <- if(!is.null(config$multivariate$run_plsda)) config$multivariate$run_plsda else FALSE

if (run_pls) {
  message("   [sPLS-DA] Identifying Global Immunological Signature (Sparse)...")
  
  tryCatch({
    X_pls <- df_global[, safe_markers]
    meta_pls <- df_global[, meta_cols, drop=FALSE]
    
    # Run Model
    pls_res <- run_splsda_model(X_pls, meta_pls, group_col = "Group", 
                                n_comp = config$multivariate$n_comp,
                                folds = config$multivariate$validation_folds)
    
    # Extract Metrics & Loadings
    top_drivers <- extract_plsda_loadings(pls_res)
    perf_metrics <- extract_plsda_performance(pls_res)
    tuning_info  <- as.data.frame(pls_res$tuning$choice.keepX)
    
    # Save Tables
    addWorksheet(wb, "Global_sPLSDA_Drivers")
    writeData(wb, "Global_sPLSDA_Drivers", top_drivers)
    addWorksheet(wb, "Global_sPLSDA_Quality")
    writeData(wb, "Global_sPLSDA_Quality", perf_metrics)
    writeData(wb, "Global_sPLSDA_Quality", tuning_info, startRow = 6, startCol = 1)
    
    # Log Success
    ber_msg <- if(is.na(perf_metrics$Overall_BER[1])) "N/A" else sprintf("%.2f%%", perf_metrics$Overall_BER[1] * 100)
    message(sprintf("   -> sPLS-DA Selected: %d markers (Comp1). BER: %s", 
                    pls_res$tuning$choice.keepX[1], ber_msg))
    
    # Prepare Loadings Colors
    pls_cols <- c(
      positive = if(!is.null(config$colors$plsda$positive)) config$colors$plsda$positive else "#CD5C5C",
      negative = if(!is.null(config$colors$plsda$negative)) config$colors$plsda$negative else "#4682B4"
    )
    
    # Save Loadings Plot (Multi-Page using New Wrapper)
    # Note: tuning_info is passed. The wrapper handles extraction safely now.
    viz_save_plsda_loadings(
      loadings_df = top_drivers,
      n_comp = pls_res$model$ncomp,
      tuning_df = tuning_info,
      file_path = file.path(out_dir, "Global_sPLSDA_Loadings.pdf"),
      bar_colors = pls_cols
    )
    
    # Save Biplot (Previously failed because viz_save_plsda_loadings crashed)
    pdf(file.path(out_dir, "Global_sPLSDA_Biplot.pdf"), width = 8, height = 7)
    mixOmics::plotIndiv(pls_res$model, comp = c(1,2), group = meta_pls$Group, 
                        ellipse = TRUE, legend = TRUE, title = "sPLS-DA: Sparse Signature",
                        star.plot = TRUE) 
    dev.off()
    message("   -> sPLS-DA Biplot saved.")
    
  }, error = function(e) {
    message(paste("   [ERROR] sPLS-DA Execution Failed:", e$message))
  })
}

# B. Local PERMANOVA + ILR Decoding
if (length(ilr_list) > 0) {
  message("\n[Stats] Running Local Analysis on Compositional Groups...")
  summary_local <- data.frame()
  
  for (grp_name in names(ilr_list)) {
    mat_ilr <- ilr_list[[grp_name]]
    
    res_perm <- test_coda_permanova(cbind(metadata[, "Group", drop=FALSE], as.data.frame(mat_ilr)), group_col = "Group", n_perm = config$stats$n_perm)
    pval <- res_perm$`Pr(>F)`[1]
    
    summary_local <- rbind(summary_local, data.frame(SubGroup = grp_name, P_Value = pval, R2_Percent = res_perm$R2[1] * 100))
    sheet_name <- substr(paste0("ILR_", grp_name), 1, 31) 
    addWorksheet(wb, sheet_name); writeData(wb, sheet_name, as.data.frame(res_perm), rowNames = TRUE)
    
    if (pval < 0.1) {
      target_mks <- config$hybrid_groups[[grp_name]]
      valid_mks <- intersect(target_mks, colnames(raw_matrix))
      if(length(valid_mks) > 1) {
        common_ids <- intersect(rownames(mat_ilr), rownames(raw_matrix))
        if (length(common_ids) >= 3) {
          decoding_table <- decode_ilr_to_clr(mat_ilr[common_ids,,drop=FALSE], raw_matrix[common_ids, valid_mks, drop=FALSE])
          if (nrow(decoding_table) > 0) {
            ds_name <- substr(paste0("Dec_", grp_name), 1, 31)
            addWorksheet(wb, ds_name); writeData(wb, ds_name, decoding_table)
          }
        }
      }
    }
  }
  addWorksheet(wb, "Summary_Local_Tests"); writeData(wb, "Summary_Local_Tests", summary_local)
}

# 6. Save Output
# ------------------------------------------------------------------------------
excel_path <- file.path(out_dir, "Statistical_Test_Results.xlsx")
saveWorkbook(wb, excel_path, overwrite = TRUE)
message(sprintf("\n[Output] Results saved to: %s", out_dir))
message("=== STEP 2 COMPLETE ===\n")