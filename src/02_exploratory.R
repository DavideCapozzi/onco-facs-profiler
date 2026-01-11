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
my_colors <- get_palette(config)
safe_markers <- DATA$hybrid_markers
raw_matrix <- DATA$raw_matrix

out_dir <- file.path(config$output_root, "02_exploratory")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(sprintf("[Data] Global Matrix: %d Samples x %d Markers (Z-Scored)", 
                nrow(df_global), ncol(df_global) - length(meta_cols)))

# 2. Data diagnostics
# ------------------------------------------------------------------------------
message("[Diagnostics] Generating unified distribution plots for RAW data...")

# A. Prepare Raw Data with Metadata
if (!all(rownames(raw_matrix) %in% metadata$Patient_ID)) {
  warning("[Diagnostics] Mismatch between Raw Matrix IDs and Metadata IDs.")
}

meta_ordered <- metadata[match(rownames(raw_matrix), metadata$Patient_ID), ]
df_raw_viz <- cbind(meta_ordered, as.data.frame(raw_matrix))

# Initialize with distinct groups (default state)
df_viz_merged <- df_raw_viz
viz_colors <- my_colors
ctrl_grp <- config$control_group 

# Check configuration for merging strategy
should_merge <- if (!is.null(config$viz$merge_case_groups)) config$viz$merge_case_groups else FALSE

if (should_merge) {
  # Combine all case group names (e.g., "GroupA + GroupB")
  case_label_str <- paste(config$case_groups, collapse = " + ")
  
  # Update Group Label: Control vs Combined Cases
  df_viz_merged$Plot_Label <- ifelse(
    df_viz_merged$Group == ctrl_grp, 
    ctrl_grp, 
    case_label_str 
  )
  df_viz_merged$Group <- df_viz_merged$Plot_Label
  
  # Update Color Palette for the new labels
  # Retrieve generic 'Case' color (fallback to first case color if generic not defined)
  c_case <- if (!is.null(config$colors$Case)) config$colors$Case else config$colors$cases[1]
  
  # Reset colors to match the new merged keys
  viz_colors <- character()
  viz_colors[[ctrl_grp]] <- config$colors$control
  viz_colors[[case_label_str]] <- c_case
  
  message(sprintf("   -> Visualization: Merged View (%s vs %s)", ctrl_grp, case_label_str))
  
} else {
  # Keep original Groups and Palette
  message(sprintf("   -> Visualization: Distinct Groups (%s)", paste(names(viz_colors), collapse = ", ")))
}

# B. Output
raw_pdf_path <- file.path(out_dir, "Distributions_Raw_Merged_Stats.pdf")
message(sprintf("   -> Saving plots to: %s", raw_pdf_path))

hl_pattern <- if(!is.null(config$viz$highlight_pattern)) config$viz$highlight_pattern else ""

pdf(raw_pdf_path, width = 8, height = 6) 

for (marker in DATA$markers) {
  if (marker %in% names(df_viz_merged)) {
    p <- plot_raw_distribution_merged(
      data_df = df_viz_merged, 
      marker_name = marker,
      colors = viz_colors,     
      highlight_pattern = hl_pattern,   
      highlight_color = "#FFD700"  
    )
    if (!is.null(p)) print(p)
  }
}

dev.off()
message("[Diagnostics] Raw distribution analysis complete.\n")

# 3. PCA
# ------------------------------------------------------------------------------
message("[PCA] Running PCA on Global Z-Scored Matrix...")

# Prepare Matrix (Exclude metadata)
pca_input <- df_global[, safe_markers]

# Run PCA (FactoMineR)
res_pca <- PCA(pca_input, scale.unit = FALSE, graph = FALSE)

highlight_map <- c("_LS" = "#FFD700")

# A. Individuals PCA plots
message("[PCA] Saving Individuals plots to single PDF...")

pdf_ind_path <- file.path(out_dir, "PCA_Global_Individuals_MultiPage.pdf")
pdf(pdf_ind_path, width = 8, height = 6)

SHOW_LABELS_PCA <- TRUE 

# Page 1: PC1 vs PC2
print(plot_pca_custom(res_pca, df_global, my_colors, dims = c(1, 2), show_labels = SHOW_LABELS_PCA, highlight_patterns = highlight_map))

# Page 2: PC1 vs PC3
print(plot_pca_custom(res_pca, df_global, my_colors, dims = c(1, 3), show_labels = SHOW_LABELS_PCA, highlight_patterns = highlight_map))

# Page 3: PC2 vs PC3
print(plot_pca_custom(res_pca, df_global, my_colors, dims = c(2, 3), show_labels = SHOW_LABELS_PCA, highlight_patterns = highlight_map))

dev.off() 


# B. Variables/Loadings plots
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


# C. Scree Plot
p_scree <- fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 50)) + 
  theme_minimal() + ggtitle("Scree Plot (Global)")
pdf(file.path(out_dir, "PCA_Global_Scree.pdf"), width = 6, height = 4)
print(p_scree)
dev.off()

# 4. Patient Stratification Heatmap
# ------------------------------------------------------------------------------
message("[Viz] Generating Patient Stratification Heatmap...")

# A. Prepare Data
mat_heatmap_z <- as.matrix(df_global[, safe_markers])
rownames(mat_heatmap_z) <- df_global$Patient_ID

# Dynamic Metadata Selection
target_cols <- config$viz$heatmap_metadata
if (is.null(target_cols)) target_cols <- c("Group")

valid_cols <- intersect(target_cols, names(df_global))
if (length(valid_cols) == 0) valid_cols <- c("Group")

meta_heatmap <- df_global[, valid_cols, drop = FALSE]

# B. Smart Color Generation
annotation_colors <- list()

# Group Colors (Base Palette)
present_groups <- unique(meta_heatmap$Group)
group_pal <- my_colors[intersect(names(my_colors), present_groups)]
annotation_colors[["Group"]] <- group_pal

# Extra Metadata Colors (Inheritance or Highlight Logic)
extra_cols <- setdiff(valid_cols, "Group")
hl_pattern <- config$viz$highlight_pattern
hl_color   <- config$colors$highlight

for (col_name in extra_cols) {
  
  # Filter out NAs to prevent errors
  raw_vals <- as.character(meta_heatmap[[col_name]])
  vals <- sort(unique(raw_vals[!is.na(raw_vals)]))
  
  # [FIX] Initialize as character() to force Named Vector, NOT List
  col_pal <- character()
  
  for (v in vals) {
    # Priority 1: Exact match to a Group Name -> Inherit Group Color
    if (v %in% names(group_pal)) {
      col_pal[[v]] <- group_pal[[v]]
      
      # Priority 2: Contains Highlight Pattern (e.g., "_LS") -> Use Highlight Color
    } else if (!is.null(hl_pattern) && grepl(hl_pattern, v)) {
      col_pal[[v]] <- hl_color
      
      # Priority 3: Subgroup containing Parent Group Name -> Inherit Parent Color
    } else {
      parent_match <- NA
      for (grp in names(group_pal)) {
        if (grepl(grp, v)) {
          parent_match <- grp
          break
        }
      }
      
      if (!is.na(parent_match)) {
        col_pal[[v]] <- group_pal[[parent_match]]
      } else {
        # Priority 4: Fallback
        col_pal[[v]] <- "#D3D3D3" # LightGray
      }
    }
  }
  annotation_colors[[col_name]] <- col_pal
}

# C. Output PDF
pdf_heat_path <- file.path(out_dir, "Heatmap_Stratification_Clustered.pdf")
pdf(pdf_heat_path, width = 10, height = 8)

tryCatch({
  if(any(is.na(meta_heatmap))) {
    warning("[Viz] Metadata contains NAs. Heatmap may show white spaces in annotation.")
  }
  
  ht_obj <- plot_stratification_heatmap(
    mat_z = mat_heatmap_z,
    metadata = meta_heatmap,
    annotation_colors_list = annotation_colors,
    title = "Global Clustering (Hybrid Z-Score)"
  )
  draw(ht_obj, merge_legend = TRUE)
  
}, error = function(e) {
  warning(paste("Heatmap generation failed:", e$message))
  plot.new()
  text(0.5, 0.5, paste("Heatmap Failed:", e$message), cex = 0.8)
})

dev.off()
message(sprintf("   -> Heatmap saved to: %s", pdf_heat_path))

# 5. Statistical Testing & Advanced Interpretation
# ------------------------------------------------------------------------------
wb <- createWorkbook()

# A. Global Analysis (PERMANOVA + PLS-DA)
message("\n[Stats] Running Global PERMANOVA & Driver Analysis...")

# 1. PERMANOVA
df_stats_global <- df_global[, c("Group", safe_markers)]
perm_global <- test_coda_permanova(
  data_input = df_stats_global, 
  group_col = "Group", 
  n_perm = config$stats$n_perm
)
addWorksheet(wb, "Global_PERMANOVA")
writeData(wb, "Global_PERMANOVA", as.data.frame(perm_global), rowNames = TRUE)

# 2. PLS-DA Driver Analysis
if (config$multivariate$run_plsda && exists("plsda")) {
  message("   [PLS-DA] Identifying Global Immunological Signature...")
  
  tryCatch({
    X_pls <- df_global[, safe_markers]
    meta_pls <- df_global[, meta_cols, drop=FALSE]
    
    # Run Model
    pls_res <- run_plsda_model(X_pls, meta_pls, group_col = "Group", n_comp = config$multivariate$n_comp)
    
    # Extract Drivers
    top_drivers <- extract_plsda_loadings(pls_res, top_n = config$multivariate$top_n_loadings)
    
    # Extract Performance Metrics (Quality Check)
    perf_metrics <- extract_plsda_performance(pls_res)
    
    # Save Drivers
    addWorksheet(wb, "Global_PLSDA_Drivers")
    writeData(wb, "Global_PLSDA_Drivers", top_drivers)
    
    # Save Performance
    addWorksheet(wb, "Global_PLSDA_Quality")
    writeData(wb, "Global_PLSDA_Quality", perf_metrics)
    
    message(sprintf("   -> Identified %d top driver markers.", nrow(top_drivers)))
    message(sprintf("   -> Model Error Rate (BER): %.2f%% (Lower is better)", perf_metrics$Overall_BER[1] * 100))
    
    # Save Plot
    pdf(file.path(out_dir, "Global_PLSDA_Biplot.pdf"), width = 8, height = 7)
    mixOmics::plotIndiv(pls_res$model, comp = c(1,2), group = meta_pls$Group, 
                        ellipse = TRUE, legend = TRUE, title = "PLS-DA: Global Signature")
    dev.off()
    
  }, error = function(e) {
    warning(paste("PLS-DA Failed:", e$message))
  })
}

# B. Local PERMANOVA + ILR Decoding
if (length(ilr_list) > 0) {
  message("\n[Stats] Running Local Analysis on Compositional Groups...")
  
  summary_local <- data.frame()
  
  for (grp_name in names(ilr_list)) {
    
    mat_ilr <- ilr_list[[grp_name]]
    
    # 1. PERMANOVA
    df_stats_local <- cbind(metadata[, "Group", drop = FALSE], as.data.frame(mat_ilr))
    res_perm <- test_coda_permanova(df_stats_local, group_col = "Group", n_perm = config$stats$n_perm)
    
    pval <- res_perm$`Pr(>F)`[1]
    r2   <- res_perm$R2[1]
    
    message(sprintf("   -> Group '%s': p = %.5f (R2=%.1f%%)", grp_name, pval, r2 * 100))
    
    summary_local <- rbind(summary_local, data.frame(
      SubGroup = grp_name,
      P_Value = pval,
      R2_Percent = r2 * 100,
      F_Model = res_perm$F[1]
    ))
    
    sheet_name <- substr(paste0("ILR_", grp_name), 1, 31) 
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, as.data.frame(res_perm), rowNames = TRUE)
    
    # 2. ILR DECODING 
    if (pval < 0.1) {
      target_mks <- config$hybrid_groups[[grp_name]]
      valid_mks <- intersect(target_mks, colnames(raw_matrix))
      
      if(length(valid_mks) > 1) {
        message(sprintf("      [Decoder] Decoding balances for '%s'...", grp_name))
        
        # A. Subset Raw Matrix for Markers
        raw_parts_subset <- raw_matrix[, valid_mks, drop = FALSE]
        
        # B. Find Common Patients (INTERSECTION) 
        ids_ilr <- rownames(mat_ilr)
        ids_raw <- rownames(raw_parts_subset)
        common_ids <- intersect(ids_ilr, ids_raw)
        
        # Check if we have enough matching patients
        if (length(common_ids) < 3) {
          warning(sprintf("      [Skip] Not enough matching patients for %s decoder.", grp_name))
          next
        }
        
        # C. Align Both Matrices perfectly
        # Subset using the safe intersection list
        mat_ilr_safe <- mat_ilr[common_ids, , drop = FALSE]
        raw_parts_safe <- raw_parts_subset[common_ids, , drop = FALSE]
        
        # D. Run Decoder
        decoding_table <- decode_ilr_to_clr(mat_ilr_safe, raw_parts_safe)
        
        if (nrow(decoding_table) > 0) {
          decode_sheet <- substr(paste0("Dec_", grp_name), 1, 31)
          addWorksheet(wb, decode_sheet)
          writeData(wb, decode_sheet, decoding_table)
        }
      } else {
        message(sprintf("      [Skip] Markers for '%s' not present in raw matrix.", grp_name))
      }
    }
  }
  
  addWorksheet(wb, "Summary_Local_Tests")
  writeData(wb, "Summary_Local_Tests", summary_local)
}

# 6. Save Output
# ------------------------------------------------------------------------------

excel_path <- file.path(out_dir, "Statistical_Test_Results.xlsx")
saveWorkbook(wb, excel_path, overwrite = TRUE)

message(sprintf("\n[Output] Results saved to: %s", out_dir))
message("=== STEP 2 COMPLETE ===\n")