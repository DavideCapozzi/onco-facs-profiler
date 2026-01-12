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
  library(circlize)
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

# --- DATA SETUP: Separation between Visualization (Viz) and Statistics (Stats) ---

# A. Data for Visualization (Keeps original groups: NSCLC, HNSCC, Healthy)
meta_viz <- DATA$metadata
safe_markers <- DATA$hybrid_markers
raw_matrix   <- DATA$raw_matrix
ilr_list     <- DATA$ilr_balances 

# Palette "Distinct" for plots
colors_viz <- get_palette(config)

# B. Data for Statistics (Applies Binary/Merge logic)
df_global <- DATA$hybrid_data_z # Starts with full data

mode <- config$analysis_design$mode
ls_pattern <- config$analysis_design$ls_pattern
subgroup_col <- config$metadata$subgroup_col

message(sprintf("[Setup] Applying Analysis Mode: '%s' for STATISTICS only.", mode))

# Merging Logic for df_global (Used for tests only)
if (mode == "binary") {
  df_global <- df_global %>%
    mutate(Group = case_when(
      Group %in% config$control_group ~ "Control",
      Group %in% config$case_groups ~ "Case",
      TRUE ~ Group
    ))
} else if (mode == "stratified_ls") {
  is_ls_row <- grepl(ls_pattern, df_global[[subgroup_col]])
  df_global <- df_global %>%
    mutate(Group = case_when(
      Group %in% config$control_group ~ "Control",
      Group %in% config$case_groups & is_ls_row ~ "Case_LS",
      Group %in% config$case_groups & !is_ls_row ~ "Case_Std",
      TRUE ~ Group
    ))
}

# Verify Groups
groups_viz  <- unique(meta_viz$Group)
groups_stat <- unique(df_global$Group)

message(sprintf("[Data] Visualization Groups: %s", paste(groups_viz, collapse=", ")))
message(sprintf("[Data] Statistical Groups:   %s", paste(groups_stat, collapse=", ")))

out_dir <- file.path(config$output_root, "02_exploratory")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# 2. Data diagnostics (Raw Data)
# ------------------------------------------------------------------------------
message("[Diagnostics] Generating unified distribution plots for RAW data...")

# Use meta_viz (Original Groups) for distribution plots
meta_ordered <- meta_viz[match(rownames(raw_matrix), meta_viz$Patient_ID), ]
df_raw_viz <- cbind(meta_ordered, as.data.frame(raw_matrix))

raw_pdf_path <- file.path(out_dir, "Distributions_Raw_Distinct.pdf")
hl_pattern <- if(!is.null(config$viz$highlight_pattern)) config$viz$highlight_pattern else ""

pdf(raw_pdf_path, width = 8, height = 6) 
for (marker in DATA$markers) {
  if (marker %in% names(df_raw_viz)) {
    # Use colors_viz to distinguish distinct tumors
    p <- plot_raw_distribution_merged(
      data_df = df_raw_viz, 
      marker_name = marker,
      colors = colors_viz,     
      highlight_pattern = hl_pattern,   
      highlight_color = "#FFD700"  
    )
    if (!is.null(p)) print(p)
  }
}
dev.off()


# 3. PCA (Exploratory) - Multi-Page
# ------------------------------------------------------------------------------
message("[PCA] Running PCA on Global Z-Scored Matrix...")

# PCA calculated on numeric matrix (independent of groups)
pca_input <- df_global[, safe_markers]
res_pca <- PCA(pca_input, scale.unit = FALSE, graph = FALSE)

# Setup Highlight (Optional)
highlight_map <- c()
if (!is.null(config$viz$highlight_pattern)) {
  highlight_map[[config$viz$highlight_pattern]] <- "#FFD700"
}

# A. Individuals PCA plots (Multi-Page)
message("[PCA] Saving Individuals plots to single PDF...")
pdf_ind_path <- file.path(out_dir, "PCA_Global_Individuals_MultiPage.pdf")
pdf(pdf_ind_path, width = 9, height = 7)

# Page 1: PC1 vs PC2
print(plot_pca_custom(res_pca, meta_viz, colors_viz, dims = c(1, 2), 
                      show_labels = TRUE, highlight_patterns = highlight_map))
# Page 2: PC1 vs PC3
print(plot_pca_custom(res_pca, meta_viz, colors_viz, dims = c(1, 3), 
                      show_labels = TRUE, highlight_patterns = highlight_map))
# Page 3: PC2 vs PC3
print(plot_pca_custom(res_pca, meta_viz, colors_viz, dims = c(2, 3), 
                      show_labels = TRUE, highlight_patterns = highlight_map))

dev.off() 

# B. Variables/Loadings plots (Multi-Page)
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
    theme_bw()
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


# 4. Patient Stratification Heatmap (Restored Logic)
# ------------------------------------------------------------------------------
message("[Viz] Generating Patient Stratification Heatmap...")

mat_heatmap_z <- as.matrix(df_global[, safe_markers])
rownames(mat_heatmap_z) <- df_global$Patient_ID

# A. Dynamic Metadata Selection (from Config)
target_cols <- config$viz$heatmap_metadata
if (is.null(target_cols)) target_cols <- c("Group")
valid_cols <- intersect(target_cols, names(meta_viz))
if (length(valid_cols) == 0) valid_cols <- c("Group")

meta_heatmap <- meta_viz[, valid_cols, drop = FALSE]

# B. Smart Color Generation (Restored from oldold)
annotation_colors <- list()

# Group Colors (Base Palette)
present_groups <- unique(meta_heatmap$Group)
group_pal <- colors_viz[intersect(names(colors_viz), present_groups)]
annotation_colors[["Group"]] <- group_pal

# Extra Metadata Colors (Inheritance or Highlight Logic)
extra_cols <- setdiff(valid_cols, "Group")
hl_pattern <- config$viz$highlight_pattern
hl_color   <- if(!is.null(config$viz$highlight_color)) config$viz$highlight_color else "#FFD700"

for (col_name in extra_cols) {
  raw_vals <- as.character(meta_heatmap[[col_name]])
  vals <- sort(unique(raw_vals[!is.na(raw_vals)]))
  col_pal <- character()
  
  for (v in vals) {
    if (v %in% names(group_pal)) {
      col_pal[[v]] <- group_pal[[v]]
    } else if (!is.null(hl_pattern) && grepl(hl_pattern, v)) {
      col_pal[[v]] <- hl_color
    } else {
      # Try to find parent match
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
        col_pal[[v]] <- "#D3D3D3" # LightGray fallback
      }
    }
  }
  annotation_colors[[col_name]] <- col_pal
}

# C. Output PDF
pdf_heat_path <- file.path(out_dir, "Heatmap_Stratification.pdf")
pdf(pdf_heat_path, width = 10, height = 8)

tryCatch({
  ht_obj <- plot_stratification_heatmap(
    mat_z = mat_heatmap_z,
    metadata = meta_heatmap,
    annotation_colors_list = annotation_colors,
    title = paste0("Clustering (Mode: ", mode, ")")
  )
  draw(ht_obj, merge_legend = TRUE)
}, error = function(e) warning(paste("Heatmap failed:", e$message)))

dev.off()


# 5. Statistical Testing & sPLS-DA (Analysis Groups)
# ------------------------------------------------------------------------------
wb <- createWorkbook()

# A. Global Analysis
message("\n[Stats] Running Global PERMANOVA & Driver Analysis...")

# 1. PERMANOVA (Uses df_global -> Merged Groups)
df_stats_global <- df_global[, c("Group", safe_markers)]
perm_global <- test_coda_permanova(
  data_input = df_stats_global, 
  group_col = "Group", 
  n_perm = config$stats$n_perm
)
addWorksheet(wb, "Global_PERMANOVA")
writeData(wb, "Global_PERMANOVA", as.data.frame(perm_global), rowNames = TRUE)

# 2. sPLS-DA (Uses df_global -> Merged Groups)
run_pls <- if(!is.null(config$multivariate$run_plsda)) config$multivariate$run_plsda else FALSE

if (run_pls) {
  message(sprintf("   [sPLS-DA] Fitting model on Statistical Groups: %s", paste(unique(df_global$Group), collapse=" vs ")))
  
  tryCatch({
    X_pls <- df_global[, safe_markers]
    # Metadata for training (contains Case/Control)
    meta_stats <- df_global[, colnames(DATA$metadata), drop=FALSE] 
    
    # Run Sparse Model
    pls_res <- run_splsda_model(X_pls, meta_stats, group_col = "Group", 
                                n_comp = config$multivariate$n_comp,
                                folds = config$multivariate$validation_folds)
    
    # Extract Results
    top_drivers <- extract_plsda_loadings(pls_res)
    perf_metrics <- extract_plsda_performance(pls_res)
    tuning_info  <- as.data.frame(pls_res$tuning$choice.keepX)
    colnames(tuning_info) <- "Selected_Features"
    
    # Save Excel
    addWorksheet(wb, "Global_sPLSDA_Drivers")
    writeData(wb, "Global_sPLSDA_Drivers", top_drivers)
    addWorksheet(wb, "Global_sPLSDA_Quality")
    writeData(wb, "Global_sPLSDA_Quality", perf_metrics)
    writeData(wb, "Global_sPLSDA_Quality", tuning_info, startRow = 6, startCol = 1)
    
    # --- PLOTS for sPLS-DA ---
    pdf(file.path(out_dir, "Global_sPLSDA_Results.pdf"), width = 11, height = 8)
    
    # Plot 1: Biplot (Visualized with distinct groups from meta_viz)
    plot_group_factor <- as.factor(meta_viz$Group)
    plot_colors <- colors_viz[levels(plot_group_factor)]
    
    # Assign to variable 'void' to prevent verbose console output
    void <- mixOmics::plotIndiv(pls_res$model, 
                                comp = c(1,2), 
                                group = plot_group_factor, # Color by NSCLC/HNSCC/Healthy
                                col.per.group = plot_colors, 
                                ellipse = TRUE, 
                                legend = TRUE, 
                                title = "sPLS-DA: Binary Model (Colored by Subtype)",
                                star.plot = TRUE)
    
    # Plot 2: Loadings (Component 1)
    p_load <- ggplot(top_drivers, aes(x = reorder(Marker, Comp1_Weight), y = Comp1_Weight, fill = Direction)) +
      geom_bar(stat = "identity", width = 0.7) +
      coord_flip() +
      scale_fill_manual(values = c("Positive_Assoc" = "#CD5C5C", "Negative_Assoc" = "#4682B4")) +
      labs(title = "sPLS-DA Selected Features (Comp 1)", 
           subtitle = "Drivers of the Statistical Separation",
           x = "Marker", y = "Weight") +
      theme_minimal()
    print(p_load)
    
    # Plot 3: Clustered Image Map (CIM)
    if (nrow(top_drivers) > 1) {
      row_cols <- colors_viz[as.character(meta_viz$Group)]
      mixOmics::cim(pls_res$model, 
                    row.sideColors = row_cols,
                    title = "Signature Heatmap (CIM)",
                    save = NULL)
    }
    
    # Plot 4: Stability
    if (!is.null(pls_res$performance$features$stable[[1]])) {
      stab_vec <- pls_res$performance$features$stable[[1]]
      df_stab <- data.frame(Marker = names(stab_vec), Stability = as.numeric(stab_vec)) %>%
        arrange(desc(Stability)) %>% head(20)
      
      p_stab <- ggplot(df_stab, aes(x = reorder(Marker, Stability), y = Stability)) +
        geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
        coord_flip() +
        labs(title = "Feature Stability (Comp 1)", y = "Frequency of Selection", x = "") +
        theme_minimal()
      print(p_stab)
    }
    
    dev.off()
    
  }, error = function(e) {
    message(paste("   [ERROR] sPLS-DA Execution Failed:", e$message))
  })
}


# B. Local PERMANOVA + ILR Decoding (Restored)
# ------------------------------------------------------------------------------
if (length(ilr_list) > 0) {
  message("\n[Stats] Running Local Analysis on Compositional Groups...")
  
  summary_local <- data.frame()
  
  # Ensure we use the statistical group (df_global$Group) for testing
  # This propagates the "Binary" or "Stratified" logic to ILR tests
  stat_groups <- df_global$Group 
  
  for (grp_name in names(ilr_list)) {
    
    mat_ilr <- ilr_list[[grp_name]]
    
    # 1. PERMANOVA on ILR balances
    # We construct a temp dataframe with the Binary/Stat group
    df_stats_local <- cbind(data.frame(Group = stat_groups), as.data.frame(mat_ilr))
    
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
    
    # Save detailed Permanova table
    sheet_name <- substr(paste0("ILR_", grp_name), 1, 31) 
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, as.data.frame(res_perm), rowNames = TRUE)
    
    # 2. ILR DECODING (Correlating abstract balances with raw markers)
    # We only decode if there is a signal (p < 0.1 relaxed threshold)
    if (pval < 0.1) {
      target_mks <- config$hybrid_groups[[grp_name]]
      valid_mks <- intersect(target_mks, colnames(raw_matrix))
      
      if(length(valid_mks) > 1) {
        # Subset Raw Matrix for Markers & Common Patients
        raw_parts_subset <- raw_matrix[, valid_mks, drop = FALSE]
        ids_ilr <- rownames(mat_ilr)
        ids_raw <- rownames(raw_parts_subset)
        common_ids <- intersect(ids_ilr, ids_raw)
        
        if (length(common_ids) >= 3) {
          mat_ilr_safe <- mat_ilr[common_ids, , drop = FALSE]
          raw_parts_safe <- raw_parts_subset[common_ids, , drop = FALSE]
          
          decoding_table <- decode_ilr_to_clr(mat_ilr_safe, raw_parts_safe)
          
          if (nrow(decoding_table) > 0) {
            decode_sheet <- substr(paste0("Dec_", grp_name), 1, 31)
            addWorksheet(wb, decode_sheet)
            writeData(wb, decode_sheet, decoding_table)
          }
        }
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