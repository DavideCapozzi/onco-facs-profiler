# src/03_statistical_analysis.R
# ==============================================================================
# STEP 03: STATISTICAL ANALYSIS
# Description: Hypothesis testing (PERMANOVA, PLS-DA, Driver Analysis).
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
  library(vegan)
  library(mixOmics)
})

source("R/utils_io.R")          
source("R/modules_hypothesis.R") 
source("R/modules_multivariate.R")   
source("R/modules_interpretation.R")
source("R/modules_viz.R") 

message("\n=== PIPELINE STEP 3: STATISTICAL ANALYSIS ===")

# 1. Load Config & Data
# ------------------------------------------------------------------------------
config <- load_config("config/global_params.yml")
input_file <- file.path(config$output_root, "01_data_processing", "data_processed.rds")

if (!file.exists(input_file)) stop("Step 01 output not found. Run src/01_data_processing.R first.")

DATA <- readRDS(input_file)

# Data Setup
df_global <- DATA$hybrid_data_z
safe_markers <- DATA$hybrid_markers
raw_matrix <- DATA$raw_matrix
ilr_list <- DATA$ilr_balances
meta_orig <- DATA$metadata

# Define Analysis Mode (Binary vs Stratified)
mode <- config$analysis_design$mode
ls_pattern <- config$analysis_design$ls_pattern
subgroup_col <- config$metadata$subgroup_col

message(sprintf("[Setup] Applying Analysis Mode: '%s'", mode))

# Apply Merging Logic to df_global
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

message(sprintf("[Data] Statistical Groups: %s", paste(unique(df_global$Group), collapse=", ")))

out_dir <- file.path(config$output_root, "03_statistics")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Initialize Excel Workbook
wb <- createWorkbook()

# 2. Global PERMANOVA
# ------------------------------------------------------------------------------
message("[Stats] Running Global PERMANOVA...")

df_stats_global <- df_global[, c("Group", safe_markers)]
perm_global <- test_coda_permanova(
  data_input = df_stats_global, 
  group_col = "Group", 
  n_perm = config$stats$n_perm
)

addWorksheet(wb, "Global_PERMANOVA")
writeData(wb, "Global_PERMANOVA", as.data.frame(perm_global), rowNames = TRUE)


# 3. sPLS-DA (Driver Analysis)
# ------------------------------------------------------------------------------
run_pls <- if(!is.null(config$multivariate$run_plsda)) config$multivariate$run_plsda else FALSE

if (run_pls) {
  message(sprintf("   [sPLS-DA] Fitting model on Statistical Groups..."))
  
  tryCatch({
    X_pls <- df_global[, safe_markers]
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
    
    # Save Excel Data
    addWorksheet(wb, "Global_sPLSDA_Drivers")
    writeData(wb, "Global_sPLSDA_Drivers", top_drivers)
    addWorksheet(wb, "Global_sPLSDA_Quality")
    writeData(wb, "Global_sPLSDA_Quality", perf_metrics)
    writeData(wb, "Global_sPLSDA_Quality", tuning_info, startRow = 6, startCol = 1)
    
    # Generate Plots using the Viz Module
    # Note: We pass original metadata/colors for Visualization (Points colored by subtype)
    colors_viz <- get_palette(config)
    
    # Identify the levels used in PLS
    pls_levels <- levels(factor(meta_stats$Group))
    # Assuming standard order (Control first, Case second) or check Config
    lbl_vec <- c("Negative" = "Reference/Control", "Positive" = "Case/Target")
    
    # Try to be more specific if possible using config
    if (length(config$case_groups) > 0) {
      lbl_vec["Positive"] <- paste(config$case_groups, collapse="+")
    }
    if (!is.null(config$control_group)) {
      lbl_vec["Negative"] <- config$control_group
    }
    
    viz_report_plsda(
      pls_res = pls_res, 
      drivers_df = top_drivers, 
      metadata_viz = meta_orig, 
      colors_viz = colors_viz, 
      out_path = file.path(out_dir, "Global_sPLSDA_Results.pdf"),
      binary_labels = lbl_vec  
    )
    
  }, error = function(e) {
    message(paste("   [ERROR] sPLS-DA Execution Failed:", e$message))
  })
}


# 4. Local Analysis (ILR & Decoding)
# ------------------------------------------------------------------------------
if (length(ilr_list) > 0) {
  message("\n[Stats] Running Local Analysis on Compositional Groups...")
  
  summary_local <- data.frame()
  stat_groups <- df_global$Group 
  
  for (grp_name in names(ilr_list)) {
    
    mat_ilr <- ilr_list[[grp_name]]
    
    # PERMANOVA on ILR balances
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
    
    sheet_name <- substr(paste0("ILR_", grp_name), 1, 31) 
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, as.data.frame(res_perm), rowNames = TRUE)
    
    # Decoding (Correlating balances with raw markers)
    if (pval < 0.1) {
      target_mks <- config$hybrid_groups[[grp_name]]
      valid_mks <- intersect(target_mks, colnames(raw_matrix))
      
      if(length(valid_mks) > 1) {
        raw_parts_subset <- raw_matrix[, valid_mks, drop = FALSE]
        common_ids <- intersect(rownames(mat_ilr), rownames(raw_parts_subset))
        
        if (length(common_ids) >= 3) {
          decoding_table <- decode_ilr_to_clr(
            mat_ilr[common_ids, , drop = FALSE], 
            raw_parts_subset[common_ids, , drop = FALSE]
          )
          
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

# 5. Save Final Excel
# ------------------------------------------------------------------------------
excel_path <- file.path(out_dir, "Statistical_Test_Results.xlsx")
saveWorkbook(wb, excel_path, overwrite = TRUE)

message(sprintf("\n[Output] Results saved to: %s", out_dir))
message("=== STEP 3 COMPLETE ===\n")