# src/03_statistical_analysis.R
# ==============================================================================
# STEP 03: STATISTICAL ANALYSIS
# Description: Hypothesis testing (PERMANOVA, sPLS-DA) and ILR Decoding.
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

# 1. Load Configuration & Data
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
meta_viz <- DATA$metadata 

# 2. Dynamic Selection Logic
# ------------------------------------------------------------------------------
strat_col <- config$stratification$column
message(sprintf("[Setup] Stratification Column: '%s'", strat_col))

if (!strat_col %in% names(df_global)) {
  stop(sprintf("[Error] Column '%s' not found in processed data.", strat_col))
}

# Overwrite 'Group' with the specific granularity requested
df_global$Group <- df_global[[strat_col]]
meta_viz$Group  <- meta_viz[[strat_col]]

# 3. Filter Targets & Dynamic Color Mapping
# ------------------------------------------------------------------------------
# Ensure targets are vectors (handles single or multiple values)
target_control <- as.vector(config$control_group)
target_cases <- as.vector(config$case_groups)
target_all <- unique(c(target_control, target_cases))

message(sprintf("[Setup] Target Groups: [%s] (Control) vs [%s] (Case)", 
                paste(target_control, collapse=", "), 
                paste(target_cases, collapse=", ")))

# Check availability
available_groups <- unique(df_global$Group)
missing <- setdiff(target_all, available_groups)
if(length(missing) > 0) {
  stop(sprintf("[Config Error] Requested groups not found in data: %s. Available: %s", 
               paste(missing, collapse=", "), paste(available_groups, collapse=", ")))
}

# Apply Filter
df_global <- df_global %>% filter(Group %in% target_all)
meta_viz <- meta_viz %>% filter(Patient_ID %in% df_global$Patient_ID)

# Ensure Factor Order (Control groups first)
df_global$Group <- factor(df_global$Group, levels = target_all)
meta_viz$Group <- factor(meta_viz$Group, levels = target_all)

# --- DYNAMIC COLOR ASSIGNMENT ---
clean_targets <- unique(as.character(na.omit(target_all)))
full_palette <- get_palette(config, match_groups = clean_targets)
current_palette <- full_palette[clean_targets]

if (any(is.na(current_palette))) {
  na_grps <- target_all[is.na(current_palette)]
  warning(paste("[Stats] Colors missing for groups:", paste(na_grps, collapse=", ")))
  current_palette[is.na(current_palette)] <- "grey50"
}

# Initialize Excel Workbook
out_dir <- file.path(config$output_root, "03_statistics")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
wb <- createWorkbook()

# 4a. Global PERMANOVA
# ------------------------------------------------------------------------------
message("[Stats] Running Global PERMANOVA (Multiclass)...")

df_stats_global <- df_global[, c("Group", safe_markers)]
perm_global <- test_coda_permanova(
  data_input = df_stats_global, 
  group_col = "Group", 
  n_perm = config$stats$n_perm
)

addWorksheet(wb, "Global_PERMANOVA")
writeData(wb, "Global_PERMANOVA", as.data.frame(perm_global), rowNames = TRUE)

# 4b. Pairwise PERMANOVA (Post-Hoc)
# ------------------------------------------------------------------------------
# Only run if we have more than 2 groups
if (length(unique(df_stats_global$Group)) > 2) {
  message("[Stats] Detecting >2 groups. Running Pairwise Post-Hoc Tests...")
  
  # We set min_n=4 to allow HNSCC_LS (n=4) strictly if needed, 
  # but n=5 is recommended for statistical stability. 
  # Given user context: HNSCC_LS is n=4. 
  # If we want to see it, we must lower min_n to 4, accepting low power.
  
  pair_res <- run_pairwise_permanova(
    data_input = df_stats_global, 
    group_col = "Group", 
    n_perm = config$stats$n_perm,
    min_n = config$stats$min_sample_size
  )
  
  if (!is.null(pair_res) && nrow(pair_res) > 0) {
    addWorksheet(wb, "Pairwise_PERMANOVA")
    writeData(wb, "Pairwise_PERMANOVA", pair_res)
    message(sprintf("   -> Pairwise results computed for %d pairs.", nrow(pair_res)))
  } else {
    message("   -> No pairs met the sample size criteria (min_n) for Pairwise PERMANOVA.")
  }
}

# 4c. Dispersion Check (With Warning)
# ------------------------------------------------------------------------------
message("[Stats] Running Global Dispersion Check...")
tryCatch({
  # --- Global Dispersion Test ---
  disp_global <- test_coda_dispersion(
    data_input = df_stats_global, 
    group_col = "Group", 
    n_perm = config$stats$n_perm
  )
  
  # Extract p-value for Heteroscedasticity check
  p_val_disp <- disp_global$anova_table["Groups", "Pr(>F)"]
  if (!is.na(p_val_disp) && p_val_disp < 0.05) {
    warning(sprintf("\n[STAT WARNING] Significant Beta-Dispersion detected (p=%.4f)!\n   -> PERMANOVA results may be driven by variance heterogeneity rather than location differences.\n   -> Interpret 'Global_PERMANOVA' with caution.", p_val_disp))
  }
  
  addWorksheet(wb, "Global_Dispersion_Test")
  writeData(wb, "Global_Dispersion_Test", disp_global$anova_table, rowNames = TRUE)
  addWorksheet(wb, "Global_Dispersion_Distances")
  writeData(wb, "Global_Dispersion_Distances", disp_global$group_distances)
  
  # --- Pairwise Dispersion Test ---
  if (length(unique(df_stats_global$Group)) > 2) {
    pair_disp_res <- run_pairwise_betadisper(
      data_input = df_stats_global,
      group_col = "Group",
      n_perm = config$stats$n_perm,
      min_n = config$stats$min_sample_size
    )
    
    if (!is.null(pair_disp_res) && nrow(pair_disp_res) > 0) {
      addWorksheet(wb, "Pairwise_Dispersion")
      writeData(wb, "Pairwise_Dispersion", pair_disp_res)
      message(sprintf("   -> Pairwise Dispersion results computed for %d pairs.", nrow(pair_disp_res)))
    }
  }
  
}, error = function(e) message(paste("   [ERROR] Beta-Dispersion failed:", e$message)))

# 5a. sPLS-DA (Multiclass / Stratified)
# ------------------------------------------------------------------------------
if (config$multivariate$run_plsda) {
  message(sprintf("   [sPLS-DA] Fitting STRATIFIED model (Multiclass)..."))
  set.seed(config$stats$seed) 
  
  tryCatch({
    X_pls <- df_global[, safe_markers]
    meta_stats <- meta_viz 
    
    n_rep <- if(!is.null(config$multivariate$n_repeat_cv)) config$multivariate$n_repeat_cv else 50
    
    pls_res <- run_splsda_model(X_pls, meta_stats, group_col = "Group", 
                                n_comp = config$multivariate$n_comp,
                                folds = config$multivariate$validation_folds,
                                n_repeat = n_rep) 
    
    top_drivers <- extract_plsda_loadings(pls_res)
    perf_metrics <- extract_plsda_performance(pls_res)
    
    addWorksheet(wb, "Global_sPLSDA_Drivers")
    writeData(wb, "Global_sPLSDA_Drivers", top_drivers)
    addWorksheet(wb, "Global_sPLSDA_Quality")
    writeData(wb, "Global_sPLSDA_Quality", perf_metrics)
    
    viz_report_plsda(
      pls_res = pls_res, 
      drivers_df = top_drivers, 
      metadata_viz = meta_viz, 
      colors_viz = current_palette,  
      out_path = file.path(out_dir, "Global_sPLSDA_Results.pdf"),
      group_col = "Group" 
    )
    
  }, error = function(e) message(paste("   [ERROR] Stratified sPLS-DA Failed:", e$message)))
}

# 5b. sPLS-DA (Binary: Control vs Case)
# ------------------------------------------------------------------------------
if (config$multivariate$run_plsda) {
  message(sprintf("   [sPLS-DA] Fitting BINARY model (Pooled Control vs Pooled Cases)..."))
  set.seed(config$stats$seed) 
  
  tryCatch({
    # Create Binary Labels
    meta_binary <- meta_viz %>% 
      mutate(Condition = ifelse(Group %in% config$control_group, "Control", "Case"))
    
    col_bin_case <- if(!is.null(config$colors$Case)) config$colors$Case else "firebrick"
    binary_palette <- c("Control" = config$colors$control, "Case" = col_bin_case)
    
    X_pls <- df_global[, safe_markers]
    
    n_rep <- if(!is.null(config$multivariate$n_repeat_cv)) config$multivariate$n_repeat_cv else 50
    
    pls_res_bin <- run_splsda_model(X_pls, meta_binary, group_col = "Condition", 
                                    n_comp = config$multivariate$n_comp,
                                    folds = config$multivariate$validation_folds,
                                    n_repeat = n_rep) 
    
    top_drivers_bin <- extract_plsda_loadings(pls_res_bin)
    perf_metrics_bin <- extract_plsda_performance(pls_res_bin)
    
    addWorksheet(wb, "Binary_sPLSDA_Drivers")
    writeData(wb, "Binary_sPLSDA_Drivers", top_drivers_bin)
    addWorksheet(wb, "Binary_sPLSDA_Quality")
    writeData(wb, "Binary_sPLSDA_Quality", perf_metrics_bin)
    
    viz_report_plsda(
      pls_res = pls_res_bin, 
      drivers_df = top_drivers_bin, 
      metadata_viz = meta_binary, 
      colors_viz = binary_palette,  
      out_path = file.path(out_dir, "Binary_sPLSDA_Results.pdf"),
      group_col = "Condition" 
    )
    
  }, error = function(e) message(paste("   [ERROR] Binary sPLS-DA Failed:", e$message)))
}

# 6. Local Analysis (ILR Test & Decoding)
# ------------------------------------------------------------------------------
if (length(ilr_list) > 0) {
  message("\n[Stats] Running Local Analysis on Compositional Groups (Test + Decoding)...")
  summary_local <- data.frame()
  
  for (grp_name in names(ilr_list)) {
    mat_ilr <- ilr_list[[grp_name]]
    
    # Identify samples present in current analysis
    common_ids <- intersect(df_global$Patient_ID, rownames(mat_ilr))
    
    if (length(common_ids) < 3) {
      message(sprintf("   [Skip] Group %s: Insufficient matching samples (%d)", grp_name, length(common_ids)))
      next 
    }
    
    # Subset ILR matrix
    mat_ilr_sub <- mat_ilr[common_ids, , drop=FALSE]
    match_idx <- match(common_ids, df_global$Patient_ID)
    group_sub <- df_global$Group[match_idx]
    
    # A) Run PERMANOVA on Balances
    df_stats_local <- cbind(data.frame(Group = group_sub), as.data.frame(mat_ilr_sub))
    
    res_perm <- tryCatch({
      test_coda_permanova(df_stats_local, group_col = "Group", n_perm = config$stats$n_perm)
    }, error = function(e) {
      message(sprintf("      [Error] PERMANOVA failed for %s: %s", grp_name, e$message))
      return(NULL)
    })
    
    if (!is.null(res_perm)) {
      summary_local <- rbind(summary_local, data.frame(
        SubGroup = grp_name,
        P_Value = res_perm$`Pr(>F)`[1],
        R2_Percent = res_perm$R2[1] * 100,
        F_Model = res_perm$F[1]
      ))
      
      # Save Test Results
      sheet_name_test <- substr(paste0("ILR_Test_", grp_name), 1, 31)
      if(!sheet_name_test %in% names(wb)) addWorksheet(wb, sheet_name_test)
      writeData(wb, sheet_name_test, as.data.frame(res_perm), rowNames = TRUE)
      
      # B) Run Decoding (Interpretation)
      # Retrieve the raw parts associated with this group from config to ensure correct decoding
      target_parts <- config$hybrid_groups[[grp_name]]
      
      if (!is.null(target_parts)) {
        # Subset raw matrix (Must match rows of ILR matrix exactly)
        raw_sub <- raw_matrix[rownames(mat_ilr_sub), target_parts, drop = FALSE]
        
        decoding_results <- decode_ilr_to_clr(mat_ilr_sub, raw_sub, p_threshold = 0.05)
        
        if (nrow(decoding_results) > 0) {
          sheet_name_interp <- substr(paste0("Interp_", grp_name), 1, 31)
          addWorksheet(wb, sheet_name_interp)
          writeData(wb, sheet_name_interp, decoding_results)
          message(sprintf("      -> Decoded %s: %d correlations found.", grp_name, nrow(decoding_results)))
        } else {
          message(sprintf("      -> Decoded %s: No significant correlations found.", grp_name))
        }
      }
    }
  }
  
  if(nrow(summary_local) > 0) {
    if(!"Summary_Local_Tests" %in% names(wb)) addWorksheet(wb, "Summary_Local_Tests")
    writeData(wb, "Summary_Local_Tests", summary_local)
  }
}

saveWorkbook(wb, file.path(out_dir, "Statistical_Test_Results.xlsx"), overwrite = TRUE)
message("=== STEP 3 COMPLETE ===\n")