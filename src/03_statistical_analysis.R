# src/03_statistical_analysis.R
# ==============================================================================
# STEP 03: STATISTICAL ANALYSIS
# Description: Hypothesis testing (PERMANOVA, sPLS-DA)
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

# --- DYNAMIC COLOR ASSIGNMENT (Vector Safe) ---
assign_color <- function(grp_name, cfg) {
  # 1. Check specific definition (Safely check names first)
  if (grp_name %in% names(cfg$colors$groups)) {
    return(cfg$colors$groups[[grp_name]])
  }
  # 2. Check if valid control (using %in% handles vectors correctly)
  if (grp_name %in% cfg$control_group) {
    return(cfg$colors$control)
  }
  # 3. Fallback
  return(NULL)
}

# Assign colors to Control groups
cols_control <- setNames(
  sapply(target_control, function(g) {
    col <- assign_color(g, config)
    if(is.null(col)) return("grey50") else return(col)
  }), 
  target_control
)

# Assign colors to Case groups
cols_cases <- c()
generic_palette <- config$colors$cases
generic_idx <- 1

for (case in target_cases) {
  specific_col <- assign_color(case, config)
  if (!is.null(specific_col)) {
    cols_cases[case] <- specific_col
  } else {
    cols_cases[case] <- generic_palette[(generic_idx - 1) %% length(generic_palette) + 1]
    generic_idx <- generic_idx + 1
  }
}

current_palette <- c(cols_control, cols_cases)

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

# 4b. Global Dispersion Check
# ------------------------------------------------------------------------------
message("[Stats] Running Global Dispersion Check...")
tryCatch({
  disp_global <- test_coda_dispersion(
    data_input = df_stats_global, 
    group_col = "Group", 
    n_perm = config$stats$n_perm
  )
  addWorksheet(wb, "Global_Dispersion_Test")
  writeData(wb, "Global_Dispersion_Test", disp_global$anova_table, rowNames = TRUE)
  addWorksheet(wb, "Global_Dispersion_Distances")
  writeData(wb, "Global_Dispersion_Distances", disp_global$group_distances)
}, error = function(e) message(paste("   [ERROR] Beta-Dispersion failed:", e$message)))

# 5. sPLS-DA (Multiclass / Stratified)
# ------------------------------------------------------------------------------
if (config$multivariate$run_plsda) {
  message(sprintf("   [sPLS-DA] Fitting STRATIFIED model (Multiclass)..."))
  set.seed(config$stats$seed) 
  
  tryCatch({
    X_pls <- df_global[, safe_markers]
    meta_stats <- meta_viz 
    
    pls_res <- run_splsda_model(X_pls, meta_stats, group_col = "Group", 
                                n_comp = config$multivariate$n_comp,
                                folds = config$multivariate$validation_folds)
    
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

# 6. sPLS-DA (Binary: Control vs Case)
# ------------------------------------------------------------------------------
if (config$multivariate$run_plsda) {
  message(sprintf("   [sPLS-DA] Fitting BINARY model (Pooled Control vs Pooled Cases)..."))
  set.seed(config$stats$seed) 
  
  tryCatch({
    # Create Binary Labels
    # Logic: If Group is ANY of the control groups -> Control, else Case
    meta_binary <- meta_viz %>% 
      mutate(Condition = ifelse(Group %in% config$control_group, "Control", "Case"))
    
    col_bin_case <- if(!is.null(config$colors$Case)) config$colors$Case else "firebrick"
    binary_palette <- c("Control" = config$colors$control, "Case" = col_bin_case)
    
    X_pls <- df_global[, safe_markers]
    
    pls_res_bin <- run_splsda_model(X_pls, meta_binary, group_col = "Condition", 
                                    n_comp = config$multivariate$n_comp,
                                    folds = config$multivariate$validation_folds)
    
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

# 7. Local Analysis (ILR)
# ------------------------------------------------------------------------------
if (length(ilr_list) > 0) {
  message("\n[Stats] Running Local Analysis on Compositional Groups...")
  summary_local <- data.frame()
  
  for (grp_name in names(ilr_list)) {
    mat_ilr <- ilr_list[[grp_name]]
    
    common_ids <- intersect(df_global$Patient_ID, rownames(mat_ilr))
    
    if (length(common_ids) < 3) {
      message(sprintf("   [Skip] Group %s: Insufficient matching samples (%d)", grp_name, length(common_ids)))
      next 
    }
    
    mat_ilr_sub <- mat_ilr[common_ids, , drop=FALSE]
    match_idx <- match(common_ids, df_global$Patient_ID)
    group_sub <- df_global$Group[match_idx]
    
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
      
      sheet_name <- substr(paste0("ILR_", grp_name), 1, 31)
      if(!sheet_name %in% names(wb)) addWorksheet(wb, sheet_name)
      writeData(wb, sheet_name, as.data.frame(res_perm), rowNames = TRUE)
    }
  }
  
  if(nrow(summary_local) > 0) {
    if(!"Summary_Local_Tests" %in% names(wb)) addWorksheet(wb, "Summary_Local_Tests")
    writeData(wb, "Summary_Local_Tests", summary_local)
  }
}

saveWorkbook(wb, file.path(out_dir, "Statistical_Test_Results.xlsx"), overwrite = TRUE)
message("=== STEP 3 COMPLETE ===\n")