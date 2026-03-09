# src/03_statistical_analysis.R
# ==============================================================================
# STEP 03: STATISTICAL ANALYSIS & SCENARIO ORCHESTRATION
# Description: Global Hypothesis testing (PERMANOVA, sPLS-DA) and 
#              Iterative Multi-Scenario execution via workflows module.
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
source("R/modules_viz.R") 
source("R/workflows.R") # Injected Orchestrator Macros

message("\n=== PIPELINE STEP 3: STATISTICAL ANALYSIS & SCENARIOS ===")

# ==============================================================================
# 1. INITIALIZATION & DATA SETUP
# ==============================================================================
config <- load_config("config/global_params.yml")
input_file <- file.path(config$output_root, "01_data_processing", "data_processed.rds")

if (!file.exists(input_file)) stop("[Fatal] Step 01 output not found. Run src/01_data_processing.R first.")
DATA <- readRDS(input_file)

# Setup Directories
results_dir <- file.path(config$output_root, "results_analysis")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

global_dir <- file.path(results_dir, "Global_Stats")
if (!dir.exists(global_dir)) dir.create(global_dir, recursive = TRUE)

# Initialize the Universal Master Workbook
wb_master <- createWorkbook()

# Data Setup
df_global <- DATA$hybrid_data_z
safe_markers <- DATA$hybrid_markers
raw_matrix <- DATA$raw_matrix
meta_viz <- DATA$metadata 

strat_col <- config$stratification$column
message(sprintf("[Setup] Stratification Column: '%s'", strat_col))

if (!strat_col %in% names(df_global)) {
  stop(sprintf("[Fatal] Column '%s' not found in processed data.", strat_col))
}

# Overwrite 'Group' with the specific granularity requested
df_global$Group <- df_global[[strat_col]]
meta_viz$Group  <- meta_viz[[strat_col]]

# ==============================================================================
# 2. GLOBAL STATISTICAL ANALYSIS
# ==============================================================================
message("\n--- RUNNING GLOBAL STATISTICS ---")

target_control <- as.vector(config$control_group)
target_cases <- as.vector(config$case_groups)
target_all <- unique(c(target_control, target_cases))

# Filter Global Data
df_stats_global <- df_global %>% filter(Group %in% target_all)
meta_stats_global <- meta_viz %>% filter(Patient_ID %in% df_stats_global$Patient_ID)

df_stats_global$Group <- factor(df_stats_global$Group, levels = target_all)
meta_stats_global$Group <- factor(meta_stats_global$Group, levels = target_all)

# Dynamic Color Assignment
clean_targets <- unique(as.character(na.omit(target_all)))
current_palette <- get_palette(config, match_groups = clean_targets)[clean_targets]
if (any(is.na(current_palette))) current_palette[is.na(current_palette)] <- "grey50"

# 2.1 Global PERMANOVA
message("   [Stats] Running Global PERMANOVA (Multiclass)...")
perm_global <- test_coda_permanova(
  data_input = df_stats_global[, c("Group", safe_markers)], 
  group_col = "Group", 
  n_perm = config$stats$n_perm
)
addWorksheet(wb_master, "Global_PERMANOVA")
writeData(wb_master, "Global_PERMANOVA", as.data.frame(perm_global), rowNames = TRUE)

# 2.2 Global Dispersion & All-Pairs Contrasts
message("   [Stats] Running Global Dispersion Check...")
tryCatch({
  # 1. Global Dispersion Test (all groups together)
  disp_global <- test_coda_dispersion(
    data_input = df_stats_global[, c("Group", safe_markers)], 
    group_col = "Group", 
    n_perm = config$stats$n_perm
  )
  p_val_disp <- disp_global$anova_table["Groups", "Pr(>F)"]
  if (!is.na(p_val_disp) && p_val_disp < 0.05) {
    warning(sprintf("\n[STAT WARNING] Significant Beta-Dispersion detected (p=%.4f)!\n   -> Interpret 'Global_PERMANOVA' with caution.", p_val_disp))
  }
  addWorksheet(wb_master, "Global_Dispersion")
  writeData(wb_master, "Global_Dispersion", disp_global$anova_table, rowNames = TRUE)
  
  # 2. Pairwise Dispersion Test (all combinations)
  message("   [Stats] Running All-Pairs Dispersion Check...")
  disp_pairwise <- run_pairwise_betadisper(
    data_input = df_stats_global[, c("Group", safe_markers)], 
    group_col = "Group", 
    metadata_cols = c("Patient_ID"),
    n_perm = config$stats$n_perm,
    min_n = if(!is.null(config$stats$min_sample_size)) config$stats$min_sample_size else 4
  )
  
  if (!is.null(disp_pairwise) && nrow(disp_pairwise) > 0) {
    addWorksheet(wb_master, "All_Pairs_Dispersion")
    writeData(wb_master, "All_Pairs_Dispersion", disp_pairwise)
  }
  
}, error = function(e) message(paste("   [ERROR] Beta-Dispersion failed:", e$message)))

# 2.3 Global sPLS-DA (Multiclass)
if (config$multivariate$run_plsda) {
  message("   [sPLS-DA] Fitting STRATIFIED model (Multiclass)...")
  set.seed(config$stats$seed) 
  tryCatch({
    pls_res <- run_splsda_model(
      data_z = df_stats_global[, safe_markers], 
      metadata = meta_stats_global, 
      group_col = "Group", 
      n_comp = config$multivariate$n_comp,
      folds = config$multivariate$validation_folds,
      n_repeat = if(!is.null(config$multivariate$n_repeat_cv)) config$multivariate$n_repeat_cv else 50
    ) 
    top_drivers <- extract_plsda_loadings(pls_res)
    addWorksheet(wb_master, "Global_sPLSDA_Drivers")
    writeData(wb_master, "Global_sPLSDA_Drivers", top_drivers)
    
    n_boxes <- if(!is.null(config$multivariate$top_n_loadings)) config$multivariate$top_n_loadings else 9
    
    viz_report_plsda(
      pls_res = pls_res, 
      drivers_df = top_drivers, 
      metadata_viz = meta_stats_global, 
      colors_viz = current_palette,  
      out_path = file.path(global_dir, "Global_sPLSDA_Results.pdf"),
      group_col = "Group", 
      n_top_boxplots = n_boxes
    )
  }, error = function(e) message(paste("   [ERROR] Stratified sPLS-DA Failed:", e$message)))
}

# ==============================================================================
# 3. SCENARIO ORCHESTRATION (Iterative Execution)
# ==============================================================================
message("\n--- RUNNING SCENARIO-SPECIFIC STATISTICS ---")

all_results <- list()

if (is.null(config$analysis_scenarios) || length(config$analysis_scenarios) == 0) {
  warning("[Stats] No 'analysis_scenarios' found in config. Skipping Scenario Loop.")
} else {
  
  full_meta <- DATA$metadata
  if (strat_col %in% colnames(full_meta)) {
    full_meta$Group <- full_meta[[strat_col]]
  }
  
  full_mat <- DATA$hybrid_data_z %>% 
    dplyr::select(dplyr::all_of(safe_markers)) %>% 
    as.matrix()
  rownames(full_mat) <- DATA$metadata$Patient_ID
  
  for (scenario in config$analysis_scenarios) {
    message(sprintf("\n   [Scenario] Executing: %s", scenario$id))
    
    scen_dir <- file.path(results_dir, scenario$id)
    if (!dir.exists(scen_dir)) dir.create(scen_dir, recursive = TRUE)
    
    tryCatch({
      # Delegate processing to module logic ensuring pure functions and clean I/O
      scen_res <- run_scenario_pipeline(
        scenario = scenario, 
        full_mat = full_mat, 
        full_meta = full_meta, 
        config = config, 
        out_dir = scen_dir
      )
      
      if (!is.null(scen_res)) {
        # Securely write valid outputs to Master Workbook
        
        # Dispersion Output
        if (!is.null(scen_res$dispersion)) {
          disp_sheet <- substr(paste0(scenario$id, "_Disp"), 1, 31)
          addWorksheet(wb_master, disp_sheet)
          writeData(wb_master, disp_sheet, scen_res$dispersion, rowNames = TRUE)
        }
        
        # PERMANOVA Output
        if (!is.null(scen_res$permanova)) {
          perm_sheet <- substr(paste0(scenario$id, "_Perm"), 1, 31)
          addWorksheet(wb_master, perm_sheet)
          writeData(wb_master, perm_sheet, scen_res$permanova, rowNames = TRUE)
        }
        
        # Drivers (sPLS-DA) Output
        if (!is.null(scen_res$drivers) && nrow(scen_res$drivers) > 0) {
          drv_sheet <- substr(paste0(scenario$id, "_Drv"), 1, 31)
          addWorksheet(wb_master, drv_sheet)
          writeData(wb_master, drv_sheet, scen_res$drivers)
        }
        
        # Collect for Consolidated Reporting
        all_results[[scenario$id]] <- list(id = scenario$id, drivers = scen_res$drivers)
      }
      
    }, error = function(e) {
      message(sprintf("!!! [Error] Scenario '%s' failed: %s", scenario$id, e$message))
    })
  }
  
  # ==============================================================================
  # 4. CONSOLIDATED REPORTING
  # ==============================================================================
  if (length(all_results) > 0) {
    message("\n   [Output] Generating Consolidated Multi-Scenario Summary...")
    
    combined_drivers_list <- lapply(names(all_results), function(scen_id) {
      res <- all_results[[scen_id]]
      if (is.null(res$drivers) || nrow(res$drivers) == 0) return(NULL)
      
      df <- res$drivers
      df$Scenario <- scen_id 
      
      df_clean <- df %>%
        dplyr::rename_with(.fn = ~ "Weight_Case_VS_Control_PC1", .cols = dplyr::matches("^Weight_.*_PC1$")) %>%
        dplyr::rename_with(.fn = ~ "Weight_Case_VS_Control_PC2", .cols = dplyr::matches("^Weight_.*_PC2$")) %>%
        dplyr::select(Scenario, Marker, dplyr::matches("^Weight_Case_VS_Control_PC[12]$"))
      
      return(df_clean)
    })
    
    combined_drivers <- dplyr::bind_rows(combined_drivers_list)
    
    if (nrow(combined_drivers) > 0) {
      addWorksheet(wb_master, "All_Drivers_Summary")
      writeData(wb_master, "All_Drivers_Summary", combined_drivers)
      
      # Safely move the summary sheet to the first position
      n_sheets <- length(names(wb_master))
      if (n_sheets > 1) {
        openxlsx::worksheetOrder(wb_master) <- c(n_sheets, 1:(n_sheets - 1))
      }
    }
  }
}

# Save the unified Master Workbook
master_report_path <- file.path(results_dir, "Multi_Scenario_Analysis_Report.xlsx")
saveWorkbook(wb_master, master_report_path, overwrite = TRUE)
message(sprintf("   [Output] Master Report saved: %s", master_report_path))

message("=== STEP 3 COMPLETE ===\n")