# R/utils_io.R
# ==============================================================================
# INFRASTRUCTURE UTILITIES
# Description: Config loading, Logging setup, and basic I/O helpers.
# ==============================================================================

library(yaml)
library(readxl)
library(dplyr)
library(openxlsx)

#' @title Load Global Configuration
#' @description Loads the YAML config and validates critical paths.
#' @param config_path Path to the YAML file.
#' @return A list containing configuration parameters.
load_config <- function(config_path = "config/global_params.yml") {
  if (!file.exists(config_path)) {
    stop(sprintf("Configuration file not found at: %s", config_path))
  }
  
  config <- yaml::read_yaml(config_path)
  
  # If project_name is defined, append it to output_root
  if (!is.null(config$project_name) && config$project_name != "") {
    config$output_root <- paste0(config$output_root, "_", config$project_name)
  }
  
  message(sprintf("[System] Configuration loaded from %s", config_path))
  
  return(config)
}

#' @title Load Raw Data (Excel)
#' @description Iterates over sheets defined in config and merges them, preserving metadata.
#' @param config The loaded configuration object.
#' @return A raw dataframe with Patient_ID, Group, and Subgroup info.
load_raw_data <- function(config) {
  input_file <- config$input_file
  subgroup_col <- config$metadata$subgroup_col # Read from config
  
  if (!file.exists(input_file)) {
    stop(sprintf("Input data file not found: %s", input_file))
  }
  
  df_list <- list()
  message("[IO] Loading Excel sheets...")
  
  for (cohort_name in names(config$cohorts)) {
    sheet_name <- config$cohorts[[cohort_name]]
    
    tryCatch({
      # Load sheet
      raw_tmp <- read_excel(input_file, sheet = sheet_name)
      
      # Standardize Patient_ID
      colnames(raw_tmp)[1] <- "Patient_ID"
      
      # 1. Identify Marker columns (exclude ID and the specific Subgroup column)
      cols_to_exclude <- c("Patient_ID")
      if (!is.null(subgroup_col) && subgroup_col %in% colnames(raw_tmp)) {
        cols_to_exclude <- c(cols_to_exclude, subgroup_col)
      }
      
      # 2. Convert ONLY marker columns to numeric
      clean_tmp <- raw_tmp %>%
        mutate(across(-all_of(cols_to_exclude), ~suppressWarnings(as.numeric(as.character(.))))) %>%
        mutate(Group = cohort_name, .after = Patient_ID)
      
      df_list[[cohort_name]] <- clean_tmp
      message(sprintf("    -> Loaded %s: %d samples", cohort_name, nrow(clean_tmp)))
      
    }, error = function(e) {
      warning(sprintf("    [WARN] Failed to load sheet '%s': %s", sheet_name, e$message))
    })
  }
  
  full_data <- bind_rows(df_list)
  return(full_data)
}

#' @title Save Quality Control Report to Excel
#' @description Generates a multi-sheet Excel report with filtering statistics and group breakdowns.
#' @param qc_list A list containing summary stats and dataframes of dropped items.
#' @param out_path Path to save the .xlsx file.
save_qc_report <- function(qc_list, out_path) {
  
  wb <- createWorkbook()
  
  # --- Sheet 1: Summary with Group Breakdown ---
  addWorksheet(wb, "Summary")
  
  # 1. Prepare Base Metrics
  metrics <- c("Initial Samples", "Initial Markers",
               "Dropped Samples (A Priori)", "Dropped Markers (A Priori)",
               "Dropped Samples (QC)", "Dropped Markers (High NA)", "Dropped Markers (Zero Var)",
               "Final Samples", "Final Markers")
  
  # 2. Extract Total Counts
  n_samples_apriori <- if(!is.null(qc_list$dropped_samples_apriori)) nrow(qc_list$dropped_samples_apriori) else 0
  n_markers_apriori <- if(!is.null(qc_list$dropped_markers_apriori)) nrow(qc_list$dropped_markers_apriori) else 0
  
  totals <- c(qc_list$n_row_init, qc_list$n_col_init,
              n_samples_apriori, n_markers_apriori,
              qc_list$n_row_dropped, qc_list$n_col_dropped, qc_list$n_col_zerovar,
              qc_list$n_row_final, qc_list$n_col_final)
  
  # 3. Build Dynamic Dataframe
  summary_df <- data.frame(Metric = metrics, Total = totals, stringsAsFactors = FALSE)
  
  # 4. Integrate Group Breakdowns (Dynamic Columns)
  # Extract unique group names from initial breakdown
  if (!is.null(qc_list$breakdown_init)) {
    group_names <- names(qc_list$breakdown_init)
    
    # Initialize group columns with NA/0
    for (g in group_names) {
      summary_df[[g]] <- NA
    }
    
    # Helper to map counts to the dataframe row
    fill_row_counts <- function(df, metric_name, counts) {
      if (is.null(counts)) return(df)
      row_idx <- which(df$Metric == metric_name)
      if (length(row_idx) == 0) return(df)
      
      for (g in names(counts)) {
        if (g %in% colnames(df)) {
          df[row_idx, g] <- as.numeric(counts[g])
        }
      }
      return(df)
    }
    
    # Fill Initial Samples
    summary_df <- fill_row_counts(summary_df, "Initial Samples", qc_list$breakdown_init)
    
    # Fill Final Samples
    summary_df <- fill_row_counts(summary_df, "Final Samples", qc_list$breakdown_final)
    
    # Calculate and Fill Dropped (QC) implicitly if needed, or leave NA for markers
    # For now, we show explicit counts where available.
  }
  
  writeData(wb, "Summary", summary_df)
  
  # --- Sheet 2: Detailed Dropped Items ---
  addWorksheet(wb, "Details_Dropped")
  
  curr_row <- 1
  
  # Write Dropped Samples (A Priori)
  if (!is.null(qc_list$dropped_samples_apriori) && nrow(qc_list$dropped_samples_apriori) > 0) {
    writeData(wb, "Details_Dropped", "Samples Excluded by Config (A Priori - Blacklist):", startRow = curr_row)
    curr_row <- curr_row + 1
    writeData(wb, "Details_Dropped", qc_list$dropped_samples_apriori, startRow = curr_row)
    curr_row <- curr_row + nrow(qc_list$dropped_samples_apriori) + 3
  }
  
  # Write Dropped Markers (A Priori)
  if (!is.null(qc_list$dropped_markers_apriori) && nrow(qc_list$dropped_markers_apriori) > 0) {
    writeData(wb, "Details_Dropped", "Markers Excluded by Config (A Priori):", startRow = curr_row)
    curr_row <- curr_row + 1
    writeData(wb, "Details_Dropped", qc_list$dropped_markers_apriori, startRow = curr_row)
    curr_row <- curr_row + nrow(qc_list$dropped_markers_apriori) + 3
  }
  
  # Write Dropped Patients (QC)
  if (nrow(qc_list$dropped_rows_detail) > 0) {
    writeData(wb, "Details_Dropped", "Dropped Patients (QC - High NA or Outliers):", startRow = curr_row)
    curr_row <- curr_row + 1
    writeData(wb, "Details_Dropped", qc_list$dropped_rows_detail, startRow = curr_row)
    curr_row <- curr_row + nrow(qc_list$dropped_rows_detail) + 3
  }
  
  # Write Dropped Markers (QC)
  if (nrow(qc_list$dropped_cols_detail) > 0) {
    writeData(wb, "Details_Dropped", "Dropped Markers (QC - > Threshold NA):", startRow = curr_row)
    curr_row <- curr_row + 1
    writeData(wb, "Details_Dropped", qc_list$dropped_cols_detail, startRow = curr_row)
  }
  
  saveWorkbook(wb, out_path, overwrite = TRUE)
  message(sprintf("[IO] QC Report saved to: %s", out_path))
}