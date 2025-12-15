# R/infrastructure.R
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
  message(sprintf("[System] Configuration loaded from %s", config_path))
  
  return(config)
}

#' @title Load Raw Data (Excel)
#' @description Iterates over sheets defined in config and merges them.
#' @param config The loaded configuration object.
#' @return A raw dataframe with Patient_ID and Group.
load_raw_data <- function(config) {
  input_file <- config$input_file
  
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
      
      # Rename first column to Patient_ID standard
      colnames(raw_tmp)[1] <- "Patient_ID"
      
      # Ensure numeric conversion (suppressing warnings for non-numeric coercion)
      clean_tmp <- raw_tmp %>%
        mutate(across(-1, ~suppressWarnings(as.numeric(as.character(.))))) %>%
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
#' @description Generates a multi-sheet Excel report with filtering statistics.
#' @param qc_list A list containing summary stats and dataframes of dropped items.
#' @param out_path Path to save the .xlsx file.
save_qc_report <- function(qc_list, out_path) {
  
  wb <- createWorkbook()
  
  # Sheet 1: General Summary
  addWorksheet(wb, "Summary")
  
  # Create a clean summary dataframe
  summary_df <- data.frame(
    Metric = c("Initial Samples", "Initial Markers",
               "Dropped Samples (High NA)", "Dropped Markers (High NA)", "Dropped Markers (Zero Var)",
               "Final Samples", "Final Markers"),
    Count = c(qc_list$n_row_init, qc_list$n_col_init,
              qc_list$n_row_dropped, qc_list$n_col_dropped, qc_list$n_col_zerovar,
              qc_list$n_row_final, qc_list$n_col_final)
  )
  writeData(wb, "Summary", summary_df)
  
  # Sheet 2: Detailed Dropped Items
  addWorksheet(wb, "Details_Dropped")
  
  curr_row <- 1
  
  # Write Dropped Patients
  if (nrow(qc_list$dropped_rows_detail) > 0) {
    writeData(wb, "Details_Dropped", "Dropped Patients (> Threshold NA):", startRow = curr_row)
    curr_row <- curr_row + 1
    writeData(wb, "Details_Dropped", qc_list$dropped_rows_detail, startRow = curr_row)
    curr_row <- curr_row + nrow(qc_list$dropped_rows_detail) + 3
  }
  
  # Write Dropped Markers
  if (nrow(qc_list$dropped_cols_detail) > 0) {
    writeData(wb, "Details_Dropped", "Dropped Markers (> Threshold NA):", startRow = curr_row)
    curr_row <- curr_row + 1
    writeData(wb, "Details_Dropped", qc_list$dropped_cols_detail, startRow = curr_row)
  }
  
  saveWorkbook(wb, out_path, overwrite = TRUE)
  message(sprintf("[IO] QC Report saved to: %s", out_path))
}