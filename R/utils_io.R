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
  
  # 1. Prepare Metrics Vectors
  metrics_samples <- c("Initial Samples", 
                       "Dropped Samples (A Priori)", 
                       "Dropped Samples (QC)", 
                       "Final Samples")
  
  metrics_markers <- c("Initial Markers", 
                       "Dropped Markers (A Priori)", 
                       "Dropped Markers (High NA)", 
                       "Dropped Markers (Zero Var)",
                       "Final Markers")
  
  # 2. Extract Counts
  n_samples_apriori <- if(!is.null(qc_list$dropped_samples_apriori)) nrow(qc_list$dropped_samples_apriori) else 0
  n_markers_apriori <- if(!is.null(qc_list$dropped_markers_apriori)) nrow(qc_list$dropped_markers_apriori) else 0
  
  totals_samples <- c(qc_list$n_row_init,
                      n_samples_apriori,
                      qc_list$n_row_dropped,
                      qc_list$n_row_final)
  
  totals_markers <- c(qc_list$n_col_init,
                      n_markers_apriori,
                      qc_list$n_col_dropped,
                      qc_list$n_col_zerovar,
                      qc_list$n_col_final)
  
  # 3. Build Dataframes
  df_samples <- data.frame(Metric = metrics_samples, Total = totals_samples, stringsAsFactors = FALSE)
  df_markers <- data.frame(Metric = metrics_markers, Total = totals_markers, stringsAsFactors = FALSE)
  
  # 4. Integrate Group Breakdowns (Dynamic Columns) 
  if (!is.null(qc_list$breakdown_init)) {
    group_names <- names(qc_list$breakdown_init)
    
    # Initialize group columns with NA
    for (g in group_names) {
      df_samples[[g]] <- NA
    }
    
    # Helper to map counts
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
    
    # Fill Sample Counts
    df_samples <- fill_row_counts(df_samples, "Initial Samples", qc_list$breakdown_init)
    df_samples <- fill_row_counts(df_samples, "Final Samples", qc_list$breakdown_final)
  }
  
  # 5. Write to Sheet (Samples Section)
  curr_row <- 1
  writeData(wb, "Summary", "SAMPLES METRICS:", startRow = curr_row)
  addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1)
  curr_row <- curr_row + 1
  
  writeData(wb, "Summary", df_samples, startRow = curr_row)
  curr_row <- curr_row + nrow(df_samples) + 3 # Add spacing (2 empty rows + header)
  
  # 6. Write to Sheet (Markers Section)
  writeData(wb, "Summary", "MARKERS METRICS:", startRow = curr_row)
  addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1)
  curr_row <- curr_row + 1
  
  # Keep track of where the markers table starts
  marker_start_row <- curr_row
  writeData(wb, "Summary", df_markers, startRow = marker_start_row)
  
  # 7. List Final Markers (Formatted 5 per row)
  if (!is.null(qc_list$final_markers_names) && length(qc_list$final_markers_names) > 0) {
    
    final_mks <- qc_list$final_markers_names
    n_per_row <- 5
    
    # Create matrix for formatted output
    n_needed <- ceiling(length(final_mks) / n_per_row) * n_per_row
    mks_padded <- c(final_mks, rep("", n_needed - length(final_mks)))
    mks_matrix <- matrix(mks_padded, ncol = n_per_row, byrow = TRUE)
    
    # Determine placement:
    # Row: The row corresponding to "Final Markers" in the df_markers table
    # Col: Immediately to the right of the "Total" column (Col 3)
    final_mk_rel_idx <- which(df_markers$Metric == "Final Markers")
    insert_row <- marker_start_row + final_mk_rel_idx # +1 for header is handled by writeData default
    insert_col <- 3 
    
    writeData(wb, "Summary", mks_matrix, startRow = insert_row, startCol = insert_col, colNames = FALSE)
  }
  
  # --- Sheet 2: Detailed Dropped Items ---
  addWorksheet(wb, "Details_Dropped")
  
  curr_row_det <- 1
  
  dropped_patients_all <- data.frame()
  
  # Combine A Priori Drops
  if (!is.null(qc_list$dropped_samples_apriori) && nrow(qc_list$dropped_samples_apriori) > 0) {
    dropped_patients_all <- rbind(dropped_patients_all, qc_list$dropped_samples_apriori)
  }
  
  # Combine QC Drops
  if (!is.null(qc_list$dropped_rows_detail) && nrow(qc_list$dropped_rows_detail) > 0) {
    dropped_patients_all <- rbind(dropped_patients_all, qc_list$dropped_rows_detail)
  }
  
  # Write the unified table
  if (nrow(dropped_patients_all) > 0) {
    cols_order <- intersect(c("Patient_ID", "NA_Percent", "Reason", "Original_Source"), names(dropped_patients_all))
    cols_order <- c(cols_order, setdiff(names(dropped_patients_all), cols_order))
    
    dropped_patients_all <- dropped_patients_all[, cols_order, drop=FALSE]
    
    writeData(wb, "Details_Dropped", "Dropped Samples (A priori & QC Filtering):", startRow = curr_row_det)
    addStyle(wb, "Details_Dropped", createStyle(textDecoration = "bold"), rows = curr_row_det, cols = 1)
    curr_row_det <- curr_row_det + 1
    writeData(wb, "Details_Dropped", dropped_patients_all, startRow = curr_row_det)
    curr_row_det <- curr_row_det + nrow(dropped_patients_all) + 3
  }
  
  # Write Dropped Markers (A Priori)
  if (!is.null(qc_list$dropped_markers_apriori) && nrow(qc_list$dropped_markers_apriori) > 0) {
    writeData(wb, "Details_Dropped", "Markers Excluded by Config (A Priori):", startRow = curr_row_det)
    addStyle(wb, "Details_Dropped", createStyle(textDecoration = "bold"), rows = curr_row_det, cols = 1)
    curr_row_det <- curr_row_det + 1
    writeData(wb, "Details_Dropped", qc_list$dropped_markers_apriori, startRow = curr_row_det)
    curr_row_det <- curr_row_det + nrow(qc_list$dropped_markers_apriori) + 3
  }
  
  # Write Dropped Markers (QC)
  if (nrow(qc_list$dropped_cols_detail) > 0) {
    writeData(wb, "Details_Dropped", "Dropped Markers (QC - > Threshold NA):", startRow = curr_row_det)
    addStyle(wb, "Details_Dropped", createStyle(textDecoration = "bold"), rows = curr_row_det, cols = 1)
    curr_row_det <- curr_row_det + 1
    writeData(wb, "Details_Dropped", qc_list$dropped_cols_detail, startRow = curr_row_det)
  }
  
  saveWorkbook(wb, out_path, overwrite = TRUE)
  message(sprintf("[IO] QC Report saved to: %s", out_path))
}