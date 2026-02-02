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
#' @description Generates a multi-sheet Excel report with filtering statistics, group breakdowns, and categorized marker lists.
#' @param qc_list A list containing summary stats, dataframes of dropped items, and group mappings.
#' @param out_path Path to save the .xlsx file.
#' @param config Configuration object (optional) containing 'qc_reporting' settings for marker categorization.
save_qc_report <- function(qc_list, out_path, config = NULL) {
  
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
  
  # 5. Calculate Totals per Macro-Group and add a "Total (by Group)" row
  if (!is.null(qc_list$group_mapping) && !is.null(qc_list$breakdown_final)) {
    
    # Create the new row structure
    total_row <- df_samples[1, ]
    total_row[,] <- NA
    total_row$Metric <- "Total (by Group)"
    total_row$Total <- df_samples[df_samples$Metric == "Final Samples", "Total"] # Grand Total
    
    # Identify unique parent groups
    unique_parents <- unique(qc_list$group_mapping$Group)
    
    for (pg in unique_parents) {
      # Find subgroups belonging to this parent
      subgroups <- qc_list$group_mapping$Subgroup[qc_list$group_mapping$Group == pg]
      
      # Intersect with columns present in the dataframe
      valid_cols <- intersect(subgroups, colnames(df_samples))
      
      if (length(valid_cols) > 0) {
        # Calculate sum of final samples for these columns
        final_row_idx <- which(df_samples$Metric == "Final Samples")
        
        # Safe sum handling NAs
        vals <- as.numeric(df_samples[final_row_idx, valid_cols])
        group_sum <- sum(vals, na.rm = TRUE)
        
        # Identify the LAST column index among the valid columns for this group
        # This places the total under the rightmost column of the group (e.g., HNSCC_LS)
        col_indices <- match(valid_cols, colnames(df_samples))
        target_col_idx <- max(col_indices)
        target_col_name <- colnames(df_samples)[target_col_idx]
        
        total_row[[target_col_name]] <- group_sum
      }
    }
    
    # Append the row
    df_samples <- rbind(df_samples, total_row)
  }
  
  # 6. Write to Sheet (Samples Section)
  curr_row <- 1
  writeData(wb, "Summary", "SAMPLES METRICS:", startRow = curr_row)
  addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1)
  curr_row <- curr_row + 1
  
  writeData(wb, "Summary", df_samples, startRow = curr_row)
  curr_row <- curr_row + nrow(df_samples) + 3 
  
  # 7. Write to Sheet (Markers Section)
  writeData(wb, "Summary", "MARKERS METRICS:", startRow = curr_row)
  addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1)
  curr_row <- curr_row + 1
  
  # Keep track of where the markers table starts
  marker_start_row <- curr_row
  writeData(wb, "Summary", df_markers, startRow = marker_start_row)
  
  # 8. List Final Markers (Organized by Biological Categories)
  if (!is.null(qc_list$final_markers_names) && length(qc_list$final_markers_names) > 0) {
    
    final_mks <- qc_list$final_markers_names
    formatted_df <- NULL
    
    # Check if category config exists
    if (!is.null(config) && !is.null(config$qc_reporting$marker_categories)) {
      
      cats <- config$qc_reporting$marker_categories
      
      # Create list of vectors for the dataframe
      df_list <- list()
      categorized_markers <- c()
      
      for (cat_name in names(cats)) {
        # Filter markers belonging to this category that exist in the data
        expected <- cats[[cat_name]]
        present <- intersect(expected, final_mks)
        
        # Keep track of used markers
        categorized_markers <- c(categorized_markers, present)
        
        if (length(present) > 0) {
          df_list[[cat_name]] <- present
        } else {
          df_list[[cat_name]] <- character(0) 
        }
      }
      
      # Handle "Others" (Markers in data but not in config)
      others <- setdiff(final_mks, categorized_markers)
      if (length(others) > 0) {
        df_list[["Others"]] <- others
      }
      
      # Convert list to dataframe (pad with empty strings)
      max_len <- max(sapply(df_list, length))
      if (max_len > 0) {
        formatted_df <- data.frame(lapply(df_list, function(x) {
          c(x, rep("", max_len - length(x)))
        }), stringsAsFactors = FALSE)
      }
      
    } else {
      # Fallback to simple matrix if no config provided
      n_per_row <- 5
      n_needed <- ceiling(length(final_mks) / n_per_row) * n_per_row
      mks_padded <- c(final_mks, rep("", n_needed - length(final_mks)))
      formatted_df <- as.data.frame(matrix(mks_padded, ncol = n_per_row, byrow = TRUE))
      colnames(formatted_df) <- paste0("Col_", 1:n_per_row)
    }
    
    if (!is.null(formatted_df)) {
      # Determine placement: Right of the "Total" column (Metrics table)
      # We insert it to the right of the Markers metrics table
      insert_row <- marker_start_row 
      insert_col <- 3 + ncol(df_markers) 
      
      writeData(wb, "Summary", "FINAL MARKERS LIST:", startRow = insert_row, startCol = insert_col)
      addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = insert_row, cols = insert_col)
      
      writeData(wb, "Summary", formatted_df, startRow = insert_row + 1, startCol = insert_col)
    }
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