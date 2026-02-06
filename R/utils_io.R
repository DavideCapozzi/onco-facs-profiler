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
  
  # 1. Data Preparation & Calculation
  
  # Extract Base Counts per Subgroup
  # Final Counts
  final_counts <- if(!is.null(qc_list$breakdown_final)) as.numeric(qc_list$breakdown_final) else numeric(0)
  if(!is.null(qc_list$breakdown_final)) names(final_counts) <- names(qc_list$breakdown_final)
  
  # Dropped QC Counts
  drop_qc_df <- qc_list$dropped_rows_detail
  drop_qc_counts <- numeric(0)
  if (!is.null(drop_qc_df) && nrow(drop_qc_df) > 0 && "Original_Source" %in% names(drop_qc_df)) {
    tbl <- table(drop_qc_df$Original_Source)
    drop_qc_counts <- as.numeric(tbl)
    names(drop_qc_counts) <- names(tbl)
  }
  
  # Dropped A Priori Counts
  drop_pre_df <- qc_list$dropped_samples_apriori
  drop_pre_counts <- numeric(0)
  if (!is.null(drop_pre_df) && nrow(drop_pre_df) > 0 && "Original_Source" %in% names(drop_pre_df)) {
    tbl <- table(drop_pre_df$Original_Source)
    drop_pre_counts <- as.numeric(tbl)
    names(drop_pre_counts) <- names(tbl)
  }
  
  # Calculate Derived Metrics per Subgroup
  all_subgroups <- unique(c(names(final_counts), names(drop_qc_counts), names(drop_pre_counts)))
  all_subgroups <- sort(all_subgroups)
  
  # Structure: list(SubgroupName = list(Init, DropPre, DropQC, Final))
  subgroup_stats <- list()
  
  for(sg in all_subgroups) {
    v_fin <- if(sg %in% names(final_counts)) final_counts[[sg]] else 0
    v_qc  <- if(sg %in% names(drop_qc_counts)) drop_qc_counts[[sg]] else 0
    v_pre <- if(sg %in% names(drop_pre_counts)) drop_pre_counts[[sg]] else 0
    v_init <- v_fin + v_qc + v_pre
    
    subgroup_stats[[sg]] <- list(Initial = v_init, DropPre = v_pre, DropQC  = v_qc, Final   = v_fin)
  }
  
  # Calculate Grand Totals
  total_final <- sum(sapply(subgroup_stats, function(x) x$Final))
  total_qc    <- sum(sapply(subgroup_stats, function(x) x$DropQC))
  total_pre   <- sum(sapply(subgroup_stats, function(x) x$DropPre))
  total_init  <- total_final + total_qc + total_pre
  
  # 2. Build Table 1: Detailed Dropped Metrics (Subgroup Level)
  
  metrics_detailed <- c("Initial Samples", 
                        "Dropped Samples (A Priori)", 
                        "Dropped Samples (QC)", 
                        "Final Samples")
  
  totals_vector <- c(total_init, -total_pre, -total_qc, total_final)
  
  df_detailed <- data.frame(Metric = metrics_detailed, Total = totals_vector, stringsAsFactors = FALSE)
  
  for (sg in all_subgroups) {
    s <- subgroup_stats[[sg]]
    df_detailed[[sg]] <- c(s$Initial, -s$DropPre, -s$DropQC, s$Final)
  }
  
  # 3. Build Table 2: By Group Dropping Metrics (Macro/Parent Level)
  
  # Identify Parent Groups
  if(!is.null(qc_list$group_mapping)) {
    unique_parents <- unique(qc_list$group_mapping$Group)
    unique_parents <- sort(unique_parents) 
  } else {
    unique_parents <- character(0)
  }
  
  # Identify Control vs Case Parents
  is_control_func <- function(g_name) grepl("Healthy|Control", g_name, ignore.case = TRUE)
  ctrl_parents <- unique_parents[sapply(unique_parents, is_control_func)]
  case_parents <- unique_parents[!sapply(unique_parents, is_control_func)]
  
  # Helper: Aggregate metrics for a list of parents
  get_aggregated_macro_stats <- function(parent_list, mapping, sub_stats) {
    sum_init <- 0
    sum_drop <- 0
    sum_fin  <- 0
    
    # Get all subgroups belonging to these parents
    relevant_subgroups <- mapping$Subgroup[mapping$Group %in% parent_list]
    relevant_subgroups <- intersect(relevant_subgroups, names(sub_stats))
    
    for(sg in relevant_subgroups) {
      s <- sub_stats[[sg]]
      sum_init <- sum_init + s$Initial
      sum_drop <- sum_drop - (s$DropPre + s$DropQC) # Sum and ensure negative
      sum_fin  <- sum_fin + s$Final
    }
    return(c(sum_init, sum_drop, sum_fin))
  }
  
  # Initialize DataFrame
  metrics_macro <- c("Initial", "Dropped Samples", "Final")
  total_dropped_combined <- -(total_pre + total_qc)
  
  df_macro <- data.frame(Metric = metrics_parent <- metrics_macro, 
                         Total = c(total_init, total_dropped_combined, total_final), 
                         stringsAsFactors = FALSE)
  
  # Add Control Column (Aggregated)
  if(length(ctrl_parents) > 0) {
    # Naming logic: Healthy_Donors_tot + ...
    col_name <- paste0(ctrl_parents, "_tot", collapse = " + ")
    df_macro[[col_name]] <- get_aggregated_macro_stats(ctrl_parents, qc_list$group_mapping, subgroup_stats)
  }
  
  # Add Case Column (Aggregated)
  if(length(case_parents) > 0) {
    # Naming logic: HNSCC_tot + NSCLC_tot
    col_name <- paste0(case_parents, "_tot", collapse = " + ")
    df_macro[[col_name]] <- get_aggregated_macro_stats(case_parents, qc_list$group_mapping, subgroup_stats)
  }
  
  # 4. Writing to Excel
  
  curr_row <- 1
  
  # Section 1: Detailed Dropped Metrics
  writeData(wb, "Summary", "DETAILED DROPPED METRICS:", startRow = curr_row)
  addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1)
  curr_row <- curr_row + 1
  
  writeData(wb, "Summary", df_detailed, startRow = curr_row)
  addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1:ncol(df_detailed), gridExpand = TRUE)
  
  curr_row <- curr_row + nrow(df_detailed) + 2 
  
  # Section 2: By Group Dropping Metrics (Macro Level)
  writeData(wb, "Summary", "BY GROUP DROPPING METRICS:", startRow = curr_row)
  addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1)
  curr_row <- curr_row + 1
  
  writeData(wb, "Summary", df_macro, startRow = curr_row)
  addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1:ncol(df_macro), gridExpand = TRUE)
  
  curr_row <- curr_row + nrow(df_macro) + 2
  
  # Section 3: Markers Metrics
  n_markers_apriori <- if(!is.null(qc_list$dropped_markers_apriori)) nrow(qc_list$dropped_markers_apriori) else 0
  totals_markers <- c(qc_list$n_col_init,
                      n_markers_apriori,
                      qc_list$n_col_dropped,
                      qc_list$n_col_zerovar,
                      qc_list$n_col_final)
  metrics_markers <- c("Initial Markers", "Dropped Markers (A Priori)", "Dropped Markers (High NA)", "Dropped Markers (Zero Var)", "Final Markers")
  df_markers <- data.frame(Metric = metrics_markers, Total = totals_markers, stringsAsFactors = FALSE)
  
  writeData(wb, "Summary", "MARKERS METRICS:", startRow = curr_row)
  addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1)
  curr_row <- curr_row + 1
  writeData(wb, "Summary", df_markers, startRow = curr_row)
  curr_row <- curr_row + nrow(df_markers) + 2 
  
  # Section 4: Final Markers List
  if (!is.null(qc_list$final_markers_names) && length(qc_list$final_markers_names) > 0) {
    final_mks <- qc_list$final_markers_names
    formatted_df <- NULL
    
    if (!is.null(config) && !is.null(config$qc_reporting$marker_categories)) {
      cats <- config$qc_reporting$marker_categories
      df_list <- list()
      categorized_markers <- c()
      
      for (cat_name in names(cats)) {
        expected <- cats[[cat_name]]
        present <- intersect(expected, final_mks)
        categorized_markers <- c(categorized_markers, present)
        df_list[[cat_name]] <- if (length(present) > 0) present else character(0)
      }
      
      others <- setdiff(final_mks, categorized_markers)
      if (length(others) > 0) df_list[["Others"]] <- others
      
      max_len <- max(sapply(df_list, length))
      if (max_len > 0) {
        formatted_df <- data.frame(lapply(df_list, function(x) {
          c(x, rep("", max_len - length(x)))
        }), stringsAsFactors = FALSE)
      }
    } else {
      n_per_row <- 5
      n_needed <- ceiling(length(final_mks) / n_per_row) * n_per_row
      mks_padded <- c(final_mks, rep("", n_needed - length(final_mks)))
      formatted_df <- as.data.frame(matrix(mks_padded, ncol = n_per_row, byrow = TRUE))
      colnames(formatted_df) <- paste0("Col_", 1:n_per_row)
    }
    
    if (!is.null(formatted_df)) {
      writeData(wb, "Summary", "FINAL MARKERS LIST:", startRow = curr_row)
      addStyle(wb, "Summary", createStyle(textDecoration = "bold"), rows = curr_row, cols = 1)
      writeData(wb, "Summary", formatted_df, startRow = curr_row + 1, colNames = FALSE)
    }
  }
  
  # Sheet 2: Details Dropped
  addWorksheet(wb, "Details_Dropped")
  curr_row_det <- 1
  dropped_patients_all <- data.frame()
  if (!is.null(qc_list$dropped_samples_apriori) && nrow(qc_list$dropped_samples_apriori) > 0) {
    dropped_patients_all <- rbind(dropped_patients_all, qc_list$dropped_samples_apriori)
  }
  if (!is.null(qc_list$dropped_rows_detail) && nrow(qc_list$dropped_rows_detail) > 0) {
    dropped_patients_all <- rbind(dropped_patients_all, qc_list$dropped_rows_detail)
  }
  
  if (nrow(dropped_patients_all) > 0) {
    cols_order <- intersect(c("Patient_ID", "NA_Percent", "Reason", "Original_Source"), names(dropped_patients_all))
    cols_order <- c(cols_order, setdiff(names(dropped_patients_all), cols_order))
    dropped_patients_all <- dropped_patients_all[, cols_order, drop=FALSE]
    
    writeData(wb, "Details_Dropped", "Dropped Samples:", startRow = curr_row_det)
    addStyle(wb, "Details_Dropped", createStyle(textDecoration = "bold"), rows = curr_row_det, cols = 1)
    writeData(wb, "Details_Dropped", dropped_patients_all, startRow = curr_row_det + 1)
    curr_row_det <- curr_row_det + nrow(dropped_patients_all) + 3
  }
  
  if (!is.null(qc_list$dropped_markers_apriori) && nrow(qc_list$dropped_markers_apriori) > 0) {
    writeData(wb, "Details_Dropped", "Markers Excluded (A Priori):", startRow = curr_row_det)
    addStyle(wb, "Details_Dropped", createStyle(textDecoration = "bold"), rows = curr_row_det, cols = 1)
    writeData(wb, "Details_Dropped", qc_list$dropped_markers_apriori, startRow = curr_row_det + 1)
    curr_row_det <- curr_row_det + nrow(qc_list$dropped_markers_apriori) + 3
  }
  
  if (nrow(qc_list$dropped_cols_detail) > 0) {
    writeData(wb, "Details_Dropped", "Dropped Markers (QC):", startRow = curr_row_det)
    addStyle(wb, "Details_Dropped", createStyle(textDecoration = "bold"), rows = curr_row_det, cols = 1)
    writeData(wb, "Details_Dropped", qc_list$dropped_cols_detail, startRow = curr_row_det + 1)
  }
  
  saveWorkbook(wb, out_path, overwrite = TRUE)
  message(sprintf("[IO] QC Report saved to: %s", out_path))
}