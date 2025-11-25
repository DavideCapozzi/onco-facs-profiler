# ==============================================================================
# FACS IMMUNOPHENOTYPING ANALYSIS PIPELINE - v3 (GROUPED PLOTS)
# ==============================================================================
# Role: Data Scientist & Immunologist
# Context: Multi-tumor comparison (HNSCC, NSCLC) vs Healthy Donors
# Updates: Explicit Logging, Dual PCA, Grouped Boxplots
# ==============================================================================

# 1. LIBRARY SETUP -------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)    
  library(readxl)       
  library(ggpubr)       
  library(rstatix)      
  library(FactoMineR)   
  library(factoextra)   
  library(ggrepel)      # Essential for labeled PCA
})

# 2. CONFIGURATION & PARAMETERS ------------------------------------------------
CONFIG <- list(
  file_path = "/mnt/c/Users/Davide/Desktop/Davide/bandi/0dottorato/projects/long survivors/DB_pulito_standardizzato.xlsx",
  
  sheets_to_read = list(
    "HNSCC" = "HNSCC", 
    "NSCLC" = "NSCLC", 
    "Healthy_Donors" = "Healthy_Donors"
  ),
  
  # Statistical Settings
  p_val_method = "wilcox.test",
  p_val_adjust = "BH",
  
  # Visualization Settings
  colors = c("Healthy_Donors" = "#2E8B57", "HNSCC" = "#CD5C5C", "NSCLC" = "#4682B4"),
  show_pca_loadings = TRUE,
  n_loadings_to_show = 8,
  
  # Data Quality Thresholds
  max_na_proportion = 0.90,      
  remove_high_na_samples = TRUE,
  
  # === NEW: MARKER GROUPS FOR BOXPLOTS ===
  # Define lists of markers to be plotted together. 
  # Ensure spelling matches Excel columns exactly.
  marker_groups = list(
    "Panel_1_General_T" = c("CD3", "CD4", "CD8TOT"),
    
    "Panel_2_Myeloid" = c("CD45", "NEUTROPHILS", "PMN", "M-MDSC", 
                          "LOX1-PMN-MDSC", "MONOCYTES"),
    
    "Panel_3_Tregs" = c("T-REG", "RESTING", "ACTIVE", "SUPPRESSIVE", "TREGCD137"),
    
    "Panel_4_Differentiation_Spectrum" = c("CM", "EFF", "EM", "NAIVE"),
    
    "Panel_5_CD137" = c("CD137TOT", "CD28", "CD137CD4", "CD137CD8", 
                  "CD137CM", "CD137EMRA", "CD137NAIVE", "CD137EFF", 
                  "CD137KI67", "CD137PD1"),
    
    "Panel_6" = c("KI67", "KI67CD4", "KI67CD8"),
    
    "Panel_7" = c("PD1", "PD1CD4", "PD1CD8", 
                  "PD1CM", "PD1EMRA", "PD1NAIVE", "PD1EFF")
  )
)


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

clean_and_convert <- function(x) {
  x_char <- as.character(x)
  x_char <- gsub(",", ".", x_char)
  x_char <- gsub("[<>]", "", x_char)
  x_char <- gsub("^(?i)(na|nd|n\\.a\\.|nan|missing)$", NA, x_char)
  suppressWarnings(as.numeric(x_char))
}

find_id_column <- function(df) {
  id_candidates <- names(df)[grepl("patient|sample|id|code|subject", names(df), ignore.case = TRUE)]
  if(length(id_candidates) == 0) return(names(df)[1])
  return(id_candidates[1])
}

# ==============================================================================
# STEP 1: IMPORT 
# ==============================================================================

import_data_multisheet <- function(path, sheet_map) {
  message("\n========================================")
  message(">> STEP 1: DATA IMPORT & HARMONIZATION")
  message("========================================\n")
  
  all_data_list <- list()
  all_columns_by_sheet <- list()
  
  for (sheet_name in names(sheet_map)) {
    group_label <- sheet_map[[sheet_name]]
    tryCatch({
      raw_df <- read_excel(path, sheet = sheet_name, col_types = "text")
      message(paste("   ✓ Loaded sheet:", sheet_name, "→", group_label, 
                    "(Rows:", nrow(raw_df), ")"))
      
      all_columns_by_sheet[[group_label]] <- names(raw_df)
      id_col <- find_id_column(raw_df)
      names(raw_df)[names(raw_df) == id_col] <- "sample_id"
      
      clean_df <- raw_df %>%
        mutate(across(-sample_id, clean_and_convert)) %>%
        mutate(group = group_label, .before = 2)
      
      all_data_list[[group_label]] <- clean_df
    }, error = function(e) warning(paste("Error reading", sheet_name, ":", e$message)))
  }
  
  common_cols <- Reduce(intersect, all_columns_by_sheet)
  common_cols <- setdiff(common_cols, find_id_column(all_data_list[[1]]))
  
  message(paste("   Shared markers found:", length(common_cols)))
  
  common_cols_with_meta <- c("sample_id", "group", common_cols)
  harmonized_list <- lapply(all_data_list, function(df) {
    cols_to_keep <- intersect(common_cols_with_meta, names(df))
    df %>% select(all_of(cols_to_keep))
  })
  
  final_df <- bind_rows(harmonized_list)
  return(final_df)
}

# ==============================================================================
# STEP 2: PREPROCESSING
# ==============================================================================

preprocess_data <- function(df) {
  message("\n========================================")
  message(">> STEP 2: PREPROCESSING & CLR")
  message("========================================\n")
  
  marker_cols <- df %>% select(where(is.numeric)) %>% colnames()
  
  # --- Filter Samples ---
  na_proportion <- rowMeans(is.na(df[marker_cols]))
  problematic_samples <- which(na_proportion > CONFIG$max_na_proportion)
  
  if(length(problematic_samples) > 0) {
    removed_ids <- df$sample_id[problematic_samples]
    message(paste("   ⚠ WARNING:", length(problematic_samples), 
                  "samples have >", CONFIG$max_na_proportion*100, "% missing"))
    
    # EXPLICIT LOGGING OF REMOVED SAMPLES
    message(paste("   -> REMOVING IDs:", paste(removed_ids, collapse = ", ")))
    
    if(CONFIG$remove_high_na_samples) {
      df <- df[-problematic_samples, ]
    }
  }
  
  # --- Imputation ---
  df_clean <- df
  for(col in marker_cols) {
    col_data <- df[[col]]
    valid_values <- col_data[col_data > 0 & !is.na(col_data)]
    
    pseudo_count <- if(length(valid_values) > 0) min(valid_values)/2 else 0.001
    
    # Log if column is entirely empty
    if(length(valid_values) == 0) {
      # message(paste("   ⚠ Marker", col, "is empty. Imputing 0.001."))
    }
    
    df_clean[[col]][is.na(df_clean[[col]]) | df_clean[[col]] == 0] <- pseudo_count
  }
  
  # --- CLR ---
  clr_data <- df_clean
  gmean <- function(x) exp(mean(log(x)))
  
  mat <- as.matrix(df_clean[, marker_cols])
  clr_mat <- t(apply(mat, 1, function(x) {
    gm <- gmean(x)
    if(is.na(gm) || gm == 0) rep(NA, length(x)) else log(x / gm)
  }))
  clr_data[, marker_cols] <- clr_mat
  
  return(list(raw = df_clean, clr = clr_data, markers = marker_cols))
}

# ==============================================================================
# STEP 3: PCA (DUAL PLOT)
# ==============================================================================

run_pca_analysis <- function(data_list) {
  message("\n========================================")
  message(">> STEP 3: PCA ANALYSIS")
  message("========================================\n")
  
  df_pca <- data_list$clr
  markers <- data_list$markers
  
  # Variance check for PCA
  col_vars <- apply(df_pca[markers], 2, var, na.rm = TRUE)
  zero_var <- names(col_vars[col_vars == 0])
  if(length(zero_var) > 0) {
    message(paste("   ⚠ Excluding zero-variance markers from PCA:", paste(zero_var, collapse=", ")))
    markers <- setdiff(markers, zero_var)
  }
  
  res.pca <- PCA(df_pca[, markers], graph = FALSE, scale.unit = FALSE)
  
  ind_coord <- as.data.frame(res.pca$ind$coord)
  ind_coord$group <- df_pca$group
  ind_coord$sample_id <- df_pca$sample_id # Needed for labels
  
  # Handle Colors dynamically
  used_groups <- unique(ind_coord$group)
  avail_colors <- CONFIG$colors[used_groups]
  if(any(is.na(names(avail_colors)))) {
    avail_colors <- setNames(rainbow(length(used_groups)), used_groups)
  }
  
  # --- PLOT 1: STANDARD (No Labels) ---
  p1 <- ggplot(ind_coord, aes(x = Dim.1, y = Dim.2, color = group)) +
    geom_point(size = 3.5, alpha = 0.8) +
    stat_ellipse(geom = "polygon", alpha = 0.1, level = 0.95, aes(fill=group), show.legend = F) +
    scale_color_manual(values = CONFIG$colors) +
    scale_fill_manual(values = CONFIG$colors) +
    labs(title = "PCA Plot (Overview)", 
         x = paste0("PC1 (", round(res.pca$eig[1,2], 1), "%)"),
         y = paste0("PC2 (", round(res.pca$eig[2,2], 1), "%)")) +
    theme_bw()
  
  # --- PLOT 2: LABELED ---
  p2 <- p1 + 
    geom_text_repel(aes(label = sample_id), size = 3, max.overlaps = 20, show.legend = FALSE) +
    labs(title = "PCA Plot (Sample IDs)")
  
  return(list(overview = p1, labeled = p2))
}

# ==============================================================================
# STEP 4: BOXPLOTS (GROUPED)
# ==============================================================================

generate_stats_plots_grouped <- function(data_list) {
  message("\n========================================")
  message(">> STEP 4: STATISTICAL BOXPLOTS (GROUPED)")
  message("========================================\n")
  
  df <- data_list$raw 
  all_markers_in_data <- data_list$markers
  
  # Variance Check (Global)
  # We identify valid markers first to avoid crashing later
  valid_markers <- df %>%
    select(all_of(all_markers_in_data)) %>%
    pivot_longer(everything(), names_to = "Marker", values_to = "Value") %>%
    group_by(Marker) %>%
    summarize(Total_Var = var(Value, na.rm = TRUE), .groups = 'drop') %>%
    filter(Total_Var > 0) %>%
    pull(Marker)
  
  skipped_markers <- setdiff(all_markers_in_data, valid_markers)
  if(length(skipped_markers) > 0) {
    message("   ⚠ SKIPPING markers with zero variance (cannot calculate stats):")
    message(paste("     ->", paste(skipped_markers, collapse = ", ")))
  }
  
  # Define Comparisons
  all_groups <- unique(df$group)
  ref_group <- "Healthy_Donors"
  if(!ref_group %in% all_groups) ref_group <- all_groups[1]
  tumor_groups <- setdiff(all_groups, ref_group)
  comparisons <- lapply(tumor_groups, function(x) c(ref_group, x))
  
  message(paste("   Comparisons defined: vs", ref_group))
  
  # === ITERATE OVER CONFIG GROUPS ===
  plot_list <- list()
  
  for (panel_name in names(CONFIG$marker_groups)) {
    
    # 1. Get markers for this panel
    target_markers <- unique(CONFIG$marker_groups[[panel_name]])
    
    # 2. Intersect with valid markers (exist in data AND have variance)
    markers_to_plot <- intersect(target_markers, valid_markers)
    
    # 3. Check if we have anything to plot
    if(length(markers_to_plot) == 0) {
      message(paste("\n   [", panel_name, "] Skipped: No valid markers found from list."))
      next
    }
    
    message(paste("\n   Processing Panel:", panel_name, "(", length(markers_to_plot), "markers )"))
    
    # 4. Prepare Data
    df_long <- df %>%
      select(sample_id, group, all_of(markers_to_plot)) %>%
      pivot_longer(cols = all_of(markers_to_plot), names_to = "Marker", values_to = "Value")
    
    # 5. Generate Plot
    # Calculate rows/cols dynamically based on number of markers to avoid tiny plots
    n_vars <- length(markers_to_plot)
    n_cols <- min(5, n_vars) 
    
    p <- ggplot(df_long, aes(x = group, y = Value, fill = group)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, size = 1, alpha = 0.4) +
      facet_wrap(~Marker, scales = "free_y", ncol = n_cols) + 
      scale_fill_manual(values = CONFIG$colors) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        legend.position = "none" # Hide legend to save space, colors are obvious
      ) +
      labs(title = paste("Panel:", panel_name), x = NULL, y = "Value")
    
    # 6. Add Statistics
    tryCatch({
      p <- p + stat_compare_means(
        comparisons = comparisons, 
        method = CONFIG$p_val_method, 
        label = "p.signif",
        hide.ns = FALSE
      )
    }, error = function(e) warning(paste("Stats failed for panel", panel_name)))
    
    plot_list[[panel_name]] <- p
  }
  
  return(plot_list)
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

message("\n╔════════════════════════════════════════╗")
message("║  FACS PIPELINE - v3 STARTED            ║")
message("╚════════════════════════════════════════╝")

tryCatch({
  # 1. Import
  raw_data <- import_data_multisheet(CONFIG$file_path, CONFIG$sheets_to_read)
  
  # 2. Preprocess
  processed_data <- preprocess_data(raw_data)
  
  # 3. PCA
  pca_results <- run_pca_analysis(processed_data)
  print(pca_results$overview)
  message("   ✓ Printed PCA Overview")
  
  print(pca_results$labeled)
  message("   ✓ Printed PCA with Labels")
  
  # 4. Grouped Boxplots
  box_plots <- generate_stats_plots_grouped(processed_data)
  
  if(length(box_plots) > 0) {
    message("\n   --- Printing Panel Plots ---")
    for(name in names(box_plots)) {
      print(box_plots[[name]])
      message(paste("   ✓ Printed:", name))
    }
  } else {
    warning("No boxplots were generated.")
  }
  
  message("\n╔════════════════════════════════════════╗")
  message("║  ✓✓✓ PIPELINE COMPLETED               ║")
  message("╚════════════════════════════════════════╝\n")
  
}, error = function(e) {
  message("\n!!! CRITICAL FAILURE !!!")
  message(e$message)
  traceback()
})
