# R/modules_viz.R
# ==============================================================================
# VISUALIZATION MODULE
# Description: Custom themes, PCA wrappers, and Plotting functions.
# Dependencies: ggplot2, FactoMineR, factoextra
# ==============================================================================

library(ggplot2)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(igraph)
library(tidygraph)
library(ggraph)
library(ComplexHeatmap)
library(circlize)

#' @title Standard Project Theme
#' @description A consistent ggplot2 theme for publication-quality figures.
#' @return A ggplot theme object.
theme_coda <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "gray95")
    )
}

#' @title Generate Dynamic Group Palette
#' @description Maps generic config colors to specific project groups defined in control/case_groups.
#' @param config The loaded configuration list.
#' @return A named vector of hex codes.
get_palette <- function(config) {
  
  # 1. Defaults if config is incomplete
  if (is.null(config$colors) || is.null(config$control_group)) {
    return(c(Control = "grey", Case = "red"))
  }
  
  # [FIX] Initialize as character vector to prevent list creation on assignment
  final_palette <- character()
  
  # 2. Assign Control Color
  ctrl_name <- config$control_group
  ctrl_color <- if(!is.null(config$colors$control)) config$colors$control else "#2E8B57"
  final_palette[[ctrl_name]] <- ctrl_color
  
  # 3. Assign Case Colors
  case_names <- config$case_groups
  case_colors <- if(!is.null(config$colors$cases)) config$colors$cases else c("#CD5C5C", "#4682B4")
  
  if (length(case_names) > 0) {
    # Interpolate if we have more groups than defined colors
    if (length(case_colors) < length(case_names)) {
      pal_func <- colorRampPalette(case_colors)
      case_colors <- pal_func(length(case_names))
    }
    
    # Map colors to Case Group Names
    for (i in seq_along(case_names)) {
      grp_name <- case_names[i]
      final_palette[[grp_name]] <- case_colors[i]
    }
  }
  
  # 4. Ensure generic "Case" key exists (for merged plots)
  if (!("Case" %in% names(final_palette))) {
    final_palette[["Case"]] <- if(!is.null(config$colors$Case)) config$colors$Case else "#CD5C5C"
  }
  
  return(final_palette)
}

#' @title Plot Merged Raw Distribution with Highlights & Stats
#' @description 
#' Visualizes raw marker distribution.
#' Features: Violin (density), Jitter (points), Mean (Dashed Red), Median (Solid Black).
#' @param data_df Dataframe containing 'Group', metadata columns, and the marker.
#' @param marker_name Name of the marker column to plot.
#' @param colors Named vector of colors for the Groups.
#' @param highlight_pattern String pattern to search in 'Original_Source' (default: "_LS").
#' @param highlight_color Color for the border of highlighted points (default: "#FFD700").
#' @return A ggplot object.
plot_raw_distribution_merged <- function(data_df, marker_name, colors, 
                                         highlight_pattern = "_LS", 
                                         highlight_color = "#FFD700") {
  
  require(ggplot2)
  require(dplyr)
  
  # 1. Validation
  if (!marker_name %in% names(data_df)) return(NULL)
  
  # 2. Prepare Data & Highlighting Logic
  src_col <- names(data_df)[grepl("Original|Source", names(data_df), ignore.case = TRUE)][1]
  
  plot_data <- data_df
  plot_data$Value <- plot_data[[marker_name]]
  
  # Define Highlight Attributes defaults
  plot_data$Highlight_Type <- "Standard"
  plot_data$Stroke_Size <- 0.2
  
  if (!is.na(src_col)) {
    is_target <- grepl(highlight_pattern, plot_data[[src_col]])
    if (any(is_target)) {
      plot_data$Highlight_Type[is_target] <- "Target"
      plot_data$Stroke_Size[is_target] <- 1.5
    }
  }
  
  # 3. Calculate Stats for Subtitle
  n_na <- sum(is.na(plot_data$Value))
  pct_na <- round((n_na / nrow(plot_data)) * 100, 1)
  
  # 4. Plotting
  p <- ggplot(plot_data, aes(x = Group, y = Value, fill = Group)) +
    
    # Layer 1: Violin
    geom_violin(alpha = 0.4, trim = FALSE, color = NA, scale = "width") +
    
    # Layer 2: Jittered Points
    geom_jitter(aes(color = Highlight_Type, stroke = Stroke_Size), 
                width = 0.2, size = 2.5, shape = 21, alpha = 0.8) +
    
    # Layer 3: Median (Solid Black Line)
    stat_summary(fun = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.5, size = 0.8, color = "black") +
    
    # Layer 4: Mean (Dashed Red Line)
    stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.5, size = 0.8, color = "darkred", linetype = "dashed") +
    
    # Colors
    scale_fill_manual(values = colors) +
    scale_color_manual(values = c("Standard" = "white", "Target" = highlight_color), 
                       guide = "none") +
    scale_continuous_identity(aesthetics = "stroke") +
    
    # Styling
    labs(
      title = paste("Raw Distribution:", marker_name),
      subtitle = sprintf("Stats: Black Line = Median | Red Dashed = Mean | Highlight: '%s'", 
                         highlight_pattern),
      y = "Raw Value",
      x = NULL
    ) +
    theme_coda() +
    theme(
      legend.position = "none", # Legend redundant as X axis labels exist
      axis.text.x = element_text(face = "bold", size = 11)
    )
  
  return(p)
}

#' @title Save Distribution Report (Multi-Page PDF)
#' @description 
#' Iterates through all markers and saves their raw distribution plots 
#' into a single PDF file. Encapsulates the looping logic.
#' 
#' @param data_df Merged dataframe with metadata and raw values.
#' @param markers Vector of marker names to plot.
#' @param file_path Output path for the PDF.
#' @param colors Named vector of colors for groups.
#' @param hl_pattern String pattern for highlighting (e.g. "_LS").
viz_save_distribution_report <- function(data_df, markers, file_path, colors, hl_pattern = "") {
  
  if (length(markers) == 0) {
    warning("[Viz] No markers provided for distribution report.")
    return(NULL)
  }
  
  message(sprintf("   [Viz] Saving distribution report to: %s", basename(file_path)))
  
  pdf(file_path, width = 8, height = 6)
  # Ensure device is closed even if the loop crashes
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)
  
  for (marker in markers) {
    if (marker %in% names(data_df)) {
      # Try-catch ensures one bad plot doesn't crash the whole PDF generation
      tryCatch({
        p <- plot_raw_distribution_merged(
          data_df = data_df, 
          marker_name = marker,
          colors = colors,     
          highlight_pattern = hl_pattern,   
          highlight_color = "#FFD700"
        )
        if (!is.null(p)) print(p)
      }, error = function(e) {
        warning(sprintf("      [WARN] Failed to plot marker '%s': %s", marker, e$message))
      })
    }
  }
  
  dev.off()
}

#' @title Run PCA on Compositional Data
#' @description Wraps FactoMineR::PCA with settings appropriate for CLR data.
#' @param data_matrix A numeric matrix (CLR transformed markers).
#' @return A PCA result object.
run_coda_pca <- function(data_matrix) {
  # Validation
  if (any(is.na(data_matrix))) stop("PCA Input contains NAs. Fix upstream.")
  
  # Run PCA
  # scale.unit = FALSE is crucial for CLR data to preserve Aitchison geometry.
  res_pca <- FactoMineR::PCA(data_matrix, scale.unit = FALSE, graph = FALSE)
  
  return(res_pca)
}

#' @title Custom PCA Biplot with Multi-Pattern Highlighting
#' @description Plots Individuals with flexible highlighting based on metadata substrings.
#' @param pca_res Result from run_coda_pca().
#' @param metadata Dataframe containing 'Patient_ID', 'Group' and the subgroup column.
#' @param colors Named vector of colors for main Groups (fill).
#' @param dims Integer vector of length 2 indicating which PCs to plot.
#' @param show_labels Logical. If TRUE, adds Patient_ID labels.
#' @param highlight_patterns Named vector (Pattern -> Color). E.g., c("_LS" = "#FFD700").
#' @param find_col_keyword Keyword to identify the metadata column to search.
#' @return A ggplot object.
plot_pca_custom <- function(pca_res, metadata, colors, dims = c(1, 2), 
                            show_labels = FALSE, 
                            highlight_patterns = NULL,
                            find_col_keyword = "Original|Source") {
  
  # 1. Setup Coordinates & Metadata
  if (length(dims) != 2) stop("dims parameter must be a vector of length 2")
  ind_coords <- as.data.frame(pca_res$ind$coord[, dims])
  colnames(ind_coords) <- c("X_Coord", "Y_Coord")
  
  if (nrow(ind_coords) != nrow(metadata)) stop("Metadata/PCA dimension mismatch.")
  plot_data <- cbind(metadata, ind_coords)
  
  # 2. Logic for Highlights
  # Default state: No Highlight
  plot_data$Highlight_Type <- "Standard"
  plot_data$Stroke_Size <- 0.5 
  
  # Identify the column to search
  target_col <- names(metadata)[grepl(find_col_keyword, names(metadata), ignore.case = TRUE)][1]
  
  has_highlights <- !is.null(highlight_patterns) && length(highlight_patterns) > 0 && !is.na(target_col)
  
  if (has_highlights) {
    target_values <- as.character(plot_data[[target_col]])
    
    for (pattern in names(highlight_patterns)) {
      matches <- grepl(pattern, target_values)
      # Only overwrite if currently Standard (prevents overwriting higher priority matches)
      to_update <- matches & (plot_data$Highlight_Type == "Standard")
      
      if (any(to_update)) {
        plot_data$Highlight_Type[to_update] <- pattern 
        plot_data$Stroke_Size[to_update] <- 2.0  
      }
    }
  }
  
  # Prepare Color Scale for Borders
  border_colors <- c("Standard" = "white")
  if (has_highlights) {
    border_colors <- c(highlight_patterns, border_colors)
  }
  
  # 3. Variance Explained
  eig_val <- get_eigenvalue(pca_res)
  var_x <- round(eig_val[dims[1], 2], 1)
  var_y <- round(eig_val[dims[2], 2], 1)
  
  # 4. Plotting
  p <- ggplot(plot_data, aes(x = X_Coord, y = Y_Coord, fill = Group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray80") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray80") +
    stat_ellipse(geom = "polygon", alpha = 0.1, show.legend = FALSE, level = 0.95) +
    
    # Points 
    geom_point(aes(color = Highlight_Type, stroke = Stroke_Size), 
               size = 3.5, shape = 21) +
    
    # Main Fill Colors (Clinical Group)
    scale_fill_manual(values = colors, guide = guide_legend(override.aes = list(shape = 21))) +
    
    # Border Colors (Highlights)
    scale_color_manual(name = "Subgroup", values = border_colors) +
    scale_continuous_identity(aesthetics = "stroke") +
    
    labs(
      x = sprintf("PC%d (%s%%)", dims[1], var_x),
      y = sprintf("PC%d (%s%%)", dims[2], var_y),
      title = sprintf("PCA (PC%d vs PC%d)", dims[1], dims[2]),
      subtitle = if(has_highlights) "Colored borders indicate specific subgroups" else NULL
    ) +
    theme_coda()
  
  # Hide legend if no highlights exist
  if (!has_highlights) {
    p <- p + guides(color = "none")
  }
  
  if (show_labels) {
    p <- p + geom_text_repel(aes(label = Patient_ID), size = 3, show.legend = FALSE, max.overlaps = 15)
  }
  
  return(p)
}

#' @title Variable Contribution Plot
#' @description Wraps fviz_pca_var for consistent styling.
plot_pca_variables <- function(pca_res) {
  fviz_pca_var(pca_res, 
               col.var = "contrib", 
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE) +
    labs(title = "Marker Contribution to Variance") +
    theme_coda()
}

#' @title Plot Stratification Heatmap (ComplexHeatmap)
#' @description 
#' Generates a clustered heatmap using Z-scored data.
#' @param mat_z A Z-scored numeric matrix (Samples x Markers).
#' @param metadata Dataframe containing patient metadata (aligned with mat_z rows).
#' @param annotation_colors_list A named list of color vectors for annotations.
#' @param title Plot title.
#' @return A ComplexHeatmap object (drawn).
plot_stratification_heatmap <- function(mat_z, metadata, annotation_colors_list, title = "Stratification") {
  
  # 1. Prepare Data
  mat_plot <- t(mat_z)
  
  if (!all(rownames(mat_z) == metadata$Patient_ID)) {
    stop("Mismatch between matrix rownames and metadata Patient_ID")
  }
  
  # [FIX] Ensure metadata is a pure data.frame (not tibble) for HeatmapAnnotation
  metadata_clean <- as.data.frame(metadata)
  
  # 2. Setup Annotations
  ha <- HeatmapAnnotation(
    df = metadata_clean,
    col = annotation_colors_list,
    simple_anno_size = unit(0.3, "cm"),
    show_annotation_name = TRUE
  )
  
  # 3. Define Color Map for Z-Scores
  col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B"))
  
  # 4. Draw Heatmap
  hm <- Heatmap(
    mat_plot,
    name = "Z-Score",
    col = col_fun,
    
    # Clustering
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_distance_columns = "euclidean",
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    
    # Annotations
    top_annotation = ha,
    
    # Aesthetics
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 6),
    show_column_names = FALSE, 
    row_dend_width = unit(2, "cm"),
    column_dend_height = unit(2, "cm"),
    
    column_title = paste0(title, " (n=", ncol(mat_plot), ")"),
    row_title = "Markers"
  )
  
  return(hm)
}

#' @title Generate Complex Heatmap Annotation Colors
#' @description 
#' Handles the logic for assigning colors to metadata columns in heatmaps.
#' Implements inheritance (Subgroup inherits from Group) and highlighting rules.
#' 
#' @param metadata Dataframe of metadata.
#' @param group_col Name of the main grouping column.
#' @param base_colors Named vector of base group colors.
#' @param hl_pattern Pattern to search for special highlighting.
#' @param hl_color Color to use for highlighted items.
#' @return A list of named vectors suitable for ComplexHeatmap.
viz_generate_complex_heatmap_colors <- function(metadata, group_col, base_colors, 
                                                hl_pattern = NULL, hl_color = "#FFD700") {
  
  annotation_colors <- list()
  
  # 1. Base Group Colors
  present_groups <- unique(metadata[[group_col]])
  # Intersect to avoid errors if config has more colors than data
  group_pal <- base_colors[intersect(names(base_colors), present_groups)]
  annotation_colors[[group_col]] <- group_pal
  
  # 2. Extra Columns Logic
  extra_cols <- setdiff(colnames(metadata), group_col)
  
  for (col_name in extra_cols) {
    # Get unique non-NA values
    raw_vals <- as.character(metadata[[col_name]])
    vals <- sort(unique(raw_vals[!is.na(raw_vals)]))
    
    col_pal <- character() # Force named vector
    
    for (v in vals) {
      # Priority 1: Exact Match to a known Group Name -> Inherit Color
      if (v %in% names(group_pal)) {
        col_pal[[v]] <- group_pal[[v]]
        
        # Priority 2: Matches Highlight Pattern -> Use Highlight Color
      } else if (!is.null(hl_pattern) && hl_pattern != "" && grepl(hl_pattern, v)) {
        col_pal[[v]] <- hl_color
        
        # Priority 3: Parent Group Containment (e.g. "NSCLC_LS" contains "NSCLC")
      } else {
        parent_match <- NA
        for (grp in names(group_pal)) {
          if (grepl(grp, v)) {
            parent_match <- grp
            break
          }
        }
        # Priority 4: Fallback Gray
        col_pal[[v]] <- if (!is.na(parent_match)) group_pal[[parent_match]] else "#D3D3D3"
      }
    }
    annotation_colors[[col_name]] <- col_pal
  }
  
  return(annotation_colors)
}

#' @title Plot Signature Boxplots (Top Drivers)
#' @description 
#' Helper function to visualize the distribution of the top discriminating markers.
#' Draws boxplots of Z-scores for the features with highest loadings.
#' 
#' @param data_matrix Numeric matrix of data used in PLS (usually Z-scored).
#' @param group_factor Factor vector defining groups.
#' @param loadings_df Dataframe of loadings from extract_plsda_loadings.
#' @param comp Integer. Which component to visualize?
#' @param n_top Integer. How many top markers to plot?
#' @param colors Named vector of colors.
#' @return A ggplot object.
viz_plot_signature_boxplots <- function(data_matrix, group_factor, loadings_df, 
                                        comp = 1, n_top = 6, colors) {
  
  col_name <- paste0("Comp", comp, "_Weight")
  if (!col_name %in% names(loadings_df)) return(NULL)
  
  # Filter top absolute weights
  top_mks <- loadings_df %>%
    arrange(desc(abs(!!sym(col_name)))) %>%
    head(n_top) %>%
    pull(Marker)
  
  if (length(top_mks) == 0) return(NULL)
  
  # Prepare long dataframe for plotting
  plot_df <- as.data.frame(data_matrix[, top_mks, drop=FALSE])
  plot_df$Group <- group_factor
  
  long_df <- plot_df %>%
    pivot_longer(-Group, names_to = "Marker", values_to = "Z_Score") %>%
    mutate(Marker = factor(Marker, levels = top_mks)) # Preserve order
  
  p <- ggplot(long_df, aes(x = Group, y = Z_Score, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
    facet_wrap(~Marker, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = colors) +
    labs(
      title = sprintf("Top %d Drivers Distribution (Comp %d)", n_top, comp),
      subtitle = "Visual validation of PLS-DA weights (Z-Scores)",
      y = "Standardized Abundance (Z)", x = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", strip.background = element_rect(fill="gray95"))
  
  return(p)
}

#' @title Report sPLS-DA Visualization 
#' @description 
#' Generates a comprehensive PDF report for sPLS-DA results.
#' Includes robust error handling for PDF device closing.
#' 
#' @param pls_res The result object from run_splsda_model.
#' @param drivers_df Dataframe of extracted loadings.
#' @param metadata_viz Dataframe for plotting.
#' @param colors_viz Named vector of colors for groups.
#' @param out_path Path to save the PDF.
#' @param group_col The name of the metadata column to use for grouping/coloring (default: "Group").
viz_report_plsda <- function(pls_res, drivers_df, metadata_viz, colors_viz, out_path, group_col = "Group") {
  
  if (is.null(pls_res$model)) return(NULL)
  
  requireNamespace("mixOmics", quietly = TRUE)
  require(ggplot2)
  require(dplyr)
  require(tidyr)
  
  message(sprintf("   [Viz] Generating sPLS-DA graphical report: %s", basename(out_path)))
  
  # Validate column existence immediately
  if (!group_col %in% names(metadata_viz)) {
    stop(sprintf("Column '%s' not found in metadata for sPLS-DA visualization.", group_col))
  }
  
  # [FIX] Open PDF and ensure it closes using on.exit immediately
  pdf(out_path, width = 11, height = 8)
  # This ensures dev.off() is called even if the code below crashes
  on.exit(try(dev.off(), silent = TRUE), add = TRUE) 
  
  # Setup data for plotting
  # Use the dynamic column specified by group_col
  plot_group_factor <- as.factor(metadata_viz[[group_col]])
  
  # Map colors safely
  # We subset colors_viz by the levels present in the data to avoid mismatches
  levels_present <- levels(plot_group_factor)
  plot_colors <- colors_viz[levels_present]
  
  # Fallback for missing colors
  if (any(is.na(plot_colors))) {
    missing_grps <- levels_present[is.na(plot_colors)]
    warning(paste("Missing colors for:", paste(missing_grps, collapse=", ")))
    # Assign gray to missing colors
    plot_colors[is.na(plot_colors)] <- "gray50"
    names(plot_colors) <- levels_present
  }
  
  n_comps <- pls_res$model$ncomp
  
  # 1. Indiv Biplot
  # Note: mixOmics plotIndiv usually expects 'group' as a vector matching X rows
  tryCatch({
    mixOmics::plotIndiv(pls_res$model, 
                        comp = c(1,2), 
                        group = plot_group_factor, 
                        col.per.group = plot_colors, 
                        ellipse = TRUE, 
                        legend = TRUE, 
                        title = "sPLS-DA: Model Separation",
                        subtitle = paste("Grouping by:", group_col),
                        star.plot = TRUE)
  }, error = function(e) {
    plot.new(); text(0.5, 0.5, paste("Biplot Error:", e$message))
  })
  
  # Iterate over components for Loadings
  for (i in 1:n_comps) {
    
    comp_col <- paste0("Comp", i, "_Weight")
    if (!comp_col %in% names(drivers_df)) next
    
    # Auto-Detection Logic for Association Direction
    variates <- pls_res$model$variates$X[, i]
    stat_groups <- plot_group_factor
    group_means <- tapply(variates, stat_groups, mean)
    
    pos_group_name <- names(group_means)[which.max(group_means)]
    neg_group_name <- names(group_means)[which.min(group_means)]
    
    label_sub <- sprintf("Direction: Positive -> %s | Negative -> %s", 
                         pos_group_name, neg_group_name)
    
    # Prepare fill colors for the bar plot
    current_fill_colors <- c()
    if(pos_group_name %in% names(plot_colors)) current_fill_colors[[pos_group_name]] <- plot_colors[[pos_group_name]]
    if(neg_group_name %in% names(plot_colors)) current_fill_colors[[neg_group_name]] <- plot_colors[[neg_group_name]]
    
    # Fill missing with gray if necessary
    if(is.null(current_fill_colors[[pos_group_name]])) current_fill_colors[[pos_group_name]] <- "gray"
    if(is.null(current_fill_colors[[neg_group_name]])) current_fill_colors[[neg_group_name]] <- "gray"
    
    df_comp <- drivers_df[abs(drivers_df[[comp_col]]) > 0, ]
    
    if (nrow(df_comp) > 0) {
      # Determine association based on weight sign
      df_comp$Association <- ifelse(df_comp[[comp_col]] > 0, pos_group_name, neg_group_name)
      
      # 2. Loading Plot 
      p_load <- ggplot(df_comp, aes(x = reorder(Marker, abs(.data[[comp_col]])), 
                                    y = .data[[comp_col]], 
                                    fill = Association)) +
        geom_bar(stat = "identity", width = 0.7) +
        coord_flip() +
        scale_fill_manual(values = current_fill_colors) +
        labs(title = sprintf("sPLS-DA Loadings (Component %d)", i), 
             subtitle = label_sub,
             x = "Marker", y = "Weight Contribution",
             fill = "Associated Group") +
        theme_coda() + theme(legend.position = "bottom")
      print(p_load)
      
      # 3. Signature Boxplots (Helper function assumed to be present)
      if(exists("viz_plot_signature_boxplots")) {
        p_box <- viz_plot_signature_boxplots(
          data_matrix = pls_res$model$X, 
          group_factor = plot_group_factor,
          loadings_df = df_comp, 
          comp = i, 
          n_top = 9, 
          colors = plot_colors
        )
        if (!is.null(p_box)) print(p_box)
      }
    }
  }
  
  # 4. Clustered Image Map (CIM)
  if (nrow(drivers_df) > 1) {
    tryCatch({
      # Map row side colors using the correctly matched vectors
      groups_vec <- as.character(metadata_viz[[group_col]])
      row_cols <- plot_colors[groups_vec]
      
      # Check for NAs in colors before plotting
      row_cols[is.na(row_cols)] <- "grey"
      
      mixOmics::cim(pls_res$model, 
                    row.sideColors = row_cols,
                    title = "Global Signature Heatmap (CIM)",
                    margins = c(7, 7),
                    save = NULL) 
    }, error = function(e) {
      plot.new(); text(0.5, 0.5, paste("CIM Error:", e$message))
    })
  }
  
  # dev.off() is handled by on.exit() automatically
}

#' @title Save sPLS-DA Loadings Report (Multi-Page) - FIXED & DECOUPLED
#' @param bar_colors Named vector with 'positive' and 'negative'.
viz_save_plsda_loadings <- function(loadings_df, n_comp, tuning_df = NULL, file_path,
                                    bar_colors = c(positive="#CD5C5C", negative="#4682B4")) {
  
  if (nrow(loadings_df) == 0) return(NULL)
  
  message(sprintf("   [Viz] Saving sPLS-DA loadings (%d comps) to: %s", n_comp, basename(file_path)))
  
  # Ensure colors are named correctly for manual scale
  fill_values <- c("Positive_Assoc" = unname(bar_colors["positive"]), 
                   "Negative_Assoc" = unname(bar_colors["negative"]))
  
  pdf(file_path, width = 8, height = 6)
  
  for (i in 1:n_comp) {
    comp_col <- paste0("Comp", i, "_Weight")
    if (!comp_col %in% names(loadings_df)) {
      warning(sprintf("      [WARN] Column '%s' not found. Skipping.", comp_col))
      next
    }
    
    df_comp <- loadings_df[abs(loadings_df[[comp_col]]) > 0, ]
    df_comp$Importance <- abs(df_comp[[comp_col]])
    df_comp <- df_comp[order(df_comp$Importance, decreasing = TRUE), ]
    df_comp$Direction <- ifelse(df_comp[[comp_col]] > 0, "Positive_Assoc", "Negative_Assoc")
    
    # [FIX] Safer extraction of selected feature count
    n_selected <- NA
    if (!is.null(tuning_df) && nrow(tuning_df) > 0) {
      # Handle if tuning_df is vector or dataframe
      if (ncol(tuning_df) >= i) {
        # Using [1, i] assumes 1 row and columns are components. 
        # Using safely [[i]] if it's a list/vector structure hidden in DF
        val <- tryCatch(as.numeric(tuning_df[1, i]), error = function(e) NA)
        if (!is.na(val)) n_selected <- val
      }
    }
    
    title_suffix <- if(!is.na(n_selected)) sprintf("(Selected Features: %d)", n_selected) else ""
    
    if (nrow(df_comp) > 0) {
      p <- ggplot(df_comp, aes(x = reorder(Marker, abs(df_comp[[comp_col]])), 
                               y = df_comp[[comp_col]], 
                               fill = Direction)) +
        geom_bar(stat = "identity", width = 0.7) +
        coord_flip() +
        scale_fill_manual(values = fill_values) +
        labs(
          title = sprintf("sPLS-DA Signature: Component %d", i),
          subtitle = title_suffix,
          x = "Marker", y = "Loading Weight"
        ) +
        theme_coda() + theme(legend.position = "bottom")
      print(p)
    } else {
      plot.new(); text(0.5, 0.5, sprintf("No features selected for Component %d", i))
    }
  }
  dev.off()
}

#' @title Plot Network Graph (ggraph)
#' @description 
#' Visualizes the inferred network.
plot_network_structure <- function(adj_mat, weight_mat, title = "Network", 
                                   layout_type = "nicely", min_cor = 0) { 
  
  # 1. Build Graph Object
  g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
  
  # Add attributes
  E(g)$weight_raw <- NA
  E(g)$sign <- NA
  E(g)$weight <- NA 
  
  # Get edges
  el <- igraph::as_data_frame(g, what = "edges")
  
  if (nrow(el) > 0) {
    weights <- numeric(nrow(el))
    signs   <- character(nrow(el))
    keep_edge <- logical(nrow(el)) 
    
    for(k in 1:nrow(el)) {
      w <- weight_mat[el[k,1], el[k,2]]
      
      if (abs(w) >= min_cor) {
        keep_edge[k] <- TRUE
        weights[k] <- abs(w)
        signs[k]   <- ifelse(w > 0, "Positive", "Negative")
      } else {
        keep_edge[k] <- FALSE
      }
    }
    
    edges_to_remove <- E(g)[!keep_edge]
    g <- delete_edges(g, edges_to_remove)
    
    if (ecount(g) > 0) {
      E(g)$weight <- weights[keep_edge]
      E(g)$sign   <- signs[keep_edge]
    }
  }
  
  V(g)$degree <- igraph::degree(g)
  
  # 2. Plotting (Standard ggraph logic)
  tg <- tidygraph::as_tbl_graph(g)
  
  p <- ggraph(tg, layout = layout_type) + 
    geom_edge_link(aes(width = weight, color = sign), alpha = 0.6) +
    scale_edge_width(range = c(0.2, 1.5), guide = "none") +
    scale_edge_color_manual(values = c("Positive" = "#00BFC4", "Negative" = "#F8766D")) +
    
    geom_node_point(aes(size = degree), color = "gray20", fill = "white", shape=21) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = 20) +
    
    scale_size_continuous(range = c(2, 8)) +
    labs(title = title, 
         subtitle = sprintf("Nodes: %d | Edges: %d (Filter: |rho| > %.2f)", 
                            vcount(g), ecount(g), min_cor),
         edge_color = "Association") +
    theme_void() +
    theme(
      plot.title = element_text(face="bold", hjust=0.5),
      plot.subtitle = element_text(hjust=0.5, color="gray50"),
      legend.position = "bottom"
    )
  
  return(p)
}
