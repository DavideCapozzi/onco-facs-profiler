# R/modules_viz.R
# ==============================================================================
# VISUALIZATION MODULE (Refactored for Configurable Colors & Bug Fixes)
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
get_palette <- function(config) {
  # 1. Defaults if config is incomplete
  if (is.null(config$colors) || is.null(config$control_group)) {
    return(c(Control = "grey", Case = "red"))
  }
  
  final_palette <- character()
  
  # 2. Assign Control Color
  ctrl_name <- config$control_group
  ctrl_color <- if(!is.null(config$colors$control)) config$colors$control else "#2E8B57"
  final_palette[[ctrl_name]] <- ctrl_color
  
  # 3. Assign Case Colors
  case_names <- config$case_groups
  case_colors <- if(!is.null(config$colors$cases)) config$colors$cases else c("#CD5C5C", "#4682B4")
  
  if (length(case_names) > 0) {
    if (length(case_colors) < length(case_names)) {
      pal_func <- colorRampPalette(case_colors)
      case_colors <- pal_func(length(case_names))
    }
    for (i in seq_along(case_names)) {
      grp_name <- case_names[i]
      final_palette[[grp_name]] <- case_colors[i]
    }
  }
  
  # 4. Ensure generic "Case" key exists
  if (!("Case" %in% names(final_palette))) {
    final_palette[["Case"]] <- if(!is.null(config$colors$Case)) config$colors$Case else "#CD5C5C"
  }
  
  return(final_palette)
}

#' @title Plot Merged Raw Distribution with Highlights & Stats
#' @param stat_colors Named vector with 'mean' and 'median' colors.
plot_raw_distribution_merged <- function(data_df, marker_name, colors, 
                                         highlight_pattern = "_LS", 
                                         highlight_color = "#FFD700",
                                         stat_colors = c(mean="darkred", median="black")) {
  require(ggplot2)
  require(dplyr)
  
  if (!marker_name %in% names(data_df)) return(NULL)
  
  # Prepare Data
  src_col <- names(data_df)[grepl("Original|Source", names(data_df), ignore.case = TRUE)][1]
  plot_data <- data_df
  plot_data$Value <- plot_data[[marker_name]]
  plot_data$Highlight_Type <- "Standard"
  plot_data$Stroke_Size <- 0.2
  
  if (!is.na(src_col)) {
    is_target <- grepl(highlight_pattern, plot_data[[src_col]])
    if (any(is_target)) {
      plot_data$Highlight_Type[is_target] <- "Target"
      plot_data$Stroke_Size[is_target] <- 1.5
    }
  }
  
  # Check colors
  c_mean <- if(!is.null(stat_colors["mean"])) stat_colors["mean"] else "darkred"
  c_med  <- if(!is.null(stat_colors["median"])) stat_colors["median"] else "black"
  
  p <- ggplot(plot_data, aes(x = Group, y = Value, fill = Group)) +
    geom_violin(alpha = 0.4, trim = FALSE, color = NA, scale = "width") +
    geom_jitter(aes(color = Highlight_Type, stroke = Stroke_Size), 
                width = 0.2, size = 2.5, shape = 21, alpha = 0.8) +
    stat_summary(fun = median, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.5, size = 0.8, color = c_med) +
    stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..),
                 width = 0.5, size = 0.8, color = c_mean, linetype = "dashed") +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = c("Standard" = "white", "Target" = highlight_color), guide = "none") +
    scale_continuous_identity(aesthetics = "stroke") +
    labs(
      title = paste("Raw Distribution:", marker_name),
      subtitle = sprintf("Stats: %s Line = Median | %s Dashed = Mean", "Solid", "Dashed"),
      y = "Raw Value", x = NULL
    ) +
    theme_coda() +
    theme(legend.position = "none", axis.text.x = element_text(face = "bold", size = 11))
  
  return(p)
}

#' @title Save Distribution Report (Multi-Page PDF)
viz_save_distribution_report <- function(data_df, markers, file_path, colors, hl_pattern = "", 
                                         stat_colors = c(mean="darkred", median="black")) {
  if (length(markers) == 0) return(NULL)
  
  message(sprintf("   [Viz] Saving distribution report to: %s", basename(file_path)))
  
  pdf(file_path, width = 8, height = 6)
  for (marker in markers) {
    if (marker %in% names(data_df)) {
      tryCatch({
        p <- plot_raw_distribution_merged(
          data_df = data_df, 
          marker_name = marker,
          colors = colors,     
          highlight_pattern = hl_pattern,   
          highlight_color = "#FFD700",
          stat_colors = stat_colors
        )
        if (!is.null(p)) print(p)
      }, error = function(e) warning(sprintf("Failed to plot marker %s: %s", marker, e$message)))
    }
  }
  dev.off()
}

#' @title Custom PCA Biplot
plot_pca_custom <- function(pca_res, metadata, colors, dims = c(1, 2), 
                            show_labels = FALSE, highlight_patterns = NULL,
                            find_col_keyword = "Original|Source") {
  # Setup Coordinates
  ind_coords <- as.data.frame(pca_res$ind$coord[, dims])
  colnames(ind_coords) <- c("X_Coord", "Y_Coord")
  
  if (nrow(ind_coords) != nrow(metadata)) stop("Metadata/PCA dimension mismatch.")
  plot_data <- cbind(metadata, ind_coords)
  
  plot_data$Highlight_Type <- "Standard"
  plot_data$Stroke_Size <- 0.5 
  
  target_col <- names(metadata)[grepl(find_col_keyword, names(metadata), ignore.case = TRUE)][1]
  has_highlights <- !is.null(highlight_patterns) && length(highlight_patterns) > 0 && !is.na(target_col)
  
  if (has_highlights) {
    target_values <- as.character(plot_data[[target_col]])
    for (pattern in names(highlight_patterns)) {
      matches <- grepl(pattern, target_values)
      to_update <- matches & (plot_data$Highlight_Type == "Standard")
      if (any(to_update)) {
        plot_data$Highlight_Type[to_update] <- pattern 
        plot_data$Stroke_Size[to_update] <- 2.0  
      }
    }
  }
  
  border_colors <- c("Standard" = "white")
  if (has_highlights) border_colors <- c(highlight_patterns, border_colors)
  
  eig_val <- get_eigenvalue(pca_res)
  var_x <- round(eig_val[dims[1], 2], 1)
  var_y <- round(eig_val[dims[2], 2], 1)
  
  p <- ggplot(plot_data, aes(x = X_Coord, y = Y_Coord, fill = Group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray80") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray80") +
    stat_ellipse(geom = "polygon", alpha = 0.1, show.legend = FALSE, level = 0.95) +
    geom_point(aes(color = Highlight_Type, stroke = Stroke_Size), size = 3.5, shape = 21) +
    scale_fill_manual(values = colors, guide = guide_legend(override.aes = list(shape = 21))) +
    scale_color_manual(name = "Subgroup", values = border_colors) +
    scale_continuous_identity(aesthetics = "stroke") +
    labs(x = sprintf("PC%d (%s%%)", dims[1], var_x),
         y = sprintf("PC%d (%s%%)", dims[2], var_y),
         title = sprintf("PCA (PC%d vs PC%d)", dims[1], dims[2])) +
    theme_coda()
  
  if (!has_highlights) p <- p + guides(color = "none")
  if (show_labels) p <- p + geom_text_repel(aes(label = Patient_ID), size = 3, show.legend = FALSE, max.overlaps = 15)
  return(p)
}

#' @title Plot Stratification Heatmap
#' @param gradient_colors Named vector with 'low', 'mid', 'high'.
plot_stratification_heatmap <- function(mat_z, metadata, annotation_colors_list, title = "Stratification",
                                        gradient_colors = c(low="#2166AC", mid="#F7F7F7", high="#B2182B")) {
  mat_plot <- t(mat_z)
  metadata_clean <- as.data.frame(metadata)
  
  ha <- HeatmapAnnotation(
    df = metadata_clean,
    col = annotation_colors_list,
    simple_anno_size = unit(0.3, "cm"),
    show_annotation_name = TRUE
  )
  
  # Use decoupled colors
  col_fun <- colorRamp2(c(-2, 0, 2), c(gradient_colors["low"], gradient_colors["mid"], gradient_colors["high"]))
  
  hm <- Heatmap(
    mat_plot,
    name = "Z-Score",
    col = col_fun,
    cluster_rows = TRUE, cluster_columns = TRUE,
    clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
    clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
    top_annotation = ha,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 6),
    show_column_names = FALSE, 
    column_title = paste0(title, " (n=", ncol(mat_plot), ")"),
    row_title = "Markers"
  )
  return(hm)
}

#' @title Generate Complex Heatmap Annotation Colors
viz_generate_complex_heatmap_colors <- function(metadata, group_col, base_colors, 
                                                hl_pattern = NULL, hl_color = "#FFD700") {
  annotation_colors <- list()
  present_groups <- unique(metadata[[group_col]])
  group_pal <- base_colors[intersect(names(base_colors), present_groups)]
  annotation_colors[[group_col]] <- group_pal
  
  extra_cols <- setdiff(colnames(metadata), group_col)
  for (col_name in extra_cols) {
    raw_vals <- as.character(metadata[[col_name]])
    vals <- sort(unique(raw_vals[!is.na(raw_vals)]))
    col_pal <- character()
    for (v in vals) {
      if (v %in% names(group_pal)) {
        col_pal[[v]] <- group_pal[[v]]
      } else if (!is.null(hl_pattern) && hl_pattern != "" && grepl(hl_pattern, v)) {
        col_pal[[v]] <- hl_color
      } else {
        parent_match <- NA
        for (grp in names(group_pal)) {
          if (grepl(grp, v)) {
            parent_match <- grp; break
          }
        }
        col_pal[[v]] <- if (!is.na(parent_match)) group_pal[[parent_match]] else "#D3D3D3"
      }
    }
    annotation_colors[[col_name]] <- col_pal
  }
  return(annotation_colors)
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

# (Other functions like plot_network_structure remain unchanged if not using colors)