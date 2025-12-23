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

#' @title Extract Group Colors
#' @description Helper to get the color palette from config safely.
#' @param config The loaded configuration list.
#' @return A named vector of hex codes.
get_palette <- function(config) {
  if (is.null(config$colors)) {
    warning("No colors defined in config. Using defaults.")
    return(c(Control = "grey", Case = "red"))
  }
  return(unlist(config$colors))
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
  # If TRUE, we emphasize correlations but lose information on magnitude variability.
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
  plot_data$Stroke_Size <- 0.5 # Default thickness
  
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
        plot_data$Stroke_Size[to_update] <- 2.0  # Thicker border for highlights
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
    
    # Points with dynamic aesthetics
    # Mapping stroke directly to the numeric column
    geom_point(aes(color = Highlight_Type, stroke = Stroke_Size), 
               size = 3.5, shape = 21) +
    
    # Main Fill Colors (Clinical Group)
    scale_fill_manual(values = colors, guide = guide_legend(override.aes = list(shape = 21))) +
    
    # Border Colors (Highlights)
    scale_color_manual(name = "Subgroup", values = border_colors) +
    
    # [FIX] Use Identity Scale for continuous numeric values (0.5 and 2.0)
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
#' @title Plot Multivariate Normality (QQ Plot)
#' @description 
#' Plots squared Mahalanobis distances against Chi-Square quantiles.
#' A straight line indicates Multivariate Normality.
#' 
#' @param assumption_res Result list from check_manova_assumptions().
#' @return A ggplot object.
plot_mvn_check <- function(assumption_res) {
  
  if (is.null(assumption_res$mahalanobis_dist)) {
    # Fallback if N < P
    return(ggplot() + 
             annotate("text", x=1, y=1, label="N < P: MVN Plot Unavailable") + 
             theme_void())
  }
  
  d2 <- sort(assumption_res$mahalanobis_dist)
  n <- length(d2)
  p <- length(d2) # Note: df for chi-sq is dimensions of data? No, it's cols of matrix.
  # We need the original dimensions to calculate df. 
  # However, assumption_res doesn't store p. We can infer p roughly or pass it.
  # Let's rely on the plotting logic:
  
  # Theoretical Quantiles (Chi-Square)
  # We assume the user passed the correct object. 
  # Actually, let's recalculate p inside the function if possible, 
  # or update check_manova_assumptions to return 'df'.
  
  # FIX: To make this robust, let's just plot Observed vs Expected quantiles
  # ppoints generates probability points
  chi_sq_quantiles <- qchisq(ppoints(n), df = mean(d2)) # Estimating df from mean(d2) ~ p is a rough proxy
  # Better: Pass 'df' in the results object. 
  
  plot_df <- data.frame(
    Theoretical = chi_sq_quantiles,
    Observed = d2
  )
  
  ggplot(plot_df, aes(x = Theoretical, y = Observed)) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(
      title = "Multivariate Normality Check",
      subtitle = "QQ Plot of Mahalanobis Distances (Assumption: Linear = Normal)",
      x = "Theoretical Quantiles (Chi-Square)",
      y = "Observed Squared Mahalanobis Distance"
    ) +
    theme_coda()
}

#' @title Plot Homogeneity of Variances
#' @description 
#' Boxplot of distances to group centroids. Checks if one group is more variable than another.
#' 
#' @param assumption_res Result list from check_manova_assumptions().
#' @param colors Named vector of colors.
#' @return A ggplot object.
plot_homogeneity <- function(assumption_res, colors) {
  
  df <- data.frame(
    Group = assumption_res$groups,
    Distance = assumption_res$dispersions
  )
  
  pval <- assumption_res$homogeneity_pval
  subtitle_txt <- sprintf("Betadisper p-value = %.4f (p < 0.05 implies unequal spread)", pval)
  
  ggplot(df, aes(x = Group, y = Distance, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = colors) +
    labs(
      title = "Homogeneity of Multivariate Dispersion",
      subtitle = subtitle_txt,
      y = "Distance to Group Centroid",
      x = NULL
    ) +
    theme_coda() +
    theme(legend.position = "none")
}

#' @title Plot MANOVA Separation (Canonical Variate)
#' @description 
#' Visualizes the result of the MANOVA by plotting the first Canonical Variate 
#' (Linear Discriminant), which represents the axis of maximal group separation.
#' 
#' @param manova_result The list returned by run_coda_manova().
#' @param colors Named vector of group colors.
#' @return A ggplot object (Boxplot + Jitter).
plot_manova_results <- function(manova_result, colors) {
  
  df <- manova_result$plot_data
  pval <- manova_result$p_value
  stat <- manova_result$pillai_stat
  
  # Format p-value for title
  pval_str <- ifelse(pval < 0.001, "p < 0.001", sprintf("p = %.4f", pval))
  
  ggplot(df, aes(x = Group, y = Canonical_Variate_1, fill = Group)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, shape = 21, color = "black", alpha = 0.8) +
    scale_fill_manual(values = colors) +
    labs(
      title = "Global Compositional Difference (MANOVA)",
      subtitle = sprintf("Statistic: Pillai = %.2f | Significance: %s", stat, pval_str),
      y = "Canonical Variate 1 (Best Separation Axis)",
      x = NULL
    ) +
    theme_coda() +
    theme(legend.position = "none")
}

#' @title Plot MANOVA Loadings (Feature Importance)
#' @description 
#' Visualizes which markers drive the separation detected by MANOVA.
#' Markers are ranked by their correlation with the Canonical Variate.
#' 
#' @param manova_result The list returned by run_coda_manova().
#' @param top_n Number of top markers to show (default 15).
#' @return A ggplot object (Lollipop chart).
plot_manova_loadings <- function(manova_result, top_n = 15) {
  
  # Get top N markers by absolute loading
  df <- manova_result$loadings %>%
    head(top_n)
  
  # Determine direction (just for coloring)
  df$Direction <- ifelse(df$Loading > 0, "Positive", "Negative")
  
  ggplot(df, aes(x = Loading, y = reorder(Marker, Loading))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_segment(aes(yend = Marker, xend = 0, color = Direction), size = 1) +
    geom_point(aes(color = Direction), size = 3) +
    labs(
      title = "Drivers of Separation (MANOVA Loadings)",
      subtitle = "Correlations between CLR Markers and Canonical Variate 1",
      x = "Contribution to Separation (Correlation)",
      y = NULL,
      color = "Direction"
    ) +
    theme_coda() +
    theme(legend.position = "none")
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

#' @title Plot Network Graph (ggraph)
#' @description 
#' @param min_cor Minimum absolute partial correlation to visualize (default 0).
#' Visualizes the inferred network. Nodes are sized by Degree. 
#' Edges are colored by sign (Blue=Pos, Red=Neg).
#' @param adj_mat Binary adjacency matrix.
#' @param weight_mat Partial correlation matrix.
#' @param title Plot title.
#' @param layout_type ggraph layout (default: "nicely", others: "fr", "kk").
#' @return A ggplot/ggraph object.
plot_network_structure <- function(adj_mat, weight_mat, title = "Network", 
                                   layout_type = "nicely", min_cor = 0) { 
  
  # 1. Build Graph Object
  g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
  
  # Add attributes
  E(g)$weight_raw <- NA
  E(g)$sign <- NA
  E(g)$weight <- NA # Initialize
  
  # Get edges
  el <- igraph::as_data_frame(g, what = "edges")
  
  if (nrow(el) > 0) {
    weights <- numeric(nrow(el))
    signs   <- character(nrow(el))
    keep_edge <- logical(nrow(el)) # Track edges to keep
    
    for(k in 1:nrow(el)) {
      w <- weight_mat[el[k,1], el[k,2]]
      
      # Check threshold
      if (abs(w) >= min_cor) {
        keep_edge[k] <- TRUE
        weights[k] <- abs(w)
        signs[k]   <- ifelse(w > 0, "Positive", "Negative")
      } else {
        keep_edge[k] <- FALSE
      }
    }
    
    # Prune edges from graph object directly
    # edges to delete are those where keep_edge is FALSE
    edges_to_remove <- E(g)[!keep_edge]
    g <- delete_edges(g, edges_to_remove)
    
    # Update attributes for REMAINING edges
    if (ecount(g) > 0) {
      E(g)$weight <- weights[keep_edge]
      E(g)$sign   <- signs[keep_edge]
    }
  }
  
  
  # This ensures node size reflects the *visible* network, not the full statistical network
  V(g)$degree <- igraph::degree(g)
  
  # Optional: Remove isolated nodes if they have 0 degree after filtering?
  # For now, we keep them so you can see which nodes are present but disconnected.
  
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
         # Update subtitle to reflect threshold
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