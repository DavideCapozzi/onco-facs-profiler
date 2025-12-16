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

#' @title Custom PCA Biplot
#' @description Plots Individuals (dots) and Variables (arrows) in the same space.
#' @param pca_res Result from run_coda_pca().
#' @param metadata Dataframe containing 'Patient_ID' and 'Group'.
#' @param colors Named vector of colors for groups.
#' @param show_labels Logical. If TRUE, adds Patient_ID labels.
#' @return A ggplot object.
plot_pca_custom <- function(pca_res, metadata, colors, show_labels = FALSE) {
  
  # 1. Extract Individual Coordinates
  ind_coords <- as.data.frame(pca_res$ind$coord)
  colnames(ind_coords)[1:2] <- c("Dim1", "Dim2")
  
  # Bind metadata (ensure order matches)
  if (nrow(ind_coords) != nrow(metadata)) stop("Metadata/PCA dimension mismatch.")
  plot_data <- cbind(metadata, ind_coords)
  
  # 2. Extract Variance Explained
  eig_val <- get_eigenvalue(pca_res)
  var_d1 <- round(eig_val[1, 2], 1)
  var_d2 <- round(eig_val[2, 2], 1)
  
  # 3. Base Plot (Individuals)
  p <- ggplot(plot_data, aes(x = Dim1, y = Dim2, color = Group, fill = Group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray80") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray80") +
    # Ellipses for groups
    stat_ellipse(geom = "polygon", alpha = 0.1, show.legend = FALSE, level = 0.95) +
    # Points
    geom_point(size = 3, alpha = 0.8, shape = 21, color = "white", stroke = 0.5) +
    aes(fill = Group) + 
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(
      x = sprintf("PC1 (%s%%)", var_d1),
      y = sprintf("PC2 (%s%%)", var_d2),
      title = "PCA: Immunological Landscape (CLR Space)",
      subtitle = "Euclidean distance on CLR = Aitchison distance on Raw"
    ) +
    theme_coda()
  
  # 4. Optional Labels
  if (show_labels) {
    p <- p + geom_text_repel(aes(label = Patient_ID), size = 3, show.legend = FALSE)
  }
  
  return(p)
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