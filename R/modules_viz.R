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