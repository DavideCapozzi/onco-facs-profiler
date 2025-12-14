# R/utils_viz.R
# ==============================================================================
# VISUALIZATION UTILITIES
# Description: Custom themes and color helpers based on config.
# ==============================================================================

library(ggplot2)
library(yaml)

#' Get Group Colors
#' @param config The loaded configuration list
#' @return Named vector of hex codes
get_group_colors <- function(config) {
  cols <- unlist(config$colors)
  return(cols)
}

#' Standard Project Theme
#' @return ggplot theme object
theme_coda <- function() {
  theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

#' PCA Biplot Helper
#' @param pca_res Result from FactoMineR::PCA
#' @param metadata Dataframe with Group and Patient_ID
#' @param colors Named vector of colors
#' @return ggplot object
plot_pca_custom <- function(pca_res, metadata, colors) {
  # Extract coordinates
  coords <- data.frame(
    Dim1 = pca_res$ind$coord[,1],
    Dim2 = pca_res$ind$coord[,2],
    Group = metadata$Group,
    ID = metadata$Patient_ID
  )
  
  # Calculate variance explained
  var_exp <- pca_res$eig[,2]
  
  ggplot(coords, aes(x = Dim1, y = Dim2, color = Group, fill = Group)) +
    geom_hline(yintercept = 0, linetype="dashed", color="gray") +
    geom_vline(xintercept = 0, linetype="dashed", color="gray") +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(geom = "polygon", alpha = 0.1, show.legend = FALSE) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    labs(
      x = sprintf("Dim 1 (%.1f%%)", var_exp[1]),
      y = sprintf("Dim 2 (%.1f%%)", var_exp[2]),
      title = "PCA: Immunological Landscape (CLR)"
    ) +
    theme_coda()
}
