# ============================================================================
# FACS Compositional Analysis - Modular Pipeline for Immunological Data
# Optimized for: NSCLC patients (n=10) vs Healthy donors (n=50)
# Best practices from: ALDEx2, CoDA, Flow Cytometry Guidelines
# ============================================================================

library(tidyverse)
library(readxl)
library(ALDEx2)
library(pheatmap)
library(RColorBrewer)
library(ggtern)      
library(patchwork)   
library(scales)

# ============================================================================
# STEP 1: IMPORT FROM EXCEL (Multi-sheet)
# ============================================================================

step1_import_excel <- function(file_path, 
                               tumor_sheet = "Tumor", 
                               healthy_sheet = "Healthy") {
  
  cat("\n===== STEP 1: Import from Excel =====\n\n")
  
  # 1. Read Tumor Data
  # ------------------
  tumor_data <- read_excel(file_path, sheet = tumor_sheet)
  tumor_data$Condition <- "Tumor"
  
  # Force conversion to numeric to catch non-numeric chars early
  tumor_data <- tumor_data %>%
    mutate(across(-any_of(c("Patient_ID", "Condition")),
                  ~ as.numeric(as.character(.)), 
                  .names = "{.col}"))
  
  # 2. Read Healthy Data
  # --------------------
  healthy_data <- read_excel(file_path, sheet = healthy_sheet)
  healthy_data$Condition <- "Healthy"
  
  healthy_data <- healthy_data %>%
    mutate(across(-any_of(c("Patient_ID", "Condition")), 
                  ~ as.numeric(as.character(.)), 
                  .names = "{.col}"))
  
  # 3. Combine and Format
  # ---------------------
  combined_data <- bind_rows(tumor_data, healthy_data)
  
  # Extract Metadata
  metadata <- combined_data %>% 
    select(Patient_ID, Condition) %>%
    mutate(Condition = factor(Condition, levels = c("Healthy", "Tumor"))) %>%
    as.data.frame() # Convert to base data.frame to safely handle rownames
  
  # Extract Abundance
  # We explicitly remove metadata columns to leave only measurements
  abundance <- combined_data %>% 
    select(-Patient_ID, -Condition) %>%
    as.data.frame() # CRITICAL: Convert tibble to data.frame for rownames support
  
  # Assign Row Names
  # This now works without warnings because we converted to data.frame above
  rownames(abundance) <- metadata$Patient_ID
  rownames(metadata) <- metadata$Patient_ID
  
  # 4. Diagnostics
  # --------------
  cat(sprintf("Total samples: %d\n", nrow(abundance)))
  cat(sprintf("Populations measured: %d\n", ncol(abundance)))
  
  # Check for non-numeric columns that might have slipped through
  non_numeric_cols <- names(abundance)[!sapply(abundance, is.numeric)]
  if(length(non_numeric_cols) > 0) {
    warning(paste("Non-numeric columns detected in abundance:", 
                  paste(non_numeric_cols, collapse=", ")))
    cat("  Attempting to remove non-numeric columns automatically...\n")
    abundance <- abundance[, sapply(abundance, is.numeric)]
  }
  
  cat(sprintf("Missing values: %.1f%%\n", 
              100 * sum(is.na(abundance)) / prod(dim(abundance))))
  
  cat("\n✓ Data import complete\n")
  
  return(list(abundance = abundance, metadata = metadata))
}

# ============================================================================
# STEP 2: EXPLORATORY DATA ANALYSIS
# ============================================================================

step2_explore_facs <- function(abundance, metadata) {
  
  cat("\n===== STEP 2: Exploratory Analysis =====\n\n")
  
  # 2.1: Missing data
  cat("2.1 Missing Data Pattern:\n")
  missing_by_pop <- colSums(is.na(abundance)) / nrow(abundance) * 100
  missing_by_sample <- rowSums(is.na(abundance)) / ncol(abundance) * 100
  
  high_miss_pops <- names(missing_by_pop[missing_by_pop > 50])
  if (length(high_miss_pops) > 0) {
    cat(sprintf("  Populations with >50%% missing (%d):\n", length(high_miss_pops)))
    for (pop in head(high_miss_pops, 10)) {
      cat(sprintf("    %s: %.1f%%\n", pop, missing_by_pop[pop]))
    }
  }
  
  # 2.2: Population statistics by condition
  cat("\n2.2 Population Abundance by Condition:\n")
  
  # Calculate median abundance per condition
  abundance_with_condition <- abundance %>%
    mutate(Condition = metadata$Condition)
  
  medians_by_condition <- abundance_with_condition %>%
    group_by(Condition) %>%
    summarise(across(where(is.numeric), ~median(.x, na.rm = TRUE)))
  
  cat("  Top 5 most abundant populations (overall median):\n")
  overall_medians <- apply(abundance, 2, median, na.rm = TRUE)
  top5 <- names(sort(overall_medians, decreasing = TRUE)[1:5])
  for (pop in top5) {
    cat(sprintf("    %s: %.2f%%\n", pop, overall_medians[pop]))
  }
  
  # 2.3: Key immunological ratios
  cat("\n2.3 Key Immunological Ratios:\n")
  
  if ("CD4" %in% names(abundance) && "CD8tot" %in% names(abundance)) {
    cd4_cd8_ratio <- abundance$CD4 / (abundance$CD8tot + 0.001)
    cat(sprintf("  CD4/CD8 ratio: %.2f ± %.2f\n", 
                mean(cd4_cd8_ratio, na.rm = TRUE),
                sd(cd4_cd8_ratio, na.rm = TRUE)))
  }
  
  if (all(c("NAIVE", "CM", "EM", "EFF") %in% names(abundance))) {
    memory_total <- rowSums(abundance[, c("CM", "EM", "EFF")], na.rm = TRUE)
    naive_memory_ratio <- abundance$NAIVE / (memory_total + 0.001)
    cat(sprintf("  Naive/Memory ratio: %.2f ± %.2f\n",
                mean(naive_memory_ratio, na.rm = TRUE),
                sd(naive_memory_ratio, na.rm = TRUE)))
  }
  
  eda_results <- list(
    missing_by_pop = missing_by_pop,
    missing_by_sample = missing_by_sample,
    medians_overall = overall_medians,
    medians_by_condition = medians_by_condition
  )
  
  cat("\n✓ Exploratory analysis complete\n")
  return(eda_results)
}

# ============================================================================
# STEP 3: QUALITY CONTROL WITH INTELLIGENT FILTERING
# ============================================================================

step3_qc_facs <- function(abundance, metadata,
                          pop_missing_threshold = 0.7,  # Keep if <70% missing
                          impute_method = "geometric_mean") {
  
  cat("\n===== STEP 3: Quality Control =====\n\n")
  
  initial_pops <- ncol(abundance)
  
  # Remove populations with excessive missingness
  missing_rate <- colSums(is.na(abundance)) / nrow(abundance)
  pops_to_keep <- names(missing_rate[missing_rate < pop_missing_threshold])
  pops_removed <- setdiff(names(abundance), pops_to_keep)
  
  if (length(pops_removed) > 0) {
    cat(sprintf("Removing %d populations with ≥%.0f%% missing:\n",
                length(pops_removed), pop_missing_threshold * 100))
    cat(sprintf("  %s\n", paste(head(pops_removed, 10), collapse = ", ")))
    if (length(pops_removed) > 10) cat("  ...\n")
  }
  
  abundance_filtered <- abundance[, pops_to_keep]
  
  # Impute remaining missing values
  if (any(is.na(abundance_filtered))) {
    cat(sprintf("\nImputing %d remaining NA values...\n", 
                sum(is.na(abundance_filtered))))
    
    for (col in names(abundance_filtered)) {
      na_idx <- is.na(abundance_filtered[[col]])
      if (any(na_idx)) {
        non_na_vals <- abundance_filtered[[col]][!na_idx]
        if (length(non_na_vals) > 0) {
          # Geometric mean for compositional data
          geom_mean <- exp(mean(log(non_na_vals + 0.001)))
          abundance_filtered[[col]][na_idx] <- geom_mean
        } else {
          # If all NA, use overall median across all populations
          abundance_filtered[[col]][na_idx] <- 0.1
        }
      }
    }
  }
  
  # Compositional closure (renormalize to 100%)
  cat("\nApplying compositional closure...\n")
  row_sums <- rowSums(abundance_filtered)
  abundance_closed <- abundance_filtered / row_sums * 100
  
  # Match metadata to remaining samples
  metadata_matched <- metadata[rownames(abundance_closed), ]
  
  cat("\n--- QC Summary ---\n")
  cat(sprintf("Populations: %d → %d\n", initial_pops, ncol(abundance_closed)))
  cat(sprintf("Samples: %d (no samples removed)\n", nrow(abundance_closed)))
  cat(sprintf("Final composition: Tumor=%d, Healthy=%d\n",
              sum(metadata_matched$Condition == "Tumor"),
              sum(metadata_matched$Condition == "Healthy")))
  
  cat("\n✓ Quality control complete\n")
  
  return(list(
    abundance = abundance_closed,
    metadata = metadata_matched,
    removed_pops = pops_removed
  ))
}

# ============================================================================
# STEP 4: CLR TRANSFORMATION
# ============================================================================

step4_clr_transform <- function(abundance, pseudocount = 0.01) {
  
  cat("\n===== STEP 4: CLR Transformation =====\n\n")
  
  # Add pseudocount
  abundance_pc <- abundance + pseudocount
  
  # CLR transformation
  clr_data <- t(apply(abundance_pc, 1, function(row) {
    gm <- exp(mean(log(row)))
    log(row / gm)
  }))
  
  colnames(clr_data) <- colnames(abundance)
  rownames(clr_data) <- rownames(abundance)
  
  cat(sprintf("CLR transformation complete\n"))
  cat(sprintf("  Pseudocount: %.3f\n", pseudocount))
  cat(sprintf("  Value range: [%.2f, %.2f]\n", 
              min(clr_data), max(clr_data)))
  
  cat("\n✓ CLR transformation complete\n")
  return(clr_data)
}

# ============================================================================
# STEP 5: PCA WITH IMPROVED VISUALIZATION
# ============================================================================

step5_pca_facs <- function(clr_data, metadata, n_top_loadings = 10) {
  
  cat("\n===== STEP 5: Principal Component Analysis =====\n\n")
  
  # Run PCA
  pca <- prcomp(clr_data, scale. = TRUE, center = TRUE)
  
  # Variance explained
  var_exp <- summary(pca)$importance[2, ]
  cum_var <- cumsum(var_exp)
  
  cat("Variance Explained:\n")
  cat(sprintf("  PC1: %.1f%%\n", var_exp[1] * 100))
  cat(sprintf("  PC2: %.1f%%\n", var_exp[2] * 100))
  cat(sprintf("  PC1+PC2: %.1f%%\n", cum_var[2] * 100))
  
  # Get top loadings for PC1 and PC2
  pc1_loadings <- abs(pca$rotation[, 1])
  pc2_loadings <- abs(pca$rotation[, 2])
  
  top_pc1 <- names(sort(pc1_loadings, decreasing = TRUE)[1:n_top_loadings])
  top_pc2 <- names(sort(pc2_loadings, decreasing = TRUE)[1:n_top_loadings])
  
  cat(sprintf("\nTop %d contributors to PC1:\n", n_top_loadings))
  for (pop in head(top_pc1, 5)) {
    cat(sprintf("  %s: %.3f\n", pop, pca$rotation[pop, 1]))
  }
  
  cat("\n✓ PCA complete\n")
  
  return(list(
    pca = pca,
    top_pc1 = top_pc1,
    top_pc2 = top_pc2
  ))
}

# ============================================================================
# STEP 6: DIFFERENTIAL ABUNDANCE (ALDEx2)
# ============================================================================

step6_differential_abundance <- function(abundance, metadata,
                                         mc_samples = 500,
                                         fdr_threshold = 0.1) {
  
  cat("\n===== STEP 6: Differential Abundance (ALDEx2) =====\n\n")
  
  # Ensure Metadata aligns with Abundance
  # -------------------------------------
  common_samples <- intersect(rownames(abundance), rownames(metadata))
  abundance <- abundance[common_samples, ]
  metadata <- metadata[common_samples, ]
  
  conditions <- metadata$Condition
  
  cat(sprintf("Comparison: %s (n=%d) vs %s (n=%d)\n",
              levels(conditions)[1], sum(conditions == levels(conditions)[1]),
              levels(conditions)[2], sum(conditions == levels(conditions)[2])))
  
  # -----------------------------------------------
  # ROBUST DATA PREPARATION
  # -----------------------------------------------
  
  # 1. Convert to Matrix
  # We use data.matrix() which handles factors better than as.matrix() 
  # just in case a factor slipped in.
  counts_matrix <- data.matrix(abundance)
  
  # 2. Handle NAs (CRITICAL FIX)
  # ALDEx2 fails with NAs. We must check if they exist.
  if (any(is.na(counts_matrix))) {
    cat("\nWARNING: NA values detected in data passed to ALDEx2.\n")
    
    # Option A: Simple Imputation (replace NA with small value/zero)
    # This prevents dropping samples, which is better for small cohorts.
    cat("  Imputing NAs with 0 for ALDEx2 analysis...\n")
    counts_matrix[is.na(counts_matrix)] <- 0
    
    # Alternative logic (commented out): Remove rows with NAs
    # na_rows <- which(!complete.cases(counts_matrix))
    # if(length(na_rows) > 0) {
    #    counts_matrix <- counts_matrix[-na_rows, ]
    #    conditions <- conditions[-na_rows]
    # }
  }
  
  # 3. Prepare for ALDEx2 (Integers)
  # ALDEx2 requires integer counts. Since FACS data is % or MFI,
  # we multiply by 100 (or 1000) and round to create pseudo-counts.
  # The error "round not meaningful for factors" is avoided because 
  # 'counts_matrix' is guaranteed to be a numeric matrix now.
  counts_integer <- round(counts_matrix * 100)
  
  # Ensure strictly integer mode
  mode(counts_integer) <- "integer"
  
  # 4. Filter Zero Variance Features
  # Remove columns that are all 0 across all samples
  keep_cols <- colSums(counts_integer) > 0
  if(sum(!keep_cols) > 0) {
    cat(sprintf("  Removing %d features with 0 total abundance.\n", sum(!keep_cols)))
    counts_integer <- counts_integer[, keep_cols]
  }
  
  # -----------------------------------------------
  # Run ALDEx2
  # -----------------------------------------------
  
  cat(sprintf("\nRunning ALDEx2 with %d Monte Carlo samples...\n", mc_samples))
  cat("(This analysis assumes inputs are compositional counts)\n")
  
  # Run the modular ALDEx2 steps manually for better error control
  # or use the wrapper 'aldex'. Here we use the wrapper.
  aldex_result <- tryCatch({
    ALDEx2::aldex(
      reads = counts_integer,
      conditions = conditions,
      test = "t",
      effect = TRUE,
      denom = "iqlr", # inter-quartile log-ratio is robust for asymmetrical data
      mc.samples = mc_samples,
      verbose = FALSE
    )
  }, error = function(e) {
    cat("\nERROR inside ALDEx2 execution:\n")
    print(e)
    return(NULL)
  })
  
  if(is.null(aldex_result)) return(NULL)
  
  # -----------------------------------------------
  # Format and Return Results
  # -----------------------------------------------
  
  # Convert results to a clean Tibble
  da_results <- aldex_result %>%
    tibble::rownames_to_column("Population") %>%
    arrange(we.eBH) %>% # Sort by Benjamini-Hochberg corrected p-value
    mutate(
      Significant = we.eBH < fdr_threshold,
      Direction = case_when(
        effect > 0 ~ paste0("UP in ", levels(conditions)[2]),
        effect < 0 ~ paste0("UP in ", levels(conditions)[1]),
        TRUE ~ "No change"
      ),
      EffectSize_Category = case_when(
        abs(effect) < 0.5 ~ "Small",
        abs(effect) < 1.0 ~ "Medium",
        TRUE ~ "Large"
      )
    )
  
  # Summary Output
  n_sig <- sum(da_results$Significant, na.rm = TRUE)
  cat(sprintf("\n--- ALDEx2 Results ---\n"))
  cat(sprintf("Significant features (FDR < %.2f): %d\n", fdr_threshold, n_sig))
  
  if (n_sig > 0) {
    print(head(filter(da_results, Significant), 10))
  } else {
    cat("No significant features found. Showing top effect sizes:\n")
    print(head(da_results %>% arrange(desc(abs(effect))), 5))
  }
  
  cat("\n✓ Differential abundance analysis complete\n")
  
  return(da_results)
}

# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================

# Categorize populations for color coding
categorize_populations <- function(pop_names) {
  
  categories <- case_when(
    pop_names %in% c("CD3", "CD4", "CD8tot") ~ "Major Lineages",
    pop_names %in% c("T.reg", "Resting", "Active", "Non-suppressive", "TregCD137") ~ "Treg States",
    pop_names %in% c("CM", "EFF", "EM", "NAIVE") ~ "Differentiation",
    grepl("^CD28", pop_names) ~ "CD28+ (Costimulation)",
    grepl("^CD137", pop_names) ~ "CD137+ (Activation)",
    grepl("^KI67|^KI657", pop_names) ~ "KI67+ (Proliferation)",
    grepl("^PD1", pop_names) ~ "PD1+ (Exhaustion)",
    grepl("M-MDSC|LOX1-PMN-MDSC", pop_names) ~ "MDSCs",
    TRUE ~ "Other"
  )
  
  return(categories)
}

# 1. Improved PCA Biplot
plot_pca_biplot_clean <- function(pca_result, metadata, 
                                  top_loadings = NULL,
                                  title = "PCA - Samples and Key Populations") {
  
  # Extract PC scores
  scores <- as.data.frame(pca_result$pca$x[, 1:2])
  scores$Condition <- metadata$Condition
  scores$Sample <- rownames(scores)
  
  # Variance explained
  var_exp <- summary(pca_result$pca)$importance[2, 1:2] * 100
  
  # Get top loadings if specified
  if (!is.null(top_loadings)) {
    loadings_to_plot <- unique(c(pca_result$top_pc1, pca_result$top_pc2))
    loadings_to_plot <- head(loadings_to_plot, top_loadings)
  } else {
    loadings_to_plot <- rownames(pca_result$pca$rotation)
  }
  
  loadings <- as.data.frame(pca_result$pca$rotation[loadings_to_plot, 1:2] * 5)
  loadings$Population <- rownames(loadings)
  loadings$Category <- categorize_populations(loadings$Population)
  
  # Plot
  p <- ggplot() +
    # Sample points
    geom_point(data = scores, 
               aes(x = PC1, y = PC2, color = Condition, shape = Condition),
               size = 4, alpha = 0.7) +
    # Loading arrows (only top contributors)
    geom_segment(data = loadings,
                 aes(x = 0, y = 0, xend = PC1, yend = PC2, color = Category),
                 arrow = arrow(length = unit(0.2, "cm")),
                 alpha = 0.6, linewidth = 0.5) +
    geom_text(data = loadings,
              aes(x = PC1, y = PC2, label = Population, color = Category),
              hjust = -0.1, vjust = -0.1, size = 2.5, fontface = "bold") +
    # Styling
    scale_color_manual(values = c("Healthy" = "#4DAF4A", "Tumor" = "#E41A1C",
                                  "Major Lineages" = "#377EB8",
                                  "Treg States" = "#FF7F00",
                                  "Differentiation" = "#984EA3",
                                  "CD28+ (Costimulation)" = "#A65628",
                                  "CD137+ (Activation)" = "#F781BF",
                                  "KI67+ (Proliferation)" = "#999999",
                                  "PD1+ (Exhaustion)" = "#000000",
                                  "MDSCs" = "#FFFF33",
                                  "Other" = "#CCCCCC")) +
    scale_shape_manual(values = c("Healthy" = 16, "Tumor" = 17)) +
    labs(x = sprintf("PC1 (%.1f%%)", var_exp[1]),
         y = sprintf("PC2 (%.1f%%)", var_exp[2]),
         title = title) +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(p)
}

# 2. Stratified boxplots for key populations
plot_stratified_boxplots <- function(abundance, metadata, 
                                    populations = NULL,
                                    ncol = 3) {
  
  if (is.null(populations)) {
    # Select top 9 most variable populations
    cv <- apply(abundance, 2, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
    populations <- names(sort(cv, decreasing = TRUE)[1:9])
  }
  
  # Reshape data
  plot_data <- abundance[, populations] %>%
    mutate(Condition = metadata$Condition) %>%
    pivot_longer(-Condition, names_to = "Population", values_to = "Abundance")
  
  # Plot
  p <- ggplot(plot_data, aes(x = Condition, y = Abundance, fill = Condition)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
    facet_wrap(~Population, scales = "free_y", ncol = ncol) +
    scale_fill_manual(values = c("Healthy" = "#4DAF4A", "Tumor" = "#E41A1C")) +
    labs(y = "Abundance (%)", title = "Key Populations by Condition") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "lightgray"),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(p)
}

# 3. Heatmap ordered by condition
plot_heatmap_ordered <- function(abundance, metadata, top_n = 30) {
  
  # Select most variable populations
  cv <- apply(abundance, 2, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
  top_pops <- names(sort(cv, decreasing = TRUE)[1:min(top_n, ncol(abundance))])
  
  abundance_subset <- abundance[, top_pops]
  
  # Order samples by condition
  sample_order <- order(metadata$Condition)
  abundance_ordered <- abundance_subset[sample_order, ]
  
  # Annotation
  annotation_row <- data.frame(
    Condition = metadata$Condition[sample_order]
  )
  rownames(annotation_row) <- rownames(abundance_ordered)
  
  ann_colors <- list(
    Condition = c(Healthy = "#4DAF4A", Tumor = "#E41A1C")
  )
  
  # Plot
  pheatmap(t(abundance_ordered),
           scale = "row",
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "euclidean",
           annotation_col = annotation_row,
           annotation_colors = ann_colors,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           main = "Population Abundance Heatmap (Top Variable)",
           fontsize_row = 8,
           fontsize_col = 8,
           show_colnames = TRUE)
}

# 4. Compositional stacked barplot
plot_composition_barplot <- function(abundance, metadata, top_n = 15) {
  
  # Select top populations by mean abundance
  means <- colMeans(abundance, na.rm = TRUE)
  top_pops <- names(sort(means, decreasing = TRUE)[1:top_n])
  
  # Create "Other" category
  abundance_plot <- abundance %>%
    mutate(Other = 100 - rowSums(select(., all_of(top_pops)))) %>%
    select(all_of(c(top_pops, "Other"))) %>%
    mutate(Sample = rownames(.), Condition = metadata$Condition) %>%
    pivot_longer(-c(Sample, Condition), names_to = "Population", values_to = "Abundance")
  
  # Order samples by condition
  sample_levels <- abundance_plot %>%
    arrange(Condition, Sample) %>%
    pull(Sample) %>%
    unique()
  
  abundance_plot$Sample <- factor(abundance_plot$Sample, levels = sample_levels)
  
  # Plot
  p <- ggplot(abundance_plot, aes(x = Sample, y = Abundance, fill = Population)) +
    geom_bar(stat = "identity", width = 0.9) +
    facet_grid(~Condition, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c(colorRampPalette(brewer.pal(12, "Set3"))(top_n), 
                                 "Other" = "gray80")) +
    labs(y = "Relative Abundance (%)", 
         title = sprintf("Compositional Profile (Top %d Populations)", top_n)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(p)
}

# 5. Volcano plot for differential abundance
plot_volcano <- function(da_results, fdr_threshold = 0.1, 
                        effect_threshold = 0.5) {
  
  da_results <- da_results %>%
    mutate(
      Significance = case_when(
        we.eBH < fdr_threshold & abs(effect) > effect_threshold ~ "Significant",
        we.eBH < fdr_threshold ~ "FDR < 0.1",
        abs(effect) > effect_threshold ~ "|Effect| > 0.5",
        TRUE ~ "Not significant"
      )
    )
  
  p <- ggplot(da_results, aes(x = effect, y = -log10(we.eBH))) +
    geom_point(aes(color = Significance), alpha = 0.6, size = 3) +
    geom_hline(yintercept = -log10(fdr_threshold), 
               linetype = "dashed", color = "red", alpha = 0.5) +
    geom_vline(xintercept = c(-effect_threshold, effect_threshold),
               linetype = "dashed", color = "blue", alpha = 0.5) +
    geom_text_repel(data = filter(da_results, Significance == "Significant"),
                    aes(label = Population), size = 3, max.overlaps = 20) +
    scale_color_manual(values = c("Significant" = "red",
                                  "FDR < 0.1" = "orange",
                                  "|Effect| > 0.5" = "steelblue",
                                  "Not significant" = "gray60")) +
    labs(x = "Effect Size", 
         y = "-log10(FDR-adjusted p-value)",
         title = "Differential Abundance - Volcano Plot") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "bottom")
  
  return(p)
}

# 6. MA plot from ALDEx2
plot_ma_aldex <- function(da_results) {
  
  p <- ggplot(da_results, aes(x = rab.all, y = diff.btw)) +
    geom_point(aes(color = Significant), alpha = 0.6, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray60")) +
    labs(x = "Log-ratio Abundance (CLR)",
         y = "Difference Between Conditions",
         title = "MA Plot - Abundance vs Difference") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(p)
}

# 7. Key immunological ratios comparison
plot_immunological_ratios <- function(abundance, metadata) {
  
  ratios <- data.frame(
    Sample = rownames(abundance),
    Condition = metadata$Condition
  )
  
  # CD4/CD8 ratio
  if ("CD4" %in% names(abundance) && "CD8tot" %in% names(abundance)) {
    ratios$CD4_CD8 <- abundance$CD4 / (abundance$CD8tot + 0.001)
  }
  
  # Naive/Memory ratio
  if (all(c("NAIVE", "CM", "EM", "EFF") %in% names(abundance))) {
    memory_total <- rowSums(abundance[, c("CM", "EM", "EFF")], na.rm = TRUE)
    ratios$Naive_Memory <- abundance$NAIVE / (memory_total + 0.001)
  }
  
  # Treg/CD4 ratio
  if ("T.reg" %in% names(abundance) && "CD4" %in% names(abundance)) {
    ratios$Treg_CD4 <- abundance$T.reg / (abundance$CD4 + 0.001)
  }
  
  # Reshape for plotting
  ratios_long <- ratios %>%
    pivot_longer(-c(Sample, Condition), names_to = "Ratio", values_to = "Value")
  
  p <- ggplot(ratios_long, aes(x = Condition, y = Value, fill = Condition)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
    facet_wrap(~Ratio, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c("Healthy" = "#4DAF4A", "Tumor" = "#E41A1C")) +
    labs(y = "Ratio", title = "Key Immunological Ratios") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(p)
}

# ============================================================================
# COMPLETE WORKFLOW
# ============================================================================

run_complete_workflow <- function(file_path,
                                 tumor_sheet = "Tumor",
                                 healthy_sheet = "Healthy",
                                 generate_plots = TRUE) {
  
  cat("\n╔══════════════════════════════════════════════════════╗\n")
  cat("║  FACS Compositional Analysis - Complete Workflow     ║\n")
  cat("╚══════════════════════════════════════════════════════╝\n")
  
  # Step 1: Import
  data <- step1_import_excel(file_path, tumor_sheet, healthy_sheet)
  
  # Step 2: Explore
  eda <- step2_explore_facs(data$abundance, data$metadata)
  
  # Step 3: QC
  qc <- step3_qc_facs(data$abundance, data$metadata)
  
  # Step 4: CLR
  clr_data <- step4_clr_transform(qc$abundance)
  
  # Step 5: PCA
  pca_result <- step5_pca_facs(clr_data, qc$metadata, n_top_loadings = 15)
  
  # Step 6: Differential Abundance
  da_results <- step6_differential_abundance(qc$abundance, qc$metadata,
                                            mc_samples = 500, fdr_threshold = 0.1)
  
  # Generate plots
  plots <- NULL
  if (generate_plots) {
    cat("\n===== Generating Visualizations =====\n\n")
    
    plots <- list(
      pca_clean = plot_pca_biplot_clean(pca_result, qc$metadata, top_loadings = 12),
      boxplots = plot_stratified_boxplots(qc$abundance, qc$metadata),
      composition = plot_composition_barplot(qc$abundance, qc$metadata, top_n = 12),
      volcano = plot_volcano(da_results, fdr_threshold = 0.1),
      ma_plot = plot_ma_aldex(da_results),
      ratios = plot_immunological_ratios(qc$abundance, qc$metadata)
    )
    
    # Heatmap (separate as it uses base graphics)
    cat("Generating heatmap...\n")
    plot_heatmap_ordered(qc$abundance, qc$metadata, top_n = 25)
    
    cat("✓ Visualizations complete\n")
  }
  
  cat("\n╔══════════════════════════════════════════════════════╗\n")
  cat("║              Analysis Complete ✓                     ║\n")
  cat("╚══════════════════════════════════════════════════════╝\n")
  
  return(list(
    abundance_clean = qc$abundance,
    metadata = qc$metadata,
    clr_data = clr_data,
    pca = pca_result,
    da_results = da_results,
    eda = eda,
    plots = plots
  ))
}

# ============================================================================
# USAGE EXAMPLE
# ============================================================================

# Run complete analysis
results <- run_complete_workflow(
  file_path = "/mnt/c/Users/Davide/Desktop/Davide/bandi/0dottorato/projects/long survivors/DB_pulito_anonimo.xlsx",
  tumor_sheet = "NSCLC OS > 18",     # or "NSCLC" or whatever your sheet is named
  healthy_sheet = "HD",  # or "Control" etc.
  generate_plots = TRUE
)

# Access results
abundance_clean <- results$abundance_clean
da_results <- results$da_results

# Display plots
print(results$plots$pca_clean)
print(results$plots$volcano)
print(results$plots$boxplots)
print(results$plots$composition)
print(results$plots$ma_plot)
print(results$plots$ratios)

# Save plots
ggsave("pca_clean.pdf", results$plots$pca_clean, width = 12, height = 8)
ggsave("volcano.pdf", results$plots$volcano, width = 10, height = 8)
ggsave("boxplots.pdf", results$plots$boxplots, width = 12, height = 10)

# Export results
write_csv(results$da_results, "differential_abundance_results.csv")
