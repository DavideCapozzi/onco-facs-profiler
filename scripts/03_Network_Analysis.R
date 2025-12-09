# ==============================================================================
# SCRIPT 03: DIFFERENTIAL NETWORK ANALYSIS
# PURPOSE: Permutation test for Spearman correlations (Tumor vs Healthy)
# INPUT: clean_data.rds
# OUTPUT: Network visualizations, edges tables (CSV + XLSX)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(qgraph)
  library(foreach)
  library(doParallel)
  library(openxlsx)
  library(reshape2)
})

# 1. CONFIGURATION -------------------------------------------------------------
CONFIG <- list(
  input_file = "/home/davidec/projects/compositional_analysis/processed_data/clean_data.rds",
  output_dir = "/home/davidec/projects/compositional_analysis/results/03_Network",
  export_dir = "/mnt/c/Users/Davide/Desktop/Davide/bandi/0dottorato/projects/long survivors/networks",
  
  # Define Comparison Groups
  group_healthy = "Healthy",
  
  # Network Parameters
  n_permutations = 5000,
  fdr_threshold = 0.05,
  min_delta_rho = 0.20,
  min_rho_single = 0.30,
  n_cores = parallel::detectCores() - 1
)

if(!dir.exists(CONFIG$output_dir)) dir.create(CONFIG$output_dir, recursive = TRUE)
if(!dir.exists(CONFIG$export_dir)) dir.create(CONFIG$export_dir, recursive = TRUE)

cat("=== PIPELINE STEP 3: NETWORK ANALYSIS ===\n")

# 2. DATA PREPARATION ----------------------------------------------------------
DATA <- readRDS(CONFIG$input_file)
df_clr <- DATA$clr_transformed
markers <- DATA$markers

# Pan-Cancer Grouping: NOT Healthy = Tumor (NSCLC + HNSCC)
df_analysis <- df_clr %>%
  dplyr::mutate(Network_Group = ifelse(Group == CONFIG$group_healthy, "Healthy", "Tumor"))

tumor_data <- df_analysis %>% 
  dplyr::filter(Network_Group == "Tumor") %>% 
  dplyr::select(all_of(markers))

healthy_data <- df_analysis %>% 
  dplyr::filter(Network_Group == "Healthy") %>% 
  dplyr::select(all_of(markers))

cat(sprintf("   -> [INFO] Grouping Strategy: Pan-Cancer Analysis\n"))
cat(sprintf("   -> Tumor (NSCLC + HNSCC): %d samples\n", nrow(tumor_data)))
cat(sprintf("   -> Healthy: %d samples\n", nrow(healthy_data)))

if(nrow(tumor_data) < 10) {
  cat("   [WARNING] Sample size for Tumor is very low (<10). Results may be unstable.\n")
}

# 3. PERMUTATION TEST ----------------------------------------------------------
run_permutation_test <- function(d1, d2, n_perms, cores) {
  
  # Observed correlations
  c1 <- cor(d1, method = "spearman")
  c2 <- cor(d2, method = "spearman")
  diff_obs <- c1 - c2
  
  # Prepare pooled data
  combined <- rbind(d1, d2)
  n1 <- nrow(d1)
  n_tot <- nrow(combined)
  
  cat(sprintf("   -> Running %d permutations on %d cores...\n", n_perms, cores))
  
  # Parallel permutation loop
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  perm_results <- foreach(i = 1:n_perms, .combine = 'c') %dopar% {
    perm_idx <- sample(1:n_tot)
    d1_p <- combined[perm_idx[1:n1], ]
    d2_p <- combined[perm_idx[(n1 + 1):n_tot], ]
    list(cor(d1_p, method = "spearman") - cor(d2_p, method = "spearman"))
  }
  
  stopCluster(cl)
  
  # Compute p-values
  p_mat <- matrix(0, nrow = nrow(diff_obs), ncol = ncol(diff_obs))
  for(r in 1:nrow(diff_obs)) {
    for(c in 1:ncol(diff_obs)) {
      obs_val <- abs(diff_obs[r, c])
      greater <- sum(sapply(perm_results, function(x) abs(x[r, c]) >= obs_val))
      p_mat[r, c] <- (greater + 1) / (n_perms + 1)
    }
  }
  
  colnames(p_mat) <- colnames(diff_obs)
  rownames(p_mat) <- rownames(diff_obs)
  
  return(list(
    diff = diff_obs, 
    p_val = p_mat, 
    r_tumor = c1, 
    r_healthy = c2
  ))
}

results <- run_permutation_test(tumor_data, healthy_data, 
                                CONFIG$n_permutations, CONFIG$n_cores)

# 4. UNIFIED EDGE TABLE CREATION -----------------------------------------------
# This function ensures consistent FDR calculation across all exports

create_edge_table <- function(corr_mat, p_mat = NULL, net_type, fdr_threshold = NULL) {
  
  # Step 1: Melt matrices
  edges <- reshape2::melt(corr_mat)
  colnames(edges) <- c("Source", "Target", "Weight")
  
  # Step 2: Add p-values if provided (before filtering to maintain alignment)
  if(!is.null(p_mat)) {
    p_vals <- reshape2::melt(p_mat)
    
    # Safety check: dimensions must match
    if(nrow(edges) != nrow(p_vals)) {
      stop("Error: Correlation and P-value matrices have different dimensions!")
    }
    
    edges$P_Value <- p_vals$value
  } else {
    edges$P_Value <- NA
  }
  
  # Step 3: Convert to character and filter for unique pairs
  edges <- edges %>%
    mutate(
      Source = as.character(Source),
      Target = as.character(Target),
      Network = net_type
    ) %>%
    filter(Source < Target)  # Upper triangle only, removes self-loops
  
  # Step 4: Calculate FDR on the filtered unique pairs (CRITICAL FIX)
  if(!is.null(p_mat)) {
    edges <- edges %>%
      mutate(
        P_FDR = p.adjust(P_Value, method = "fdr"),
        Interaction_Class = case_when(
          P_FDR < fdr_threshold & Weight > 0 ~ "Gain_in_Tumor",
          P_FDR < fdr_threshold & Weight < 0 ~ "Gain_in_Healthy",
          TRUE ~ "Non_Significant"
        )
      )
  } else {
    # For single networks without p-values
    edges <- edges %>%
      mutate(
        P_FDR = NA,
        Interaction_Class = ifelse(Weight > 0, "Positive_Correlation", "Negative_Correlation")
      )
  }
  
  return(edges)
}

# 5. EXPORT DIFFERENTIAL NETWORK -----------------------------------------------
cat("[4] Exporting Differential Network Results...\n")

# Create edge table with unified FDR calculation
edges_diff_all <- create_edge_table(
  corr_mat = results$diff,
  p_mat = results$p_val,
  net_type = "Differential_Tumor_Minus_Healthy",
  fdr_threshold = CONFIG$fdr_threshold
)

# Filter significant edges only
edges_diff_sig <- edges_diff_all %>%
  dplyr::filter(P_FDR < CONFIG$fdr_threshold & abs(Weight) > CONFIG$min_delta_rho) %>%
  dplyr::arrange(P_FDR) %>%
  dplyr::select(Source, Target, Delta_Rho = Weight, P_Value, P_FDR, Category = Interaction_Class)

# Export significant edges
write.xlsx(edges_diff_sig, 
           file.path(CONFIG$output_dir, "Significant_Edges.xlsx"))

write.csv(edges_diff_sig, 
          file.path(CONFIG$output_dir, "Significant_Edges.csv"), 
          row.names = FALSE, quote = FALSE)

cat(sprintf("   -> [INFO] Found %d significant edges (FDR < %.2f, |Delta_Rho| > %.2f)\n", 
            nrow(edges_diff_sig), CONFIG$fdr_threshold, CONFIG$min_delta_rho))

# Export full differential network (for Cytoscape/ggraph)
edges_diff_export <- edges_diff_all %>%
  dplyr::select(Source, Target, Weight, P_Value, Network, P_FDR, Interaction_Class)

write.csv(edges_diff_export, 
          file.path(CONFIG$export_dir, "Network_Differential.csv"), 
          row.names = FALSE, quote = FALSE)

write.xlsx(edges_diff_export, 
           file.path(CONFIG$export_dir, "Network_Differential.xlsx"))

cat(sprintf("   -> Exported full differential network (%d edges) to export directory\n", 
            nrow(edges_diff_export)))

# 6. EXPORT SINGLE NETWORKS ----------------------------------------------------
cat("[5] Exporting Individual Network Results...\n")

# Tumor Network
edges_tumor <- create_edge_table(
  corr_mat = results$r_tumor,
  p_mat = NULL,
  net_type = "Tumor_Only"
) %>%
  dplyr::select(Source, Target, Weight, Network, Interaction_Class)

write.csv(edges_tumor, 
          file.path(CONFIG$export_dir, "Network_Tumor.csv"), 
          row.names = FALSE, quote = FALSE)

write.xlsx(edges_tumor, 
           file.path(CONFIG$export_dir, "Network_Tumor.xlsx"))

# Healthy Network
edges_healthy <- create_edge_table(
  corr_mat = results$r_healthy,
  p_mat = NULL,
  net_type = "Healthy_Only"
) %>%
  dplyr::select(Source, Target, Weight, Network, Interaction_Class)

write.csv(edges_healthy, 
          file.path(CONFIG$export_dir, "Network_Healthy.csv"), 
          row.names = FALSE, quote = FALSE)

write.xlsx(edges_healthy, 
           file.path(CONFIG$export_dir, "Network_Healthy.xlsx"))

cat(sprintf("   -> Exported Tumor network (%d edges)\n", nrow(edges_tumor)))
cat(sprintf("   -> Exported Healthy network (%d edges)\n", nrow(edges_healthy)))

# 7. EXPORT NODE TABLE ---------------------------------------------------------
# Master list of all unique nodes
all_nodes <- unique(c(as.character(edges_diff_export$Source), 
                      as.character(edges_diff_export$Target)))

node_table <- data.frame(
  ID = all_nodes,
  Label = all_nodes,
  stringsAsFactors = FALSE
)

write.csv(node_table, 
          file.path(CONFIG$export_dir, "Node_Table.csv"), 
          row.names = FALSE, quote = FALSE)

write.xlsx(node_table, 
           file.path(CONFIG$export_dir, "Node_Table.xlsx"))

cat(sprintf("   -> Exported node table (%d nodes)\n", nrow(node_table)))

# 8. NETWORK VISUALIZATION -----------------------------------------------------
cat("[6] Generating Network Visualizations...\n")

# Calculate layout based on Healthy network (most stable reference)
L <- qgraph(results$r_healthy, DoNotPlot = TRUE)$layout

pdf(file.path(CONFIG$output_dir, "Differential_Network_PanCancer.pdf"), 
    width = 14, height = 8)

# Page 1: Differential Network
qgraph(results$diff, 
       layout = L, 
       graph = "default", 
       diag = FALSE,
       title = "Differential Network (Tumor [NSCLC+HNSCC] - Healthy)\nRed = Stronger in Tumor | Blue = Stronger in Healthy",
       minimum = CONFIG$min_delta_rho, 
       cut = CONFIG$min_delta_rho,
       posCol = "darkred", 
       negCol = "darkblue", 
       vsize = 4, 
       labels = colnames(results$diff), 
       label.cex = 0.9)

# Page 2: Side-by-side Single Networks
par(mfrow = c(1, 2))

# Tumor Network
qgraph(results$r_tumor, 
       layout = L, 
       graph = "cor", 
       diag = FALSE,
       title = "Tumor Network (Spearman)",
       minimum = CONFIG$min_rho_single, 
       cut = CONFIG$min_rho_single,
       posCol = "darkred", 
       negCol = "darkblue", 
       vsize = 4, 
       labels = colnames(results$diff), 
       label.cex = 0.8)

# Healthy Network
qgraph(results$r_healthy, 
       layout = L, 
       graph = "cor", 
       diag = FALSE,
       title = "Healthy Network (Spearman)",
       minimum = CONFIG$min_rho_single, 
       cut = CONFIG$min_rho_single,
       posCol = "darkred", 
       negCol = "darkblue", 
       vsize = 4, 
       labels = colnames(results$diff), 
       label.cex = 0.8)

dev.off()

cat(sprintf("   -> Network plots saved to: %s\n", CONFIG$output_dir))

# 9. VALIDATION SUMMARY --------------------------------------------------------
cat("\n=== EXPORT VALIDATION ===\n")
cat(sprintf("Significant edges (FDR < %.2f, |Delta| > %.2f): %d\n", 
            CONFIG$fdr_threshold, CONFIG$min_delta_rho, nrow(edges_diff_sig)))

cat(sprintf("Full differential network edges: %d\n", nrow(edges_diff_export)))
cat(sprintf("Tumor network edges: %d\n", nrow(edges_tumor)))
cat(sprintf("Healthy network edges: %d\n", nrow(edges_healthy)))

# Verify consistency
n_sig_in_full <- sum(edges_diff_export$Interaction_Class != "Non_Significant")
if(n_sig_in_full == nrow(edges_diff_sig)) {
  cat(sprintf("   -> [PASS] Consistency check: %d significant edges match across exports\n", 
              n_sig_in_full))
} else {
  cat(sprintf("   -> [WARNING] Inconsistency detected: %d vs %d significant edges\n", 
              n_sig_in_full, nrow(edges_diff_sig)))
}

cat(sprintf("\n[SUCCESS] All results saved to:\n"))
cat(sprintf("   - Main output: %s\n", CONFIG$output_dir))
cat(sprintf("   - Export files: %s\n", CONFIG$export_dir))
cat("=== PIPELINE STEP 3 COMPLETED ===\n")

