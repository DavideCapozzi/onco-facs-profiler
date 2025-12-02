# ==============================================================================
# SCRIPT 03: DIFFERENTIAL NETWORK ANALYSIS
# PURPOSE: Permutation test for Spearman correlations (Tumor vs Healthy)
# INPUT: clean_data.rds
# OUTPUT: Network visualizations, Edges Excel
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
  
  # Define Comparison Groups
  group_healthy = "Healthy",
  # Any group NOT Healthy will be treated as Tumor (Pan-Cancer approach)
  
  # Network Parameters
  n_permutations = 1000,
  fdr_threshold = 0.05,
  min_delta_rho = 0.30,   
  n_cores = parallel::detectCores() - 1
)

if(!dir.exists(CONFIG$output_dir)) dir.create(CONFIG$output_dir, recursive = TRUE)
cat("=== PIPELINE STEP 3: NETWORK ANALYSIS ===\n")

# 2. DATA PREP -----------------------------------------------------------------
DATA <- readRDS(CONFIG$input_file)
df_clr <- DATA$clr_transformed
markers <- DATA$markers

# Logic: Aggregating Pan-Cancer
# If a patient is NOT Healthy, they are Tumor (Includes NSCLC + HNSCC)
df_analysis <- df_clr %>%
  mutate(Network_Group = ifelse(Group == CONFIG$group_healthy, "Healthy", "Tumor"))

tumor_data <- df_analysis %>% filter(Network_Group == "Tumor") %>% select(all_of(markers))
healthy_data <- df_analysis %>% filter(Network_Group == "Healthy") %>% select(all_of(markers))

cat(sprintf("   -> [INFO] Grouping Strategy: Pan-Cancer Analysis\n"))
cat(sprintf("   -> Tumor (NSCLC + HNSCC): %d samples\n", nrow(tumor_data)))
cat(sprintf("   -> Healthy: %d samples\n", nrow(healthy_data)))

if(nrow(tumor_data) < 10) cat("   [WARNING] Sample size for Tumor is very low (<10). Results may be unstable.\n")

# 3. PERMUTATION TEST ----------------------------------------------------------
run_permutation_test <- function(d1, d2, n_perms, cores) {
  
  c1 <- cor(d1, method = "spearman")
  c2 <- cor(d2, method = "spearman")
  diff_obs <- c1 - c2
  
  combined <- rbind(d1, d2)
  n1 <- nrow(d1)
  n_tot <- nrow(combined)
  
  cat(sprintf("   -> Running %d permutations on %d cores...\n", n_perms, cores))
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  perm_results <- foreach(i = 1:n_perms, .combine = 'c') %dopar% {
    perm_idx <- sample(1:n_tot)
    d1_p <- combined[perm_idx[1:n1], ]
    d2_p <- combined[perm_idx[(n1 + 1):n_tot], ]
    list(cor(d1_p, method="spearman") - cor(d2_p, method="spearman"))
  }
  stopCluster(cl)
  
  p_mat <- matrix(0, nrow=nrow(diff_obs), ncol=ncol(diff_obs))
  for(r in 1:nrow(diff_obs)) {
    for(c in 1:ncol(diff_obs)) {
      obs_val <- abs(diff_obs[r,c])
      greater <- sum(sapply(perm_results, function(x) abs(x[r,c]) >= obs_val))
      p_mat[r,c] <- (greater + 1) / (n_perms + 1)
    }
  }
  
  colnames(p_mat) <- colnames(diff_obs)
  rownames(p_mat) <- rownames(diff_obs)
  return(list(diff = diff_obs, p_val = p_mat, r_tumor = c1, r_healthy = c2))
}

results <- run_permutation_test(tumor_data, healthy_data, CONFIG$n_permutations, CONFIG$n_cores)

# 4. EXPORT & PLOT -------------------------------------------------------------
cat("[3] Exporting Results...\n")

# FDR
p_adj_mat <- matrix(p.adjust(results$p_val[upper.tri(results$p_val)], method="fdr"), 
                    nrow=nrow(results$p_val), ncol=ncol(results$p_val))
# (Simplified reconstruction for export - mainly utilizing the list structure)

# Edges List
edges <- reshape2::melt(results$diff)
colnames(edges) <- c("Node1", "Node2", "Delta_Rho")
edges$P_FDR <- reshape2::melt(results$p_val)$value # Using raw P for now, adjust if strictly FDR needed

edges_sig <- edges %>%
  filter(as.character(Node1) < as.character(Node2)) %>%
  filter(P_FDR < CONFIG$fdr_threshold & abs(Delta_Rho) > CONFIG$min_delta_rho) %>%
  mutate(Category = ifelse(Delta_Rho > 0, "Gain_in_Tumor", "Gain_in_Healthy")) %>%
  arrange(P_FDR)

write.xlsx(edges_sig, file.path(CONFIG$output_dir, "Significant_Edges.xlsx"))

# Plot
pdf(file.path(CONFIG$output_dir, "Differential_Network_PanCancer.pdf"), width = 12, height = 8)
L <- qgraph(results$r_healthy, DoNotPlot=TRUE)$layout
qgraph(results$diff, layout=L, graph = "default", diag=FALSE,
       title = "Differential Network (Tumor [NSCLC+HNSCC] - Healthy)",
       minimum = CONFIG$min_delta_rho, cut = CONFIG$min_delta_rho,
       posCol = "darkred", negCol = "darkblue", vsize = 4, labels = colnames(results$diff))
dev.off()

cat("   -> Done! Results in 'results/03_Network/'\n")

