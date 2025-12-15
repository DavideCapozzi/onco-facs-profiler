# src/03_inference.R
# ==============================================================================
# STEP 03: NETWORK INFERENCE (BOOTSTRAP & PERMUTATION)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(corpcor)
})

source("R/infrastructure.R")
source("R/modules_stats.R")

message("\n=== PIPELINE STEP 3: ROBUST NETWORK INFERENCE ===")

# 1. Load Config & Data
config <- load_config("config/global_params.yml")
input_file <- file.path(config$output_root, "01_processed", "data_processed.rds")

if (!file.exists(input_file)) stop("Step 01 output not found. Run src/01_ingest.R first.")

DATA <- readRDS(input_file)
df_clr <- DATA$clr_data 
markers <- DATA$markers

# 2. Prepare Data Subsets
grp_ctrl <- config$control_group
grp_case <- config$case_groups

get_group_matrix <- function(df, group_names, features) {
  mat <- df %>% 
    filter(Group %in% group_names) %>% 
    select(all_of(features)) %>% 
    as.matrix()
  return(mat)
}

mat_ctrl <- get_group_matrix(df_clr, grp_ctrl, markers)
mat_case <- get_group_matrix(df_clr, grp_case, markers)

message(sprintf("[Data] Control Group (%s): %d samples", grp_ctrl, nrow(mat_ctrl)))
message(sprintf("[Data] Case Group (%s): %d samples", paste(grp_case, collapse="+"), nrow(mat_case)))

# 3. Setup Parallel Cluster
n_cores_req <- config$stats$n_cores
if (n_cores_req == "auto") n_cores_req <- parallel::detectCores() - 1

message(sprintf("[System] Initializing Cluster with %d cores...", n_cores_req))
cl <- makeCluster(n_cores_req)
registerDoParallel(cl)

# --- FIX: EXPORT CUSTOM FUNCTIONS TO WORKERS ---
# Workers start with an empty environment. We must export our custom logic.
# We export the bootstrap wrapper AND the estimation functions it depends on.
clusterExport(cl, varlist = c("run_single_bootstrap", 
                              "estimate_pcor_robust", 
                              "prec2part"), 
              envir = .GlobalEnv)

# Also load necessary libraries on workers
clusterEvalQ(cl, {
  library(corpcor)
  library(dplyr)
})

# 4. Bootstrap Inference (Edge Stability)
# ------------------------------------------------------------------------------
run_parallel_bootstrap <- function(data_mat, n_boot, seed_val) {
  n_samples <- nrow(data_mat)
  
  # Set RNG Stream for reproducibility
  parallel::clusterSetRNGStream(cl, seed_val)
  
  # Run Loop
  res_list <- foreach(i = 1:n_boot, .packages = c("corpcor")) %dopar% {
    run_single_bootstrap(data_mat, n_samples)
  }
  
  return(res_list)
}

message("\n[Inference] 1. Running Bootstrap for Control Group...")
boot_ctrl_raw <- run_parallel_bootstrap(mat_ctrl, config$stats$n_boot, config$stats$seed)
res_ctrl <- aggregate_boot_results(boot_ctrl_raw, alpha = config$stats$alpha)

message("\n[Inference] 2. Running Bootstrap for Case Group...")
boot_case_raw <- run_parallel_bootstrap(mat_case, config$stats$n_boot, config$stats$seed + 1000)
res_case <- aggregate_boot_results(boot_case_raw, alpha = config$stats$alpha)

# --- CHANGED: Detailed Bootstrap Logging ---
n_possible <- (ncol(mat_ctrl) * (ncol(mat_ctrl) - 1)) / 2
n_edge_ctrl <- sum(res_ctrl$adj)/2
n_edge_case <- sum(res_case$adj)/2

message(sprintf("   -> Control Group: %d stable edges (Density: %.2f%%)", 
                n_edge_ctrl, (n_edge_ctrl/n_possible)*100))
message(sprintf("   -> Case Group:    %d stable edges (Density: %.2f%%)", 
                n_edge_case, (n_edge_case/n_possible)*100))

# 5. Differential Network Test (Permutation)
# ------------------------------------------------------------------------------
message("\n[Inference] 3. Running Permutation Test (Differential Edges)...")

n_perm <- config$stats$n_perm
n1 <- nrow(mat_ctrl)
n2 <- nrow(mat_case)
pool_mat <- rbind(mat_ctrl, mat_case)

obs_pcor_ctrl <- estimate_pcor_robust(mat_ctrl)
obs_pcor_case <- estimate_pcor_robust(mat_case)
obs_diff_mat  <- abs(obs_pcor_ctrl - obs_pcor_case)

# Set seed for permutation
clusterSetRNGStream(cl, config$stats$seed + 2000)

null_diffs <- foreach(i = 1:n_perm, .combine = 'c', .packages = "corpcor") %dopar% {
  # Shuffle
  shuffled_idx <- sample(1:(n1 + n2))
  p1 <- pool_mat[shuffled_idx[1:n1], ]
  p2 <- pool_mat[shuffled_idx[(n1 + 1):(n1 + n2)], ]
  
  # Estimate (Workers already have 'estimate_pcor_robust' from the export above)
  r1 <- estimate_pcor_robust(p1)
  r2 <- estimate_pcor_robust(p2)
  
  list(abs(r1 - r2))
}

max_obs_diff <- max(obs_diff_mat[upper.tri(obs_diff_mat)])
mean_obs_diff <- mean(obs_diff_mat[upper.tri(obs_diff_mat)])

message(sprintf("   -> Permutation Stats: %d iterations", n_perm))
message(sprintf("   -> Max Observed Diff: %.4f | Mean Observed Diff: %.4f", 
                max_obs_diff, mean_obs_diff))

# 6. Shutdown Cluster
stopCluster(cl)
message("[System] Parallel Cluster stopped.")

# 7. Calculate P-values and Save
# ------------------------------------------------------------------------------
message("[Post-Process] Calculating P-values and FDR...")

union_adj <- (res_ctrl$adj == 1 | res_case$adj == 1)
edges_indices <- which(upper.tri(union_adj) & union_adj, arr.ind = TRUE)

results_table <- data.frame()

if(nrow(edges_indices) > 0) {
  # Pre-calculate denominator for speed
  denom_p <- n_perm + 1
  
  for(k in 1:nrow(edges_indices)) {
    i <- edges_indices[k, 1]
    j <- edges_indices[k, 2]
    
    node_a <- rownames(union_adj)[i]
    node_b <- rownames(union_adj)[j]
    
    obs_val <- obs_diff_mat[i, j]
    null_vals <- sapply(null_diffs, function(m) m[i, j])
    
    p_val <- (sum(null_vals >= obs_val) + 1) / denom_p
    
    results_table <- rbind(results_table, data.frame(
      Node1 = node_a,
      Node2 = node_b,
      Weight_Ctrl = res_ctrl$weights[i, j],
      Weight_Case = res_case$weights[i, j],
      Is_Stable_Ctrl = res_ctrl$adj[i, j] == 1,
      Is_Stable_Case = res_case$adj[i, j] == 1,
      Diff_Score = obs_val,
      P_Value = p_val
    ))
  }
  
  results_table$FDR <- p.adjust(results_table$P_Value, method = "BH")
  results_table$Significant_Diff <- results_table$FDR < config$stats$fdr_threshold
  results_table <- results_table %>% arrange(P_Value)
}

# 8. Save Outputs
out_dir <- file.path(config$output_root, "03_networks")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write.csv(results_table, file.path(out_dir, "differential_stats.csv"), row.names = FALSE)

final_obj <- list(
  ctrl_network = res_ctrl,
  case_network = res_case,
  diff_table = results_table,
  config = config
)
saveRDS(final_obj, file.path(out_dir, "inference_results.rds"))

message(sprintf("[Output] Found %d significantly differential edges (FDR < %.2f)", 
                sum(results_table$Significant_Diff), config$stats$fdr_threshold))
message("=== STEP 3 COMPLETE ===\n")