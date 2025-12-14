# scripts/03_inference.R
# ==============================================================================
# STEP 03: ROBUST NETWORK INFERENCE (Bootstrap & Permutation)
# ==============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
  library(parallel)
  library(doParallel)
  library(igraph)
  library(conflicted)
  library(corpcor) # Explicitly load corpcor for the main thread
})

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("count", "dplyr")

# Load Config & Stats Utils
config <- read_yaml("config/global_params.yml")

# Ensure utils are loaded
if(file.exists("R/utils_stats.R")) {
  source("R/utils_stats.R")
} else {
  stop("Critical: R/utils_stats.R not found.")
}

# I/O Paths
input_rds <- file.path("data", "processed", "clean_data.rds")
output_dir <- file.path(config$output_root, "03_Networks")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

message("=== PIPELINE STEP 3: NETWORK INFERENCE ===")

# ------------------------------------------------------------------------------
# 1. Load Data & Prepare Matrices
# ------------------------------------------------------------------------------
if(!file.exists(input_rds)) stop("Run Step 01 first!")
DATA <- readRDS(input_rds)
df_clr <- DATA$clr
markers <- DATA$markers

# Define Groups
grp_ctrl <- config$control_group
grp_case <- config$case_groups

# Extract Matrices
get_mat <- function(g) df_clr %>% filter(Group %in% g) %>% select(all_of(markers)) %>% as.matrix()

mat_ctrl <- get_mat(grp_ctrl)
mat_case <- get_mat(grp_case)

message(sprintf("   -> Control Group (%s): %d samples", grp_ctrl, nrow(mat_ctrl)))
message(sprintf("   -> Case Group (%s): %d samples", paste(grp_case, collapse="+"), nrow(mat_case)))

# ------------------------------------------------------------------------------
# 2. Setup Parallel Cluster (Global)
# ------------------------------------------------------------------------------
n_cores <- if(config$stats$n_cores == "auto") detectCores() - 1 else config$stats$n_cores
message(sprintf("[System] Initializing Cluster with %d cores...", n_cores))

cl <- makeCluster(n_cores)
registerDoParallel(cl)
# Pre-load libraries/functions on workers
clusterEvalQ(cl, { 
  library(corpcor) 
})

# ------------------------------------------------------------------------------
# 3. Bootstrap Inference (Edge Stability)
# ------------------------------------------------------------------------------
message("\n[1] Running Bootstrap Inference (Individual Networks)...")

# Control Network
message(sprintf("   -> Bootstrapping %s (B=%d)...", grp_ctrl, config$stats$n_boot))
set.seed(config$stats$seed)
res_ctrl <- run_bootstrap_inference(mat_ctrl, B = config$stats$n_boot, cl = cl)

# Case Network
message(sprintf("   -> Bootstrapping Case Group (B=%d)...", config$stats$n_boot))
set.seed(config$stats$seed + 111)
res_case <- run_bootstrap_inference(mat_case, B = config$stats$n_boot, cl = cl)

# ------------------------------------------------------------------------------
# 4. Permutation Test (Differential Network)
# ------------------------------------------------------------------------------
message("\n[2] Running Permutation Test (Group Comparison)...")

# --- CRITICAL FIX START ---
# We must calculate the OBSERVED difference using single point estimates (raw data),
# NOT the bootstrap means. This ensures consistency with the null distribution
# generated in the permutation loop below (which uses single estimates).

message("   -> Calculating Observed Difference (Point Estimates)...")
pcor_obs_ctrl <- estimate_pcor_robust(mat_ctrl)
pcor_obs_case <- estimate_pcor_robust(mat_case)
obs_diff <- abs(pcor_obs_ctrl - pcor_obs_case)
# --- CRITICAL FIX END ---

# Prepare for Permutation
n1 <- nrow(mat_ctrl); n2 <- nrow(mat_case)
pool <- rbind(mat_ctrl, mat_case)

message(sprintf("   -> Permuting %d times...", config$stats$n_perm))
set.seed(config$stats$seed + 222)

# Export helper function to cluster manually if needed, 
# or rely on simple corpcor calls inside the loop to be safe.
clusterExport(cl, varlist = c("prec2part"), envir = environment())

null_diffs <- foreach(i = 1:config$stats$n_perm, .combine='c', .packages='corpcor') %dopar% {
  # Shuffle labels
  idx <- sample(1:(n1+n2))
  p1 <- pool[idx[1:n1], ]
  p2 <- pool[idx[(n1+1):(n1+n2)], ]
  
  # Fast estimation using shrinkage (Must match estimate_pcor_robust logic)
  # We inline the logic to avoid dependency issues on workers
  c1 <- corpcor::invcov.shrink(p1, verbose=FALSE)
  c2 <- corpcor::invcov.shrink(p2, verbose=FALSE)
  
  # Convert to Partial Correlation
  # Note: using prec2part from utils_stats (exported above) or inline math
  d1 <- diag(c1); r1 <- -as.matrix(c1)/sqrt(outer(d1,d1)); diag(r1)<-1
  d2 <- diag(c2); r2 <- -as.matrix(c2)/sqrt(outer(d2,d2)); diag(r2)<-1
  
  list(abs(r1 - r2))
}

# ------------------------------------------------------------------------------
# 5. Shutdown Cluster
# ------------------------------------------------------------------------------
stopCluster(cl)
message("[System] Cluster stopped.")

# ------------------------------------------------------------------------------
# 6. Calculate P-values & FDR
# ------------------------------------------------------------------------------
message("[3] Calculating Statistics (FDR on Union Network)...")

adj_ctrl <- res_ctrl$adj
adj_case <- res_case$adj
# Define Universe: Edges present in at least one condition (Union)
union_adj <- (adj_ctrl == 1 | adj_case == 1)

edges_to_test <- which(upper.tri(union_adj) & union_adj, arr.ind = TRUE)

message(sprintf("   -> Testing %d edges (Union of Bootstrap Networks)...", nrow(edges_to_test)))

if(nrow(edges_to_test) == 0) {
  warning("No stable edges found in either group. No statistics to compute.")
  edge_stats <- data.frame()
} else {
  
  edge_stats <- data.frame()
  
  for(k in 1:nrow(edges_to_test)) {
    i <- edges_to_test[k, 1]
    j <- edges_to_test[k, 2]
    
    # Observed Difference (from Point Estimates)
    obs_val <- obs_diff[i, j]
    
    # Null Distribution for this specific edge
    null_vals <- sapply(null_diffs, function(m) m[i, j])
    
    # P-value calculation
    # (Count how many null differences are >= observed difference)
    p_val <- (sum(null_vals >= obs_val) + 1) / (config$stats$n_perm + 1)
    
    # Store stats
    # Note: We store Bootstrap Means for visualization (Weight_Ctrl/Case),
    # but the P-value is derived from the Point Estimate comparison.
    edge_stats <- rbind(edge_stats, data.frame(
      Node1 = colnames(mat_ctrl)[i],
      Node2 = colnames(mat_ctrl)[j],
      Weight_Ctrl_BootMean = res_ctrl$pcor_mean[i, j],
      Weight_Case_BootMean = res_case$pcor_mean[i, j],
      Weight_Ctrl_Obs = pcor_obs_ctrl[i, j], # Added for clarity
      Weight_Case_Obs = pcor_obs_case[i, j], # Added for clarity
      Edge_In_Ctrl = adj_ctrl[i, j] == 1,
      Edge_In_Case = adj_case[i, j] == 1,
      Diff_Obs = obs_val,
      P_Value = p_val
    ))
  }
  
  # Apply FDR Correction (Benjamini-Hochberg)
  edge_stats$Q_Value <- p.adjust(edge_stats$P_Value, method = "BH")
  edge_stats$Significant <- edge_stats$Q_Value < config$stats$fdr_threshold
  
  # Save Results
  write.csv(edge_stats, file.path(output_dir, "Differential_Network_Stats.csv"), row.names = FALSE)
  
  n_sig <- sum(edge_stats$Significant)
  message(sprintf("   -> Found %d Significant Differential Edges (FDR < %.2f)", 
                  n_sig, config$stats$fdr_threshold))
}

# ------------------------------------------------------------------------------
# 7. Save Results
# ------------------------------------------------------------------------------
final_res <- list(
  ctrl = res_ctrl,
  case = res_case,
  obs_matrices = list(ctrl = pcor_obs_ctrl, case = pcor_obs_case),
  diff_stats = list(obs_diff=obs_diff, null_diffs=null_diffs)
)
saveRDS(final_res, file.path(output_dir, "inference_results.rds"))

message("=== PIPELINE STEP 3 COMPLETE ===")