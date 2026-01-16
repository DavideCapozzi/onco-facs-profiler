# src/04_network_inference.R
# ==============================================================================
# STEP 04: NETWORK INFERENCE (BOOTSTRAP & PERMUTATION)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(corpcor)
  library(openxlsx) 
})

source("R/utils_io.R")          
source("R/modules_network.R")   

message("\n=== PIPELINE STEP 4: ROBUST NETWORK INFERENCE ===")

# 1. Load Config & Data
config <- load_config("config/global_params.yml")
input_file <- file.path(config$output_root, "01_data_processing", "data_processed.rds")
if (!file.exists(input_file)) stop("Step 01 output not found.")

DATA <- readRDS(input_file)
df_input <- DATA$hybrid_data_z 
markers <- DATA$hybrid_markers

# 2. Dynamic Stratification
# ------------------------------------------------------------------------------
strat_mode <- config$stratification$mode
message(sprintf("[Setup] Stratification Mode: '%s'", strat_mode))

if (strat_mode == "stratified_ls") {
  strat_conf <- config$stratification
  df_input <- df_input %>%
    mutate(Group = case_when(
      Group %in% config$control_group ~ Group, 
      grepl(strat_conf$pattern, .data[[strat_conf$column]]) ~ paste0(Group, strat_conf$suffix_match),
      TRUE ~ paste0(Group, strat_conf$suffix_no_match)
    ))
} else if (strat_mode == "binary") {
  df_input <- df_input %>%
    mutate(Group = case_when(
      Group %in% config$control_group ~ "Control",
      TRUE ~ "Case"
    ))
}

# 3. Target Selection
# ------------------------------------------------------------------------------
grp_ctrl_name <- config$control_group
grp_case_names <- config$case_groups

message(sprintf("[Config] Network Contrast: %s (Control) vs %s (Case)", 
                grp_ctrl_name, paste(grp_case_names, collapse="+")))

# Validation
if (!grp_ctrl_name %in% df_input$Group) stop(paste("Reference group not found:", grp_ctrl_name))
missing_cases <- setdiff(grp_case_names, df_input$Group)
if (length(missing_cases) > 0) stop(paste("Test groups not found:", paste(missing_cases, collapse=", ")))

# 4. Matrix Preparation
# ------------------------------------------------------------------------------
get_group_matrix <- function(df, group_names, features) {
  df %>% 
    filter(Group %in% group_names) %>% 
    select(all_of(features)) %>% 
    as.matrix()
}

mat_ctrl <- get_group_matrix(df_input, grp_ctrl_name, markers)
mat_case <- get_group_matrix(df_input, grp_case_names, markers)

message(sprintf("[Data] Control Matrix: %d samples | Case Matrix: %d samples", nrow(mat_ctrl), nrow(mat_case)))
if (nrow(mat_ctrl) < 5 || nrow(mat_case) < 5) stop("Insufficient sample size (<5) for network inference.")

# 5. Global Lambda (For Bootstrap Stability)
# ------------------------------------------------------------------------------
# We keep fixed lambda for Bootstrap to ensure we test stability of specific model
lambda_ctrl <- corpcor::estimate.lambda(mat_ctrl, verbose = FALSE)
lambda_case <- corpcor::estimate.lambda(mat_case, verbose = FALSE)
message(sprintf("[Config] Lambda (Stability): Control=%.4f | Case=%.4f", lambda_ctrl, lambda_case))

n_cores_req <- if(config$stats$n_cores == "auto") parallel::detectCores() - 1 else config$stats$n_cores
cl <- makeCluster(n_cores_req)
registerDoParallel(cl)
clusterExport(cl, varlist = c("boot_worker_pcor", "infer_network_pcor"), envir = .GlobalEnv)
clusterEvalQ(cl, { library(corpcor); library(dplyr) })

# 6. Bootstrap Inference (Stability)
# ------------------------------------------------------------------------------
run_parallel_bootstrap <- function(data_mat, n_boot, seed_val, fixed_lambda) {
  n_samples <- nrow(data_mat)
  parallel::clusterSetRNGStream(cl, seed_val)
  foreach(i = 1:n_boot, .packages = c("corpcor")) %dopar% {
    boot_worker_pcor(data_mat, n_samples, lambda_val = fixed_lambda)
  }
}

message("[Inference] Bootstrapping Control...")
boot_ctrl <- run_parallel_bootstrap(mat_ctrl, config$stats$n_boot, config$stats$seed, lambda_ctrl)
res_ctrl <- aggregate_boot_results(boot_ctrl, alpha = config$stats$alpha)

message("[Inference] Bootstrapping Case...")
boot_case <- run_parallel_bootstrap(mat_case, config$stats$n_boot, config$stats$seed + 1000, lambda_case)
res_case <- aggregate_boot_results(boot_case, alpha = config$stats$alpha)

# 7. Permutation Test (Differential)
# ------------------------------------------------------------------------------
message("[Inference] Running Permutation Test...")
n_perm <- config$stats$n_perm
pool_mat <- rbind(mat_ctrl, mat_case)
n1 <- nrow(mat_ctrl); n2 <- nrow(mat_case)

# [FIX] Use Dynamic Lambda (NULL) for Observed Difference to match Permutations logic
# This reproduces the "Old" behavior where lambda was re-estimated
obs_pcor_ctrl <- infer_network_pcor(mat_ctrl, fixed_lambda = NULL)
obs_pcor_case <- infer_network_pcor(mat_case, fixed_lambda = NULL)
obs_diff_mat  <- abs(obs_pcor_ctrl - obs_pcor_case)

clusterSetRNGStream(cl, config$stats$seed + 2000)

# Return a LIST of matrices, not a flattened vector
null_diffs <- foreach(i = 1:n_perm, .packages = "corpcor") %dopar% {
  shuffled <- sample(1:(n1 + n2))
  p1 <- pool_mat[shuffled[1:n1], ]
  p2 <- pool_mat[shuffled[(n1 + 1):(n1 + n2)], ]
  
  # Dynamic lambda for null distribution (Standard)
  r1 <- infer_network_pcor(p1, fixed_lambda = NULL)
  r2 <- infer_network_pcor(p2, fixed_lambda = NULL)
  abs(r1 - r2) 
}
stopCluster(cl)

# 8. FDR & Export
# ------------------------------------------------------------------------------
union_adj <- (res_ctrl$adj == 1 | res_case$adj == 1)
edges_idx <- which(upper.tri(union_adj) & union_adj, arr.ind = TRUE)
results_table <- data.frame()

if(nrow(edges_idx) > 0) {
  denom <- n_perm + 1
  
  # [FIX] Handle list structure of null_diffs correctly
  for(k in 1:nrow(edges_idx)) {
    i <- edges_idx[k,1]; j <- edges_idx[k,2]
    obs_val <- obs_diff_mat[i,j]
    
    # Extract value from each matrix in the list
    null_vals <- sapply(null_diffs, function(m) m[i, j])
    
    p_val <- (sum(null_vals >= obs_val) + 1) / denom
    
    results_table <- rbind(results_table, data.frame(
      Node1 = rownames(union_adj)[i], Node2 = rownames(union_adj)[j],
      Weight_Ctrl = res_ctrl$weights[i,j], Weight_Case = res_case$weights[i,j],
      Is_Stable_Ctrl = res_ctrl$adj[i,j] == 1,
      Is_Stable_Case = res_case$adj[i,j] == 1,
      Diff_Score = obs_val, P_Value = p_val
    ))
  }
  
  if(nrow(results_table) > 0) {
    results_table$FDR <- p.adjust(results_table$P_Value, method = "BH")
    results_table$Significant_Diff <- results_table$FDR < config$stats$fdr_threshold
    results_table <- results_table %>% arrange(P_Value)
  }
}

out_dir <- file.path(config$output_root, "04_network_inference")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

wb_diff <- createWorkbook()
addWorksheet(wb_diff, "Differential_Edges")
if (nrow(results_table) > 0) {
  writeData(wb_diff, "Differential_Edges", results_table)
} else {
  writeData(wb_diff, "Differential_Edges", data.frame(Message = "No significant edges found"))
}
saveWorkbook(wb_diff, file.path(out_dir, "Differential_Stats.xlsx"), overwrite = TRUE)

final_obj <- list(ctrl_network = res_ctrl, case_network = res_case, 
                  diff_table = results_table, config = config)
saveRDS(final_obj, file.path(out_dir, "inference_results.rds"))

message("=== STEP 4 COMPLETE ===\n")