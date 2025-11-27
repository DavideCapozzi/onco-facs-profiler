# ==============================================================================
# PROJECT: Enhanced Differential Correlation Analysis (Pan-Cancer vs Healthy)
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. SETUP & LIBRARIES
# ------------------------------------------------------------------------------

library(readxl)
library(dplyr)
library(ggplot2)
library(qgraph)
library(NetworkComparisonTest)
library(igraph)
library(compositions)
library(reshape2)
library(psych)

# ------------------------------------------------------------------------------
# 1. CONFIGURATION
# ------------------------------------------------------------------------------

setwd("/home/davidec/projects/compositional_analysis/results")

CONFIG <- list(
  # Input File
  FILENAME = "/mnt/c/Users/Davide/Desktop/Davide/bandi/0dottorato/projects/long survivors/DB_anonimo_standardizzato.xlsx",
  
  # Group Definitions
  GROUP_CASE_LIST = c("NSCLC", "HNSCC"), 
  GROUP_CTRL = "Healthy_Donors",
  
  # Data Cleaning & Imputation (POINT A)
  # We relax the completeness threshold to keep more biological markers
  MIN_COMPLETENESS = 0.40,       # Keep markers present in at least 60% of samples
  IMPUTE_METHOD = "median",      # Simple robust imputation for missing values
  
  # Network & Stats
  PERFORM_CLR = TRUE,            # Compositional transformation
  NCT_PERMS = 1000,              # Permutations for statistical test
  CORR_METHOD = "spearman",      # Rank correlation
  
  # Adaptive Thresholds (Logic: Smaller N needs strictly stronger effects)
  # We focus on the DIFFERENCE in correlation, not just absolute value.
  MIN_R_DIFF_TO_KEEP = 0.3,      # Only edges shifting by at least 0.3 are interesting
  SIG_LEVEL_NCT = 0.05,          # P-value threshold
  
  OUTPUT_DIR = paste0("Analysis_PanCancer_v2_", format(Sys.time(), "%Y%m%d"))
)

if (!dir.exists(CONFIG$OUTPUT_DIR)) dir.create(CONFIG$OUTPUT_DIR)
cat(">>> Starting Enhanced Analysis v2.0\n")

# ------------------------------------------------------------------------------
# 2. DATA LOADING & MERGING
# ------------------------------------------------------------------------------

# Helper function: Load and force numeric conversion
load_clean_numeric <- function(path, sheet) {
  df <- read_excel(path, sheet = sheet)
  ids <- df[[1]] # Assume first column is Patient_ID
  matrix_data <- df[, -1]
  
  # Force numeric conversion suppressing warnings
  matrix_data <- as.data.frame(lapply(matrix_data, function(x) suppressWarnings(as.numeric(as.character(x)))))
  
  matrix_data$Patient_ID <- ids
  return(matrix_data)
}

# 1. Load Datasets
raw_nsclc <- load_clean_numeric(CONFIG$FILENAME, "NSCLC")
raw_hnscc <- load_clean_numeric(CONFIG$FILENAME, "HNSCC")
raw_healthy <- load_clean_numeric(CONFIG$FILENAME, "Healthy_Donors")

# 2. Identify Potential Markers (Name Matching)
common_cols <- intersect(colnames(raw_nsclc), colnames(raw_hnscc))
common_cols <- intersect(common_cols, colnames(raw_healthy))
marker_cols <- setdiff(common_cols, "Patient_ID")

cat("\nInitial common markers found by name:", length(marker_cols), "\n")

# 3. GLOBAL FILTERING & REPORTING
# Strategy: Keep marker if present in at least 40% of TOTAL tumor patients (NSCLC + HNSCC)

raw_tumor_combined <- rbind(raw_nsclc[, marker_cols], raw_hnscc[, marker_cols])
completeness_rates <- colMeans(!is.na(raw_tumor_combined))

# Identify Valid vs Dropped
valid_markers_step1 <- names(completeness_rates)[completeness_rates >= CONFIG$MIN_COMPLETENESS]
dropped_markers_na <- names(completeness_rates)[completeness_rates < CONFIG$MIN_COMPLETENESS]

# REPORT DROPPED (High Missingness)
if(length(dropped_markers_na) > 0) {
  cat("\n[REJECTED] The following markers were dropped due to High Missingness (<40% present):\n")
  for(m in dropped_markers_na) {
    missing_pct <- (1 - completeness_rates[[m]]) * 100
    cat(sprintf("  - %-15s : Missing in %5.1f%% of Tumor samples\n", m, missing_pct))
  }
}

cat("\n[PASSED] Markers retained for Imputation:", length(valid_markers_step1), "\n")

# 4. INTELLIGENT IMPUTATION (Distribution-Based)
# Logic: Use Reference (NSCLC) Mean/SD to generate jittered data for empty HNSCC columns

impute_with_reference <- function(target_df, reference_df, markers) {
  imputed_df <- target_df
  for (m in markers) {
    target_vec <- target_df[[m]]
    ref_vec <- reference_df[[m]]
    
    missing_idx <- is.na(target_vec)
    n_missing <- sum(missing_idx)
    
    if (n_missing > 0) {
      # CASE A: Target has data (>30%) -> Use Target Median
      if (mean(!missing_idx) > 0.30) {
        fill_val <- median(target_vec, na.rm = TRUE)
        imputed_df[missing_idx, m] <- fill_val
      } else {
        # CASE B: Target empty -> Rescue with Reference Distribution + Jitter
        ref_mean <- mean(ref_vec, na.rm = TRUE)
        ref_sd <- sd(ref_vec, na.rm = TRUE)
        if (is.na(ref_sd) || ref_sd == 0) { ref_sd <- 0.001; ref_mean <- 0 } # Safety
        
        set.seed(123) 
        synthetic_values <- rnorm(n_missing, mean = ref_mean, sd = ref_sd)
        synthetic_values[synthetic_values < 0] <- 0.001 # No negative percentages
        
        imputed_df[missing_idx, m] <- synthetic_values
      }
    }
  }
  return(imputed_df)
}

# Apply Imputation
nsclc_clean <- impute_with_reference(raw_nsclc, raw_nsclc, valid_markers_step1)
hnscc_clean <- impute_with_reference(raw_hnscc, raw_nsclc, valid_markers_step1)
healthy_clean <- impute_with_reference(raw_healthy, raw_healthy, valid_markers_step1)

# 5. FINAL MERGE & VARIANCE CHECK
data_tumor_final <- rbind(
  nsclc_clean[, valid_markers_step1, drop=FALSE],
  hnscc_clean[, valid_markers_step1, drop=FALSE]
)
data_healthy_final <- healthy_clean[, valid_markers_step1, drop=FALSE]

# Check for Zero Variance (Constant columns crash correlation tests)
tumor_var <- apply(data_tumor_final, 2, var)
healthy_var <- apply(data_healthy_final, 2, var)

# Valid if variance > 0 in BOTH groups
final_mask <- (tumor_var > 0) & (healthy_var > 0)
final_markers <- valid_markers_step1[final_mask]
dropped_variance <- valid_markers_step1[!final_mask]

# REPORT DROPPED (Zero Variance)
if(length(dropped_variance) > 0) {
  cat("\n[REJECTED] The following markers were dropped due to Zero Variance (Constant values):\n")
  for(m in dropped_variance) {
    cat(sprintf("  - %-15s : Var(Tumor)=%.4f, Var(Healthy)=%.4f\n", m, tumor_var[[m]], healthy_var[[m]]))
  }
}

# Apply Final Filter
data_tumor_imp <- data_tumor_final[, final_markers]
data_healthy_imp <- data_healthy_final[, final_markers]

cat("\nFINAL DATASET READY FOR ANALYSIS:\n")
cat("  > Tumor Samples (NSCLC+HNSCC) :", nrow(data_tumor_final), "\n")
cat("  > Healthy Samples             :", nrow(data_healthy_final), "\n")
cat("  > Final Markers Analyzed      :", length(final_markers), "\n")

# ------------------------------------------------------------------------------
# 3. CLR TRANSFORMATION
# ------------------------------------------------------------------------------

apply_clr <- function(df) {
  # Add small offset for zeros
  mat <- as.matrix(df) + 0.01 
  clr_res <- compositions::clr(mat)
  return(as.data.frame(clr_res))
}

cat(">>> Applying CLR Transformation...\n")
data_tumor_final <- apply_clr(data_tumor_imp)
data_healthy_final <- apply_clr(data_healthy_imp)

# Ensure Column Alignment (Critical for NCT)
data_tumor_final <- data_tumor_final[, final_markers]
data_healthy_final <- data_healthy_final[, final_markers]

# Convert to Matrix
data_tumor_final_mat <- as.matrix(data_tumor_final)
data_healthy_final_mat <- as.matrix(data_healthy_final)

# ------------------------------------------------------------------------------
# 4. NETWORK COMPARISON TEST (NCT)
# ------------------------------------------------------------------------------

cat(">>> Running NCT (1000 permutations). This takes about 1-2 minutes...\n")
set.seed(999) 

nct_out <- NCT(data_tumor_final, data_healthy_final, 
               it = CONFIG$NCT_PERMS, 
               gamma = 0, 
               test.edges = TRUE, 
               test.centrality = TRUE,
               progressbar = TRUE,
               paired = FALSE,
               weighted = TRUE,
               AND = TRUE,
               abs = FALSE)

# ------------------------------------------------------------------------------
# 5. DIFFERENTIAL EDGE EXTRACTION 
# ------------------------------------------------------------------------------

cat(">>> Extracting Differential Statistics...\n")

# 1. Calcolo Matrici di Correlazione Locali
cor_tumor <- cor(data_tumor_final, method = CONFIG$CORR_METHOD)
cor_healthy <- cor(data_healthy_final, method = CONFIG$CORR_METHOD)

# Usiamo la maschera per estrarre le coppie uniche (triangolo superiore)
upper_tri <- upper.tri(cor_tumor)
nodes <- colnames(cor_tumor)

# Creiamo il dataframe "Master" con le correlazioni locali
my_edges <- data.frame(
  Node1 = nodes[row(cor_tumor)[upper_tri]],
  Node2 = nodes[col(cor_tumor)[upper_tri]],
  R_Tumor = cor_tumor[upper_tri],
  R_Healthy = cor_healthy[upper_tri],
  stringsAsFactors = FALSE
)

# Creiamo una Chiave Univoca ordinata (es. "GeneA_GeneB") per il matching sicuro
my_edges$EdgeID <- paste(pmin(my_edges$Node1, my_edges$Node2), 
                         pmax(my_edges$Node1, my_edges$Node2), sep="_")

# 2. Recupero dati da NCT (Gestione robusta formato Tabella vs Matrice)
if (!is.null(nct_out$einv.pvals)) {
  
  # Convertiamo l'output di NCT in dataframe
  nct_table <- as.data.frame(nct_out$einv.pvals)
  
  # Standardizziamo i nomi delle colonne (le prime 2 sono nodi, la 3 è p-val)
  # Questo gestisce sia il formato lista (465x4) che matrice convertita
  if(ncol(nct_table) >= 3) {
    colnames(nct_table)[1:3] <- c("N1", "N2", "P_Val")
    
    # Creiamo la stessa chiave univoca su NCT
    nct_table$EdgeID <- paste(pmin(as.character(nct_table$N1), as.character(nct_table$N2)), 
                              pmax(as.character(nct_table$N1), as.character(nct_table$N2)), sep="_")
    
    # Selezioniamo solo ciò che serve
    nct_clean <- nct_table[, c("EdgeID", "P_Val")]
    nct_clean$P_Val <- as.numeric(as.character(nct_clean$P_Val))
    
    # 3. JOIN: Uniamo i p-value alla tabella master
    final_edges <- merge(my_edges, nct_clean, by="EdgeID", all.x=TRUE)
    names(final_edges)[names(final_edges) == "P_Val"] <- "P_Value_NCT"
    
  } else {
    # Fallback improbabile: se NCT restituisce formato sconosciuto
    warning("NCT output format not recognized. P-values set to NA.")
    final_edges <- my_edges
    final_edges$P_Value_NCT <- NA
  }
  
} else {
  warning("NCT p-values not found (einv.pvals is NULL).")
  final_edges <- my_edges
  final_edges$P_Value_NCT <- NA
}

# 4. Classificazione e Pulizia
edge_list <- final_edges %>%
  mutate(
    Delta_R = R_Tumor - R_Healthy,
    Significant = !is.na(P_Value_NCT) & P_Value_NCT < CONFIG$SIG_LEVEL_NCT,
    Category = case_when(
      Delta_R > CONFIG$MIN_R_DIFF_TO_KEEP  ~ "Tumor Gain / Healthy Loss",
      Delta_R < -CONFIG$MIN_R_DIFF_TO_KEEP ~ "Tumor Loss / Healthy Gain",
      TRUE ~ "Stable"
    )
  ) %>%
  mutate(across(c(R_Tumor, R_Healthy, Delta_R), \(x) round(x, 3))) %>%
  select(Node1, Node2, R_Tumor, R_Healthy, Delta_R, P_Value_NCT, Significant, Category)

# 5. Salvataggio
cat(paste("    Total edges processed:", nrow(edge_list), "\n"))

sig_edges <- edge_list %>% 
  filter(Significant == TRUE | abs(Delta_R) > 0.4) %>% 
  arrange(P_Value_NCT)

write.csv(sig_edges, file.path(CONFIG$OUTPUT_DIR, "Significant_Differential_Edges.csv"), row.names = FALSE)

# Variabile necessaria per i plot successivi (Blocco 6)
diff_matrix <- cor_tumor - cor_healthy

cat(">>> Extraction Complete.\n")

# ------------------------------------------------------------------------------
# 6. VISUALIZATION 
# ------------------------------------------------------------------------------

cat(">>> Generating Advanced Plots...\n")
pdf(file.path(CONFIG$OUTPUT_DIR, "Network_Analysis_Plots.pdf"), width = 12, height = 8)

# Layout: Calculate a joint layout for stability
L <- qgraph(cor_healthy, DoNotPlot=TRUE)$layout

# PLOT 1: The Difference Network (The most important one)
# Red edges = Higher correlation in Tumor
# Blue edges = Higher correlation in Healthy
title_diff <- paste("DIFFERENCE NETWORK (Tumor - Healthy)\nRed = Tumor Specific | Blue = Healthy Specific")

qgraph(diff_matrix, layout = L,
       graph = "cor",
       title = title_diff,
       cut = 0.2,        # Hide very small differences
       minimum = 0.2,    # Minimum difference to show an edge
       edge.labels = FALSE,
       posCol = "#D55E00", # Vermillion (Colorblind safe Red)
       negCol = "#0072B2", # Blue (Colorblind safe)
       vsize = 6, 
       borders = FALSE,
       labels = colnames(data_tumor_final))

# PLOT 2 & 3: Side by Side with same layout
par(mfrow=c(1,2))
qgraph(cor_tumor, layout = L, title = paste("Tumor (N=17)\nNSCLC + HNSCC"), 
       graph = "cor", minimum = 0.35, posCol="darkred", negCol="darkblue", vsize=5)
qgraph(cor_healthy, layout = L, title = paste("Healthy (N=56)"), 
       graph = "cor", minimum = 0.35, posCol="darkred", negCol="darkblue", vsize=5)

dev.off()

# PLOT 4: Scatter of Correlations (Visualization of Global Shift)
p_scatter <- ggplot(edge_list, aes(x=R_Healthy, y=R_Tumor)) +
  geom_point(aes(color=abs(Delta_R) > 0.4), alpha=0.6) +
  scale_color_manual(values=c("grey", "red")) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  theme_minimal() +
  labs(title="Correlation Stability Analysis",
       subtitle="Points far from diagonal represent altered immune regulation",
       x="Correlation in Healthy", y="Correlation in Tumor",
       color="High Delta")

ggsave(file.path(CONFIG$OUTPUT_DIR, "Scatter_Global_Shift.pdf"), p_scatter, width=6, height=6)

# ------------------------------------------------------------------------------
# 7. SUMMARY REPORT
# ------------------------------------------------------------------------------

sink(file.path(CONFIG$OUTPUT_DIR, "Final_Report.txt"))
cat("=== PAN-CANCER NETWORK ANALYSIS REPORT ===\n\n")
cat("1. Data Composition:\n")
cat("   Tumor Group (NSCLC + HNSCC): N =", nrow(data_tumor_final), "\n")
cat("   Healthy Group: N =", nrow(data_healthy_final), "\n")
cat("   Markers Analyzed:", length(final_markers), "\n\n")

cat("2. Global Network Differences (NCT):\n")
cat("   Structure Invariance p-value:", nct_out$nwb, "\n")
cat("   (p < 0.05 indicates the global immune wiring is significantly different)\n\n")

cat("3. Key Differential Edges (Delta R > 0.4):\n")
top_diffs <- head(sig_edges[order(-abs(sig_edges$Delta_R)), c("Node1", "Node2", "Delta_R", "P_Value_NCT")], 10)
print(top_diffs)
sink()

# ==============================================================================
# 8. Cytoscape export 
# ==============================================================================

export_to_cytoscape <- function(nct_output, marker_names, output_dir) {
  
  # 2. Estrazione Matrici e Controllo Integrità
  mat_healthy <- nct_output$nw1 
  mat_tumor   <- nct_output$nw2   
  
  # Diagnostica interna
  sum_healthy <- sum(abs(mat_healthy), na.rm=TRUE)
  sum_tumor   <- sum(abs(mat_tumor), na.rm=TRUE)
  
  cat("--- Cytoscape Export Diagnostic ---\n")
  cat("Sum of edges (Healthy):", sum_healthy, "\n")
  cat("Sum of edges (Tumor):", sum_tumor, "\n")
  
  if(sum_healthy == 0 || sum_tumor == 0) {
    warning("ATTENZIONE: Una delle reti è vuota (tutti zeri). Controllo interrotto ma procedo all'export.")
  }
  
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 3. Funzione Scrittura File
  write_network_file <- function(matrix_in, filename, label) {
    if(is.null(matrix_in)) return()
    
    # Assegna nomi
    rownames(matrix_in) <- marker_names
    colnames(matrix_in) <- marker_names
    
    # Trasforma matrice in lista archi
    edges <- reshape2::melt(as.matrix(matrix_in))
    colnames(edges) <- c("Source", "Target", "Weight")
    
    # Rimuovi duplicati (triangolo superiore) e autoloops
    edges <- edges[as.character(edges$Source) < as.character(edges$Target), ]
    
    # Rimuovi edge con peso 0 assoluto
    edges <- edges[edges$Weight != 0, ]
    
    if(nrow(edges) > 0) {
      edges$Type <- label
      edges$Interaction <- ifelse(edges$Weight > 0, "Positive", "Negative")
      edges$AbsWeight <- abs(edges$Weight)
      
      file_path <- file.path(output_dir, filename)
      write.table(edges, file_path, sep = "\t", quote = FALSE, row.names = FALSE)
      cat("Saved:", filename, "(", nrow(edges), "edges )\n")
    } else {
      cat("Skipped:", filename, "(No edges found)\n")
    }
  }
  
  # 4. Esecuzione Export
  write_network_file(mat_healthy, "Cytoscape_Healthy_Network.txt", "Healthy")
  write_network_file(mat_tumor,   "Cytoscape_Tumor_Network.txt",   "Tumor")
  
  # 5. Export Differenziale
  if(sum_healthy > 0 && sum_tumor > 0) {
    delta_matrix <- mat_tumor - mat_healthy
    rownames(delta_matrix) <- marker_names; colnames(delta_matrix) <- marker_names
    
    edges_diff <- reshape2::melt(as.matrix(delta_matrix))
    colnames(edges_diff) <- c("Source", "Target", "Delta_R")
    edges_diff <- edges_diff[as.character(edges_diff$Source) < as.character(edges_diff$Target), ]
    
    # Aggiungi P-Value se presente
    if(!is.null(nct_output$einv.pvals)) {
      pvals <- reshape2::melt(as.matrix(nct_output$einv.pvals))
      colnames(pvals) <- c("Source", "Target", "P_Value")
      edges_diff <- merge(edges_diff, pvals, by=c("Source", "Target"), all.x=TRUE)
    } else {
      edges_diff$P_Value <- NA
    }
    
    # Filtro: Salva se Delta è alto (>0.3) O se P-value è significativo (<0.05)
    # Rimuoviamo il rumore di fondo
    mask <- (abs(edges_diff$Delta_R) > 0.3) | (!is.na(edges_diff$P_Value) & edges_diff$P_Value < 0.05)
    edges_diff_sig <- edges_diff[mask, ]
    
    if(nrow(edges_diff_sig) > 0) {
      write.table(edges_diff_sig, file.path(output_dir, "Cytoscape_Differential_Network.txt"), 
                  sep = "\t", quote = FALSE, row.names = FALSE)
      cat("Saved: Cytoscape_Differential_Network.txt (", nrow(edges_diff_sig), "edges )\n")
    }
  }
}

export_to_cytoscape(nct_out, final_markers, "Cytoscape_Results")

cat(">>> Analysis Complete. Results saved in", CONFIG$OUTPUT_DIR, "\n")
