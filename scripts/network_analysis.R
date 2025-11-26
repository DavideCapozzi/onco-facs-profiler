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

setwd("/home/davidec/projects/compositional_analysis/")

CONFIG <- list(
  # Input File
  FILENAME = "/mnt/c/Users/Davide/Desktop/Davide/bandi/0dottorato/projects/long survivors/DB_anonimo_standardizzato.xlsx",
  
  # Group Definitions
  GROUP_CASE_LIST = c("NSCLC", "HNSCC"), 
  GROUP_CTRL = "Healthy_Donors",
  
  # Data Cleaning & Imputation (POINT A)
  # We relax the completeness threshold to keep more biological markers
  MIN_COMPLETENESS = 0.60,       # Keep markers present in at least 60% of samples (was 70%)
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
# 2. DATA LOADING & MERGING (Point A - Implementation)
# ------------------------------------------------------------------------------

# Function to safely load and force numeric conversion
load_and_clean_sheet <- function(file, sheet) {
  cat(paste("    Loading sheet:", sheet, "\n"))
  
  # 1. Read Excel treating various "NA" strings as missing values
  df <- tryCatch(
    read_excel(file, sheet = sheet, na = c("", "NA", "na", "NaN", "nd", "n.d.")), 
    error = function(e) { cat("Error reading sheet\n"); return(NULL) }
  )
  
  if(is.null(df)) return(NULL)
  
  # 2. Clean Column Names (remove spaces/hidden chars)
  names(df) <- trimws(names(df))
  
  # 3. Force Numeric Conversion for Marker Columns (exclude Patient_ID)
  # This fixes the "need numeric data" error
  cols_to_fix <- setdiff(names(df), "Patient_ID")
  
  for(col in cols_to_fix) {
    # Suppress warnings about NAs introduced by coercion
    df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
  }
  
  return(df)
}

# --- A. Load Datasets ---
cat(">>> Loading Data...\n")
df_nsclc <- load_and_clean_sheet(CONFIG$FILENAME, "NSCLC")
df_hnscc <- load_and_clean_sheet(CONFIG$FILENAME, "HNSCC")
df_healthy <- load_and_clean_sheet(CONFIG$FILENAME, CONFIG$GROUP_CTRL)

# --- B. Merge Tumor Groups ---
# Identify common columns first to avoid mismatch errors
common_cols_tumor <- intersect(names(df_nsclc), names(df_hnscc))
cat(paste("    Common columns between NSCLC and HNSCC:", length(common_cols_tumor), "\n"))

df_tumor_raw <- rbind(
  df_nsclc[, common_cols_tumor],
  df_hnscc[, common_cols_tumor]
)

cat(paste(">>> Tumor Group Merged (N=", nrow(df_tumor_raw), ")\n"))

# --- C. Feature Selection 

# 1. Define universe of markers (Intersection of Tumor and Healthy columns)
all_potential_markers <- intersect(names(df_tumor_raw), names(df_healthy))
all_potential_markers <- setdiff(all_potential_markers, "Patient_ID")

cat("    Potential common markers found:", length(all_potential_markers), "\n")

# 2. Check Completeness with a Lower Threshold (0.5 instead of 0.6)
# Since HNSCC is sparse, we must be more tolerant to avoid dropping everything.
TOLERANCE_THRESHOLD <- 0.50 

check_completeness <- function(df, markers, threshold, group_name) {
  valid_list <- c()
  for(m in markers) {
    comp <- mean(!is.na(df[[m]]))
    if(comp >= threshold) {
      valid_list <- c(valid_list, m)
    } else {
      # Optional: Print what is being dropped for debug
      # cat(paste("      Dropping", m, "in", group_name, "- Completeness:", round(comp*100,1), "%\n"))
    }
  }
  return(valid_list)
}

cat(">>> Filtering Markers (Threshold:", TOLERANCE_THRESHOLD*100, "%)...\n")
valid_tumor <- check_completeness(df_tumor_raw, all_potential_markers, TOLERANCE_THRESHOLD, "Tumor")
valid_healthy <- check_completeness(df_healthy, all_potential_markers, TOLERANCE_THRESHOLD, "Healthy")

final_markers <- intersect(valid_tumor, valid_healthy)

cat(paste(">>> FINAL MARKERS RETAINED:", length(final_markers), "\n"))
print(final_markers) # Stampa i nomi per verifica

if(length(final_markers) < 10) {
  warning("Warning: Number of markers is still low. Consider lowering TOLERANCE_THRESHOLD to 0.4")
}

# --- D. Imputation (Safe Version) ---

impute_data <- function(df, markers) {
  d <- df[, markers, drop=FALSE] # drop=FALSE ensures it stays a dataframe
  
  for(i in 1:ncol(d)) {
    vals <- d[[i]] # Extract column vector
    
    if(all(is.na(vals))) {
      # If a column is ALL NA (shouldn't happen with filter, but just in case)
      d[[i]] <- 0 
    } else if(any(is.na(vals))) {
      # Impute with median
      med_val <- median(vals, na.rm = TRUE)
      d[is.na(vals), i] <- med_val
    }
  }
  return(d)
}

cat(">>> Applying Median Imputation...\n")
data_tumor_imp <- impute_data(df_tumor_raw, final_markers)
data_healthy_imp <- impute_data(df_healthy, final_markers)

cat(">>> Data Preparation Complete.\n")

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

cat(">>> Analysis Complete. Results saved in", CONFIG$OUTPUT_DIR, "\n")
