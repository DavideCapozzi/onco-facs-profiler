# src/02_exploratory.R
# ==============================================================================
# STEP 02: EXPLORATORY DATA ANALYSIS & ROBUST HYPOTHESIS TESTING
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(FactoMineR)
  library(factoextra)
  library(openxlsx)
  library(vegan)
})

# Source modules
source("R/utils_io.R")          
source("R/modules_hypothesis.R") 
source("R/modules_viz.R")

message("\n=== PIPELINE STEP 2: EXPLORATORY & STATS (Hybrid/Z-Score) ===")

# 1. Load Config & Data
config <- load_config("config/global_params.yml")
input_file <- file.path(config$output_root, "01_QC", "data_processed.rds")

if (!file.exists(input_file)) stop("Step 01 output not found. Run src/01_ingest.R first.")

DATA <- readRDS(input_file)

# MAPPING NEW DATA STRUCTURE (from Step 1)
# ------------------------------------------------------------------------------
# hybrid_data_z: Matrice Z-scored completa (Markers funzionali + CLR).
#                Include Patient_ID e Group.
df_global <- DATA$hybrid_data_z 

# ilr_balances: Lista di matrici ILR per i sottogruppi (es. Differentiation).
#               NOTA: Queste matrici non hanno metadati dentro (sono numeric matrix),
#               dobbiamo ri-associarli per i test.
ilr_list <- DATA$ilr_balances 

meta_cols <- c("Patient_ID", "Group")
metadata  <- DATA$metadata
my_colors <- get_palette(config)

out_dir <- file.path(config$output_root, "02_exploratory")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

message(sprintf("[Data] Global Matrix: %d Samples x %d Markers (Z-Scored)", 
                nrow(df_global), ncol(df_global) - length(meta_cols)))

# ==============================================================================
# PART A: PCA (Global Hybrid View)
# ==============================================================================
message("[PCA] Running PCA on Global Z-Scored Matrix...")

# Prepare Matrix (Exclude metadata)
pca_input <- df_global %>% select(-all_of(meta_cols)) %>% as.matrix()

# Run PCA (FactoMineR)
# scale.unit = FALSE perché i dati sono già Z-scored, ma TRUE non fa danni (1/1=1)
res_pca <- PCA(pca_input, scale.unit = FALSE, graph = FALSE)

# Plotting
p_pca <- plot_pca_custom(res_pca, df_global, my_colors, show_labels = FALSE) # Set TRUE if few samples
ggsave(file.path(out_dir, "PCA_Global_Score.pdf"), p_pca, width = 8, height = 6)

p_scree <- fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 50)) + 
  theme_minimal() + ggtitle("Scree Plot (Global)")
ggsave(file.path(out_dir, "PCA_Global_Scree.pdf"), p_scree, width = 6, height = 4)

p_vars <- fviz_pca_var(res_pca, col.var = "contrib", 
                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                       repel = TRUE) + ggtitle("Variables Contribution")
ggsave(file.path(out_dir, "PCA_Global_Variables.pdf"), p_vars, width = 8, height = 6)


# ==============================================================================
# PART B: STATISTICAL TESTING (PERMANOVA)
# ==============================================================================
# Inizializziamo un Workbook Excel per salvare tutti i risultati
wb <- createWorkbook()

# 1. GLOBAL PERMANOVA (Is there a generic difference?)
# ------------------------------------------------------------------------------
message("\n[Stats] Running Global PERMANOVA (All Markers)...")
perm_global <- test_coda_permanova(df_global, group_col = "Group", 
                                   metadata_cols = "Patient_ID", 
                                   n_perm = config$stats$n_perm)

message(sprintf("   -> Global p-value: %.5f (R2: %.2f%%)", 
                perm_global$`Pr(>F)`[1], perm_global$R2[1] * 100))

addWorksheet(wb, "Global_PERMANOVA")
writeData(wb, "Global_PERMANOVA", as.data.frame(perm_global), rowNames = TRUE)


# 2. LOCAL PERMANOVA (Sub-groups ILR Balances)
# ------------------------------------------------------------------------------
# Testiamo se specifici assi composizionali (es. Naive vs Memory) sono diversi.
# Questo sostituisce la MANOVA.

if (length(ilr_list) > 0) {
  message("\n[Stats] Running Local PERMANOVA on Compositional Groups (ILR)...")
  
  summary_local <- data.frame()
  
  for (grp_name in names(ilr_list)) {
    # Recupera la matrice ILR
    mat_ilr <- ilr_list[[grp_name]]
    
    # Ricostruisci un df con i metadati per la funzione permanova
    # Assumiamo che l'ordine delle righe sia preservato (Step 1 garantisce questo)
    df_test <- cbind(metadata, as.data.frame(mat_ilr))
    
    # Run Test
    res_perm <- test_coda_permanova(df_test, group_col = "Group", 
                                    metadata_cols = "Patient_ID", 
                                    n_perm = config$stats$n_perm)
    
    # Store Summary
    pval <- res_perm$`Pr(>F)`[1]
    r2   <- res_perm$R2[1]
    
    message(sprintf("   -> Group '%s': p = %.5f", grp_name, pval))
    
    summary_local <- rbind(summary_local, data.frame(
      SubGroup = grp_name,
      P_Value = pval,
      R2_Percent = r2 * 100,
      F_Model = res_perm$F[1]
    ))
    
    # Save detailed result to Excel sheet
    sheet_name <- substr(paste0("ILR_", grp_name), 1, 31) # Excel limit 31 chars
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, as.data.frame(res_perm), rowNames = TRUE)
  }
  
  # Save Summary Table
  addWorksheet(wb, "Summary_Local_Tests")
  writeData(wb, "Summary_Local_Tests", summary_local)
}

# ==============================================================================
# PART C: SAVE OUTPUT
# ==============================================================================

excel_path <- file.path(out_dir, "Statistical_Test_Results.xlsx")
saveWorkbook(wb, excel_path, overwrite = TRUE)

message(sprintf("\n[Output] Results saved to: %s", out_dir))
message("=== STEP 2 COMPLETE ===\n")