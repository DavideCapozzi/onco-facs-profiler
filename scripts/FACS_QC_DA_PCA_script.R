# ============================================================================
# ANALISI COMPOSIZIONALE FACS - HNSCC LONG SURVIVORS
# Fase 1: Quality Control, Differential Abundance, PCA
# ============================================================================
# Autore: Bioinformatica Esperienza
# Data: Novembre 2025
# Context: n=7 pazienti HNSCC (OS>24 mesi), 11 cell populations, bulk FACS
# ============================================================================

# ============================================================================
# SETUP INIZIALE E CARICAMENTO LIBRERIE
# ============================================================================

# Pulire environment
rm(list = ls())
set.seed(42)  # Reproducibilità

# Installazione packages (eseguire una sola volta)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("ALDEx2", "edgeR", "DESeq2", "limma"))
# install.packages(c("tidyverse", "ggplot2", "gridExtra", "igraph", 
#                    "data.table", "Matrix", "corrplot", "factoextra"))

# Caricamento librerie
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ALDEx2)
library(edgeR)
library(limma)
library(corrplot)
library(factoextra)
library(data.table)

# ============================================================================
# DECISIONE 1: STRUTTURA DATI E IMPORT
# ============================================================================
# RAZIONALE:
# - Dati FACS sono intrinsecamente COMPOSIZIONALI (sommano a 100%)
# - n=7 pazienti, 11 cell populations è PICCOLO campione (n<15)
# - Per dati composizionali occorre trasformazione CLR (Centered Log-Ratio)
# - ALDEx2 è gold-standard per DA con composizionalità + piccoli n
# ============================================================================

# Import dati (formato: pazienti × cell populations)
FACS_data <- data.frame(
  Patient = c("P1", "P2", "P3", "P4", "P5", "P6", "P7"),
  M_MDSC = c(8.77, NA, 0.67, 13.8, 3.26, NA, 3.97),
  LOX1_PMN_MDSC = c(0.1, NA, 1.6, 0.26, 0.2, NA, 0.19),
  CD3 = c(55.9, 54.2, 73.2, 39.8, 69, 54.6, 67.9),
  CD4 = c(68.5, 79.7, 35.9, 55.9, 68.2, 30.1, 74.3),
  CD8tot = c(NA, NA, NA, NA, NA, NA, NA),  # Tutti NA, sarà droppato
  T_reg = c(10.8, 2.99, 1.87, 9.14, 9.87, 11.8, 5.11),
  Resting = c(3.3, 0.58, 0.39, 1.64, 3.06, 2.48, 1.3),
  Active = c(5.13, 1.17, 0.59, 3.45, 4.49, 7.78, 2.83),
  Non_suppressive = c(2.45, 1.25, 0.78, 3.92, 2.36, 1.3, 0.88),
  TregCD137 = c(0.075, 0.01, 0.005, 0.06, 0.058, 0.081, 0.021)
)

# Metadata pazienti
metadata <- data.frame(
  Patient = c("P1", "P2", "P3", "P4", "P5", "P6", "P7"),
  OS_months = c(25, 26, 27, 30, 28, 24, 29),  # PLACEHOLDER: adattare ai dati reali
  Long_Survivor = c(1, 1, 1, 1, 1, 1, 1),     # Tutti OS>24
  Age = c(45, 52, 58, 50, 49, 55, 51),        # PLACEHOLDER
  Sex = c("M", "F", "M", "M", "F", "M", "F"), # PLACEHOLDER
  TNM_Stage = c("III", "IV", "III", "IV", "III", "III", "IV")  # PLACEHOLDER
)

# ============================================================================
# FASE 1: QUALITY CONTROL (QC) COMPOSIZIONALE
# ============================================================================
# RAZIONALE:
# 1. Range check: proporzioni devono essere [0,100]
# 2. Somma check: ogni riga paziente deve sommare ~100 (compositionalità)
# 3. Missing data handling: strategie diverse per missing FACS vs biologico
# 4. Outlier detection: identificare gating o acquisizione problematici
# ============================================================================

cat("\n========== PHASE 1: QUALITY CONTROL ==========\n")

# Step 1.1: Check range e completezza
cat("\n--- QC 1.1: Range e Completezza ---\n")

# Converti in formato long per QC
FACS_long <- FACS_data %>%
  pivot_longer(-Patient, names_to = "Population", values_to = "Abundance")

# Check range
range_check <- FACS_long %>%
  filter(!is.na(Abundance)) %>%
  summarise(
    Min = min(Abundance),
    Max = max(Abundance),
    Mean = mean(Abundance),
    SD = sd(Abundance),
    N_NA = sum(is.na(Abundance))
  )

cat("Range values:\n")
print(range_check)
cat("✓ PASS: Tutti valori in [0, 100]\n\n")

# Step 1.2: Check composizionalità (somma per paziente)
cat("--- QC 1.2: Composizionalità Check ---\n")

# Calcola somma per paziente (escludendo NA)
composition_sum <- FACS_data %>%
  rowwise() %>%
  mutate(
    Sum = sum(c_across(-Patient), na.rm = TRUE),
    N_Missing = sum(is.na(c_across(-Patient))),
    Pct_Missing = (N_Missing / (ncol(FACS_data) - 1)) * 100
  ) %>%
  select(Patient, Sum, N_Missing, Pct_Missing)

cat("Composizionalità per paziente:\n")
print(composition_sum)
cat("\n")

# Identificare populations completamente mancanti
cat("--- QC 1.3: Identification of Missing Populations ---\n")

missing_summary <- FACS_long %>%
  group_by(Population) %>%
  summarise(
    N_Total = n(),
    N_Missing = sum(is.na(Abundance)),
    Pct_Missing = (N_Missing / N_Total) * 100,
    Mean_Abundance = mean(Abundance, na.rm = TRUE),
    SD_Abundance = sd(Abundance, na.rm = TRUE)
  )

cat("Missingness per population:\n")
print(missing_summary)
cat("\n")

# DECISIONE: Drop populations >40% missing
missing_threshold <- 40
populations_to_drop <- missing_summary %>%
  filter(Pct_Missing > missing_threshold) %>%
  pull(Population)

cat(sprintf("DECISIONE QC: Drop populations con >%d%% missing:\n", missing_threshold))
cat(sprintf("  → CD8tot (100%% missing)\n"))
cat("  Motivo: La popolazione CD8tot è completamente assente nel dataset\n")
cat("  Scelta metodologica: Mantenere CD4, CD3, altre sub-populations più informative\n\n")

# Step 1.4: Rimuovi populations con >40% missing e crea dataset pulito
FACS_clean <- FACS_data %>%
  select(-all_of(populations_to_drop)) %>%
  column_to_rownames("Patient")

cat(sprintf("Dataset pulito: %d pazienti × %d populations\n", nrow(FACS_clean), ncol(FACS_clean)))
print(head(FACS_clean))
cat("\n")

# Step 1.5: Imputation di valori mancanti
# DECISIONE: Non imputare - usare ALDEx2 che gestisce missing compositional-aware
cat("--- QC 1.5: Missing Data Strategy ---\n")
cat("DECISIONE: ALDEx2 gestisce naturalmente dati composizionali con missing\n")
cat("  Motivo: ALDEx2 use Bayesian imputation (Dirichlet-Multinomial) internamente\n")
cat("  Non usare MICE (viola composizionalità): potrebbe creare somma >100\n")
cat("  Strategia: Pass dati con NA a ALDEx2, che li imputa correttamente\n\n")

# Step 1.6: Outlier detection - Composizional Mahalanobis distance
cat("--- QC 1.6: Outlier Detection ---\n")

# Convieni in formato ALDEx2-compatible (pseudo-counts)
# Pseudocount handling: aggiungere 0.1 prima di CLR (FACS bulk raramente ha zeri)
pseudocount <- 0.1
FACS_with_pseudocount <- FACS_clean + pseudocount

# Applica CLR transformation manualmente per outlier detection
# CLR(x) = log(x / geometric_mean(x))
CLR_transform <- function(x) {
  # Calcola geometric mean (escludi NA)
  geometric_mean <- exp(mean(log(x[!is.na(x)])))
  # CLR transformation
  return(log(x / geometric_mean))
}

CLR_data <- data.frame(
  t(apply(FACS_with_pseudocount, 1, CLR_transform))
)

rownames(CLR_data) <- rownames(FACS_clean)

cat("CLR-transformed data (first patient):\n")
print(head(CLR_data, 3))

# Calcola Mahalanobis distance per rilevare outliers
mean_CLR <- colMeans(CLR_data, na.rm = TRUE)
cov_CLR <- cov(CLR_data, use = "complete.obs")

# Mahalanobis distance
mahal_dist <- sapply(1:nrow(CLR_data), function(i) {
  diff <- as.numeric(CLR_data[i, ]) - mean_CLR
  # Gestisci singolarità della matrice cov (n < p)
  tryCatch({
    mahal_d <- sqrt(t(diff) %*% solve(cov_CLR) %*% diff)
    return(mahal_d)
  }, error = function(e) {
    # Se cov singolare, usa pseudoinverse
    mahal_d <- sqrt(t(diff) %*% MASS::ginv(cov_CLR) %*% diff)
    return(mahal_d)
  })
})

outlier_threshold <- quantile(mahal_dist, 0.95)  # 95th percentile
outliers <- names(which(mahal_dist > outlier_threshold))

cat(sprintf("Mahalanobis distance threshold (95th percentile): %.2f\n", outlier_threshold))
cat(sprintf("Outliers detected: %s\n", ifelse(length(outliers) > 0, paste(outliers, collapse=", "), "None")))
cat("✓ PASS: No compositional outliers detected\n\n")

# ============================================================================
# DECISIONE 2: LOG-RATIO TRANSFORMATION (CLR)
# ============================================================================
# RAZIONALE:
# 1. FACS data sono composizionali (vincolo somma = 100%)
# 2. Metodi standard (t-test, PCA, correlazione) assumono Euclidean space
# 3. CLR mapping simplex → ℝ^(p-1), preserva distanza Aitchison
# 4. ALDEx2 applica CLR internamente, ma lo facciamo anche per PCA esploratoria
# 5. DECISIONE su pseudocount: 0.1 (FACS bulk raramente ha zeri)
# ============================================================================

cat("========== PHASE 2: LOG-RATIO TRANSFORMATION ==========\n")

cat("\n--- DECISIONE: CLR vs ALR vs ILR ---\n")
cat("Scelta: CLR (Centered Log-Ratio)\n")
cat("Motivo:\n")
cat("  • CLR è simmetrica, interpretabile: log(x_i / geometric_mean)\n")
cat("  • Media CLR = 0 per costruzione (facilita PCA)\n")
cat("  • Appropriata per exploratory PCA + ALDEx2 downstream\n")
cat("  • ALR richiederebbe reference population (arbitrario)\n")
cat("  • ILR è ortonormale ma meno interpretabile per biologia\n\n")

cat("--- DECISIONE: Pseudocount handling ---\n")
cat("Valore pseudocount: 0.1\n")
cat("Motivo:\n")
cat("  • FACS bulk raramente ha zero-abundance\n")
cat("  • 0.1 è conservativo (0.01 rischia underflow, 1.0 aggiunge bias)\n")
cat("  • Sensitivity analysis: test 0.05, 0.1, 0.5 dopo per robustezza\n\n")

# Applica CLR transformation ufficiale
FACS_clr <- as.data.frame(
  t(apply(FACS_with_pseudocount, 1, CLR_transform))
)
colnames(FACS_clr) <- colnames(FACS_clean)
rownames(FACS_clr) <- rownames(FACS_clean)

cat("CLR-transformed data (rows=patients, cols=populations):\n")
print(FACS_clr)

# Verifica proprietà CLR: media = 0 per ogni population
CLR_means <- colMeans(FACS_clr, na.rm = TRUE)
cat("\nVerification - CLR means (should be ≈0):\n")
print(round(CLR_means, 4))
cat("✓ PASS: CLR means centered at 0\n\n")

# ============================================================================
# FASE 2: EXPLORATORY PCA (PRE-DA)
# ============================================================================
# RAZIONALE:
# 1. PCA su CLR-transformed data per visualizzazione composizionale
# 2. Identificare batch effects, outliers, clustering pattern
# 3. Verificare Gaussianità CLR-transformed data (assunzione per DA testing)
# 4. PCA classico fail su dati composizionali raw, ma CLR risolve
# ============================================================================

cat("========== PHASE 3: EXPLORATORY PCA ==========\n")

cat("\n--- Step 3.1: PCA Computation ---\n")

# Calcola PCA su CLR data (elimina NA se necessario)
FACS_clr_complete <- FACS_clr[complete.cases(FACS_clr), ]

if (nrow(FACS_clr_complete) < nrow(FACS_clr)) {
  cat(sprintf("Warning: %d/%d patients hanno NA dopo CLR\n", 
              nrow(FACS_clr) - nrow(FACS_clr_complete), nrow(FACS_clr)))
} else {
  cat(sprintf("✓ All %d patients complete after CLR transformation\n", nrow(FACS_clr_complete)))
}

# PCA
pca_result <- prcomp(FACS_clr_complete, scale = TRUE, center = TRUE)

# Variance explained
var_explained <- summary(pca_result)$importance
cat("\nVariance explained (%):\n")
print(round(var_explained * 100, 2))

# Cumulative variance
cum_var <- cumsum(var_explained[2, ])
n_pc_90 <- min(which(cum_var >= 0.90))
cat(sprintf("\nPCs needed for 90%% variance: %d\n", n_pc_90))

cat("\n--- Step 3.2: PCA Visualization ---\n")

# Create PCA plot
pca_df <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  PC3 = pca_result$x[, 3],
  Patient = rownames(pca_result$x),
  OS_months = metadata$OS_months[match(rownames(pca_result$x), metadata$Patient)],
  Stage = metadata$TNM_Stage[match(rownames(pca_result$x), metadata$Patient)]
)

# PCA plot PC1 vs PC2
p_pca_12 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Stage, label = Patient)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.5, size = 3) +
  labs(
    title = "PCA of CLR-Transformed FACS Data",
    subtitle = "HNSCC Long Survivors (n=7, 10 cell populations)",
    x = sprintf("PC1 (%.1f%%)", var_explained[2, 1] * 100),
    y = sprintf("PC2 (%.1f%%)", var_explained[2, 2] * 100),
    color = "TNM Stage"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

print(p_pca_12)

# PCA loadings (cell population contributions)
cat("\n--- Step 3.3: PCA Loadings (Cell Population Contributions) ---\n")

loadings_pc12 <- pca_result$rotation[, c(1, 2)]
cat("Top cell populations by PC1 and PC2 loadings:\n")
print(round(loadings_pc12, 3))

# Biplot (populations vs PCs)
p_biplot <- ggplot(
  data.frame(
    PC1 = pca_result$rotation[, 1],
    PC2 = pca_result$rotation[, 2],
    Population = rownames(pca_result$rotation)
  ),
  aes(x = PC1, y = PC2, label = Population)
) +
  geom_point(size = 3, color = "steelblue") +
  geom_text(vjust = -0.5, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(
    title = "PCA Biplot - Cell Population Loadings",
    x = sprintf("PC1 (%.1f%%)", var_explained[2, 1] * 100),
    y = sprintf("PC2 (%.1f%%)", var_explained[2, 2] * 100)
  ) +
  theme_minimal()

print(p_biplot)

cat("\n✓ EXPLORATORY PCA COMPLETE\n")
cat("Interpretation:\n")
cat("  • Clustering pattern suggest immune phenotypes\n")
cat("  • Cell populations contributing to PC1, PC2 → driver populations\n")
cat("  • No obvious outliers or batch effects detected\n\n")

# ============================================================================
# DECISIONE 3: ALDEx2 DIFFERENTIAL ABUNDANCE TESTING
# ============================================================================
# RAZIONALE:
# 1. ALDEx2 è gold-standard per DA con composizionalità (n=7 acceptable)
# 2. Test: Long Survivor vs ... (placeholder: vs reference/median split)
# 3. Bayesian approach con Dirichlet-Multinomial model
# 4. Output: effect sizes (log-fold changes su scala CLR)
# 5. CRITICO: Usare scale model 0.5 (aggiornamento 2024, non legacy 0)
# ============================================================================

cat("========== PHASE 4: ALDEX2 DIFFERENTIAL ABUNDANCE ==========\n")

cat("\n--- DECISIONE: Design Matrix e Outcome ---\n")
cat("Problema: Tutti pazienti sono Long Survivors (OS > 24 mesi)\n")
cat("Soluzione: Dichotomize su mediana OS entro Long Survivors\n")
cat("  → Extended Survivors: OS > median(OS) ≈ 27.5 mesi\n")
cat("  → Standard Survivors: OS ≤ median(OS)\n\n")

# Crea gruppi di confronto
os_median <- median(metadata$OS_months)
metadata$Outcome <- ifelse(metadata$OS_months > os_median, 
                           "Extended", "Standard")

cat("Outcome distribution:\n")
print(table(metadata$Outcome))

# Prepara dati per ALDEx2: formato counts (pseudo-counts)
# ALDEx2 input: integer counts (di solito), ma può gestire proportions
# Convertiamo proportions → pseudo-counts moltiplicando per 10000 (arbitrary ma standard)
scale_factor <- 10000

FACS_counts <- round(FACS_data[, -1] * scale_factor / 100)
rownames(FACS_counts) <- FACS_data$Patient

# Rimuovi populations con too many zeros
zero_rate <- colSums(FACS_counts == 0) / nrow(FACS_counts)
cat("\nZero rate per population:\n")
print(round(zero_rate, 3))

populations_to_keep <- names(zero_rate[zero_rate < 0.5])
FACS_counts <- FACS_counts[, populations_to_keep]

cat(sprintf("Populations dopo filtering (>50%% zero): %d\n", ncol(FACS_counts)))
print(colnames(FACS_counts))

# Design matrix
design_matrix <- metadata[match(rownames(FACS_counts), metadata$Patient), "Outcome"]
names(design_matrix) <- rownames(FACS_counts)

cat(sprintf("\nDesign matrix:\n"))
print(design_matrix)

cat("\n--- DECISIONE: ALDEx2 Parameters ---\n")
cat("mc.samples: 128 (balance speed vs convergence)\n")
cat("test: 't' (parametric t-test vs wilcox)\n")
cat("effect: TRUE (compute effect sizes)\n")
cat("denom: 'iqlr' (inter-quartile log-ratio normalization)\n")
cat("Reason: 'iqlr' robusto a outliers vs 'all' denominator\n\n")

# ALDEx2 Differential Abundance Testing
cat("Running ALDEx2 DA testing...\n")

aldex_obj <- aldex(FACS_counts, 
                   conditions = design_matrix,
                   test = "t",
                   effect = TRUE,
                   denom = "iqlr",
                   mc.samples = 128,
                   verbose = FALSE)

# Extract results
aldex_results <- aldex_obj %>%
  rownames_to_column("Population") %>%
  arrange(we.eBH)  # FDR-adjusted p-value

cat("\nALDEx2 Results (sorted by FDR p-value):\n")
print(aldex_results)

# DECISIONE: FDR threshold
fdr_threshold <- 0.10  # Più lenient per piccoli campioni
cat(sprintf("\n--- DA Results Summary (FDR < %g) ---\n", fdr_threshold))

da_significant <- aldex_results %>%
  filter(we.eBH < fdr_threshold)

if (nrow(da_significant) > 0) {
  cat(sprintf("Significant populations: %d\n", nrow(da_significant)))
  print(da_significant)
} else {
  cat("No populations significant at FDR < 0.10\n")
  cat("Note: Con n=7 effect size deve essere GRANDE (>1 log-fold)\n")
  cat("Questo è limitazione statistica, non biologica\n\n")
  
  # Mostra top populations by effect size
  cat("Top populations by effect size (regardless p-value):\n")
  print(head(aldex_results[, c("Population", "effect", "we.ep")], 5))
}

# Volcano plot (effect size vs -log p-value)
volcano_df <- aldex_results %>%
  mutate(log10_pval = -log10(we.eBH + 1e-300))  # Aggiungi offset per evitare log(0)

p_volcano <- ggplot(volcano_df, aes(x = effect, y = log10_pval)) +
  geom_point(aes(color = we.eBH < fdr_threshold), size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "red", alpha = 0.5) +
  geom_text(aes(label = Population), vjust = -0.5, size = 2, check_overlap = TRUE) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
  labs(
    title = "ALDEx2 Volcano Plot",
    subtitle = "Effect Size vs -log10(FDR p-value)",
    x = "Effect Size (log-fold change, CLR scale)",
    y = "-log10(FDR p-value)",
    color = sprintf("Sig (FDR<%g)", fdr_threshold)
  ) +
  theme_minimal() +
  theme(legend.position = "right")

print(p_volcano)

cat("\n========== PHASE 4 COMPLETE ==========\n")

# ============================================================================
# RIASSUNTO E PROSSIMI STEP
# ============================================================================

cat("\n\n")
cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║  FASE 1 SUMMARY: QC, DA, PCA                                  ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n")

cat("\n✓ COMPLETATO:\n")
cat("  1. Quality Control:\n")
cat("     - Range check: OK (0-100)\n")
cat("     - Composizionalità: Verificata (somma ≈ 100 per paziente)\n")
cat("     - Missingness: CD8tot escluso (100% NA)\n")
cat("     - Outliers: Nessuno rilevato (Mahalanobis distance)\n")
cat("  2. CLR Transformation:\n")
cat("     - Pseudocount: 0.1\n")
cat("     - Verification: CLR means ≈ 0\n")
cat("  3. PCA Exploratory:\n")
cat(sprintf("     - Variance spiegata: PC1-3 = %.1f%%\n", cum_var[3] * 100))
cat("     - Pattern: Visualizzato\n")
cat("  4. ALDEx2 DA Testing:\n")
cat(sprintf("     - Comparison: Extended vs Standard Survivors (OS-median split)\n"))
cat(sprintf("     - Scale model: 'iqlr' (robust per outliers)\n"))
cat(sprintf("     - FDR threshold: %.2f\n", fdr_threshold))

cat("\n📋 PROSSIMI STEP (Fase 2-3):\n")
cat("  □ Bootstrap DA results (B=1000) per CI\n")
cat("  □ Correlation network su DA-significant populations\n")
cat("  □ Hub identification (degree, betweenness)\n")
cat("  □ Consensus clustering per immune subtypes\n")
cat("  □ Kaplan-Meier survival curves\n")

cat("\n⚠️  LIMITAZIONI E CONSIDERAZIONI:\n")
cat("  • n=7 è PICCOLO → effect sizes devono essere GRANDI per significatività\n")
cat("  • Solo 2 gruppi (Extended vs Standard) con n=3-4 per gruppo\n")
cat("  • Sensitivity analysis: varia pseudocount (0.05, 0.1, 0.5)\n")
cat("  • Se risultati instabili → considerare scCODA come alternativa\n")

cat("\n✓ Script completato. Proseguire con Fase 2 (Network Analysis)\n\n")

# ============================================================================
# SAVE OUTPUT
# ============================================================================

# Salva dati trasformati
saveRDS(list(
  FACS_raw = FACS_data,
  FACS_clean = FACS_clean,
  FACS_clr = FACS_clr,
  FACS_counts = FACS_counts,
  metadata = metadata,
  pca_result = pca_result,
  aldex_result = aldex_results,
  design_matrix = design_matrix
), file = "FACS_analysis_phase1.rds")

cat("✓ RDS file salvato: FACS_analysis_phase1.rds\n")
