#!/bin/bash
set -e

# === SETUP RENVIRON LOCALE AL PROGETTO ===
# Questo script genera un file .Renviron specifico per questo progetto
# basandosi sull'ambiente micromamba ATTIVO.

# 1. Verifica che l'ambiente sia attivo
if [ -z "$CONDA_PREFIX" ]; then
    echo "❌ Errore: Nessun ambiente Conda/Micromamba attivo."
    echo "👉 Attiva prima l'ambiente con: micromamba activate facs-comp-analysis"
    exit 1
fi

# 2. Definisci il file target (nella root del progetto corrente)
# Usa $(pwd) per assicurarsi che finisca nella cartella da cui lanci lo script (la root del progetto)
PROJECT_ROOT=$(pwd)
PROJECT_RENVIRON="$PROJECT_ROOT/.Renviron"

echo "🔧 Configurazione ambiente R per il progetto in: $PROJECT_ROOT"
echo "   Ambiente sorgente: $CONDA_PREFIX"

# 3. Backup se esiste già
if [ -f "$PROJECT_RENVIRON" ]; then
    echo "⚠️  Trovato .Renviron esistente. Backup in .Renviron.bak"
    cp "$PROJECT_RENVIRON" "${PROJECT_RENVIRON}.bak"
fi

# 4. Generazione del file .Renviron locale
cat <<EOF > "$PROJECT_RENVIRON"
# --- Generato automaticamente da setup/renviron.sh ---
# Progetto: $(basename "$PROJECT_ROOT")
# Ambiente Micromamba: $(basename "$CONDA_PREFIX")

# PATH: Priorità ai binari dell'ambiente
PATH="${CONDA_PREFIX}/bin:\${PATH}"

# LIBRERIE SISTEMA: Fix per CXXABI e compilatori
LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:\${LD_LIBRARY_PATH}"

# R LIBS: Percorso pacchetti specifico dell'ambiente
R_LIBS_USER="${CONDA_PREFIX}/lib/R/library"

# R HOME: Eseguibile R dell'ambiente
R_HOME="${CONDA_PREFIX}/lib/R"
EOF

echo "✓ File .Renviron creato nel progetto."
echo "🚀 Riavvia la sessione R (Session -> Restart R) per applicare."
