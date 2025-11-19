#!/bin/bash
# Setup script per RStudio Server con ambiente micromamba
# Risolve il problema CXXABI_1.3.15 configurando systemd

set -e  # Exit on error

CONDA_ENV_NAME="facs-comp-analysis"
CONDA_PREFIX="$HOME/micromamba/envs/$CONDA_ENV_NAME"

echo "=== Setup RStudio Server per ambiente $CONDA_ENV_NAME ==="

# === STEP 1: Verifica che l'ambiente esista ===
if [ ! -d "$CONDA_PREFIX" ]; then
    echo "❌ Errore: Ambiente $CONDA_ENV_NAME non trovato"
    exit 1
fi

echo "✓ Ambiente trovato: $CONDA_PREFIX"

# === STEP 2: Verifica versione libstdc++ nell'ambiente ===
LIBSTDCXX_PATH="$CONDA_PREFIX/lib/libstdc++.so.6"
if [ -f "$LIBSTDCXX_PATH" ]; then
    echo "✓ libstdc++.so.6 trovato in conda"
    strings "$LIBSTDCXX_PATH" | grep "CXXABI_1.3.15" && echo "✓ CXXABI_1.3.15 disponibile" || echo "⚠️  CXXABI_1.3.15 NON trovato"
else
    echo "❌ Errore: libstdc++ non trovato nell'ambiente"
    exit 1
fi

# === STEP 3: Crea file di configurazione RStudio Server systemd override ===
SYSTEMD_OVERRIDE_DIR="/etc/systemd/system/rstudio-server.service.d"
SYSTEMD_OVERRIDE_FILE="$SYSTEMD_OVERRIDE_DIR/conda-env.conf"

echo ""
echo "=== Configurazione RStudio Server systemd ==="
echo "Questo richiede privilegi sudo per modificare il servizio systemd"
echo ""

# Crea la directory override
sudo mkdir -p "$SYSTEMD_OVERRIDE_DIR"

# Crea il file di configurazione override
sudo tee "$SYSTEMD_OVERRIDE_FILE" > /dev/null <<EOF
[Service]
# Forza l'uso delle librerie conda
Environment="LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$CONDA_PREFIX/lib/R/lib"
Environment="R_HOME=$CONDA_PREFIX/lib/R"
Environment="R_LIBS_USER=$CONDA_PREFIX/lib/R/library"

# Path aggiornato per trovare binari R
Environment="PATH=$CONDA_PREFIX/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

# Configurazioni aggiuntive per debugging
Environment="R_LD_LIBRARY_PATH=$CONDA_PREFIX/lib"
EOF

echo "✓ File override creato: $SYSTEMD_OVERRIDE_FILE"

# === STEP 4: Reload systemd e restart RStudio Server ===
echo ""
echo "=== Riavvio RStudio Server ==="
sudo systemctl daemon-reload
sudo systemctl restart rstudio-server
sudo systemctl status rstudio-server --no-pager || true

echo ""
echo "✓ RStudio Server riconfigurato"

# === STEP 5: Verifica configurazione ===
echo ""
echo "=== Test di verifica ==="
echo "Eseguendo test R per verificare il caricamento librerie..."

# Crea script R temporaneo per test
TEST_SCRIPT=$(mktemp --suffix=.R)
cat > "$TEST_SCRIPT" <<'REOF'
# Test caricamento librerie
cat("=== Verifica ambiente R ===\n")
cat("R.home():", R.home(), "\n")
cat(".libPaths():", .libPaths(), "\n")

# Test caricamento pacchetti critici
cat("\n=== Test caricamento pacchetti ===\n")
tryCatch({
  library(latticeExtra)
  cat("✓ latticeExtra caricato correttamente\n")
}, error = function(e) {
  cat("❌ Errore latticeExtra:", conditionMessage(e), "\n")
})

tryCatch({
  library(ALDEx2)
  cat("✓ ALDEx2 caricato correttamente\n")
}, error = function(e) {
  cat("❌ Errore ALDEx2:", conditionMessage(e), "\n")
})

# Verifica librerie dinamiche
cat("\n=== Librerie dinamiche caricate ===\n")
loaded <- getLoadedDLLs()
if ("interp" %in% names(loaded)) {
  cat("✓ interp.so caricato correttamente\n")
  cat("  Path:", loaded$interp[["path"]], "\n")
}

cat("\n=== Test completato ===\n")
REOF

# Esegui test con R dell'ambiente
"$CONDA_PREFIX/bin/R" --vanilla --slave < "$TEST_SCRIPT"
rm "$TEST_SCRIPT"

echo ""
echo "=== SETUP COMPLETATO ==="
echo ""
echo "📋 PROSSIMI PASSI:"
echo "1. Accedi a RStudio Server via browser"
echo "2. In una sessione R, verifica con:"
echo "   > library(ALDEx2)"
echo "   > sessionInfo()"
echo ""
echo "3. Se vedi ancora errori, esegui in R:"
echo "   > Sys.getenv('LD_LIBRARY_PATH')"
echo "   > .libPaths()"
echo ""
echo "🔧 TROUBLESHOOTING:"
echo "- Log RStudio Server: sudo journalctl -u rstudio-server -n 50"
echo "- Verifica override: systemctl cat rstudio-server"
echo "- Test manuale R: $CONDA_PREFIX/bin/R"
