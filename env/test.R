# File: test_fixed.R
library(ALDEx2)
library(compositions)

# 1. Simulazione di dati di CONTEGGIO GREZZI
# Simula dati a virgola mobile
data_raw_float <- matrix(runif(100), nrow=10, ncol=10)
# Moltiplica e arrotonda per simulare CONTEGGI INTERI
# Il fattore 1000 simula un "profondità di sequenziamento"
data_counts <- round(data_raw_float * 1000)
print("Matrice di Conteggi Interi:")
print(data_counts)

# --- Test CLR transformation (con il pacchetto compositions) ---
# Se si vuole usare il clr da 'compositions', si usano le proporzioni o i conteggi
data_proportions <- data_counts / rowSums(data_counts) # Normalizza a proporzioni
print("Matrice di Proporzioni:")
print(data_proportions)
clr_test <- clr(data_proportions)
print("Risultato CLR (da compositions):")
print(clr_test)
print("CLR transformation successful!")

# --- Test ALDEx2 ---
# ALDEx2 deve usare i dati di CONTEGGIO GREZZI (data_counts)
# e internamente eseguirà il suo proprio CLR e procedura Monte Carlo.
print("Inizio ALDEx2...")
aldex_test <- aldex(data_counts, rep(c("A","B"), each=5), mc.samples=128)
print("Risultato ALDEx2 (parziale):")
print(head(aldex_test))
print("ALDEx2 funzionale!")