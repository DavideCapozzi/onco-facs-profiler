# ==============================================================================
# IMMUNOLOGY DATA PIPELINE: CLEAN, STANDARDIZE, ANONYMIZE
# ==============================================================================

library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

# ==============================================================================
# 1. CONFIGURATION
# ==============================================================================

# Percorsi file (Verifica che siano corretti per il tuo PC)
file_path   <- "C:/Users/Davide/Desktop/Davide/bandi/0dottorato/projects/long survivors/DB_pulito.xlsx"
output_path <- "C:/Users/Davide/Desktop/Davide/bandi/0dottorato/projects/long survivors/DB_anonimo_standardizzato.xlsx"

# Configurazione Fogli
sheet_config <- list(
  "HNSCC OS>24"    = list(type = "standard",      out_name = "HNSCC", include = TRUE),
  "NSCLC OS > 18"  = list(type = "standard",      out_name = "NSCLC", include = TRUE),
  "GB OS>18"       = list(type = "double_header", out_name = "GB",    include = FALSE), 
  "HD"             = list(type = "simple_hd",     out_name = "Healthy_Donors", include = TRUE)
)

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

get_id_prefix <- function(sheet_name) {
  str_extract(sheet_name, "^\\S+") 
}

clean_name_keep_hyphen <- function(x) {
  if (is.na(x)) return("")
  x %>%
    str_trim() %>%
    str_replace_all("\\s+", "-") %>%
    str_replace_all("[^A-Za-z0-9-]", "") %>%
    str_to_upper()
}

get_bio_root <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("[^a-z0-9]", "") %>%
    str_remove_all("(tot|total|abs|absolute|count|freq|val|value)$")
}

# ==============================================================================
# 3. PARSING LOGIC
# ==============================================================================

process_standard <- function(path, sheet) {
  df <- read.xlsx(path, sheet = sheet)
  names(df) <- map_chr(names(df), clean_name_keep_hyphen)
  
  prefix <- get_id_prefix(sheet)
  df <- df %>%
    rename(Patient_ID = 1) %>%
    mutate(Patient_ID = paste0(prefix, row_number()))
  return(df)
}

process_double_header <- function(path, sheet) {
  raw_headers <- read.xlsx(path, sheet = sheet, rows = 1:2, colNames = FALSE)
  h1 <- as.character(raw_headers[1, ])
  h2 <- as.character(raw_headers[2, ])
  
  filled_h1 <- tibble(val = h1) %>% fill(val, .direction = "down") %>% pull(val)
  
  new_names <- map2_chr(filled_h1, h2, function(top, sub) {
    c_top <- clean_name_keep_hyphen(top)
    c_sub <- clean_name_keep_hyphen(sub)
    if (c_sub == "") return("REMOVE")
    if (c_top == "" || c_top == "NA") return(c_sub)
    return(paste(c_sub, c_top, sep = "-")) 
  })
  
  df <- read.xlsx(path, sheet = sheet, startRow = 3, colNames = FALSE)
  names(df) <- new_names
  prefix <- get_id_prefix(sheet)
  
  df <- df %>%
    select(-any_of("REMOVE")) %>%
    rename(Patient_ID = 1) %>%
    mutate(Patient_ID = paste0(prefix, row_number()))
  return(df)
}

process_simple_hd <- function(path, sheet) {
  raw <- read.xlsx(path, sheet = sheet, colNames = FALSE)
  headers <- raw[1, ] %>% 
    unlist() %>%
    str_remove_all("\\.fcs$") %>% 
    map_chr(clean_name_keep_hyphen)
  
  df <- raw[-1, ]
  names(df) <- headers
  prefix <- get_id_prefix(sheet)
  
  df <- type.convert(df, as.is = TRUE) %>%
    rename(Patient_ID = 1) %>%
    mutate(Patient_ID = paste0(prefix, row_number()))
  return(df)
}

# ==============================================================================
# 4. EXECUTION: LOAD & PARSE
# ==============================================================================

processed_data <- list()
original_sheets <- getSheetNames(file_path)

message("--- 1. READING DATA ---")

for (orig_name in names(sheet_config)) {
  cfg <- sheet_config[[orig_name]]
  if (!cfg$include) next
  
  if (orig_name %in% original_sheets) {
    if (cfg$type == "standard") {
      processed_data[[cfg$out_name]] <- process_standard(file_path, orig_name)
    } else if (cfg$type == "double_header") {
      processed_data[[cfg$out_name]] <- process_double_header(file_path, orig_name)
    } else if (cfg$type == "simple_hd") {
      processed_data[[cfg$out_name]] <- process_simple_hd(file_path, orig_name)
    }
  } else {
    warning(paste("Sheet not found:", orig_name))
  }
}

# ==============================================================================
# 5. ANALYSIS & STANDARDIZATION
# ==============================================================================

message("--- 2. ANALYZING COLUMN SIMILARITIES ---")

all_cols <- map_dfr(names(processed_data), function(s) {
  tibble(Sheet = s, Current_Name = names(processed_data[[s]]))
}) %>% filter(Current_Name != "Patient_ID")

renaming_logic <- all_cols %>%
  mutate(
    Bio_Root = map_chr(Current_Name, get_bio_root),
    Clean_Base = Current_Name 
  ) %>%
  group_by(Bio_Root) %>%
  mutate(
    # AUTOMATISMO: Sceglie il nome più corto (es. vince "CD8" su "CD8TOT")
    Final_Proposal = Clean_Base[which.min(nchar(Clean_Base))] 
  ) %>%
  ungroup() 

# Log Conflicts
conflicts_view <- renaming_logic %>%
  group_by(Bio_Root) %>%
  filter(n_distinct(Current_Name) > 1) %>%
  ungroup() %>%
  select(Bio_Root, Sheet, Current_Name, Final_Proposal) %>%
  pivot_wider(names_from = Sheet, values_from = Current_Name) %>%
  relocate(Final_Proposal, .after = last_col())

if (nrow(conflicts_view) > 0) {
  message("--- CONFLICTS DETECTED (Will be auto-corrected to Final_Proposal) ---")
  print(as.data.frame(conflicts_view))
} else {
  message("No naming conflicts found.")
}

# Crea la mappa automatica: VECCHIO_NOME -> NUOVO_NOME
auto_map <- setNames(renaming_logic$Final_Proposal, renaming_logic$Current_Name)

# ==============================================================================
# 6. MANUAL OVERRIDES (OPZIONALE)
# ==============================================================================

message("\n--- 3. APPLYING RENAMING ---")

# Se vuoi forzare un nome diverso da quello automatico, scommenta e aggiungi qui sotto.
# Altrimenti lascia NULL e lo script userà la "Final_Proposal" automatica.

manual_overrides <- c(
  # "VECCHIO_NOME_NEL_DB" = "NOME_CHE_VUOI_TU",
  # Es: "LOX1-PMN" = "LOX1-MDSC",
  NULL 
)

final_map <- auto_map
if (!is.null(manual_overrides)) {
  # Sovrascrive la logica automatica solo per le voci specificate
  final_map[names(manual_overrides)] <- manual_overrides
}

# ==============================================================================
# 7. UPDATE DATAFRAMES (CORRETTO)
# ==============================================================================

final_datasets <- map(processed_data, function(df) {
  
  # 1. Identifica le colonne presenti
  cols_to_rename <- intersect(names(df), names(final_map))
  
  if (length(cols_to_rename) > 0) {
    # 2. Prepara la sotto-mappa
    current_selection <- final_map[cols_to_rename]
    
    # 3. CRUCIALE: Inverte la mappa per dplyr::rename
    # La mappa originale è:   c("CD8TOT" = "CD8")   (VECCHIO = NUOVO)
    # rename vuole:           c("CD8" = "CD8TOT")   (NUOVO = VECCHIO)
    rename_vec <- setNames(names(current_selection), current_selection)
    
    # 4. Esegue il rinominamento
    df <- df %>% rename(any_of(rename_vec))
  }
  return(df)
})

# ==============================================================================
# 8. EXPORT
# ==============================================================================

write.xlsx(final_datasets, file = output_path, overwrite = TRUE)
message(paste("Success! File saved to:", output_path))
