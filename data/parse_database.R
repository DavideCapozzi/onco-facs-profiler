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

file_path   <- "path/to/excel"
output_path <- "output/path/excel"

# Sheet Configuration
# include: TRUE to process, FALSE to skip.
# out_name: The name of the sheet in the final Excel file.
sheet_config <- list(
  #Input sheet name 
  "HNSCC OS>24"    = list(type = "standard",      out_name = "HNSCC", include = TRUE),
  "NSCLC OS > 18"  = list(type = "standard",      out_name = "NSCLC", include = TRUE),
  "GB OS>18"       = list(type = "double_header", out_name = "GB",    include = FALSE), 
  "HD"             = list(type = "simple_hd",     out_name = "Healthy_Donors", include = TRUE)
)

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

#' Generates ID prefix from sheet name.
#' Logic: Takes the first word (everything before the first space).
#' Example: "HNSCC OS>24" -> "HNSCC"
get_id_prefix <- function(sheet_name) {
  str_extract(sheet_name, "^\\S+") 
}

#' Cleans column names while PRESERVING HYPHENS.
#' Converts spaces to hyphens, removes other symbols, converts to Uppercase.
clean_name_keep_hyphen <- function(x) {
  if (is.na(x)) return("")
  x %>%
    str_trim() %>%
    str_replace_all("\\s+", "-") %>%        # Space -> Hyphen
    str_replace_all("[^A-Za-z0-9-]", "") %>% # Keep only alphanumeric and hyphens
    str_to_upper()
}

#' Extract 'Biological Root' for matching logic.
#' Removes separators and technical suffixes to find similarities.
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
  
  # ID Logic: First word of sheet name + row number
  prefix <- get_id_prefix(sheet)
  
  df <- df %>%
    rename(Patient_ID = 1) %>%
    mutate(Patient_ID = paste0(prefix, row_number()))
  
  return(df)
}

process_double_header <- function(path, sheet) {
  # Read headers (first 2 rows)
  raw_headers <- read.xlsx(path, sheet = sheet, rows = 1:2, colNames = FALSE)
  h1 <- as.character(raw_headers[1, ])
  h2 <- as.character(raw_headers[2, ])
  
  # Fill merged cells logic
  filled_h1 <- tibble(val = h1) %>% fill(val, .direction = "down") %>% pull(val)
  
  # Merge Top and Sub headers
  new_names <- map2_chr(filled_h1, h2, function(top, sub) {
    c_top <- clean_name_keep_hyphen(top)
    c_sub <- clean_name_keep_hyphen(sub)
    
    if (c_sub == "") return("REMOVE")
    if (c_top == "" || c_top == "NA") return(c_sub)
    return(paste(c_sub, c_top, sep = "-")) 
  })
  
  # Read Data
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
  
  # Manual header parsing
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

# A. Extract all column names across all processed sheets
all_cols <- map_dfr(names(processed_data), function(s) {
  tibble(Sheet = s, Current_Name = names(processed_data[[s]]))
}) %>% filter(Current_Name != "Patient_ID")

# B. Calculate Proposed Names
renaming_logic <- all_cols %>%
  mutate(
    Bio_Root = map_chr(Current_Name, get_bio_root),
    Clean_Base = Current_Name 
  ) %>%
  group_by(Bio_Root) %>%
  mutate(
    # Heuristic: Choose the name that appears most often, or the shortest one
    Final_Proposal = Clean_Base[which.min(nchar(Clean_Base))] 
  ) %>%
  ungroup() 

# C. CREATE PIVOT VIEW (Wide Format)
# This transforms the data so each Sheet is a column
wide_comparison <- renaming_logic %>%
  select(Bio_Root, Sheet, Current_Name, Final_Proposal) %>%
  pivot_wider(
    names_from = Sheet, 
    values_from = Current_Name
  ) %>%
  # Reorder columns to put Final_Proposal at the end
  relocate(Final_Proposal, .after = last_col()) %>% 
  arrange(Bio_Root)

# Print the FULL comparison table
message("\n--- FULL COLUMN COMPARISON ---")
print(as.data.frame(wide_comparison))

# D. IDENTIFY AND PRINT CONFLICTS ONLY
# Filter logic: Group by Bio_Root and check if there is more than 1 distinct name used.
conflicts_view <- renaming_logic %>%
  group_by(Bio_Root) %>%
  filter(n_distinct(Current_Name) > 1) %>% # Keep only if names differ across sheets
  ungroup() %>%
  select(Bio_Root, Sheet, Current_Name, Final_Proposal) %>%
  pivot_wider(
    names_from = Sheet, 
    values_from = Current_Name
  ) %>%
  relocate(Final_Proposal, .after = last_col())

message("\n--- CONFLICTS DETECTED (Differences between sheets) ---")
# Only print if there are actual conflicts
if (nrow(conflicts_view) > 0) {
  print(as.data.frame(conflicts_view))
} else {
  message("No naming conflicts found across sheets.")
}

# D. Create Automatic Map
auto_map <- setNames(renaming_logic$Final_Proposal, renaming_logic$Current_Name)

# ==============================================================================
# 6. MANUAL OVERRIDES
# ==============================================================================

message("\n--- 3. APPLYING RENAMING ---")

# DEFINE MANUAL CHANGES HERE
# Use the printed table above. If 'Proposed_Name' is not what you want,
# add an entry here: "Current_Name_Seen_In_Table" = "Your_Desired_Name"

manual_overrides <- c(
  # "LOX1-PMN" = "LOX1-PMN-MDSC", 
  # "CD8-TOT"  = "CD8",
  NULL # Placeholder
)

# Apply overrides to the map
final_map <- auto_map
if (!is.null(manual_overrides)) {
  final_map[names(manual_overrides)] <- manual_overrides
}

# Update Dataframes
final_datasets <- map(processed_data, function(df) {
  cols_to_rename <- intersect(names(df), names(final_map))
  if (length(cols_to_rename) > 0) {
    df <- df %>% rename(any_of(final_map))
  }
  return(df)
})

# ==============================================================================
# 7. EXPORT
# ==============================================================================

write.xlsx(final_datasets, file = output_path, overwrite = TRUE)
message(paste("Success! File saved to:", output_path))
