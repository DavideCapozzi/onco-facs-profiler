# test/test_logit_behavior.R
# ==============================================================================
# DIAGNOSTIC MODULE: LOGIT TRANSFORMATION CHECK
# Description: Validates the mathematical behavior of the Logit function used
#              in the pipeline and visualizes patient data distribution on the curve.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(gridExtra)
  library(reshape2)
})

# Load Project Modules
source("R/utils_io.R")
source("R/modules_coda.R")

message("\n=== DIAGNOSTIC: LOGIT BEHAVIOR & DATA MAPPING ===")

# 1. Setup & Configuration
# ------------------------------------------------------------------------------
config <- load_config("config/global_params.yml")
output_dir <- file.path(config$output_root, "diagnostics")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 2. Load Real Patient Data (from Step 01)
# ------------------------------------------------------------------------------
data_path <- file.path(config$output_root, "01_data_processing", "data_processed.rds")

if (!file.exists(data_path)) {
  stop("[Error] Processed data not found. Please run src/01_data_processing.R first.")
}

DATA <- readRDS(data_path)
mat_raw <- DATA$raw_matrix

message(sprintf("[Data] Loaded %d samples x %d markers.", nrow(mat_raw), ncol(mat_raw)))

# 3. Mathematical Stress Test
# ------------------------------------------------------------------------------
message("[Check] Running mathematical stress tests on 'coda_transform_logit'...")

# Define edge cases: Exact 0, Exact 100, infinitesimal, normal range
test_vec <- c(0, 0.0001, 0.1, 1, 50, 99, 99.9, 99.9999, 100)
test_mat <- matrix(test_vec, ncol = 1)
colnames(test_mat) <- "Test_Value"

# Run transformation using the pipeline's exact function
# We use a small epsilon to verify clamping works (default in module is usually 1e-6)
trans_mat <- tryCatch({
  coda_transform_logit(test_mat, input_type = "percentage")
}, error = function(e) stop("Transformation function failed: ", e$message))

# Check for Infinite or NaN values
n_inf <- sum(is.infinite(trans_mat))
n_nan <- sum(is.nan(trans_mat))

if (n_inf > 0 || n_nan > 0) {
  stop(sprintf("[FAIL] Transformation generated %d Infinite and %d NaN values.", n_inf, n_nan))
} else {
  message("   -> [PASS] No Infinite or NaN artifacts generated at boundaries.")
}

# Check Monotonicity (Order preservation)
if (all(order(test_mat) == order(trans_mat))) {
  message("   -> [PASS] Monotonicity preserved (Order of values is unchanged).")
} else {
  warning("   -> [WARN] Monotonicity violated. Check epsilon handling.")
}

# 4. Generate Visualization Data
# ------------------------------------------------------------------------------

# A. Theoretical Curve Data
# Generate a dense sequence from 0 to 100 to draw the smooth curve
x_theo <- seq(0, 100, length.out = 2000)
mat_theo <- matrix(x_theo, ncol = 1)
y_theo <- coda_transform_logit(mat_theo, input_type = "percentage")

df_curve <- data.frame(
  Raw_Percentage = x_theo,
  Logit_Value = as.vector(y_theo),
  Type = "Theoretical Function"
)

# B. Real Data Mapping
# Flatten the user matrix to plot individual data points
# We sub-sample if data is too large to keep plot rendering efficient
set.seed(123)
n_points <- min(10000, length(mat_raw))
sample_indices <- sample(length(mat_raw), n_points)

raw_sample <- as.matrix(mat_raw)[sample_indices]
# Transform sample using exactly the same logic
mat_sample <- matrix(raw_sample, ncol = 1)
logit_sample <- coda_transform_logit(mat_sample, input_type = "percentage")

df_points <- data.frame(
  Raw_Percentage = raw_sample,
  Logit_Value = as.vector(logit_sample),
  Type = "Observed Data"
)

# 5. Plotting
# ------------------------------------------------------------------------------
message("[Plot] Generating diagnostic plots...")

# PLOT 1: The Theoretical Curve (Full Range 0-100)
p1 <- ggplot(df_curve, aes(x = Raw_Percentage, y = Logit_Value)) +
  geom_line(color = "#2c3e50", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 50, linetype = "dashed", color = "gray60") +
  annotate("text", x = 10, y = 3, label = "Expansion\nRegion", color = "firebrick", size=3) +
  annotate("text", x = 90, y = -3, label = "Expansion\nRegion", color = "firebrick", size=3) +
  annotate("text", x = 50, y = 0.5, label = "Linear-like\nRegion", size=3) +
  labs(
    title = "1. Theoretical Logit Transformation Curve",
    subtitle = "Visualizing how Input % (X) maps to Logit Space (Y)",
    x = "Input Percentage (0-100)",
    y = "Logit Value (Real Space)"
  ) +
  theme_bw()

# PLOT 2: The Tail Behavior (Zoom in 0-1%)
# Focus on where the 'magic' happens for rare cells
df_curve_zoom <- df_curve %>% filter(Raw_Percentage <= 1)
p2 <- ggplot(df_curve_zoom, aes(x = Raw_Percentage, y = Logit_Value)) +
  geom_line(color = "firebrick", size = 1) +
  labs(
    title = "2. Zoom: Rare Population Expansion (0-1%)",
    subtitle = "Note the steep slope: small raw changes create large logit diffs",
    x = "Input Percentage (0-1%)",
    y = "Logit Value"
  ) +
  theme_bw()

# PLOT 3: User Data Overlay
# Overlay real data points on the curve to see distribution
p3 <- ggplot() +
  # Background Curve
  geom_line(data = df_curve, aes(x = Raw_Percentage, y = Logit_Value), 
            color = "gray70", alpha = 0.5, size = 0.8) +
  # User Data Points
  geom_point(data = df_points, aes(x = Raw_Percentage, y = Logit_Value), 
             color = "#2980b9", alpha = 0.4, size = 1.5, shape = 16) +
  labs(
    title = "3. Your Data Distribution on the Curve",
    subtitle = "Blue dots = Your actual data samples overlaid on the theoretical function",
    x = "Raw Percentage",
    y = "Logit Value"
  ) +
  theme_bw()

# PLOT 4: Density Comparison (Raw vs Logit)
# Shows how the distribution shape changes
p4_raw <- ggplot(df_points, aes(x = Raw_Percentage)) +
  geom_histogram(fill = "gray40", bins = 50) +
  labs(title = "Distribution: Raw %", y = "Count", x = NULL) + theme_minimal()

p4_logit <- ggplot(df_points, aes(x = Logit_Value)) +
  geom_histogram(fill = "#2980b9", bins = 50) +
  labs(title = "Distribution: Logit Transformed", y = "Count", x = NULL) + theme_minimal()

p4 <- gridExtra::arrangeGrob(p4_raw, p4_logit, ncol = 2)

# 6. Save Report
# ------------------------------------------------------------------------------
pdf_file <- file.path(output_dir, "Diagnostic_Logit_Transform.pdf")
pdf(pdf_file, width = 10, height = 12)
grid.arrange(p1, p2, p3, p4, ncol = 1, heights = c(1, 1, 1.2, 1))
dev.off()

message(sprintf("[Output] Diagnostic report saved to: %s", pdf_file))
message("=== DIAGNOSTIC COMPLETE ===\n")