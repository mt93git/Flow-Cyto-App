# tests/verify_logic.R
# Standalone automated integration test for Flow-Cyto App pipeline

library(readxl)
library(dplyr)
library(stringr)

cat("====================================================\n")
cat("--- 1. Testing Environment & Parser Ingestion ---\n")
cat("====================================================\n")

# Path discovery
base_dir <- getwd()
facs_file <- file.path("..", "demo_data", "demo_facs_data.xls")
meta_file <- file.path("..", "demo_data", "demo_metadata.xlsx")
marker_file <- file.path("..", "demo_data", "demo_marker_map.csv")

if(!file.exists(facs_file)) {
    facs_file <- file.path("demo_data", "demo_facs_data.xls")
    meta_file <- file.path("demo_data", "demo_metadata.xlsx")
    marker_file <- file.path("demo_data", "demo_marker_map.csv")
    source("R/modules/mod_data_loader.R")
} else {
    source("../R/modules/mod_data_loader.R")
}

if(!file.exists(facs_file)) stop("CRITICAL: Demo FACS file not found.")

cat("Loading and executing complete end-to-end flow pipeline...\n")
res <- execute_flow_pipeline(
    facs_path = facs_file,
    meta_path = meta_file,
    marker_path = marker_file,
    bead_pattern = "CountBright|Bead",
    add_log = function(m) cat("  [LOG]", m, "\n")
)

cat("\n====================================================\n")
cat("--- 2. Validating Pipeline Outputs ---\n")
cat("====================================================\n")

stopifnot(!is.null(res$data))
stopifnot(!is.null(res$count_data))
stopifnot(nrow(res$data) > 0)
stopifnot(nrow(res$count_data) > 0)

cat("✓ MFI Processed Rows:", nrow(res$data), "\n")
cat("✓ Absolute Counts Rows:", nrow(res$count_data), "\n")
cat("✓ Transformed Values Range:", round(range(res$data$Value_Transformed, na.rm=TRUE), 3), "\n")
cat("✓ Z-Score Mean (should be ~0):", round(mean(res$data$Z_Score, na.rm=TRUE), 4), "\n")

cat("\n✓ Sample Preview of Normalized Data:\n")
print(head(res$data %>% select(SampleID, marker, condition, Value_Transformed, Z_Score), 5))

cat("\n====================================================\n")
cat(">>> ALL INTEGRATION TESTS PASSED SUCCESSFULLY! <<<\n")
cat("====================================================\n")
