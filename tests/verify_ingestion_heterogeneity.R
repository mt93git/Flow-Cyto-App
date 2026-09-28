
# tests/verify_ingestion_heterogeneity.R
# Purpose: Verify that the refactored 'parse_legacy_facs' handles both Baseline and Secondary datasets correctly.

suppressPackageStartupMessages({
    library(readxl)
    library(stringr)
    library(dplyr)
    library(testthat)
})

# Source the module file to access the internal function
path_mod <- if(file.exists("R/modules/mod_data_loader.R")) "R/modules/mod_data_loader.R" else "../R/modules/mod_data_loader.R"
source(path_mod)

test_that("Baseline Data Ingestion", {
    cat("\n--- Testing Baseline Data ---\n")
    path_baseline <- "../demo_data/demo_facs_data.xls"
    if(!file.exists(path_baseline)) path_baseline <- "demo_data/demo_facs_data.xls"
    if(!file.exists(path_baseline)) skip("Demo FACS file not found")
    
    # Run Parser
    res <- parse_legacy_facs(path_baseline)
    
    # Assertions
    expect_true(nrow(res$mfi) > 0)
    expect_true(nrow(res$counts) > 0)
    
    # Check Sample Pattern (legacy format check)
    # Expected: "BM_Neu_001.fcs" -> ID 1
    sample_ex <- head(res$mfi$Sample, 1)
    id_ex <- head(res$mfi$SampleID, 1)
    
    cat("Baseline Example Sample:", sample_ex, "-> ID:", id_ex, "\n")
    expect_true(!is.na(id_ex))
    
    # Check Population Names
    pops <- unique(res$counts$Population)
    cat("Baseline Pops:", paste(head(pops, 5), collapse=", "), "\n")
    expect_true("CD45+" %in% pops)
})

test_that("Dynamic Population Remapping & Polymorphism", {
    cat("\n--- Testing Population Remapping on Standard Benchmark ---\n")
    path_data <- "demo_data/demo_facs_data.xls"
    if(!file.exists(path_data)) path_data <- "../demo_data/demo_facs_data.xls"
    if(!file.exists(path_data)) skip("Demo FACS file not found")
    
    # Polymorphic Map: Rename 'Cells' to 'Total_Viable_Cells'
    poly_map <- list("Cells" = "Total_Viable_Cells")
    
    # Run Parser with Mapping
    res <- parse_legacy_facs(path_data, population_mapping = poly_map)
    
    # Assertions
    expect_true(nrow(res$mfi) > 0)
    expect_true(nrow(res$counts) > 0)
    
    pops <- unique(res$counts$Population)
    cat("Mapped Pops Preview:", paste(head(pops, 5), collapse=", "), "\n")
    
    # Check that 'Total_Viable_Cells' is present and 'Cells' is replaced
    expect_true("Total_Viable_Cells" %in% pops)
    expect_false("Cells" %in% pops)
})

cat("\nVerification Complete\n")
