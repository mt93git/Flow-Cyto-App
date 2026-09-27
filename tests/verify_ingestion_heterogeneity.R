
# tests/verify_ingestion_heterogeneity.R
# Purpose: Verify that the refactored 'parse_legacy_facs' handles both Baseline and Secondary datasets correctly.

suppressPackageStartupMessages({
    library(readxl)
    library(stringr)
    library(dplyr)
    library(testthat)
})

# Source the module file to access the internal function
# Note: In a package we'd export it or use :::, but here we source
source("../R/modules/mod_data_loader.R")

test_that("Baseline Data Ingestion", {
    cat("\n--- Testing Baseline Data ---\n")
    path_baseline <- "../28-Oct-2025_exp6_7_LPS_PDAC.xls"
    
    if(!file.exists(path_baseline)) skip("Baseline file not found")
    
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

test_that("Secondary Data Ingestion (Heterogeneous)", {
    cat("\n--- Testing Secondary Data ---\n")
    path_secondary <- "../New files for Heat map data-2/20251105 and 20251107 and 20251113 BoneMarrow FloJo Gating Analysis.xls"
    
    if(!file.exists(path_secondary)) skip("Secondary file not found")
    
    # Polymorphic Map (Duplicate from Server logic)
    poly_map <- list("immune" = "CD45+")
    
    # Run Parser with Mapping
    res <- parse_legacy_facs(path_secondary, population_mapping = poly_map)
    
    # Assertions
    expect_true(nrow(res$mfi) > 0)
    expect_true(nrow(res$counts) > 0)
    
    # Check MFI Capture (was failing before due to negative values/format?)
    # "Geometric Mean : Comp-APC-A = -136.833755493"
    cat("Secondary MFI Rows:", nrow(res$mfi), "\n")
    if(nrow(res$mfi) > 0) {
        cat("Secondary MFI Example:", res$mfi$Parameter[1], "=", res$mfi$Value[1], "\n")
    }
    
    # Check Population Mapping
    # "immune" should be mapped to "CD45+"
    pops <- unique(res$counts$Population)
    cat("Secondary Pops (Mapped):", paste(head(pops, 5), collapse=", "), "\n")
    
    expect_true("CD45+" %in% pops)
    expect_false("immune" %in% pops) # Should be replaced
})

cat("\nVerification Complete\n")
