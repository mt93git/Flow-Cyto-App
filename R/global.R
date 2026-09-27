# R/global.R
# Global configuration and library loading

# --- 1. Path Management (Anti-renv Strategy) ---
# Explicitly set libPaths to ensure system libraries are accessible
# This prevents "package not found" errors in restricted environments
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

# --- 2. Load Libraries ---
suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(DT)
  library(ComplexHeatmap)
  library(circlize)
  library(shinyjqui) # For resizable elements
  library(colourpicker)
  library(tibble)
  library(sortable)
})

# --- 3. Global Constants ---
COFACTOR <- 150 # for asinh transformation

# --- 4. Helper Functions ---
# (Will add shared helpers here if needed)
