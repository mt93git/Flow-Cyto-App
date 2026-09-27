# setup_dependencies.R

# 1. Install BiocManager if missing
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# 2. Check and Install ComplexHeatmap
if (!require("ComplexHeatmap", quietly = TRUE)) {
    message("Installing ComplexHeatmap from Bioconductor...")
    BiocManager::install("ComplexHeatmap", ask = FALSE)
}

# 3. Check other key dependencies
libs <- c("shiny", "bslib", "readxl", "dplyr", "tidyr", "stringr", "DT", "circlize", "shinyjqui", "colourpicker", "tibble", "sortable")
missing_libs <- libs[!(libs %in% installed.packages()[,"Package"])]
if(length(missing_libs)) {
  message("Installing missing CRAN packages: ", paste(missing_libs, collapse=", "))
  install.packages(missing_libs)
}

message("Setup Complete! You can now run the app.")
