# app.R
# Main Entry Point

# Source global config
source("R/global.R")

# Source modules
source("R/modules/mod_data_loader.R")
source("R/modules/mod_heatmap.R")

# Define UI
ui <- page_fluid(
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#2C3E50"),
  
  titlePanel("Flow Cytometry Analysis App v4.7"),
  
  layout_sidebar(
    sidebar = sidebar(
      title = "Data Controls",
      width = 350,
      mod_data_loader_ui("loader")
    ),
    
    # Main Content Area
    navset_card_underline(
      nav_panel("Heatmap", mod_heatmap_ui("heatmap")),
      nav_panel("Data Table", DTOutput("debug_table")),
      nav_panel("Methods", 
                div(class="card", style="margin: 20px; padding: 20px;",
                    markdown("
### **Methods & Protocols**

#### **1. Absolute Count Normalization Strategy**
Absolute cell quantification is derived using exogenous reference standard (CountBright, Absolute Counting Beads) spiked at a known concentration; Allowing to mitigate volumetric variability during acquisition.

**Formula Used:**
AbsCount = (CellCount_target / BeadCount_ref) * BeadsInput * (Vol_total / Vol_sample)

**Parameter Definitions:**
*   **`CellCount`**: Raw event frequency of the target biological population.
*   **`BeadCount`**: Raw event frequency of the reference coordinate, identified via the *Control Regex* mechanism.
*   **`BeadsInput`**: Initial spike-in quantity (Metadata: `beads_input`).
*   **`VolTotal`**: Resuspension volume (Metadata: `volume_total`).
*   **`VolSample`**: Acquisition aliquot volume.

#### **2. Control Population Identification (Regex Mechanism)**
The system employs **Regular Expressions (Regex)** to programmatically segregate reference bead events from biological populations within the input matrix.

*   **Mechanism**: A pattern matching algorithm scans the `Name` column of the input dataset.
*   **Syntax**: The query strings utilize standard POSIX Extended Regular Expression syntax.
    *   **Logical Disjunction (`|`)**: The pipe operator functions as a logical `OR`.
        *   *Usage*: The default pattern `CountBright|Bead` targets any identifier containing the substring \"CountBright\" **OR** \"Bead\".
    *   **Application**: Rows validating TRUE against this pattern are designated as `Reference` and sequestered from downstream Clustering/Heatmap analysis, serving exclusively as the normalization denominator.

#### **3. Data Transformation & Variance Stabilization**
*   **Arcsinh Transformation**: Raw fluorescence intensities (defined as x) undergo an inverse hyperbolic sine transformation: f(x) = asinh(x / 150). This linearizes data in the lower range while matching log behavior in the upper range.
*   **Z-Score Standardization**: For visualization, marker expression is standardized across the sample cohort: Z = (x - mean) / standard_deviation. This centers the distribution at 0 with unit variance.

---
**Author**: Maxence Tricaud for Libreros Lab
**Date**: Dec 18 2025
**Contact**: maxence.benjamin@gmail.com
                    ")
                )
      )
    )
  )
)

# Define Server
server <- function(input, output, session) {
  
  # Module: Data Loader
  # Returns: list(data = reactiveVal, log = reactiveVal)
  loader_res <- mod_data_loader_server("loader")
  
  # Extract Data for Heatmap
  mfi_data <- reactive({
    res <- loader_res$data()
    req(res)
    res$data
  })
  
  # Render Log
  output$system_log <- renderText({
     loader_res$log()
  })
  
  # Module: Heatmap
  # Consumes the processed MFI data
  mod_heatmap_server("heatmap", mfi_data)
  
  # Debug Table
  output$debug_table <- renderDT({
    req(mfi_data())
    datatable(mfi_data(), options = list(scrollX = TRUE, pageLength = 10))
  })
}

shinyApp(ui, server)
