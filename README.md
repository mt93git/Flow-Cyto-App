# Flow-Cyto App (v4.7)

[![R Version](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![Shiny](https://img.shields.io/badge/Framework-Shiny%20%7C%20bslib-2C3E50.svg)](https://shiny.posit.co/)
[![Bioconductor](https://img.shields.io/badge/Bioc-ComplexHeatmap-brightgreen.svg)](https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **Interactive R/Shiny Platform for Multiparametric Flow Cytometry (FACS) Exploration, Absolute Cell Quantification, and Complex Phenotypic Profiling.**  
> *Developed by Maxence Tricaud (Libreros Lab — Yale University / Université Laval).*

---

## 1. Overview & Biological Purpose

**Flow-Cyto App** is a specialized computational environment engineered to process, normalize, and visualize high-dimensional flow cytometry datasets (conventional and spectral FACS). The platform solves key analytical bottlenecks in cytometry pipelines:
1. **Volumetric Acquisition Bias:** Standardizes cell event numbers into absolute cell counts using exogenous reference counting beads (CountBright™ standard).
2. **Dynamic Gating Hierarchy Ingestion:** Recursively parses multi-tier hierarchical population trees exported from FlowJo™ / BD FACSDiva™ reports.
3. **Variance Stabilization & Cross-Marker Scaling:** Couples cofactor-adjusted arcsinh transformation with sample-level Z-score standardization for high-contrast biomarker discovery.
4. **Hierarchical Multi-Track Visualization:** Leverages a `ComplexHeatmap` visualization engine with real-time reactive reordering by experimental cohort, tissue origin, or activation state.

---

## 2. ⚡ Quick Start: 1-Click Evaluation (Built-in Demo)

The repository includes a curated, fully functional preclinical demonstration dataset (`demo_data/`) comparing murine bone marrow leukocyte activation under LPS challenge vs. steady-state conditions.

### Option A: Launch in RStudio
1. Open `FlowCytometryAnalysisApp.Rproj` in RStudio.
2. Open `app.R` and click **Run App** (or execute `shiny::runApp()`).
3. In the sidebar, simply click the green button:  
   👉 **`📊 Load Built-in Demo Dataset`**
4. The complete interactive Heatmap, Z-scores, Absolute Counts, and data tables will populate dynamically!

### Option B: Terminal Launch
```bash
# Clone the repository
git clone https://github.com/mt93git/Flow-Cyto-App.git
cd Flow-Cyto-App

# Install dependencies (once)
Rscript setup_dependencies.R

# Launch the interactive application
Rscript -e "shiny::runApp(port = 3838, launch.browser = TRUE)"
```

---

## 3. Mathematical & Algorithmic Methods

### A. Absolute Cell Count Normalization Strategy
To convert raw event counts into true volumetric concentrations without volumetric fluidics error, the engine isolates reference microspheres spiked at known concentrations:

$$\text{AbsCount} = \left( \frac{\text{CellCount}_{\text{target}}}{\text{BeadCount}_{\text{ref}}} \right) \times \text{BeadsInput} \times \left( \frac{\text{Vol}_{\text{total}}}{\text{Vol}_{\text{sample}}} \right)$$

* **`CellCount`**: Raw event frequency of the gated biological phenotype.
* **`BeadCount`**: Reference microsphere event frequency, targeted via dynamic POSIX regular expressions (`CountBright|Bead`).
* **`BeadsInput`**: Initial reference bead spike quantity.
* **`VolTotal` / `VolSample`**: Resuspension volume and acquired aliquot volume.

### B. Variance Stabilization & Transformation
Raw Geometric Mean Fluorescence Intensities (MFI, denoted $x$) span multiple orders of magnitude with near-zero baseline noise. Values undergo an inverse hyperbolic sine transformation:

$$f(x) = \operatorname{asinh}\left( \frac{x}{\text{cofactor}} \right), \quad \text{with } \text{cofactor} = 150$$

### C. Cohort Standardization
For comparative cross-panel clustering, transformed values are standardized per marker and cell population across all biological replicates:

$$Z = \frac{f(x) - \mu}{\sigma}$$

---

## 4. Repository Structure

```
Flow-Cyto-App/
├── app.R                       # Main Shiny application entry point
├── FlowCytometryAnalysisApp.Rproj # RStudio Project configuration
├── LICENSE                     # Formal MIT Open-Source License
├── README.md                   # Technical documentation
├── setup_dependencies.R        # Automated CRAN/Bioconductor installer
├── R/
│   ├── global.R                # Path management and global constants
│   └── modules/
│       ├── mod_data_loader.R   # Parsing engine, normalization & demo handler
│       └── mod_heatmap.R       # ComplexHeatmap reactive visualization UI/server
├── demo_data/                  # Built-in demonstration benchmark
│   ├── demo_facs_data.xls      # Multi-population hierarchical FACS report
│   ├── demo_metadata.xlsx      # Biological metadata (treatment, sex, volumes)
│   └── demo_marker_map.csv     # Fluorophore-to-marker channel mapping
└── tests/
    ├── verify_logic.R          # Standalone end-to-end integration test
    └── verify_ingestion_heterogeneity.R # Edge-case hierarchy test suite
```

---

## 5. Input Data Format

When uploading your own datasets, the application expects three complementary files:
1. **FACS Export (`.xls`):** Hierarchical population report containing columns `Name`, `Statistic`, and `#Cells`.
2. **Metadata Table (`.xlsx`):** Must contain `sample_no`, experimental condition (`condition`), and bead volumetric columns (`volume_total`, `volume_sample_Neu`, `beads_input`).
3. **Marker Map (`.csv`):** Mapping table with columns `fluor` (e.g. `APC-A`, `BUV395-A`) and `marker` (e.g. `CD117`, `Sca1`).

---

## 6. License & Citation

Distributed under the **MIT License**. See `LICENSE` for details.

If utilizing **Flow-Cyto App** in academic studies, please cite:
> Tricaud M. *Flow-Cyto App: Interactive Flow Cytometry Analysis and Absolute Count Normalization Suite*. Yale University / Université Laval (2025-2026). GitHub: `https://github.com/mt93git/Flow-Cyto-App`.
