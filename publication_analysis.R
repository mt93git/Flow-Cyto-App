suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ggpubr))

source("R/modules/mod_data_loader.R")

cat("Initializing Data Processing for Publication & Reproducibility...\n")

# -- FILE PATHS (Curated Anonymized Demo Data) --
meta_path <- "./demo_data/demo_metadata.xlsx"
marker_path <- "./demo_data/demo_marker_map.csv"
facs_path <- "./demo_data/demo_facs_data.xls"

if (!file.exists(meta_path) || !file.exists(marker_path) || !file.exists(facs_path)) {
  stop("Curated demo datasets not found in ./demo_data/.")
}

# -- 1. DATA INGESTION --
meta <- read_excel(meta_path)
meta$sample_no <- as.character(meta$sample_no)

markers <- read.csv(marker_path)
parsed <- parse_legacy_facs(facs_path)

# Prepare Output Report File
stat_report <- file("Publication_Statistical_Report.txt", open="wt")
writeLines("Publication Statistical & Phenotypic Quality Report\n===================================================\n", stat_report)

# -- 2. MFI MATRIX PREPARATION & SCALING --
cat("Processing MFI Matrix...\n")
mfi_df <- parsed$mfi
mfi_combined <- mfi_df %>%
  mutate(SampleID = as.character(SampleID)) %>%
  inner_join(meta, by = c("SampleID" = "sample_no")) %>%
  mutate(CleanFluor = trimws(gsub("Comp-", "", Parameter))) %>%
  left_join(markers, by = c("CleanFluor" = "fluor")) %>%
  filter(!is.na(marker)) %>%
  mutate(Value_Transformed = asinh(as.numeric(Value) / 150)) %>%
  group_by(marker, Population) %>%
  mutate(Z_Score = (Value_Transformed - mean(Value_Transformed, na.rm=TRUE)) / sd(Value_Transformed, na.rm=TRUE)) %>%
  ungroup()

mature_mfi <- mfi_combined %>% filter(Population == "MatureNeu" | Population == "mature neutrophils" | Population == "TotalNeutrophil", is.finite(Z_Score))
if (nrow(mature_mfi) == 0) {
  mature_mfi <- mfi_combined %>% filter(is.finite(Z_Score))
}

mat_df <- mature_mfi %>%
  select(marker, SampleID, Z_Score) %>%
  group_by(marker, SampleID) %>%
  summarise(val = mean(Z_Score, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from = SampleID, values_from = val) %>%
  tibble::column_to_rownames("marker")

mat <- as.matrix(mat_df)
mat_meta <- meta %>% filter(sample_no %in% colnames(mat))
mat <- mat[, mat_meta$sample_no]

# -- 3. PCA ANALYSIS --
cat("Computing PCA...\n")
mat_complete <- t(na.omit(mat))
pca_res <- prcomp(mat_complete, scale. = FALSE)
pca_data <- as.data.frame(pca_res$x)
pca_data$sample_no <- rownames(pca_data)
pca_data <- pca_data %>% left_join(meta, by="sample_no")

variance_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

p_pca <- ggplot(pca_data, aes(x=PC1, y=PC2, color=condition, shape=bone_type)) +
  geom_point(size=4, alpha=0.85) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCA of High-Dimensional Cytometry Profiles",
    x = paste0("PC1 (", variance_explained[1], "% variance)"),
    y = paste0("PC2 (", variance_explained[2], "% variance)")
  ) +
  theme(legend.position = "right", panel.grid.minor = element_blank())

ggsave("Publication_Fig1_PCA.pdf", p_pca, width=8, height=5.5)
cat("Saved Publication_Fig1_PCA.pdf\n")

# -- 4. HIGH-DIMENSIONAL PHENOTYPIC PROFILING (HEATMAP) --
cat("Generating ComplexHeatmap...\n")
ha <- HeatmapAnnotation(
  Condition = mat_meta$condition,
  Compartment = mat_meta$bone_type,
  col = list(
    Condition = c("Healthy_Control" = "#2ecc71", "Inflammatory_Cohort" = "#e74c3c", "Oncology_Cohort" = "#3498db"),
    Compartment = c("Peripheral_Blood" = "#e6550d", "Bone_Marrow_Aspirate" = "#31a354", "Tissue_Compartment_A" = "#756bb1", "Tissue_Compartment_B" = "#636363")
  )
)

col_fun <- colorRamp2(c(-2, 0, 2), c("#2166ac", "#f7f7f7", "#b2182b"))
ht <- Heatmap(mat,
        name = "Z-Score\n(MFI)",
        col = col_fun,
        top_annotation = ha,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 10, fontface = "italic"),
        column_title = "Multiparametric Spectral Phenotypic Profile"
)

pdf("Publication_Fig2_Activation_Heatmap.pdf", width=10, height=7)
draw(ht)
dev.off()
cat("Saved Publication_Fig2_Activation_Heatmap.pdf\n")

# -- 5. ABSOLUTE COUNTS KINETICS & STATISTICS --
cat("Calculating Counts and Verification Statistics...\n")
counts_df <- parsed$counts
bead_df <- counts_df %>%
  filter(grepl("bead", Population, ignore.case=TRUE)) %>%
  mutate(SampleID = as.character(SampleID)) %>%
  select(SampleID, BeadCount = Count) %>%
  group_by(SampleID) %>%
  summarise(BeadCount = mean(BeadCount, na.rm=TRUE), .groups="drop")

counts_normalized <- counts_df %>%
  filter(!grepl("bead", Population, ignore.case=TRUE)) %>%
  mutate(SampleID = as.character(SampleID)) %>%
  left_join(bead_df, by="SampleID") %>%
  left_join(meta, by=c("SampleID" = "sample_no")) %>%
  mutate(AbsCount = (Count / BeadCount) * beads_input * (volume_total / volume_sample_Neu))

writeLines(paste0("Total Analyzed Samples: ", nrow(meta), "\n",
                  "MFI Feature Dimensions: ", nrow(mat), " markers x ", ncol(mat), " samples\n",
                  "PCA Variance Explained: PC1=", variance_explained[1], "%, PC2=", variance_explained[2], "%\n",
                  "Absolute Quantification Normalization: Complete\n"), stat_report)

close(stat_report)
cat("Saved Publication_Statistical_Report.txt\n")
cat("Publication pipeline successfully completed.\n")
