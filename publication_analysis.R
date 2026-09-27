suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ggpubr))

source("R/modules/mod_data_loader.R")

cat("Initializing Data Processing for Publication...\n")

# -- FILE PATHS --
meta_path <- "./Bone_Data/20251105 and 20251107 and 20251113 WT 5LO 15LO Naive and Zymozan Metadata.xlsx"
marker_path <- "./Bone_Data/20251105 and 20251107 and 20251113 Flow Panel Marker Log.csv"
facs_path <- "./Bone_Data/20251105 and 20251107 and 20251113 BoneMarrow FloJo Gating Analysis.xls"

# -- 1. DATA INGESTION --
meta <- read_excel(meta_path)
meta$sample_no <- as.character(meta$sample_no)
meta <- meta %>% mutate(
  Time = case_when(
    grepl("naive", condition, ignore.case=TRUE) ~ "0h",
    grepl("4", condition) ~ "4h",
    grepl("24", condition) ~ "24h",
    grepl("72", condition) ~ "72h",
    TRUE ~ "Unknown"
  ),
  Time = factor(Time, levels=c("0h", "4h", "24h", "72h")),
  Genotype = case_when(
    grepl("WT|Wt", condition, ignore.case=TRUE) ~ "WT",
    grepl("5LO", condition, ignore.case=TRUE) ~ "5LO",
    grepl("15LO", condition, ignore.case=TRUE) ~ "15LO",
    TRUE ~ "Unknown"
  ),
  Genotype = factor(Genotype, levels=c("WT", "5LO", "15LO"))
)

markers <- read.csv(marker_path)
parsed <- parse_legacy_facs(facs_path)

# Prepare Output Report File
stat_report <- file("Publication_Statistical_Report.txt", open="wt")
writeLines("Nature-Grade Publication Statistical Report\n============================================\n", stat_report)


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

mature_mfi <- mfi_combined %>% filter(Population == "mature neutrophils", is.finite(Z_Score))

mat_df <- mature_mfi %>%
  select(marker, SampleID, Z_Score) %>%
  group_by(marker, SampleID) %>%
  summarise(val = mean(Z_Score, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from = SampleID, values_from = val) %>%
  tibble::column_to_rownames("marker")

mat <- as.matrix(mat_df)
# Sync metadata with matrix columns
mat_meta <- meta %>% filter(sample_no %in% colnames(mat))
# Reorder matrix columns to match metadata order
mat <- mat[, mat_meta$sample_no]


# -- 3. PCA ANALYSIS --
cat("Computing PCA...\n")
# For PCA, we need samples as rows, features as columns.
# We also drop rows with NA to run prcomp.
mat_complete <- t(na.omit(mat))
pca_res <- prcomp(mat_complete, scale. = FALSE) # already scaled via Z-score
pca_data <- as.data.frame(pca_res$x)
pca_data$sample_no <- rownames(pca_data)
pca_data <- pca_data %>% left_join(meta, by="sample_no")

variance_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

p_pca <- ggplot(pca_data, aes(x=PC1, y=PC2, color=Genotype, shape=Time)) +
  geom_point(size=4, alpha=0.8) +
  theme_bw(base_size = 14) +
  labs(
    title = "PCA of Mature Neutrophil Phenotypes",
    x = paste0("PC1 (", variance_explained[1], "% variance)"),
    y = paste0("PC2 (", variance_explained[2], "% variance)")
  ) +
  scale_color_manual(values = c("WT" = "#1f77b4", "5LO" = "#d62728", "15LO" = "#2ca02c")) +
  theme(legend.position = "right", panel.grid.minor = element_blank())
  
if(nrow(pca_data) > 10) { # Add ellipses if enough data
    p_pca <- p_pca + stat_ellipse(aes(group=Genotype), type="t", linetype=2, alpha=0.5)
}

ggsave("Publication_Fig1_PCA.pdf", p_pca, width=7, height=5)
cat("Saved Publication_Fig1_PCA.pdf\n")


# -- 4. HIGH-DIMENSIONAL PHENOTYPIC PROFILING (HEATMAP) --
cat("Generating ComplexHeatmap...\n")
ha <- HeatmapAnnotation(
  Genotype = mat_meta$Genotype,
  Time = mat_meta$Time,
  Bone = mat_meta$bone_type,
  col = list(
    Genotype = c("WT" = "#1f77b4", "5LO" = "#d62728", "15LO" = "#2ca02c"),
    Time = c("0h"="#fcfbfd", "4h"="#bcbddc", "24h"="#807dba", "72h"="#3f007d"),
    Bone = c("Femur"="#e6550d", "Calvaria"="#31a354", "Lumbar"="#756bb1", "Sternum"="#636363")
  )
)

col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
ht <- Heatmap(mat,
        name = "Z-Score\n(MFI)",
        col = col_fun,
        top_annotation = ha,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_column_names = FALSE,
        row_names_gp = gpar(fontsize = 10, fontface = "italic"),
        column_title = "Mature Neutrophils Activation Profile (All Niches)"
)

pdf("Publication_Fig2_Activation_Heatmap.pdf", width=10, height=7)
draw(ht)
dev.off()
cat("Saved Publication_Fig2_Activation_Heatmap.pdf\n")


# -- 5. SPATIOTEMPORAL KINETICS & STATISTICS --
cat("Calculating Kinetics and Statistics...\n")
counts_df <- parsed$counts
bead_df <- counts_df %>%
  filter(grepl("single beads", Population, ignore.case=TRUE)) %>%
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

neu_kinetics <- counts_normalized %>%
  filter(Population %in% c("immature nuetrophils", "mature neutrophils"), Genotype %in% c("WT", "5LO"))

# Statistical Test: 5LO vs WT at 72h (Femur)
stat_data <- neu_kinetics %>% filter(Time == "72h", bone_type == "Femur", Population == "immature nuetrophils")
if(nrow(stat_data %>% filter(Genotype=="WT")) >= 2 && nrow(stat_data %>% filter(Genotype=="5LO")) >= 2) {
    res_wilcox <- wilcox.test(AbsCount ~ Genotype, data = stat_data, exact=FALSE)
    writeLines(paste0("\nStatistical Test: Wilcoxon Rank Sum\n",
                      "Comparison: Immature Neutrophils Accumulation in Femur at 72h (WT vs 5LO)\n",
                      "W-Statistic: ", res_wilcox$statistic, "\n",
                      "P-value: ", signif(res_wilcox$p.value, 4), "\n",
                      "Conclusion: ", ifelse(res_wilcox$p.value < 0.05, "Significant Difference", "Not Significant")), 
               stat_report)
}

# Plot: Immature & Mature Neu Kinetics (Femur & Calvaria)
plot_data <- neu_kinetics %>%
  filter(bone_type %in% c("Femur", "Calvaria")) %>%
  mutate(Population = ifelse(Population == "immature nuetrophils", "Immature Neu.", "Mature Neu.")) %>%
  group_by(Genotype, Time, bone_type, Population) %>%
  summarise(MeanAbs = mean(AbsCount, na.rm=TRUE), SE = sd(AbsCount, na.rm=TRUE)/sqrt(max(1, n()-1)), .groups="drop")

p_kinetics <- ggplot(plot_data, aes(x=Time, y=MeanAbs, color=Genotype, group=Genotype)) +
  geom_line(linewidth=1.2) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=pmax(0, MeanAbs-SE), ymax=MeanAbs+SE), width=0.2, linewidth=0.8) +
  facet_grid(Population ~ bone_type, scales="free_y") +
  theme_bw(base_size=14) +
  scale_color_manual(values = c("WT" = "#1f77b4", "5LO" = "#d62728")) +
  theme(strip.background = element_rect(fill="#f0f0f0", color="black"),
        strip.text = element_text(face="bold")) +
  labs(
    title = "Neutrophil Accumulation Kinetics (Mean ± SEM)",
    y = "Absolute Cell Count",
    x = "Time post-Zymosan (Hours)"
  ) 

ggsave("Publication_Fig3_Neutrophil_Kinetics.pdf", p_kinetics, width=10, height=7)
cat("Saved Publication_Fig3_Neutrophil_Kinetics.pdf\n")

# Stat Report for Calvaria CD106 vs Femur CD106 (WT 4h)
cd106_stat <- mfi_combined %>% filter(Population == "mature neutrophils", Genotype == "WT", Time == "4h", marker == "CD106", bone_type %in% c("Femur", "Calvaria"))
if(nrow(cd106_stat %>% filter(bone_type=="Femur")) >= 2 && nrow(cd106_stat %>% filter(bone_type=="Calvaria")) >= 2) {
    res_wilcox2 <- wilcox.test(Z_Score ~ bone_type, data = cd106_stat, exact=FALSE)
    writeLines(paste0("\nStatistical Test: Wilcoxon Rank Sum\n",
                      "Comparison: CD106 Activation State (Z-Score) in WT at 4h (Femur vs Calvaria)\n",
                      "W-Statistic: ", res_wilcox2$statistic, "\n",
                      "P-value: ", signif(res_wilcox2$p.value, 4), "\n",
                      "Conclusion: ", ifelse(res_wilcox2$p.value < 0.05, "Significant Spatial Heterogeneity", "Not Significant")), 
               stat_report)
}

close(stat_report)
cat("Saved Publication_Statistical_Report.txt\n")
cat("Publication pipeline successfully completed.\n")
