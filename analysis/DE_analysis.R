# =============================================================================
# Differential Expression Analysis — Cas13 RIP-Seq
# =============================================================================
# Paper: CRISPR-Cas13 Localises Within RNA-Rich Condensates To Drive
#        RNA Recognition and Cleavage
# Authors: Gill Jagjeet Singh GK, Hu W, Shembrey C, Hodel A,
#          Viskoboinik I, McMillan P, Trapani J, Fareh M
# Peter MacCallum Cancer Centre, Melbourne, Australia
#
# Description:
#   Differential expression analysis of RIP-Seq data across four fractions:
#     - Total    : Total cytoplasmic RNA
#     - CCEF     : Condensate-enriched fraction (CCEF)
#     - Total_IP : Cytoplasmic Cas13 immunoprecipitation
#     - CCEF_IP  : Condensate-enriched Cas13 immunoprecipitation
#
#   Note: R-safe names use underscores (Total_IP, CCEF_IP) internally.
#         Output CSV files use display names (Total-IP, CCEF-IP).
#
#   12 samples total: 3 biological replicates per condition
#   Pipeline: Galaxy (STAR + featureCounts, June 2024) -> R (limma-voom)
#
# R version: 4.3.3
# limma version: 3.58.1
# edgeR version: 4.0.16
# =============================================================================


# 1. Load libraries -----------------------------------------------------------

library(limma)
library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)


# 2. Load data ----------------------------------------------------------------

# Count matrix from featureCounts (run via Galaxy)
# Rows = Entrez Gene IDs, Columns = samples 1-12
# Update file path to match your local directory
count_matrix <- read_excel("COUNTDATA.xlsx", sheet = "CountData")
count_matrix <- as.data.frame(count_matrix)
rownames(count_matrix) <- count_matrix$`Entrez Gene ID`
count_matrix <- count_matrix[, -1]  # Remove Entrez Gene ID column

# Sample metadata
# SampleName: 1-12
# SampleType: Total (n=3), CCEF (n=3), Total_IP (n=3), CCEF_IP (n=3)
FactorData <- read_excel("COUNTDATA.xlsx", sheet = "Factor Data")
FactorData$SampleName <- as.character(FactorData$SampleName)

# Note: hyphens replaced with underscores for R compatibility
# Ensure FactorData SampleType column uses: Total, CCEF, Total_IP, CCEF_IP
FactorData$SampleType <- factor(
  FactorData$SampleType,
  levels = c("Total", "CCEF", "Total_IP", "CCEF_IP")
)


# 3. Create DGEList and filter low-count genes --------------------------------

dge <- DGEList(counts = count_matrix, group = FactorData$SampleType)
dge <- calcNormFactors(dge)

# Retain genes with CPM >= 0.5 in at least 2 samples
cpm_values <- cpm(count_matrix)
pass_filter <- rowSums(cpm_values >= 0.5) >= 2
dge_filtered <- dge[pass_filter, , keep.lib.sizes = FALSE]

cat("Genes after filtering:", nrow(dge_filtered), "\n")


# 4. Quality control ----------------------------------------------------------

# MDS plot - visualise sample clustering
plotMDS(
  dge_filtered,
  col = as.numeric(FactorData$SampleType),
  labels = FactorData$SampleName,
  main = "MDS Plot of RNA-seq Samples"
)
legend("topright",
       legend = levels(FactorData$SampleType),
       col = 1:length(levels(FactorData$SampleType)),
       pch = 16, cex = 0.8)


# 5. Voom transformation and model fitting ------------------------------------

design <- model.matrix(~0 + SampleType, data = FactorData)
colnames(design) <- levels(FactorData$SampleType)

# Voom transformation (mean-variance trend plot shown)
v <- voom(dge_filtered, design, plot = TRUE)

# Fit linear model
fit <- lmFit(v, design)


# 6. Define contrasts and run differential expression ------------------------
#
# Five comparisons:
#   1. Total_IP vs Total      (cytoplasmic Cas13 IP vs total cytoplasmic RNA)
#   2. CCEF_IP vs Total       (condensate IP vs total cytoplasmic RNA)
#   3. CCEF_IP vs CCEF        (condensate IP vs condensate fraction)
#   4. CCEF_IP vs Total_IP    (condensate IP vs cytoplasmic IP)
#   5. CCEF vs Total          (condensate fraction vs total cytoplasmic RNA)

contrast_matrix <- makeContrasts(
  TotalIP_vs_Total  = Total_IP - Total,
  CCEFIP_vs_Total   = CCEF_IP  - Total,
  CCEFIP_vs_CCEF    = CCEF_IP  - CCEF,
  CCEFIP_vs_TotalIP = CCEF_IP  - Total_IP,
  CCEF_vs_Total     = CCEF     - Total,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)


# 7. Extract results and annotate with gene symbols --------------------------

# Helper function: extract and annotate top table for a given contrast
get_annotated_results <- function(fit, contrast_name) {
  tt <- topTable(fit, coef = contrast_name, adjust = "BH", number = Inf)
  tt$EntrezGeneID <- rownames(tt)

  annotations <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = tt$EntrezGeneID,
    columns = c("SYMBOL", "GENENAME"),
    keytype = "ENTREZID"
  )

  merge(tt, annotations, by.x = "EntrezGeneID", by.y = "ENTREZID", all.x = TRUE)
}

# Extract results for all five contrasts
results_TotalIP_Total  <- get_annotated_results(fit2, "TotalIP_vs_Total")
results_CCEFIP_Total   <- get_annotated_results(fit2, "CCEFIP_vs_Total")
results_CCEFIP_CCEF    <- get_annotated_results(fit2, "CCEFIP_vs_CCEF")
results_CCEFIP_TotalIP <- get_annotated_results(fit2, "CCEFIP_vs_TotalIP")
results_CCEF_Total     <- get_annotated_results(fit2, "CCEF_vs_Total")

# Create output directories
dir.create("data",    showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Save DE results as CSV
write.csv(results_TotalIP_Total,  "data/Total-IP_Total.csv",   row.names = FALSE)
write.csv(results_CCEFIP_Total,   "data/CCEF-IP_Total.csv",    row.names = FALSE)
write.csv(results_CCEFIP_CCEF,    "data/CCEF-IP_CCEF.csv",     row.names = FALSE)
write.csv(results_CCEFIP_TotalIP, "data/CCEF-IP_Total-IP.csv", row.names = FALSE)
write.csv(results_CCEF_Total,     "data/CCEF_Total.csv",       row.names = FALSE)


# 8. Volcano plots ------------------------------------------------------------
#
# Significance thresholds: |log2FC| > 1.5, FDR < 0.05
# Top 10 enriched and depleted genes are labelled per plot.
# Plots are saved to the figures/ directory at 300 dpi.
#
# Note: Interactive versions of all volcano plots are available at
#       https://farehlab.github.io/Cas13IP-RNA-Seq-Profiles/
# Note: Raw count heatmaps (pheatmap) are in the paper supplementary figures.
#       Interactive heatmaps are available via the link above.

plot_volcano <- function(results, title) {

  df <- results %>%
    mutate(
      log2fc = logFC,
      pvalue = P.Value,
      fdr    = adj.P.Val,
      labels = SYMBOL,
      sig = case_when(
        fdr < 0.05 & log2fc >  1.5 ~ "Enriched",
        fdr < 0.05 & log2fc < -1.5 ~ "Depleted",
        TRUE ~ "Not Significant"
      )
    )

  # Top 10 enriched and depleted genes for labelling
  top_genes <- bind_rows(
    df %>% filter(sig == "Enriched") %>% arrange(desc(log2fc)) %>% slice_head(n = 10),
    df %>% filter(sig == "Depleted") %>% arrange(log2fc)       %>% slice_head(n = 10)
  )

  ggplot(df, aes(x = log2fc, y = -log10(pvalue))) +
    geom_vline(xintercept = c(-1.5, 1.5), col = "gray", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05),  col = "gray", linetype = "dashed") +
    geom_point(aes(colour = sig), size = 0.1) +
    scale_color_manual(
      values = c("Depleted" = "#00AFBB", "Not Significant" = "grey", "Enriched" = "#bb0c00")
    ) +
    coord_cartesian(xlim = c(-6, 6), ylim = c(0, 12)) +
    geom_label_repel(
      data = df %>% filter(labels %in% top_genes$SYMBOL),
      aes(label = labels),
      size = 2, fontface = "bold.italic",
      nudge_y = 0.5, max.overlaps = Inf,
      segment.color = "black", segment.size = 0.5,
      box.padding = 0.005, point.padding = 0.005,
      force = 1, force_pull = 0.5, min.segment.length = 0,
      show.legend = FALSE, na.rm = TRUE
    ) +
    labs(title = title, x = "log2(fold change)", y = "-log10(p-value)", colour = "") +
    theme_classic(base_size = 11) +
    theme(
      aspect.ratio = 1,
      axis.title.y = element_text(face = "bold", margin = margin(0, 20, 0, 0), size = rel(1.0), color = "black"),
      axis.title.x = element_text(face = "bold", margin = margin(20, 5, 5, 5), size = rel(1.1), color = "black"),
      plot.title   = element_text(hjust = 0.5)
    )
}

# Generate and save all five volcano plots
plots <- list(
  "Total-IP_vs_Total"   = plot_volcano(results_TotalIP_Total,  "Total-IP vs Total"),
  "CCEF-IP_vs_Total"    = plot_volcano(results_CCEFIP_Total,   "CCEF-IP vs Total"),
  "CCEF-IP_vs_CCEF"     = plot_volcano(results_CCEFIP_CCEF,    "CCEF-IP vs CCEF"),
  "CCEF-IP_vs_Total-IP" = plot_volcano(results_CCEFIP_TotalIP, "CCEF-IP vs Total-IP"),
  "CCEF_vs_Total"       = plot_volcano(results_CCEF_Total,     "CCEF vs Total")
)

for (name in names(plots)) {
  ggsave(paste0("figures/", name, ".png"),
         plot = plots[[name]], width = 7, height = 7, dpi = 300)
}

cat("Volcano plots saved to figures/\n")


# 9. Session info -------------------------------------------------------------

sessionInfo()

