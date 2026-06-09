## =========================================================
## 0. Environment setup
## =========================================================
suppressPackageStartupMessages({
  library(SeuratObject)
  library(Seurat)
  library(DoubletFinder)
  library(ggsci)
  library(Matrix)
  library(ggpubr)
  library(cowplot)
  library(gridExtra)
  library(ggnewscale)
  library(tidyverse)
  library(patchwork)
  library(tibble)
  library(ggrepel)
  library(EnhancedVolcano)
  library(igraph)
})

setwd("/Users/kuo/Desktop/TYC/20260423")

dir.create("results", showWarnings = FALSE)
dir.create("results/qc", showWarnings = FALSE)
dir.create("results/umap", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE)
dir.create("results/rds", showWarnings = FALSE)
dir.create("results/logs", showWarnings = FALSE)

sink("results/logs/sessionInfo.txt")
sessionInfo()
sink()

## ---------------------------------------------------------
## Optional: consistent plotting theme
## ---------------------------------------------------------
theme_set(
  theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
)

## =========================================================
## 1. Read 10X Genomics data
## =========================================================
Ctrl <- Read10X_h5("/Users/kuo/Desktop/TYC/ctrl sample_filtered_feature_bc_matrix.h5")
COPD_quit <- Read10X_h5("/Users/kuo/Desktop/TYC/COPD quit sample_filtered_feature_bc_matrix.h5")
COPD_active <- Read10X_h5("/Users/kuo/Desktop/TYC/COPD avtive sample_filtered_feature_bc_matrix.h5")

raw_dims <- tibble(
  sample = c("Control", "COPD-quit", "COPD-active"),
  n_genes = c(nrow(Ctrl), nrow(COPD_quit), nrow(COPD_active)),
  n_cells = c(ncol(Ctrl), ncol(COPD_quit), ncol(COPD_active))
)

write.csv(raw_dims, "results/tables/raw_matrix_dimensions.csv", row.names = FALSE)

## =========================================================
## 2. Create Seurat objects
## =========================================================
Ctrl_seurat <- CreateSeuratObject(
  counts = Ctrl,
  project = "Ctrl",
  min.cells = 3,
  min.features = 200
)

COPD_quit_seurat <- CreateSeuratObject(
  counts = COPD_quit,
  project = "COPD_quit",
  min.cells = 3,
  min.features = 200
)

COPD_active_seurat <- CreateSeuratObject(
  counts = COPD_active,
  project = "COPD_active",
  min.cells = 3,
  min.features = 200
)

Ctrl_seurat$group <- "Control"
COPD_quit_seurat$group <- "COPD-quit"
COPD_active_seurat$group <- "COPD-active"

pre_qc_counts <- tibble(
  group = c("Control", "COPD-quit", "COPD-active"),
  n_cells = c(ncol(Ctrl_seurat), ncol(COPD_quit_seurat), ncol(COPD_active_seurat))
)

write.csv(pre_qc_counts, "results/tables/cell_counts_pre_qc.csv", row.names = FALSE)

## =========================================================
## 3. QC metrics and filtering
## =========================================================
Ctrl_seurat[["percent.mt"]] <- PercentageFeatureSet(Ctrl_seurat, pattern = "^MT-")
COPD_quit_seurat[["percent.mt"]] <- PercentageFeatureSet(COPD_quit_seurat, pattern = "^MT-")
COPD_active_seurat[["percent.mt"]] <- PercentageFeatureSet(COPD_active_seurat, pattern = "^MT-")

Ctrl_seurat <- subset(
  Ctrl_seurat,
  subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 10
)

COPD_quit_seurat <- subset(
  COPD_quit_seurat,
  subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 10
)

COPD_active_seurat <- subset(
  COPD_active_seurat,
  subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 10
)

post_qc_counts <- tibble(
  group = c("Control", "COPD-quit", "COPD-active"),
  n_cells = c(ncol(Ctrl_seurat), ncol(COPD_quit_seurat), ncol(COPD_active_seurat))
)

write.csv(post_qc_counts, "results/tables/cell_counts_post_qc.csv", row.names = FALSE)

## =========================================================
## 4. Merge objects
## =========================================================
balmerged_singleR_obj <- merge(
  x = Ctrl_seurat,
  y = list(COPD_quit_seurat, COPD_active_seurat),
  add.cell.ids = c("Ctrl", "COPD_quit", "COPD_active"),
  project = "MergedProject"
)

balmerged_singleR_obj$group <- factor(
  balmerged_singleR_obj$group,
  levels = c("Control", "COPD-quit", "COPD-active")
)

DefaultAssay(balmerged_singleR_obj) <- "RNA"

merged_counts <- as.data.frame(table(balmerged_singleR_obj$group))
colnames(merged_counts) <- c("group", "n_cells")
write.csv(merged_counts, "results/tables/cell_counts_merged_post_qc.csv", row.names = FALSE)

## Save merged object after QC
saveRDS(
  balmerged_singleR_obj,
  "results/rds/BALF_scRNAseq_postQC_merged.rds"
)

# =========================================================
# paper_theme.R
# Unified publication theme system for the whole paper
# Suitable for volcano / violin / UMAP / dotplot / barplot
# =========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})
# =========================================================
# 1) Global color palettes
# =========================================================
group_colors <- c(
  "Control"     = "#A0A0A0",  
  "COPD-quit"   = "#E5C494",  
  "COPD-active" = "#B82E2E"   
)
# =========================================================
# 2) Main paper theme: with border
#    Use for volcano / barplot / dotplot / general figures
# =========================================================
paper_theme <- function(base_size = 14, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) +
    theme(
      # -------- plot title / subtitle --------
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      plot.subtitle = element_text(
        size = 12,
        face = "plain",
        hjust = 0.5,
        color = "black"
      ),
      
      # -------- axis --------
      axis.title = element_text(
        size = 14,
        face = "bold",
        color = "black"
      ),
      axis.text = element_text(
        size = 12,
        color = "black"
      ),
      axis.line = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.7,
        color = "black"
      ),
      axis.ticks.length = unit(0.18, "cm"),
      
      # -------- legend --------
      legend.title = element_text(
        size = 12,
        face = "bold",
        color = "black"
      ),
      legend.text = element_text(
        size = 11,
        color = "black"
      ),
      legend.key = element_blank(),
      legend.background = element_blank(),
      
      # -------- panel --------
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(
        fill = NA,
        color = "black",
        linewidth = 0.8
      ),
      
      # -------- facet strip --------
      strip.background = element_rect(
        fill = "white",
        color = "black",
        linewidth = 0.8
      ),
      strip.text = element_text(
        size = 12,
        face = "bold",
        color = "black"
      ),
      
      # -------- overall --------
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# =========================================================
# 3) Axes-only theme: no border
#    Use for violin / boxplot / simple scatter
# =========================================================
paper_theme_axes <- function(base_size = 14, base_family = "") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      plot.subtitle = element_text(
        size = 12,
        face = "plain",
        hjust = 0.5,
        color = "black"
      ),
      axis.title = element_text(
        size = 14,
        face = "bold",
        color = "black"
      ),
      axis.text = element_text(
        size = 12,
        color = "black"
      ),
      axis.line = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.7,
        color = "black"
      ),
      axis.ticks.length = unit(0.18, "cm"),
      legend.title = element_text(
        size = 12,
        face = "bold",
        color = "black"
      ),
      legend.text = element_text(
        size = 11,
        color = "black"
      ),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# =========================================================
# 4) Specialized theme: volcano
# =========================================================
theme_volcano_paper <- function() {
  paper_theme(base_size = 14) +
    theme(
      legend.position = "none",
      plot.margin = margin(8, 18, 8, 8)
    )
}

# =========================================================
# 5) Specialized theme: violin
# =========================================================
theme_violin_paper <- function() {
  paper_theme_axes(base_size = 14) +
    theme(
      legend.position = "none",
      plot.margin = margin(12, 14, 20, 12)
    )
}

# =========================================================
# 6) Specialized theme: UMAP / trajectory
# =========================================================
theme_umap_paper <- function() {
  theme_void(base_size = 14) +
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      legend.title = element_text(
        size = 12,
        face = "bold",
        color = "black"
      ),
      legend.text = element_text(
        size = 11,
        color = "black"
      ),
      plot.margin = margin(8, 8, 8, 8)
    )
}

# =========================================================
# 7) Specialized theme: dotplot
# =========================================================
theme_dotplot_paper <- function() {
  paper_theme(base_size = 14) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 11,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 12,
        color = "black"
      ),
      legend.position = "right"
    )
}

# =========================================================
# 8) Manual scales
# =========================================================
scale_fill_paper_groups <- function(drop = FALSE) {
  scale_fill_manual(values = group_colors, drop = drop)
}

scale_color_paper_groups <- function(drop = FALSE) {
  scale_color_manual(values = group_colors, drop = drop)
}

scale_fill_paper_macrophage <- function(drop = FALSE) {
  scale_fill_manual(values = mac_colors, drop = drop)
}

scale_color_paper_macrophage <- function(drop = FALSE) {
  scale_color_manual(values = mac_colors, drop = drop)
}

# =========================================================
# 9) Optional helper: TIFF / PDF / PNG saving
# =========================================================
save_plot_paper <- function(plot_obj, filename_base, outdir = ".",
                            width = 4.8, height = 5.2, dpi = 600) {
  
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  
  # PDF
  ggsave(
    filename = file.path(outdir, paste0(filename_base, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
  
  # PNG
  ggsave(
    filename = file.path(outdir, paste0(filename_base, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  
  # TIFF
  ggsave(
    filename = file.path(outdir, paste0(filename_base, ".tiff")),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    compression = "lzw",
    bg = "white"
  )
}

# =========================================================
# 10) Set global theme once
# =========================================================
theme_set(paper_theme())


## =========================================================
## 5. QC visualization
## =========================================================

qc_df <- balmerged_singleR_obj@meta.data %>%
  dplyr::mutate(cell = rownames(.))

write.csv(
  qc_df,
  "results/tables/qc_metadata_merged.csv",
  row.names = FALSE
)

## =========================================================
## QC violin plot
## =========================================================

p_vln <- VlnPlot(
  balmerged_singleR_obj,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  group.by = "group",
  pt.size = 0
)

p_vln <- p_vln &
  scale_fill_manual(values = group_colors) &
  theme_violin_paper() &
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 12,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 14,
      face = "bold"
    ),
    axis.title.y = element_text(
      size = 14,
      face = "bold"
    ),
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 38,
      l = 10
    )
  )

## =========================================================
## Histogram
## =========================================================

p_hist <- ggplot(
  qc_df,
  aes(x = nFeature_RNA)
) +
  geom_histogram(
    bins = 100,
    fill = "grey70",
    color = "black"
  ) +
  geom_vline(
    xintercept = c(200, 5000),
    linetype = "dashed",
    color = "red"
  ) +
  paper_theme() +
  labs(
    title = "Distribution of nFeature_RNA",
    x = "nFeature_RNA",
    y = "Number of cells"
  )

## =========================================================
## Scatter plot
## =========================================================

p_scatter <- ggplot(
  qc_df,
  aes(
    x = nCount_RNA,
    y = nFeature_RNA,
    color = group
  )
) +
  geom_point(
    size = 0.4,
    alpha = 0.4
  ) +
  geom_hline(
    yintercept = c(200, 5000),
    linetype = "dashed",
    color = "red"
  ) +
  scale_color_manual(values = group_colors) +
  paper_theme() +
  labs(
    title = "nCount_RNA vs nFeature_RNA",
    x = "nCount_RNA",
    y = "nFeature_RNA"
  )

## =========================================================
## Save figures
## =========================================================

ggsave(
  "results/qc/QC_violinplot_group.png",
  p_vln,
  width = 10,
  height = 5.4,
  dpi = 300,
  bg = "white"
)

ggsave(
  "results/qc/QC_violinplot_group.tiff",
  p_vln,
  width = 10,
  height = 5.4,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  "results/qc/QC_nFeature_RNA_histogram.png",
  p_hist,
  width = 7.2,
  height = 5.2,
  dpi = 300,
  bg = "white"
)

ggsave(
  "results/qc/QC_nFeature_RNA_histogram.tiff",
  p_hist,
  width = 7.2,
  height = 5.2,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  "results/qc/QC_nCount_vs_nFeature_scatter.png",
  p_scatter,
  width = 7.2,
  height = 5.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  "results/qc/QC_nCount_vs_nFeature_scatter.tiff",
  p_scatter,
  width = 7.2,
  height = 5.5,
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

## =========================================================
## 6. Check RNA layers for integration
## =========================================================
DefaultAssay(balmerged_singleR_obj) <- "RNA"

layer_info <- tibble(layer = Layers(balmerged_singleR_obj[["RNA"]]))
print(layer_info)
write.csv(layer_info, "results/tables/rna_layers_before_integration.csv", row.names = FALSE)

## =========================================================
## 7. Preprocessing before integration
## =========================================================
balmerged_singleR_obj <- NormalizeData(
  balmerged_singleR_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

balmerged_singleR_obj <- FindVariableFeatures(
  balmerged_singleR_obj,
  selection.method = "vst",
  nfeatures = 2000
)

balmerged_singleR_obj <- ScaleData(balmerged_singleR_obj)

balmerged_singleR_obj <- RunPCA(balmerged_singleR_obj)

p_elbow <- ElbowPlot(balmerged_singleR_obj)
ggsave("results/qc/ElbowPlot_pre_integration.png", p_elbow, width = 6.5, height = 4.5, dpi = 300)

## =========================================================
## 8. Layer integration
## =========================================================
balmerged_singleR_obj <- IntegrateLayers(
  object = balmerged_singleR_obj,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)

saveRDS(
  balmerged_singleR_obj,
  "results/rds/BALF_scRNAseq_postIntegration.rds"
)

## =========================================================
balmerged_singleR_obj <- FindNeighbors(
  balmerged_singleR_obj,
  reduction = "integrated.cca",
  dims = 1:20
)

balmerged_singleR_obj <- FindClusters(
  balmerged_singleR_obj,
  resolution = 0.5
)

balmerged_singleR_obj <- RunUMAP(
  balmerged_singleR_obj,
  reduction = "integrated.cca",
  dims = 1:20
)

cluster_counts_before <- as.data.frame(table(balmerged_singleR_obj$seurat_clusters))
colnames(cluster_counts_before) <- c("cluster", "n_cells")
write.csv(cluster_counts_before, "results/tables/cluster_counts_before_doublet_removal.csv", row.names = FALSE)

p_umap_group_before <- DimPlot(
  balmerged_singleR_obj,
  reduction = "umap",
  group.by = "group"
) + ggtitle("Group")

p_umap_split_before <- DimPlot(
  balmerged_singleR_obj,
  reduction = "umap",
  split.by = "group",
  label = TRUE
) + ggtitle("Split by group")

p_umap_cluster_before <- DimPlot(
  balmerged_singleR_obj,
  reduction = "umap",
  label = TRUE
) + ggtitle("Clusters")

ggsave("results/umap/UMAP_Group_before_doublet_removal.png", p_umap_group_before, width = 8, height = 6, dpi = 300)
ggsave("results/umap/UMAP_Split_before_doublet_removal.png", p_umap_split_before, width = 10, height = 6, dpi = 300)
ggsave("results/umap/UMAP_Clusters_before_doublet_removal.png", p_umap_cluster_before, width = 8, height = 6, dpi = 300)

## =========================================================
## 10. Prepare output folders and DoubletFinder input
## =========================================================
dir.create("results", showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/umap", recursive = TRUE, showWarnings = FALSE)
dir.create("results/rds", recursive = TRUE, showWarnings = FALSE)

DefaultAssay(balmerged_singleR_obj) <- "RNA"
balmerged_singleR_obj[["RNA"]] <- JoinLayers(balmerged_singleR_obj[["RNA"]])

## Ensure meta.data is a base data.frame to avoid downstream issues caused by tibble objects
balmerged_singleR_obj@meta.data <- as.data.frame(
  balmerged_singleR_obj@meta.data,
  stringsAsFactors = FALSE
)
rownames(balmerged_singleR_obj@meta.data) <- colnames(balmerged_singleR_obj)


old_df_cols <- grep(
  "^(DF.classifications_|pANN_)",
  colnames(balmerged_singleR_obj@meta.data),
  value = TRUE
)

if (length(old_df_cols) > 0) {
  message("Removing old DoubletFinder columns: ", paste(old_df_cols, collapse = ", "))
  balmerged_singleR_obj@meta.data <- balmerged_singleR_obj@meta.data[
    , !colnames(balmerged_singleR_obj@meta.data) %in% old_df_cols,
    drop = FALSE
  ]
}

## =========================================================
## 11. DoubletFinder parameter sweep
## =========================================================
sweep.res.list <- paramSweep(
  balmerged_singleR_obj,
  PCs = 1:20,
  sct = FALSE
)

sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)

if (!"BCmetric" %in% colnames(sweep.stats)) {
  sweep.stats$BCmetric <- sweep.stats$BCreal
}

write.csv(
  sweep.stats,
  "results/tables/DoubletFinder_paramSweep_stats.csv",
  row.names = FALSE
)

bcmvn <- sweep.stats %>%
  dplyr::group_by(pK) %>%
  dplyr::summarise(BCmetric = mean(BCmetric), .groups = "drop") %>%
  dplyr::arrange(desc(BCmetric))

write.csv(
  bcmvn,
  "results/tables/DoubletFinder_pK_summary.csv",
  row.names = FALSE
)

pK <- as.numeric(as.character(bcmvn$pK[1]))

selected_pk <- data.frame(selected_pK = pK)
write.csv(
  selected_pk,
  "results/tables/DoubletFinder_selected_pK.csv",
  row.names = FALSE
)

## =========================================================
## 12. Estimate expected doublet number
## =========================================================
if (!"seurat_clusters" %in% colnames(balmerged_singleR_obj@meta.data)) {
  stop("seurat_clusters not found in meta.data. Please run FindClusters() before DoubletFinder.")
}

annotations <- balmerged_singleR_obj$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)

nExp_poi <- round(0.08 * ncol(balmerged_singleR_obj))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

doublet_estimates <- data.frame(
  homotypic.prop = homotypic.prop,
  nExp_poi = nExp_poi,
  nExp_poi.adj = nExp_poi.adj
)

write.csv(
  doublet_estimates,
  "results/tables/DoubletFinder_expected_doublets.csv",
  row.names = FALSE
)

cat("Selected pK =", pK, "\n")
cat("Homotypic proportion =", homotypic.prop, "\n")
cat("Expected doublets (raw) =", nExp_poi, "\n")
cat("Expected doublets (adjusted) =", nExp_poi.adj, "\n")

## =========================================================
## 13. Custom DoubletFinder function (based on your old code)
## =========================================================
doubletFinder_custom <- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, sct = FALSE) {
  require(Seurat)
  require(fields)
  
  real.cells <- colnames(seu)
  
  data <- GetAssayData(
    object = seu,
    assay = "RNA",
    layer = "counts"
  )[, real.cells, drop = FALSE]
  
  n_real.cells <- length(real.cells)
  n_doublets <- round(n_real.cells / (1 - pN) - n_real.cells)
  
  message("Creating ", n_doublets, " artificial doublets...")
  
  set.seed(123)
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  
  doublets <- (data[, real.cells1, drop = FALSE] + data[, real.cells2, drop = FALSE]) / 2
  colnames(doublets) <- paste0("X", seq_len(n_doublets))
  
  data_wdoublets <- cbind(data, doublets)
  
  seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
  DefaultAssay(seu_wdoublets) <- "RNA"
  
  seu_wdoublets <- NormalizeData(seu_wdoublets)
  seu_wdoublets <- FindVariableFeatures(seu_wdoublets)
  seu_wdoublets <- ScaleData(seu_wdoublets)
  seu_wdoublets <- RunPCA(seu_wdoublets, npcs = max(PCs))
  
  pca.coord <- Embeddings(seu_wdoublets, "pca")[, PCs, drop = FALSE]
  
  ## cell-cell distance matrix
  dist.mat <- fields::rdist(pca.coord)
  
  ## Determine the number of neighbors based on the number of cells, not the number of PCs
  k <- round(nrow(pca.coord) * pK)
  k <- max(k, 5)
  
  pANN <- rep(0, n_real.cells)
  names(pANN) <- real.cells
  
  for (i in seq_len(n_real.cells)) {
    neighbors <- order(dist.mat[, i])[2:(k + 1)]
    pANN[i] <- sum(neighbors > n_real.cells) / k
  }
  
  classifications <- rep("Singlet", n_real.cells)
  classifications[order(pANN, decreasing = TRUE)[seq_len(min(nExp, n_real.cells))]] <- "Doublet"
  
  seu@meta.data[[paste0("pANN_", pN, "_", pK, "_", nExp)]] <- pANN[colnames(seu)]
  seu@meta.data[[paste0("DF.classifications_", pN, "_", pK, "_", nExp)]] <- classifications
  
  return(seu)
}

## =========================================================
## 14. Run custom DoubletFinder
## =========================================================
pK <- as.numeric(pK)[1]
nExp_poi.adj <- as.integer(round(as.numeric(nExp_poi.adj)[1]))

stopifnot(!is.na(pK), length(pK) == 1)
stopifnot(!is.na(nExp_poi.adj), length(nExp_poi.adj) == 1)

balmerged_singleR_obj <- doubletFinder_custom(
  seu = balmerged_singleR_obj,
  PCs = 1:20,
  pN = 0.25,
  pK = pK,
  nExp = nExp_poi.adj,
  reuse.pANN = FALSE,
  sct = FALSE
)

df_col <- grep(
  "^DF\\.classifications",
  colnames(balmerged_singleR_obj@meta.data),
  value = TRUE
)

pann_col <- grep(
  "^pANN",
  colnames(balmerged_singleR_obj@meta.data),
  value = TRUE
)

if (length(df_col) == 0) stop("No DF.classifications column found after DoubletFinder.")
if (length(pann_col) == 0) stop("No pANN column found after DoubletFinder.")

df_col <- tail(df_col, 1)
pann_col <- tail(pann_col, 1)

message("Using DF column: ", df_col)
message("Using pANN column: ", pann_col)

print(head(balmerged_singleR_obj@meta.data[, df_col, drop = FALSE]))

## Export metadata
meta_keep <- c("group", "seurat_clusters")
meta_keep <- meta_keep[meta_keep %in% colnames(balmerged_singleR_obj@meta.data)]

doublet_meta <- balmerged_singleR_obj@meta.data %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  dplyr::mutate(cell = rownames(.)) %>%
  dplyr::select(cell, dplyr::all_of(meta_keep), dplyr::all_of(c(df_col, pann_col)))

write.csv(
  doublet_meta,
  "results/tables/DoubletFinder_classification_metadata.csv",
  row.names = FALSE
)

## Export counts
doublet_counts <- as.data.frame(table(balmerged_singleR_obj@meta.data[[df_col]]))
colnames(doublet_counts) <- c("classification", "n_cells")

write.csv(
  doublet_counts,
  "results/tables/DoubletFinder_classification_counts.csv",
  row.names = FALSE
)

## =========================================================
## 15. Visualize DoubletFinder results
## =========================================================

if (!"umap" %in% names(balmerged_singleR_obj@reductions)) {
  if ("integrated.cca" %in% names(balmerged_singleR_obj@reductions)) {
    balmerged_singleR_obj <- RunUMAP(
      balmerged_singleR_obj,
      reduction = "integrated.cca",
      dims = 1:20
    )
  } else {
    balmerged_singleR_obj <- RunUMAP(
      balmerged_singleR_obj,
      reduction = "pca",
      dims = 1:20
    )
  }
}

p_doublet_class <- DimPlot(
  balmerged_singleR_obj,
  group.by = df_col,
  reduction = "umap"
) + ggtitle("DoubletFinder Classification")

ggsave(
  "results/umap/UMAP_DoubletFinder_classification.png",
  plot = p_doublet_class,
  width = 8,
  height = 6,
  dpi = 300
)

## =========================================================
## 16. Keep singlets only
## =========================================================
singlet_cells <- rownames(
  balmerged_singleR_obj@meta.data[
    balmerged_singleR_obj@meta.data[[df_col]] == "Singlet",
    ,
    drop = FALSE
  ]
)

balmerged_singleR_obj <- subset(
  balmerged_singleR_obj,
  cells = singlet_cells
)

cell_counts_after_doublet <- as.data.frame(table(balmerged_singleR_obj$group))
colnames(cell_counts_after_doublet) <- c("group", "n_cells")

write.csv(
  cell_counts_after_doublet,
  "results/tables/cell_counts_after_doublet_removal.csv",
  row.names = FALSE
)

saveRDS(
  balmerged_singleR_obj,
  "results/rds/BALF_scRNAseq_singlets.rds"
)

## =========================================================
## 17. Recompute graph / clustering / UMAP after doublet removal
## =========================================================
if ("integrated.cca" %in% names(balmerged_singleR_obj@reductions)) {
  balmerged_singleR_obj <- FindNeighbors(
    balmerged_singleR_obj,
    reduction = "integrated.cca",
    dims = 1:20
  )
  
  balmerged_singleR_obj <- FindClusters(
    balmerged_singleR_obj,
    resolution = 0.5
  )
  
  balmerged_singleR_obj <- RunUMAP(
    balmerged_singleR_obj,
    reduction = "integrated.cca",
    dims = 1:20
  )
} else {
  balmerged_singleR_obj <- NormalizeData(balmerged_singleR_obj)
  balmerged_singleR_obj <- FindVariableFeatures(balmerged_singleR_obj)
  balmerged_singleR_obj <- ScaleData(balmerged_singleR_obj)
  balmerged_singleR_obj <- RunPCA(balmerged_singleR_obj, npcs = 30)
  
  balmerged_singleR_obj <- FindNeighbors(
    balmerged_singleR_obj,
    reduction = "pca",
    dims = 1:20
  )
  
  balmerged_singleR_obj <- FindClusters(
    balmerged_singleR_obj,
    resolution = 0.5
  )
  
  balmerged_singleR_obj <- RunUMAP(
    balmerged_singleR_obj,
    reduction = "pca",
    dims = 1:20
  )
}

cluster_counts_after <- as.data.frame(table(balmerged_singleR_obj$seurat_clusters))
colnames(cluster_counts_after) <- c("cluster", "n_cells")

write.csv(
  cluster_counts_after,
  "results/tables/cluster_counts_after_doublet_removal.csv",
  row.names = FALSE
)

p_umap_group_after <- DimPlot(
  balmerged_singleR_obj,
  group.by = "group",
  reduction = "umap"
) + ggtitle("After Doublet Removal")

p_umap_split_after <- DimPlot(
  balmerged_singleR_obj,
  reduction = "umap",
  split.by = "group",
  label = TRUE
) + ggtitle("Split by Group")

p_umap_cluster_after <- DimPlot(
  balmerged_singleR_obj,
  reduction = "umap",
  label = TRUE
) + ggtitle("Clusters after Doublet Removal")

ggsave(
  "results/umap/UMAP_Group_after_doublet_removal.png",
  plot = p_umap_group_after,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "results/umap/UMAP_Split_after_doublet_removal.png",
  plot = p_umap_split_after,
  width = 10,
  height = 6,
  dpi = 300
)

ggsave(
  "results/umap/UMAP_Clusters_after_doublet_removal.png",
  plot = p_umap_cluster_after,
  width = 8,
  height = 6,
  dpi = 300
)

table(balmerged_singleR_obj$group)
saveRDS(
  balmerged_singleR_obj,
  "results/rds/BALF_scRNAseq_final_annotated.rds"
)

## =========================================================
## 0. Setup
## =========================================================

## install BiocManager first
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c(
      "SingleR", "celldex", "scran",
      "SummarizedExperiment", "SingleCellExperiment"
    )) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

req_pkgs <- c(
  "SingleR", "celldex", "Seurat", "SeuratObject",
  "SummarizedExperiment", "SingleCellExperiment",
  "ggplot2", "scran", "dplyr", "patchwork", "ggpubr",
  "ggrepel", "scales", "Matrix", "tibble"
)

invisible(lapply(req_pkgs, install_if_missing))
invisible(lapply(req_pkgs, library, character.only = TRUE))

options(future.globals.maxSize = 5 * 1024^3)

## =========================================================
## Global publication-style theme
## =========================================================

theme_set(
  theme_classic(base_size = 14) +
    theme(
      plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title   = element_text(size = 14, face = "bold"),
      axis.text    = element_text(size = 12),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text  = element_text(size = 11),
      strip.text   = element_text(size = 13, face = "bold")
    )
)

theme_umap_paper <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14, face = "bold", color = "black"),
      axis.text = element_text(size = 12, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 12),
      panel.border = element_blank(),
      panel.grid = element_blank()
    )
}

save_plot <- function(filename, plot, width, height, dpi = 300) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
}

save_plot_paper <- function(
    plot_obj,
    filename_base,
    outdir = "results",
    width = 6,
    height = 5,
    dpi = 600
) {
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  
  ggsave(
    filename = file.path(outdir, paste0(filename_base, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi
  )
  
  ggsave(
    filename = file.path(outdir, paste0(filename_base, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height
  )
}

sessionInfo()

## =========================================================
## Load packages
## =========================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(SingleR)
library(celldex)

## =========================================================
## Output directories
## =========================================================

outdir_umap <- "results/umap"
outdir_barplot <- "results/barplot"

dir.create(outdir_umap, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_barplot, recursive = TRUE, showWarnings = FALSE)

cat("Current working directory:\n")
print(getwd())

## =========================================================
## 1. SingleR annotation
## =========================================================

DefaultAssay(balmerged_singleR_obj) <- "RNA"

ref_grch38 <- celldex::BlueprintEncodeData()

expr <- GetAssayData(
  balmerged_singleR_obj,
  assay = "RNA",
  layer = "data"
)

clusters <- Idents(balmerged_singleR_obj)

singleR_results <- SingleR(
  test = expr,
  ref = ref_grch38,
  labels = ref_grch38$label.main,
  de.method = "classic",
  fine.tune = TRUE,
  quantile = 0.8
)

cluster_results <- SingleR(
  test = expr,
  ref = ref_grch38,
  labels = ref_grch38$label.main,
  clusters = clusters
)

rownames(singleR_results) <- colnames(balmerged_singleR_obj)
balmerged_singleR_obj$SingleR_labels <- singleR_results$labels

## ---------------------------------------------------------
## 1A. SingleR barplot by Seurat clusters
## ---------------------------------------------------------

p_singler_bar <- ggplot(
  balmerged_singleR_obj@meta.data,
  aes(x = seurat_clusters, fill = SingleR_labels)
) +
  geom_bar() +
  labs(
    title = "SingleR annotations by clusters",
    x = "Seurat clusters",
    y = "Cell count",
    fill = "SingleR labels"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.text = element_text(size = 13)
  )

ggsave(
  filename = file.path(outdir_umap, "SingleR_Annotations_by_Clusters.png"),
  plot = p_singler_bar,
  width = 8,
  height = 6,
  dpi = 600
)

## ---------------------------------------------------------
## 1B. UMAP colored by SingleR labels
## ---------------------------------------------------------

label_counts <- table(balmerged_singleR_obj$SingleR_labels)

balmerged_singleR_obj$SingleR_labels_with_counts <- paste0(
  balmerged_singleR_obj$SingleR_labels,
  " (",
  label_counts[balmerged_singleR_obj$SingleR_labels],
  ")"
)

p_umap_singler <- DimPlot(
  balmerged_singleR_obj,
  reduction = "umap",
  group.by = "SingleR_labels_with_counts",
  label = TRUE,
  repel = TRUE,
  label.size = 4.5,
  pt.size = 0.5
) +
  ggtitle("UMAP colored by SingleR labels") +
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold", color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key = element_blank()
  )

ggsave(
  filename = file.path(outdir_umap, "UMAP_SingleR_Labels.png"),
  plot = p_umap_singler,
  width = 6.5,
  height = 5.5,
  dpi = 600
)

## =========================================================
## 2. Merge SingleR labels into broader categories
## =========================================================

balmerged_singleR_obj$Merged_Labels <- dplyr::recode(
  balmerged_singleR_obj$SingleR_labels,
  
  "CD8+ T-cells"     = "CD8+ T-cells",
  "CD4+ T-cells"     = "CD4+ T-cells",
  "B-cell"           = "B-cells",
  "NK cell"          = "NK cells",
  
  "Macrophage"       = "Macrophages",
  "Monocytes"        = "Monocytes",
  "Neutrophils"      = "Neutrophils",
  "Eosinophils"      = "Eosinophils",
  "DC"               = "DC",
  
  "Keratinocytes"    = "Epithelial cells",
  "Epithelial cells" = "Epithelial cells",
  "Fibroblasts"      = "Fibroblasts",
  "Mesangial cells"  = "Mesangial cells",
  "Adipocytes"       = "Adipocytes"
)

celltype_levels <- c(
  "Macrophages",
  "Monocytes",
  "DC",
  "Neutrophils",
  "Eosinophils",
  "CD4+ T-cells",
  "CD8+ T-cells",
  "NK cells",
  "B-cells",
  "Epithelial cells",
  "Fibroblasts",
  "Mesangial cells",
  "Adipocytes"
)

balmerged_singleR_obj$Merged_Labels <- factor(
  balmerged_singleR_obj$Merged_Labels,
  levels = celltype_levels
)

print(table(balmerged_singleR_obj$Merged_Labels, useNA = "ifany"))

Idents(balmerged_singleR_obj) <- "Merged_Labels"

## =========================================================
## Custom colors
## =========================================================

custom_colors <- c(
  "Macrophages"      = "#3B6FB6",
  "Monocytes"        = "#5FA8A6",
  "DC"               = "#9DC3E6",
  "Neutrophils"      = "#E08D3C",
  "Eosinophils"      = "#C96B3B",
  
  "CD4+ T-cells"     = "#A66BBE",
  "CD8+ T-cells"     = "#7BC96F",
  "NK cells"         = "#2F8F83",
  "B-cells"          = "#D98BB9",
  
  "Epithelial cells" = "#B33C4A",
  "Fibroblasts"      = "#8C6D5A",
  "Mesangial cells"  = "#B8B0A5",
  "Adipocytes"       = "#7A5C58"
)

present_levels <- celltype_levels[
  celltype_levels %in%
    unique(as.character(balmerged_singleR_obj$Merged_Labels))
]

main_cols_use <- custom_colors[present_levels]

print(main_cols_use)

## =========================================================
## 2A. UMAP with merged labels
## =========================================================

p_umap_merged <- DimPlot(
  balmerged_singleR_obj,
  reduction = "umap",
  group.by = "Merged_Labels",
  cols = main_cols_use,
  label = FALSE,
  repel = TRUE,
  pt.size = 0.5
) +
  ggtitle("Cell type annotation") +
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = 16,
        size = 4
      )
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold", color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13),
    legend.key = element_blank()
  )

ggsave(
  filename = file.path(outdir_umap, "DimPlot_Merged_Labels.png"),
  plot = p_umap_merged,
  width = 8,
  height = 6,
  dpi = 600
)

## =========================================================
## 2B. Split UMAP by group
## =========================================================

balmerged_singleR_obj$group <- factor(
  balmerged_singleR_obj$group,
  levels = c("Control", "COPD-quit", "COPD-active")
)

p_umap_merged_split <- DimPlot(
  balmerged_singleR_obj,
  reduction = "umap",
  group.by = "Merged_Labels",
  split.by = "group",
  cols = main_cols_use,
  label = FALSE,
  repel = TRUE,
  pt.size = 0.35,
  ncol = 3
) +
  ggtitle("Cell type distribution by group") +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = 16,
        size = 4
      )
    )
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12, face = "bold", color = "black")
  )

ggsave(
  filename = file.path(outdir_umap, "DimPlot_Merged_Labels_Split_by_Group.png"),
  plot = p_umap_merged_split,
  width = 12,
  height = 5,
  dpi = 600
)

## =========================================================
## 2C. Cell proportion barplot
## =========================================================

df_main <- balmerged_singleR_obj@meta.data %>%
  dplyr::filter(
    !is.na(group),
    !is.na(Merged_Labels)
  ) %>%
  dplyr::mutate(
    Condition = factor(
      group,
      levels = c("Control", "COPD-quit", "COPD-active")
    ),
    Cell_Type = factor(
      Merged_Labels,
      levels = celltype_levels
    )
  ) %>%
  dplyr::filter(
    !is.na(Condition),
    !is.na(Cell_Type)
  ) %>%
  dplyr::count(Condition, Cell_Type, name = "CellNumber") %>%
  dplyr::group_by(Condition) %>%
  dplyr::mutate(
    Percentage = CellNumber / sum(CellNumber) * 100
  ) %>%
  dplyr::ungroup()

group_colors <- c(
  "Control"     = "#A0A0A0",
  "COPD-quit"   = "#E5C494",
  "COPD-active" = "#B82E2E"
)

p_main_prop <- ggplot(
  df_main,
  aes(
    x = Cell_Type,
    y = Percentage,
    fill = Condition
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    colour = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = group_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "Cell proportion (%)",
    x = NULL,
    y = "Cell proportion (%)",
    fill = "Group"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 13),
    plot.margin = margin(10, 10, 35, 10)
  )

ggsave(
  filename = file.path(outdir_barplot, "Cell_Proportion_Main_Celltypes.png"),
  plot = p_main_prop,
  width = 8,
  height = 5.5,
  dpi = 600
)

## =========================================================
## Check saved files
## =========================================================

cat("\nSaved files:\n")
print(list.files("results", recursive = TRUE))



## ---------------------------------------------------------
## 2D. Save merged object and annotation
## ---------------------------------------------------------

saveRDS(
  balmerged_singleR_obj,
  file = "balmerged_singleR_final_obj.rds"
)

annot_df <- data.frame(
  CellID = colnames(balmerged_singleR_obj),
  CellType = balmerged_singleR_obj$Merged_Labels
)

write.csv(
  annot_df,
  file = "celltype_annotation.csv",
  row.names = FALSE
)


## =========================================================
## 3. Extract macrophages
## =========================================================
Idents(balmerged_singleR_obj) <- "Merged_Labels"
macrophage_obj <- subset(balmerged_singleR_obj, idents = "Macrophages")

if (!"group" %in% colnames(macrophage_obj@meta.data)) {
  stop("Error: 'group' column does not exist in meta.data.")
}

macrophage_obj$group <- factor(
  macrophage_obj$group,
  levels = c("Control", "COPD-quit", "COPD-active")
)

## ---------------------------------------------------------
## 3A. Reprocessing macrophage subset
## ---------------------------------------------------------
DefaultAssay(macrophage_obj) <- "RNA"

macrophage_obj <- NormalizeData(macrophage_obj)
macrophage_obj <- FindVariableFeatures(macrophage_obj, selection.method = "vst", nfeatures = 2000)
macrophage_obj <- ScaleData(macrophage_obj)
macrophage_obj <- RunPCA(macrophage_obj, npcs = 20)

p_elbow_mac <- ElbowPlot(macrophage_obj) + ggtitle("Elbow plot: macrophage subset")
save_plot("ElbowPlot_Macrophage.png", p_elbow_mac, 6, 4.5)

macrophage_obj <- FindNeighbors(macrophage_obj, dims = 1:15)
macrophage_obj <- FindClusters(macrophage_obj, resolution = 0.7)
macrophage_obj <- RunUMAP(macrophage_obj, dims = 1:15)

## ---------------------------------------------------------
## 3B. Macrophage UMAP by group
## ---------------------------------------------------------

macrophage_obj$group <- factor(
  macrophage_obj$group,
  levels = c("Control", "COPD-quit", "COPD-active")
)

p_mac_group <- DimPlot(
  macrophage_obj,
  reduction = "umap",
  group.by = "group",
  cols = group_colors,
  pt.size = 1
) +
  ggtitle("Macrophage UMAP by group") +
  theme_umap_paper()

p_mac_split <- DimPlot(
  macrophage_obj,
  reduction = "umap",
  group.by = "group",
  split.by = "group",
  cols = group_colors,
  pt.size = 0.8
) +
  ggtitle("Macrophage UMAP split by group") +
  theme_umap_paper() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12)
  )

save_plot_paper(p_mac_group, "Macrophage_UMAP_by_Group", width = 8, height = 6)
save_plot_paper(p_mac_split, "Macrophage_UMAP_Split", width = 10, height = 4)

## =========================================================
## 4. Marker QC for macrophage subset
## =========================================================

Idents(macrophage_obj) <- macrophage_obj$seurat_clusters

p_tcell <- VlnPlot(
  macrophage_obj,
  features = c("TRAC", "CD3D", "CD3E", "IL32", "MS4A1", "CD79A"),
  group.by = "seurat_clusters",
  pt.size = 0,
  stack = TRUE,
  flip = TRUE
) +
  ggtitle("T/B cell markers") +
  theme_violin_paper()

p_epi <- VlnPlot(
  macrophage_obj,
  features = c("SCGB1A1", "CFAP157", "WFDC2", "EPCAM", "KRT8", "KRT13",
               "SCEL", "SPRR3", "TMPRSS11B", "KRT4"),
  group.by = "seurat_clusters",
  pt.size = 0,
  stack = TRUE,
  flip = TRUE
) +
  ggtitle("Epithelial markers") +
  theme_violin_paper()

p_macro <- VlnPlot(
  macrophage_obj,
  features = c("CD68", "CSF1R", "MARCO", "LYZ", "MSR1", "MRC1"),
  group.by = "seurat_clusters",
  pt.size = 0,
  stack = TRUE,
  flip = TRUE
) +
  ggtitle("Macrophage markers") +
  theme_violin_paper()

save_plot_paper(p_tcell,  "Violin_TBcell_Markers_Macrophage_QC", width = 10, height = 8)
save_plot_paper(p_epi,    "Violin_Epithelial_Markers_Macrophage_QC", width = 10, height = 12)
save_plot_paper(p_macro,  "Violin_Macrophage_Markers_Macrophage_QC", width = 10, height = 8)

## ---------------------------------------------------------
## Marker discovery
## ---------------------------------------------------------
all_markers <- FindAllMarkers(
  macrophage_obj,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.1,
  return.thresh = 0.05
)

fc_col <- if ("avg_log2FC" %in% colnames(all_markers)) "avg_log2FC" else "avg_logFC"

matop10_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = .data[[fc_col]], n = 10)

write.csv(all_markers, "Macrophage_AllMarkers.csv", row.names = FALSE)
write.csv(matop10_markers, "Macrophage_Top10_Markers.csv", row.names = FALSE)

Idents(macrophage_obj) <- "seurat_clusters"

table(Idents(macrophage_obj))

## =========================================================
## 4D. DotPlot for QC markers across macrophage clusters
## =========================================================

qc_genes <- c(
  "MARCO", "PPARG", "FABP4", "CD68", 
  "CD14", "VCAN", "FCN1", "CCR2",
  "MPO", "FCGR3B", "CSF3R",
  "EPCAM", "KRT8", "KRT18", "SFTPA1", "SFTPB", "SFTPC", "SCGB1A1",
  "TRAC", "CD3D", "CD3E",
  "MS4A1", "CD79A",
  "HBB", "HBA1", "HBA2"
)

qc_genes_use <- qc_genes[qc_genes %in% rownames(macrophage_obj)]
print(qc_genes_use)


macrophage_obj$seurat_clusters <- factor(
  macrophage_obj$seurat_clusters,
  levels = sort(unique(as.numeric(as.character(macrophage_obj$seurat_clusters))))
)

Idents(macrophage_obj) <- macrophage_obj$seurat_clusters

p_dot <- DotPlot(
  macrophage_obj,
  features = qc_genes_use,
  group.by = "seurat_clusters",
  cols = c("white", "#D55E00"),
  dot.scale = 6
) +
  RotatedAxis() +
  labs(
    title = "QC marker expression across macrophage clusters",
    x = NULL,
    y = "Seurat clusters"
  ) +
  theme_dotplot_paper() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 10
    ),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

save_plot_paper(
  p_dot,
  "QC_marker_dotplot_macrophage_clusters",
  width = 15,
  height = 8,
  dpi = 600
)


## ---------------------------------------------------------
## 4A. Remove contaminating clusters
## ---------------------------------------------------------


Idents(macrophage_obj) <- macrophage_obj$seurat_clusters


clusters_to_remove <- c("9", "10", "11", "12", "14")

macrophage_obj <- subset(
  macrophage_obj,
  idents = clusters_to_remove,
  invert = TRUE
)


Idents(macrophage_obj) <- macrophage_obj$seurat_clusters


table(Idents(macrophage_obj))

## ---------------------------------------------------------
## 4B. QC metrics across retained macrophage clusters
## ---------------------------------------------------------
p_qc_cluster_metrics <- VlnPlot(
  macrophage_obj,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  group.by = "seurat_clusters",
  pt.size = 0,
  stack = FALSE
) +
  ggtitle("QC metrics across macrophage clusters") +
  theme_violin_paper()

save_plot_paper(
  p_qc_cluster_metrics,
  "QC_Macrophage_Cluster_Metrics",
  width = 10,
  height = 6
)

## ---------------------------------------------------------
## 4C. Feature check for contamination markers
## ---------------------------------------------------------
feature_check_genes <- c(
  "LYZ", "CSF1R", "CD68",
  "TRAC", "CD3D", "CD3E",
  "EPCAM", "KRT19", "SCGB1A1",
  "HBB"
)

feature_check_genes_use <- feature_check_genes[
  feature_check_genes %in% rownames(macrophage_obj)
]

p_feature_check <- FeaturePlot(
  macrophage_obj,
  features = feature_check_genes_use,
  reduction = "umap",
  ncol = 4,
  pt.size = 0.4,
  order = TRUE
) &
  theme_umap_paper()

save_plot_paper(
  p_feature_check,
  "FeaturePlot_Macrophage_Contamination_Check",
  width = 14,
  height = 10
)
## =========================================================
## 4B. Reprocess macrophage object after removing contaminants
## =========================================================
DefaultAssay(macrophage_obj) <- "RNA"
Idents(macrophage_obj) <- "seurat_clusters"

macrophage_obj <- NormalizeData(macrophage_obj, verbose = FALSE)
macrophage_obj <- FindVariableFeatures(macrophage_obj, verbose = FALSE)
macrophage_obj <- ScaleData(
  macrophage_obj,
  features = VariableFeatures(macrophage_obj),
  verbose = FALSE
)
macrophage_obj <- RunPCA(
  macrophage_obj,
  features = VariableFeatures(macrophage_obj),
  verbose = FALSE
)

macrophage_obj <- FindNeighbors(macrophage_obj, dims = 1:15, verbose = FALSE)
macrophage_obj <- FindClusters(macrophage_obj, resolution = 0.7, verbose = FALSE)
macrophage_obj <- RunUMAP(macrophage_obj, dims = 1:15, verbose = FALSE)

Idents(macrophage_obj) <- "seurat_clusters"


cluster_order <- sort(as.character(unique(Idents(macrophage_obj))))
cluster_order <- cluster_order[order(as.numeric(cluster_order))]

macrophage_obj$seurat_clusters <- factor(
  macrophage_obj$seurat_clusters,
  levels = cluster_order
)

Idents(macrophage_obj) <- "seurat_clusters"

print(levels(Idents(macrophage_obj)))
print(table(Idents(macrophage_obj)))

## =========================================================
## Contamination marker check
## =========================================================

feature_check_genes <- c(
  ## macrophage markers
  "LYZ", "MARCO", "CD68",
  
  ## T cell contaminants
  "TRAC", "CD3D", "CD3E",
  
  ## epithelial contaminants
  "EPCAM", "KRT19", "SCGB1A1"
)

feature_check_genes_use <- feature_check_genes[
  feature_check_genes %in% rownames(macrophage_obj)
]

feature_check_genes_missing <- setdiff(
  feature_check_genes,
  feature_check_genes_use
)

if (length(feature_check_genes_missing) > 0) {
  message(
    "Genes not found in macrophage_obj: ",
    paste(feature_check_genes_missing, collapse = ", ")
  )
}

if (length(feature_check_genes_use) > 0) {
  
  p_feature_check <- FeaturePlot(
    macrophage_obj,
    features = feature_check_genes_use,
    reduction = "umap",
    ncol = 3,
    pt.size = 2,
    order = TRUE,
    raster = TRUE,
    cols = c("grey85", "#0047CC"),
    min.cutoff = "q10",
    max.cutoff = "q95"
  ) &
    theme_umap_paper() &
    theme(
      plot.title = element_text(
        size = 12,
        face = "bold"
      ),
      
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.4
      )
    )
  
  save_plot_paper(
    p_feature_check,
    "FeaturePlot_Macrophage_Contamination_Check",
    width = 12,
    height = 10
  )
  
} else {
  message("No contamination check genes found in macrophage_obj.")
}
## =========================================================
## 4C. DotPlot of macrophage programs
## =========================================================
gene_list <- unique(c(
  # 1. Resident-like AMs
  "FABP4", "PPARG", "MARCO", "MCEMP1",
  
  # 2. Transitional Macrophages
  "SPP1", "CD9", "LGALS3", "VCAN",
  
  # 3. Monocyte-derived M1-like AMs
  "CD14", "FCN1", "CD86", "IL1B", "TREM1",
  
  # 4. Monocyte-derived M2-like AMs
  "CD163", "MRC1", "CCL18", "TGM2",
  
  # 5. Foamy Macrophages
  "PLIN3", "APOE", "TREM2", "GPNMB",
  
  # 6. Oxidative Foamy Macrophages
  "CYBB", "NCF2", "HMOX1", "G6PD",
  
  # 7. Proliferating Macrophages
  "MKI67", "TOP2A", "CDK1"
))

gene_list_use <- gene_list[gene_list %in% rownames(macrophage_obj)]
print(gene_list_use)

p_dot_mac_program <- DotPlot(
  macrophage_obj,
  features = gene_list_use,
  group.by = "seurat_clusters",
  cols = c("lightgrey", "#1F77B4"),
  dot.scale = 6
) +
  RotatedAxis() +
  ggtitle("Macrophage marker and functional program genes") +
  theme_dotplot_paper() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      face = "italic",
      size = 11
    ),
    axis.text.y = element_text(size = 11)
  )

save_plot_paper(
  p_dot_mac_program,
  "DotPlot_Macrophage_Programs",
  width = 14,
  height = 8
)

## =========================================================
## 4D. UMAP of reclustered macrophages
## =========================================================
p_umap_clusters <- DimPlot(
  macrophage_obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.4
) +
  ggtitle("UMAP of macrophage clusters") +
  theme_umap_paper()

save_plot_paper(
  p_umap_clusters,
  "UMAP_Macrophage_Clusters",
  width = 8,
  height = 6
)

## =========================================================
## 4E. FindAllMarkers
## =========================================================
Idents(macrophage_obj) <- "seurat_clusters"

markers_mac <- FindAllMarkers(
  object = macrophage_obj,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.2,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  return.thresh = 0.05
)

fc_col <- if ("avg_log2FC" %in% colnames(markers_mac)) "avg_log2FC" else "avg_logFC"

top20_markers <- markers_mac %>%
  group_by(cluster) %>%
  filter(pct.1 >= 0.25) %>%
  arrange(desc(.data[[fc_col]]), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()

write.csv(markers_mac, "FindAllMarkers_macrophage_all.csv", row.names = FALSE)
write.csv(top20_markers, "FindAllMarkers_macrophage_top20.csv", row.names = FALSE)
## =========================================================
## =========================================================
## 0. Libraries
## =========================================================
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(Matrix)

## =========================================================
## 1. Macrophage subtype annotation
## =========================================================
Idents(macrophage_obj) <- "seurat_clusters"

macrophage_obj$Macrophage_Type <- dplyr::recode(
  as.character(macrophage_obj$seurat_clusters),
  "0"  = "Oxidative Foamy Macrophages",
  "1"  = "Oxidative Foamy Macrophages",
  "2"  = "Monocyte-derived M2-like AMs",
  "3"  = "Foamy Macrophages",
  "4"  = "Transitional LAMs",
  "5"  = "Foamy Macrophages",
  "6"  = "Resident-like AMs",
  "7"  = "Oxidative Foamy Macrophages",
  "8"  = "Monocyte-derived M1-like AMs",
  "9"  = "Oxidative Foamy Macrophages",
  "10" = "Oxidative Foamy Macrophages",
  "11" = "Transitional LAMs",
  "12" = "Proliferating macrophages",
  "13" = "Monocyte-derived M1-like AMs"
)

mac_order <- c(
  "Proliferating macrophages",
  "Monocyte-derived M1-like AMs",
  "Resident-like AMs",
  "Transitional LAMs",
  "Monocyte-derived M2-like AMs",
  "Foamy Macrophages",
  "Oxidative Foamy Macrophages"
)

macrophage_obj$Macrophage_Type <- factor(
  macrophage_obj$Macrophage_Type,
  levels = mac_order
)

cat("Macrophage_Type counts:\n")
print(table(macrophage_obj$Macrophage_Type, useNA = "ifany"))

## =========================================================
## 2. Publication theme / colors
## =========================================================
group_colors <- c(
  "Control"     = "#A0A0A0",  
  "COPD-quit"   = "#E5C494",  
  "COPD-active" = "#B82E2E"   
)


mac_colors <- c(
  "Foamy Macrophages"            = "#4E79A7",  
  "Monocyte-derived M1-like AMs" = "#F28E2B",  
  "Monocyte-derived M2-like AMs" = "#59A14F",  
  "Oxidative Foamy Macrophages"  = "#D3A2C7",  
  "Proliferating macrophages"    = "#B07AA1",  
  "Resident-like AMs"            = "#9C755F",  
  "Transitional LAMs"            = "#76B7B2"   
)
paper_theme <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title      = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title      = element_text(size = 14, face = "bold"),
      axis.text.x     = element_text(size = 12, color = "black"),
      axis.text.y     = element_text(size = 12, color = "black"),
      legend.title    = element_text(size = 12, face = "bold"),
      legend.text     = element_text(size = 11, color = "black"),
      panel.grid      = element_blank(),
      panel.border    = element_rect(linewidth = 0.8, color = "black")
    )
}

theme_umap_paper <- function(base_size = 14) {
  paper_theme(base_size) +
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank()
    )
}

theme_dotplot_paper <- function(base_size = 14) {
  paper_theme(base_size) +
    theme(
      axis.text.x = element_text(
        angle = 90, vjust = 0.5, hjust = 1,
        size = 11, color = "black"
      ),
      axis.text.y = element_text(size = 12, color = "black")
    )
}

save_plot_paper <- function(plot, filename, width = 8, height = 6, dpi = 600) {
  ggsave(
    paste0(filename, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
  ggsave(
    paste0(filename, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
  ggsave(
    paste0(filename, ".tiff"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    compression = "lzw"
  )
}

## =========================================================
## 3. Prepare object for plotting
## =========================================================
DefaultAssay(macrophage_obj) <- "RNA"

# Seurat v5 merged-layer safety
macrophage_obj <- JoinLayers(macrophage_obj, assay = "RNA")

# keep only cells with valid subtype label
macrophage_obj_plot <- subset(
  macrophage_obj,
  cells = colnames(macrophage_obj)[!is.na(macrophage_obj$Macrophage_Type)]
)

macrophage_obj_plot$Macrophage_Type <- factor(
  as.character(macrophage_obj_plot$Macrophage_Type),
  levels = mac_order
)
macrophage_obj_plot$Macrophage_Type <- droplevels(macrophage_obj_plot$Macrophage_Type)

Idents(macrophage_obj_plot) <- "Macrophage_Type"

cols_use <- mac_colors[levels(macrophage_obj_plot$Macrophage_Type)]

## =========================================================
## 4. Marker genes
## =========================================================
gene_list <- unique(c(
  # 1. Resident-like AMs
  "FABP4", "PPARG", "MARCO", "MCEMP1",
  
  # 2. Transitional Macrophages
  "SPP1", "CD9", "LGALS3", "VCAN",
  
  # 3. Monocyte-derived M1-like AMs
  "CD14", "FCN1", "CD86", "IL1B", "TREM1",
  
  # 4. Monocyte-derived M2-like AMs
  "CD163", "MRC1", "CCL18", "TGM2",
  
  # 5. Foamy Macrophages
  "PLIN3", "APOE", "TREM2", "GPNMB",
  
  # 6. Oxidative Foamy Macrophages
  "CYBB", "NCF2", "HMOX1", "G6PD",
  
  # 7. Proliferating Macrophages
  "MKI67", "TOP2A", "CDK1"
))

expr_mat <- GetAssayData(macrophage_obj_plot, assay = "RNA", layer = "data")
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]

gene_list_filtered <- gene_list[gene_list %in% rownames(expr_mat)]
gene_list_filtered <- gene_list_filtered[
  Matrix::rowSums(expr_mat[gene_list_filtered, , drop = FALSE]) > 0
]

cat("Genes used in dot plot:\n")
print(gene_list_filtered)

# add line breaks every 3rd gene for readability
x_labels <- gene_list_filtered
if (length(x_labels) >= 3) {
  spacing_idx <- seq(3, length(x_labels), by = 3)
  x_labels[spacing_idx] <- paste0(x_labels[spacing_idx], "\n")
}

## =========================================================
## 5. Custom DotPlot (robust replacement of Seurat::DotPlot)
## =========================================================
meta_df <- macrophage_obj_plot@meta.data %>%
  dplyr::select(Macrophage_Type) %>%
  tibble::rownames_to_column("cell")

meta_df$Macrophage_Type <- factor(meta_df$Macrophage_Type, levels = mac_order)

# average expression per subtype
avg_exp <- sapply(levels(meta_df$Macrophage_Type), function(ct) {
  cells_use <- meta_df$cell[meta_df$Macrophage_Type == ct]
  if (length(cells_use) == 0) {
    rep(NA_real_, length(gene_list_filtered))
  } else {
    Matrix::rowMeans(expr_mat[gene_list_filtered, cells_use, drop = FALSE])
  }
})

avg_exp <- as.data.frame(avg_exp)
colnames(avg_exp) <- levels(meta_df$Macrophage_Type)
avg_exp$gene <- gene_list_filtered

avg_long <- avg_exp %>%
  pivot_longer(
    cols = -gene,
    names_to = "Macrophage_Type",
    values_to = "avg_exp"
  )

# percent expressing per subtype
pct_exp <- sapply(levels(meta_df$Macrophage_Type), function(ct) {
  cells_use <- meta_df$cell[meta_df$Macrophage_Type == ct]
  if (length(cells_use) == 0) {
    rep(NA_real_, length(gene_list_filtered))
  } else {
    Matrix::rowMeans(expr_mat[gene_list_filtered, cells_use, drop = FALSE] > 0) * 100
  }
})

pct_exp <- as.data.frame(pct_exp)
colnames(pct_exp) <- levels(meta_df$Macrophage_Type)
pct_exp$gene <- gene_list_filtered

pct_long <- pct_exp %>%
  pivot_longer(
    cols = -gene,
    names_to = "Macrophage_Type",
    values_to = "pct_exp"
  )

dot_df <- avg_long %>%
  left_join(pct_long, by = c("gene", "Macrophage_Type")) %>%
  group_by(gene) %>%
  mutate(avg_exp_scaled = as.numeric(scale(avg_exp))) %>%
  ungroup()

# replace NA scale values (e.g. zero variance genes within all groups)
dot_df$avg_exp_scaled[is.na(dot_df$avg_exp_scaled)] <- 0

dot_df$gene <- factor(dot_df$gene, levels = gene_list_filtered)
dot_df$Macrophage_Type <- factor(dot_df$Macrophage_Type, levels = rev(mac_order))

p_dot_macrophage_markers <- ggplot(
  dot_df,
  aes(x = gene, y = Macrophage_Type)
) +
  geom_point(aes(size = pct_exp, color = avg_exp_scaled)) +
  scale_size(range = c(1, 8), name = "Pct. expressing") +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Scaled expression"
  ) +
  scale_x_discrete(labels = x_labels) +
  labs(
    title = "Marker expression across macrophage subtypes",
    x = NULL,
    y = NULL
  ) +
  theme_dotplot_paper() +
  theme(
    legend.position = "right"
  )

save_plot_paper(
  plot = p_dot_macrophage_markers,
  filename = "Fig_Macrophage_DotPlot_Markers",
  width = 12,
  height = 6,
  dpi = 600
)

## =========================================================
## 6. UMAP of macrophage subtypes
## =========================================================
p_mac_umap_pub <- DimPlot(
  macrophage_obj_plot,
  reduction = "umap",
  group.by = "Macrophage_Type",
  cols = cols_use,
  label = FALSE,
  pt.size = 0.45
) +
  ggtitle("Macrophage subtype annotation") +
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_umap_paper() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold", color = "black"),
    
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.key = element_blank()
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = 21,
        size = 4,
        fill = unname(cols_use),
        color = "black",
        stroke = 0.6
      ),
      ncol = 1
    )
  )

save_plot_paper(
  plot = p_mac_umap_pub,
  filename = "UMAP_Macrophage_Type_publication",
  width = 10,
  height = 7,
  dpi = 600
)

## =========================================================
## 6B. UMAP colored by group
## =========================================================

macrophage_obj_plot$group <- factor(
  macrophage_obj_plot$group,
  levels = c("Control", "COPD-quit", "COPD-active")
)

group_colors <- c(
  "Control"     = "#A0A0A0", 
  "COPD-quit"   = "#E5C494", 
  "COPD-active" = "#B82E2E"   
)


p_mac_umap_group <- DimPlot(
  macrophage_obj_plot,
  reduction = "umap",
  group.by = "group",
  cols = group_colors,
  label = FALSE,
  pt.size = 0.45
) +
  ggtitle("Macrophage UMAP by group") +
  labs(
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_umap_paper() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 13, face = "bold", color = "black"),
    
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.key = element_blank()
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = 21,
        size = 4,
        fill = unname(group_colors),
        color = "black",
        stroke = 0.6
      ),
      ncol = 1
    )
  )
save_plot_paper(
  plot = p_mac_umap_group,
  filename = "UMAP_Macrophage_Group_publication",
  width = 8,
  height = 6,
  dpi = 600
)
## =========================================================
## 7. Print plots in session
## =========================================================
print(p_dot_macrophage_markers)
print(p_mac_umap_pub)

macrophage_obj$group <- factor(
  macrophage_obj$group,
  levels = c("Control", "COPD-quit", "COPD-active")
)

## =========================================================
## 7A. Macrophage subtype UMAP split by group
## =========================================================
p_mac_split <- DimPlot(
  macrophage_obj,
  reduction = "umap",
  group.by = "Macrophage_Type",
  split.by = "group",
  label = FALSE,
  cols = cols_use,
  pt.size = 0.4
) +
  ggtitle("Macrophage subtype UMAP split by group") +
  theme_umap_paper() +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(size = 13)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        shape = 21,
        size = 4,
        fill = unname(cols_use),
        color = "black",
        stroke = 0.6
      )
    )
  )

save_plot_paper(
  p_mac_split,
  "Macrophage_UMAP_Split_by_Group",
  width = 14,
  height = 6
)

## =========================================================
## 7B. Macrophage subtype UMAP by group with shared legend
## =========================================================
library(patchwork)

groups <- c("Control", "COPD-quit", "COPD-active")

plot_list <- vector("list", length(groups))
names(plot_list) <- groups

for (i in seq_along(groups)) {
  
  grp <- groups[i]
  
  obj_sub <- subset(
    macrophage_obj,
    subset = group == grp
  )
  
  plot_list[[i]] <- DimPlot(
    obj_sub,
    reduction = "umap",
    group.by = "Macrophage_Type",
    label = TRUE,
    repel = TRUE,
    label.size = 4.5,
    cols = cols_use,
    pt.size = 0.45
  ) +
    ggtitle(grp) +
    theme_umap_paper() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      legend.text = element_text(size = 13)
    ) +
    guides(
      color = guide_legend(
        override.aes = list(
          shape = 21,
          size = 4,
          fill = unname(cols_use),
          color = "black",
          stroke = 0.6
        )
      )
    )
}

combined_plot <- (plot_list[[1]] | plot_list[[2]] | plot_list[[3]]) +
  patchwork::plot_layout(guides = "collect") &
  theme(
    legend.position = "right"
  )

save_plot_paper(
  combined_plot,
  "UMAP_Macrophage_ByGroup_SharedLegend",
  width = 21,
  height = 7
)
## =========================================================
## Macrophage subcluster proportion by group (VERTICAL)
## =========================================================

cell_counts_sub <- macrophage_obj@meta.data %>%
  transmute(
    group = factor(
      group,
      levels = c("Control", "COPD-quit", "COPD-active")
    ),
    subcluster = factor(
      as.character(Macrophage_Type),
      levels = mac_order
    )
  ) %>%
  filter(!is.na(group), !is.na(subcluster)) %>%
  dplyr::count(group, subcluster, name = "n")

p_subcluster_prop <- ggplot(
  cell_counts_sub,
  aes(x = group, y = n, fill = subcluster)
) +
  geom_col(
    position = "fill",
    width = 0.72,
    color = "white",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = mac_colors,
    name = "Macrophage subtype",
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Macrophage subtype proportion by group",
    x = NULL,
    y = "Percentage"
  ) +
  theme_bw(base_size = 12, base_family = "") +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    
    axis.text.x = element_text(size = 13, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 13, color = "black"),
    legend.key.size = unit(0.6, "cm"),
    legend.position = "right"
  )

save_plot_paper(
  p_subcluster_prop,
  "Macrophage_Subclusters_Proportion_by_Group_VERTICAL",
  width = 8,
  height = 6
)

## =========================================================
## Stacked horizontal barplot of macrophage subtype cell numbers
## =========================================================

library(dplyr)

bar_data <- macrophage_obj@meta.data %>%
  transmute(
    group = factor(
      group,
      levels = c("Control", "COPD-quit", "COPD-active")
    ),
    Macrophage_Type = factor(
      as.character(Macrophage_Type),
      levels = levels(macrophage_obj$Macrophage_Type)
    )
  ) %>%
  filter(!is.na(group), !is.na(Macrophage_Type)) %>%
  dplyr::count(Macrophage_Type, group, name = "CellNumber")

## publication-style group colors
group_colors <- c(
  "Control"     = "#A0A0A0",  
  "COPD-quit"   = "#E5C494",  
  "COPD-active" = "#B82E2E"  
)

bar_plot_h <- ggplot(
  bar_data,
  aes(
    x = Macrophage_Type,
    y = CellNumber,
    fill = group
  )
) +
  geom_col(
    width = 0.72,
    color = "white",
    linewidth = 0.25
  ) +
  coord_flip() +
  scale_fill_manual(
    values = group_colors,
    breaks = c("Control", "COPD-quit", "COPD-active")
  ) +
  
  labs(
    title = "Cell numbers of macrophage subtypes",
    x = "Macrophage subtype",
    y = "Cell number",
    fill = "Group"
  ) +
  
  theme_bw(base_size = 12) +
  
  theme(
    ## remove grids
    panel.grid = element_blank(),
    
    ## panel border
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    ## axis text
    axis.text.x = element_text(
      size = 12,
      face = "bold",
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 12,
      face = "bold",
      color = "black"
    ),
    
    ## axis titles
    axis.title.x = element_text(
      size = 14,
      face = "bold",
      color = "black"
    ),
    
    axis.title.y = element_text(
      size = 14,
      face = "bold",
      color = "black"
    ),
    
    ## title
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    
    ## legend
    legend.title = element_text(
      size = 13,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 13,
      color = "black"
    ),
    
    legend.position = "right",
    legend.key.size = unit(0.6, "cm")
  )

save_plot_paper(
  bar_plot_h,
  "BarChart_Macrophage_Subtypes_Stacked_Horizontal",
  width = 10,
  height = 6
)

saveRDS(
  balmerged_singleR_obj,
  "results/rds/BALF_scRNAseq_final_annotated.rds"
)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(rstatix)
  library(ggpubr)
  library(tibble)
  library(Matrix)
  library(grDevices)
})

# =========================================================
# =========================================================
FIG_WIDTH  <- 3.3
FIG_HEIGHT <- 3
FIG_DPI    <- 600

SAVE_PDF   <- TRUE
SAVE_PNG   <- TRUE
SAVE_TIFF  <- TRUE

TIFF_COMPRESSION <- "lzw"
BG_COLOR         <- "white"

theme_project_paper <- function(base_size = 15) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.7),
      axis.text.x = element_text(
        color = "black",
        size = 9,
        angle = 0,
        vjust = 1,
        hjust = 0.5
      ),
      axis.text.y = element_text(color = "black", size = 12),
      axis.title = element_text(face = "bold", color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 15, face = "bold", color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      legend.position = "none",
      plot.margin = margin(10, 15, 10, 10)
    )
}

theme_set(theme_project_paper(base_size = 15))

# =========================================================
# 1) Parameter settings
# =========================================================
genes_want <- c(
  "MARCO", "MRC1",
  "LPL", "LIPA", "PLIN2",
  "ABCG1", "CYBB",
  "TREM2", "GPX4", "IL1B", "TREM1"
)

group_col   <- "group"
group_order <- c("Control", "COPD-quit", "COPD-active")
outdir      <- "Violin_byGene_Group_paper"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

group_colors_all <- c(
  "Control"     = "#A0A0A0",  
  "COPD-quit"   = "#E5C494",  
  "COPD-active" = "#B82E2E"   
)



macrophage_obj@meta.data[[group_col]] <- as.character(macrophage_obj@meta.data[[group_col]])

macrophage_obj@meta.data[[group_col]] <- dplyr::recode(
  macrophage_obj@meta.data[[group_col]],
  "COPD_active" = "COPD-active",
  "COPD_quit"   = "COPD-quit",
  "COPD active" = "COPD-active",
  "COPD quit"   = "COPD-quit"
)

macrophage_obj@meta.data[[group_col]] <- factor(
  macrophage_obj@meta.data[[group_col]],
  levels = group_order
)

message("Group counts:")
print(table(macrophage_obj@meta.data[[group_col]], useNA = "ifany"))

# =========================================================
# =========================================================
expr_mat <- tryCatch(
  {
    GetAssayData(
      macrophage_obj,
      assay = DefaultAssay(macrophage_obj),
      layer = "data"
    )
  },
  error = function(e) {
    GetAssayData(
      macrophage_obj,
      assay = DefaultAssay(macrophage_obj),
      slot = "data"
    )
  }
)

genes_use <- intersect(genes_want, rownames(expr_mat))

if (length(genes_use) == 0) {
  stop("None of the genes in genes_want are present in the rownames of the current assay.")
}

genes_use <- genes_use[
  Matrix::rowSums(expr_mat[genes_use, , drop = FALSE]) > 0
]

if (length(genes_use) == 0) {
  stop("None of the candidate genes are expressed in the data/layer of the current assay.")
}

message("Genes to plot: ", paste(genes_use, collapse = ", "))

# =========================================================
# 4) Prepare data
# =========================================================
meta_df <- macrophage_obj@meta.data %>%
  rownames_to_column("cell") %>%
  dplyr::select(cell, all_of(group_col))

colnames(meta_df)[2] <- "group"

expr_df <- FetchData(macrophage_obj, vars = genes_use) %>%
  rownames_to_column("cell")

plot_df <- meta_df %>%
  left_join(expr_df, by = "cell") %>%
  mutate(group = as.character(group)) %>%
  filter(!is.na(group), group %in% group_order) %>%
  mutate(group = factor(group, levels = group_order))

group_present <- intersect(group_order, unique(as.character(plot_df$group)))

if (length(group_present) < 2) {
  stop("Fewer than two available groups; comparison cannot be performed.")
}

plot_df <- plot_df %>%
  mutate(group = factor(group, levels = group_present))

group_colors <- group_colors_all[group_present]

# =========================================================
# 5) Statistical annotation function
#    Use * for A > C and A > Q
#    Use # for Q > C
# =========================================================
make_custom_signif_label <- function(comparison, p_adj) {
  dplyr::case_when(
    comparison %in% c(
      "COPD-active vs Control",
      "COPD-active vs COPD-quit"
    ) & p_adj < 0.0001 ~ "****",
    
    comparison %in% c(
      "COPD-active vs Control",
      "COPD-active vs COPD-quit"
    ) & p_adj < 0.001 ~ "***",
    
    comparison %in% c(
      "COPD-active vs Control",
      "COPD-active vs COPD-quit"
    ) & p_adj < 0.01 ~ "**",
    
    comparison %in% c(
      "COPD-active vs Control",
      "COPD-active vs COPD-quit"
    ) & p_adj < 0.05 ~ "*",
    
    comparison == "COPD-quit vs Control" & p_adj < 0.0001 ~ "####",
    comparison == "COPD-quit vs Control" & p_adj < 0.001  ~ "###",
    comparison == "COPD-quit vs Control" & p_adj < 0.01   ~ "##",
    comparison == "COPD-quit vs Control" & p_adj < 0.05   ~ "#",
    
    TRUE ~ "ns"
  )
}

# =========================================================
# 6) Single-gene violin plot function
# =========================================================
plot_violin_gene_paper <- function(
    df,
    gene,
    p_adjust_method = "BH",
    group_colors,
    mode = c("auto", "clipped", "full"),
    base_size = 15
) {
  mode <- match.arg(mode)
  
  if (!(gene %in% colnames(df))) {
    stop("Gene column not found: ", gene)
  }
  
  dat <- df %>%
    dplyr::select(group, all_of(gene)) %>%
    dplyr::filter(!is.na(group))
  
  colnames(dat)[2] <- "expr"
  
  dat <- dat %>%
    dplyr::filter(!is.na(expr), is.finite(expr)) %>%
    dplyr::mutate(
      group = droplevels(factor(group, levels = names(group_colors)))
    )
  
  if (nlevels(dat$group) < 2) {
    stop("Gene ", gene, " has fewer than two valid groups.")
  }
  
  # -----------------------------
  # 1. Automatically determine the y-axis range
  # -----------------------------
  q02 <- as.numeric(stats::quantile(dat$expr, 0.02, na.rm = TRUE))
  q90 <- as.numeric(stats::quantile(dat$expr, 0.90, na.rm = TRUE))
  q98 <- as.numeric(stats::quantile(dat$expr, 0.98, na.rm = TRUE))
  q99 <- as.numeric(stats::quantile(dat$expr, 0.99, na.rm = TRUE))
  
  expr_min <- min(dat$expr, na.rm = TRUE)
  expr_max <- max(dat$expr, na.rm = TRUE)
  
  auto_clip <- is.finite(q90) &&
    is.finite(q99) &&
    q90 > 0 &&
    (q99 / q90 > 2.5)
  
  use_clip <- switch(
    mode,
    auto    = auto_clip,
    clipped = TRUE,
    full    = FALSE
  )
  
  if (use_clip) {
    y_min <- min(0, q02, na.rm = TRUE)
    y_max <- q98
  } else {
    y_min <- min(0, expr_min, na.rm = TRUE)
    y_max <- expr_max
  }
  
  forced_ymax <- c(
    "ABCA1" = 3.5,
    "PLIN2" = 3.5
  )
  
  if (gene %in% names(forced_ymax)) {
    y_max <- forced_ymax[[gene]]
    y_min <- min(y_min, 0)
  }
  
  if (!is.finite(y_min) || !is.finite(y_max) || y_max <= y_min) {
    y_min <- 0
    y_max <- max(1, expr_max, na.rm = TRUE)
  }
  
  y_rng <- y_max - y_min
  
  # -----------------------------
  # 2. Dunn’s test
  # -----------------------------
  dunn <- data.frame()
  
  unique_expr_n <- length(unique(dat$expr[is.finite(dat$expr)]))
  
  if (unique_expr_n > 1) {
    dunn <- dat %>%
      rstatix::dunn_test(
        expr ~ group,
        p.adjust.method = p_adjust_method
      ) %>%
      dplyr::mutate(
        group1 = as.character(group1),
        group2 = as.character(group2),
        
        Comparison = dplyr::case_when(
          group1 == "Control" & group2 == "COPD-active" ~ "COPD-active vs Control",
          group1 == "COPD-quit" & group2 == "COPD-active" ~ "COPD-active vs COPD-quit",
          group1 == "Control" & group2 == "COPD-quit" ~ "COPD-quit vs Control",
          TRUE ~ paste(group2, "vs", group1)
        ),
        
        Comparison = factor(
          Comparison,
          levels = c(
            "COPD-active vs Control",
            "COPD-active vs COPD-quit",
            "COPD-quit vs Control"
          )
        ),
        
        p.signif = make_custom_signif_label(
          comparison = as.character(Comparison),
          p_adj = p.adj
        )
      ) %>%
      dplyr::filter(!is.na(Comparison), p.adj < 0.05) %>%
      dplyr::arrange(Comparison)
  }
  
  # -----------------------------
  # 3. Fix the bracket height order
  # -----------------------------
  if (nrow(dunn) > 0) {
    dunn <- dunn %>%
      dplyr::mutate(
        comparison_id = dplyr::row_number(),
        y.position = y_max + y_rng * (0.14 + 0.12 * (comparison_id - 1))
      )
    
    y_upper <- max(dunn$y.position, na.rm = TRUE) + y_rng * 0.12
  } else {
    y_upper <- y_max + y_rng * 0.20
  }
  
  if (!is.finite(y_upper) || y_upper <= y_max) {
    y_upper <- y_max + y_rng * 0.25
  }
  
  # -----------------------------
  # 4. Generate plot
  # -----------------------------
  p <- ggplot(dat, aes(x = group, y = expr, fill = group)) +
    geom_violin(
      trim = TRUE,
      scale = "width",
      width = 0.78,
      linewidth = 0.5,
      color = "black",
      alpha = 0.95
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 21,
      size = 2.6,
      fill = "white",
      color = "black",
      stroke = 0.8
    ) +
    stat_summary(
      fun.data = function(x) {
        m <- mean(x, na.rm = TRUE)
        s <- stats::sd(x, na.rm = TRUE)
        data.frame(
          y = m,
          ymin = m - s,
          ymax = m + s
        )
      },
      geom = "errorbar",
      width = 0.12,
      linewidth = 0.55,
      color = "black"
    ) +
    scale_fill_manual(values = group_colors, drop = FALSE) +
    coord_cartesian(
      ylim = c(y_min, y_upper),
      clip = "off"
    ) +
    labs(
      title = gene,
      x = NULL,
      y = "Expression (log-normalized)"
    ) +
    theme_project_paper(base_size = base_size) +
    theme(
      plot.margin = margin(8, 8, 8, 8),
      legend.position = "none"
    )
  
  # -----------------------------
  # 5. Add statistical symbols
  # -----------------------------
  if (nrow(dunn) > 0) {
    p <- p +
      ggpubr::stat_pvalue_manual(
        dunn,
        label = "p.signif",
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        tip.length = 0.01,
        size = 5.2,
        bracket.size = 0.65,
        step.increase = 0,
        hide.ns = TRUE
      )
  }
  
  return(p)
}

# =========================================================
# 7) Plot saving function
# =========================================================
save_plot_paper <- function(
    plot_obj,
    filename_base,
    outdir,
    width = 3.3,
    height = 3.8,
    dpi = 600,
    save_pdf = TRUE,
    save_png = TRUE,
    save_tiff = TRUE,
    bg_color = "white",
    tiff_compression = "lzw"
) {
  if (!dir.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  }
  
  if (!is.numeric(width) || length(width) != 1 || is.na(width)) {
    stop("width must be a single numeric value, for example 3.3")
  }
  if (!is.numeric(height) || length(height) != 1 || is.na(height)) {
    stop("height must be a single numeric value, for example 3.8")
  }
  if (!is.numeric(dpi) || length(dpi) != 1 || is.na(dpi)) {
    stop("dpi must be a single numeric value, for example 600")
  }
  
  if (isTRUE(save_pdf)) {
    ggsave(
      filename = file.path(outdir, paste0(filename_base, ".pdf")),
      plot = plot_obj,
      width = width,
      height = height,
      units = "in",
      device = grDevices::cairo_pdf,
      bg = bg_color
    )
  }
  
  if (isTRUE(save_png)) {
    ggsave(
      filename = file.path(outdir, paste0(filename_base, ".png")),
      plot = plot_obj,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = bg_color
    )
  }
  
  if (isTRUE(save_tiff)) {
    ggsave(
      filename = file.path(outdir, paste0(filename_base, ".tiff")),
      plot = plot_obj,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      compression = tiff_compression,
      bg = bg_color
    )
  }
}

# =========================================================
# 8) Export MAIN violin plots
# =========================================================
failed_main <- c()

for (g in genes_use) {
  message("Processing: ", g)
  
  tryCatch({
    p_main <- plot_violin_gene_paper(
      df = plot_df,
      gene = g,
      p_adjust_method = "BH",
      group_colors = group_colors,
      mode = "auto"
    )
    
    save_plot_paper(
      plot_obj = p_main,
      filename_base = paste0("MAIN_", g),
      outdir = outdir,
      width = FIG_WIDTH,
      height = FIG_HEIGHT,
      dpi = FIG_DPI,
      save_pdf = SAVE_PDF,
      save_png = SAVE_PNG,
      save_tiff = SAVE_TIFF,
      bg_color = BG_COLOR,
      tiff_compression = TIFF_COMPRESSION
    )
  }, error = function(e) {
    message("Failed MAIN: ", g, " --> ", e$message)
    failed_main <<- c(failed_main, g)
  })
}

# =========================================================
# 9) Result summary
# =========================================================
message("Done. Output folder: ", normalizePath(outdir))

if (length(failed_main) > 0) {
  message("Failed MAIN genes: ", paste(failed_main, collapse = ", "))
} else {
  message("All MAIN plots saved successfully.")
}

# =========================================================
# 9) Result summary
# =========================================================
message("Done. Output folder: ", normalizePath(outdir))

if (length(failed_main) > 0) {
  message("Failed BOX genes: ", paste(failed_main, collapse = ", "))
} else {
  message("All BOX plots saved successfully.")
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(Matrix)
})

# =========================================================
# 0) GLOBAL SETTINGS
# =========================================================
options(stringsAsFactors = FALSE)

BASE_SIZE <- 14

group_order <- c("Control", "COPD-quit", "COPD-active")

mac_order <- c(
  "Proliferating macrophages",
  "Monocyte-derived M1-like AMs",
  "Resident-like AMs",
  "Transitional LAMs",
  "Monocyte-derived M2-like AMs",
  "Foamy Macrophages",
  "Oxidative Foamy Macrophages"
)

comp_colors <- c(
  "COPD-active vs Control"   = "#E64B35",
  "COPD-quit vs Control"     = "#4DBBD5",
  "COPD-active vs COPD-quit" = "#00A087"
)

rename_map <- c(
  "0"  = "Oxidative Foamy Macrophages",
  "1"  = "Oxidative Foamy Macrophages",
  "2"  = "Monocyte-derived M2-like AMs",
  "3"  = "Foamy Macrophages",
  "4"  = "Transitional LAMs",
  "5"  = "Foamy Macrophages",
  "6"  = "Resident-like AMs",
  "7"  = "Oxidative Foamy Macrophages",
  "8"  = "Monocyte-derived M1-like AMs",
  "9"  = "Oxidative Foamy Macrophages",
  "10" = "Oxidative Foamy Macrophages",
  "11" = "Transitional LAMs",
  "12" = "Proliferating macrophages",
  "13" = "Monocyte-derived M1-like AMs"
)

# =========================================================
# 1) PAPER THEME
# =========================================================
paper_theme <- function(base_size = BASE_SIZE) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(
        size = base_size + 2,
        face = "bold",
        hjust = 0.5
      ),
      axis.title = element_text(
        size = base_size,
        face = "bold"
      ),
      axis.text.x = element_text(
        size = base_size - 2,
        color = "black"
      ),
      axis.text.y = element_text(
        size = base_size - 2,
        color = "black"
      ),
      legend.title = element_text(
        size = base_size - 1,
        face = "bold"
      ),
      legend.text = element_text(
        size = base_size - 2
      ),
      strip.text = element_text(
        size = base_size - 1,
        face = "bold"
      ),
      panel.grid = element_blank(),
      panel.border = element_rect(
        linewidth = 0.8,
        color = "black"
      ),
      plot.margin = margin(5, 5, 5, 5)
    )
}

theme_set(paper_theme())

# =========================================================
# 2) SAVE FUNCTION
# =========================================================
save_plot_paper <- function(plot, filename, width, height, dpi = 600) {
  ggsave(
    paste0(filename, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf
  )
  ggsave(
    paste0(filename, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  ggsave(
    paste0(filename, ".tiff"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    compression = "lzw",
    bg = "white"
  )
}

# =========================================================
# 3) HELPER FUNCTIONS
# =========================================================
exclude_markers <- function(deg_df) {
  deg_df <- deg_df %>% rownames_to_column("gene")
  
  is_ig <- grepl("^IGH|^IGK|^IGL", deg_df$gene)
  
  sex_genes <- c("KDM5D", "DDX3Y", "UTY", "EIF1AY")
  is_sex <- deg_df$gene %in% sex_genes
  
  epi_patterns <- c(
    "^KRT",
    "^SCGB",
    "^MUC",
    "^BPIF",
    "^EPCAM$",
    "^CLDN",
    "^OCLN$",
    "^DSP$",
    "^DSG", "^DSC",
    "^FOXJ1$",
    "^TP63$",
    "^SPRR"
  )
  is_epi <- grepl(paste(epi_patterns, collapse = "|"), deg_df$gene)
  
  lymph_patterns <- c(
    "^CD3D$", "^CD3E$", "^CD3G$",
    "^TRAC$", "^TRBC",
    "^CD4$", "^CD8A$", "^CD8B$",
    "^MS4A1$", "^CD79A$", "^CD79B$",
    "^BANK1$",
    "^NKG7$", "^GNLY$", "^KLRD1$", "^KLRB1$", "^KLRK1$",
    "^GZMB$", "^GZMH$", "^PRF1$"
  )
  is_lymph <- grepl(paste(lymph_patterns, collapse = "|"), deg_df$gene)
  
  is_rbc <- grepl("^HB[ABDEGMQ]", deg_df$gene)
  
  deg_df %>%
    filter(!(is_ig | is_sex | is_epi | is_lymph | is_rbc)) %>%
    column_to_rownames("gene")
}

get_top_genes <- function(deg_results, n = 10) {
  deg_results <- deg_results %>% rownames_to_column("gene")
  
  top_up <- deg_results %>%
    filter(avg_log2FC > 0) %>%
    arrange(p_val_adj, desc(avg_log2FC)) %>%
    slice_head(n = n) %>%
    pull(gene)
  
  top_down <- deg_results %>%
    filter(avg_log2FC < 0) %>%
    arrange(p_val_adj, avg_log2FC) %>%
    slice_head(n = n) %>%
    pull(gene)
  
  list(up = top_up, down = top_down)
}

make_deg_df <- function(genes, deg_df, comparison_label) {
  genes <- unique(genes)
  genes <- genes[genes %in% rownames(deg_df)]
  
  data.frame(
    Gene = genes,
    Comparison = comparison_label,
    logFC = deg_df[genes, "avg_log2FC"],
    stringsAsFactors = FALSE
  )
}

# =========================================================
# 4) SUBSET MACROPHAGES
# =========================================================
Idents(balmerged_singleR_obj) <- "Merged_Labels"
macrophage_all <- subset(balmerged_singleR_obj, idents = "Macrophages")

Idents(macrophage_all) <- macrophage_all$seurat_clusters
rename_map_use <- rename_map[names(rename_map) %in% levels(Idents(macrophage_all))]
macrophage_all <- RenameIdents(macrophage_all, rename_map_use)

macrophage_all$Macrophage_Type <- as.character(Idents(macrophage_all))

macrophage_obj <- subset(
  macrophage_all,
  subset = Macrophage_Type %in% mac_order
)

macrophage_obj$Macrophage_Type <- factor(
  macrophage_obj$Macrophage_Type,
  levels = mac_order
)

cat("\n=== Macrophage subtype counts ===\n")
print(table(macrophage_obj$Macrophage_Type, useNA = "ifany"))

cat("\n=== Group counts ===\n")
print(table(macrophage_obj$group, useNA = "ifany"))

# =========================================================
# 5) DIFFERENTIAL EXPRESSION BY GROUP
# =========================================================
macrophage_obj$group <- factor(
  as.character(macrophage_obj$group),
  levels = group_order
)
Idents(macrophage_obj) <- "group"

if (!all(group_order %in% levels(Idents(macrophage_obj)))) {
  stop("group levels are missing Control / COPD-quit / COPD-active. Please check macrophage_obj$group first.")
}

de_genes_active <- FindMarkers(
  macrophage_obj,
  ident.1 = "COPD-active",
  ident.2 = "Control",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  assay = "RNA"
)

de_genes_quit <- FindMarkers(
  macrophage_obj,
  ident.1 = "COPD-quit",
  ident.2 = "Control",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  assay = "RNA"
)

de_genes_active_quit <- FindMarkers(
  macrophage_obj,
  ident.1 = "COPD-active",
  ident.2 = "COPD-quit",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  assay = "RNA"
)

cat("\n=== DEG dimensions ===\n")
print(dim(de_genes_active))
print(dim(de_genes_quit))
print(dim(de_genes_active_quit))

# =========================================================
# 6) FILTER UNWANTED GENES
# =========================================================
de_genes_active      <- exclude_markers(de_genes_active)
de_genes_quit        <- exclude_markers(de_genes_quit)
de_genes_active_quit <- exclude_markers(de_genes_active_quit)

# =========================================================
# 7) EXTRACT TOP GENES
# =========================================================
top10_active      <- get_top_genes(de_genes_active, 10)
top10_quit        <- get_top_genes(de_genes_quit, 10)
top10_active_quit <- get_top_genes(de_genes_active_quit, 10)

top10_active_genes      <- c(top10_active$up, top10_active$down)
top10_quit_genes        <- c(top10_quit$up, top10_quit$down)
top10_active_quit_genes <- c(top10_active_quit$up, top10_active_quit$down)

cat("\n=== Top gene vector lengths ===\n")
print(length(top10_active_genes))
print(length(top10_quit_genes))
print(length(top10_active_quit_genes))

# =========================================================
# 8) BUILD BARPLOT DATA
# =========================================================
top_genes_df <- bind_rows(
  make_deg_df(top10_active_genes,      de_genes_active,      "COPD-active vs Control"),
  make_deg_df(top10_quit_genes,        de_genes_quit,        "COPD-quit vs Control"),
  make_deg_df(top10_active_quit_genes, de_genes_active_quit, "COPD-active vs COPD-quit")
) %>%
  filter(!is.na(logFC), is.finite(logFC)) %>%
  mutate(
    Comparison = factor(
      Comparison,
      levels = c(
        "COPD-active vs Control",
        "COPD-quit vs Control",
        "COPD-active vs COPD-quit"
      )
    )
  ) %>%
  group_by(Comparison) %>%
  mutate(Gene_plot = factor(Gene, levels = Gene[order(logFC)])) %>%
  ungroup()

# =========================================================
# 9) PLOT: FACETED BARPLOT
# =========================================================
p_bar <- ggplot(top_genes_df, aes(x = Gene_plot, y = logFC, fill = Comparison)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.2) +
  coord_flip() +
  facet_wrap(~Comparison, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = comp_colors) +
  labs(
    title = "Top differentially expressed genes in macrophages",
    x = NULL,
    y = "Average log2 fold change"
  ) +
  paper_theme() +
  theme(
    axis.text.y = element_text(
      face = "italic",
      size = BASE_SIZE - 4
    ),
    axis.text.x = element_text(
      size = BASE_SIZE - 3
    ),
    legend.position = "none"
  )

save_plot_paper(
  plot = p_bar,
  filename = "Fig_Top10_DEGs_Faceted_filtered",
  width = 12,
  height = 5,
  dpi = 600
)


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(Matrix)
  library(ggrepel)
})

# =========================================================
# 0) GLOBAL SETTINGS
# =========================================================
options(stringsAsFactors = FALSE)

BASE_SIZE <- 14

group_order <- c("Control", "COPD-quit", "COPD-active")

mac_order <- c(
  "Proliferating macrophages",
  "Monocyte-derived M1-like AMs",
  "Resident-like AMs",
  "Transitional LAMs",
  "Monocyte-derived M2-like AMs",
  "Foamy Macrophages",
  "Oxidative Foamy Macrophages"
)

comp_levels <- c(
  "COPD-active vs Control",
  "COPD-quit vs Control",
  "COPD-active vs COPD-quit"
)

comp_colors <- c(
  "COPD-active vs Control"   = "#E64B35",
  "COPD-quit vs Control"     = "#4DBBD5",
  "COPD-active vs COPD-quit" = "#00A087"
)

rename_map <- c(
  "0"  = "Oxidative Foamy Macrophages",
  "1"  = "Oxidative Foamy Macrophages",
  "2"  = "Monocyte-derived M2-like AMs",
  "3"  = "Foamy Macrophages",
  "4"  = "Transitional LAMs",
  "5"  = "Foamy Macrophages",
  "6"  = "Resident-like AMs",
  "7"  = "Oxidative Foamy Macrophages",
  "8"  = "Monocyte-derived M1-like AMs",
  "9"  = "Oxidative Foamy Macrophages",
  "10" = "Oxidative Foamy Macrophages",
  "11" = "Transitional LAMs",
  "12" = "Proliferating macrophages",
  "13" = "Monocyte-derived M1-like AMs"
)

# =========================================================
# 1) UNIFIED PAPER THEME
# =========================================================
paper_theme <- function(base_size = BASE_SIZE) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(
        size = base_size + 2,
        face = "bold",
        hjust = 0.5
      ),
      axis.title = element_text(
        size = base_size,
        face = "bold"
      ),
      axis.text.x = element_text(
        size = base_size - 2,
        color = "black"
      ),
      axis.text.y = element_text(
        size = base_size - 2,
        color = "black"
      ),
      legend.title = element_text(
        size = base_size - 1,
        face = "bold"
      ),
      legend.text = element_text(
        size = base_size - 2
      ),
      strip.text = element_text(
        size = base_size - 1,
        face = "bold"
      ),
      panel.grid = element_blank(),
      panel.border = element_rect(
        linewidth = 0.8,
        color = "black"
      ),
      plot.margin = margin(6, 8, 6, 8)
    )
}

theme_set(paper_theme())

# =========================================================
# 2) UNIFIED SAVE FUNCTION
# =========================================================
save_plot_paper <- function(plot, filename, width, height, dpi = 600) {
  ggsave(
    paste0(filename, ".pdf"),
    plot = plot,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white"
  )
  ggsave(
    paste0(filename, ".png"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  ggsave(
    paste0(filename, ".tiff"),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    compression = "lzw",
    bg = "white"
  )
}

# =========================================================
# 3) HELPER FUNCTIONS
# =========================================================
exclude_markers <- function(deg_df) {
  deg_df <- deg_df %>% tibble::rownames_to_column("gene")
  
  is_ig <- grepl("^IGH|^IGK|^IGL", deg_df$gene)
  
  sex_genes <- c("KDM5D", "DDX3Y", "UTY", "EIF1AY")
  is_sex <- deg_df$gene %in% sex_genes
  
  epi_patterns <- c(
    "^KRT",
    "^SCGB",
    "^MUC",
    "^BPIF",
    "^EPCAM$",
    "^CLDN",
    "^OCLN$",
    "^DSP$",
    "^DSG", "^DSC",
    "^FOXJ1$",
    "^TP63$",
    "^SPRR"
  )
  is_epi <- grepl(paste(epi_patterns, collapse = "|"), deg_df$gene)
  
  lymph_patterns <- c(
    "^CD3D$", "^CD3E$", "^CD3G$",
    "^TRAC$", "^TRBC",
    "^CD4$", "^CD8A$", "^CD8B$",
    "^MS4A1$", "^CD79A$", "^CD79B$",
    "^BANK1$",
    "^NKG7$", "^GNLY$", "^KLRD1$", "^KLRB1$", "^KLRK1$",
    "^GZMB$", "^GZMH$", "^PRF1$"
  )
  is_lymph <- grepl(paste(lymph_patterns, collapse = "|"), deg_df$gene)
  
  is_rbc <- grepl("^HB[ABDEGMQ]", deg_df$gene)
  
  deg_df %>%
    filter(!(is_ig | is_sex | is_epi | is_lymph | is_rbc)) %>%
    tibble::column_to_rownames("gene")
}

get_top_genes <- function(deg_results, n = 10) {
  deg_results <- deg_results %>% tibble::rownames_to_column("gene")
  
  top_up <- deg_results %>%
    filter(avg_log2FC > 0) %>%
    arrange(p_val_adj, desc(avg_log2FC)) %>%
    slice_head(n = n) %>%
    pull(gene)
  
  top_down <- deg_results %>%
    filter(avg_log2FC < 0) %>%
    arrange(p_val_adj, avg_log2FC) %>%
    slice_head(n = n) %>%
    pull(gene)
  
  list(up = top_up, down = top_down)
}

make_deg_df <- function(genes, deg_df, comparison_label) {
  genes <- unique(genes)
  genes <- genes[genes %in% rownames(deg_df)]
  
  data.frame(
    Gene = genes,
    Comparison = comparison_label,
    logFC = deg_df[genes, "avg_log2FC"],
    stringsAsFactors = FALSE
  )
}

draw_volcano_pub_auto <- function(
    deg,
    title,
    lfc_cut = 0.5,
    padj_cut = 0.05,
    n_up = 5,
    n_down = 5,
    highlight_genes = NULL,
    x_q = 0.995,
    axis_expand = 0.08,
    y_max = 400,
    base_size = BASE_SIZE
) {
  df <- deg %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene")
  
  lfc_col <- dplyr::case_when(
    "avg_log2FC" %in% colnames(df) ~ "avg_log2FC",
    "avg_logFC"  %in% colnames(df) ~ "avg_logFC",
    TRUE ~ NA_character_
  )
  
  p_col <- dplyr::case_when(
    "p_val_adj" %in% colnames(df) ~ "p_val_adj",
    "p_adj"     %in% colnames(df) ~ "p_adj",
    TRUE ~ NA_character_
  )
  
  if (is.na(lfc_col) || is.na(p_col)) {
    stop("Cannot find logFC or adjusted p-value column in DEG table.")
  }
  
  if (is.null(highlight_genes)) {
    highlight_genes <- character(0)
  } else {
    highlight_genes <- as.character(highlight_genes)
  }
  
  df <- df %>%
    mutate(
      logFC = as.numeric(.data[[lfc_col]]),
      p_adj = as.numeric(.data[[p_col]]),
      p_adj = ifelse(is.na(p_adj), 1, p_adj),
      p_adj_plot = pmax(p_adj, 1e-300),
      neglog10_p_raw = -log10(p_adj_plot),
      neglog10_p = neglog10_p_raw,
      sig = case_when(
        p_adj < padj_cut & logFC >=  lfc_cut ~ "Up",
        p_adj < padj_cut & logFC <= -lfc_cut ~ "Down",
        TRUE ~ "NS"
      ),
      rank_score = abs(logFC) * neglog10_p_raw,
      is_manual = gene %in% highlight_genes
    )
  
  auto_up <- df %>%
    filter(sig == "Up", !is_manual) %>%
    arrange(desc(rank_score)) %>%
    slice_head(n = n_up)
  
  auto_down <- df %>%
    filter(sig == "Down", !is_manual) %>%
    arrange(desc(rank_score)) %>%
    slice_head(n = n_down)
  
  auto_genes <- c(auto_up$gene, auto_down$gene)
  
  df <- df %>%
    mutate(
      is_key = gene %in% c(highlight_genes, auto_genes),
      plot_class = case_when(
        is_key & sig == "Up"   ~ "KeyUp",
        is_key & sig == "Down" ~ "KeyDown",
        is_key & sig == "NS"   ~ "KeyNS",
        TRUE ~ sig
      )
    )
  
  x_abs <- abs(df$logFC[is.finite(df$logFC)])
  x_lim <- unname(stats::quantile(x_abs, probs = x_q, na.rm = TRUE))
  if (!is.finite(x_lim) || x_lim <= 0) x_lim <- max(x_abs, na.rm = TRUE)
  if (!is.finite(x_lim) || x_lim <= 0) x_lim <- 1
  x_lim <- x_lim * (1 + axis_expand)
  
  x_pretty <- pretty(c(-x_lim, x_lim), n = 5)
  x_lim_final <- max(abs(x_pretty), na.rm = TRUE)
  xlim_use <- c(-x_lim_final, x_lim_final)
  x_breaks <- pretty(xlim_use, n = 5)
  
  lab_df <- df %>%
    filter(is_key) %>%
    distinct(gene, .keep_all = TRUE)
  
  ggplot(df, aes(x = logFC, y = neglog10_p)) +
    geom_point(
      data = df %>% filter(!is_key, neglog10_p <= y_max),
      aes(color = plot_class),
      size = 0.6,
      alpha = 0.6
    ) +
    geom_point(
      data = df %>% filter(is_key, neglog10_p <= y_max),
      aes(color = plot_class),
      size = 1.4,
      alpha = 0.95
    ) +
    geom_vline(
      xintercept = c(-lfc_cut, lfc_cut),
      linetype = "dashed",
      linewidth = 0.5
    ) +
    geom_hline(
      yintercept = -log10(padj_cut),
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggrepel::geom_text_repel(
      data = lab_df %>% filter(neglog10_p <= y_max),
      aes(label = gene),
      size = 2.4,
      box.padding = 0.45,
      point.padding = 0.22,
      force = 1.5,
      force_pull = 0.6,
      segment.size = 0.2,
      segment.alpha = 0.8,
      min.segment.length = 0,
      max.overlaps = 40,
      seed = 123
    ) +
    scale_color_manual(
      values = c(
        Up      = "#D55E00",
        Down    = "#0072B2",
        NS      = "grey82",
        KeyUp   = "#B2182B",
        KeyDown = "#1B7837",
        KeyNS   = "#7B3294"
      )
    ) +
    scale_x_continuous(
      limits = xlim_use,
      breaks = x_breaks
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      breaks = seq(0, y_max, 100),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = title,
      x = "Average log2 fold change",
      y = expression(-log[10]("adjusted P value"))
    ) +
    paper_theme(base_size = base_size) +
    theme(
      legend.position = "none",
      axis.line = element_line(linewidth = 0.7, color = "black"),
      axis.ticks = element_line(linewidth = 0.7, color = "black"),
      plot.margin = margin(8, 24, 8, 8)
    )
}

# =========================================================
# 4) SUBSET MACROPHAGES
# =========================================================
Idents(balmerged_singleR_obj) <- "Merged_Labels"
macrophage_all <- subset(balmerged_singleR_obj, idents = "Macrophages")

Idents(macrophage_all) <- macrophage_all$seurat_clusters
rename_map_use <- rename_map[names(rename_map) %in% levels(Idents(macrophage_all))]
macrophage_all <- RenameIdents(macrophage_all, rename_map_use)

macrophage_all$Macrophage_Type <- as.character(Idents(macrophage_all))

macrophage_obj <- subset(
  macrophage_all,
  subset = Macrophage_Type %in% mac_order
)

macrophage_obj$Macrophage_Type <- factor(
  macrophage_obj$Macrophage_Type,
  levels = mac_order
)

cat("\n=== Macrophage subtype counts ===\n")
print(table(macrophage_obj$Macrophage_Type, useNA = "ifany"))

cat("\n=== Group counts ===\n")
print(table(macrophage_obj$group, useNA = "ifany"))

# =========================================================
# 5) DIFFERENTIAL EXPRESSION BY GROUP
# =========================================================
macrophage_obj$group <- factor(
  as.character(macrophage_obj$group),
  levels = group_order
)
Idents(macrophage_obj) <- "group"

if (!all(group_order %in% levels(Idents(macrophage_obj)))) {
  stop("group levels are missing Control / COPD-quit / COPD-active. Please check macrophage_obj$group first.")
}

de_genes_active <- FindMarkers(
  macrophage_obj,
  ident.1 = "COPD-active",
  ident.2 = "Control",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  assay = "RNA"
)

de_genes_quit <- FindMarkers(
  macrophage_obj,
  ident.1 = "COPD-quit",
  ident.2 = "Control",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  assay = "RNA"
)

de_genes_active_quit <- FindMarkers(
  macrophage_obj,
  ident.1 = "COPD-active",
  ident.2 = "COPD-quit",
  logfc.threshold = 0.25,
  min.pct = 0.1,
  assay = "RNA"
)

cat("\n=== DEG dimensions ===\n")
print(dim(de_genes_active))
print(dim(de_genes_quit))
print(dim(de_genes_active_quit))

# =========================================================
# 6) FILTER UNWANTED GENES
# =========================================================
de_genes_active      <- exclude_markers(de_genes_active)
de_genes_quit        <- exclude_markers(de_genes_quit)
de_genes_active_quit <- exclude_markers(de_genes_active_quit)

# =========================================================
# 7) TOP GENES FOR BARPLOT
# =========================================================
top10_active      <- get_top_genes(de_genes_active, 10)
top10_quit        <- get_top_genes(de_genes_quit, 10)
top10_active_quit <- get_top_genes(de_genes_active_quit, 10)

top10_active_genes      <- c(top10_active$up, top10_active$down)
top10_quit_genes        <- c(top10_quit$up, top10_quit$down)
top10_active_quit_genes <- c(top10_active_quit$up, top10_active_quit$down)

cat("\n=== Top gene vector lengths ===\n")
print(length(top10_active_genes))
print(length(top10_quit_genes))
print(length(top10_active_quit_genes))

# =========================================================
# 8) BUILD BARPLOT DATA
# =========================================================
top_genes_df <- bind_rows(
  make_deg_df(top10_active_genes,      de_genes_active,      "COPD-active vs Control"),
  make_deg_df(top10_quit_genes,        de_genes_quit,        "COPD-quit vs Control"),
  make_deg_df(top10_active_quit_genes, de_genes_active_quit, "COPD-active vs COPD-quit")
) %>%
  filter(!is.na(logFC), is.finite(logFC)) %>%
  mutate(
    Comparison = factor(Comparison, levels = comp_levels)
  ) %>%
  group_by(Comparison) %>%
  mutate(Gene_plot = factor(Gene, levels = Gene[order(logFC)])) %>%
  ungroup()

# =========================================================
# 9) PLOT: FACETED BARPLOT
# =========================================================
p_bar <- ggplot(top_genes_df, aes(x = Gene_plot, y = logFC, fill = Comparison)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.2) +
  coord_flip() +
  facet_wrap(~Comparison, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = comp_colors) +
  labs(
    title = "Top differentially expressed genes in macrophages",
    x = NULL,
    y = "Average log2 fold change"
  ) +
  paper_theme() +
  theme(
    axis.text.y = element_text(
      face = "italic",
      size = BASE_SIZE - 4
    ),
    axis.text.x = element_text(
      size = BASE_SIZE - 3
    ),
    legend.position = "none"
  )

save_plot_paper(
  plot = p_bar,
  filename = "Fig_Top10_DEGs_Faceted_filtered",
  width = 12,
  height = 5,
  dpi = 600
)

# =========================================================
# 10) OPTIONAL HIGHLIGHT GENES FOR VOLCANO
#    define these before plotting if needed
# =========================================================
# key_genes_active_vs_ctrl  <- c("GENE1", "GENE2")
# key_genes_quit_vs_ctrl    <- c("GENE3", "GENE4")
# key_genes_active_vs_quit  <- c("GENE5", "GENE6")

if (!exists("key_genes_active_vs_ctrl"))  key_genes_active_vs_ctrl  <- NULL
if (!exists("key_genes_quit_vs_ctrl"))    key_genes_quit_vs_ctrl    <- NULL
if (!exists("key_genes_active_vs_quit"))  key_genes_active_vs_quit  <- NULL

# =========================================================
# 11) VOLCANO PLOTS
# =========================================================
p_volcano_active_vs_ctrl <- draw_volcano_pub_auto(
  deg = de_genes_active,
  title = "COPD-active vs Control",
  highlight_genes = key_genes_active_vs_ctrl,
  y_max = 400,
  base_size = 10
)

p_volcano_quit_vs_ctrl <- draw_volcano_pub_auto(
  deg = de_genes_quit,
  title = "COPD-quit vs Control",
  highlight_genes = key_genes_quit_vs_ctrl,
  y_max = 400,
  base_size = 10
)

p_volcano_active_vs_quit <- draw_volcano_pub_auto(
  deg = de_genes_active_quit,
  title = "COPD-active vs COPD-quit",
  highlight_genes = key_genes_active_vs_quit,
  y_max = 400,
  base_size = 10
)

save_plot_paper(
  plot = p_volcano_active_vs_ctrl,
  filename = "Fig_volcano_active_vs_control",
  width = 3.35,
  height = 3.2,
  dpi = 600
)

save_plot_paper(
  plot = p_volcano_quit_vs_ctrl,
  filename = "Fig_volcano_quit_vs_control",
  width = 3.35,
  height = 3.2,
  dpi = 600
)

save_plot_paper(
  plot = p_volcano_active_vs_quit,
  filename = "Fig_volcano_active_vs_quit",
  width = 3.35,
  height = 3.2,
  dpi = 600
)

##--------------------------------------------------------------------------------------
# GSEA (Gene Set Enrichment Analysis) dot plots using clusterProfiler

library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# Check group info
table(Idents(macrophage_obj))
table(macrophage_obj@meta.data$group)
Idents(macrophage_obj) <- macrophage_obj@meta.data$group

# DEG (A vs Q)
# Compute DEGs
de_genes <- FindMarkers(macrophage_obj, ident.1 = "COPD-active", ident.2 = "COPD-quit", logfc.threshold = 0, min.pct = 0.1)

# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(rownames(de_genes), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_entrez <- gene_entrez[!duplicated(gene_entrez$SYMBOL), ]

# Merge with DEG data
de_genes$gene <- rownames(de_genes)
gene_list <- merge(gene_entrez, de_genes, by.x = "SYMBOL", by.y = "gene")

# Prepare gene ranking for GSEA
gene_ranking <- sort(jitter(gene_list$avg_log2FC), decreasing = TRUE)
names(gene_ranking) <- gene_list$ENTREZID
gene_ranking <- gene_ranking[!is.na(gene_ranking) & gene_ranking != 0]

# Run KEGG GSEA
gsea_kegg <- gseKEGG(geneList = gene_ranking, organism = "hsa", pvalueCutoff = 0.05, eps = 0)

# Run GO BP GSEA
gsea_go <- gseGO(geneList = gene_ranking, 
                 OrgDb = org.Hs.eg.db, 
                 ont = "BP", 
                 pvalueCutoff = 0.05, 
                 eps = 0,
                 nPermSimple = 10000,
                 minGSSize = 10,
                 maxGSSize = 500)

# Plot GSEA results
dotplot(gsea_go, showCategory = 10) +
  ggtitle("GO Biological Process GSEA (COPD-active vs COPD-quit)")

dotplot(gsea_kegg, showCategory = 10) +
  ggtitle("KEGG Pathway GSEA (COPD-active vs COPD-quit)")

# Custom dot plot
gsea_results <- as.data.frame(gsea_kegg@result)
top_pathways <- gsea_results %>%
  filter(p.adjust < 0.05) %>%
  arrange(NES) %>%
  top_n(15, wt = abs(NES)) %>%
  mutate(Count = setSize)

ggplot(top_pathways, aes(x = NES, y = reorder(Description, NES))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal() +
  xlab("Normalized Enrichment Score") +
  ylab("Pathway") +
  ggtitle("GSEA Enrichment Analysis (COPD-active vs COPD-quit)") +
  theme(legend.position = "right")

ggsave("A_Q_GSEA_Enrichment_Analysis.png", width = 8, height = 6, dpi = 300)
##--------------------------------------------------------------------------------------
# GSEA Analysis: COPD-quit vs Control (Q vs C)

# Compute DEGs
de_genes <- FindMarkers(macrophage_obj, ident.1 = "COPD-quit", ident.2 = "Control", logfc.threshold = 0, min.pct = 0.1)

# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(rownames(de_genes), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_entrez <- gene_entrez[!duplicated(gene_entrez$SYMBOL), ]

# Merge with DEG data
de_genes$gene <- rownames(de_genes)
gene_list <- merge(gene_entrez, de_genes, by.x = "SYMBOL", by.y = "gene")

# Prepare gene ranking for GSEA
gene_ranking <- sort(jitter(gene_list$avg_log2FC), decreasing = TRUE)
names(gene_ranking) <- gene_list$ENTREZID
gene_ranking <- gene_ranking[!is.na(gene_ranking) & gene_ranking != 0]

# Run KEGG GSEA
gsea_kegg <- gseKEGG(geneList = gene_ranking, organism = "hsa", pvalueCutoff = 0.05, eps = 0)

# Run GO BP GSEA
gsea_go <- gseGO(geneList = gene_ranking, 
                 OrgDb = org.Hs.eg.db, 
                 ont = "BP", 
                 pvalueCutoff = 0.05, 
                 eps = 0,
                 nPermSimple = 10000,
                 minGSSize = 10,
                 maxGSSize = 500)

# Plot GSEA results
dotplot(gsea_go, showCategory = 10) +
  ggtitle("GO Biological Process GSEA (COPD-quit vs Control)")

dotplot(gsea_kegg, showCategory = 10) +
  ggtitle("KEGG Pathway GSEA (COPD-quit vs Control)")

# Custom dot plot
gsea_results <- as.data.frame(gsea_kegg@result)
top_pathways <- gsea_results %>%
  filter(p.adjust < 0.05) %>%
  arrange(NES) %>%
  top_n(15, wt = abs(NES)) %>%
  mutate(Count = setSize)

ggplot(top_pathways, aes(x = NES, y = reorder(Description, NES))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal() +
  xlab("Normalized Enrichment Score") +
  ylab("Pathway") +
  ggtitle("GSEA Enrichment Analysis (COPD-quit vs Control)") +
  theme(legend.position = "right")

ggsave("Q_C_GSEA_Enrichment_Analysis.png", width = 8, height = 6, dpi = 300)
##--------------------------------------------------------------------------------------

# GSEA Analysis: COPD-active vs Control (A vs C)

# Compute DEGs
de_genes <- FindMarkers(macrophage_obj, ident.1 = "COPD-active", ident.2 = "Control", logfc.threshold = 0, min.pct = 0.1)

# Convert gene symbols to Entrez IDs
gene_entrez <- bitr(rownames(de_genes), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_entrez <- gene_entrez[!duplicated(gene_entrez$SYMBOL), ]

# Merge with DEG data
de_genes$gene <- rownames(de_genes)
gene_list <- merge(gene_entrez, de_genes, by.x = "SYMBOL", by.y = "gene")

# Prepare gene ranking for GSEA
gene_ranking <- sort(jitter(gene_list$avg_log2FC), decreasing = TRUE)
names(gene_ranking) <- gene_list$ENTREZID
gene_ranking <- gene_ranking[!is.na(gene_ranking) & gene_ranking != 0]

# Run KEGG GSEA
gsea_kegg <- gseKEGG(geneList = gene_ranking, organism = "hsa", pvalueCutoff = 0.05, eps = 0)

# Run GO BP GSEA
gsea_go <- gseGO(geneList = gene_ranking, 
                 OrgDb = org.Hs.eg.db, 
                 ont = "BP", 
                 pvalueCutoff = 0.05, 
                 eps = 0,
                 nPermSimple = 10000,
                 minGSSize = 10,
                 maxGSSize = 500)

# Plot GSEA results
dotplot(gsea_go, showCategory = 10) +
  ggtitle("GO Biological Process GSEA (COPD-active vs Control)")

dotplot(gsea_kegg, showCategory = 10) +
  ggtitle("KEGG Pathway GSEA (COPD-active vs Control)")

# Custom dot plot
gsea_results <- as.data.frame(gsea_kegg@result)
top_pathways <- gsea_results %>%
  filter(p.adjust < 0.05) %>%
  arrange(NES) %>%
  top_n(15, wt = abs(NES)) %>%
  mutate(Count = setSize)

ggplot(top_pathways, aes(x = NES, y = reorder(Description, NES))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_minimal() +
  xlab("Normalized Enrichment Score") +
  ylab("Pathway") +
  ggtitle("GSEA Enrichment Analysis (COPD-active vs Control)") +
  theme(legend.position = "right")

ggsave("A_C_GSEA_Enrichment_Analysis.png", width = 8, height = 6, dpi = 300)

library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(msigdbr)
library(stringr)

# Optional; add this later if msigdbdf is actually needed
if (!requireNamespace("msigdbdf", quietly = TRUE)) {
  install.packages("msigdbdf", repos = "https://igordot.r-universe.dev")
}

output_dir <- "GSEA_Dotplots"
if (!dir.exists(output_dir)) dir.create(output_dir)

gsea_msigdb_analysis <- function(
    group1, group2,
    group_col = "group",
    categories = c("H", "C2", "C7"),
    out_prefix = "GSEA"
){
  Idents(macrophage_obj) <- macrophage_obj[[group_col]][,1]
  
  ## ① DEGs
  degs <- FindMarkers(
    macrophage_obj,
    ident.1 = group1, ident.2 = group2,
    logfc.threshold = 0, min.pct = 0.1
  )
  degs$gene <- rownames(degs)
  
  # Symbol → Entrez
  gene_entrez <- suppressWarnings(
    bitr(degs$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  )
  gene_entrez <- gene_entrez[!duplicated(gene_entrez$SYMBOL), ]
  gene_list <- merge(gene_entrez, degs, by.x = "SYMBOL", by.y = "gene")
  
  if (nrow(gene_list) == 0) {
    warning(paste("❗ No DE genes mapped to ENTREZID for:", group1, "vs", group2))
    return(NULL)
  }
  
  gene_ranking <- sort(jitter(gene_list$avg_log2FC), decreasing = TRUE)
  names(gene_ranking) <- gene_list$ENTREZID
  gene_ranking <- gene_ranking[!is.na(gene_ranking) & gene_ranking != 0]
  
  ## ② Loop over MSigDB categories
  for (cat in categories) {
    message("🔬 Running MSigDB GSEA for category: ", cat)
    
    msigdb <- tryCatch({
      msigdbr(species = "Homo sapiens", category = cat)
    }, error = function(e) {
      warning(paste("❗ Failed to retrieve MSigDB category:", cat, "-", e$message))
      return(NULL)
    })
    if (is.null(msigdb)) next
    
    msigdb_t2g <- msigdb[, c("gs_name", "entrez_gene")]
    
    gsea_res <- GSEA(
      geneList     = gene_ranking,
      TERM2GENE    = msigdb_t2g,
      pvalueCutoff = 0.05,
      verbose      = FALSE,
      eps          = 0
    )
    
    if (nrow(gsea_res@result) == 0) {
      warning(paste("⚠ No significant enrichment for:", cat, group1, "vs", group2))
      next
    }
    
    ## ③ Convert to df & clean pathway names
    gdf <- as.data.frame(gsea_res@result)
    
    gdf <- gdf %>%
      filter(!is.na(NES), NES > 0) %>%   # Plot NES > 0 only
      mutate(
        Description = toupper(Description),  # Convert to uppercase
        Description = gsub("^(REACTOME_|KEGG_|WP_|GO_BP_|GO_CC_|GO_MF_)", "", Description),
        Description = gsub("_", " ", Description),
        Description = stringr::str_wrap(Description, width = 50)
      ) %>%
      group_by(Description) %>%                 # Keep the pathway entry with the highest NES for duplicated pathway names
      slice_max(NES, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      arrange(desc(NES)) %>%
      slice_head(n = 15)                        # Keep the top 15 entries
    
    if (nrow(gdf) == 0) {
      message("⚠ No positive NES pathways remain after cleaning for category ", cat)
      next
    }
    
    gdf <- gdf %>%
      arrange(NES) %>%
      mutate(Description = factor(Description, levels = Description))
    
    ## 4. File name (defined before saving)
    out_name <- paste0(out_prefix, "_", group1, "_vs_", group2, "_", cat)
    
    # Save result table
    write.csv(
      gdf,
      file.path(output_dir, paste0(out_name, "_results.csv")),
      row.names = FALSE
    )
    
    ## 5. Large NES bubble plot
    p <- ggplot(gdf, aes(x = NES, y = Description)) +
      geom_point(aes(size = setSize, color = p.adjust), alpha = 0.92) +
      scale_size(range = c(4, 12), name = "Count") +
      scale_color_gradient(low = "red", high = "blue", name = "p.adjust") +
      labs(
        title    = paste0("GSEA Enrichment Analysis (", group1, " vs ", group2, ")"),
        subtitle = paste0("MSigDB category: ", cat),
        x        = "Normalized Enrichment Score",
        y        = NULL
      ) +
      theme_minimal(base_size = 16) +
      theme(
        plot.title    = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 13),
        axis.text.y   = element_text(size = 10, lineheight = 1.2),
        axis.text.x   = element_text(size = 11),
        legend.title  = element_text(size = 12),
        legend.text   = element_text(size = 10),
        panel.grid.minor = element_blank(),
        plot.margin   = margin(10, 30, 10, 15),
        legend.position = "right"
      )
    
    ggsave(
      file.path(output_dir, paste0(out_name, "_bubble.png")),
      p, width = 11, height = 8, dpi = 300
    )
  }
}

# 🚀 Run GSEA for 3 comparisons
gsea_msigdb_analysis("COPD-active", "COPD-quit", out_prefix = "A_Q")
gsea_msigdb_analysis("COPD-quit", "Control",        out_prefix = "Q_C")
gsea_msigdb_analysis("COPD-active", "Control",      out_prefix = "A_C")

## Combine results from the three comparisons into a facetable long-format data frame
gsea_msigdb_analysis <- function(
    group1, group2,
    group_col = "group",
    categories = c("H"),   # Main figures usually use Hallmark only
    out_prefix = "GSEA"
){
  Idents(macrophage_obj) <- macrophage_obj[[group_col]][,1]
  
  ## 1. DEGs (same as the current workflow)
  degs <- FindMarkers(
    macrophage_obj,
    ident.1 = group1, ident.2 = group2,
    logfc.threshold = 0, min.pct = 0.1
  )
  degs$gene <- rownames(degs)
  
  gene_entrez <- suppressWarnings(
    bitr(degs$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  )
  gene_entrez <- gene_entrez[!duplicated(gene_entrez$SYMBOL), ]
  gene_list <- merge(gene_entrez, degs, by.x = "SYMBOL", by.y = "gene")
  
  gene_ranking <- sort(jitter(gene_list$avg_log2FC), decreasing = TRUE)
  names(gene_ranking) <- gene_list$ENTREZID
  gene_ranking <- gene_ranking[!is.na(gene_ranking) & gene_ranking != 0]
  
  ## New: collect results from all categories
  all_res <- list()
  
  for (cat in categories) {
    
    msigdb <- msigdbr(species = "Homo sapiens", category = cat)
    msigdb_t2g <- msigdb[, c("gs_name", "entrez_gene")]
    
    gsea_res <- GSEA(
      geneList     = gene_ranking,
      TERM2GENE    = msigdb_t2g,
      pvalueCutoff = 0.05,
      eps          = 0,
      verbose      = FALSE
    )
    
    if (nrow(gsea_res@result) == 0) next
    
    gdf <- as.data.frame(gsea_res@result) %>%
      filter(!is.na(NES)) %>%
      mutate(
        Description = toupper(Description),
        Description = gsub("^HALLMARK_", "", Description),
        Description = gsub("_", " ", Description),
        comparison  = paste(group1, "vs", group2)
      )
    
    all_res[[cat]] <- gdf
  }
  
  ## Key step: return a long-format dataframe
  bind_rows(all_res)
}

df_A_C <- gsea_msigdb_analysis("COPD-active", "Control")
df_Q_C <- gsea_msigdb_analysis("COPD-quit",   "Control")
df_A_Q <- gsea_msigdb_analysis("COPD-active", "COPD-quit")
gsea_long <- bind_rows(df_A_C, df_Q_C, df_A_Q)


## =========================================================
## ADD-ON: Supplementary running enrichment curves
##         + Leading-edge genes heatmap
## Add this AFTER:
## gsea_long <- bind_rows(df_A_C, df_Q_C, df_A_Q)
## =========================================================

library(tidyr)
library(pheatmap)

supp_dir <- file.path(output_dir, "Supplementary_GSEA")
curve_dir <- file.path(supp_dir, "Running_enrichment_curves")
heatmap_dir <- file.path(supp_dir, "Leading_edge_heatmap")

dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(curve_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(heatmap_dir, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------
## 1) Pathways to show in supplementary running curves
##    Consistent with the main figure
## ---------------------------------------------------------

target_hallmark_terms <- c(
  "HALLMARK_BILE_ACID_METABOLISM",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_APICAL_SURFACE",
  "HALLMARK_COAGULATION",
  "HALLMARK_XENOBIOTIC_METABOLISM",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_PEROXISOME",
  "HALLMARK_ADIPOGENESIS",
  "HALLMARK_MITOTIC_SPINDLE"
)

nice_term <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  x
}

safe_term <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x
}

## ---------------------------------------------------------
## 2) Re-run Hallmark GSEA using your same logic,
##    but store gsea object + gene ranking for curves
## ---------------------------------------------------------

run_gsea_for_supp <- function(group1, group2, group_col = "group") {
  
  Idents(macrophage_obj) <- macrophage_obj[[group_col]][, 1]
  
  degs <- FindMarkers(
    macrophage_obj,
    ident.1 = group1,
    ident.2 = group2,
    logfc.threshold = 0,
    min.pct = 0.1
  )
  degs$gene <- rownames(degs)
  
  gene_entrez <- suppressWarnings(
    bitr(
      degs$gene,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )
  
  gene_entrez <- gene_entrez[!duplicated(gene_entrez$SYMBOL), ]
  
  gene_list <- merge(
    gene_entrez,
    degs,
    by.x = "SYMBOL",
    by.y = "gene"
  )
  
  gene_ranking <- sort(jitter(gene_list$avg_log2FC), decreasing = TRUE)
  names(gene_ranking) <- gene_list$ENTREZID
  gene_ranking <- gene_ranking[!is.na(gene_ranking) & gene_ranking != 0]
  gene_ranking <- gene_ranking[!duplicated(names(gene_ranking))]
  
  msigdb <- msigdbr(
    species = "Homo sapiens",
    category = "H"
  )
  
  msigdb_t2g <- msigdb[, c("gs_name", "entrez_gene")]
  
  gsea_res <- GSEA(
    geneList = gene_ranking,
    TERM2GENE = msigdb_t2g,
    pvalueCutoff = 1,
    eps = 0,
    verbose = FALSE
  )
  
  list(
    gsea = gsea_res,
    ranking = gene_ranking,
    msigdb_t2g = msigdb_t2g
  )
}

supp_A_C <- run_gsea_for_supp("COPD-active", "Control")
supp_A_Q <- run_gsea_for_supp("COPD-active", "COPD-quit")
supp_Q_C <- run_gsea_for_supp("COPD-quit", "Control")

supp_list <- list(
  A_C = supp_A_C,
  A_Q = supp_A_Q,
  Q_C = supp_Q_C
)

comparison_label <- c(
  A_C = "COPD-active vs Control",
  A_Q = "COPD-active vs COPD-quit",
  Q_C = "COPD-quit vs Control"
)

## ---------------------------------------------------------
## 3) Manual running enrichment curve
## ---------------------------------------------------------

plot_manual_gsea_curve <- function(gene_ranking, gene_set, title, nes, fdr, file_png, file_pdf) {
  
  genes <- names(gene_ranking)
  hits <- genes %in% gene_set
  
  N <- length(gene_ranking)
  Nh <- sum(hits)
  Nm <- N - Nh
  
  if (Nh == 0 || Nm == 0) return(NULL)
  
  stat_abs <- abs(gene_ranking)
  
  hit_score <- ifelse(hits, stat_abs / sum(stat_abs[hits]), 0)
  miss_score <- ifelse(!hits, 1 / Nm, 0)
  running_ES <- cumsum(hit_score - miss_score)
  
  plot_df <- data.frame(
    rank = seq_along(running_ES),
    ES = running_ES,
    hit = hits
  )
  
  hit_df <- plot_df %>% filter(hit)
  
  p <- ggplot(plot_df, aes(x = rank, y = ES)) +
    geom_hline(yintercept = 0, linewidth = 0.45, color = "grey45") +
    geom_vline(
      data = hit_df,
      aes(xintercept = rank),
      linewidth = 0.22,
      color = "firebrick",
      alpha = 0.65
    ) +
    geom_line(linewidth = 1.2, color = "black") +
    theme_bw(base_size = 16) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      plot.subtitle = element_text(size = 15, hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "grey25", fill = NA, linewidth = 0.8)
    ) +
    labs(
      title = title,
      subtitle = paste0("NES = ", round(nes, 2), ", FDR = ", signif(fdr, 3)),
      x = "Rank in ordered gene list",
      y = "Running enrichment score"
    )
  
  ggsave(file_png, p, width = 8, height = 6, dpi = 600, bg = "white")
  ggsave(file_pdf, p, width = 8, height = 6, bg = "white")
  
  return(p)
}

## ---------------------------------------------------------
## Make sure output folders exist
## ---------------------------------------------------------

output_dir <- "GSEA_Dotplots"

supp_dir <- file.path(output_dir, "Supplementary_GSEA")
curve_dir <- file.path(supp_dir, "Running_enrichment_curves")
heatmap_dir <- file.path(supp_dir, "Leading_edge_heatmap")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(curve_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(heatmap_dir, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------
## 4) Output running enrichment curves
## ---------------------------------------------------------

for (cmp in names(supp_list)) {
  
  gsea_obj <- supp_list[[cmp]]$gsea
  gene_ranking <- supp_list[[cmp]]$ranking
  msigdb_t2g <- supp_list[[cmp]]$msigdb_t2g
  
  res_df <- as.data.frame(gsea_obj@result)
  
  selected_df <- res_df %>%
    filter(ID %in% target_hallmark_terms) %>%
    filter(!is.na(NES), !is.na(p.adjust)) %>%
    filter(p.adjust < 0.05)
  
  write.csv(
    selected_df,
    file.path(curve_dir, paste0(cmp, "_selected_running_curve_pathways.csv")),
    row.names = FALSE
  )
  
  for (i in seq_len(nrow(selected_df))) {
    
    term_id <- selected_df$ID[i]
    
    term_genes <- msigdb_t2g %>%
      filter(gs_name == term_id) %>%
      pull(entrez_gene) %>%
      as.character()
    
    file_base <- paste0(cmp, "_H_", safe_term(term_id))
    
    plot_manual_gsea_curve(
      gene_ranking = gene_ranking,
      gene_set = term_genes,
      title = paste0(nice_term(term_id), "\n", comparison_label[cmp]),
      nes = selected_df$NES[i],
      fdr = selected_df$p.adjust[i],
      file_png = file.path(curve_dir, paste0(file_base, ".png")),
      file_pdf = file.path(curve_dir, paste0(file_base, ".pdf"))
    )
  }
}

## ---------------------------------------------------------
## 5) Leading-edge genes table
## ---------------------------------------------------------

leading_edge_all <- list()

for (cmp in names(supp_list)) {
  
  res_df <- as.data.frame(supp_list[[cmp]]$gsea@result)
  
  tmp <- res_df %>%
    dplyr::filter(ID %in% target_hallmark_terms) %>%
    dplyr::filter(!is.na(core_enrichment)) %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::select(
      ID,
      Description,
      NES,
      p.adjust,
      setSize,
      core_enrichment
    ) %>%
    tidyr::separate_rows(core_enrichment, sep = "/") %>%
    dplyr::rename(ENTREZID = core_enrichment) %>%
    dplyr::mutate(
      comparison_code = cmp,
      comparison = comparison_label[cmp]
    )
  
  leading_edge_all[[cmp]] <- tmp
}

leading_edge_df <- bind_rows(leading_edge_all)

entrez_symbol <- suppressWarnings(
  bitr(
    unique(leading_edge_df$ENTREZID),
    fromType = "ENTREZID",
    toType = "SYMBOL",
    OrgDb = org.Hs.eg.db
  )
)

leading_edge_df <- leading_edge_df %>%
  left_join(entrez_symbol, by = "ENTREZID") %>%
  relocate(comparison_code, comparison, ID, Description, SYMBOL, ENTREZID)

write.csv(
  leading_edge_df,
  file.path(supp_dir, "Supplementary_Table_leading_edge_genes.csv"),
  row.names = FALSE
)

## ---------------------------------------------------------
## 6) Leading-edge heatmap
##    Default pathway: INFLAMMATORY RESPONSE
## ---------------------------------------------------------

heatmap_term <- "HALLMARK_INFLAMMATORY_RESPONSE"

heatmap_genes <- leading_edge_df %>%
  filter(ID == heatmap_term) %>%
  filter(!is.na(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()

heatmap_genes <- heatmap_genes[heatmap_genes %in% rownames(macrophage_obj)]

cat("Leading-edge genes used in heatmap:", length(heatmap_genes), "\n")

if (length(heatmap_genes) >= 5) {
  
  avg_expr <- AverageExpression(
    macrophage_obj,
    features = heatmap_genes,
    group.by = "group",
    assays = DefaultAssay(macrophage_obj),
    slot = "data"
  )[[DefaultAssay(macrophage_obj)]]
  
  group_order <- c("Control", "COPD-quit", "COPD-active")
  avg_expr <- avg_expr[, intersect(group_order, colnames(avg_expr)), drop = FALSE]
  
  mat_z <- t(scale(t(avg_expr)))
  mat_z[is.na(mat_z)] <- 0
  
  png(
    filename = file.path(
      heatmap_dir,
      "Leading_edge_heatmap_HALLMARK_INFLAMMATORY_RESPONSE.png"
    ),
    width = 6,
    height = 8,
    units = "in",
    res = 600
  )
  
  pheatmap(
    mat_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize = 10,
    fontsize_row = 7,
    fontsize_col = 12,
    main = "Inflammatory response leading-edge genes",
    border_color = NA
  )
  
  dev.off()
  
  pdf(
    file = file.path(
      heatmap_dir,
      "Leading_edge_heatmap_HALLMARK_INFLAMMATORY_RESPONSE.pdf"
    ),
    width = 6,
    height = 8
  )
  
  pheatmap(
    mat_z,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize = 10,
    fontsize_row = 7,
    fontsize_col = 12,
    main = "Inflammatory response leading-edge genes",
    border_color = NA
  )
  
  dev.off()
}

cat("\n[DONE] Supplementary GSEA outputs generated:\n")
cat(normalizePath(supp_dir), "\n")

## =========================================================
## ADD-ON: Combine selected 5 running curves into 3 grouped panels
## One figure per comparison
## =========================================================

library(patchwork)
library(png)
library(grid)

combined_dir <- file.path(supp_dir, "Combined_5_running_curves")
dir.create(combined_dir, showWarnings = FALSE, recursive = TRUE)

selected_5_terms <- c(
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_XENOBIOTIC_METABOLISM"
)

make_img_plot <- function(img_file) {
  img <- png::readPNG(img_file)
  grid::rasterGrob(img, interpolate = TRUE) |>
    ggplotify::as.ggplot()
}

if (!requireNamespace("ggplotify", quietly = TRUE)) {
  install.packages("ggplotify")
}
library(ggplotify)

for (cmp in names(supp_list)) {
  
  plot_list <- list()
  
  for (term_id in selected_5_terms) {
    
    file_base <- paste0(cmp, "_H_", safe_term(term_id), ".png")
    img_file <- file.path(curve_dir, file_base)
    
    if (file.exists(img_file)) {
      plot_list[[term_id]] <- make_img_plot(img_file)
    } else {
      message("Missing file: ", img_file)
    }
  }
  
  if (length(plot_list) > 0) {
    
    combined_plot <- wrap_plots(plot_list, ncol = 2) +
      plot_annotation(
        title = paste0("Selected Hallmark GSEA running enrichment curves: ", comparison_label[cmp]),
        theme = theme(
          plot.title = element_text(
            hjust = 0.5,
            face = "bold",
            size = 18
          )
        )
      )
    
    ggsave(
      file.path(combined_dir, paste0(cmp, "_selected_5_running_curves.png")),
      combined_plot,
      width = 12,
      height = 14,
      dpi = 600,
      bg = "white"
    )
    
    ggsave(
      file.path(combined_dir, paste0(cmp, "_selected_5_running_curves.pdf")),
      combined_plot,
      width = 12,
      height = 14,
      bg = "white"
    )
  }
}


## =========================================================
## Three major pathway integrated leading-edge heatmap
## Lipid / Inflammatory / Oxidative programs
## =========================================================

library(dplyr)
library(pheatmap)
library(RColorBrewer)

integrated_dir <- file.path(supp_dir, "Integrated_three_pathway_heatmap")
dir.create(integrated_dir, showWarnings = FALSE, recursive = TRUE)

pathway_modules <- data.frame(
  ID = c(
    "HALLMARK_FATTY_ACID_METABOLISM",
    "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_PEROXISOME"
  ),
  Module = c(
    "Lipid metabolism",
    "Lipid metabolism",
    "Inflammatory program",
    "Oxidative / mitochondrial stress",
    "Oxidative / mitochondrial stress"
  )
)

## Extract leading-edge genes from the three major pathways
integrated_genes_df <- leading_edge_df %>%
  dplyr::filter(ID %in% pathway_modules$ID) %>%
  dplyr::filter(!is.na(SYMBOL)) %>%
  dplyr::left_join(pathway_modules, by = "ID") %>%
  dplyr::distinct(SYMBOL, Module, ID)

## Avoid duplicate genes when the same gene appears in multiple modules
integrated_genes_df <- integrated_genes_df %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

heatmap_genes <- integrated_genes_df$SYMBOL
heatmap_genes <- heatmap_genes[heatmap_genes %in% rownames(macrophage_obj)]

integrated_genes_df <- integrated_genes_df %>%
  dplyr::filter(SYMBOL %in% heatmap_genes)

cat("Genes used in integrated heatmap:", length(heatmap_genes), "\n")

## Average expression
avg_expr <- AverageExpression(
  macrophage_obj,
  features = heatmap_genes,
  group.by = "group",
  assays = DefaultAssay(macrophage_obj),
  slot = "data"
)[[DefaultAssay(macrophage_obj)]]

group_order <- c("Control", "COPD-quit", "COPD-active")
avg_expr <- avg_expr[, intersect(group_order, colnames(avg_expr)), drop = FALSE]

## Z-score
mat_z <- t(scale(t(avg_expr)))
mat_z[is.na(mat_z)] <- 0

## row annotation
annotation_row <- integrated_genes_df %>%
  dplyr::select(SYMBOL, Module) %>%
  as.data.frame()

rownames(annotation_row) <- annotation_row$SYMBOL
annotation_row$SYMBOL <- NULL

annotation_row <- annotation_row[rownames(mat_z), , drop = FALSE]

ann_colors <- list(
  Module = c(
    "Lipid metabolism" = "#E69F00",
    "Inflammatory program" = "#D55E00",
    "Oxidative / mitochondrial stress" = "#0072B2"
  )
)

## Save gene list
write.csv(
  integrated_genes_df,
  file.path(integrated_dir, "Integrated_three_pathway_leading_edge_genes.csv"),
  row.names = FALSE
)

## Plot PNG
png(
  filename = file.path(
    integrated_dir,
    "Integrated_three_pathway_leading_edge_heatmap.png"
  ),
  width = 7,
  height = 10,
  units = "in",
  res = 600
)

pheatmap(
  mat_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  fontsize = 10,
  fontsize_row = 7,
  fontsize_col = 12,
  main = "Leading-edge genes of key Hallmark pathways",
  border_color = NA,
  angle_col = 45,
  color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100)
)

dev.off()

## Plot PDF
pdf(
  file = file.path(
    integrated_dir,
    "Integrated_three_pathway_leading_edge_heatmap.pdf"
  ),
  width = 7,
  height = 10
)

pheatmap(
  mat_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  fontsize = 10,
  fontsize_row = 7,
  fontsize_col = 12,
  main = "Leading-edge genes of key Hallmark pathways",
  border_color = NA,
  angle_col = 45,
  color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100)
)

dev.off()

cat("[DONE] Integrated three-pathway heatmap saved to:\n")
cat(normalizePath(integrated_dir), "\n")

## =========================================================
## Split integrated leading-edge heatmap into 3 separate figures
## Remove Module annotation
## 1) Lipid metabolism
## 2) Inflammatory program
## 3) Oxidative / mitochondrial stress
## =========================================================

split_dir <- file.path(supp_dir, "Split_three_pathway_heatmaps_no_module")
dir.create(split_dir, showWarnings = FALSE, recursive = TRUE)

plot_one_module_heatmap <- function(module_name, file_prefix, top_n = 25) {
  
  genes_use <- integrated_genes_df %>%
    dplyr::filter(Module == module_name) %>%
    dplyr::pull(SYMBOL) %>%
    unique()
  
  genes_use <- genes_use[genes_use %in% rownames(mat_z)]
  
  ## optional: limit number of genes to avoid overcrowding
  genes_use <- genes_use[1:min(top_n, length(genes_use))]
  
  cat(module_name, "genes used:", length(genes_use), "\n")
  
  if (length(genes_use) < 3) {
    warning(paste("Too few genes for", module_name))
    return(NULL)
  }
  
  mat_sub <- mat_z[genes_use, , drop = FALSE]
  
  fig_height <- max(4.5, length(genes_use) * 0.22)
  
  ## ---------------- PNG ----------------
  png(
    filename = file.path(split_dir, paste0(file_prefix, ".png")),
    width = 5.2,
    height = fig_height,
    units = "in",
    res = 600
  )
  
  pheatmap(
    mat_sub,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 11,
    main = module_name,
    border_color = NA,
    angle_col = 45,
    color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100)
  )
  
  dev.off()
  
  ## ---------------- PDF ----------------
  pdf(
    file = file.path(split_dir, paste0(file_prefix, ".pdf")),
    width = 5.2,
    height = fig_height
  )
  
  pheatmap(
    mat_sub,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 11,
    main = module_name,
    border_color = NA,
    angle_col = 45,
    color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100)
  )
  
  dev.off()
}

plot_one_module_heatmap(
  module_name = "Lipid metabolism",
  file_prefix = "A_Lipid_metabolism_leading_edge_heatmap",
  top_n = 25
)

plot_one_module_heatmap(
  module_name = "Inflammatory program",
  file_prefix = "B_Inflammatory_program_leading_edge_heatmap",
  top_n = 25
)

plot_one_module_heatmap(
  module_name = "Oxidative / mitochondrial stress",
  file_prefix = "C_Oxidative_mitochondrial_stress_leading_edge_heatmap",
  top_n = 25
)

cat("[DONE] Three separated pathway heatmaps without module annotation saved to:\n")
cat(normalizePath(split_dir), "\n")
