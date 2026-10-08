# cytokine_heatmap.R
library(Seurat)
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
library(grid)

setwd("D:/OneDrive/OneDrive - The Mount Sinai Hospital/RG_sc_files/20250707_JAK1_GoF_Reanalysis")
obj <- readRDS("D:/OneDrive/OneDrive - The Mount Sinai Hospital/Ji_Lab_2023/RG_MERFISH_files/20250710_WT_JAK1_only/20260622_sc_merfish_int_full_obj_annotated_full.rds")

DefaultAssay(obj) <- "RNA"
Idents(obj) <- "assay"
obj <- subset(obj, idents = "scRNA")
obj <- JoinLayers(obj)

obj$condition <- as.character(obj$condition)
obj$condition[obj$condition == "Normal"] <- "WT"
obj$condition <- factor(obj$condition, levels = c("WT", "JAK1GoF"))

condition_cols <- c(WT = "grey50", JAK1GoF = "#7A3DB8")

cytokine_genes <- c(
  "Il4","Il10","Il13","Il5","Ifng","Il17a","Tslp",
  "Il12b","Il12a","Il1a","Il23a","Lif","Il33","Il1b",
  "Il6","Osm","Tgfb1","Tnf"
)

pb <- AggregateExpression(
  obj,
  assays = "RNA",
  features = cytokine_genes,
  group.by = "condition",
  slot = "counts",
  return.seurat = FALSE,
  verbose = FALSE
)

pb_mat <- as.matrix(pb$RNA)
pb_mat <- pb_mat[, c("WT", "JAK1GoF"), drop = FALSE]

pb_log <- log1p(pb_mat)

# Calculate log2FC (JAK1GoF vs WT)
log2fc <- log2((pb_mat[, "JAK1GoF"] + 1) /
                 (pb_mat[, "WT"] + 1))

# Append log2FC to cytokine names
rownames(pb_log) <- paste0(
  rownames(pb_log),
  " (",
  sprintf("%.1f", log2fc[rownames(pb_log)]),
  ")"
)

row_cor <- cor(t(pb_log), method = "pearson")
row_cor[is.na(row_cor)] <- 0
row_dist <- as.dist(1 - row_cor)
# Order rows by highest JAK1GoF expression to lowest
row_order <- order(pb_log[, "JAK1GoF"], decreasing = TRUE)
pb_plot <- pb_log[row_order, , drop = FALSE]

# column annotation bar
ha_col <- HeatmapAnnotation(
  condition = colnames(pb_plot),
  col = list(condition = condition_cols),
  show_annotation_name = FALSE
)

# heatmap colors
col_fun <- colorRamp2(
  c(min(pb_plot), median(pb_plot), max(pb_plot)),
  c("navy", "white", "firebrick")
)

cell_size <- unit(6, "mm")

p1 <- Heatmap(
  pb_plot,
  name = "Log Pseudobulk Expression",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  top_annotation = ha_col,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 12),
  rect_gp = gpar(col = "black", lwd = 1),
  width = ncol(pb_plot) * cell_size,
  height = nrow(pb_plot) * cell_size
)

pdf("20260505_Cytokine_log_pseudobulk_counts.pdf", width = 6, height = 12)
draw(p1)
dev.off()