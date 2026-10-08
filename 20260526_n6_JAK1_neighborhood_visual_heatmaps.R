library(Seurat)
library(tidyverse)
library(Matrix)
library(reticulate)
library(rhdf5)
library(anndata)
library(qs)
library(readxl)
library(ggplot2)
library(dplyr)
library(scales)
library(pheatmap)
library(grid)

base_out <- "D:/OneDrive/OneDrive - The Mount Sinai Hospital/RG_sc_files/20250707_JAK1_GoF_Reanalysis/Neighborhoods"
dir.create(base_out, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(
  "D:/OneDrive/OneDrive - The Mount Sinai Hospital/Ji_Lab_2023/RG_MERFISH_files/20250710_WT_JAK1_only/20260527_sc_merfish_int_full_obj_annotated_full.rds"
)

Idents(obj) <- "assay"
obj <- subset(obj, idents = "MERFISH")

nhood <- readxl::read_excel(
  "D:/OneDrive/OneDrive - The Mount Sinai Hospital/Ji_Lab_2023/RG_MERFISH_files/20250710_WT_JAK1_only/20250710_MERFISH_WT_JAK1.cellcharter_cluster_assignments.xlsx",
  sheet = 1
)

nhood$cell_barcode <- as.character(nhood$cell_barcode)
rownames(nhood) <- nhood$cell_barcode

obj$cell_barcode <- Cells(obj)
idx <- match(obj$cell_barcode, nhood$cell_barcode)

obj$nhood6 <- as.character(nhood$cellcharter_cluster.k_6[idx])

# annotate nhood6
obj$nhood6_annotated <- dplyr::recode(
  as.character(obj$nhood6),
  "0" = "N0: Subcutis",
  "1" = "N1: DSJ",
  "2" = "N2: HF_DP",
  "3" = "N3: Upper Dermis",
  "4" = "N4: Chondrocyte",
  "5" = "N5: Epidermis",
  .default = NA_character_
)

nhood_levels <- c(
  "N5: Epidermis",
  "N3: Upper Dermis",
  "N1: DSJ",
  "N0: Subcutis",
  "N2: HF_DP",
  "N4: Chondrocyte"
)

obj$nhood6_annotated <- factor(obj$nhood6_annotated, levels = rev(nhood_levels))

nhood_cols <- c(
  "N0: Subcutis" = "#8ac926",
  "N1: DSJ" = "#ff595e",
  "N2: HF_DP" = "#6a4c93",
  "N3: Upper Dermis" = "#1982c4",
  "N4: Chondrocyte" = "#ffca3a",
  "N5: Epidermis" = "#ff924c"
)

# quick summaries
print(round(prop.table(table(obj$nhood6_annotated, obj$site), margin = 1), 2))
print(round(prop.table(table(obj$nhood6_annotated, obj$condition), margin = 2), 4))

# spatial reduction from center_x / center_y
coords <- obj@meta.data[, c("center_x", "center_y")]
coords$center_x <- as.numeric(coords$center_x)
coords$center_y <- as.numeric(coords$center_y)
coords <- as.matrix(coords)
rownames(coords) <- Cells(obj)
colnames(coords) <- c("spatial_1", "spatial_2")

obj[["spatial"]] <- CreateDimReducObject(
  embeddings = coords,
  key = "spatial_",
  assay = DefaultAssay(obj)
)

# main spatial plot
p0 <- DimPlot(
  obj,
  reduction = "spatial",
  group.by = "nhood6_annotated",
  cols = nhood_cols,
  split.by = "condition",
  raster = FALSE,
  pt.size = 0.15
)

ggsave(
  file.path(base_out, "20260108_JAK1_nhood6.pdf"),
  p0,
  height = 24,
  width = 48,
  units = "in"
)


Idents(obj) <- 'celltype_broad'
levels(obj) <- sort(unique(obj$celltype_subtype))
# main spatial plot
p0 <- DimPlot(
  obj,
  reduction = "spatial",
  group.by = "celltype_subtype",
  cols = detailed.cell_type.colors,
  split.by = "condition",
  raster = FALSE,
  pt.size = 0.15
)

ggsave(
  file.path(base_out, "20260108_JAK1_celltypes.pdf"),
  p0,
  height = 24,
  width = 48,
  units = "in"
)


# highlight each cluster one at a time
spatial_df <- as.data.frame(obj@meta.data)
spatial_df$cluster <- as.character(spatial_df$nhood6_annotated)
spatial_df$center_x <- as.numeric(spatial_df$center_x)
spatial_df$center_y <- as.numeric(spatial_df$center_y)

all_clus <- na.omit(unique(spatial_df$cluster))

for (clus in all_clus) {
  cols <- rep("lightgrey", length(all_clus))
  names(cols) <- all_clus
  cols[clus] <- "red"
  
  g <- ggplot(spatial_df, aes(x = center_x, y = center_y, color = cluster)) +
    facet_wrap(~ site, scales = "free", nrow = 2) +
    geom_point(size = 0.6) +
    scale_color_manual(values = cols) +
    ggtitle("nhood6") +
    theme(legend.position = "none") +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  ggsave(
    file.path(base_out, paste0("spatial_merged_nhood6_highlight_", clus, "_red.png")),
    plot = g,
    width = 8,
    height = 12
  )
}

print(table(obj$nhood6_annotated))

# Cytokine / receptor heatmaps
df <- FetchData(
  obj,
  vars = c(
    "Osm",
    "Il6",
    "Osmr",
    "Il4",
    "Il13",
    "Il4ra",
    "Il13ra1",
    "Il6st",
    "Lifr",
    "nhood6_annotated"
  )
)

avg_exp <- df %>%
  group_by(nhood6_annotated) %>%
  summarise(
    Osm = mean(Osm, na.rm = TRUE),
    Il6 = mean(Il6, na.rm = TRUE),
    Osmr = mean(Osmr, na.rm = TRUE),
    Il4 = mean(Il4, na.rm = TRUE),
    Il13 = mean(Il13, na.rm = TRUE),
    Il4ra = mean(Il4ra, na.rm = TRUE),
    Il13ra1 = mean(Il13ra1, na.rm = TRUE),
    Il6st = mean(Il6st, na.rm = TRUE),
    Lifr = mean(Lifr, na.rm = TRUE),
    .groups = "drop"
  )

cytokine_mat <- as.matrix(avg_exp[, c("Osm", "Il6", "Il4", "Il13")])
rownames(cytokine_mat) <- as.character(avg_exp$nhood6_annotated)

receptor_mat <- as.matrix(avg_exp[, c("Osmr", "Il4ra", "Il13ra1", "Il6st", "Lifr")])
rownames(receptor_mat) <- as.character(avg_exp$nhood6_annotated)

nhood_order <- c(
  "N5: Epidermis",
  "N3: Upper Dermis",
  "N1: DSJ",
  "N0: Subcutis",
  "N2: HF_DP",
  "N4: Chondrocyte"
)

cytokine_mat <- cytokine_mat[match(nhood_order, rownames(cytokine_mat)), , drop = FALSE]
receptor_mat <- receptor_mat[match(nhood_order, rownames(receptor_mat)), , drop = FALSE]

hm_cols <- colorRampPalette(c("blue", "white", "firebrick"))(100)
# row annotation for neighborhoods
ann_row <- data.frame(
  Neighborhood = factor(
    rownames(cytokine_mat),
    levels = c(
      "N5: Epidermis",
      "N3: Upper Dermis",
      "N1: DSJ",
      "N0: Subcutis",
      "N2: HF_DP",
      "N4: Chondrocyte"
    )
  )
)
rownames(ann_row) <- rownames(cytokine_mat)

ann_colors <- list(
  Neighborhood = c(
    "N0: Subcutis" = "#8ac926",
    "N1: DSJ" = "#ff595e",
    "N2: HF_DP" = "#6a4c93",
    "N3: Upper Dermis" = "#1982c4",
    "N4: Chondrocyte" = "#ffca3a",
    "N5: Epidermis" = "#ff924c"
  )
)

p1 <- pheatmap(
  cytokine_mat,
  scale = "column",
  color = hm_cols,
  border_color = "black",
  cellwidth = 25,
  cellheight = 25,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  treeheight_row = 0,
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  silent = TRUE
)

p2 <- pheatmap(
  receptor_mat,
  scale = "column",
  color = hm_cols,
  border_color = "black",
  cellwidth = 25,
  cellheight = 25,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  treeheight_row = 0,
  annotation_row = ann_row,
  annotation_colors = ann_colors,
  silent = TRUE
)

pdf(
  file.path(base_out, "20260108_cytokine_receptor_heatmaps_nhood6.pdf"),
  width = 12,
  height = 5
)

grid.arrange(
  p1$gtable,
  p2$gtable,
  ncol = 2
)

dev.off()


# tissue barplot colored by neighborhood
tab_site <- as.data.frame(table(obj$site, obj$nhood6_annotated))
colnames(tab_site) <- c("Tissue", "Neighborhood", "Freq")

tab_site$Neighborhood <- factor(
  tab_site$Neighborhood,
  levels = c(
    "N5: Epidermis",
    "N3: Upper Dermis",
    "N1: DSJ",
    "N0: Subcutis",
    "N2: HF_DP",
    "N4: Chondrocyte"
  )
)

p_site <- ggplot(tab_site, aes(x = Tissue, y = Freq, fill = Neighborhood)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = nhood_cols) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  ylab("Proportion") +
  xlab("Tissue")

ggsave("20260108_nhood6_site_barplot.pdf", p_site, width = 4, height = 3)

###### Immune Cell types in neighborhoods ############
library(dplyr)
library(pheatmap)

df <- FetchData(
  obj,
  vars = c("celltype_subtype", "nhood6_annotated", "condition")
) %>%
  filter(condition == "JAK1GoF") %>%
  filter(celltype_subtype %in% c("Basophil", "ILC2", "DC2")) %>%
  filter(!is.na(nhood6_annotated))

tab <- table(df$celltype_subtype, df$nhood6_annotated)

# fraction of each cell type across neighborhoods
prop_mat <- prop.table(tab, margin = 1)

# neighborhood order
nhood_order <- c(
  "N5: Epidermis",
  "N3: Upper Dermis",
  "N1: DSJ",
  "N0: Subcutis",
  "N2: HF_DP",
  "N4: Chondrocyte"
)

# reorder columns
prop_mat <- prop_mat[, intersect(nhood_order, colnames(prop_mat)), drop = FALSE]
tab <- tab[, colnames(prop_mat), drop = FALSE]

ann_col <- data.frame(Neighborhood = colnames(prop_mat))
rownames(ann_col) <- colnames(prop_mat)

ann_colors <- list(
  Neighborhood = c(
    "N0: Subcutis" = "#8ac926",
    "N1: DSJ" = "#ff595e",
    "N2: HF_DP" = "#6a4c93",
    "N3: Upper Dermis" = "#1982c4",
    "N4: Chondrocyte" = "#ffca3a",
    "N5: Epidermis" = "#ff924c"
  )
)

hm_cols <- colorRampPalette(c("white", "firebrick"))(100)

pdf(
  "20260108_JAK1GoF_Basophil_ILC2_DC2_neighborhood_heatmap.pdf",
  width = 6,
  height = 4
)

pheatmap(
  prop_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = hm_cols,
  border_color = "black",
  cellwidth = 35,
  cellheight = 35,
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  number_format = "%.0f",
  main = "JAK1 GoF localization of Basophil / ILC2 / DC2"
)

dev.off()


########Heatmap celltypes across all neighborhoods ##########
setwd(base_out)

detailed.cell_type.colors <- c("Adipocyte"="#FBDB79", 
                               "Basophil" ='#EE1289',
                               "Bulge" = "lightblue",
                               "Chondrocyte" ="#D9A5BA",
                               "DC1" ="chartreuse1",
                               "DC2" ="chartreuse4",
                               "FIB"="#D96F91",
                               "FIB_Subcutis" = '#FF8C00',
                               "FIB_Reticular"="red3",
                               "FIB_Subfascia" = "#0000CD",
                               "FIB_Papillary" ="#BA55D3",
                               "FIB_DP"="#FFD700",
                               "FIB_DS"="#9ACD35",
                               "HF_IRS" ='orange4',
                               "HF_ORS"='tan3',
                               "KC_Bas" ="#358A8E",
                               "KC_IFE"  ="#044D6E",
                               "KC_Spn" ="#BBE4D7" ,
                               "LC" = "#556B2F",
                               "LEC"= "#8B7355",
                               "TC_gdt"='firebrick4',
                               "TC_Treg"='lightgoldenrod',
                               "TC_Th17"='lightpink',
                               "TC_Cd8"='chocolate',
                               "NK"='steelblue',
                               "ILC2"='purple',
                               "TC_Cyc"='wheat2',
                               "MAC"="#00BFFF",
                               "Mast"="#8B4C39",
                               "Neutrophil"="#8B668B",
                               "Pericyte"="#708090",
                               "Schwann"="#EE8066",
                               "SG"="#B5B6D7",
                               "SMC"='#F82705',
                               "VEC" = "#D2B48C")

library(dplyr)
library(pheatmap)

# pull metadata
df <- FetchData(
  obj,
  vars = c("celltype_subtype", "nhood6_annotated", "condition")
) %>%
  filter(condition == "JAK1GoF") %>%
  filter(!is.na(nhood6_annotated)) %>%
  filter(!is.na(celltype_subtype))

# count cells
tab <- table(df$celltype_subtype, df$nhood6_annotated)

# row-normalize:
# fraction of each cell type across neighborhoods
prop_mat <- prop.table(tab, margin = 1)

# desired neighborhood order
nhood_order <- c(
  "N5: Epidermis",
  "N3: Upper Dermis",
  "N1: DSJ",
  "N0: Subcutis",
  "N2: HF_DP",
  "N4: Chondrocyte"
)

prop_mat <- prop_mat[, intersect(nhood_order, colnames(prop_mat)), drop = FALSE]
tab <- tab[, colnames(prop_mat), drop = FALSE]

# neighborhood annotation colors
ann_col <- data.frame(
  Neighborhood = factor(colnames(prop_mat), levels = nhood_order)
)
rownames(ann_col) <- colnames(prop_mat)

ann_colors <- list(
  Neighborhood = c(
    "N0: Subcutis" = "#8ac926",
    "N1: DSJ" = "#ff595e",
    "N2: HF_DP" = "#6a4c93",
    "N3: Upper Dermis" = "#1982c4",
    "N4: Chondrocyte" = "#ffca3a",
    "N5: Epidermis" = "#ff924c"
  )
)

# row annotation for celltype_subtype colors
row_ann <- data.frame(
  CellType = rownames(prop_mat)
)
rownames(row_ann) <- rownames(prop_mat)

# keep only colors that are actually used
celltype_cols_use <- detailed.cell_type.colors[names(detailed.cell_type.colors) %in% rownames(prop_mat)]

# optional: warn if any row types have no color assigned
missing_ct_cols <- setdiff(rownames(prop_mat), names(detailed.cell_type.colors))
if (length(missing_ct_cols) > 0) {
  message("No color assigned for: ", paste(missing_ct_cols, collapse = ", "))
}

ann_colors$row_ann <- NULL
ann_colors <- list(
  Neighborhood = ann_colors$Neighborhood,
  CellType = detailed.cell_type.colors
)

# heatmap palette
hm_cols <- colorRampPalette(c("white", "firebrick"))(100)

# save PDF
pdf(
  "20260108_JAK1GoF_all_celltypes_neighborhood_heatmap.pdf",
  width = 8,
  height = 12
)

pheatmap(
  prop_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = hm_cols,
  border_color = "black",
  annotation_col = ann_col,
  annotation_row = row_ann,
  annotation_colors = ann_colors,
  display_numbers = FALSE,
  cellwidth = 15,
  cellheight = 15,
  fontsize_row = 10,
  breaks = seq(0, 1, length.out = 101),
  legend_breaks = c(0, 0.5, 1),
  legend_labels = c("0", "0.5", "1"),
  main = "JAK1 GoF cell type localization across neighborhoods"
)

dev.off()

############


library(ggrepel)
library(patchwork)

# keep only cells with a neighborhood and the two conditions of interest
obj_de <- subset(
  obj,
  subset = !is.na(nhood6_annotated) & condition %in% c("WT", "JAK1GoF")
)

# helper to make one volcano plot
make_volcano <- function(
    seu,
    nhood_name,
    min_pct = 0.1,
    fc_cut = 0.25,
    padj_cut = 0.05
) {
  
  sub <- subset(seu, subset = nhood6_annotated == nhood_name)
  
  tab_cond <- table(sub$condition)
  if (length(tab_cond) < 2 || any(tab_cond == 0)) {
    message("Skipping ", nhood_name, " because one condition is missing.")
    return(NULL)
  }
  
  Idents(sub) <- "condition"
  
  deg <- FindMarkers(
    sub,
    ident.1 = "JAK1GoF",
    ident.2 = "WT",
    logfc.threshold = 0,
    min.pct = min_pct,
    test.use = "wilcox"
  )
  
  deg$gene <- rownames(deg)
  deg$padj <- ifelse(is.na(deg$p_val_adj), 1, deg$p_val_adj)
  deg$neglog10_padj <- -log10(deg$padj + 1e-300)
  
  deg$sig <- "NS"
  deg$sig[deg$avg_log2FC >= fc_cut & deg$padj < padj_cut] <- "Up in JAK1GoF"
  deg$sig[deg$avg_log2FC <= -fc_cut & deg$padj < padj_cut] <- "Up in WT"
  
  deg_label <- deg %>%
    filter(padj < padj_cut) %>%
    arrange(padj) %>%
    slice_head(n = 25)
  
  p <- ggplot(deg, aes(avg_log2FC, neglog10_padj)) +
    geom_point(aes(color = sig), size = 2, alpha = 0.8) +
    scale_color_manual(
      values = c(
        "NS" = "grey70",
        "Up in JAK1GoF" = "#7A3DB8",
        "Up in WT" = "#1F78B4"
      )
    ) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    geom_text_repel(
      data = deg_label,
      aes(label = gene),
      size = 6,
      max.overlaps = Inf,
      segment.color = NA
    ) +
    theme_classic(base_size = 16) +
    labs(
      title = nhood_name,
      x = "avg log2FC (JAK1GoF vs WT)",
      y = "-log10(adj p-value)"
    ) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_blank(),
      legend.text = element_text(size = 14)
    )
  
  ggsave(
    paste0(
      "volcano_",
      gsub("[^A-Za-z0-9_]+", "_", nhood_name),
      ".pdf"
    ),
    plot = p,
    width = 8,
    height = 6
  )
  
  return(deg)
}


deg_list <- lapply(
  levels(droplevels(obj_de$nhood6_annotated)),
  function(nh) {
    message("Running ", nh)
    
    res <- make_volcano(obj_de, nh)
    
    if (!is.null(res)) {
      res$gene <- rownames(res)
      res$neighborhood <- nh
    }
    
    res
  }
)

names(deg_list) <- levels(droplevels(obj_de$nhood6_annotated))

all_deg <- dplyr::bind_rows(deg_list)

sig_deg <- all_deg %>%
  filter(
    p_val_adj < 0.05,
    abs(avg_log2FC) > 0.25
  ) %>%
  arrange(neighborhood, desc(avg_log2FC))

sig_deg %>%
  group_by(neighborhood) %>%
  slice_max(avg_log2FC, n = 20)
