library(Seurat)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

obj <- readRDS("D:/OneDrive/OneDrive - The Mount Sinai Hospital/Ji_Lab_2023/RG_MERFISH_files/20250710_WT_JAK1_only/20260527_sc_merfish_int_full_obj_annotated_full.rds")

Idents(obj) <- "assay"
obj <- subset(obj, idents = "scRNA")
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
obj <- NormalizeData(obj)
obj <- ScaleData(obj)

obj$condition <- as.character(obj$condition)
obj$condition[obj$condition == "Normal"] <- "WT"
obj$condition <- factor(obj$condition, levels = c("WT", "JAK1GoF"))

obj@meta.data <- obj@meta.data %>%
  mutate(celltype_grouped = case_when(
    grepl("^FIB", celltype_subtype) ~ "FIB",
    grepl("^KC", celltype_subtype) ~ "KC",
    grepl("^HF", celltype_subtype) ~ "HF",
    grepl("^DC", celltype_subtype) ~ "DC",
    TRUE ~ celltype_subtype
  ))

receptor_genes <- c("Il4ra", "Il13ra1", "Il6ra", "Il6st", "Osmr","Il13ra2","Lifr")

avg <- AverageExpression(
  obj,
  assays = "RNA",
  features = receptor_genes,
  group.by = c("celltype_grouped", "condition"),
  slot = "data"
)

mat <- as.matrix(avg$RNA)

# order columns by cell type, then condition
celltypes <- unique(obj$celltype_grouped)
col_order <- unlist(lapply(celltypes, function(ct) {
  paste0(ct, c("_WT", "_JAK1GoF"))
}))
col_order <- intersect(col_order, colnames(mat))
mat <- mat[, col_order, drop = FALSE]

# annotation
annotation_col <- data.frame(
  celltype = sub("_(WT|JAK1GoF)$", "", colnames(mat)),
  condition = ifelse(grepl("WT$", colnames(mat)), "WT", "JAK1GoF")
)
rownames(annotation_col) <- colnames(mat)

# gaps between cell types
gaps_col <- which(annotation_col$celltype[-1] != annotation_col$celltype[-nrow(annotation_col)])

# normalize within each column
mat_norm <- scale(mat)
mat_norm[is.na(mat_norm)] <- 0

pheatmap(
  mat_norm,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  gaps_col = gaps_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = annotation_col,
  annotation_colors = list(
    celltype = setNames(
      colorRampPalette(brewer.pal(8, "Set3"))(length(unique(annotation_col$celltype))),
      unique(annotation_col$celltype)
    ),
    condition = c(WT = "grey70", JAK1GoF = "#7A3DB8")
  ),
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
)



library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

genes <- c("Il4ra", "Il13ra1", "Il6ra", "Il6st", "Osmr","Il13ra2","Lifr")

# keep only the cell types you want
obj_sub <- subset(obj, subset = celltype_grouped %in% c("KC", "FIB", "Pericyte","Chondrocyte"))

# make sure condition order is correct
obj_sub$condition <- as.character(obj_sub$condition)
obj_sub$condition[obj_sub$condition == "Normal"] <- "WT"
obj_sub$condition <- factor(obj_sub$condition, levels = c("WT", "JAK1GoF"))

obj_sub$celltype_grouped <- factor(obj_sub$celltype_grouped, levels = c("KC", "FIB", "Pericyte","Chondrocyte"))

out_pdf <- "20260505_Receptor_boxplots_KC_FIB_Pericyte.pdf"
pdf(out_pdf, width = 10, height = 4)

for (gene in genes) {
  df <- FetchData(obj_sub, vars = c(gene, "celltype_grouped", "condition"))
  colnames(df)[1] <- "expr"
  
  p <- ggplot(df, aes(x = condition, y = expr, fill = condition)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(aes(color = condition), width = 0.12, size = 1.1, alpha = 0.5) +
    facet_wrap(~celltype_grouped, nrow = 1, scales = "free_y") +
    stat_compare_means(
      comparisons = list(c("WT", "JAK1GoF")),
      method = "wilcox.test",
      label = "p.signif",
      label.y.npc = 0.96,
      size = 3
    ) +
    scale_fill_manual(values = c("WT" = "grey70", "JAK1GoF" = "#7A3DB8")) +
    scale_color_manual(values = c("WT" = "grey45", "JAK1GoF" = "#5A2C91")) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      plot.margin = margin(8, 8, 12, 8)
    ) +
    labs(title = gene, x = NULL, y = "Expression")
  
  print(p)
}

dev.off()

########## Separate Plots ##########
setwd("D:/OneDrive/OneDrive - The Mount Sinai Hospital/RG_sc_files/20250707_JAK1_GoF_Reanalysis/Receptor_Analysis")

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

genes <- c("Il4ra", "Il13ra1", "Il6ra", "Il6st", "Osmr", "Il13ra2", "Lifr","Flg","Il1rn")
celltypes_to_plot <- c("KC", "FIB", "Pericyte")

# make sure condition order is correct
obj$condition <- as.character(obj$condition)
obj$condition[obj$condition == "Normal"] <- "WT"
obj$condition <- factor(obj$condition, levels = c("WT", "JAK1GoF"))

for (ct in celltypes_to_plot) {
  obj_sub <- subset(obj, subset = celltype_grouped == ct)
  
  for (gene in genes) {
    if (!gene %in% rownames(obj_sub)) next
    
    df <- FetchData(obj_sub, vars = c(gene, "condition"))
    colnames(df)[1] <- "expr"
    
    ymax <- max(df$expr, na.rm = TRUE)
    yrng <- diff(range(df$expr, na.rm = TRUE))
    ypos <- ymax + 0.20 * yrng
    
    p <- ggplot(df, aes(x = condition, y = expr, fill = condition)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.9) +
  geom_jitter(aes(color = condition), width = 0.12, size = 1.1, alpha = 0.5) +
  stat_compare_means(
    comparisons = list(c("WT", "JAK1GoF")),
    method = "wilcox.test",
    label = "p.signif",
    label.y = ypos,
    size = 6,
    bracket.size = 0.7,
    tip.length = 0.05
  ) +
      stat_summary(
        fun = mean,
        geom = "crossbar",
        color = "black",
        width = 0.3
      ) +
  scale_fill_manual(values = c("WT" = "grey70", "JAK1GoF" = "#7A3DB8")) +
  scale_color_manual(values = c("WT" = "grey45", "JAK1GoF" = "#5A2C91")) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) + 
      coord_cartesian(clip = "off") + 
      theme_classic(base_size = 14) + 
      theme( legend.position = "none", 
             strip.background = element_blank(), 
             strip.text = element_text(face = "bold"), 
             axis.text.x = element_text(face = "bold"), 
             plot.margin = margin(8, 8, 12, 8) ) + 
      labs( title = paste(ct, gene), x = NULL, y = "Expression")
    out_pdf <- paste0("20260505_Receptor_vlnplot_", ct, "_", gene, ".pdf")
    ggsave(out_pdf, plot = p, width = 3, height = 4)
  }
}
