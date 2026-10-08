setwd("D:/OneDrive/OneDrive - The Mount Sinai Hospital/Ji_Lab_2023/RG_MERFISH_files/20250710_WT_JAK1_only")
library(Seurat)
library(ggpubr)

# obj <- readRDS("20241119_sc_merfish_int_full_obj.rds")
# 
# DefaultAssay(obj) <- "integrated"
# 
# obj <- ScaleData(obj, verbose = TRUE, vars.to.regress = c("nCount_RNA"))
# obj <- RunPCA(obj, npcs = 30, verbose = TRUE)
# 
# obj <- FindNeighbors(obj, dims = 1:30, verbose = T)
# obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)
# obj <- FindClusters(obj, resolution = 0.6, verbose = FALSE)
# obj <- FindClusters(obj, resolution = 0.8, verbose = FALSE)
# obj <- FindClusters(obj, resolution = 1, verbose = FALSE)
# obj <- FindClusters(obj, resolution = 1.1, verbose = FALSE)
# obj <- FindClusters(obj, resolution = 1.2, verbose = FALSE)
# obj <- FindClusters(obj, resolution = 1.4, verbose = FALSE)
# obj <- FindClusters(obj, resolution = 1.5, verbose = FALSE)
# obj <- FindClusters(obj, resolution = 1.6, verbose = FALSE)
# 
# 
 Idents(obj) <- 'condition'
 obj <- RenameIdents(obj, 
                     "Normal" = "WT",
                     "JAK1" = "JAK1GoF")
 obj$condition <- Idents(obj)
# 
# obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)
# DimPlot(obj, group.by = "assay")
# DimPlot(obj, split.by = "assay")
# 
# saveRDS(obj, "20241119_sc_merfish_int_full_obj.rds")



obj <- readRDS("20241119_sc_merfish_int_full_obj.rds")

Idents(obj) <- 'integrated_snn_res.1'
DefaultAssay(obj) <- 'RNA'
#markers <- FindAllMarkers(obj, only.pos = T)
#write.csv(markers, "20260622_Res1_Cluster_markers.csv")
obj <- RenameIdents(obj,
                    "0"="FIB",
                    "1"="FIB",
                    "2"="KC_IFE",
                    "3"="Chondrocyte",
                    "4"="KC_Spn",
                    "5"="KC_IFE",
                    "6"="SG",
                    "7"="SMC",
                    "8"="Myeloid",
                    "9"="HF_IRS",
                    "10"="VEC",
                    "11"="Adipocyte",
                    "12"="Bulge",
                    "13"="Myeloid",
                    "14"="FIB",
                    "15"="KC_Bas",
                    "16"="Schwann",
                    "17"="SG",
                    "18"="FIB",
                    "19"="HF_ORS",
                    "20"="Pericyte",
                    "21"="Lymphocyte",
                    "22"="Myeloid",
                    "23"="LEC",
                    "24"="Mast",
                    "25"="Myeloid",
                    "26"="Lymphocyte",
                    "27"="FIB",
                    "28"="Myeloid"
)

obj$celltype_broad <- as.character(Idents(obj))

celltype_colors <- list("Adipocyte"="#FBDB79",
                        "Bulge" = "lightblue",
                        "Chondrocyte" ="#D9A5BA",
                        "FIB"="#D96F91",
                        "HF_IRS" ='orange4',
                        "HF_ORS"='tan3',
                        "KC_Bas" ="#358A8E",
                        "KC_IFE"  ="#044D6E",
                        "KC_Spn" ="#BBE4D7" ,
                        "LEC"= "burlywood4",
                        "Lymphocyte"="plum",
                        "Mast"="salmon4",
                        "Myeloid"="green4",
                        "Pericyte"="slategrey",
                        "Schwann"="#EE8066",
                        "SG"="#B5B6D7",
                        "SMC"='#F82705',
                        "VEC" = "tan") 


detailed.cell_type.colors <- c("Adipocyte"="#FBDB79", 
                               "Basophil" ='#EE1289',
                               "Bulge" = "lightblue",
                               "Chondrocyte" ="#D9A5BA",
                               "DC_cDC1" ="#41AB5D",
                               "DC_cDC2" ="#A1D99B",
                               "DC_Mig"="#556B2F",
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
                               "LC" = "#00441B",
                               "LEC"= "#8B7355",
                               "TC_gdt"='firebrick4',
                               "TC_Treg"='lightgoldenrod',
                               "TC_Th17"='lightpink',
                               "NK_Cd8"='chocolate',
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


colors <- unname(unlist(celltype_colors))

Idents(obj) <- 'celltype_broad'
levels(obj) <- sort(unique(obj@meta.data$celltype_broad))
p1 <- DimPlot(obj, cols = colors, label = T, label.box = T, repel = T) & NoAxes() & NoLegend()
png("Res1_broad_labeled_Dimplot.png", res = 800, width = 4.75, height = 4.75, units = 'in')
print(p1)
dev.off()
p2 <- DimPlot(obj, cols = colors, split.by = 'assay') & NoAxes() & NoLegend() 
p2 <- LabelClusters(plot = p2, id = 'ident', fontface = 'bold', repel = T, size = 5)
png("Res1_labeled_Dimplot_splitby_assay_12_6.png", res = 800, width = 9.5, height = 4.75, units = 'in')
print(p2)
dev.off()

##Myeloid Recluster
Idents(obj) <- 'celltype_broad'
myl <- subset(obj, idents = c("Myeloid"))
DefaultAssay(myl) <- 'integrated'
myl <- ScaleData(myl, verbose = TRUE, vars.to.regress = c("nCount_RNA"))
myl <- RunPCA(myl, npcs = 30, verbose = TRUE)

myl <- FindNeighbors(myl, dims = 1:30, verbose = T)
myl <- FindClusters(myl, resolution = 0.2, verbose = FALSE)
myl <- FindClusters(myl, resolution = 0.3, verbose = FALSE)
myl <- FindClusters(myl, resolution = 0.4, verbose = FALSE)
myl <- FindClusters(myl, resolution = 0.6, verbose = FALSE)
myl <- FindClusters(myl, resolution = 0.8, verbose = FALSE)
myl <- FindClusters(myl, resolution = 1, verbose = FALSE)


myl <- RunUMAP(myl, reduction = "pca", dims = 1:30)
DimPlot(myl, group.by = "assay")
p2 <- DimPlot(myl, split.by = "assay", group.by = 'integrated_snn_res.0.3', label = T)
p2

Idents(myl)<- 'integrated_snn_res.0.3'
myl.markers <- FindAllMarkers(myl, only.pos = T)

DimPlot(myl, label = T, group.by = 'celltype_broad')
DimPlot(myl, label = T)


myl <- RenameIdents(myl, "0" = "MAC",
                    "1" = "MAC",
                    "2" = "MAC",
                    "3" = "MAC",
                    "4"= "LC",
                    "5"= "DC_cDC1",
                    "6" = "Basophil",
                    "7" = "MAC",
                    "8" = "Neutrophil",
                    "9" = "DC_Mig")

DimPlot(myl, label = T)

myl$celltype_subtype <- as.character(Idents(myl))
Idents(myl)<- 'celltype_subtype'
DefaultAssay(myl) <- 'RNA'
myl.markers <- FindAllMarkers(myl, only.pos = T)



myl_markers <- list("Basophil" = c("Mcpt8","Hdc"),
                    "DC_cDC1" = c("Clec9a","Xcr1"),
                    "DC_Mig" = c("Ccr7","Ccl22"),
                    "LC" = c("Cd207","Hpgds"),
                    "MAC" = c("Mrc1","Cd163"),
                    "Neutrophil" = c("Acod1","S100a8"))

myl_markers_rev <- (myl_markers)

Idents(myl) <- 'celltype_subtype'

DefaultAssay(myl) <- 'RNA'
levels(myl) <- rev(sort(as.character(unique(Idents(myl)))))

p1 <- DotPlot(myl, features = myl_markers_rev, scale.max = 50) + theme_bw() + theme(panel.spacing = unit(0.2, "lines")) + rotate_x_text(45) +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=0))+
  theme(strip.text = element_blank())+
  guides(color = guide_colorbar(title = "Average\nExpression", frame.colour = 'black'),
         size = guide_legend(title = "Percent\nExpressed")) +
  theme(legend.position = 'bottom',
        panel.grid.minor = element_line(color = "black", size = 0.5))

ggsave("myl_DotPlot_labeled.pdf",p1, width = 4, height = 4, units = 'in')







obj$celltype_subtype <- as.character(obj$celltype_broad)

obj@meta.data[rownames(myl@meta.data),]$celltype_subtype <- myl@meta.data$celltype_subtype

obj$myl_sub_res_0.3 <- "NA"
obj@meta.data[rownames(myl@meta.data),]$myl_sub_res_0.3 <- myl@meta.data$integrated_snn_res.0.3


DimPlot(obj, group.by = 'celltype_subtype', split.by = 'condition', label = T, repel = T)

########FIB Recluster##############
Idents(obj) <- 'celltype_broad'
fib <- subset(obj, idents = c("FIB"))
DefaultAssay(fib) <- 'integrated'
fib <- ScaleData(fib, verbose = TRUE, vars.to.regress = c("nCount_RNA"))
fib <- RunPCA(fib, npcs = 30, verbose = TRUE)

fib <- FindNeighbors(fib, dims = 1:30, verbose = T)
fib <- FindClusters(fib, resolution = 0.2, verbose = FALSE)
fib <- FindClusters(fib, resolution = 0.3, verbose = FALSE)
fib <- FindClusters(fib, resolution = 0.4, verbose = FALSE)
fib <- FindClusters(fib, resolution = 0.5, verbose = FALSE)
fib <- FindClusters(fib, resolution = 0.6, verbose = FALSE)
fib <- FindClusters(fib, resolution = 0.5, verbose = FALSE)
fib <- FindClusters(fib, resolution = 0.7, verbose = FALSE)
fib <- FindClusters(fib, resolution = 0.8, verbose = FALSE)
fib <- FindClusters(fib, resolution = 1, verbose = FALSE)


fib <- RunUMAP(fib, reduction = "pca", dims = 1:30)
DimPlot(fib, group.by = "assay")
p2 <- DimPlot(fib, split.by = "assay", group.by = 'integrated_snn_res.0.5', label = T)
p2

Idents(fib)<- 'integrated_snn_res.0.5'
DefaultAssay(fib) <- 'RNA'
fib.markers <- FindAllMarkers(fib, only.pos = T)

DimPlot(fib, label = T, group.by = 'celltype_broad')

Idents(fib) <- 'integrated_snn_res.0.5'
DimPlot(fib, label = T, repel = T)

Idents(fib) <- 'integrated_snn_res.0.5'
fib <- RenameIdents(fib, "0" = "FIB",
                    "1" = "FIB_Reticular",
                    "2" = "FIB_Papillary",
                    "3" = "FIB_Subcutis",
                    "4"= "FIB_DS",
                    "5"= "FIB_Subfascia",
                    "6" = "FIB",
                    "7" = "FIB_Reticular",
                    "8" = "FIB_DP",
                    "9" = "FIB_Reticular",
                    "10" = "SMC",
                    "11" = "Doublet")

DimPlot(fib, label = T, repel = T, split.by = 'assay')

fib$celltype_subtype <- as.character(Idents(fib))

DefaultAssay(fib) <- 'RNA'
fib.markers <- FindAllMarkers(fib, only.pos = T)

sort(unique(fib$celltype_subtype))

fib_markers <- list("SMC" = c("Des","Sfrp2"),
                    "FIB_Subfascia" = c("Plau","Ifit3"),
                    "FIB_Subcutis" = c("F3","Pi16"),
                    "FIB_Reticular" = c("Il1r2","Ccl11"),
                    "FIB_Papillary" = c("Ccbe1","Tnxb"),
                    "DS" = c("Col11a1","Lrrc15"),
                    "DP"= c("Crabp1","Serpine2"),
                    "FIB"= c("Col1a1","Col4a1"),
                    "Doublet" = c("H2-K1"))
                    
fib_markers_rev <- rev(fib_markers)

Idents(fib) <- 'celltype_subtype'
levels(fib) <- rev(sort(as.character(unique(Idents(fib)))))

p1 <- DotPlot(fib, features = fib_markers_rev, scale.max = 50) + theme_bw() + theme(panel.spacing = unit(0.2, "lines")) + rotate_x_text(45) +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=0))+
  theme(strip.text = element_blank())+
  guides(color = guide_colorbar(title = "Average\nExpression", frame.colour = 'black'),
         size = guide_legend(title = "Percent\nExpressed")) +
  theme(legend.position = 'bottom',
        panel.grid.minor = element_line(color = "black", size = 0.5))
ggsave("fib_DotPlot_labeled.pdf",p1, width = 5.4, height = 4, units = 'in')




obj@meta.data[rownames(fib@meta.data),]$celltype_subtype <- fib@meta.data$celltype_subtype

obj@meta.data$fib_sub_res_0.5 <- "NA"
obj@meta.data[rownames(fib@meta.data),]$fib_sub_res_0.5 <- as.character(fib@meta.data$integrated_snn_res.0.5)

DimPlot(obj, group.by = 'celltype_subtype', split.by = 'condition', label = T, repel = T)
DimPlot(
  fib,
  reduction = "spatial",
  cells.highlight = WhichCells(fib, idents = 3),
  cols.highlight = 'red',
  cols = 'grey',
  sizes.highlight = 0.1,
  pt.size = 0.1, 
  raster = F
)

############## Lymph Recluster ########
Idents(obj) <- 'celltype_broad'
lymph <- subset(obj, idents = c("Lymphocyte"))
DefaultAssay(lymph) <- 'integrated'
lymph <- ScaleData(lymph, verbose = TRUE, vars.to.regress = c("nCount_RNA"))
lymph <- RunPCA(lymph, npcs = 30, verbose = TRUE)

lymph <- FindNeighbors(lymph, dims = 1:30, verbose = T)
lymph <- FindClusters(lymph, resolution = 0.2, verbose = FALSE)
lymph <- FindClusters(lymph, resolution = 0.3, verbose = FALSE)
lymph <- FindClusters(lymph, resolution = 0.4, verbose = FALSE)
lymph <- FindClusters(lymph, resolution = 0.5, verbose = FALSE)
lymph <- FindClusters(lymph, resolution = 0.6, verbose = FALSE)
lymph <- FindClusters(lymph, resolution = 0.8, verbose = FALSE)
lymph <- FindClusters(lymph, resolution = 1, verbose = FALSE)


lymph <- RunUMAP(lymph, reduction = "pca", dims = 1:30)
DimPlot(lymph, group.by = "assay")
p2 <- DimPlot(lymph, split.by = "assay", group.by = 'integrated_snn_res.0.6', label = T)
p2

Idents(lymph)<- 'integrated_snn_res.0.6'
DefaultAssay(lymph) <- 'RNA'
lymph.markers <- FindAllMarkers(lymph, only.pos = T)

immune_features <- list( Lymphocyte = c("Cd3d","Cd3g","Il2rg","Cd2","Ptprc"),
                         TC_gdt =  c("Cd3e","Trdc","Kcnma1"),
                         TC_Th17 = c("Il23r","Rorc","Cxcr6"),
                         TC_Treg = c("Foxp3","Ctla4"),
                         TC_Cd4 = c("Cd40lg","Cd4","Ccr4"),
                         TC_Cd8 = c("Cd8a","Csf1","Ccl4"),
                         TC_Cyc = c("Mki67","Top2a"),
                         ILC2 = c("Gata3","Il7r","Kit"),
                         NK = c("Ncr1","Gzma","Irf8"))

Idents(lymph) <- 'integrated_snn_res.0.6'
lymph <- RenameIdents(lymph, "0" = "TC_gdt",
                      "1" = "TC_Cyc",
                      "2" = "Doublet",
                      "3" = "TC_Th17",
                      "4"= "ILC2",
                      "5"= "NK_Cd8",
                      "6" = "Doublet2",
                      "7"= "DC_cDC2",
                      "8"= "TC_Treg")

DimPlot(lymph, label = T, split.by = 'assay')

lymph$celltype_subtype <- as.character(Idents(lymph))
Idents(lymph) <- 'celltype_subtype'

DefaultAssay(lymph) <- 'RNA'
lymph.markers <- FindAllMarkers(lymph, only.pos = T)

levels(lymph) <- rev(sort(as.character(unique(Idents(lymph)))))

p1 <- DotPlot(lymph, features = lymph_markers_rev, scale.max = 50) + theme_bw() + theme(panel.spacing = unit(0.2, "lines")) + rotate_x_text(45) +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=0))+
  theme(strip.text = element_blank())+
  guides(color = guide_colorbar(title = "Average\nExpression", frame.colour = 'black'),
         size = guide_legend(title = "Percent\nExpressed")) +
  theme(legend.position = 'bottom',
        panel.grid.minor = element_line(color = "black", size = 0.5))


immune_features <- list( "DC_cDC2" = c("Flt3","Cd209a","Sirpa"),
                         "Doublet"= c("Col17a1","Dst"),
                         "Doublet2" = c("Col6a1","C3"),
                         "ILC2" = c("Kit","Il13","Il9r"),
                         "NK_Cd8"= c("Ncr1","Cd8a","Gzmb"),
                         "TC_Cyc"= c("Mki67","Top2a"),
                         "TC_gdt" =  c("Trdc","Kcnma1"),
                         "TC_Th17" = c("Il23r","Rorc","Cxcr6"),
                         "TC_Treg" = c("Foxp3","Ctla4"),
                         "TC" = c("Cd3d","Cd3g","Il2rg","Cd2"))
                         

order <- names(immune_features)
lymph_markers_rev <- rev(immune_features)

Idents(lymph) <- 'celltype_subtype'
levels(lymph) <- sort(as.character(unique(Idents(lymph))))

p1 <- DotPlot(lymph, features = lymph_markers_rev, scale.max = 50) + theme_bw() + theme(panel.spacing = unit(0.2, "lines")) + rotate_x_text(45) +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=0))+
  theme(strip.text = element_blank())+
  guides(color = guide_colorbar(title = "Average\nExpression", frame.colour = 'black'),
         size = guide_legend(title = "Percent\nExpressed")) +
  theme(legend.position = 'bottom',
        panel.grid.minor = element_line(color = "black", size = 0.5))
ggsave("Lymph_DotPlot_labeled.pdf",p1, width = 6.4, height = 4, units = 'in')




obj@meta.data[rownames(lymph@meta.data),]$celltype_subtype <- lymph@meta.data$celltype_subtype

obj$lymph_sub_res_0.6 <- "NA"
obj@meta.data[rownames(lymph@meta.data),]$lymph_sub_res_0.6 <- lymph@meta.data$integrated_snn_res.0.6


DimPlot(obj, group.by = 'celltype_subtype', split.by = 'condition', label = T, repel = T)


saveRDS(fib, "20260622_sc_merfish_int_fibs_annotated.rds")
saveRDS(myl, "20260622_sc_merfish_int_myl_annotated.rds")
saveRDS(lymph, "20260622_sc_merfish_int_lymph_annotated.rds")


Idents(obj) <- 'celltype_subtype'
sort(unique(obj$celltype_subtype))
obj <- subset(obj, idents = c('Doublet','Doublet2'), invert = T)

test <- names(detailed.cell_type.colors)
test %in% unique(obj$celltype_subtype)
unique(obj$celltype_subtype) %in%  test 

saveRDS(obj,"20260622_sc_merfish_int_full_obj_annotated_full.rds")

obj <- readRDS("20260622_sc_merfish_int_full_obj_annotated_full.rds")

Idents(obj) <- 'celltype_broad'
levels(obj) <- sort(unique(obj@meta.data$celltype_broad))
p2 <- DimPlot(obj, cols = colors, split.by = 'assay') & NoAxes() & NoLegend() 
p2 <- LabelClusters(plot = p2, id = 'ident', fontface = 'bold', repel = T, size = 5)
png("Full_labeled_Dimplot_splitby_assay_12_6.png", res = 800, width = 10, height = 5, units = 'in')
print(p2)
dev.off()


Idents(obj) <- 'celltype_subtype'
levels(obj) <- sort(unique(obj@meta.data$celltype_subtype))
p1 <- DimPlot(obj, cols = detailed.cell_type.colors, label = T, label.box = T, repel = T, raster = F) & NoAxes() & NoLegend()
png("Full_labeled_Dimplot.png", res = 800, width = 5.3, height = 5.3, units = 'in')
print(p1)
dev.off()



p2 <- DimPlot(obj, cols = detailed.cell_type.colors, split.by = 'assay') & NoAxes() & NoLegend() 
p2 <- LabelClusters(plot = p2, id = 'ident', fontface = 'bold', repel = T, size = 5)
png("Full_labeled_Dimplot_splitby_assay_12_6.png", res = 800, width = 12, height = 6, units = 'in')
print(p2)
dev.off()



p2 <- DimPlot(obj, cols = detailed.cell_type.colors) & NoAxes()
png("Full_labeled_legend_splitby_assay_12_6.png", res = 800, width = 12, height = 6, units = 'in')
print(p2)
dev.off()

p2 <- DimPlot(obj, cols = detailed.cell_type.colors) & NoAxes() & NoLegend() 
p2 <- LabelClusters(plot = p2, id = 'ident', fontface = 'bold', repel = T, size = 5)
png("Full_labeled_Dimplot_6_6.png", res = 800, width = 6, height = 6, units = 'in')
print(p2)
dev.off()



# create spatial coordinate matrix
coords <- cbind(
  as.numeric(obj@meta.data$center_x),
  as.numeric(obj@meta.data$center_y)
)

colnames(coords) <- c("SPATIAL_1", "SPATIAL_2")
rownames(coords) <- rownames(obj@meta.data)

obj[["spatial"]] <- CreateDimReducObject(
  embeddings = coords,
  key = "SPATIAL_",
  assay = DefaultAssay(obj)
)

Idents(obj) <- 'celltype_subtype'
levels(obj)  <- sort(unique(obj$celltype_subtype)) 
p3 <- DimPlot(obj, red = 'spatial', cols = detailed.cell_type.colors,pt.size = 0.001, raster = F, split.by = 'condition') + coord_fixed() &NoAxes() &NoLegend()
png("Full_labeled_SPATIAL_splitby_assay_12_6.png", res = 1200, width = 24, height = 12, units = 'in')
print(p3)
dev.off()

pdf("Full_labeled_SPATIAL_splitby_assay_12_6.pdf", width = 24, height = 12)   # inches
print(p3)
dev.off()


p1 <- DimPlot(obj, cols = detailed.cell_type.colors, label = T, label.box = T, repel = T, raster = F) & NoAxes() 
pdf("Dimplot_legend_subtype.pdf", width = 12, height = 5)
print(p1)
dev.off()



setwd("D:/OneDrive/OneDrive - The Mount Sinai Hospital/Ji_Lab_2023/RG_MERFISH_files/20250710_WT_JAK1_only/labeled_intres1.1_spatial_highlights_labeled/celltype_plots")

outdir <- "celltype_plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
Idents(obj) <- "celltype_subtype"
# Colors you can change
bg_grey <- "grey80"
highlight_col <- "#B2182B"  # dark red

celltypes <- unique(obj$celltype_subtype)

for (ct in celltypes) {
  safe_name <- gsub("[^A-Za-z0-9_.-]", "_", ct)
  cells_use <- WhichCells(obj, expression = celltype_subtype == ct)
  
  if (length(cells_use) == 0) {
    message(sprintf("Skipping '%s' (no cells found).", ct))
    next
  }
  
  # --- get spatial coordinates ---
  coords <- tryCatch(
    Embeddings(obj, "spatial"),
    error = function(e) NULL
  )
  
  if (is.null(coords)) {
    stop("Cannot find spatial embeddings via Embeddings(obj, 'spatial'). Please adapt to where your spatial coordinates live.")
  }
  
  # build plotting dataframe
  df <- as.data.frame(coords)
  df$cell <- rownames(df)
  xcol <- colnames(df)[1]
  ycol <- colnames(df)[2]
  
  # robustly get metadata and align rows to df$cell
  meta <- obj@meta.data
  if (all(df$cell %in% rownames(meta))) {
    meta_rows <- match(df$cell, rownames(meta))
  } else if (all(df$cell %in% colnames(obj))) {
    meta_rows <- match(df$cell, colnames(obj))
  } else {
    stop("Cannot align cell names between embeddings and metadata. Check your object cell names.")
  }
  
  # check condition column
  if (!"condition" %in% colnames(meta)) {
    stop("No 'condition' column found in obj@meta.data. Change 'condition' in the code to your column name.")
  }
  # keep factor levels from the full metadata so facets preserve order
  condition_levels <- if (is.factor(meta[ , "condition"])) levels(meta[ , "condition"]) else unique(as.character(meta[ , "condition"]))
  df$condition <- factor(as.character(meta[meta_rows, "condition"]), levels = condition_levels)
  
  # create group column and ensure ordering: other first, highlight last
  df$group <- ifelse(df$cell %in% cells_use, "highlight", "other")
  df$group <- factor(df$group, levels = c("other", "highlight"))
  
  # split data to avoid NSE issues and guarantee draw order
  df_bg <- df[df$group == "other", , drop = FALSE]
  df_hi <- df[df$group == "highlight", , drop = FALSE]
  
  # plotting parameters (tweak to taste)
  point_size <- 0.5     # ggplot2 point size
  bg_alpha <- 1         # 1 = opaque; set <1 to fade background
  
  # ggplot using .data pronoun for tidy-eval
  p <- ggplot() +
    geom_point(
      data = df_bg,
      mapping = aes(x = .data[[xcol]], y = .data[[ycol]]),
      colour = bg_grey,
      size = point_size,
      alpha = bg_alpha
    ) +
    geom_point(
      data = df_hi,
      mapping = aes(x = .data[[xcol]], y = .data[[ycol]]),
      colour = highlight_col,
      size = point_size,
      alpha = 1
    ) +
    facet_wrap(~ condition) +
    coord_fixed(expand = FALSE) +
    theme_void() +
    theme(
      strip.text = element_text(size = 12),
      panel.spacing = unit(0.5, "lines")
    )
  
  # save
  png_file <- file.path(outdir, paste0(safe_name, ".png"))
  ggsave(png_file, plot = p, device = "png", width = 24, height = 12, dpi = 800)
}




p1 <- DimPlot(myl, cols = detailed.cell_type.colors, repel = T) & NoAxes() 
p1 <- LabelClusters(plot = p1, id = 'ident', fontface = 'bold', repel = T, size = 5)

p2 <- DimPlot(fib, cols = detailed.cell_type.colors, repel = T) & NoAxes()  
p2 <- LabelClusters(plot = p2, id = 'ident', fontface = 'bold', repel = T, size = 5)

p3 <- DimPlot(lymph, cols = detailed.cell_type.colors, repel = T) & NoAxes()
p3 <- LabelClusters(plot = p3, id = 'ident', fontface = 'bold', repel = T, size = 5)

library(patchwork)

combined <- p2 | p1 | p3

combined

ggsave(
  filename = "20260622_MERFISH_Subcluster_UMAPs.pdf",
  plot = combined,
  width = 15,
  height = 4,
  units = "in"
)
