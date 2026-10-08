library(Seurat)
library(DoubletFinder)

#import sc rds
setwd("D:/OneDrive/OneDrive - The Mount Sinai Hospital/Ji_Lab_2023/RG_MERFISH_files/20250710_WT_JAK1_only")

##Import sc data
II1.data <- Read10X_h5("D:/OneDrive/OneDrive - The Mount Sinai Hospital/RG_sc_files/20240216_II_Mouse_AD/Raw_data/II1_outs/II1/outs/filtered_feature_bc_matrix.h5")
# Initialize the Seurat object with the raw (non-normalized data).
II1 <- CreateSeuratObject(counts = II1.data, project = "H_MC903_OXA", min.cells = 3)

II2.data <- Read10X_h5("D:/OneDrive/OneDrive - The Mount Sinai Hospital/RG_sc_files/20240216_II_Mouse_AD/Raw_data/II2/raw_feature_bc_matrix.h5")
# Initialize the Seurat object with the raw (non-normalized data).
II2 <- CreateSeuratObject(counts = II2.data, project = "OXA_H_JAK1", min.cells = 3)

II1@meta.data$condition <- "WT"
II2@meta.data$condition <- "JAK1GoF"

sc <- merge(II1, y = c(II2), add.cell.ids = c("WT","JAK1GoF"), project = "mouse_AD")

##QC 
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^mt-")
# Visualize QC metrics as a violin plot
#VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'condition')
sc <- subset(sc, subset = nFeature_RNA > 200 & percent.mt)
VlnPlot(sc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'condition')

do.seurat <- function(obj) {
  require(Seurat)
  obj <- JoinLayers(obj)
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.1, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.2, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:30, verbose = TRUE)
}

sc <- do.seurat(sc)

Idents(sc) <- 'RNA_snn_res.0.1'
DimPlot(sc, label = T)
#markers <- FindAllMarkers(sc, only.pos = T)
sc <- subset(sc, idents = "0", invert = T) #remove low transcript/feature cluster, keep cluster 16 for part of mast cell

library(scDblFinder)
sce <- as.SingleCellExperiment(sc)
sce <- scDblFinder(sce)
table(sce$scDblFinder.class)

# Transfer singlet/doublet calls back to Seurat metadata
sc$scDblFinder.class <- sce$scDblFinder.class
sc$scDblFinder.score <- sce$scDblFinder.score  # optional: keep scores too

# Now filter singlets directly in Seurat
sc <- subset(sc, subset = scDblFinder.class == "singlet")
sc <- do.seurat(sc)

FeaturePlot(sc, features = 'nFeature_RNA', split.by = 'condition', max.cutoff = 2500) & theme(legend.position = "right")
FeaturePlot(sc, features = 'nCount_RNA', split.by = 'condition', max.cutoff = 2500) & theme(legend.position = "right")

table(sc@meta.data$condition)

FeaturePlot(sc, features = c("Col1a1","Dcn","Lum","Col11a1"))
Idents(sc) <- 'RNA_snn_res.0.1'
DimPlot(sc, label = T, split.by = 'condition')
FeaturePlot(sc, features = c("Col1a1","Dcn","Lum","Col11a1"))
fibs <- subset(sc, idents = c("1"))
fibs <- do.seurat(fibs)
Idents(fibs) <- 'RNA_snn_res.0.2'
DimPlot(fibs, label = T, split.by = 'condition')

fib.markers <- FindAllMarkers(fibs, only.pos = T)
FeaturePlot(fibs, features = c("Col1a1","Dcn","Lum","Col11a1"))

Idents(fibs)<- 'condition'
fibs_deg <- FindAllMarkers(fibs, only.pos = T, min.pct = 0.10)

DotPlot(fibs, features = c("Tnc","Postn","Cxcl1","Ccl7","Ccl4","Cxcl14"), group.by = 'condition')


p1 <- DimPlot(fibs, label = F, cols = c("grey65","darkorchid")) & NoAxes()
png("2025117_WT_JAK1_fib_sc_Only_subcluster_dimplot.png", height = 3, width = 3.5, units = 'in', res = 600)
print(p1)
dev.off()

saveRDS(fibs, "20251117_sc_only_WT_JAK1_fibroblast.rds")

DotPlot(fibs, features = c("Tnc","Postn","Cxcl1","Ccl7", "Cxcl14","Ccl4"), group.by = 'condition')

p1 <- DotPlot(fibs, features = c("Tnc","Postn","Cxcl1","Ccl7","Ccl4", "Cxcl14"), scale.max = 80) + theme_bw() + theme(panel.spacing = unit(0.2, "lines")) + rotate_x_text(45) +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=0))
pdf("WT_JAK1_sc_fibs_DotPlot_unlabeled.pdf", width = 4, height = 2)
print(p1)
dev.off()


merfish <- readRDS("20250710_MERFISH_JAK1_WT_baysor_noint.rds")
do.seurat <- function(obj) {
  require(Seurat)
  obj <- JoinLayers(obj)
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.1, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.2, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:30, verbose = TRUE)
}

merfish <- FindClusters(merfish, resolution = 1, verbose = FALSE)
DimPlot(merfish, label = T)
celltype_features <- list( Adipo =  c("Adipoq","Ccdc80"),
                           Baso = c("Mcpt8"),
                           DC = c("Ccr7","Cd68","Ccr2","H2-Aa"),
                           DP = c("Crabp1","Dkk2"),
                           Fib = c("Col1a1","Col3a1","Tnc","Cxcl5"),
                           HF = c("Sox9","Sfrp1"),
                           KC = c("Krt1","Krt6a"),
                           LEC = c("Cldn5","Lyve1"),
                           Mac = c("Mrc1","Cd163"),
                           Mast = c("Mcpt4","Cpa3"),
                           Neutrophil = c("Cxcr2","S100a8"),
                           Pericyte = c("Rgs5","Acta2"),
                           Schwann = c("Mpz","Sox10"),
                           SM = c("Des","Vegfa"),
                           Tc = c("Cd2","Cd3g"),
                           VEC = c("Sele","Pecam1"))

DotPlot(merfish, features = celltype_features)

fibs <- subset(merfish, idents = c("1","2","9","25","23","24"))
fibs <- do.seurat(fibs)
fibs <- FindClusters(fibs, resolution = 0.6, verbose = FALSE)
DimPlot(fibs, label = T, group.by = 'RNA_snn_res.0.6')
DotPlot(fibs, group.by = 'RNA_snn_res.0.6', features = celltype_features)

##remove mac doublet cluster
fibs <- subset(fibs, idents = c("10"), invert = T)
fibs <- do.seurat(fibs)
fibs <- FindClusters(fibs, resolution = 0.5, verbose = FALSE)
DimPlot(fibs, label = T, group.by = 'RNA_snn_res.0.5')
DotPlot(fibs, group.by = 'RNA_snn_res.0.5', features = celltype_features)

Idents(fibs) <- 'condition'
p1 <- DotPlot(fibs, features = c("Tnc","Postn","Cxcl1","Ccl7","Ccl8","Ccl4", "Cxcl14"), scale.max = 10) + theme_bw() + theme(panel.spacing = unit(0.2, "lines")) + rotate_x_text(45) +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=0))

p1 <- DimPlot(fibs, label = F, cols = c("grey65","darkorchid")) & NoAxes()
png("2025117_WT_JAK1_fib_merfish_Only_subcluster_dimplot.png", height = 3, width = 3.5, units = 'in', res = 600)
print(p1)
dev.off()

Idents(fibs) <- 'condition'
fibs <- RenameIdents(fibs, "JAK1"= "JAK1GoF",
                     "Normal" = "WT")
fibs$condition <- Idents(fibs)

p1 <- DimPlot(fibs, label = F, cols = c("darkorchid","grey65")) & NoAxes()
png("2025117_WT_JAK1_fib_merfish_Only_subcluster_dimplot.png", height = 3, width = 3.5, units = 'in', res = 600)
print(p1)
dev.off()

fib_marker_table <- FindAllMarkers(fibs, only.pos = T)
levels(fibs) <- c("WT","JAK1GoF")
p1 <- DotPlot(fibs, features = c("Tnc","Postn","Cxcl1","Ccl8","Ccl4", "Cxcl14"), scale.max = 15) + theme_bw() + theme(panel.spacing = unit(0.2, "lines")) + rotate_x_text(45) +
  theme(axis.text.x = element_text(size=12, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title = element_text(size=0))
pdf("WT_JAK1_merfish_fibs_DotPlot_unlabeled.pdf", width = 4, height = 2)
print(p1)
dev.off()