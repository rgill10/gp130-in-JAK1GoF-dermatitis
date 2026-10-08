library(Seurat)
library(future)

# Joint embedding UMAP of MERFISH + scRNA-seq of OXA

#import sc rds
setwd("D:/OneDrive/OneDrive - The Mount Sinai Hospital/Ji_Lab_2023/RG_MERFISH_files/20250710_WT_JAK1_only")

merfish <- readRDS("20250710_MERFISH_JAK1_WT_baysor_noint.rds")

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
markers <- FindAllMarkers(sc, only.pos = T)
sc <- subset(sc, idents = "0", invert = T) #remove low transcript/feature cluster, keep cluster 16 for part of mast cell

sc <- do.seurat(sc)

Idents(sc) <- 'RNA_snn_res.0.1'
DimPlot(sc, label = T, split.by = 'condition')

FeaturePlot(sc, features = 'nFeature_RNA', split.by = 'condition', max.cutoff = 2500) & theme(legend.position = "right")
FeaturePlot(sc, features = 'nCount_RNA', split.by = 'condition', max.cutoff = 2500) & theme(legend.position = "right")

table(merfish@meta.data$condition)
table(sc@meta.data$condition)

sc@meta.data$assay <- "scRNA"
merfish@meta.data$assay <- "MERFISH"

# Correct gene names in merfish data
merfish_new = merfish@assays$RNA$counts
rownames(merfish_new) = gsub("H2.Aa","H2-Aa",rownames(merfish_new))
rownames(merfish_new) = gsub("H2.K1","H2-K1",rownames(merfish_new))
merfish_new = CreateSeuratObject(merfish_new) # new MERFISH object
merfish_new = AddMetaData(merfish_new, metadata = merfish@meta.data)

merfish_genes <- rownames(merfish_new@assays$RNA$counts)
sc_genes <- rownames(sc@assays$RNA$counts)
overlap_genes = intersect(rownames(merfish_new@assays$RNA$counts),rownames(sc@assays$RNA$counts))
missing <- merfish_genes [!merfish_genes  %in% overlap_genes]

merged <- merge(merfish, sc)
saveRDS(merged, file = "20250711_merfish_sc_merged_noint.rds")

merged <- readRDS("20250711_merfish_sc_merged_noint.rds")
plan()
plan("multisession", workers = 4)
options(future.globals.maxSize = 18 * 1024^3)
plan()

table(merged@meta.data$assay)
colnames(merged@meta.data)

Idents(merged) <- "assay"
DefaultAssay(merged) <- 'RNA'
merged <- JoinLayers(merged)
merged <- NormalizeData(merged)
merged <- ScaleData(merged, vars.to.regress = "nCount_RNA")

Idents(merged) <- 'assay'
merfish <- subset(merged, idents = 'MERFISH')
merfish_new = merfish@assays$RNA$counts
rownames(merfish_new) = gsub("H2.Aa","H2-Aa",rownames(merfish_new))
rownames(merfish_new) = gsub("H2.K1","H2-K1",rownames(merfish_new))
merfish_new = CreateSeuratObject(merfish_new) # new MERFISH object
merfish_new = AddMetaData(merfish_new, metadata = merfish@meta.data)
merfish_genes <- rownames(merfish_new@assays$RNA$counts)

ns.list <- SplitObject(merged, split.by = "assay")
rm(merged)

saveRDS(ns.list,"20250714_objectlist_norm.rds")

ns.list[[1]] <- NormalizeData(ns.list[[1]], verbose = T)
ns.list[[1]] <- FindVariableFeatures(ns.list[[1]], selection.method = "vst",
                                     nfeatures = 498, verbose = T)
ns.list[[2]] <- NormalizeData(ns.list[[2]], verbose = T)
ns.list[[2]] <- FindVariableFeatures(ns.list[[2]], selection.method = "vst",
                                     nfeatures = 498, verbose = T)


ns.anchors <- FindIntegrationAnchors(object.list = ns.list, dims = 1:30, anchor.features = merfish_genes, verbose = T)

saveRDS(ns.anchors, file = "20241119_int_anchors.rds")
sc_merfish <- IntegrateData(anchorset = ns.anchors, dims = 1:30)
saveRDS(sc_merfish, file = "20241119_sc_merfish_int_full_obj.rds")

DefaultAssay(sc_merfish) <- "integrated"
rm(ns.list)
rm(ns.anchors)

sc_merfish <- ScaleData(sc_merfish, verbose = TRUE, vars.to.regress = c("nCount_RNA"))
sc_merfish <- RunPCA(sc_merfish, npcs = 30, verbose = TRUE)

sc_merfish <- FindNeighbors(sc_merfish, dims = 1:30, verbose = T)
sc_merfish <- FindClusters(sc_merfish, resolution = 0.4, verbose = FALSE)
sc_merfish <- FindClusters(sc_merfish, resolution = 0.6, verbose = FALSE)
sc_merfish <- FindClusters(sc_merfish, resolution = 0.8, verbose = FALSE)
sc_merfish <- FindClusters(sc_merfish, resolution = 1, verbose = FALSE)
sc_merfish <- FindClusters(sc_merfish, resolution = 1.1, verbose = FALSE)
sc_merfish <- FindClusters(sc_merfish, resolution = 1.2, verbose = FALSE)

sc_merfish <- RunUMAP(sc_merfish, reduction = "pca", dims = 1:30)
DimPlot(sc_merfish, group.by = "assay")
DimPlot(sc_merfish, split.by = "assay")

saveRDS(sc_merfish, file = "20240711_sc_merfish_int_full_obj.rds")

sc_merfish <- readRDS("20241119_sc_merfish_int_full_obj.rds")
DefaultAssay(sc_merfish) <- 'RNA'
sc_merfish <- JoinLayers(sc_merfish)
Idents(sc_merfish) <- 'integrated_snn_res.1'
markers <- FindAllMarkers(sc_merfish, only.pos = T)
write.csv(markers, file = "Int_res.1_cluster_markers.csv")
Idents(sc_merfish) <- 'integrated_snn_res.1.1'
markers <- FindAllMarkers(sc_merfish, only.pos = T)
write.csv(markers, file = "Int_res.1.1_cluster_markers.csv")
Idents(sc_merfish) <- 'integrated_snn_res.1.2'
markers <- FindAllMarkers(sc_merfish, only.pos = T)
write.csv(markers, file = "Int_res.1.2_cluster_markers.csv")


