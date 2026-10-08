library(Seurat)

obj <- readRDS("D:/OneDrive/OneDrive - The Mount Sinai Hospital/Ji_Lab_2023/RG_MERFISH_files/20241111_mouse_merfish_merge_WT_OXA_MC903_timecourse/More_rds/20241111_MC903_cKO_MAR1_timecourse_WT_MC903_only.rds")

unique(obj@meta.data$sample)

Idents(obj) <- 'sample'
wt <- subset(obj, idents = c("d0_WT_neck"))
unique(wt@meta.data$orig.ident)
wt@meta.data$orig.ident <- "R1_WT_Neck"
rm(obj)

dirs <- list.dirs(path= 'D:/OneDrive/OneDrive - The Mount Sinai Hospital/Ji_Lab_2023/RG_MERFISH_files/20250710_WT_JAK1_only/MERFISH_baysor_outs', recursive = FALSE, full.names=FALSE)
for (x in dirs) {
  name <- x
  counts <- read.table(paste0("D:/OneDrive/OneDrive - The Mount Sinai Hospital/Ji_Lab_2023/RG_MERFISH_files/20250710_WT_JAK1_only/MERFISH_baysor_outs/",x,"/baysor_cell_by_gene.csv"),sep = ",",header=T,row.names=1,colClasses = "character")
  metadata <- read.table(paste0("D:/OneDrive/OneDrive - The Mount Sinai Hospital/Ji_Lab_2023/RG_MERFISH_files/20250710_WT_JAK1_only/MERFISH_baysor_outs/",x,"/baysor_cell_metadata.csv"),sep = ",",header=T,row.names=1,colClasses = "character")
  cell_metadata_filt = metadata[rownames(counts),]
  assign(name, CreateSeuratObject(counts = t(counts[,1:500]),meta.data = cell_metadata_filt, project = name))
  rm(counts)
  rm(metadata)
  rm(cell_metadata_filt)
}

meta <- colnames(wt@meta.data)
select <- meta[1:15]
wt@meta.data <- wt@meta.data[,select]


obj <- merge(wt, c(R0_JAK1_Ear, R1_JAK1_Inguinal,R3_WT_Ear,R4_WT_Neck, R6_JAK1_Neck, R8_WT_Flank))


do.seurat <- function(obj) {
  require(Seurat)
  obj <- JoinLayers(obj)
  DefaultAssay(obj) <- "RNA"
  obj <- subset(obj, subset = nCount_RNA > 10 & volume > 100)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nrow(obj))
  all.genes <- rownames(obj)
  obj <- ScaleData(obj, features = all.genes, vars.to.regress = c("nCount_RNA"))
  obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 30)
  obj <- FindNeighbors(obj, dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.6, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.8, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 1.5, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 1.8, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 2, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:30, verbose = TRUE)
}

obj <- do.seurat(obj)

setwd("D:/OneDrive/OneDrive - The Mount Sinai Hospital/Ji_Lab_2023/RG_MERFISH_files/20250710_WT_JAK1_only")
saveRDS(obj, "20250710_MERFISH_JAK1_WT_baysor_noint.rds")

DimPlot(obj, group.by = 'orig.ident')


##Add site and condition columns
obj@meta.data$site <- obj@meta.data$orig.ident
unique(obj@meta.data$site)
obj$site <- gsub("R0_", "", obj$site)
obj$site <- gsub("R1_JAK1", "JAK1", obj$site)
obj$site <- gsub("R3_", "", obj$site)
obj$site <- gsub("R4_", "", obj$site)
obj$site <- gsub("R6_", "", obj$site)
obj$site <- gsub("R8_", "", obj$site)
obj$site <- gsub("R1_WT_Neck", "Normal_Neck2", obj$site)
obj$site <- gsub("WT", "Normal", obj$site)
unique(obj@meta.data$site)
obj$condition <- sub("_.*", "", obj$site)
unique(obj@meta.data$condition)

saveRDS(obj, "20250710_MERFISH_JAK1_WT_baysor_noint.rds")





