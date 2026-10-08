library(Seurat)
library(CellChat)

setwd("/sc/arion/projects/ji_lab/USER/gillr10/Analysis/JAK1_cellchat")

obj <- readRDS("/sc/arion/projects/ji_lab/USER/gillr10/Analysis/JAK1_cellchat/20260107_sc_merfish_int_full_obj_annotated_full.rds")

Idents(obj) <- 'assay'
sc <- subset(obj, idents = 'scRNA') 
rm(obj)
DefaultAssay(sc) <- 'RNA'
sc <- NormalizeData(sc)

Idents(sc) <- 'celltype_subtype'

sc <- RenameIdents(sc, 
  # KC populations
  "KC_Bas" = "KC",
  "KC_Grn" = "KC",
  "KC_IFE" = "KC",
  "KC_Spn" = "KC",
  
  # FIB populations
  "FIB"           = "FIB",
  "FIB_Cyc"       = "FIB",
  "FIB_DP"        = "FIB",
  "FIB_Papillary" = "FIB",
  "FIB_Reticular" = "FIB",
  "FIB_Subcutis"  = "FIB",
  "FIB_Subfascia" = "FIB",
  
  # HF populations
  "HF_IRS"        = "HF",
  "HF_ORS_Bulge"  = "HF"
)

sc$celltype_cellchat <- as.character(Idents(sc))

sc <- SetIdent(sc, value = sc$condition)
sc.Normal <- subset(sc, idents = "Normal")
sc.Normal <- NormalizeData(sc.Normal)
sc.JAK1 <- subset(sc, idents = "JAK1GoF")
sc.JAK1 <- NormalizeData(sc.JAK1)
sc.Normal <- SetIdent(sc.Normal, value = sc.Normal@meta.data$celltype_cellchat)
sc.JAK1 <- SetIdent(sc.JAK1, value = sc.JAK1@meta.data$celltype_cellchat)
sc_samples <- c("sc.Normal","sc.JAK1")
sc_samples.cellchat <- c("sc.Normal.cellchat","sc.JAK1.cellchat")
sc.Normal.cellchat <- createCellChat(object = sc.Normal, group.by = 'celltype_cellchat', assay = 'RNA')
sc.JAK1.cellchat <- createCellChat(object = sc.JAK1, group.by = 'celltype_cellchat', assay = 'RNA')


#add db info
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB)

sc.Normal.cellchat@DB <- CellChatDB.use
sc.JAK1.cellchat@DB <- CellChatDB.use

future::plan("multisession", workers = 10) # do parallel
options(future.globals.maxSize = 500000 * 1024^2)

########################################## run cellchat on each condition ########################################## 

sc.Normal.cellchat <- subsetData(sc.Normal.cellchat)
sc.Normal.cellchat <- identifyOverExpressedGenes(sc.Normal.cellchat)
sc.Normal.cellchat <- identifyOverExpressedInteractions(sc.Normal.cellchat)
sc.Normal.cellchat <- computeCommunProb(sc.Normal.cellchat, type =  "truncatedMean", trim = 0.05, population.size = F)
sc.Normal.cellchat <- computeCommunProbPathway(sc.Normal.cellchat)
sc.Normal.cellchat <- aggregateNet(sc.Normal.cellchat)
sc.Normal.cellchat <- netAnalysis_computeCentrality(sc.Normal.cellchat, slot.name = "netP")
saveRDS(sc.Normal.cellchat, file = '20260109_Normal_sc_only_cellchat_v2.rds')
