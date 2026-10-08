library(CellChat)
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)

setwd("D:/OneDrive/OneDrive - The Mount Sinai Hospital/RG_sc_files/20250707_JAK1_GoF_Reanalysis/Cellchat")

jak <- readRDS("D:/OneDrive/OneDrive - The Mount Sinai Hospital/RG_sc_files/20250707_JAK1_GoF_Reanalysis/20260109_JAK1_sc_only_cellchat_v2.rds")
unique(jak@idents)
wt <- readRDS("D:/OneDrive/OneDrive - The Mount Sinai Hospital/RG_sc_files/20250707_JAK1_GoF_Reanalysis/20260109_Normal_sc_only_cellchat_v2.rds")
rankNet(jak, mode = "single")

object.list <- list(WT = wt, JAK1 = jak)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")


num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)


colors <- c("grey50","#7A3DB8")

detailed.celltype.colors <- c(
 "Basophil" ='#EE1289',
 "Chondrocyte" ="#D9A5BA",
 "DC1" ="chartreuse1",
 "DC2" ="chartreuse4",
 "FIB"="navy",
 "HF"='tan3',
 "ILC2"='purple',
 "KC" ="#358A8E",
 "LC" = "#556B2F",
 "LEC"= "#8B7355",
 "MAC"="#00BFFF",
 "Mast"="#8B4C39",
 "Neutrophil"="#8B668B",
 "NK"='steelblue',
 "Pericyte"="#708090",
 "Schwann"="#EE8066",
 "SMC"='#F82705',
 "TC_gdt"='firebrick4',
 "TC_Treg"='lightgoldenrod',
 "TC_Th17"='lightpink',
 "TC_Cd8"='chocolate',
 "TC_Cyc"='wheat2',
 "VEC" = "#D2B48C")


detailed.celltype.colors <- as.character(c(
  Basophil    = "#EE1289",
  Chondrocyte = "#D9A5BA",
  DC1         = grDevices::rgb(t(col2rgb("chartreuse1")), maxColorValue = 255),
  DC2         = grDevices::rgb(t(col2rgb("chartreuse4")), maxColorValue = 255),
  FIB         = grDevices::rgb(t(col2rgb("navy")), maxColorValue = 255),
  HF          = grDevices::rgb(t(col2rgb("tan3")), maxColorValue = 255),
  ILC2        = grDevices::rgb(t(col2rgb("purple")), maxColorValue = 255),
  KC          = "#358A8E",
  LC          = grDevices::rgb(t(col2rgb("darkolivegreen")), maxColorValue = 255),
  LEC         = grDevices::rgb(t(col2rgb("tan4")), maxColorValue = 255),
  MAC         = "#00BFFF",
  Mast        = grDevices::rgb(t(col2rgb("#8B4C39")), maxColorValue = 255),
  Neutrophil  = grDevices::rgb(t(col2rgb("#8B668B")), maxColorValue = 255),
  NK          = grDevices::rgb(t(col2rgb("steelblue")), maxColorValue = 255),
  Pericyte    = grDevices::rgb(t(col2rgb("slategray")), maxColorValue = 255),
  Schwann     = grDevices::rgb(t(col2rgb("#EE8066")), maxColorValue = 255),
  SMC         = "#F82705",
  TC_gdt      = grDevices::rgb(t(col2rgb("firebrick4")), maxColorValue = 255),
  TC_Treg     = grDevices::rgb(t(col2rgb("lightgoldenrod")), maxColorValue = 255),
  TC_Th17     = grDevices::rgb(t(col2rgb("lightpink")), maxColorValue = 255),
  TC_Cd8      = grDevices::rgb(t(col2rgb("chocolate")), maxColorValue = 255),
  TC_Cyc      = grDevices::rgb(t(col2rgb("wheat2")), maxColorValue = 255),
  VEC         = "#D2B48C"
))

##across whole dataset
for (i in 1:length(object.list)) {
  test=rankNet(object.list[[i]], mode = "single", stacked = T,return.data = T,do.stat = T, signaling.type = 'Secreted Signaling')
  test=test$signaling.contribution
  test <- test %>%arrange(desc(contribution.scaled))
  test <- head(test, 38)
  p<-ggplot(data=test, aes(x=name, y=contribution.scaled)) +
    geom_bar(stat="identity", fill = colors[i]) + coord_flip() + ggtitle(paste0(names(object.list)[i])) +
  pdf(paste0("20250108_",names(object.list)[i],"_pathways_rank_sc.pdf"), height =8.25, width = 2.4)
  print(p)
  dev.off()
}

#fib
for (i in 1:length(object.list)) {
  test=rankNet(object.list[[i]], mode = "single", stacked = T,return.data = T,do.stat = T, signaling.type = 'Secreted Signaling', 
               sources.use = c("MAC","Basophil","Mast","Neutrophil","TC_gdt","TC_Th17","DC1","DC2","ILC2","NK","TC_Cd8","TC_Treg"),
targets.use = c("FIB"))
  test=test$signaling.contribution
  test <- test %>%arrange(desc(contribution.scaled))
  test <- head(test, 20)
  p<-ggplot(data=test, aes(x=name, y=contribution.scaled)) +
    geom_bar(stat="identity", fill = colors[i]) + coord_flip() + ggtitle(paste0(names(object.list)[i])) + theme_classic()
  pdf(paste0("20250108_",names(object.list)[i],"_pathways_rank_sc_ImmunetoFIB.pdf"), height =5, width = 2.4)
  print(p)
  dev.off()
}

#kc
for (i in 1:length(object.list)) {
  colors <- c("grey50","#7A3DB8")
  test=rankNet(object.list[[i]], mode = "single", stacked = T,return.data = T,do.stat = T, signaling.type = 'Secreted Signaling',
sources.use = c("MAC","Basophil","Mast","Neutrophil","TC_gdt","TC_Th17","DC1","DC2","ILC2","NK","TC_Cd8","TC_Treg"),
targets.use = c("KC"))
  test=test$signaling.contribution
  test <- test %>%arrange(desc(contribution.scaled))
  test <- head(test, 20)
  p<-ggplot(data=test, aes(x=name, y=contribution.scaled)) +
    geom_bar(stat="identity", fill = colors[i]) + coord_flip() + ggtitle(paste0(names(object.list)[i])) + theme_classic()
  pdf(paste0("20250108_",names(object.list)[i],"_pathways_rank_sc_ImmunetoKC.pdf"), height =5, width = 2.4)
  print(p)
  dev.off()
}

#Immune
for (i in 1:length(object.list)) {
  colors <- c("grey50","#7A3DB8")
  test=rankNet(object.list[[i]], mode = "single", stacked = T,return.data = T,do.stat = T, signaling.type = 'Secreted Signaling',
sources.use = c("MAC","Basophil","Mast","Neutrophil","TC_gdt","TC_Th17","DC1","DC2","ILC2","NK","TC_Cd8","TC_Treg"),
targets.use = c("MAC","Basophil","Mast","Neutrophil","TC_gdt","TC_Th17","DC1","DC2","ILC2","NK","TC_Cd8","TC_Treg"))
  test=test$signaling.contribution
  test <- test %>%arrange(desc(contribution.scaled))
  test <- head(test, 20)
  p<-ggplot(data=test, aes(x=name, y=contribution.scaled)) +
    geom_bar(stat="identity", fill = colors[i]) + coord_flip() + ggtitle(paste0(names(object.list)[i])) + theme_classic()
  pdf(paste0("20250108_",names(object.list)[i],"_pathways_rank_sc_ImmunetoImmune.pdf"), height =5, width = 2.4)
  print(p)
  dev.off()
}


df=rankNet(object.list[[1]], mode = "single", stacked = T,return.data = T,do.stat = T, signaling.type = 'Secreted Signaling')
df=df$signaling.contribution

df=rankNet(object.list[[2]], mode = "single", stacked = T,return.data = T,do.stat = T, signaling.type = 'Secreted Signaling')
df=df$signaling.contribution

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, signaling.type = 'Secreted Signaling',signaling = c("OSM","IL4","IL6","IL1","TGFb","TNF","SPP1"),
color.use = c("grey50","#7A3DB8"))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, signaling.type = 'Secreted Signaling',signaling = c("OSM","IL4","IL6","IL1","TGFb","TNF","SPP1"),
color.use = c("grey50","#7A3DB8"))
gg1 + gg2

png("20250730_ranked_pathways.png", width = 5.5, height = 3, units = 'in', res = 600)
gg1 + gg2
dev.off()

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, signaling.type = 'Secreted Signaling')

color_map_list <- c(unname(detailed.celltype.colors))


ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "outgoing", signaling = c("OSM","IL4","IL6","IL1","TGFb","TNF","SPP1"), color.use = detailed.celltype.colors, title = names(object.list)[1], width = 5, height = 3)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "outgoing", signaling = c("OSM","IL4","IL6","IL1","TGFb","TNF","SPP1"),color.use = detailed.celltype.colors, title = names(object.list)[2], width = 5, height = 3)
png("20250730_WT_JAK1_cellchat_outgoingsignaling.png", width = 10, height = 3, units = 'in', res = 600)
print(draw(ht1 + ht2, ht_gap = unit(0.5, "cm")))
dev.off()

ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "outgoing", signaling = c("OSM","IL4","IL6"),color.use = detailed.celltype.colors, title = names(object.list)[2], width = 6, height = 1, color.heatmap = "OrRd")
png("20250730_JAK1_cellchat_outgoingsignaling.png", width = 4, height = 4, units = 'in', res = 600)
print(draw(ht2, ht_gap = unit(0.5, "cm")))
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming", signaling = c("OSM","IL4","IL6","IL1","TGFb","TNF","SPP1"), color.use = detailed.celltype.colors, title = names(object.list)[1], width = 5, height = 3)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming", signaling = c("OSM","IL4","IL6","IL1","TGFb","TNF","SPP1"), color.use = detailed.celltype.colors, title = names(object.list)[2], width = 5, height = 3)
png("20250730_WT_JAK1_cellchat_incomingsignaling.png", width = 10, height = 3, units = 'in', res = 600)
print(draw(ht1 + ht2, ht_gap = unit(0.5, "cm")))
dev.off()

ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming", signaling = c("OSM","IL4","IL6"), color.use = detailed.celltype.colors, title = names(object.list)[2], width = 6, height = 1, color.heatmap = "OrRd")
png("20250730_JAK1_cellchat_incomingsignaling.png", width = 4, height = 4, units = 'in', res = 600)
print(draw(ht2, ht_gap = unit(0.5, "cm")))
dev.off()


######### PDF ###########
############################
# WT outgoing
############################

ht_wt_outgoing <- netAnalysis_signalingRole_heatmap(
  object.list[[1]],
  pattern = "outgoing",
  signaling = c("OSM","IL4","IL6","IL1","TGFb","TNF","SPP1"),
  color.use = detailed.celltype.colors,
  title = names(object.list)[1],
  width = 5,
  height = 3
)

pdf(
  "20250730_WT_cellchat_outgoingsignaling.pdf",
  width = 5,
  height = 3
)

draw(ht_wt_outgoing)

dev.off()


############################
# JAK1 outgoing
############################

ht_jak_outgoing <- netAnalysis_signalingRole_heatmap(
  object.list[[2]],
  pattern = "outgoing",
  signaling = c("OSM","IL4","IL6","IL1","TGFb","TNF","SPP1"),
  color.use = detailed.celltype.colors,
  title = names(object.list)[2],
  width = 5,
  height = 3
)

pdf(
  "20250730_JAK1_cellchat_outgoingsignaling.pdf",
  width = 5,
  height = 3
)

draw(ht_jak_outgoing)

dev.off()


############################
# JAK1 focused outgoing
############################

ht_jak_outgoing_focus <- netAnalysis_signalingRole_heatmap(
  object.list[[2]],
  pattern = "outgoing",
  signaling = c("OSM","IL4","IL6"),
  color.use = detailed.celltype.colors,
  title = names(object.list)[2],
  width = 6,
  height = 1,
  color.heatmap = "OrRd"
)

pdf(
  "20250730_JAK1_cellchat_outgoingsignaling_focus.pdf",
  width = 4,
  height = 4
)

draw(ht_jak_outgoing_focus)

dev.off()


ht_jak_incoming_focus <- netAnalysis_signalingRole_heatmap(
  object.list[[2]],
  pattern = "incoming",
  signaling = c("OSM","IL4","IL6"),
  color.use = detailed.celltype.colors,
  title = names(object.list)[2],
  width = 6,
  height = 1,
  color.heatmap = "OrRd"
)

pdf(
  "20250730_JAK1_cellchat_incomingsignaling_focus.pdf",
  width = 4,
  height = 4
)

draw(ht_jak_incoming_focus)

dev.off()

############################
# WT incoming
############################

ht_wt_incoming <- netAnalysis_signalingRole_heatmap(
  object.list[[1]],
  pattern = "incoming",
  signaling = c("OSM","IL4","IL6","IL1","TGFb","TNF","SPP1"),
  color.use = detailed.celltype.colors,
  title = names(object.list)[1],
  width = 5,
  height = 3
)

pdf(
  "20250730_WT_cellchat_incomingsignaling.pdf",
  width = 5,
  height = 3
)

draw(ht_wt_incoming)

dev.off()


############################
# JAK1 incoming
############################

ht_jak_incoming <- netAnalysis_signalingRole_heatmap(
  object.list[[2]],
  pattern = "incoming",
  signaling = c("OSM","IL4","IL6","IL1","TGFb","TNF","SPP1"),
  color.use = detailed.celltype.colors,
  title = names(object.list)[2],
  width = 5,
  height = 3
)

pdf(
  "20250730_JAK1_cellchat_incomingsignaling.pdf",
  width = 5,
  height = 3
)

draw(ht_jak_incoming)

dev.off()

#########################

netVisual_bubble(cellchat, sources.use = 1, targets.use = c(5,15),  comparison = c(2,1), angle.x = 45, signaling = c("OSM","IL6","IL4"))

pdf("20260109_JAK1_LR_bubble_OSM_IL6.pdf", width = 5, height = 2)
print(
  netVisual_bubble(
    jak,
    sources.use = c("Basophil","Neutrophil","MAC","Mast","ILC2"),
    targets.use = c("KC","FIB","Pericyte"),
    signaling = c("OSM", "IL6"),
    remove.isolate = TRUE,
    sort.by.source = TRUE,
    sort.by.target = TRUE,
    angle.x = 45
  )
)
dev.off()


