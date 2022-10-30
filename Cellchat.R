library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(tidyr)
library(do)

hm12_expr<-readRDS("hm12_expr.rds")
hm12_suv<-readRDS("hm12_suv.rds")
all(rownames(hm12_suv)==colnames(hm12_expr))
hm12_suv$Global_Cluster<-hm12_suv$group
hm12_suv<-hm12_suv[,c(15,17)]

scRNA<-readRDS("E:/R/UCSC-CRC/singlecell/T-cell/scRNA.rds")
clinical<-readRDS("E:/R/UCSC-CRC/singlecell/T-cell/clinical.rds")
CRC_T<- scRNA@assays[["RNA"]]@data %>% as.data.frame()
clinical<-clinical[clinical$Global_Cluster %in% "CD8 T cell",]
clinical<-clinical[,c(15,17)]
clinical$Global_Cluster<-clinical$Sub_ClusterID
CRC_T<-CRC_T[,rownames(clinical)]


##Cellchat

data.input<-cbind(hm12_expr,CRC_T) %>% as.matrix()
identity<-rbind(hm12_suv,clinical)
all(rownames(identity)==colnames(data.input))

identity$Global_Cluster<-Replace(data=identity$Global_Cluster,pattern = c("hT12:hT12_TN","hT13:hT13_TCM","hT14:hT14_TEMRA/TEFF","hT15:hT15_TEM","hT16:hT16_TRM","hT17:hT17_IEl","hT18:hT18_TEX"))
identity$Global_Cluster<-factor(identity$Global_Cluster,levels=c("hM12_H","hM12_M","hM12_L","hT12_TN","hT13_TCM","hT14_TEMRA/TEFF","hT15_TEM","hT16_TRM","hT17_IEl","hT18_TEX"))

cellchat <- createCellChat(object = data.input)
cellchat <- addMeta(cellchat, meta = identity)
cellchat <- setIdent(cellchat, ident.use = "Global_Cluster")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human 
colnames(CellChatDB$interaction) 
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
showDatabaseCategory(CellChatDB)
unique(CellChatDB$interaction$annotation)
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat,raw.use=FALSE,population.size=TRUE)
df.net<-subsetCommunication(cellchat)
write.csv(df.net,"net_lr.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "net_pathway.csv")

#画图
cellchat <- aggregateNet(cellchat)
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Number of interactions",
                 color.use=c("#BF6734","#F2D0A7","#5D768C","#5F94D9","#7FA644","#D95E32","#724399","#D9A036","#DAA2F2","#8A6B4C"))
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T,
                 label.edge= F, title.name = "Interaction weights/strength",
                 color.use=c("#BF6734","#F2D0A7","#5D768C","#5F94D9","#7FA644","#D95E32","#724399","#D9A036","#DAA2F2","#8A6B4C"))

mat <- as.data.frame(cellchat@netP[["pathways"]],)
mat <- as.matrix(mat[mat$`cellchat@netP[["pathways"]]` %in% c("JAM","TNF","ICOS","HGF","THBS","SPP1"),])
pathways.show<-mat[1,]
netVisual_aggregate(cellchat,layout = 'hierarchy',vertex.receiver = c(1,2,3),signaling = pathways.show,
                    color.use=c("#BF6734","#F2D0A7","#5D768C","#5F94D9","#7FA644","#D95E32","#724399","#D9A036","#DAA2F2","#8A6B4C"))
netAnalysis_contribution(cellchat, signaling = pathways.show,vertex.receiver = c(1,2,3))
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = pairLR[1,],sources.use = c("hM12_H","hM12_M","hM12_L") ,layout = 'circle',
                     color.use=c("#BF6734","#F2D0A7","#5D768C","#5F94D9","#7FA644","#D95E32","#724399","#D9A036","#DAA2F2","#8A6B4C"))

