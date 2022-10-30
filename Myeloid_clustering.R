library(Seurat)
library(ggplot2)
library(dplyr)
library(ggsci)

clinical <- read.table("E:/R/UCSC-CRC/singlecell/CRC.Leukocyte.10x.Metadata.txt",sep="\t",header = T)
CRC_M<-readRDS("E:/R/UCSC-CRC/singlecell/CRC_M.rds")
CRC_M<-CRC_M[,substr(colnames(CRC_M),3,3)=="T"]
CRC_M<-log2(CRC_M+1)

rownames(clinical)<-clinical$CellName
clinical<-clinical[colnames(CRC_M),]
all(rownames(clinical)==colnames(CRC_M))

##数据清理

table(clinical$Sub_ClusterID)
clinical<-clinical[!clinical$Sub_ClusterID %in% c("hM07","hM09"),]
table(clinical$Sample)
clinical<-clinical[!clinical$Sample %in% "P0305",]

CRC_M<-CRC_M[,rownames(clinical)]


##Seurat

scRNA = CreateSeuratObject(counts=CRC_M,meta.data = clinical,project = "Allcells",min.cells = 0, min.features = 0)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000) 
scRNA <- ScaleData(scRNA, features = (rownames(scRNA)))
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 
scRNA <- JackStraw(scRNA,num.replicate = 100)
scRNA <- ScoreJackStraw(scRNA,dims = 1:20)

scRNA1 <- FindNeighbors(scRNA, dims = 1:15)
scRNA1 <- FindClusters(scRNA1, resolution = 1)
scRNA1 = RunTSNE(scRNA1, dims =1:15)

#注释
scRNA2<-scRNA1
cell.embeddings<-cbind(scRNA2@meta.data[["Global_tSNE_1"]],scRNA2@meta.data[["Global_tSNE_2"]])
rownames(cell.embeddings)<-scRNA1@meta.data[["CellName"]]
colnames(cell.embeddings)<-c("tSNE_1","tSNE_2")
scRNA2@reductions[["tsne"]]@cell.embeddings=cell.embeddings
Sub_Cluster<-gsub("hM01","Mast",scRNA2@meta.data[["Sub_ClusterID"]])
Sub_Cluster[which(Sub_Cluster %in% c("hM02","hM03","hM04"))] <-"DC"
Sub_Cluster[which(Sub_Cluster %in% c("hM05","hM06"))] <-"Mono"
Sub_Cluster[which(Sub_Cluster %in% c("hM08","hM10","hM11","hM12","hM13"))] <-"Macro"
scRNA2@meta.data[["seurat_clusters"]]=Sub_Cluster

saveRDS(scRNA2,"scRNA.rds")
saveRDS(clinical,"clinical.rds")


#tSNE图
DimPlot(scRNA2, reduction = "tsne", group.by = "Sub_ClusterID", 
        pt.size=0.5, label = TRUE,repel = TRUE,cols = c("#6587AB","#2D1040","#923003","#B16D1A","#466333","#0E3226","#D9B779","#5C554B","#A65A49","#F2856D","#8C3D35"))+
  ggtitle("Clusters of Myeloid Cells")

#比例图
cell_pro <- as.data.frame(table(scRNA2@meta.data[["Sub_Cluster"]],scRNA2@meta.data[["seurat_clusters"]]))
colnames(cell_pro) <- c('Subcelltype','Cell','proportion')
ggplot(cell_pro,aes(x=Cell,y=proportion,fill=Subcelltype))+
  geom_bar(stat = 'identity',position = 'fill')+
  scale_fill_manual(values = c("#6587AB","#2D1040","#923003","#B16D1A","#466333","#0E3226","#D9B779","#5C554B","#A65A49","#F2856D","#8C3D35"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text = element_text(face='bold',size = 8,angle = 45,family = 'sans'))

#小提琴图
VlnPlot(scRNA2,
        features = c("HLA-DMA","HOXC8","HAMP","APOE","SIRPB2"),
        stack = TRUE,
        group.by = 'Sub_Cluster',
        cols = c("#6587AB","#2D1040","#923003","#B16D1A","#466333","#0E3226","#D9B779","#5C554B","#A65A49","#F2856D","#8C3D35"),
        split.by = 'Sub_Cluster',
        split.plot = F,
        same.y.lims = F)+
  theme(axis.text.x = element_text(hjust = 0,size = 8,family = "sans",face = 'bold'),
        axis.text.y = element_text(hjust = 1,vjust = 0,size = 9,family = "sans",face = 'bold'))+
  theme(legend.position = 'none')+
  xlab('')+ylab('')

