library(Seurat)
library(ggplot2)
library(dplyr)
library(ggsci)

clinical <- read.table("E:/R/UCSC-CRC/singlecell/CRC.Leukocyte.10x.Metadata.txt",sep="\t",header = T)
CRC.T<-readRDS("E:/R/UCSC-CRC/singlecell/CRC.T.rds")
CRC.T<-log2(CRC.T+1)

rownames(clinical)<-clinical$CellName
clinical<-clinical[colnames(CRC.T),]
all(rownames(clinical)==colnames(CRC.T))

##数据清理

table(clinical$Sub_ClusterID)
clinical<-clinical[!clinical$Sub_ClusterID %in% c("hI04","hM07","hM09"),]
table(clinical$Sample)
clinical<-clinical[!clinical$Sample %in% "P0305",]

CRC.T<-CRC.T[,rownames(clinical)]


##Seurat

scRNA = CreateSeuratObject(counts=CRC.T,meta.data = clinical,project = "Allcells",min.cells = 0, min.features = 0)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000) 
scRNA <- ScaleData(scRNA, features = (rownames(scRNA)))
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA)) 
scRNA <- JackStraw(scRNA,num.replicate = 100)
scRNA <- ScoreJackStraw(scRNA,dims = 1:20)

scRNA1 <- FindNeighbors(scRNA, dims = 1:20)
scRNA1 <- FindClusters(scRNA1, resolution = 1.5)
scRNA1 = RunTSNE(scRNA1, dims =1:20)

#注释
scRNA2<-scRNA1
cell.embeddings<-cbind(scRNA2@meta.data[["Global_tSNE_1"]],scRNA2@meta.data[["Global_tSNE_2"]])
rownames(cell.embeddings)<-scRNA1@meta.data[["CellName"]]
colnames(cell.embeddings)<-c("tSNE_1","tSNE_2")
scRNA2@reductions[["tsne"]]@cell.embeddings=cell.embeddings
scRNA2@meta.data[["seurat_clusters"]]=scRNA2@meta.data[["Global_Cluster"]]
saveRDS(scRNA2,"scRNA.rds")
saveRDS(clinical,"clinical.rds")

#tSNE图
DimPlot(scRNA2, reduction = "tsne", group.by = "seurat_clusters", 
        pt.size=0.5, label = TRUE,repel = TRUE,cols = c("#F44336","#FCE5CD","#FFE599","#93C47D","#134F5C"))+
  ggtitle("Clusters of all Immune Cells")

#比例图
cell_pro <- as.data.frame(table(scRNA2@meta.data[["seurat_clusters"]],scRNA2@meta.data[["Sample"]]))
colnames(cell_pro) <- c('celltype','ID','proportion')
ggplot(cell_pro,aes(x=ID,y=proportion,fill=celltype))+
  geom_bar(stat = 'identity',position = 'fill')+
  scale_fill_manual(values = c("#F44336","#FCE5CD","#FFE599","#93C47D","#134F5C"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text = element_text(face='bold',size = 8,angle = 45,family = 'sans'))

#小提琴图
VlnPlot(scRNA2,
        features = c("HLA-DMA","HOXC8","HAMP","APOE","SIRPB2"),
        stack = TRUE,
        group.by = 'seurat_clusters',
        cols = c("#F44336","#FCE5CD","#FFE599","#93C47D","#134F5C"),
        split.by = 'seurat_clusters',
        split.plot = F,
        same.y.lims = F)+
  theme(axis.text.x = element_text(hjust = 0,size = 8,family = "sans",face = 'bold'),
        axis.text.y = element_text(hjust = 1,vjust = 0,size = 9,family = "sans",face = 'bold'))+
  theme(legend.position = 'none')+
  xlab('')+ylab('')
