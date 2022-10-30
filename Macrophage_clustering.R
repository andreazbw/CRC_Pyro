library(dplyr)
library(Seurat)
library(ggplot2)
library(ggsci)
library(GSVA)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(psych)
library(rstatix)
library(ggpubr)
library(corrplot)

scRNA<-readRDS("scRNA.rds")
clinical<-readRDS("clinical.rds")

clinical<-clinical[clinical$Sub_ClusterID %in% c("hM08","hM10","hM11","hM12","hM13"),]
clinical$Sub_ClusterID<-as.factor(clinical$Sub_ClusterID)
CRC_M<- scRNA@assays[["RNA"]]@data %>% as.data.frame()
CRC_M<-CRC_M[,rownames(clinical)]


##seurat

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
scRNA2@meta.data[["seurat_clusters"]]=scRNA2@meta.data[["Sub_ClusterID"]]

saveRDS(scRNA2,"scRNA2.rds")
saveRDS(clinical,"clinical2.rds")

#画图
VlnPlot(scRNA2,
        features = c('GSDMD','CASP1','NLRP3','PYCARD','IL1B','IL18','TREM2','APIP','BRD4'),
        stack = TRUE,
        group.by = 'Sub_ClusterID',
        cols = c("#D9B779","#5C554B","#A65A49","#F2856D","#8C3D35"),
        split.by = 'Sub_ClusterID',
        split.plot = F,
        same.y.lims = F)+
  theme(axis.text.x = element_text(hjust = 0,size = 8,family = "sans",face = 'bold'),
        axis.text.y = element_text(hjust = 1,vjust = 0,size = 9,family = "sans",face = 'bold'))+
  theme(legend.position = 'none')+
  xlab('')+ylab('')


##anova

anova_dat <- as.data.frame(t(CRC_M["APOE",]))
anova_dat$cluster <- clinical$Sub_ClusterID

anova_dat2<-anova_dat
anova_dat2$cluster<-factor(anova_dat2$cluster, levels = c("hM08", "hM10", "hM11","hM12","hM13"))
anova_pic<- anova_dat2 %>% tukey_hsd(APOE ~ cluster) %>% add_xy_position(x = "cluster")

#画图
ggboxplot(anova_dat2, x = "cluster", y = "APOE",fill = "cluster",alpha=0.7)+
  scale_y_continuous(name = "APOE expression")+
  scale_x_discrete(name = "") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face =  "bold",family = 'sans',hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11,family = 'sans',face = 'bold'),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
  scale_fill_manual(values = c("#D9B779","#5C554B","#A65A49","#F2856D","#8C3D35"))+  
  stat_pvalue_manual(anova_pic, hide.ns = TRUE) +
  labs(subtitle = get_test_label(anova_dat2 %>% anova_test(APOE ~ cluster), detailed = TRUE),
       caption = get_pwc_label(anova_pic))


##相关性

data_gene<-as.data.frame(t(CRC_M['APOE',]))
data_p<-as.data.frame(t(CRC_M[c('GSDMD','CASP1','NLRP3','PYCARD','IL1B','IL18','TREM2','APIP','BRD4'),]))
all(rownames(data_gene)==rownames(data_p))

c.mat.all<-c()
p.mat.all<-c()
for (i in levels(clinical$Sub_ClusterID)) {
  data_gene_sub <- data_gene[clinical$Sub_ClusterID==i,]
  data_p_sub <- data_p[clinical$Sub_ClusterID==i,]
  corr_sub<-corr.test(data_gene_sub,data_p_sub, method = "pearson",adjust= "fdr")
  c.mat.sub <-corr_sub$r
  p.mat.sub <-corr_sub$p.adj

  c.mat.all <- rbind(c.mat.all,c.mat.sub)
  p.mat.all <- rbind(p.mat.all,p.mat.sub)
}

rownames(c.mat.all) <- c('hm8','hm10','hm11','hm12','hm13')
rownames(p.mat.all) <- c('hm8','hm10','hm11','hm12','hm13')
dim(c.mat.all)

#画图
corrplot(c.mat.all,tl.col = 'black',
         cl.ratio = 0.2,tl.cex = 0.9,cl.cex = 0.8,col.lim = c(-0.5,0.5),is.corr = FALSE,method = 'color',
         tl.offset = 0.5,cl.align.text = 'l',outline="white",
         col = c(colorRampPalette(colors = c("#1E63A7","white","#AD152A"))(100)),
         p.mat  = p.mat.all,insig = 'label_sig',
         sig.level = c(.001, .01, .05), pch.cex = 1.2,  pch.col = 'black')
