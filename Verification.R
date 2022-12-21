
clinical <- read.table("CRC.Leukocyte.Smart-seq2.Metadata.txt",sep="\t",header = T)
CRC_tpm <- read.table("CRC.Leukocyte.Smart-seq2.TPM.txt",sep=" ",header = T)
all(clinical$CellName==colnames(CRC_tpm))

clinical<-clinical[clinical$Tissue %in% "T",]
clinical<-clinical[clinical$Sub_ClusterID %in% "hM12",]
rownames(clinical)<-clinical$CellName
clinical$Sub_ClusterID<-as.factor(clinical$Sub_ClusterID)

CRC_M<-CRC_tpm[,rownames(clinical)]
all(rownames(clinical)==colnames(CRC_M))

CRC_M<-log2(CRC_M+1)


library(dplyr)
library(tidyr)
library(ggpubr)
library(ggthemes)
library(GSVA)
library(reshape2)
library(pheatmap)
library(Seurat)
library(psych)
library(plyr)

test_suv<-cbind(clinical,t(CRC_M)[,c("APOE","NLRP3","CASP1","GSDMD","IL1B","IL18")])
test_suv$group1<-ifelse(test_suv$NLRP3>mean(test_suv$NLRP3),"High","Low")
test_suv$group2<-ifelse(test_suv$CASP1>mean(test_suv$CASP1),"High","Low")
test_suv$group3<-ifelse(test_suv$GSDMD>mean(test_suv$GSDMD),"High","Low")
test_suv<-unite(test_suv,"group",group1,group2,group3)

test_suv$group[which(test_suv$group %in% c("High_High_High"))] <-"hM12_H"
test_suv$group[which(test_suv$group %in% c("Low_High_High","High_Low_High","High_High_Low","Low_Low_High"))] <-"hM12_M"
test_suv$group[which(test_suv$group %in% c("High_Low_Low","Low_High_Low","Low_Low_Low"))] <-"hM12_L"

boxplot<-test_suv
boxplot$id<-rownames(boxplot)
boxplot<-boxplot[,24:31]
boxplot<-boxplot %>% gather(key="gene",value="expression",1:6) %>% 
  dplyr::select(id,gene,expression,group)
boxplot$group<-factor(boxplot$group,levels=c("hM12_H","hM12_M","hM12_L"))
my_comparisons <- list( c("hM12_H", "hM12_M"), c("hM12_H", "hM12_L"), c("hM12_M", "hM12_L") )
ggboxplot(boxplot, x = "group", y = "expression",
          color = "black",fill="group",facet.by = "gene",palette = "nejm",
          xlab = "", ylab = "The expression of Genes(log2 TPM)") +
  scale_fill_manual(values=c("#BF6734","#F2D0A7","#5D768C"))+
  stat_compare_means(hide.ns = TRUE,comparisons = my_comparisons,label.y = c(2.5, 3, 3.5),method = "wilcox.test") +
  stat_compare_means(label.y = 5,method = "anova") +
  theme_few()+
  theme(plot.title = element_text(size = 10, face =  "bold",family = 'sans',hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,angle = 0,family = 'sans',vjust = 1,hjust = 1,face = 'bold'),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())


##GSVA

geneset <- read.csv('Gsva.csv',row.names = 1,header = FALSE,na.strings = '') %>% t()
datExpr <- as.matrix(CRC_M)

GSVA_gs <- melt(geneset,measure.vars = colnames(geneset),variable.name = 'Description',value.name = 'gene') %>% na.omit()
GSVA_gs <- split(GSVA_gs$gene, GSVA_gs$Var2)

gsva.es <- gsva(datExpr, GSVA_gs,verbose=TRUE,method="gsva",kcdf="Gaussian",min.sz=1,max.sz=1000)

#取平均
clinical2<-test_suv

clinical2$Global_ID<-factor(clinical2$group,levels=c("hM12_H","hM12_M","hM12_L"))
gsva.es_avg<-data.frame(gsva.es[,rownames(clinical2)])

for (i in levels(clinical2$Global_ID)) {
  gsva.es_sub <- gsva.es[,clinical2$Global_ID==i]
  gsva.es_sub <- rowMeans(gsva.es_sub)
  gsva.es_avg <- cbind(gsva.es_avg,gsva.es_sub)
}
gsva.es_avg<-gsva.es_avg[,266:268]
colnames(gsva.es_avg) <- levels(clinical2$Global_ID)
gsva.es_avg <- gsva.es_avg %>% t() %>% scale() %>% t()
summary(as.numeric(gsva.es_avg))

#画热图
heatmap<-gsva.es_avg
pheatmap(heatmap,border_color = 'white',
         annotation_names_row = FALSE,
         show_rownames = T,show_colnames = T,
         angle_col = 0,
         cellwidth = 40,
         cellheight = 15,
         cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c("#1E63A7","white","#AD152A"))(100))


