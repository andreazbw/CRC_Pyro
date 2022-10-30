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

scRNA<-readRDS("scRNA2.rds")
clinical<-readRDS("clinical2.rds")
CRC_M<- scRNA@assays[["RNA"]]@data %>% as.data.frame()

clinical<-clinical[clinical$Sub_ClusterID %in% "hM12",]
clinical$Sub_ClusterID<-as.factor(clinical$Sub_ClusterID)
CRC_M<-CRC_M[,rownames(clinical)]
all(rownames(clinical)==colnames(CRC_M))

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
  stat_compare_means(hide.ns = TRUE,comparisons = my_comparisons,label.y = c(2.5, 3, 3.5),method = "t.test") +
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

clinical2<-test_suv

clinical2$Global_ID<-factor(clinical2$group,levels=c("hM12_H","hM12_M","hM12_L"))
gsva.es_avg<-data.frame(gsva.es[,rownames(clinical2)])

for (i in levels(clinical2$Global_ID)) {
  gsva.es_sub <- gsva.es[,clinical2$Global_ID==i]
  gsva.es_sub <- rowMeans(gsva.es_sub)
  gsva.es_avg <- cbind(gsva.es_avg,gsva.es_sub)
}
gsva.es_avg<-gsva.es_avg[,550:552]
colnames(gsva.es_avg) <- levels(clinical2$Global_ID)
gsva.es_avg <- gsva.es_avg %>% t() %>% scale() %>% t()
summary(as.numeric(gsva.es_avg))

#热图
heatmap<-gsva.es_avg
pheatmap(heatmap,border_color = 'white',
         annotation_names_row = FALSE,
         show_rownames = T,show_colnames = T,
         angle_col = 0,
         cellwidth = 40,
         cellheight = 15,
         cluster_cols = F,cluster_rows = T,
         color = colorRampPalette(c("#1E63A7","white","#AD152A"))(100))

saveRDS(datExpr,"hm12_expr.rds")
saveRDS(test_suv,"hm12_suv.rds")

#比例图
cell_pro <- clinical2[,c('Sample','group')]
cell_pro$Sample <- mapvalues(cell_pro$Sample,from = unique(cell_pro$Sample),
                             to = c('Stage II','Stage III','Stage III','Stage II','Stage III',
                                    'Stage I','Stage III','Stage II','Stage II'))
cell_pro <- as.data.frame(table(cell_pro$Sample,cell_pro$group))
colnames(cell_pro) <- c('stage','celltype','proportion')
cell_pro <- dcast(cell_pro,stage~celltype,)
rownames(cell_pro) <- cell_pro[,1]
cell_pro <- cell_pro[,-1]
cell_pro <- cell_pro/rowSums(cell_pro)
cell_pro <- melt(cell_pro)
cell_pro$stage <- rep(c('Stage I','Stage II','Stage III'),times=3)
colnames(cell_pro) <- c('celltype','proportion','stage')
cell_pro$celltype <- factor(cell_pro$celltype,levels = c('hM12_H','hM12_M','hM12_L'))

ggplot(cell_pro,aes(x=celltype,y=proportion,fill=stage))+
  geom_bar(stat = 'identity',position = 'fill')+
  scale_fill_manual(values=c("#E7C0B1","#7E81A5","#44302E"))+
  theme_classic()+
  theme(axis.text = element_text(face='bold',size = 10,angle = 0,family = 'sans',colour = 'black'),
        axis.text.x =element_text(vjust = 0.5),
        plot.title = element_text(hjust = 0.5,size = 20))+
  guides(fill = guide_legend(title = NULL))+
  ggtitle('Cancer stage proportion of TAMs in each pyroptosis state')+
  xlab('Celltypes')+ylab('Proportion')

