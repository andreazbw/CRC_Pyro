library(monocle)
library(survival)
library(survminer)
library(DESeq2)
library(dplyr)
library(tidyverse)
library(ggthemes)
library(pheatmap)

tumor_count<-readRDS("E:/R/UCSC-CRC/download/tumor_count.rds")
tumor_tpm<-readRDS("E:/R/UCSC-CRC/download/tumor_tpm.rds")
survival_os<-readRDS("E:/R/UCSC-CRC/download/survival_os.rds")
pyro_all<- read.csv('genecards_a1.csv',header = F)

tumor_count<-tumor_count[,intersect(colnames(tumor_count),rownames(survival_os))]
tumor_tpm<-tumor_tpm[,intersect(colnames(tumor_tpm),rownames(survival_os))]
all(colnames(tumor_count)==rownames(survival_os))
all(colnames(tumor_tpm)==rownames(survival_os))
tumor_pyro<-tumor_count[intersect(pyro_all$V1,rownames(tumor_count)),]


##monocle聚类

ct<-as.matrix(tumor_pyro)
gene_ann <- data.frame(
  gene_short_name = row.names(ct),
  row.names = row.names(ct)
)

pd <- new("AnnotatedDataFrame",
          data=survival_os)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)

cds <- newCellDataSet(
  ct,
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0)

cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)

plot_pc_variance_explained(cds, 
                           return_all = F,
                           max_components = 20)
cds <- reduceDimension(cds, 
                       max_components = 2, 
                       num_dim = 6,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 3)

plot_cell_clusters(cds,cell_size = 2.5)+
  scale_y_continuous(name = "tSNE_2")+
  scale_x_discrete(name = "tSNE_1") +
  scale_color_manual(values=c("#043F8C","#F2AC29"))

#画图
datSurv <- data.frame(survival_os[,1:3],cds@phenoData@data[["Cluster"]])
heatmap_anno<-arrange(datSurv,cds.phenoData.data...Cluster...)
heatmap_anno<-heatmap_anno[,c(2,4)]
heatmap_anno$OS<-as.factor(heatmap_anno$OS)
colnames(heatmap_anno)=c("OS","cluster")
heatmap<-tumor_tpm[,rownames(heatmap_anno)]
heatmap<-heatmap[pyro_all$V1,]
all(colnames(heatmap)==rownames(heatmap_anno))
bk <- c(seq(-4,-0.01,by=0.01),seq(0,4,by=0.01))
annotation_color <- list(OS=c("0"='#BFBFBF',"1"='#363432'),
                         cluster=c("2"="#F2AC29","1"="#043F8C"))
pheatmap(heatmap, 
         color=c(colorRampPalette(colors = c("#043F8C","white"))(length(bk)/2),colorRampPalette(colors = c("white","#F2AC29"))(length(bk)/2)),
         breaks = bk,
         scale = "row",
         cluster_rows=T,
         cluster_cols = F,
         cellwidth = 2,
         cellheight = 7,
         fontsize = 8,
         annotation_col = heatmap_anno,
         annotation_colors = annotation_color,
         show_colnames = FALSE,
         show_rownames = TRUE
)


##生存分析

datSurv$cds.phenoData.data...Cluster...<-as.numeric(datSurv$cds.phenoData.data...Cluster...)
datSurv$cds.phenoData.data...Cluster...[which(datSurv$cds.phenoData.data...Cluster... %in% c("1"))] <-"Pyroptosis_L"
datSurv$cds.phenoData.data...Cluster...[which(datSurv$cds.phenoData.data...Cluster... %in% c("2"))] <-"Pyroptosis_H"
group<-datSurv$cds.phenoData.data...Cluster...

fit_km <- surv_fit(Surv(OS.time,OS) ~ group, data = datSurv)
ggsurvplot(fit_km,
           legend.title='Pyroptosis Groups',
           conf.int = F,
           pval = T,
           pval.method = T,
           combine = T,
           pval.size = 8,
           risk.table = T,
           risk.table.title = "",
           palette = c("#F2AC29","#043F8C"),
           legend.labs = c("Pyroptosis_H","Pyroptosis_L"),
           tables.height = 0.2,
           xlab = "Time In Days",
           ylab = "Overall Survival",
           surv.median.line = "hv",
           ncensor.plot = F)


##差异分析

des_count<-tumor_count[,rownames(anno)]
all(colnames(des_count)==rownames(anno))

group_list=anno$cds.phenoData.data...Cluster...
group_list <- factor(group_list,levels = c('high','low'))
coldata <- data.frame(row.names=colnames(des_count),tissue=group_list)
dds <- DESeqDataSetFromMatrix(countData = des_count, 
                              colData = coldata, 
                              design = ~tissue)
dds$tissue<- relevel(dds$tissue, ref = "low")
dds <- DESeq(dds)

DESeq2_result <- as.data.frame(results(dds))
DESeq2<-subset(DESeq2_result,DESeq2_result$padj< 0.05 & abs(DESeq2_result$log2FoldChange)>1)

saveRDS(DESeq2,"DESeq2.rds")
saveRDS(datSurv,"datSurv.rds")

#画图
deseq_reu <- DESeq2_result
deseq_reu <- deseq_reu[c('NLRP1','NLRP3','AIM2','GSDMC','GSDMD','DFNA5','CASP1','CASP4','CASP5','IL1B','IL6','IL18'),]

datbar <- deseq_reu[,c(2,6)]
datbar$group <- as.factor(rep(1:4,each=3))
datbar <- datbar[order(datbar$group,-datbar$log2FoldChange),]
datbar$gene <- factor(rownames(datbar),levels = rownames(datbar))

star <- data.frame(group=datbar$gene,value=datbar$log2FoldChange+0.1)
colorboard <- ifelse(datbar$padj<0.05,'red','grey')


ggplot(data=datbar,mapping=aes(x=gene,y=log2FoldChange,fill=group))+
  geom_bar(stat = 'identity')+
  geom_text(star,mapping = aes(group,value),label='*',color=colorboard,size=7)+
  scale_fill_manual(values = c('#E47112','#FFB522','#FFE759','#FFE695',rep('white',12)))+
  theme_classic()+
  ggtitle('')+
  xlab('')+ylab('log2(fold change)')+
  theme(axis.text = element_text(face='bold',size = 10,angle = 0,family = 'sans',colour = 'black'),
        axis.text.x =element_text(vjust = 0.5,angle = 0),
        plot.title = element_text(hjust = 0.5,size = 6),
        legend.position='none')
