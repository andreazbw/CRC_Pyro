library(GSVA)
library(dplyr)
library(tidyverse)
library(survival)
library(survminer)
library(pheatmap)
library(estimate)
library(ggthemes)
library(DESeq2)
library(patchwork)

tumor_tpm<-readRDS("E:/R/UCSC-CRC/download/tumor_tpm.rds")
tumor_count<-readRDS("E:/R/UCSC-CRC/download/tumor_count.rds")
survival_os<-readRDS("E:/R/UCSC-CRC/download/survival_os.rds")

tumor_tpm<-tumor_tpm[,intersect(colnames(tumor_tpm),rownames(survival_os))]
all(colnames(tumor_tpm)==rownames(survival_os))
tumor_count<-tumor_count[,intersect(colnames(tumor_count),rownames(survival_os))]
all(colnames(tumor_count)==rownames(survival_os))


##ssgsea

tisidb <- read.csv("29Reports.txt",header = F,sep = "\t",)
tisidb <- tisidb %>%column_to_rownames("V1")%>%t()
a <- tisidb
a <- a[1:nrow(a),]
set <- colnames(a)
l <- list()
#i <- "Activated CD8 T cell"
for (i in set) {
  x <-  as.character(a[,i])
  x <- x[nchar(x)!=0]
  x <-  as.character(x)
  l[[i]] <-x
}

tpm_immune<- gsva(tumor_tpm, l,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)


##层次聚类

df <- scale(t(tpm_immune))
result <- dist(df, method = "euclidean")
result_hc <- hclust(d = result, method = "ward.D")
cluster<-cutree(result_hc,k=2)
table(cluster)


##两组生存分析

datSurv <- data.frame(survival_os[,1:3],cluster)
datSurv$cluster[which(datSurv$cluster %in% c("1"))] <-"Immunity_L"
datSurv$cluster[which(datSurv$cluster %in% c("2"))] <-"Immunity_H"
group<-datSurv$cluster
fit_km <- surv_fit(Surv(OS.time,OS) ~ group, data = datSurv)
ggsurvplot(fit_km,
           legend.title='Immunity Groups',
           conf.int = F,
           pval = T,
           pval.method = T,
           combine = T,
           pval.size = 8,
           risk.table = T,
           risk.table.title = "",
           palette = c("#F26E50","#3C92A6"),
           legend.labs = c("Immunity_H","Immunity_L"),
           tables.height = 0.2,
           xlab = "Time In Days",
           ylab = "Overall Survival",
           surv.median.line = "hv",
           ncensor.plot = F)

saveRDS(datSurv,"immune_group.rds")


##Estimate评分

estimate <- function(tumor_tpm,pro){
  input.f=paste0(pro,'_estimate_input.txt')
  output.f=paste0(pro,'_estimate_gene.gct')
  output.ds=paste0(pro,'_estimate_score.gct')
  write.table(tumor_tpm,file = input.f,sep = '\t',quote = F)
  filterCommonGenes(input.f=input.f,
                    output.f=output.f ,
                    id="GeneSymbol")
  estimateScore(input.ds = output.f,
                output.ds=output.ds,
                platform="affymetrix")
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}
pro='CRC'
scores=estimate(tumor_tpm,pro)
a<-gsub("[.]","-",rownames(scores))
rownames(scores)<-a

write.csv(scores,"scores.csv")


##画图

#热图
all(rownames(datSurv)==rownames(scores))
anno<-cbind(datSurv,scores)
anno<-arrange(anno,cluster)
anno<-as.data.frame(anno[,2:8])
anno<-anno[,-2]
anno<-anno[,c(2,1,3:6)]
anno$OS<-as.factor(anno$OS)

heatmap<-tpm_immune[,rownames(anno)]
all(colnames(heatmap)==rownames(anno))
bk <- c(seq(-5,-0.01,by=0.01),seq(0,5,by=0.01))

annotation_color <- list(ImmuneScore=c('#F2F2F2','#EF6024'),
                         ESTIMATEScore=c('#F2F2F2','#F0941F'),
                         StromalScore=c('#F2F2F2','#90A19D'),
                         TumorPurity=c('#F2F2F2','#196774'),
                         OS=c("0"='#BFBFBF',"1"='#363432'),
                         cluster=c(Immunity_H="#F26E50",Immunity_L="#3C92A6"))
pheatmap(heatmap, 
         color=c(colorRampPalette(colors = c("#3C92A6","white"))(length(bk)/2),colorRampPalette(colors = c("white","#F26E50"))(length(bk)/2)),
         breaks = bk,
         scale = "row",
         cluster_rows=T,
         cluster_cols = F,
         cellwidth = 1,
         cellheight = 20,
         fontsize = 8,
         annotation_col = anno,
         annotation_colors = annotation_color,
         show_colnames = FALSE,
         show_rownames = TRUE
)

#四种评分箱线图
boxplot<-anno
boxplot$id<-rownames(boxplot)
boxplot<-boxplot %>% gather(key="score",value="expression",3) %>% 
  dplyr::select(id,score,expression,cluster) 

ggplot(aes(x = cluster, y = expression,fill = cluster),data = boxplot) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "StromalScore")+
  scale_x_discrete(name = "") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face =  "bold",family = 'sans',hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11,family = 'sans',face = 'bold'),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
  scale_fill_manual(values = c("#F26E50","#3C92A6"))+
  stat_compare_means(method = 'wilcox.test',label = 'p.signif',hide.ns = TRUE,label.x = 1.5)

#HLA基因箱线图
boxplot<-as.data.frame(tumor_tpm)
boxplot<-t(boxplot[c("HLA-G","HLA-F","HLA-E","HLA-DRB5","HLA-DRB1","HLA-DRA",
                   "HLA-DQB2","HLA-DQB1","HLA-DQA2","HLA-DQA1","HLA-DPB1",
                   "HLA-DPA1","HLA-DOB","HLA-DOA","HLA-DMB","HLA-DMA",
                   "HLA-C","HLA-B","HLA-A"),])
boxplot<-boxplot[rownames(anno),]
all(rownames(boxplot)==rownames(anno))
boxplot<-cbind(boxplot,anno)
boxplot$id<-rownames(boxplot)
boxplot<-boxplot %>% gather(key="gene",value="expression",1:19) %>% 
  dplyr::select(id,gene,expression,cluster) 

ggplot(aes(x = gene, y = expression,fill = cluster),data = boxplot) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Gene expression")+
  scale_x_discrete(name = "") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face =  "bold",family = 'sans',hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.y=element_text(size = 11,family = 'sans',face = 'bold'),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
  scale_fill_manual(values = c("#F26E50","#3C92A6"))+
  stat_compare_means(method = 'wilcox.test',label = 'p.signif',hide.ns = TRUE,label.x = 1.5)+
  coord_flip()

#免疫检查点基因箱线图
boxplot<-as.data.frame(tumor_tpm)
boxplot<-t(boxplot[c('PDCD1','CD274','PDCD1LG2','CTLA4','CD80','CD86'),])
boxplot<-boxplot[rownames(anno),]
all(rownames(boxplot)==rownames(anno))
boxplot<-cbind(boxplot,anno)
boxplot$id<-rownames(boxplot)
boxplot<-boxplot %>% gather(key="gene",value="expression",6) %>% 
  dplyr::select(id,gene,expression,cluster) 

ggplot(aes(x = cluster, y = expression,fill = cluster),data = boxplot) +
  geom_violin(aes(fill = cluster),alpha=0.7)+
  geom_boxplot(alpha=0) +
  scale_y_continuous(name = "CD86 expression")+
  scale_x_discrete(name = "") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face =  "bold",family = 'sans',hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.y=element_text(size = 11,family = 'sans',face = 'bold'),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
  scale_fill_manual(values = c("#F26E50","#3C92A6"))+
  stat_compare_means(method = 'wilcox.test',label = 'p.signif',hide.ns = TRUE,label.x = 1.5)


##差异分析

des_count<-tumor_count[,rownames(anno)]
all(colnames(des_count)==rownames(anno))

group_list<-anno$cluster
group_list[which(group_list %in% c("Immunity_L"))] <-"Low"
group_list[which(group_list %in% c("Immunity_H"))] <-"High"
group_list <- factor(group_list,levels = c('High','Low'))
coldata <- data.frame(row.names=colnames(des_count),tissue=group_list)
dds <- DESeqDataSetFromMatrix(countData = des_count, 
                              colData = coldata, 
                              design = ~tissue)
dds$tissue<- relevel(dds$tissue, ref = "Low")
dds <- DESeq(dds)

DESeq2_result <- as.data.frame(results(dds))
DESeq2<-subset(DESeq2_result,DESeq2_result$padj< 0.05 & abs(DESeq2_result$log2FoldChange)>1)

saveRDS(DESeq2,"DESeq2.rds")

