library(GSVA)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(psych)
library(ggplot2)
library(rstatix)
library(ggpubr)

scRNA<-readRDS("scRNA2.rds")
clinical<-readRDS("clinical2.rds")


##Gsva

geneset <- read.csv('Gsva.csv',row.names = 1,header = FALSE,na.strings = '') %>% t()
datExpr <- scRNA@assays[["RNA"]]@data %>% as.matrix()

GSVA_gs <- melt(geneset,measure.vars = colnames(geneset),variable.name = 'Description',value.name = 'gene') %>% na.omit()
GSVA_gs <- split(GSVA_gs$gene, GSVA_gs$Var2)

gsva.es <- gsva(datExpr, GSVA_gs,verbose=TRUE,method="gsva",kcdf="Gaussian",min.sz=1,max.sz=1000)

clinical$Global_ID<-as.factor(clinical$Sub_ClusterID)
gsva.es_avg<-data.frame(gsva.es[,rownames(clinical)])

for (i in levels(clinical$Global_ID)) {
  gsva.es_sub <- gsva.es[,clinical$Global_ID==i]
  gsva.es_sub <- rowMeans(gsva.es_sub)
  gsva.es_avg <- cbind(gsva.es_avg,gsva.es_sub)
}
gsva.es_avg<-gsva.es_avg[,1784:1788]
colnames(gsva.es_avg) <- levels(clinical$Global_ID)
gsva.es_avg <- gsva.es_avg %>% t() %>% scale() %>% t()
summary(as.numeric(gsva.es_avg))

#画图
heatmap<-gsva.es_avg
range(heatmap)
bk <- c(seq(-1,-0.01,by=0.01),seq(0,1,by=0.01))
pheatmap(heatmap, 
         color = colorRampPalette(c("#1E63A7", "white", "#AD152A"))(100),
         border_color = 'white',
         cluster_rows= T,
         cluster_cols = F,
         cellwidth = 29,
         cellheight = 11,
         fontsize = 9,
         angle_col = 45,
         show_colnames = TRUE,
         show_rownames = TRUE
)


##相关性

data_gsva<-t(gsva.es[,rownames(clinical)])[,c("PYROPTOSIS UP","PYROPTOSIS DOWN")]
data_gene<-as.matrix(t(datExpr[,rownames(clinical)])[,'APOE'])
colnames(data_gene)<-'APOE'
all(rownames(data_gene)==rownames(data_gsva))

c.mat.all<-c()
p.mat.all<-c()
for (i in levels(clinical$Global_ID)) {
  data_gene_sub <- data_gene[clinical$Global_ID==i,]
  data_gsva_sub <- data_gsva[clinical$Global_ID==i,]
  corr_sub<-corr.test(data_gsva_sub,data_gene_sub, method = "pearson",adjust= "fdr")
  c.mat.sub <-corr_sub$r
  p.mat.sub <-corr_sub$p.adj
  if(!is.null(p.mat.sub)){
    sssmt <- p.mat.sub< 0.001
    p.mat.sub[sssmt] <- '***'
    ssmt <- p.mat.sub > 0.001& p.mat.sub< 0.01
    p.mat.sub[ssmt] <- '**'
    smt <- p.mat.sub > 0.01& p.mat.sub < 0.05
    p.mat.sub[smt] <- '*'
    p.mat.sub[!sssmt&!ssmt&!smt]<- 'ns'
  } else{
    p.mat.sub <- F
  }
  c.mat.all <- rbind(c.mat.all,c.mat.sub)
  p.mat.all <- rbind(p.mat.all,p.mat.sub)
}

c.mat.fil <- data.frame(row.names = c('hm8','hm10','hm11','hm12','hm13'),up=c.mat.all[c(1,3,5,7,9),],down=c.mat.all[c(2,4,6,8,10)]) %>% as.matrix()
p.mat.fil <- data.frame(row.names = c('hm8','hm10','hm11','hm12','hm13'),up=p.mat.all[c(1,3,5,7,9),],down=p.mat.all[c(2,4,6,8,10)]) %>% as.matrix()
dim(c.mat.fil)

range(c.mat.fil)

datBar <- cbind(melt(c.mat.fil),melt(p.mat.fil))[,c(1,2,3,6)]
colnames(datBar) <- c('cell','group','cor','p')
datBar$cell <- c('hm08','hm10','hm11','hm12','hm13','hm08.','hm10.','hm11.','hm12.','hm13.')
ggplot(datBar,aes(x=cell,y=cor,group=group,fill=group))+
  geom_hline(yintercept = 0,color='black',linetype=c('dashed'))+
  geom_segment(aes(x=cell,xend=cell,y=0,yend=cor),size=1,color='darkgrey',linetype="solid")+
  geom_point(size=12,shape=21,color=ifelse(datBar$group=='up','#AD152A','#1E63A7'),
             fill=ifelse(datBar$p=='ns','white',ifelse(datBar$group=='up','#AD152A','#1E63A7')))+
  ylab('Correlation')+xlab('Cluster')+
  theme_classic()+
  theme(axis.text.x=element_text(angle = 0,family = 'sans',vjust = 1.5,hjust = -0.01,face = 'bold'),)+
  ggtitle("Correlation between APOE and Pyroptosis") 

