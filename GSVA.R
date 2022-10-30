library(GSVA)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(limma)
library(psych)
library(corrplot)

survival_mix<-readRDS("E:/R/UCSC-CRC/intersect/survival_mix.rds")
tumor_tpm<-readRDS("E:/R/UCSC-CRC/download/tumor_tpm.rds")
tumor_tpm<-tumor_tpm[,rownames(survival_mix)]
all(colnames(tumor_tpm)==rownames(survival_mix))


##Gsva

geneset <- read.csv('Gsva.csv',row.names = 1,header = FALSE,na.strings = '') %>% t()

GSVA_gs <- melt(geneset,measure.vars = colnames(geneset),variable.name = 'Description',value.name = 'gene') %>% na.omit()
GSVA_gs <- split(GSVA_gs$gene, GSVA_gs$Var2)

gsva.es <- gsva(tumor_tpm, GSVA_gs,verbose=TRUE,method="gsva",kcdf="Gaussian",min.sz=1,max.sz=1000)


##差异分析

survival_os<-survival_mix
survival_os$num_group[which(survival_os$num_group %in% "Pyroptosis_H&Immunity_H")] <-"high"
survival_os$num_group[which(survival_os$num_group %in% "Pyroptosis_L&Immunity_L")] <-"low"

group_list=survival_os$num_group
group_list <- factor(group_list,levels = c('high','low'))

design_matrix <- model.matrix(~0+group_list)
colnames(design_matrix) <- levels(group_list)
rownames(design_matrix) <- rownames(survival_os)

contracts = paste(levels(group_list),collapse = "-")
contracts_matrix <- makeContrasts(contracts,levels=design_matrix)

fit <- lmFit(gsva.es, design_matrix)
fit <- contrasts.fit(fit,contracts_matrix)
fit <- eBayes(fit)
DEG <- topTable(fit,coef=contracts,n=Inf)
DEG = na.omit(DEG)

#画图
datBar<- data.frame(DEG$t,DEG$P.Value,rownames(DEG))
colnames(datBar) <- c('t','padj','pathway')
rownames(datBar) <- rownames(DEG)
datBar <- datBar %>% mutate(.,group=as.factor(ifelse(padj>=0.05,'ns',ifelse(t>0,'enriched in High P&I group','enriched in Low P&I group'))))
datBar <- datBar %>% mutate(.,pathway=factor(pathway,levels = pathway[order(t,decreasing = F)]))

ggplot(data=datBar,mapping=aes(x=t,y=pathway,fill=group))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = c('#A65866','#9CBCD9','#CCCCCC'))+
  geom_vline(xintercept = c(-2,2),linetype="dashed")+
  theme_test()+
  ggtitle('')+
  xlab('t value of GSVA score,High P&I group versus Low P&I group')+ylab('')+
  theme(axis.text = element_text(face='bold',size = 10,angle = 0,family = 'sans',colour = 'black'),
        axis.text.x =element_text(vjust = 0.5,angle = 0),
        plot.title = element_text(hjust = 0.5,size = 6))


##相关性

#基因与通路
data_gsva<-t(as.matrix(gsva.es))
data_gene<-t(tumor_tpm[c("HLA-DMA","HOXC8","HAMP","APOE","SIRPB2"),])
all(rownames(data_gene)==rownames(data_gsva))

corr<-corr.test(data_gsva,data_gene,method = "pearson",adjust= "fdr")
c.mat <-corr$r
p.mat <-corr$p.adj

corrplot(c.mat,method = "square",tl.col = 'black',col.lim = c(-1,1),
         cl.ratio = 0.3,tl.cex = 0.7,cl.cex = 0.7,
         col = c(colorRampPalette(colors = c("#1E63A7","white"))(100),
                 colorRampPalette(colors = c("white","#AD152A"))(100)),
         p.mat  = p.mat,insig = 'label_sig',
         sig.level = c(.001, .01, .05), pch.cex = .9,  pch.col = 'black')

#焦亡通路与其它通路
data_gsva<-t(as.matrix(gsva.es[!rownames(gsva.es) %in% c("PYROPTOSIS","FERROPTOSIS","NECROPTOSIS","APOPTOSIS"),]))
data_gene<-as.matrix(gsva.es["PYROPTOSIS",])
colnames(data_gene)<-"PYROPTOSIS"
all(rownames(data_gene)==rownames(data_gsva))

corr<-corr.test(data_gene,data_gsva, method = "pearson",adjust= "fdr")
c.mat <-corr$r
p.mat <-corr$p.adj

datBar<-as.data.frame(t(c.mat))
datBar<-cbind(datBar,t(p.mat))
datBar$Pathways<-rownames(datBar)
colnames(datBar)<-c("Correlation","Pvalues","Pathways")
datBar <- datBar %>% mutate(.,Pathways=factor(Pathways,levels = Pathways[order(Correlation,decreasing = T)]))

ggplot(datBar,aes(x=Pathways,y=Correlation))+
  geom_hline(yintercept = c(-1,0,1),color=c('grey','black','grey'),linetype=c('dashed'))+
  geom_segment(aes(x=Pathways,xend=Pathways,y=0,yend=Correlation),size=1,color='darkgrey',linetype="solid")+
  geom_point(size=8,color=ifelse(datBar$Correlation>0,'#A65866','#9CBCD9'),
             fill=ifelse(datBar$Pvalues>0.05,'white',ifelse(datBar$Correlation>0,'#A65866','#9CBCD9')),
             shape=21)+
  ylab('Correlation')+xlab('Pathways')+
  theme_classic()+
  theme(axis.text.x=element_text(angle = -75,family = 'sans',vjust = 1.5,hjust = -0.01,face = 'bold'),)+
  ggtitle("Correlation between Pyroptosis and other pathways") 



