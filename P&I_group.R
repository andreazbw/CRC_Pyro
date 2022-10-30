library(survival)
library(survminer)
library(dplyr)
library(tidyverse)
library(ggthemes)
library(DESeq2)
library(ggrepel)
library(ggpubr)
library(VennDiagram)
library(patchwork)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(ggforce)
library(reshape2)
library(plyr)

datSurv<-readRDS("E:/R/UCSC-CRC/monocle/datSurv.rds")
estimate<-readRDS("E:/R/UCSC-CRC/estimate/immune_group.rds")
estimate<-estimate[rownames(datSurv),]
all(rownames(estimate)==rownames(datSurv))
tumor_tpm<-readRDS("E:/R/UCSC-CRC/download/tumor_tpm.rds")
tumor_count<-readRDS("E:/R/UCSC-CRC/download/tumor_count.rds")


##取双定义组

survival_mix<-estimate
survival_mix$num<-datSurv$cds.phenoData.data...Cluster...
survival_mix<-unite(survival_mix,"num_group",num,cluster,sep = "&")
write.csv(survival_mix,"crc_sankey.csv")
survival_mix$num_group[which(survival_mix$num_group %in% "Pyroptosis_L&Immunity_H")] <-"mix"
survival_mix$num_group[which(survival_mix$num_group %in% "Pyroptosis_H&Immunity_L")] <-"mix"
survival_mix<-survival_mix[!survival_mix$num_group %in% "mix",]
sfit<-survfit(Surv(survival_mix$OS.time,survival_mix$OS)~num_group,data = survival_mix)

ggsurvplot(sfit,
           legend.title='Pyroptosis&Immunity Groups',
           conf.int = F,
           pval = T,
           pval.method = T,
           combine = T,
           pval.size = 8,
           risk.table = T,
           risk.table.title = "",
           palette = c("#A65866","#9CBCD9"),
           legend.labs = c("Pyroptosis_H&Immunity_H","Pyroptosis_L&Immunity_L"),
           tables.height = 0.2,
           xlab = "Time In Days",
           ylab = "Overall Survival",
           surv.median.line = "hv",
           ncensor.plot = F)

saveRDS(survival_mix,"survival_mix.rds")


##sankey图

crc_sankey <- read.csv('crc_sankey.csv',header = T,row.names = 1,check.names = F)
sfit<-survfit(Surv(crc_sankey$OS.time,crc_sankey$OS)~num_group,data = crc_sankey)
ggsurvplot(sfit,
           legend.title='Pyroptosis&Immunity Groups',
           conf.int = F,
           pval = T,
           pval.method = T,
           combine = T,
           pval.size = 8,
           risk.table = T,
           risk.table.title = "",
           palette = c("#A65866","#76A660","#D98841","#9CBCD9"),
           legend.labs = c("Pyroptosis_H&Immunity_H","Pyroptosis_H&Immunity_L","Pyroptosis_L&Immunity_H","Pyroptosis_L&Immunity_L"),
           tables.height = 0.2,
           xlab = "Time In Days",
           ylab = "Overall Survival",
           surv.median.line = "hv",
           ncensor.plot = F)

crc_sankey <- table(crc_sankey$num,crc_sankey$num_group,crc_sankey$group) %>% 
  melt()
colnames(crc_sankey) <- c('Pyroptosis group','P&I group','Immunity group','value')

crc_sankey <- gather_set_data(crc_sankey,1:3)
crc_sankey$x <- factor(crc_sankey$x,levels = c('Pyroptosis group','P&I group','Immunity group'))

ggplot(crc_sankey, aes(x, id = id, split = y, value = value)) +
  geom_parallel_sets(aes(fill = `P&I group`), alpha = 0.7, axis.width = 0.1) +
  geom_parallel_sets_axes(axis.width = 0.2,fill="#D9D8D7",color="white") +
  geom_parallel_sets_labels(colour = 'black',angle = 0) +
  scale_fill_manual(values = c("#A65866","#76A660","#D98841","#9CBCD9"))+
  scale_x_discrete(name = "")+
  scale_y_continuous(name = "Number of patients")+
  theme_bw()+
  theme(plot.title = element_text(size = 20, face =  "bold",family = 'sans',hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11,family = 'sans',face = 'bold'),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())

#卡方检验
crc_sankey <- read.csv('crc_sankey.csv',header = T,row.names = 1,check.names = F)
crc_chisq <- table(crc_sankey$num,crc_sankey$group)
result <- chisq.test(crc_chisq)
result[["p.value"]]

#分组与肿瘤进展的关系
crc_sankey <- read.csv('crc_sankey.csv',header = T,row.names = 1,check.names = F)
clinical<-readRDS("E:/R/UCSC-CRC/download/clinical.rds")

pat_id <- intersect(rownames(clinical),rownames(crc_sankey))
crc_cli <- crc_sankey[pat_id,] %>% 
  subset(.,num_group%in%c('Pyroptosis_L&Immunity_L','Pyroptosis_H&Immunity_H')) 
crc_cli <- cbind(crc_cli,clinical[rownames(crc_cli),'pathologic_stage'])          
crc_cli <- crc_cli[,c(4,7)]
colnames(crc_cli) <- c('group1','group2')
crc_cli <- na.omit(crc_cli)
crc_cli <- crc_cli[crc_cli$group2!='[Discrepancy]',]
crc_cli$group2 <- mapvalues(crc_cli$group2,from = unique(crc_cli$group2),
                            to = c('Stage III','Stage I','Stage II','Stage III',
                                   'Stage IV','Stage II','Stage IV','Stage II',
                                   'Stage IV','Stage III','Stage II','Stage III','Stage I'))
crc_chisq <- table(crc_cli$group1,crc_cli$group2)
result <- chisq.test(crc_chisq)
result[["p.value"]]

#折线图
datline <- as.data.frame(table(crc_cli$group1,crc_cli$group2)) %>% spread(.,key = Var1,value = Freq)
datline[,-1] <- (datline[,-1]/rowSums(datline[,-1])) 
datline <- melt(datline,id.vars = 'Var2')
colnames(datline) <- c('stage','group','proportion')
ggplot(datline , mapping = aes(x = stage, y = proportion,group=group,color=group)) + geom_line(cex=1)+
  scale_color_manual(values = c("#A65866","#9CBCD9"))+
  geom_point(pch=20,cex=5)+
  theme_classic()+
  theme(axis.text = element_text(face='bold',size = 10,angle = 0,family = 'sans',colour = 'black'),
        axis.text.x =element_text(vjust = 0.5),
        plot.title = element_text(hjust = 0.5,size = 15))+
  xlab('')+ylab('Proportion')+
  ggtitle('Changes in different status of patients at each tumor stage')


##焦亡核心基因两组表达情况

boxplot<-t(tumor_tpm[,rownames(survival_mix)])
all(rownames(boxplot)==rownames(survival_mix))

boxplot<-as.data.frame(boxplot[,c('NLRP1','NLRP3','AIM2','GSDMC','GSDMD','DFNA5','CASP1','CASP4','CASP5','IL1B','IL6','IL18')])
boxplot$group<-survival_mix$num_group
boxplot$id<-rownames(boxplot)
boxplot<-boxplot %>% gather(key="gene",value="expression",1:12) %>% 
  dplyr::select(id,gene,expression,group) 

ggplot(aes(x = gene, y = expression,fill = group),data = boxplot) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Gene expression")+
  scale_x_discrete(name = "") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face =  "bold",family = 'sans',hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.y=element_text(size = 11,family = 'sans',face = 'bold'),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
  scale_fill_manual(values = c("#A65866","#9CBCD9"))+
  stat_compare_means(method = 'wilcox.test',label = 'p.signif',hide.ns = TRUE,label.x = 1.5)
  

##差异分析

des_count<-tumor_count[,rownames(survival_mix)]
all(colnames(des_count)==rownames(survival_mix))

group_list=survival_mix$num_group
group_list[which(group_list %in% c("Pyroptosis_L&Immunity_L"))] <-"Low"
group_list[which(group_list %in% c("Pyroptosis_H&Immunity_H"))] <-"High"
group_list <- factor(group_list,levels = c('High','Low'))
coldata <- data.frame(row.names=colnames(des_count),tissue=group_list)
dds <- DESeqDataSetFromMatrix(countData = des_count, 
                              colData = coldata, 
                              design = ~tissue)
dds$tissue<- relevel(dds$tissue, ref = "Low")
dds <- DESeq(dds)

DESeq2_result <- as.data.frame(results(dds))
DESeq2<-subset(DESeq2_result,DESeq2_result$padj< 0.001 & abs(DESeq2_result$log2FoldChange)>1)

saveRDS(DESeq2,"DESeq2.rds")

DESeq2_immune<-readRDS("E:/R/UCSC-CRC/estimate/DESeq2.rds")
DESeq2_immune<-subset(DESeq2_immune,DESeq2_immune$padj< 0.001 & abs(DESeq2_immune$log2FoldChange)>1)
DESeq2_pyro<-readRDS("E:/R/UCSC-CRC/monocle/DESeq2.rds")
DESeq2_pyro<-subset(DESeq2_pyro,DESeq2_pyro$padj< 0.001 & abs(DESeq2_pyro$log2FoldChange)>1)

DEG.pyro<-DESeq2[intersect(rownames(DESeq2_immune),rownames(DESeq2)),]
DEG.pyro<-DEG.pyro[intersect(rownames(DESeq2_pyro),rownames(DEG.pyro)),]

saveRDS(DEG.pyro,"DEG.pyro.rds")

#韦恩图
A<-rownames(DESeq2_immune)
B<-rownames(DESeq2_pyro)
C<-rownames(DESeq2)
D<-venn.diagram(x= list(A = A,B = B,C = C), 
             category.names = c("Immunity","Pyroptosis","Pyroptosis&Immunity"),
             filename = NULL, height = 6000, width = 6000,resolution =3000,
             col="white",fill=c("#F26E50","#F2AC29","#A65866"),alpha = 0.50, cex=0.45, cat.cex=0.45)

pdf("1-E.pdf",width = 2,height = 2)
grid.draw(D)
dev.off()

#火山图
DESeq2_vol<-DESeq2_result
DESeq2_vol$logP <- -log10(DESeq2_vol$padj)
DESeq2_vol$Group <- 'non-significant'
DESeq2_vol$Group[match(rownames(DEG.pyro[DEG.pyro$log2FoldChange>0,]),rownames(DESeq2_vol))] <- 'up-regulated'
DESeq2_vol$Group[match(rownames(DEG.pyro[DEG.pyro$log2FoldChange<0,]),rownames(DESeq2_vol))] <- 'down-regulated'
table(DESeq2_vol$Group)

DESeq2_vol$Label <- ''
DESeq2_vol$Label[match(rownames(DEG.pyro),rownames(DESeq2_vol))] <- rownames(DEG.pyro)

ggscatter(DESeq2_vol,
          x='log2FoldChange',
          y='logP',
          color = 'Group',
          palette = c('#9CBCD9','#BBBBBB','#A65866'),
          size = 1,
          #label = 'Label',
          repel = F,
          font.label = 8)+
  theme_base()+
  geom_hline(yintercept = 0,linetype='dashed')+
  geom_vline(xintercept = c(-1,1),linetype='dashed')

#焦亡相关基因箱线图
boxplot<-t(tumor_tpm[,rownames(survival_mix)])
all(rownames(boxplot)==rownames(survival_mix))
boxplot<-as.data.frame(boxplot[,c("HLA-DMA","HOXC8","HAMP","APOE","SIRPB2")])
boxplot$`P&I group`<-survival_mix$num_group
boxplot$id<-rownames(boxplot)
boxplot<-boxplot %>% gather(key="gene",value="expression",5) %>% 
  dplyr::select(id,gene,expression,`P&I group`) 

ggplot(aes(x = gene, y = expression,fill = `P&I group`),data = boxplot) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "SIRPB2")+
  scale_x_discrete(name = "") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face =  "bold",family = 'sans',hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.y=element_text(size = 11,family = 'sans',face = 'bold'),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
  scale_fill_manual(values = c("#A65866","#9CBCD9"))+
  stat_compare_means(method = 'wilcox.test',label = 'p.signif',hide.ns = TRUE,label.x = 1.5)


##生存分析

survival_os<-survival_mix[survival_mix$num_group=="Pyroptosis_H&Immunity_H",]
survival_os$OS=as.numeric(survival_os$OS)
survival_os$OS.time=as.numeric(survival_os$OS.time)
tpm_pyro<-tumor_tpm[intersect(rownames(DEG.pyro),rownames(tumor_tpm)),]
tpm_pyro<-tpm_pyro[,intersect(colnames(tpm_pyro),rownames(survival_os))]
all(colnames(tpm_pyro)==rownames(survival_os))

my.surv <- Surv(survival_os$OS.time,survival_os$OS)
cox_results <- apply(tpm_pyro, 1, function(gene){
  survival_os$gene<-gene
  survival_os$group=ifelse(gene>median(gene),'high','low')
  m=coxph(my.surv ~ group,data=survival_os)
  
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  
  tmp <- round(cbind(coef = beta, se = se, z = beta/se, p = 1 - pchisq((beta/se)^2, 1),
                     HR = HR, HRse = HRse,
                     HRz = (HR - 1) / HRse, HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
                     HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
                     HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
  return(tmp['grouplow',])
  
})
cox_results=t(cox_results)
cox_results0.05<-cox_results[cox_results[,4]<0.05,]

#画图
survival_mix2<-cbind(survival_os,t(tpm_pyro))
gene<-"HLA-DMA"
group<-ifelse(survival_mix2[,gene]>median(survival_mix2[,gene]),"High","Low")
sfit2<-survfit(Surv(survival_mix2$OS.time,survival_mix2$OS)~group,data = survival_mix2)
ggsurvplot(sfit2,
           legend.title='HLA-DMA',
           conf.int = F,
           pval = T,
           pval.method = T,
           combine = T,
           pval.size = 8,
           risk.table = T,
           risk.table.title = "",
           palette = c("#D75929","#403C2B"),
           legend.labs = c("High","Low"),
           tables.height = 0.2,
           xlab = "Time In Days",
           ylab = "Overall Survival",
           surv.median.line = "hv",
           ncensor.plot = F)


##GO和KEGG富集分析

write.csv(hg,'DEG.pyro.csv')

gomf <- read.table('GOBP.txt',sep = '\t',header = T,check.names = F)
gomf <- gomf %>% subset(.,FDR<0.05)
gomf <- gomf[order(gomf$Count,decreasing = T),]
gomf <- gomf[1:10,]

datbar <- gomf[,c(2,3,13)]
datbar$Term <- factor(datbar$Term,levels = rev(datbar$Term))
ggplot(datbar,aes(x=Count,y=Term,fill=FDR))+
  geom_bar(stat = 'identity',mapping = aes(fill=FDR))+
  theme_test()+
  scale_colour_gradient(low="#8C423B",high="#F2C094",aesthetics = "fill")+
  theme(axis.text = element_text(face='bold',size = 8,angle = 0,family = 'sans',colour = 'black'),
        legend.position = c(0.9,0.5))+
  xlab('Gene Counts')+ylab('')
  
