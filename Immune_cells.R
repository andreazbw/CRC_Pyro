library(tibble)
library(GSVA)
library(dplyr)
library(pheatmap)
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(survival)
library(survminer)
library(psych)
library(reshape2)
library(corrplot)
library(reshape2)

survival_mix<-readRDS("E:/R/UCSC-CRC/intersect/survival_mix.rds")
tumor_tpm<-readRDS("E:/R/UCSC-CRC/download/tumor_tpm.rds")
tumor_tpm<-tumor_tpm[,rownames(survival_mix)]
all(colnames(tumor_tpm)==rownames(survival_mix))


##ssGSEA

tisidb <- read.csv("CellReports_Innate.txt",header = F,sep = "\t",)
tisidb <- tisidb %>%column_to_rownames("V1")%>%t()
a <- tisidb
a <- a[1:nrow(a),]
set <- colnames(a)
l <- list()
for (i in set) {
  x <-  as.character(a[,i])
  x <- x[nchar(x)!=0]
  x <-  as.character(x)
  l[[i]] <-x
}

tpm_immune<- gsva(tumor_tpm, l,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
tisidb2 <- read.csv("Gsva.csv",header = F,row.names = 1,check.names = F)['PYROPTOSIS',] %>% t()
tisidb2 <- list(as.vector(tisidb2$PYROPTOSIS))
pyroptosis <- gsva(tumor_tpm,tisidb2,method="ssgsea",kcdf="Gaussian",abs.ranking=TRUE)
rownames(pyroptosis)<-"pyroptosis"
pyroptosis<-t(pyroptosis)

anno<-arrange(survival_mix,num_group)
anno<-as.data.frame(anno[,c(2,4)])
anno$OS<-as.factor(anno$OS)
colnames(anno)<-c("OS","P&I group")

#热图
heatmap<-tpm_immune[,rownames(anno)]
all(colnames(heatmap)==rownames(anno))
annotation_color <- list(OS=c("0"='#BFBFBF',"1"='#363432'),
                         `P&I group`=c("Pyroptosis_H&Immunity_H"="#A65866","Pyroptosis_L&Immunity_L"="#9CBCD9"))
pheatmap(heatmap, 
         color = colorRampPalette(c("#9CBCD9", "white", "#A65866"))(100),
         scale = "row",
         cluster_rows=T,
         cluster_cols = F,
         cellwidth = 1,
         cellheight = 30,
         fontsize = 8,
         annotation_col = anno,
         annotation_colors = annotation_color,
         show_colnames = FALSE,
         show_rownames = TRUE
)

#箱线图
boxplot<-as.data.frame(t(tpm_immune))
boxplot<-boxplot[rownames(anno),]
all(rownames(boxplot)==rownames(anno))
boxplot$group<-anno$`P&I group`
boxplot$id<-rownames(boxplot)
boxplot<-boxplot %>% gather(key="cell",value="expression",1:14) %>% 
  dplyr::select(id,cell,expression,group) 
ggplot(aes(x = cell, y = expression,fill = group),data = boxplot) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "The expression of Cells(log2 TPM)")+
  scale_x_discrete(name = "Cell types") +
  ggtitle("Different Cell types of CRC patients") +
  theme_bw() +
  theme(plot.title = element_text(size = 20, face =  "bold",family = 'sans',hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(angle = -75,family = 'sans',vjust = 1.5,hjust = -0.01,face = 'bold'),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
  scale_fill_manual(values = c("#A65866", "#9CBCD9"))+
  stat_compare_means(method = 'wilcox.test',label = 'p.signif',hide.ns = TRUE)
ggsave(paste0("boxplot",".pdf"),width = 20,height = 10)

cell<-tpm_immune[!rownames(tpm_immune) %in% "CD56bright natural killer cell",]
saveRDS(cell,"tumor_immune.rds")


##相关性

tpm_pyro<-t(tumor_tpm[c("HLA-DMA","HOXC8","HAMP","APOE","SIRPB2",'NLRP3','CASP1','GSDMD'),])
all(rownames(tpm_pyro)==rownames(t(cell)))
tpm_cell<-t(cell)
all(rownames(tpm_pyro)==rownames(tpm_cell))

corr<-corr.test(tpm_cell,tpm_pyro, method = "pearson",adjust= "fdr")
c.mat <-corr$r
p.mat <-corr$p.adj

corheatmap <-melt(c.mat,value.name= "cor")
corheatmap$pvalue <- as.vector(p.mat)
corheatmap<-corheatmap[which(corheatmap$pvalue<0.05),]
corheatmap<-corheatmap[which(corheatmap$cor>0.7),]
table(corheatmap$Var2)

corrplot(c.mat,tl.col = 'black',col.lim = c(-1,1),
         cl.ratio = 0.2,tl.cex = 0.7,cl.cex = 0.7,
         col = c(colorRampPalette(colors = c("#1E63A7","white"))(100),
                 colorRampPalette(colors = c("white","#AD152A"))(100)),
         p.mat  = p.mat,insig = 'label_sig',
         sig.level = c(.001, .01, .05), pch.cex = .9,  pch.col = 'black')

datbox<-as.data.frame(cell["Macrophage",])
colnames(datbox)<-"Macrophage"
all(rownames(datbox)==rownames(pyroptosis))
datbox<-cbind(datbox,pyroptosis)
ggscatter(datbox,x = "pyroptosis", y = 'Macrophage',size = 2,
         xlab = 'Pyroptosis',ylab = 'Macrophage',color = '#A65866',
         add = "reg.line", conf.int = TRUE,    
         add.params = list(fill = "darkgray"))+
  stat_cor(method = "pearson")

