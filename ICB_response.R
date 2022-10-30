library(easier)
library(maftools)
library(ggpubr)
library(ggstatsplot)
library(reshape2)
library(tidyverse)
library(psych)
library(corrplot)

crc_counts <- readRDS('tumor_count.rds')
crc_logtpm <- readRDS('tumor_tpm.rds')
crc_status <- readRDS('survival_mix.rds')

coad_maf <- read.maf('TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.gz',isTCGA = T)
read_maf <- read.maf('TCGA.READ.mutect.faa5f62a-2731-4867-a264-0e85b7074e87.DR-10.0.somatic.maf.gz',isTCGA = T)


##数据处理

coad_tmb <- tmb(coad_maf)
read_tmb <- tmb(read_maf)
crc_tmb <- rbind(coad_tmb,read_tmb)

TMB <- crc_tmb$total_perMB
names(TMB) <- crc_tmb$Tumor_Sample_Barcode

colnames(crc_counts) <- str_sub(colnames(crc_counts),1,12)
crc_counts <- crc_counts[,unique(colnames(crc_counts))]

colnames(crc_logtpm) <- str_sub(colnames(crc_logtpm),1,12)
crc_tpm <- crc_logtpm[,unique(colnames(crc_logtpm))]

all(colnames(crc_logtpm)==colnames(crc_counts))

TMB <- TMB[intersect(names(TMB),crc_status$X_PATIENT)]
crc_status <- crc_status[crc_status$X_PATIENT%in%names(TMB),]
rownames(crc_status) <- crc_status$X_PATIENT
crc_status <- crc_status[names(TMB),]
all(rownames(crc_status)%in%names(TMB))

crc_counts <- crc_counts[,names(TMB)]
crc_logtpm <- crc_logtpm[,names(TMB)]

crc_tpm <- 2^crc_logtpm-1

RNA_counts <- crc_counts
RNA_tpm <- crc_tpm


##Easier

hallmarks_of_immune_response <- c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")
immune_response_scores <- compute_scores_immune_response(RNA_tpm = RNA_tpm, 
                                                         selected_scores = hallmarks_of_immune_response)#tpm数据
cell_fractions <- compute_cell_fractions(RNA_tpm = RNA_tpm)
pathway_activities <- compute_pathway_activity(RNA_counts = RNA_counts,
                                               remove_sig_genes_immune_response = TRUE)
tf_activities <- compute_TF_activity(RNA_tpm = RNA_tpm)
lrpair_weights <- compute_LR_pairs(RNA_tpm = RNA_tpm,
                                   cancer_type = "CRC")
ccpair_scores <- compute_CC_pairs(lrpairs = lrpair_weights, 
                                  cancer_type = "CRC")

predictions <- predict_immune_response(pathways = pathway_activities,
                                       immunecells = cell_fractions,
                                       tfs = tf_activities,
                                       lrpairs = lrpair_weights,
                                       ccpairs = ccpair_scores,
                                       cancer_type = 'CRC', 
                                       verbose = TRUE)

easier_derived_scores <- retrieve_easier_score(predictions_immune_response = predictions,
                                               TMB_values = TMB,
                                               easier_with_TMB = c("weighted_average", 
                                                                   "penalized_score"),
                                               weight_penalty = 0.5)


#画图
datbox <- data.frame(row.names = names(TMB),likelihood=easier_derived_scores$w_avg_score,group=crc_status$num_group)
datbox$group[which(datbox$group %in% "high_High")] <-"Pyroptosis_H&Immunity_H"
datbox$group[which(datbox$group %in% "low_Low")] <-"Pyroptosis_L&Immunity_L"

ggplot(aes(x = group, y = likelihood,fill = group),data = datbox) +
          stat_boxplot(geom = 'errorbar',width=0.3,cex=1)+
          geom_boxplot(outlier.size = -1,width=0.5) +
          scale_y_continuous(name = "Response score",)+
          scale_x_discrete(name = "") +
          ggtitle("Likehood of response to immune therapy") +
          theme_classic()+
          scale_fill_manual(values = c("#A65866","#9CBCD9"))+
          theme(plot.title = element_text(size = 10, face =  "bold",family = 'sans',hjust = 0.5),
                text = element_text(size = 12),
                axis.title = element_text(face="bold"),
                axis.text.x=element_text(size = 10,angle = 45,family = 'sans',vjust = 1,hjust = 1,face = 'bold'),
                panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
          geom_signif(comparisons = list(c(1,2)),test="wilcox.test",annotations = '****',
                      tip_length = 0.02,size = 0.8,textsize = 6,y_position = 0.45)


