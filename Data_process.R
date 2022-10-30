library("stringr")
library("rtracklayer")
library(dplyr)


##注释文件

annotation = import('gencode.v22.annotation.gtf')
annotation = as.data.frame(annotation)
annotation = annotation[,c("type","gene_id","gene_type","gene_name")]
annotation_all<-annotation[!duplicated(annotation$gene_id),]
annotation_mRNA<-annotation_all[annotation_all$gene_type=="protein_coding",]

#除去小数点后数字
annotation_mRNA$gene_name=unlist(str_split(annotation_mRNA$gene_name,"[.]",simplify=T))[,1]
annotation_all$gene_name=unlist(str_split(annotation_all$gene_name,"[.]",simplify=T))[,1]
annotation_mRNA<-annotation_mRNA[!duplicated(annotation_mRNA$gene_name),]
annotation_all<-annotation_all[!duplicated(annotation_all$gene_name),]

#log2(count+1)
coad_count<-read.table("TCGA-COAD.htseq_counts.tsv",header = T,row.names = 1)
read_count<-read.table("TCGA-READ.htseq_counts.tsv",header = T,row.names = 1)
colnames(coad_count)<-gsub(pattern = "[.]",replacement = "-",x = colnames(coad_count))
colnames(read_count)<-gsub(pattern = "[.]",replacement = "-",x = colnames(read_count))
all(rownames(coad_count)==rownames(read_count))
CRC_count<-cbind(coad_count,read_count)

#注释
CRC_count<-CRC_count[match(annotation_mRNA$gene_id,rownames(CRC_count)),]
all(rownames(CRC_count)==annotation_mRNA$gene_id)
rownames(CRC_count)<-annotation_mRNA$gene_name
colnames(CRC_count)<-substr(colnames(CRC_count),1,15)
CRC_count<-round(2^CRC_count-1,digits = 0)

#log2(fpkm+1)
coad_fpkm<-read.table("TCGA-COAD.htseq_fpkm.tsv",header = T,row.names = 1)
read_fpkm<-read.table("TCGA-READ.htseq_fpkm.tsv",header = T,row.names = 1)
colnames(coad_fpkm)<-gsub(pattern = "[.]",replacement = "-",x = colnames(coad_fpkm))
colnames(read_fpkm)<-gsub(pattern = "[.]",replacement = "-",x = colnames(read_fpkm))
all(rownames(coad_fpkm)==rownames(read_fpkm))
CRC_fpkm<-cbind(coad_fpkm,read_fpkm)
colnames(CRC_fpkm)<-substr(colnames(CRC_fpkm),1,15)
CRC_tpm <- 2^CRC_fpkm-1
fpkmToTpm <- function(fpkm)
{exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
CRC_tpm <- apply(CRC_tpm,2,fpkmToTpm)
colSums(CRC_tpm)
CRC_tpm <- log2(CRC_tpm+1)

#注释
CRC_tpm_all<-CRC_tpm[match(annotation_all$gene_id,rownames(CRC_tpm)),]
all(rownames(CRC_tpm_all)==annotation_all$gene_id)
rownames(CRC_tpm_all)<-annotation_all$gene_name
CRC_tpm<-CRC_tpm[match(annotation_mRNA$gene_id,rownames(CRC_tpm)),]
all(rownames(CRC_tpm)==annotation_mRNA$gene_id)
rownames(CRC_tpm)<-annotation_mRNA$gene_name


##临床信息处理

clinical<-read.table("TCGA.COADREAD.sampleMap_COADREAD_clinicalMatrix",header = T,row.names = 1,na.strings = "",sep = "\t")
survival<-read.table("survival_COADREAD_survival.txt",header = T,row.names = 1,na.strings = "",sep = "\t")
survival<-survival[!duplicated(survival$X_PATIENT),]
clinical<-clinical[intersect(rownames(clinical),rownames(survival)),]
all(rownames(clinical)==rownames(survival))
all(match(colnames(CRC_count),colnames(CRC_tpm)))

CRC_tpm <- CRC_tpm[,intersect(colnames(CRC_tpm),rownames(survival))]
CRC_tpm <- CRC_tpm[,intersect(colnames(CRC_tpm),rownames(clinical))]
clinical <- clinical[intersect(colnames(CRC_tpm),rownames(clinical)),]
survival <- survival[intersect(colnames(CRC_tpm),rownames(survival)),]
all(rownames(clinical)==rownames(survival))
all(colnames(CRC_tpm)==rownames(survival))
table(clinical$sample_type)
CRC_tpm_all <- CRC_tpm_all[,intersect(colnames(CRC_tpm),colnames(CRC_tpm_all))]
CRC_count <- CRC_count[,intersect(colnames(CRC_tpm),colnames(CRC_count))]

#癌症样本
tumor_tpm<-CRC_tpm[,intersect(colnames(CRC_tpm),
                              rownames(clinical[which(clinical$sample_type %in% c('Metastatic','Primary Tumor','Recurrent Tumor')),]))]
table(substr(colnames(tumor_tpm),14,15))
tumor_tpm_all<-CRC_tpm_all[,intersect(colnames(CRC_tpm_all),
                              rownames(clinical[which(clinical$sample_type %in% c('Metastatic','Primary Tumor','Recurrent Tumor')),]))]
tumor_count<-CRC_count[,intersect(colnames(CRC_count),
                              rownames(clinical[which(clinical$sample_type %in% c('Metastatic','Primary Tumor','Recurrent Tumor')),]))]

#生存时间os
survival_os<-survival[!survival$OS.time=="0",]
survival_os<-survival_os[!is.na(survival_os$OS.time),]
survival_os<-survival_os[!is.na(survival_os$OS),]
survival_os <- survival_os[intersect(colnames(tumor_tpm),rownames(survival_os)),]


##单细胞数据处理

clinical <- read.table("CRC.Leukocyte.10x.Metadata.txt",sep="\t",header = T)
CRC_tpm <- read.table("CRC.Leukocyte.10x.TPM.txt",sep=" ",header = T)

CRC.T<-CRC_tpm[,substr(clinical$CellName,3,3)=="T"]
CRC_M<-CRC_tpm[,substr(clinical$Sub_ClusterID,2,2)=="M"]
CRC_T<-CRC_tpm[,substr(clinical$Sub_ClusterID,2,2)=="T"]


##保存数据

saveRDS(tumor_tpm,"tumor_tpm.rds")
saveRDS(tumor_count,"tumor_count.rds")
saveRDS(clinical,"clinical.rds")
saveRDS(survival_os,"survival_os.rds")
saveRDS(CRC.T,"CRC.T.rds")
saveRDS(CRC_M,"CRC_M.rds")
saveRDS(CRC_T,"CRC_T.rds")
