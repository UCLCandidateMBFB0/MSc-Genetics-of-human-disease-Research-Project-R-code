### This script is used for retrieving copy number data, merging CN and SNV data and calculating CCF values##
library(TCGAbiolinks)
library(maftools)
library(dplyr)
library(mixtools)
library(ggplot2)
library(hrbrthemes)
library(tidyr)
library(pastecs)
library(viridis)

purity <- read.delim("TCGA_mastercalls.abs_tables_JSedit.fixed.txt")
quiescence <-read.csv("PancancerQuiescenceGroups.csv")

############################# Cancer types ################################
# Download copy number data for LUSC:
query.maf.hg38 <- GDCquery(project = "TCGA-LUSC", 
                           data.category = "Copy Number Variation", 
                           data.type = "Copy Number Segment",
                           access = "open", 
                           legacy = FALSE)
GDCdownload(query.maf.hg38)
cnv <- GDCprepare(query.maf.hg38)
df.cnv <- data.frame(cnv)
# derive copy number value from segment mean:
df.cnv$CopyNumber <- 2*(2^cnv$Segment_Mean)
# add "chr" to chromosome numbers (will be needed later):
df.cnv$Chromosome <- paste0("chr",df.cnv$Chromosome)
#save(df.cnv, file="df.cnv.RData")

# Download MAF data for LUSC:
maf <- GDCquery_Maf("LUSC", pipelines = "mutect2") %>% read.maf
quiescence<-quiescence[quiescence$CancerType=="LUSC",]

# extract data for some patient IDs:
substitutionStats = titv(maf = maf, plot = FALSE, useSyn = TRUE)
maf_barcode<-substitutionStats$fraction.contribution[,1]
maf_barcode$Tumor_Sample_Barcode<-as.character(maf_barcode$Tumor_Sample_Barcode)
snvs_dormant <- data.frame(subsetMaf(maf = maf, 
                                     tsb =maf_barcode$Tumor_Sample_Barcode, 
                                     mafObj = FALSE))

#########
# Next, we need to be able to match TCGA IDs from the MAF and CN tables 
# - this is done via the first three bits of the TCGA ID.

# extract this information for both tables and store it in a new column:
snvs_dormant$Case <- sapply(snvs_dormant$Tumor_Sample_Barcode,
                            function(x) paste(strsplit(as.character(x),"-")[[1]][1:4],collapse="-"))
df.cnv$Case <- sapply(df.cnv$Sample,
                      function(x) paste(strsplit(x,"-")[[1]][1:4],collapse="-"))

purity$Case<- sapply(purity$sample,
                     function(x) paste(strsplit(x,"-")[[1]][1:4],collapse="-"))

snvs_dormant<- merge(snvs_dormant, purity,
                     by.x="Case", by.y="Case",
                     all.x=FALSE, all.y=FALSE)

quiescence$Case<- sapply(quiescence$Barcode,
                         function(x) paste(strsplit(x,"-")[[1]][1:4],collapse="-"))

#attach the quiescence data to the original dataframe
snvs_dormant <- merge(snvs_dormant, quiescence,
                      by.x="Case", by.y="Case",
                      all.x=FALSE, all.y=FALSE)

uni_patient<-unique(quiescence$Case)

###### for every patients in each cancer type #####
## Select one random patient to work on (this will need to be repeated for every patient):
patient=uni_patient[1]
snvs_dormant.selected <- snvs_dormant[which(snvs_dormant$Case==patient),]
df.cnv.selected <- df.cnv[which(df.cnv$Case == patient),]

#### Employ the genomic ranges package to intersect genomic intervals:
library(GenomicRanges)
# build a genomic ranges object for SNVs:
gr.snvs <- makeGRangesFromDataFrame(snvs_dormant.selected, keep.extra.columns = TRUE,
                                    seqnames.field = c("Chromosome"), 
                                    start.field = "Start_Position", 
                                    end.field = "End_Position")
# build a genomic ranges object for CN:
gr.cnvs <- makeGRangesFromDataFrame(df.cnv.selected, keep.extra.columns = TRUE,
                                    seqnames.field = c("Chromosome"), 
                                    start.field = "Start", 
                                    end.field = "End")
# intersect the two:
gr.overlap <- data.frame(findOverlaps(gr.snvs, gr.cnvs))

# merge the relevant columns from the intersection:
gr.VAFandCN <- cbind(data.frame(gr.snvs)[gr.overlap$queryHits,],
                     data.frame(gr.cnvs)[gr.overlap$subjectHits,])
Store_VAF_CN<-gr.VAFandCN  #run it for the first gr.VAFandCN

######################### data store for each sample #########################

for (i in 2:length(uni_patient)) {
  patient=uni_patient[i]
  snvs_dormant.selected <- snvs_dormant[which(snvs_dormant$Case==patient),]
  df.cnv.selected <- df.cnv[which(df.cnv$Case == patient),]
  
  #### Employ the genomic ranges package to intersect genomic intervals:
  gr.snvs <- makeGRangesFromDataFrame(snvs_dormant.selected, keep.extra.columns = TRUE,
                                      seqnames.field = c("Chromosome"), 
                                      start.field = "Start_Position", 
                                      end.field = "End_Position")
  # build a genomic ranges object for CN:
  gr.cnvs <- makeGRangesFromDataFrame(df.cnv.selected, keep.extra.columns = TRUE,
                                      seqnames.field = c("Chromosome"), 
                                      start.field = "Start", 
                                      end.field = "End")
  # intersect the two:
  gr.overlap <- data.frame(findOverlaps(gr.snvs, gr.cnvs))
  
  # merge the relevant columns from the intersection:
  gr.VAFandCN <- cbind(data.frame(gr.snvs)[gr.overlap$queryHits,],
                       data.frame(gr.cnvs)[gr.overlap$subjectHits,])
  
  Store_VAF_CN<-rbind(Store_VAF_CN,gr.VAFandCN)
}


Store_VAF_CN$VAF<-Store_VAF_CN$t_alt_count/(Store_VAF_CN$t_alt_count+Store_VAF_CN$t_ref_count)
Store_VAF_CN$CCF<-(2+Store_VAF_CN$purity*(Store_VAF_CN$CopyNumber-2))*Store_VAF_CN$VAF/Store_VAF_CN$purity
Store_VAF_CN<-Store_VAF_CN[!is.na(Store_VAF_CN$CCF),]
Store_VAF_CN[Store_VAF_CN$CCF>1,]$CCF<-1
rownames(Store_VAF_CN) <- NULL
Store_VAF_CN<-as.data.frame(Store_VAF_CN)
CCF_dormant<-Store_VAF_CN[Store_VAF_CN$Group=="Highly Quiescent",]
dormant_plot= subset(CCF_dormant, select = -c(Case,end,seqnames,start,strand,width))
CCF_cycling<-Store_VAF_CN[Store_VAF_CN$Group=="Fast Cycling",]
cycling_plot= subset(CCF_cycling, select = -c(Case,end,seqnames,start,strand,width))
all_plot=rbind(dormant_plot,cycling_plot)
Gaussian_dormant<-dormant_plot[dormant_plot$CCF<1,]
Gaussian_cycling<-cycling_plot[cycling_plot$CCF<1,]


###visualisation for CCF
ggplot(dormant_plot, aes(CCF,group=Tumor_Sample_Barcode,fill=Tumor_Sample_Barcode)) +
  geom_density(adjust=1.5,alpha=0.6)+
  theme_ipsum()+ggtitle("Density of CCF in dormant tumours in LUSC")+ theme(legend.position = "none")

ggplot(cycling_plot, aes(CCF, group=Tumor_Sample_Barcode,fill=Tumor_Sample_Barcode)) +
  geom_density(adjust=1.5,alpha=0.6)+
  theme_ipsum()+ggtitle("Density of CCF in fast cycling tumours in LUSC")+ theme(legend.position = "none")

ggplot(all_plot, aes(CCF, group=Group,fill=Group)) +
  geom_density(adjust=1.5,alpha=0.3)+
  theme_ipsum()+ggtitle("Density of CCF in different types of tumours in LUSC")
