### This script downloads and analyses mutation data for one cancer type from the TCGA project.

if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("maftools")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(maftools)
library(dplyr)
library(reshape)
library(ggpubr)
library(ggplot2)

##There are totally 31 cancer types, see details in Table 1 of the dissertation: 
#"ACC"   "BLCA"  "BRCA"     "CESC"   "CHOL" "COAD" "ESCA" 
#"GBM"   "HNSC"   "KICH"     "KIRC"   "KIRP"    "LGG"    "LIHC" "LUAD"
#"LUSC”  "MESO"   "OV"       "PAAD"  "PCPG" "PRAD" "READ" "SARC" 
#"SKCM"  "STAD"  "TGCT"    "THCA" "THYM" "UCEC"  "UCS"   "UVM"

#Here we use ESCA as an example##
?GDCquery_Maf
maf <- GDCquery_Maf("ESCA", pipelines = "mutect2") %>% read.maf
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)

## Get the number of substitution counts per sample:
substitutionStats = titv(maf = maf, plot = FALSE, useSyn = TRUE)
?titv

#Plot transitions/transversion graphs:
pdf("titv.ESCA.pdf")
plotTiTv(res = substitutionStats)
dev.off()

# Convert the substitution object to data frame in preparation for next step:
df.subs1 <- substitutionStats$fraction.contribution
df.subs2 <- substitutionStats$raw.counts # not informative for our analysis
df.subs3 <- substitutionStats$TiTv.fractions

df.subs.merged <- merge(df.subs1,df.subs3,
                         by.x="Tumor_Sample_Barcode", by.y="Tumor_Sample_Barcode",
                        all.x=FALSE, all.y=FALSE)


## Next, read the pan-cancer quiescence scores:
quiescence <- read.csv("PancancerQuiescenceGroups.csv",
                       header=TRUE, stringsAsFactors = FALSE)
quiescence$PatientID <- sapply(quiescence$Barcode,
                              function(x) paste(strsplit(x,"-")[[1]][1:4],collapse="-"))
df.subs.merged$PatientID <- sapply(as.character(df.subs.merged$Tumor_Sample_Barcode),
                           function(x) paste(strsplit(x,"-")[[1]][1:4],collapse="-"))

## Merge quiescence data frame with the substitutions data frame:
df.merged <- merge(quiescence, df.subs.merged,
                   by.x="PatientID", by.y="PatientID",
                   all.x=FALSE, all.y=FALSE)

##Visualisation
df.merged$`C>A` <- as.numeric(df.merged$`C>A`)
my_comparisons <- list(c("Fast Cycling","Highly Quiescent"))
ggboxplot(df.merged, x = "Group", y = "C>A",
          color = "Group", 
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)

df.melt <- melt(df.merged, id.vars=c("PatientID","CancerType",
                                     "Barcode","QuiescenceScore",
                                     "Group","Tumor_Sample_Barcode"))
df.melt$value <- as.numeric(df.melt$value)


pdf("substitutionsCompared.ESCA.pdf",w=10,h=8)
ggboxplot(df.melt, x = "Group", y = "value",
          color = "Group", 
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+ # Add pairwise comparisons p-value
  facet_wrap(~variable, scales = "free",nrow=2)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())+
  ylab("% substitutions")+
  xlab("")
dev.off()

