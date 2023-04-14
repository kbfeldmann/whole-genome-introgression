##Modified from Kathryn Grabenstein
##Plot whole genome data for BCCH/CACH Cognition Data

#install.packages("devtools")
#install.packages("tidyverse")
#install.packages("adegenet")
#install.packages("ggplot2")
library(devtools)
library(tidyverse)
library(adegenet)
library(ggplot2)
#install_github("zhengxwen/gdsfmt")
#if(!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("SNPRelate")
library(gdsfmt)
library(SNPRelate)

##PCA for Whole Genome
#######################

#INPUT VCF TABLE
print("input VCF table")
#vcf.fn <- paste("/data5/K_Feldmann_data/filtered_snps_3_0.8_0.05.vcf")
vcf.fn <- paste("/data5/K_Feldmann_data/PCA/filtered_snps_3_0.8_0.05_DROPPED_SIB.vcf")

#REFORMAT VCF
print("reformat VCF")
#snpgdsVCF2GDS(vcf.fn, "/data5/K_Feldmann_data/filtered_snps_3_0.8_0.05.gds", method="biallelic.only")
snpgdsVCF2GDS(vcf.fn, "/data5/K_Feldmann_data/PCA/filtered_snps_3_0.8_0.05_DROPPED_SIB.gds", method="biallelic.only")

#SUMMARY OF VCF
#snpgdsSummary("/data5/K_Feldmann_data/filtered_snps_3_0.8_0.05.gds")
#snpgdsSummary("/data5/K_Feldmann_data/PCA/filtered_snps_3_0.8_0.05_DROPPED.gds")

#OPEN THE GDS FILE
#print("open the GDS file")
#genofile <- snpgdsOpen("/data5/K_Feldmann_data/filtered_snps_3_0.8_0.05.gds")
#genofile <- snpgdsOpen("/data5/K_Feldmann_data/PCA/filtered_snps_3_0.8_0.05_DROPPED.gds")

#RUN PCA
#print("Run pca")
#pca <- snpgdsPCA(genofile, num.thread=4, autosome.only = FALSE)

#CALCULATE: percent of variation that is accounted for by top PCA components:
#print("calculate percent of variation that is accounted for by top PCA components")
#pc.percent <- pca$varprop*100
#head(round(pc.percent, 2))

#SAMPLE ID AND POPULATION INFORMATION
#print("get sample ID")
#sample.id.pop <- read.delim("/data5/K_Feldmann_data/sample_pop.txt", sep = ",", header = TRUE)
#sample.id.pop <- read.delim("/data5/K_Feldmann_data/PCA/sample_pop_DROPPED.txt", sep = ",", header = TRUE)
#chickadeeInfo <- read.table(pipe("pbpaste"), sep = "\t", header = T, na.strings = c("","NA"))
#nrow(sample.id.pop)

#print("make new dataframe with all variables for plotting")
#tab <- data.frame(Sample_ID = sample.id.pop$Sample_ID,
                  #pop = sample.id.pop$Population,
                  #EV1 = pca$eigenvect[,1],    # the first eigenvector
                  #EV2 = pca$eigenvect[,2],    # the second eigenvector
                  #stringsAsFactors = FALSE)

#tab <- merge(tab, chickadeeInfo, by = "Sample_ID")
#tab$Location <- factor(tab$Location, levels = c("Ithaca, New York","DeRuyter, New York","Hickory Run State Park","Lehigh campus","Jacobsburg State Park","DeSales University","East Carolina West Research Campus","LSU"))
#head(tab)

#PLOT PCA
#print("plot in ggplot")
#PCA.Whole.Genome <- ggplot(tab, aes(x=EV1, y=EV2)) + geom_point(aes(color=Location), size=3, position = "jitter") +
#geom_text(aes(label=Sample_ID), position = position_jitter(width = 0.02, height = 0.02)) +
  #theme_bw() +
  #theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#PCA.Whole.Genome + scale_color_manual(values=c("dodgerblue4","dodgerblue2","dodgerblue","darkorchid4","darkorchid3","darkorchid1","firebrick4","firebrick1", "darkgrey"))
#PCA.Whole.Genome + scale_shape_manual(c(21,21,21,23,23,23,22,22,8))
#PCA.Whole.Genome + scale_color_manual(values=c("#909BFF", "#FFFB00", "#a134eb"))

#print("saving as pdf")
#pdf("filtered_snps_3_0.8_0.05_DROPPED.pdf")
#print(PCA.Whole.Genome)
#dev.off()

#########################



