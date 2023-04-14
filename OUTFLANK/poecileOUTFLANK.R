#install.packages("devtools")
#library(devtools)
#if(!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("qvalue")
#install_github("whitlock/OutFLANK")
#install.packages("vcfR")
#if(!("bigstatsr" %in% installed.packages())){install.packages("bigstatsr")}
#if(!("bigsnpr" %in% installed.packages())){devtools::install_github("privefl/bigsnpr")}
#library(bigstatsr)
#library(bigsnpr)   # package for SNP trimming
#library(OutFLANK)
#library(vcfR)
#library(stringr)

############
# ORIGINAL #
############

#Define Populations
#sample_pop <- read.csv("/data5/K_Feldmann_data/sample_pop.txt", header = TRUE)
#poplist <- sample_pop$Population

###### Convert VCF to .geno-esque file using this code:
#obj.vcfR <- read.vcfR("/data5/K_Feldmann_data/filtered_snps_3_0.8_0.05.vcf")

#geno <- extract.gt(obj.vcfR) # Character matrix containing the genotypes
#positions <- getPOS(obj.vcfR) # Positions in bp
#chromosome_unedited <- getCHROM(obj.vcfR) # Chromosome information
#locinames <- paste(chromosome_unedited, positions, sep = ".")
#str_chromosome <- data.frame(sapply(str_split(chromosome_unedited, "-"), `[`, 3))
#colnames(str_chromosome) <- "Chromosome"
#numeric_chromosome <- read.csv(file = "/data5/K_Feldmann_data/numeric_chromosome.csv", header = T)
#chrom_data <- merge(str_chromosome, numeric_chromosome, by = "Chromosome")
#chromosomes <- chrom_data[,"Numeric_Chromosome"]

#G <- matrix(NA, nrow = nrow(geno), ncol = ncol(geno))
#G[geno %in% c("0/0", "0|0")] <- 0
#G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
#G[geno %in% c("1/1", "1|1")] <- 2
#G[geno %in% c("./.", ".|.")] <- 9
#G[is.na(G)] <- 9

#table(as.vector(G))

#save(G, file = "G.RData")
#load("G.RData")

#Trim SNPs
#G_na <- G
#G_na[G_na == 9] <- NA
#G_na <- cbind(G_na,positions,chromosomes)
#G_na <- G_na[complete.cases(G_na),]
#G_na <- G_na[order(G_na[,"chromosomes"],G_na[,"positions"]),]

#poecile <- list()
#poecile$G <- G_na[,!colnames(G_na) %in% c("positions","chromosomes")]
#poecile$positions <- G_na[,"positions"]
#poecile$chromosomes <- G_na[,"chromosomes"]

#G <- add_code256(big_copy(t(poecile$G), type = "raw"), code = bigsnpr:::CODE_012)
#newpc <- snp_autoSVD(G = G, infos.chr = poecile$chromosomes, infos.pos = poecile$positions, size = ncol(G), roll.size = 0)
#which_pruned <- attr(newpc, which = "subset") # Indexes of remaining SNPS after pruning
#save(which_pruned, file = "which_pruned.RData")
#load("which_pruned.RData")

#poecileOutFLANK <- MakeDiploidFSTMat(SNPmat = t(G), locusNames = locinames, popNames = poplist)

#Eliminate rows with NA
#poecileOutFLANK.edited <- poecileOutFLANK[complete.cases(poecileOutFLANK),]

#write.csv(poecileOutFLANK.edited, file = "/data5/K_Feldmann_data/OUTFLANK/outflankInput.csv")

#poecileOutFLANK.edited <- read.csv("/data5/K_Feldmann_data/OUTFLANK/outflankInput.csv", header = T)
#filtered_LD <- data.frame(read.table("/data5/K_Feldmann_data/filtered_snps_LD.prune.in", header = F))
#colnames(filtered_LD) <- "LocusName"
#poecileOutFLANK.edited.LD <- merge(filtered_LD, poecileOutFLANK.edited, by = "LocusName")

#pdf("outflankInput.pdf")
#plot(poecileOutFLANK.edited$He, poecileOutFLANK.edited$FST)

#plot(poecileOutFLANK.edited$FST, poecileOutFLANK.edited$FSTNoCorr, xlim = c(-0.1, 1), ylim = c(-0.1, 1), pch = 20)
#abline(0, 1) # Checking the effect of sample size on Fst since FSTNoCorr will be used in the follow

#plot(poecileOutFLANK.edited$He, poecileOutFLANK.edited$FSTNoCorr, pch=20, col="grey")

#hist(poecileOutFLANK.edited$FSTNoCorr)
#dev.off()

# Call Outliers

#outlierResults <- OutFLANK(poecileOutFLANK.edited.LD, NumberOfSamples = 8, RightTrimFraction = 0.35, qthreshold = 0.05, Hmin = 0.1)
#P1 <- pOutlierFinderChiSqNoCorr(poecileOutFLANK.edited, Fstbar = outlierResults$FSTNoCorrbar, dfInferred = outlierResults$dfInferred, qthreshold = 0.05, Hmin = 0.1)
#my_out <- P1$OutlierFlag==TRUE

outlierResults <- read.csv("/data5/K_Feldmann_data/OUTFLANK/outflankResults.csv", header = T)
outlierFST <- outlierResults[which(outlierResults$results.OutlierFlag == TRUE),"results.FST"]

# Plot the ditribution of Fst with the chi squared distribution
pdf("outflankResults.pdf")
#OutFLANKResultsPlotter(outlierResults, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.005, Zoom = FALSE, titletext = NULL)
#OutFLANKResultsPlotter(outlierResults, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.005, Zoom = TRUE, RightZoomFraction = 0.1, titletext = NULL)
#hist(outlierResults$results$pvaluesRightTail)
#plot(P1$He, P1$FST, pch=19, col=rgb(0,0,0,0.1))
#points(P1$He[my_out], P1$FST[my_out], col="blue")
#hist(P1$pvaluesRightTail)
#plot(P1$LocusName[P1$He>0.1], P1$FST[P1$He>0.1], xlab="Position", ylab="FST", col=rgb(0,0,0,0.2))
#points(P1$LocusName[my_out], P1$FST[my_out], col="magenta", pch=20)
hist(outlierFST)
dev.off()

# Identify Outliers
#outlierTable <- outlierResults$results$LocusName[outlierResults$results$OutlierFlag == TRUE]

#write.table(outlierTable, file = "/data5/K_Feldmann_data/OUTFLANK/outlierTable.txt", sep = "\t")
#write.csv(outlierResults, file = "/data5/K_Feldmann_data/OUTFLANK/outflankResults.csv")

###############
# ALTERNATIVE #
###############

#print("read genotypes file")
#genotypes <- read.table("poecileOUTFLANK.012")
#print("read loci file")
#loci <- read.table("poecileOUTFLANK.pos")
#print("read populations file")
#populations <- read.table("poecileOUTFLANK.pop")

#print("create OUTFLANK input file")
#outflank_input <- MakeDiploidFSTMat(SNPmat = genotypes, locusNames = loci, popNames = populations)
#outflank_input.edited <- outflank_input[complete.cases(outflank_input),]
#write.csv(outflank_input.edited, file = "/data5/K_Feldmann_data/outflankInput_alternative.csv")

#outlierResults <- OutFLANK(outflank_input.edited, NumberOfSamples = 8, RightTrimFraction = 0.05, qthreshold = 0.05, Hmin = 0.1)

#pdf("outflankResults_alternative.pdf")
#OutFLANKResultsPlotter(outlierResults, withOutliers = TRUE, NoCorr = TRUE, Hmin = 0.1, binwidth = 0.005, Zoom = FALSE, RightZoomFraction = 0.05, titletext = NULL)
#dev.off()

#outlierTable <- outlierResults$results$LocusName[outlierResults$results$OutlierFlag == TRUE]
