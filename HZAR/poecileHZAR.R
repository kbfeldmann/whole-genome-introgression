library(tidyr)
library(stringr)
library(maditr)

## CREATE INPUT FILE ##

#LSU <- read.table(file = "LSUfrequencies.frq", sep = "\t", header = T, row.names = NULL)
#LSU$Location <- "LSU"
#ECW <- read.table(file = "ECWfrequencies.frq", sep = "\t", header = T, row.names = NULL)
#ECW$Location <- "ECW"
#DSU <- read.table(file = "DSUfrequencies.frq", sep = "\t", header = T, row.names = NULL)
#DSU$Location <- "DSU"
#LHU <- read.table(file = "LHUfrequencies.frq", sep = "\t", header = T, row.names = NULL)
#LHU$Location <- "LHU"
#JSP <- read.table(file = "JSPfrequencies.frq", sep = "\t", header = T, row.names = NULL)
#JSP$Location <- "JSP"
#HRP <- read.table(file = "HRPfrequencies.frq", sep = "\t", header = T, row.names = NULL)
#HRP$Location <- "HRP"
#INY <- read.table(file = "INYfrequencies.frq", sep = "\t", header = T, row.names = NULL)
#INY$Location <- "INY"
#DNY <- read.table(file = "DNYfrequencies.frq", sep = "\t", header = T, row.names = NULL)
#DNY$Location <- "DNY"

#print("combine input data")
#LCTN <- rbind(LSU, ECW, DSU, LHU, JSP, HRP, INY, DNY)

#write.table(LCTN, file = "combinedHZAR.txt")

#LCTN <- read.table(file = "combinedHZAR.txt", header = T)
#colnames(LCTN) <- c("Scaffold","Position","NumAlleles","NumChr","AlleleFreq_1","AlleleFreq_2","Location")

#print("Entire Dataframe Processing")
#LCTN <- LCTN %>% separate(AlleleFreq_1, sep = ":", c("Allele1","AlleleFrequency1"))
#LCTN <- LCTN %>% separate(AlleleFreq_2, sep = ":", c("Allele2","AlleleFrequency2"))
#LCTN$Scaffold <- paste(LCTN$Scaffold, LCTN$Position, sep = ".")
#LCTN$AlleleFrequency1 <- as.numeric(LCTN$AlleleFrequency1)
#LCTN$AlleleFrequency2 <- as.numeric(LCTN$AlleleFrequency2)

#print("Split Data to Restructure")
#alleleInfo <- LCTN[,c("Scaffold","NumChr","Location")]
#allele1LCTN <- LCTN[,c("Scaffold","Allele1","AlleleFrequency1","Location")]
#allele2LCTN <- LCTN[,c("Scaffold","Allele2","AlleleFrequency2","Location")]

#alleleInfo$Scaffold <- paste(alleleInfo$Scaffold, "samples",sep = ".")
#allele1LCTN$Scaffold <- paste(LCTN$Scaffold, LCTN$Allele1, sep = ".")
#allele1LCTN$Allele1 <- NULL
#allele2LCTN$Scaffold <- paste(LCTN$Scaffold, LCTN$Allele2, sep = ".")
#allele2LCTN$Allele2 <- NULL

#print("Restructure Data")
#alleleInfo_dcast <- alleleInfo %>% dcast(Location ~ Scaffold, value.var = "NumChr")
#allele1LCTN_dcast <- allele1LCTN %>% dcast(Location ~ Scaffold, value.var = "AlleleFrequency1")
#allele2LCTN_dcast <- allele2LCTN %>% dcast(Location ~ Scaffold, value.var = "AlleleFrequency2")
#HZARinput <- merge(allele1LCTN_dcast, allele2LCTN_dcast, by = "Location")
#HZARinput <- merge(HZARinput, alleleInfo_dcast, by = "Location")

#print("Add Distance Data")
#distances <- data.frame(Location = c("LSU","ECW","DSU","LHU","JSP","HRP","INY","DNY"), Distance_km = c(0,580,1125,1132,1154,1180,1337,1372))
#HZARinput <- merge(HZARinput, distances, by = "Location")
#write.csv(HZARinput, file = "inputHZAR.csv")

## SELECT LOCI TO PLOT ##
num_SNPS <- 10

print("Read OutFLANK Data")
outflank_results <- read.csv("/data5/K_Feldmann_data/outflankResults.csv", header = T)
pvalue_results <- outflank_results[order(outflank_results$results.pvalues, na.last = T),c("results.LocusName","results.pvalues")]
pvalue_subset <- pvalue_results[1:num_SNPS,]

print("Read HZAR Input Data")
HZARinput <- read.csv("/data5/K_Feldmann_data/HZAR/inputHZAR.csv", header = T)
save(HZARinput, file = "inputHZAR.RData")
#load("inputHZAR.RData")

print("Isolate Signficant Loci in HZAR Data")
plotHZAR <- HZARinput[,c("Location","distance",grep(pattern = pvalue_subset$results.LocusName, colnames(HZARinput)))]
print(head(plotHZAR))

