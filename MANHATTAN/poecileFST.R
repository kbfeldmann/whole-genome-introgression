#library(devtools)
#install_github("spflanagan/gwscaR")
#library(gwscaR)
#install.packages("qqman")
library(qqman)

data <- read.table("Fst_OUTFLANK_plot.txt", header = T)
#data <- read.table("Fst_VCFTOOLS_plot.txt", header = T)

data[!complete.cases(data),"FST"] <- 0

print("Creating Plot")
SNP <- c(1:(nrow(data)))
mydf <- data.frame(SNP,data)
write.table(mydf, file = "Fst_OUTFLANK_edited.txt")
pdf(file = "manhattanOUTFLANK.pdf")
plot <- manhattan(mydf, chr = "CHROM", bp = "POS", p = "FST", snp = "SNP", col = c("darkblue","orange3"), cex = 0.5, logp = FALSE, ylab = "Weir and Cockerham Fst")
dev.off()

