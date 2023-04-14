#install.packages("devtools")
#library(devtools)
#devtools::install_github("ribailey/gghybrid")
library("gghybrid")

# Modified 'plot_clinecurve' function to reduce line widths (lwd): https://rdrr.io/github/ribailey/gghybrid/src/R/plot_clinecurve.r
#------------------------------------------------------------------------------------#
plot_clinecurve = function(ggcline.object,data.prep.object,esth.object,test.subject = "INDLABEL",include.Source = FALSE,cline.locus,locus.column="locus",cline.col="black",null.line.locus=NULL,null.line.col="black",plot.data=NULL,data.col=NULL,plot.genotype=FALSE,PLOIDY,cline.centre.line=NULL,cline.centre.col="black",...){

    if(!include.Source){
        setkey(data.prep.object,Source);
        data.prep.object = data.prep.object[Source == "TEST"];
        esth.object = esth.object[Source == "TEST"];
        };

    setkeyv(data.prep.object,test.subject);setkeyv(esth.object,test.subject);
    prep.hi = data.prep.object[esth.object];

    setkeyv(prep.hi,locus.column);
    setkeyv(ggcline.object,locus.column);

    v = numeric(1);
    u = numeric(1);
    S1.prop_1 = numeric(1);
    S0.prop_1 = numeric(1);

    plot(prep.hi[cline.locus[1],Source_allele]~prep.hi[cline.locus[1],qlogis(h_posterior_mode)],type="n",xlim=c(0,1),ann=FALSE,...);

    if(is.null(null.line.locus)==FALSE){
		if(null.line.col=="black"){null.line.col=rep("black",length(null.line.locus))};
        for(i in 1:length(null.line.locus)){
            S1.prop_1 = unique(prep.hi[null.line.locus[i],S1.prop_1]);
            S0.prop_1 = unique(prep.hi[null.line.locus[i],S0.prop_1]);
            abline(a=S0.prop_1,b=(S1.prop_1 - S0.prop_1),lty=2,lwd=2,col=null.line.col[i]);
            };
    };

    if(is.null(cline.centre.line)==FALSE){
        centre=ggcline.object[cline.centre.line,centre_mean];
        abline(v=centre,col=cline.centre.col,lty=2,lwd=2);
        setkey(ggcline.object,NULL);
    };

    setkeyv(ggcline.object,locus.column);
    if(cline.col=="black"){cline.col=rep("black",length(cline.locus))};
    for(i in 1:length(cline.locus)){
        v = ggcline.object[cline.locus[i],v_mean];
        u = qlogis(ggcline.object[cline.locus[i],centre_mean])*ggcline.object[cline.locus[i],v_mean];
        S1.prop_1 = unique(prep.hi[cline.locus[i],S1.prop_1]);
        S0.prop_1 = unique(prep.hi[cline.locus[i],S0.prop_1]);

        par(new=T);

        curve(S0.prop_1 + (x^v/(x^v + (1 - x)^v*exp(u)))*(S1.prop_1 - S0.prop_1),from=0,to=1,axes=F,xlab="",ylab="", col=cline.col[i],lwd=0.5,ylim=c(0,1)); # MODIFIED

        if(i == length(cline.locus)){par(new=F)};
        };

    if(is.null(plot.data)==FALSE){

        setnames(prep.hi,locus.column,"loc");

        if(plot.genotype==TRUE){
            setnames(prep.hi,test.subject,"INDLABEL");
            prep.hi.geno= melt(dcast(prep.hi, INDLABEL + h_posterior_mode~loc,value.var="Source_allele",fun=sum),id = c(test.subject,"h_posterior_mode"),variable.name = locus.column,value.name = "genotype");
            prep.hi.geno[,genotype:=genotype/PLOIDY];
		
            setkeyv(prep.hi.geno,locus.column);###was setkey(prep.hi.geno,loc)###

            for(i in 1:length(plot.data)){
                points(prep.hi.geno[cline.locus[i],genotype]~prep.hi.geno[cline.locus[i],h_posterior_mode],col=data.col[i],...);
            };
										 
            setnames(prep.hi,"INDLABEL",test.subject);									 
        }else{
            setkey(prep.hi,loc);
            for(i in 1:length(plot.data)){
                points(prep.hi[cline.locus[i],Source_allele]~prep.hi[cline.locus[i],h_posterior_mode],col=data.col[i],...);
            };
        };
    };
}

plot_h = function(data,subset.by=NULL,test.subject="INDLABEL",POPID.name="POPID",
    mean.h.by=NULL,sort.by="h_posterior_mode",
    col.group=NULL,group.sep=NULL,fill.source=FALSE,
    basic.lines=TRUE,source.col=NULL,source.limits=NULL,
    custom.abline=NULL,...){

    setkey(data,NULL);

    if(is.null(col.group)==FALSE){
	col.Dark2 = c("#666666","orchid1","#E6AB02","#66A61E","#E7298A","darkorchid4","darkorchid3","#1B9E77")
        #col.Dark2 = c("navy","dodgerblue4","dodgerblue","darkorchid4","darkorchid3","orchid1","firebrick4","firebrick1")
        col.Dark2 = data.table(col.Dark2);	

        setkeyv(data,col.group);
        pop = data.table(unique(data[,col.group,with=F]));
        pop.col = cbind(pop,col.Dark2);

        setkeyv(pop.col,col.group);
        data=data[pop.col];

        setkey(data,NULL);setkey(data,col.Dark2);
        col.summary = data[,head(.SD,1),by=col.Dark2]
                                 };

    setkey(data,NULL);

    if(is.na(chmatch("Source",sort.by))==FALSE){
        setkey(data,Source);
        data["S0",Source2:="A"];
        data["TEST",Source2:="B"];
        data["S1",Source2:="C"];
        sort.by[chmatch("Source",sort.by)]="Source2";
                                               };

    setkey(data,NULL);
    setkeyv(data,test.subject);

    if(is.null(mean.h.by)){
        setkeyv(data,sort.by);
        data[,rn:=row(data)[,1]];
		                  }else{
        setkeyv(data,mean.h.by);
        data[,mean_h:=mean(h_posterior_mode),by=mean.h.by];
        setkeyv(data,sort.by);
        data[,rn:=row(data)[,1]];
                               };

    setkey(data,NULL);

    plot(data[,h_posterior_mode]~data[,rn],type="n",
        main="Hybrid index and 95% C.I.",
        xlab="Individual",ylab="Hybrid index",...);

    if(is.null(group.sep)==FALSE){
        setkeyv(data,group.sep);
        abline(v=0.5,col="grey");
        abline(v=data[,(max(rn)+0.5),by=group.sep]$V1,col="grey");
                                 };

    setkey(data,NULL);

	setnames(data,POPID.name,"POPID");

    if(fill.source==TRUE){
        setkey(data,Source);

    for(i in 1:(data[Source=="S0",length(unique(POPID))])){
        rect(data[Source=="S0",(min(rn)),by=POPID]$V1[i] - 0.5,
            par("usr")[3],
            data[Source=="S0",(max(rn)),by=POPID]$V1[i] + 0.5,
            par("usr")[4], col = "grey95");
                                                          };

    for(i in 1:(data[Source=="S1",length(unique(POPID))])){
        rect(data[Source=="S1",(min(rn)),by=POPID]$V1[i] - 0.5,
            par("usr")[3],
            data[Source=="S1",(max(rn)),by=POPID]$V1[i] + 0.5,
            par("usr")[4], col = "grey95");
                                                          };
                         };

    setkey(data,NULL);

    if(basic.lines==TRUE){
        abline(h=c(0,0.5,1),col="grey",lty=2,lwd=2);
                         };

    if(is.null(source.col)==FALSE){
        setkey(data,Source);

        if(is.null(col.group)){
            data[,col.Dark2:="black"];
                              };
        data["S0",col.Dark2:=source.col[1]];
        data["S1",col.Dark2:=source.col[2]];

        setkey(data,NULL);setkey(data,col.Dark2);
        col.summary = data[,head(.SD,1),by=col.Dark2];
        setkey(col.summary,Source);
        col.summary["S0",POPID:="S0"];
        col.summary["S1",POPID:="S1"];
                                  };

    setkey(data,NULL);

    if(is.null(col.group) & is.null(source.col)){
        points(data[,h_posterior_mode]~data[,rn],...);
        arrows(data[,
            rn],data[,h_cred_int_lower],data[,
            rn],data[,h_cred_int_upper],angle=90,code=3,length=0);
                                                }else{
        points(data[,h_posterior_mode]~data[,rn],
            col=data[,col.Dark2],...);
        arrows(data[,
            rn],data[,h_cred_int_lower],data[,
            rn],data[,h_cred_int_upper],angle=90,code=3,length=0,
            col=data[,col.Dark2]);
                                                     };

    if(is.null(source.limits)==FALSE){
        abline(h=data[Source=="S0",max(h_cred_int_upper)],
            col=source.limits[1],lty=2,lwd=2);
        abline(h=data[Source=="S1",min(h_cred_int_lower)],
            col=source.limits[2],lty=2,lwd=2);
                                     };

    if(is.null(custom.abline)==FALSE){
        custom.abline;
                                     };

	setnames(col.summary,"POPID",POPID.name);

    if(is.null(mean.h.by)){
        col.summary=col.summary[,c(test.subject,"Source","col.Dark2","rn"),with=F];
                          }else{
        col.summary=col.summary[,c(test.subject,mean.h.by,"Source","mean_h","col.Dark2","rn"),with=F];
                               };

	return(col.summary)
}
#------------------------------------------------------------------------------------#

#print("read data")
#ggdata <- read.data(file = "gghybridFst0.6.str", nprecol = 2, MISSINGVAL = -9)

#save(ggdata, file = "ggdata.RData")
load("/data5/K_Feldmann_data/GGHYBRID/OUTLIERS/ggdata.RData")

#print("data preparation and filtering")
#prepdata <- data.prep(data=ggdata$data, loci=ggdata$loci, alleles=ggdata$alleles, S0=c("HRP","INY","DNY"), S1=c("LSU","ECW"), precols=ggdata$precols, max.S.MAF = 0.1, return.genotype.table=T, return.locus.table=T)

#save(prepdata, file = "prepdata.RData")
load("/data5/K_Feldmann_data/GGHYBRID/OUTLIERS/prepdata.RData")

#hindlabel <- esth(data.prep.object = prepdata$data.prep, read.data.precols = ggdata$precols, include.Source = TRUE, plot.col = c("dodgerblue4","darkorchid4","darkorchid3","orchid1","firebrick3"), nitt=10000, burnin=5000)

#save(hindlabel, file = "hindlabel_new.RData")
load("/data5/K_Feldmann_data/GGHYBRID/OUTLIERS/hindlabel.RData")

#setkey(hindlabel$hi,POPID)

# Order populations by latitudinal distance
#hindlabel$hi$OrderLetter <- hindlabel$hi$POPID
#setattr(hindlabel$hi$OrderLetter,"levels",c(levels(hindlabel$hi$OrderLetter),"A","B","C","D","E","F","G","H"))
#hindlabel$hi$OrderLetter[which(hindlabel$hi$OrderLetter == "DNY")] <- "A"
#hindlabel$hi$OrderLetter[which(hindlabel$hi$OrderLetter == "INY")] <- "B"
#hindlabel$hi$OrderLetter[which(hindlabel$hi$OrderLetter == "HRP")] <- "C"
#hindlabel$hi$OrderLetter[which(hindlabel$hi$OrderLetter == "JSP")] <- "D"
#hindlabel$hi$OrderLetter[which(hindlabel$hi$OrderLetter == "LHU")] <- "E"
#hindlabel$hi$OrderLetter[which(hindlabel$hi$OrderLetter == "DSU")] <- "F"
#hindlabel$hi$OrderLetter[which(hindlabel$hi$OrderLetter == "ECW")] <- "G"
#hindlabel$hi$OrderLetter[which(hindlabel$hi$OrderLetter == "LSU")] <- "H"

# Plot hybrid index by individual
#pdf("hybridIndex.pdf")
#abc <- plot_h(data=hindlabel$hi[c("LSU","ECW","DSU","LHU","JSP","HRP","INY","DNY")],test.subject=hindlabel$test.subject,mean.h.by="POPID",sort.by=c("OrderLetter","h_posterior_mode"),col.group="POPID",group.sep="POPID",fill.source=TRUE,basic.lines=FALSE,source.col=c("dodgerblue4","firebrick3"),cex=1,pch=16,cex.lab=1.5,cex.main=1.5,ylim=c(0,1),asp=30)
#setkey(abc,rn)
#legend("topleft",abc[,POPID],bg="white",text.col=c("black"),pch=22,col=abc[,col.Dark2],pt.bg=abc[,col.Dark2],ncol=2,cex=1, pt.cex=1)
#dev.off()

#print("run Bayesian MCMC")
#gc1 <- ggcline(data.prep.object=prepdata$data.prep,esth.object=hindlabel,read.data.precols=dat$precols,nitt = 10000, burnin = 5000, print.k = 50)

#save(gc1, file = "gc1.RData")
load("/data5/K_Feldmann_data/GGHYBRID/OUTLIERS/gc1.RData")

print(head(gc1))

# Plot histogram of hybrid index values
#pdf("hybridindexHistogram.pdf")
#hist(hindlabel$hi$h_posterior_mode, main = "", xlab = "Hybrid Index", ylab = "Number of Individuals")
#dev.off()


# PLOT GENOMIC CLINE #
#-------------------------------------------------------------------------#
num_SNPs <- 20600
#-------------------------------------------------------------------------#

#####################
# P-Value Heat Map #
#####################

# GGHYBRID #
#gghybrid_pvalue <- gc1$gc[,c("locus","centre_pvalue")]
#pvalue_table <- gghybrid_pvalue[order(gghybrid_pvalue$centre_pvalue, na.last = T),]
#colnames(pvalue_table) <- c("locus","pvalue")

# OutFLANK #
#positions <- read.table("positions_Fst0.6.txt") # For FST Data
#positions <- read.table("/data5/K_Feldmann_data/OUTFLANK/filtered_snps_outliers.txt") # For Outlier Data
#positions <- paste(positions$V1, positions$V2, sep = ".")
#snp_positions <- cbind(positions,ggdata$loci)
#colnames(snp_positions) <- c("results.LocusName","SNP_Number")

#outflank_data <- read.csv("/data5/K_Feldmann_data/OUTFLANK/outflankResults.csv", header  = T)
#outflank_pvalue <- outflank_data[,c("results.LocusName","results.pvalues")]
#outflank_pvalue <- merge(snp_positions, outflank_pvalue, by = "results.LocusName")
#outflank_pvalue <- outflank_pvalue[which(gc1$gc$locus %in% outflank_pvalue$SNP_Number),]
#pvalue_table <- outflank_pvalue[order(outflank_pvalue$results.pvalues, na.last = T),c("SNP_Number","results.pvalues")]
#colnames(pvalue_table) <- c("locus","pvalue")

## Plot ##

# Group SNPs into equal sized bins based on p-value
#lociA <- pvalue_table[1:(num_SNPs/5),]
#lociB <- pvalue_table[(num_SNPs/5+1):(2*num_SNPs/5),]
#lociC <- pvalue_table[(2*num_SNPs/5+1):(3*num_SNPs/5),]
#lociD <- pvalue_table[(3*num_SNPs/5+1):(4*num_SNPs/5),]
#lociE <- pvalue_table[(4*num_SNPs/5+1):num_SNPs,]

# Program output describing the range of p-values in each color bin
#print(paste("Total Number of SNPs Plotted:", num_SNPs, sep = " "))
#print(paste("RED:", min(lociA$pvalue), "-", max(lociA$pvalue), ":", nrow(lociA), "SNPs", sep = " "))
#print(paste("ORANGE:", min(lociB$pvalue), "-", max(lociB$pvalue), ":", nrow(lociB), "SNPs", sep = " "))
#print(paste("GREEN:", min(lociC$pvalue), "-", max(lociC$pvalue), ":", nrow(lociC), "SNPs", sep = " "))
#print(paste("BLUE:", min(lociD$pvalue), "-", max(lociD$pvalue), ":", nrow(lociD), "SNPs", sep = " "))
#print(paste("BLACK:", min(lociE$pvalue), "-", max(lociE$pvalue), ":", nrow(lociE), "SNPs", sep = " "))

# Create input for genomic cline plot
#loci <- c(lociE$locus,lociD$locus,lociC$locus,lociB$locus,lociA$locus)
#colors <- c(rep("gray0",length(lociE$locus)),rep("blue",length(lociD$locus)),rep("green",length(lociC$locus)),rep("orange",length(lociB$locus)),rep("red",length(lociA$locus)))

########################
# Color by Chromosome #
########################
library(stringr)

# GGHYBRID #
#positions <- read.table("positions_Fst0.6.txt") # For FST Data
positions <- read.table("/data5/K_Feldmann_data/OUTFLANK/filtered_snps_outliers.txt") # For Outlier Data
gghybrid_pvalue <- gc1$gc[,c("locus","centre_pvalue")]
pvalue_table <- gghybrid_pvalue[order(gghybrid_pvalue$centre_pvalue, na.last = T),]
colnames(pvalue_table) <- c("locus","pvalue")
pvalue_table <- pvalue_table[1:num_SNPs,]

# OutFLANK P-Value #
#positions <- read.table("positions_Fst0.6.txt") # For FST Data
#positions <- read.table("/data5/K_Feldmann_data/OUTFLANK/filtered_snps_outliers.txt") # For Outlier Data
#positions.combined <- paste(positions$V1, positions$V2, sep = ".")
#snp_positions <- cbind(positions.combined,ggdata$loci)
#colnames(snp_positions) <- c("results.LocusName","locus")

#P-value
#outflank_data <- read.csv("/data5/K_Feldmann_data/OUTFLANK/outflankResults.csv", header  = T)
#outflank_pvalue <- outflank_data[,c("results.LocusName","results.pvalues")]
#outflank_pvalue <- merge(snp_positions, outflank_pvalue, by = "results.LocusName")
#outflank_pvalue <- outflank_pvalue[which(gc1$gc$locus %in% outflank_pvalue$locus),]
#pvalue_table <- outflank_pvalue[order(outflank_pvalue$results.pvalues, na.last = T),]
#pvalue_table <- pvalue_table[1:num_SNPs,]

#FST
#outflank_data <- read.csv("/data5/K_Feldmann_data/OUTFLANK/outflankResults.csv", header  = T)
#outflank_fst <- outflank_data[,c("results.LocusName","results.FST")]
#outflank_fst <- merge(snp_positions, outflank_fst, by = "results.LocusName")
#outflank_fst <- outflank_fst[which(gc1$gc$locus %in% outflank_fst$locus),]
#pvalue_table <- outflank_fst[order(outflank_fst$results.FST, na.last = T, decreasing = T),]
#pvalue_table <- pvalue_table[1:num_SNPs,]

#He
#outflank_data <- read.csv("/data5/K_Feldmann_data/OUTFLANK/outflankResults.csv", header  = T)
#outflank_he <- outflank_data[,c("results.LocusName","results.He")]
#outflank_he <- merge(snp_positions, outflank_he, by = "results.LocusName")
#outflank_he <- outflank_he[which(gc1$gc$locus %in% outflank_he$locus),]
#pvalue_table <- outflank_he[order(outflank_he$results.He, na.last = T, decreasing = F),]
#pvalue_table <- pvalue_table[1:num_SNPs,]

#print(max(pvalue_table$results.FST))
#print(min(pvalue_table$results.FST))

# OutFLANK Outlier #
#positions <- read.table("/data5/K_Feldmann_data/positions_Fst0.6.txt")
#positions.combined <- paste(positions$V1, positions$V2, sep = ".")
#snp_positions <- cbind(positions.combined,ggdata$loci)
#colnames(snp_positions) <- c("results.LocusName","locus")

#outflank_data <- read.csv("/data5/K_Feldmann_data/OUTFLANK/outflankResults.csv", header  = T)
#outflank_pvalue <- outflank_data[,c("results.LocusName","results.pvalues")]
#outflank_pvalue <- merge(snp_positions, outflank_pvalue, by = "results.LocusName")
#outflank_pvalue <- outflank_pvalue[which(gc1$gc$locus %in% outflank_pvalue$locus),]
#colnames(outflank_pvalue) <- c("results.LocusName","locus","results.pvalues")

#print(nrow(outflank_pvalue))

#outliers <- read.table("/data5/K_Feldmann_data/OUTFLANK/outlierTable.txt", header = T)
#colnames(outliers) <- "results.LocusName"
#pvalue_table <- merge(outliers, outflank_pvalue, by = "results.LocusName")

## Plot ##

positions <- positions$V1
positions <- str_split(positions, "-")
chromosome <- sapply(positions, "[[", 3)
snp_positions <- cbind(chromosome,ggdata$loci)
colnames(snp_positions) <- c("Chromosome","locus")
chromData <- merge(pvalue_table, snp_positions, by = "locus")
chromData <- chromData[order(chromData$pvalue, decreasing = T, na.last = T),]

# Assign a color to every possible chromosome
chromData$color <- ""
chromData[which(chromData$Chromosome == "CHR_1"),"color"] <- "purple"
chromData[which(chromData$Chromosome == "CHR_1A"),"color"] <- "green"
chromData[which(chromData$Chromosome == "CHR_1B"),"color"] <- "cornsilk"
chromData[which(chromData$Chromosome == "CHR_2"),"color"] <- "firebrick4"
chromData[which(chromData$Chromosome == "CHR_3"),"color"] <- "deeppink"
chromData[which(chromData$Chromosome == "CHR_4"),"color"] <- "darkviolet"
chromData[which(chromData$Chromosome == "CHR_4A"),"color"] <- "burlywood4"
chromData[which(chromData$Chromosome == "CHR_5"),"color"] <- "deepskyblue"
chromData[which(chromData$Chromosome == "CHR_6"),"color"] <- "blue"
chromData[which(chromData$Chromosome == "CHR_7"),"color"] <- "cornflowerblue"
chromData[which(chromData$Chromosome == "CHR_8"),"color"] <- "coral"
chromData[which(chromData$Chromosome == "CHR_9"),"color"] <- "chocolate"
chromData[which(chromData$Chromosome == "CHR_10"),"color"] <- "darkcyan"
chromData[which(chromData$Chromosome == "CHR_11"),"color"] <- "gold"
chromData[which(chromData$Chromosome == "CHR_12"),"color"] <- "darkolivegreen1"
chromData[which(chromData$Chromosome == "CHR_13"),"color"] <- "chartreuse"
chromData[which(chromData$Chromosome == "CHR_14"),"color"] <- "darkorange"
chromData[which(chromData$Chromosome == "CHR_15"),"color"] <- "aquamarine4"
chromData[which(chromData$Chromosome == "CHR_17"),"color"] <- "blueviolet"
chromData[which(chromData$Chromosome == "CHR_18"),"color"] <- "darkkhaki"
chromData[which(chromData$Chromosome == "CHR_19"),"color"] <- "dodgerblue4"
chromData[which(chromData$Chromosome == "CHR_20"),"color"] <- "purple"
chromData[which(chromData$Chromosome == "CHR_21"),"color"] <- "darkslategray1"
chromData[which(chromData$Chromosome == "CHR_22"),"color"] <- "darkmagenta"
chromData[which(chromData$Chromosome == "CHR_23"),"color"] <- "deeppink4"
chromData[which(chromData$Chromosome == "CHR_24"),"color"] <- "goldenrod1"
chromData[which(chromData$Chromosome == "CHR_25"),"color"] <- "aquamarine"
chromData[which(chromData$Chromosome == "CHR_26"),"color"] <- "darksalmon"
chromData[which(chromData$Chromosome == "CHR_27"),"color"] <- "cyan"
chromData[which(chromData$Chromosome == "CHR_28"),"color"] <- "darkgreen"
chromData[which(chromData$Chromosome == "CHR_MT"),"color"] <- "black"
chromData[which(chromData$Chromosome == "CHR_Z"),"color"] <- "darkblue"
chromData[which(chromData$Chromosome == "CHR_UNK"),"color"] <- "black"
chromData[which(chromData$Chromosome == "LGE22"),"color"] <- "darkslateblue"

# Isolate by chromosome
#chromData <- chromData[which(chromData$Chromosome == "CHR_1" | chromData$Chromosome == "CHR_Z" | chromData$Chromosome == "CHR_1A" | chromData$Chromosome == "CHR_3" | chromData$Chromosome == "CHR_5" | chromData$Chromosome == "CHR_2"),]
chromData <- chromData[which(chromData$Chromosome == "CHR_MT"),]

# Output describing the chromosome-color combinations
print(unique(chromData[,c("Chromosome","color")]))	
print(table(chromData$Chromosome))
print(head(chromData$pvalue, n=100))

# Randomly mix lines so chromosomes are not plotted in layers - plotting in layers could give false information on what chromosomes are primarily being plotted
set.seed(42)
rows <- sample(nrow(chromData))
chromData <- chromData[rows,]

# Create input for genomic cline plot
loci <- as.character(chromData$locus)
colors <- as.character(chromData$color)

#########################
# GGHYBRID vs. OutFLANK #
#########################
#'%notin%' <- Negate('%in%')

# GGHYBRID #
#gghybrid_pvalue <- gc1$gc[,c("locus","centre_pvalue")]
#gghybrid_pvalue <- gghybrid_pvalue[order(gghybrid_pvalue$centre_pvalue),]
#gghybrid_subset <- gghybrid_pvalue[1:num_SNPs,]

# OutFLANK #
#positions <- read.table("positions_Fst0.6.txt")
#positions <- paste(positions$V1, positions$V2, sep = ".")
#snp_positions <- cbind(positions,ggdata$loci)
#colnames(snp_positions) <- c("results.LocusName","SNP_Number")

#outflank_data <- read.csv("outflankResults.csv", header  = T)
#outflank_pvalue <- outflank_data[,c("results.LocusName","results.pvalues")]
#outflank_pvalue <- merge(snp_positions, outflank_pvalue, by = "results.LocusName")
#outflank_pvalue$results.LocusName <- NULL
#outflank_pvalue <- outflank_pvalue[which(gc1$gc$locus %in% outflank_pvalue$SNP_Number),]

#colnames(outflank_pvalue) <- c("locus","results.pvalues")
#outflank_pvalue <- outflank_pvalue[order(outflank_pvalue$results.pvalues),]
#outflank_subset <- outflank_pvalue[1:num_SNPs,]

## Plot ##

#print(paste("GGHYBRID:", min(gghybrid_subset$centre_pvalue), "-", max(gghybrid_subset$centre_pvalue), sep = " "))
#print(paste("OutFLANK:", min(outflank_subset$results.pvalues), "-", max(outflank_subset$results.pvalues), sep = " "))

#combined_loci <- merge(outflank_subset, gghybrid_subset, by = "locus")
#gghybrid_loci <- gghybrid_subset[which(gghybrid_subset$locus %notin% combined_loci$locus),]
#outflank_loci <- outflank_subset[which(outflank_subset$locus %notin% combined_loci$locus),]

#colnames(gghybrid_loci) <- c("locus","pvalue")
#colnames(outflank_loci) <- c("locus","pvalue")
#unique_loci <- rbind(gghybrid_loci, outflank_loci)

#set.seed(42)
#rows <- sample(nrow(unique_loci))
#unique_loci <- unique_loci[rows,]

# Create input for genomic cline plot
#loci <- c(unique_loci$locus,combined_loci$locus)
#colors <- c(rep("grey",length(unique_loci$locus)),rep("black",length(combined_loci$locus)))

#####################
# OutFLANK Outliers #
#####################

#positions <- read.table("positions_Fst0.6.txt")
#positions <- paste(positions$V1, positions$V2, sep = ".")
#snp_positions <- cbind(positions,ggdata$loci)
#colnames(snp_positions) <- c("results.LocusName","SNP_Number")

#outflank_data <- read.csv("outflankResults.csv", header  = T)
#outflank_pvalue <- outflank_data[,c("results.LocusName","results.pvalues")]
#outflank_pvalue <- merge(snp_positions, outflank_pvalue, by = "results.LocusName")
#outflank_pvalue <- outflank_pvalue[which(gc1$gc$locus %in% outflank_pvalue$SNP_Number),]

#colnames(outflank_pvalue) <- c("results.LocusName","locus","results.pvalues")
#outflank_pvalue <- outflank_pvalue[order(outflank_pvalue$results.pvalues),]
#outflank_subset <- outflank_pvalue[1:num_SNPs,]

#outliers <- read.table("outlierTable.txt", header = T)
#colnames(outliers) <- "results.LocusName"
#outlier_merge <- merge(outliers, outflank_pvalue, by = "results.LocusName")

# Create input for genomic cline plot
#loci <- c(as.character(outflank_subset$locus),as.character(outlier_merge$locus))
#colors <- c(rep("grey",length(outflank_subset$locus)),rep("red",length(outlier_merge$locus)))

# Plot genomic clines
pdf("genomicCline.pdf")
gghybrid_plot <- plot_clinecurve(ggcline.object=gc1$gc,data.prep.object=prepdata$data.prep,esth.object=hindlabel$hi,cline.locus=loci,cline.col=colors,plot.genotype=TRUE,PLOIDY=2)
title(xlab="Hybrid index",ylab="Genotype frequency")
dev.off()
