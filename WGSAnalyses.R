#library(devtools)
#install_github("dinmatias/reconproGS")
#devtools::install_github("ribailey/gghybrid")
#install.packages("doMC")
#install.packages("MCMCpack")
#isntall.packages("hzar")
{library(reshape2)
library(ggplot2)
library(adegenet)
library(reconproGS)
library(gghybrid)
library(stringr)
library(MCMCpack)
library(hzar)
library(stringr)
library(tidyr)}

setwd("/Users/katherinefeldmann/Desktop/CU_Boulder/Taylor_Lab/Masters Thesis/Code")

###############
## ADMIXTURE ##
###############
admixtureTable <- read.table("poecileADMIXTURE.2.Q")
chickadeeInfo <- read.table(pipe("pbpaste"), sep = "\t", header = T, na.strings = c("","NA"))
samples <- read.csv("samplesID_pop.txt")
samples$Population <- NULL
admixtureTable <- cbind(admixtureTable, samples)
admixtureTable <- reshape2::melt(admixtureTable, id = "Sample_ID")
admixtureTable <- merge(admixtureTable, chickadeeInfo, by = "Sample_ID")
colnames(admixtureTable)[1:3] <- c("Individual","Population","Admixture")

orderNames <- c("B_25801","B_25902","B_25903","B_25904","BCCH1","BCCH10","BCCH1A","BCCH2","BCCH30","BCCH4","BCCH5","BCCH6","BCCH7","BCCH8","BCCH3","E458","E461","E463","E464","E466","E635","E809","E980","E983","E639","E454-2","E457","E941","E975","E442","E927","E699","E641","E441-2","E992","E452","E935","E447","E446","E942","E700","E696","E698","E451","E938","E990","E453","E688","E392","E346","E035","E031","E137","E082","E030","E128","E816","E350","E070","E414","E232","E959","E450","E039","E136","E304","E012","E085","E415","E987","E456-3","E247","E945","E271","E924","E060","E042","E592","E328","E103","E240","E836","E076","E175","E104","E986","E111","E958","E790","E533","E448","E239","E110","E095","E081","E630","E578","E143","E054","E583","E336","E162","E132","E058","E563","E564","CACH32.1","CACH5.1","E1032","E1033","E1034","E1035","S_77292","S_77293","S_77294","S_77296","S_77297")

admixtureTable$Individual <- factor(admixtureTable$Individual, levels = orderNames)

ggplot(admixtureTable, aes(x = Individual, y = Admixture, fill = Population))+
  geom_bar(position="fill",stat="identity")+
  scale_y_continuous(expand = c(0,0))+
  theme_bw() +
  theme(axis.text.x = element_text(size = 8, angle = 90), legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())+
  scale_fill_manual(values = c("firebrick1","dodgerblue3"))

##########
## HZAR ##
##########
## Quantitative Trail Model Input File ##
#------------------------------------------------------------------------------#
hzar <- read.table("inputHZAR.txt")
Location_ID <- c("DNY","DSU","ECW","HRP","INY","JSP","LHU","LSU")
meanInput <- hzar[,c("Location","EV1")]
meanAgg <- aggregate(EV1 ~ Location, data = meanInput, mean)
colnames(meanAgg) <- c("Location","pcaMean")
varInput <- hzar[,c("Location","EV1")]
varAgg <- aggregate(EV1 ~ Location, data = varInput, var)
colnames(varAgg) <- c("Location","pcaVar")
mean_var <- merge(meanAgg, varAgg, by = "Location")
nSamples <- data.frame(table(hzar$Location))
colnames(nSamples) <- c("Location","nSamples")
merged <- merge(mean_var, nSamples, by = "Location")
inputHzar <- cbind(Location_ID, merged)
inputHzar[which(inputHzar$Location == "LSU"),"Location"] <- "Louisiana State University"
inputHzar[which(inputHzar$Location == "Lehigh campus"),"Location"] <- "Lehigh University"
distances <- data.frame(Location_ID = c("LSU","ECW","DSU","LHU","JSP","HRP","INY","DNY"), Distance_km = c(0,580,1125,1132,1154,1180,1337,1372))
inputHzar <- merge(inputHzar, distances, by = "Location_ID")

write.table(inputHzar, file = "inputHZAR_updated.txt")

## Run HZAR: Modified from Georgy Semenov ##
data <- read.table("inputHZAR_updated.txt", header=T)

hzar.face<-list()
hzar.face$obs.list<-list()
hzar.face$model.list<-list()
hzar.face$fitRs.list<-list()
hzar.face$runs.list<-list()
hzar.face$analysis.list<-list()

# Create lists for each data object
hzar.face$obs.list$face <- hzar.doNormalData1DPops(data$Distance_km, data$pcaMean, data$pcaVar, data$nSamples, data$Location_ID)
hzar.plot.obsData(hzar.face$obs.list$face)

chainLength=1e5
mainSeed=
  list(
    A=c(978,544,99,596,528,124),
    B=c(544,99,596,528,124,978),
    C=c(99,596,528,124,978,544))

# If you have doMC, use foreach in parallel mode to speed up computation.
if(require(doMC)){
  registerDoMC()
} else { # Use foreach in sequential mode
  registerDoSEQ();
}

# Load Models
# what we want is for model list to be a list of loci names
# within each locus name are all the models
# want to apply the load model function to each 

hzar.face$model.list<-list()
for(i in 1:length(hzar.face$obs.list)){
  tmp<-hzar.face$obs.list[[i]]
  free_none <- hzar.makeCline1DNormal(tmp,  "none")
  free_both<- hzar.makeCline1DNormal(tmp, "both")
  free_right<- hzar.makeCline1DNormal(tmp,  "right")
  free_left<- hzar.makeCline1DNormal(tmp,  "left")
  free_mirror<- hzar.makeCline1DNormal(tmp, "mirror")
  mod.list<-list(free_none, free_both, free_right, free_left, free_mirror)
  names(mod.list)<-c("free_none", "free_both", "free_right", "free_left", "free_mirror")
  name<-names(hzar.face$obs.list[i])
  hzar.face$model.list[[name]]<-mod.list
}
hzar.face$obs.list$face$frame
hzar.face$model.list$face$free_none

# Modify all models to focus on observed region
for(i in 1:length(hzar.face$obs.list)){
  tmp<-hzar.face$model.list[[i]]
  s<-sapply(tmp,
            hzar.model.addBoxReq,
            -30 , 1400,
            simplify=FALSE)
  hzar.face$model.list[[i]]<-s 
}
print(hzar.face$model.list$face$free_none)

# Compile each of the models to prepare for fitting
hzar.face$fitRs.list$init <- list()
for(i in 1:length(hzar.face$obs.list)){
  tmp<-hzar.face$model.list[[i]]
  obs<-hzar.face$obs.list[[i]]
  s<-sapply(tmp,
            hzar.first.fitRequest.gC,
            obsData=obs,
            verbose=FALSE,
            simplify=FALSE)
  name<-names(hzar.face$model.list[i])
  hzar.face$fitRs.list$init[[name]]<-s
}
hzar.face$fitRs.list$init$face$free_none

# Run initial chains simultaneously using lapply
hzar.face$runs.list$init<-list()
for(i in 1:length(hzar.face$model.list)){
  tmp<-hzar.face$fitRs.list$init[[i]]
  mod<-lapply(tmp, function(x) hzar.doFit(x))
  name<-names(hzar.face$model.list[i])
  hzar.face$runs.list$init[[name]]<-mod
}
plot(hzar.mcmc.bindLL(hzar.face$runs.list$init$face$free_none))

# New Fit Requests Based on Initial Chains: Compile a new set of fit requests using the initial chains
hzar.face$fitRs.list$chains<-list()
for(i in 1:length(hzar.face$runs.list$init)){
  tmp<-hzar.face$runs.list$init[[i]]
  fitr.mod<-lapply(tmp,hzar.next.fitRequest)
  name<-names(hzar.face$model.list[i])
  hzar.face$fitRs.list$chains[[name]]<-fitr.mod
}
hzar.face$fitRs.list$chains$face$free_none$modelParam

# Replicate each fit request 3 times, keeping the original seeds while switching to a new seed channel.
for(i in 1:length(hzar.face$obs.list)){
  mod<-hzar.multiFitRequest(hzar.face$fitRs.list$chains[[i]], 
                            each=3, baseSeed=NULL)
  name<-names(hzar.face$model.list[i])
  hzar.face$fitRs.list$chains[[name]]<-mod
}

# Run a chain of 3 runs for every fit request
hzar.face$runs.list$chains<-list()
for(i in 1:length(hzar.face$model.list)){
  tmp<-hzar.face$fitRs.list$chains[[i]]
  long.chains<-hzar.doChain.multi(tmp,doPar=TRUE,
                                  inOrder=FALSE,
                                  count=3)
  name<-names(hzar.face$model.list[i])
  hzar.face$runs.list$chains[[name]]<-long.chains  
}

#Did the model converge?
summary(do.call(mcmc.list,
                lapply(hzar.face$runs.list$chains$face[1:3],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )
summary(do.call(mcmc.list,
                lapply(hzar.face$runs.list$chains$face[4:6],
                       function(x) hzar.mcmc.bindLL(x[[3]]) )) )

# Analysis: Start aggregation of data for analysis
# Create a model data group (hzar.dataGroup object) for each model from the initial runs.
hzar.face$analysis.list<-list()
for(i in 1:length(hzar.face$model.list)){
  mods.list<-hzar.face$runs.list$init[[i]]
  mod.runs<-lapply(mods.list, function(x) hzar.dataGroup.add(x))
  mod.names<-names(mods.list)
  name<-names(hzar.face$model.list[i])
  hzar.face$analysis.list[[name]]$initDGs[mod.names]<-mod.runs
}

# Create a hzar.obsDataGroup object from the four hzar.dataGroup just created
for(i in 1:length(hzar.face$obs.list)){
  tmp<-hzar.face$analysis.list[[i]]$initDGs
  obj<-hzar.make.obsDataGroup(tmp)
  hzar.face$analysis.list[[i]]$oDG <-obj
  hzar.face$analysis.list[[i]]$oDG<-hzar.copyModelLabels(hzar.face$analysis.list[[i]]$initDGs,
                                                         hzar.face$analysis.list[[i]]$oDG)
}  

# Convert all 27 runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
for(i in 1:length(hzar.face$obs.list)){
  chains<-hzar.face$runs.list$chains[[i]]
  odgs<-hzar.face$analysis.list[[i]]$oDG
  obj<-hzar.make.obsDataGroup(lapply(chains,hzar.dataGroup.add),odgs)
  hzar.face$analysis.list[[i]]$oDG<-obj                       
}

# Check data groups
print(summary(hzar.face$analysis.list$face$oDG$data.groups))

# Compare the cline models to the null model graphically
par(mfrow=c(2,1))
for(i in 1:length(hzar.face$obs.list)){
  hzar.plot.cline(hzar.face$analysis.list[[i]]$oDG)
}

# Do model selection based on the AICc scores
for(i in 1:length(hzar.face$obs.list)){
  AIC.tab<-hzar.AICc.hzar.obsDataGroup(hzar.face$analysis.list[[i]]$oDG)
  hzar.face$analysis.list[[i]]$AICcTable<-AIC.tab
}

# Print out the model with the minimum AICc score: free_none
for(i in 1:length(hzar.face$obs.list)){
  print(hzar.face$analysis.list[[i]]$model.name<-
          rownames(hzar.face$analysis.list[[i]]$AICcTable)[[which.min(hzar.face$analysis.list[[i]]$AICcTable$AICc)]])
}
hzar.face$analysis.list$face$AICcTable

# Extract the hzar.dataGroup object for the selected model
for(i in 1:length(hzar.face$obs.list)){
  hzar.face$analysis.list[[i]]$model.selected<-
    hzar.face$analysis.list[[i]]$oDG$data.groups[[hzar.face$analysis.list[[i]]$model.name]]
}

# Look at the variation in parameters for the selected model, this function will fail when null model is best model
for(i in 1:length(hzar.face$obs.list)){
  mod.sel<-hzar.face$analysis.list[[i]]$model.selected
  params<-hzar.face$analysis.list[[i]]$model.selected$data.param
  ml.params<-hzar.getLLCutParam(mod.sel, names(params))
  print(ml.params)
}

# Print the maximum likelihood cline for the selected model
for(i in 1:length(hzar.face$obs.list)){
  print(hzar.get.ML.cline(hzar.face$analysis.list[[i]]$model.selected))
}
hzar.get.ML.cline(hzar.face$analysis.list$face$model.selected)

# Plot the maximum likelihood cline for the selected model
for(i in 1:length(hzar.face$obs.list)){
  hzar.plot.cline(hzar.face$analysis.list[[i]]$model.selected)
  #par(new=TRUE)
}

# Plot the 95% credible cline region for the selected model also won't work for null model
hzar.plot.fzCline(hzar.face$analysis.list$face$model.selected, pch=3, fzCol="pink")
#------------------------------------------------------------------------------#

## Genetic Model Input File ##
#------------------------------------------------------------------------------#
# Create Input File
setwd("/Users/katherinefeldmann/Desktop/")
print("read frequency data")
LSU <- read.table(file = "headLSU.frq", sep = "\t", header = T, row.names = NULL)
LSU$Location <- "LSU"
ECW <- read.table(file = "headECW.frq", sep = "\t", header = T, row.names = NULL)
ECW$Location <- "ECW"
DSU <- read.table(file = "headDSU.frq", sep = "\t", header = T, row.names = NULL)
DSU$Location <- "DSU"
LHU <- read.table(file = "headLHU.frq", sep = "\t", header = T, row.names = NULL)
LHU$Location <- "LHU"
JSP <- read.table(file = "headJSP.frq", sep = "\t", header = T, row.names = NULL)
JSP$Location <- "JSP"
HRP <- read.table(file = "headHRP.frq", sep = "\t", header = T, row.names = NULL)
HRP$Location <- "HRP"
INY <- read.table(file = "headINY.frq", sep = "\t", header = T, row.names = NULL)
INY$Location <- "INY"
DNY <- read.table(file = "headDNY.frq", sep = "\t", header = T, row.names = NULL)
DNY$Location <- "DNY"

print("combine input data")
LCTN <- rbind(LSU, ECW, DSU, LHU, JSP, HRP, INY, DNY)

print("data processing")
colnames(LCTN) <- c("Scaffold","Position","NumAlleles","NumChr","AlleleFreq_1","AlleleFreq_2","Location")
LCTN <- LCTN %>% separate(AlleleFreq_1, sep = ":", c("Allele1","AlleleFrequency1"))
LCTN <- LCTN %>% separate(AlleleFreq_2, sep = ":", c("Allele2","AlleleFrequency2"))
LCTN$Scaffold <- paste(str_match(LCTN$Scaffold,"scaffold[0-9]"), LCTN$Position, sep = ".")
LCTN$AlleleFrequency1 <- as.numeric(LCTN$AlleleFrequency1)
LCTN$AlleleFrequency2 <- as.numeric(LCTN$AlleleFrequency2)
alleleInfo <- LCTN[,c("Scaffold","NumChr","Location")]
allele1LCTN <- LCTN[,c("Scaffold","Allele1","AlleleFrequency1","Location")]
allele2LCTN <- LCTN[,c("Scaffold","Allele2","AlleleFrequency2","Location")]
alleleInfo$Scaffold <- paste(alleleInfo$Scaffold,"samples",sep = ".")
alleleInfo <- dcast(alleleInfo, Location ~ Scaffold, value.var = "NumChr")
allele1LCTN <- dcast(allele1LCTN, Location ~ Scaffold + Allele1, value.var = "AlleleFrequency1")
allele2LCTN <- dcast(allele2LCTN, Location ~ Scaffold + Allele2, value.var = "AlleleFrequency2")
HZARinput <- merge(allele1LCTN, allele2LCTN, by = "Location")
HZARinput <- merge(HZARinput, alleleInfo, by = "Location")
HZARinput <- HZARinput[,order(colnames(HZARinput))] 

distances <- data.frame(Location = c("LSU","ECW","DSU","LHU","JSP","HRP","INY","DNY"), distance = c(0,580,1125,1132,1154,1180,1337,1372))
HZARinput <- merge(HZARinput, distances, by = "Location")

# Plot HZAR

## A typical chain length.  This value is the default setting in the package.
chainLength=1e5;                       

## Make each model run off a separate seed
mainSeed=list(A=c(596,528,124,978,544,99),B=c(528,124,978,544,99,596),C=c(124,978,544,99,596,528))

if(require(doMC)){
  ## If you have doMC, use foreach in parallel mode
  ## to speed up computation.
  registerDoMC()
} else {
  ## Use foreach in sequential mode
  registerDoSEQ();
}

## Molecular Analysis

allele <- "scaffold1.11711_C"
numSamples <- "scaffold1.11711.samples"

## Blank out space in memory to hold molecular analysis
if(length(apropos("^poecile$",ignore.case=FALSE)) == 0 ||
   !is.list(poecile) ) poecile <- list()

for(){
  ## We are doing just the one allele at one locus, but it is good to stay organized.
  poecile$Allele <- list();
  ## Space to hold the observed data
  poecile$Allele$obs <- list();
  ## Space to hold the models to fit
  poecile$Allele$models <- list();
  ## Space to hold the compiled fit requests
  poecile$Allele$fitRs <- list();
  ## Space to hold the output data chains
  poecile$Allele$runs <- list();
  ## Space to hold the analysed data
  poecile$Allele$analysis <- list();
  
  ## Locus Ada, Allele A from Brumfield et al 2001
  poecile$Allele$obs <- hzar.doMolecularData1DPops(HZARinput[,"distance"], HZARinput[,allele], HZARinput[,numSamples]);
  
  ## Look at a graph of the observed data
  #hzar.plot.obsData(poecile$Allele$obs);
  
  #### Show some plot modification commands for use by new users
  
  ## Make a helper function
  poecile.loadAdaAmodel <- function(scaling,tails,id=paste(scaling,tails,sep="."))
    poecile$Allele$models[[id]] <<- hzar.makeCline1DFreq(poecile$Allele$obs, scaling, tails)
  
  poecile.loadAdaAmodel("fixed","none","modelI");
  poecile.loadAdaAmodel("free" ,"none","modelII");
  poecile.loadAdaAmodel("free" ,"both","modelIII");
  
  ## Modify all models to focus on the region where the observed
  ## data were collected.
  ## Observations were between 0 and 570 km.
  poecile$Allele$models <- sapply(poecile$Allele$models, hzar.model.addBoxReq, -30 , 1400, simplify=FALSE)
  
  ## Compile each of the models to prepare for fitting
  poecile$Allele$fitRs$init <- sapply(poecile$Allele$models, hzar.first.fitRequest.old.ML, obsData=poecile$Allele$obs, verbose=FALSE, simplify=FALSE)
  ## Update the settings for the fitter if desired.
  poecile$Allele$fitRs$init$modelI$mcmcParam$chainLength <- chainLength;
  poecile$Allele$fitRs$init$modelI$mcmcParam$burnin <- chainLength %/% 10;
  poecile$Allele$fitRs$init$modelI$mcmcParam$seed[[1]] <- mainSeed$A
  
  poecile$Allele$fitRs$init$modelII$mcmcParam$chainLength <- chainLength;
  poecile$Allele$fitRs$init$modelII$mcmcParam$burnin <- chainLength %/% 10;
  poecile$Allele$fitRs$init$modelII$mcmcParam$seed[[1]] <- mainSeed$B 
  
  poecile$Allele$fitRs$init$modelIII$mcmcParam$chainLength <- chainLength;
  poecile$Allele$fitRs$init$modelIII$mcmcParam$burnin <- chainLength %/% 10;
  poecile$Allele$fitRs$init$modelIII$mcmcParam$seed[[1]] <- mainSeed$C 
  
  ## Run just one of the models for an initial chain
  poecile$Allele$runs$init <- list()
  poecile$Allele$runs$init$modelI <- hzar.doFit(poecile$Allele$fitRs$init$modelI)
  
  ## Plot the trace
  #plot(hzar.mcmc.bindLL(poecile$Allele$runs$init$modelI))
  
  ## Run another model for an initial chain
  poecile$Allele$runs$init$modelII <- hzar.doFit(poecile$Allele$fitRs$init$modelII)
  
  ## Plot the trace
  #plot(hzar.mcmc.bindLL(poecile$Allele$runs$init$modelII))
  
  ## Run another model for an initial chain
  poecile$Allele$runs$init$modelIII <- hzar.doFit(poecile$Allele$fitRs$init$modelIII)
  
  ## Plot the trace
  #plot(hzar.mcmc.bindLL(poecile$Allele$runs$init$modelIII))
  
  ## Compile a new set of fit requests using the initial chains 
  poecile$Allele$fitRs$chains <- lapply(poecile$Allele$runs$init, hzar.next.fitRequest)
  
  ## Replicate each fit request 3 times, keeping the original
  ## seeds while switching to a new seed channel.
  poecile$Allele$fitRs$chains <- hzar.multiFitRequest(poecile$Allele$fitRs$chains, each=3, baseSeed=NULL)
  
  ## Just to be thorough, randomize the initial value for each fit
  
  ## runif(9,-30,600) center for modelI, modelII, modelIII
  poecile$Allele$fitRs$chains[[1]]$modelParam$init["center"]= 120.08256
  poecile$Allele$fitRs$chains[[2]]$modelParam$init["center"]= 345.55072
  poecile$Allele$fitRs$chains[[3]]$modelParam$init["center"]= 158.34218
  poecile$Allele$fitRs$chains[[4]]$modelParam$init["center"]= 208.49748
  poecile$Allele$fitRs$chains[[5]]$modelParam$init["center"]= 89.08333
  poecile$Allele$fitRs$chains[[6]]$modelParam$init["center"]= 286.89100
  poecile$Allele$fitRs$chains[[7]]$modelParam$init["center"]= 529.46032
  poecile$Allele$fitRs$chains[[8]]$modelParam$init["center"]= 67.07035
  poecile$Allele$fitRs$chains[[9]]$modelParam$init["center"]= 513.06504
  
  ## runif(9,0,630) width for modelI, modelII, modelIII
  poecile$Allele$fitRs$chains[[1]]$modelParam$init["width"]= 283.09967 
  poecile$Allele$fitRs$chains[[2]]$modelParam$init["width"]= 436.41024  
  poecile$Allele$fitRs$chains[[3]]$modelParam$init["width"]= 84.82553 
  poecile$Allele$fitRs$chains[[4]]$modelParam$init["width"]= 385.90734 
  poecile$Allele$fitRs$chains[[5]]$modelParam$init["width"]= 279.79723  
  poecile$Allele$fitRs$chains[[6]]$modelParam$init["width"]= 38.92028  
  poecile$Allele$fitRs$chains[[7]]$modelParam$init["width"]= 42.87985
  poecile$Allele$fitRs$chains[[8]]$modelParam$init["width"]=  98.21601 
  poecile$Allele$fitRs$chains[[9]]$modelParam$init["width"]= 230.45149
  
  ## runif(6,0,1) pMin for modelII, modelIII
  poecile$Allele$fitRs$chains[[4]]$modelParam$init["pMin"]= 0.6778914 
  poecile$Allele$fitRs$chains[[5]]$modelParam$init["pMin"]= 0.4866557
  poecile$Allele$fitRs$chains[[6]]$modelParam$init["pMin"]=  0.4565415 
  poecile$Allele$fitRs$chains[[7]]$modelParam$init["pMin"]= 0.8985962 
  poecile$Allele$fitRs$chains[[8]]$modelParam$init["pMin"]= 0.3901906 
  poecile$Allele$fitRs$chains[[9]]$modelParam$init["pMin"]= 0.1531809
  
  ## runif(6,0,1) pMax for modelII, modelIII 
  poecile$Allele$fitRs$chains[[4]]$modelParam$init["pMax"]= 0.4567925
  poecile$Allele$fitRs$chains[[5]]$modelParam$init["pMax"]= 0.7801083 
  poecile$Allele$fitRs$chains[[6]]$modelParam$init["pMax"]= 0.3347732
  poecile$Allele$fitRs$chains[[7]]$modelParam$init["pMax"]=  0.8740998
  poecile$Allele$fitRs$chains[[8]]$modelParam$init["pMax"]=  0.4463791
  poecile$Allele$fitRs$chains[[9]]$modelParam$init["pMax"]=  0.1640979
  
  ## runif(3,0,630) deltaL for modelIII
  poecile$Allele$fitRs$chains[[7]]$modelParam$init["deltaL"]= 290.2459 
  poecile$Allele$fitRs$chains[[8]]$modelParam$init["deltaL"]= 242.7785 
  poecile$Allele$fitRs$chains[[9]]$modelParam$init["deltaL"]= 351.2703
  
  ## runif(3,0,1) tauL for modelIII 
  poecile$Allele$fitRs$chains[[7]]$modelParam$init["tauL"]= 0.3205238
  poecile$Allele$fitRs$chains[[8]]$modelParam$init["tauL"]= 0.9736836 
  poecile$Allele$fitRs$chains[[9]]$modelParam$init["tauL"]= 0.6674259
  
  ## runif(3,0,630) deltaR for modelIII
  poecile$Allele$fitRs$chains[[7]]$modelParam$init["deltaR"]=472.6632 
  poecile$Allele$fitRs$chains[[8]]$modelParam$init["deltaR"]=390.6439 
  poecile$Allele$fitRs$chains[[9]]$modelParam$init["deltaR"]=545.4608
  
  ## runif(3,0,1) tauR for modelIII
  poecile$Allele$fitRs$chains[[7]]$modelParam$init["tauR"]=0.9146836 
  poecile$Allele$fitRs$chains[[8]]$modelParam$init["tauR"]=0.2064311
  poecile$Allele$fitRs$chains[[9]]$modelParam$init["tauR"]= 0.2435238
  
  ## Go ahead and run a chain of 3 runs for every fit request
  poecile$Allele$runs$chains <-  hzar.doChain.multi(poecile$Allele$fitRs$chains, doPar=TRUE, inOrder=FALSE, count=3)
  
  ## Did modelI converge?
  #summary(do.call(mcmc.list, lapply(poecile$Allele$runs$chains[1:3], function(x) hzar.mcmc.bindLL(x[[3]]) )) )
  
  ## Did modelII converge?
  #summary(do.call(mcmc.list, lapply(poecile$Allele$runs$chains[4:6], function(x) hzar.mcmc.bindLL(x[[3]]) )) )
  
  ## Did modelIII converge?
  #summary(do.call(mcmc.list, lapply(poecile$Allele$runs$chains[7:9], function(x) hzar.mcmc.bindLL(x[[3]]) )) )
  
  ## All three models have three convergent chains, so additional runs are not needed.
  
  ## Just to show how, an additonal run is requested for each chain of modelI.  These three runs will not be included in the analysis stage.
  
  poecile$Allele$fitRs$extra.modelI <- lapply(lapply(poecile$Allele$runs$chains[1:3], function(x) x[[3]]), hzar.next.fitRequest)             
  
  ## Only do a single run.
  poecile$Allele$runs$extra.modelI <- hzar.doFit.multi(poecile$Allele$fitRs$extra.modelI, doPar=TRUE, inOrder=FALSE)
  
  ## Double check, are these three runs also convergent?
  summary(do.call(mcmc.list, lapply(poecile$Allele$runs$extra.modelI, hzar.mcmc.bindLL )))
  
  ## Start aggregation of data for analysis
  
  ## Create a model data group for the null model (expected allele frequency independent of distance along cline) to include in analysis.
  poecile$Allele$analysis$initDGs <- list(nullModel =  hzar.dataGroup.null(poecile$Allele$obs))
  
  ## Create a model data group (hzar.dataGroup object) for each model from the initial runs.
  poecile$Allele$analysis$initDGs$modelI <- hzar.dataGroup.add(poecile$Allele$runs$init$modelI)
  poecile$Allele$analysis$initDGs$modelII <- hzar.dataGroup.add(poecile$Allele$runs$init$modelII)
  poecile$Allele$analysis$initDGs$modelIII <- hzar.dataGroup.add(poecile$Allele$runs$init$modelIII)
  
  ## Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, copying the naming scheme (nullModel, modelI, modelII, modelIII).
  poecile$Allele$analysis$oDG <- hzar.make.obsDataGroup(poecile$Allele$analysis$initDGs)
  poecile$Allele$analysis$oDG <- hzar.copyModelLabels(poecile$Allele$analysis$initDGs, poecile$Allele$analysis$oDG)
  
  ## Convert all 27 runs to hzar.dataGroup objects, adding them to the hzar.obsDataGroup object.
  poecile$Allele$analysis$oDG <- hzar.make.obsDataGroup(lapply(poecile$Allele$runs$chains, hzar.dataGroup.add), poecile$Allele$analysis$oDG);
  
  ## Check to make sure that there are only four hzar.dataGroup objects named nullModel, modelI, modelII, and modelIII in the hzar.obsDataGroup object.
  #print(summary(poecile$Allele$analysis$oDG$data.groups))
  
  ## Compare the 3 cline models to the null model graphically
  #hzar.plot.cline(poecile$Allele$analysis$oDG);
  
  ## Do model selection based on the AICc scores
  #print(poecile$Allele$analysis$AICcTable <- hzar.AICc.hzar.obsDataGroup(poecile$Allele$analysis$oDG));
  
  ## Print out the model with the minimum AICc score
  #print(poecile$Allele$analysis$model.name <- rownames(poecile$Allele$analysis$AICcTable)[[ which.min(poecile$Allele$analysis$AICcTable$AICc )]])
  
  ## Extract the hzar.dataGroup object for the selected model
  poecile$Allele$analysis$model.selected <- poecile$Allele$analysis$oDG$data.groups[[poecile$Allele$analysis$model.name]]
  
  ## Look at the variation in parameters for the selected model
  #print(hzar.getLLCutParam(poecile$Allele$analysis$model.selected, names(poecile$Allele$analysis$model.selected$data.param)));
  print(hzar.getLLCutParam(poecile$Allele$analysis$model.selected, c("center","width")))
  
  ## Print the maximum likelihood cline for the selected model
  #print(hzar.get.ML.cline(poecile$Allele$analysis$model.selected))
  
  ## Plot the maximum likelihood cline for the selected model
  hzar.plot.cline(poecile$Allele$analysis$model.selected);
  
  ## Plot the 95% credible cline region for the selected model
  #hzar.plot.fzCline(poecile$Allele$analysis$model.selected);
}

#########
## MAP ##
#########
{library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggmap)
library(maps)
library(Hmisc)
library(rgeos)
library(ebirdst)
library(raster)
library(smoothr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(rgdal)}
select <- dplyr::select

bkcchi <- readOGR(dsn = "/Users/katherinefeldmann/Desktop/CU_Boulder/Taylor_Lab/Masters Thesis/Code/bkcchi/bkcchi-range-2019.gpkg")
carchi <- readOGR(dsn = "/Users/katherinefeldmann/Desktop/CU_Boulder/Taylor_Lab/Masters Thesis/Code/carchi/carchi-range-2019.gpkg")

world <- st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))
states <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))

locations <- read.csv("/Users/katherinefeldmann/Desktop/CU_Boulder/Taylor_Lab/Masters Thesis/Code/Sample_Locations.csv")
locations$Location <- factor(locations$Location, levels = c("DNY","INY","HRP","JSP","LHU","DSU","ECW","LSU"))
locations$Order <- rep(NA, length(locations$Location))
locations[which(locations$Location == "Lehigh University, Pennsylvania" | locations$Location == "Jacobsburg State Park, Pennsylvania" | locations$Location == "DeSales University, Pennsylvania"),"Order"] <- 1
locations[which(locations$Location == "Ithaca, New York" | locations$Location == "DeRuyter, New York" | locations$Location == "Hickory Run State Park, Pennsylvania"),"Order"] <- 2
locations[which(locations$Location == "East Carolina University, North Carolina" | locations$Location == "Louisiana State University, Louisiana"),"Order"] <- 3
locations <- locations[order(locations$Order),]

# Large Map #
ggplot()+
  geom_sf(data = world, fill = NA, lwd = 0.5)+
  geom_sf(data = states, fill = NA, lwd = 0.3)+
  geom_polygon(data = carchi, aes(x = long, y = lat, group = group), fill = "firebrick3", alpha = 0.7)+
  geom_polygon(data = bkcchi, aes(x = long, y = lat, group = group), fill = "dodgerblue4", alpha = 0.7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  coord_sf(xlim = c(-164,-57), ylim = c(26,70))
  #geom_point(data = locations, aes(x = Longitude, y = Latitude, fill = Location, shape = Location), size = 4)+
  #scale_fill_manual(values=c("navy","dodgerblue4","dodgerblue","darkorchid4","darkorchid3","orchid1","firebrick4","firebrick1"))+
  #scale_shape_manual(values=c(21,21,21,23,23,23,22,22))

# Small Map #
ggplot()+
  geom_sf(data = world, fill = NA, lwd = 0.5)+
  geom_sf(data = states, fill = NA, lwd = 0.3)+
  geom_polygon(data = carchi, aes(x = long, y = lat, group = group), fill = "firebrick3", alpha = 0.7)+
  geom_polygon(data = bkcchi, aes(x = long, y = lat, group = group), fill = "dodgerblue4", alpha = 0.7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  coord_sf(xlim = c(-77.5,-73.5), ylim = c(39,44), expand = FALSE)+
  geom_point(data = locations, aes(x = Longitude, y = Latitude, fill = Location, shape = Location), size = 10)+
  scale_fill_manual(values=c("navy","dodgerblue4","dodgerblue","darkorchid4","darkorchid3","orchid1","firebrick4","firebrick1"))+
  scale_shape_manual(values=c(21,21,21,23,23,23,22,22))

