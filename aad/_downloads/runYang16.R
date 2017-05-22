#!/usr/bin/Rscript
# source("https://bioconductor.org/biocLite.R")
# biocLite("hgu133plus2.db")
# phenotype

#1. control - 74
#2. allergy only (no asthma) - 70
#3. asthma only (no allergy) - 24
#4. allergy and asthma - 82 

#library(affy)
library("limma")
library("Biobase")
library("gtools")

###############################################################################
args <-commandArgs(TRUE)
if (length(args) == 0){
    outliers = FALSE
}else{
    outliers = TRUE
}

## set dirs
gseId <-  'yang16'
experimentId <- gseId
scriptDir <- file.path(path.expand("~"),"documents","ucdenver-docs","ivana","u19-asthma","expression")
dataDir <- file.path(path.expand("~"),"data","ivana","u19-asthma","Expression")
agilentDir <- file.path(dataDir,"Agilent_txt_files")
scriptDir <- file.path(path.expand("~"), "documents","ucdenver-docs","samir")

## variables
if (outliers == TRUE){
    logFile = file.path(scriptDir,"results",paste(gseId,"-o-summary.log",sep=""))
}else{
    logFile = file.path(scriptDir,"results",paste(gseId,"-summary.log",sep=""))
}
cat(paste(format(Sys.time(), "%a %b %d %X %Y"),"\n",sep=""),file=logFile)
print(paste("saving...",logFile,sep=" "))

## load saved targets
targetsFile <- file.path(dataDir,"agilent-covariates.csv")
covs <- readTargets(targetsFile,sep=",")

dxList = as.character(covs$Dx)
normalInds <- which(dxList=="Normal")
asthmaInds <- which(dxList=="Asthma")

tissueList = as.character(covs$Tissue)
phenotypeList = as.character(covs$category)
pbmcInds <- which(tissueList=="PBMC")
nasalInds <- which(tissueList=="Nasal")
bronchInds <- which(tissueList=="Bronchial")

tokeep <- bronchInds
targets <- data.frame("dx"=as.character(covs$Dx)[tokeep],
                      "subject_id"=as.character(covs$ID)[tokeep],
                      "FileName"=as.character(rownames(covs))[tokeep],
                      "age" = as.numeric(as.character(covs$Age)[tokeep]),
                      "gender" = as.character(covs$Gender)[tokeep],
                      "tissue" = as.character(covs$Tissue)[tokeep],
                      "batch" = as.character(covs$Subarray)[tokeep])          

targets$FileName <- sapply(targets$FileName,function(x) paste(x,".txt",sep=""))
dxList = as.character(targets$dx)
normalInds <- which(dxList=="Normal")
asthmaInds <- which(dxList=="Asthma")
targets$dx = as.character(targets$dx)
targets$dx[normalInds] = 'control'
targets$dx[asthmaInds] = 'asthma'

normalization <- 'quantile'
drops <- c("Description")
threshold <- 4
exprFile <- file.path(dataDir,paste("expr-",normalization,".csv",sep=""))

## read in img data
print("loading data and saveing expr file...")
x <- read.maimages(targets, source="agilent", green.only=TRUE,path=agilentDir)

## background correct and normalize
y <- backgroundCorrect(x,method="normexp")
y <- normalizeBetweenArrays(y,method=normalization)

## filter out low expressed probes
neg95 <- apply(y$E[y$genes$ControlType==-1,],2,function(x) quantile(x,p=0.99))
cutoff <- matrix(1.1*neg95,nrow(y),ncol(y),byrow=TRUE)
isexpr <- rowSums(y$E > cutoff) >= threshold
y0 <- y[y$genes$ControlType==0 & isexpr,]

## write the expression values to file
write.csv(as.data.frame(y0),file=exprFile,quote=FALSE)

## read raw data
print("...reading raw data")
exprs <- read.csv(exprFile,header=TRUE, sep=",",row.names=1)
probeIds <- exprs[,"ProbeName"]
geneIds <- exprs[,"SystematicName"]
systematicName <- exprs[,"SystematicName"]
exprs <- as.matrix(exprs[,6:dim(exprs)[2]])
rownames(targets) = sapply(targets$FileName,function(x) gsub(".txt","",x))

print(paste("targets dimensions", dim(targets)[1],dim(targets)[2],sep=" "))
print(paste("features dimensions",dim(exprs)[1],dim(exprs)[2],sep=" "))

## error checking
if (dim(targets)[1] != dim(exprs)[2]){
    print("ERROR: dimension mismatch")
    q()
}

## sanity check
if (all(rownames(targets)==colnames(exprs)) == FALSE){
    print("ERROR: ordering mismatch")
    q()
}

#############################################################################
## load data and specify covariates
fData <- AnnotatedDataFrame(data.frame(probes=probeIds,genes=geneIds))
pData <- AnnotatedDataFrame(targets)
eset <- ExpressionSet(assayData=exprs,phenoData=pData,featureData=fData)

############################################################################
## specify covariates
dx <- factor(targets$dx,levels=c("control","asthma"))
age <- as.numeric(targets$age)
gender <- as.factor(targets$gender)
batch <- as.factor(targets$batch)
subject <- as.factor(targets$subject)

outlierData <- read.csv("sample-outliers.csv",stringsAsFactors=FALSE)
outlierInds <- strsplit(outlierData$pca.outliers[which(outlierData$gseId==sub("-1","",gseId))],";")
outlierInds <- as.numeric(unlist(outlierInds))

if (outliers == TRUE){
    toKeep <- setdiff(1:length(age),outlierInds)
    experimentId <- paste(experimentId,"-o",sep="")

    dx <- dx[toKeep]
    age <- age[toKeep]
    gender <-  gender[toKeep]
    batch <- batch[toKeep]
    subject <- subject[toKeep]
    eset <- eset[,toKeep]
    targets <- targets[toKeep,]
}

## write the data to a file
cat(paste("probes",dim(eset)[1],"\n"),file=logFile,append=TRUE)
cat(paste("samples",dim(eset)[2],"\n"),file=logFile,append=TRUE)
towriteEset <- as.data.frame(exprs(eset))
towriteEset <- cbind(probes=probeIds,towriteEset)

#towriteEset$probes <- probeIds
write.csv(towriteEset,file=file.path(scriptDir,"data",paste(experimentId,"-exp.csv",sep="")),quote=FALSE,row.names=FALSE)
write.csv(as.data.frame(targets),file=file.path(scriptDir,"data",paste(experimentId,"-cov.csv",sep="")),quote=FALSE)
    
## plots
png(file=file.path(scriptDir,'figs',paste(experimentId,'-densities.png',sep="")),width=8,height=8,unit='in',res=300)
plotDensities(eset,legend=FALSE,main="normalized array densities")
dev.off()

png(file=file.path(scriptDir,'figs',paste(experimentId,'-limma.png',sep="")),width=8,height=4,unit='in',res=300)
par(mfrow = c(1,2))

#probeIds <- rownames(toWrite)

############################################################################
## probe level fit
design <- model.matrix(~0+dx+age+gender+batch)
cat(paste("model:","~0+dx+age+gender+batch","\n",sep=" "),file=logFile,append=TRUE)
print(colnames(design))
contMatrix <- makeContrasts(asthmaVsControl=dxasthma-dxcontrol,levels=design)
fitProbes <- lmFit(eset,design)
fitProbes2 <- contrasts.fit(fitProbes, contMatrix)
fitProbes2 <- eBayes(fitProbes2,trend=TRUE)

## volcano plot
print("...making volcano plot")
plotSA(fitProbes2, main="expression - probes")
volcanoplot(fitProbes2,highlight=8,names=probeIds,main="volcano - probes")
x <- dev.off()

## filter on gene names
print(names(as.data.frame(fitProbes2)))
print("...filtering")
if ('genes.genes' %in% names(as.data.frame(fitProbes2)) == FALSE){
    print("ERROR: need to specify correct gene id")
    q()
}

geneIds = as.character(as.data.frame(fitProbes2)$genes.genes)
a <- which(geneIds != "")
b <- which(is.na(geneIds)==FALSE)
hasGene <- intersect(a,b)
cat(paste("gene filter ",length(hasGene),"/",length(geneIds),"\n",sep=""),file=logFile,append=TRUE)
fitProbes2 <- fitProbes2[hasGene,]
print(paste("...filtering on gene names ",length(hasGene),"/",length(geneIds),sep=""))

## filter out specific columns
results <- topTable(fitProbes2,n=nrow(eset),sort.by="none",coef="asthmaVsControl",adjust.method="BH")
#toRemove = c("Representative.Public.ID","Sequence Type","Sequence Source","Species.Scientific.Name","Annotation.Date","Sequence.Type","Gene.Title","Target.Description")
#toRemove = c(toRemove,c("Gene.Ontology.Biological.Process","Gene.Ontology.Cellular.Component","Gene.Ontology.Molecular.Function","RefSeq.Transcript.ID"))
#toRemove = c(toRemove,c("GO_ID","TIGR_ID","ACCESSION_STRING","CHROMOSOMAL_LOCATION","CYTOBAND","SEQUENCE","ENSEMBL_ID","DESCRIPTION","UNIGENE_ID"))

print("configuring results...")
results <- as.data.frame(results)#[,!(names(results) %in% toRemove)]
results$systematicName <- systematicName
results$probeId <- probeIds
orderedInds <- order(results$adj.P.Val)

## write results
print("writing resutls...")
fileName <- file.path(scriptDir,'results',paste(experimentId,"-probes.tab",sep=""))
write.table(results[orderedInds,], file=fileName, sep="\t",row.names=FALSE)


## create a heatmap
print("...making heatmap")
results <- read.csv(fileName,header=TRUE,sep="\t")
orderedInds <- order(results$adj.P.Val)
sigInds <- orderedInds[which(results$adj.P.Val[orderedInds] < 0.05)]

print(paste("There were ", length(sigInds), " significant probes",sep=""))
cat(paste("There were ", length(sigInds), " significant probes",sep=""),file=logFile,append=TRUE)
png(file=file.path(scriptDir,'figs',paste(experimentId,'-heatmap.png',sep="")),width=8,height=4,unit='in',res=300)
    
color_map <- function(diag) { if (diag=="asthma") "#FF0000" else "#0000FF" }
patientcolors <- unlist(lapply(dx, color_map))

heatmap(exprs(eset[orderedInds[1:50],]), col=topo.colors(100),ColSideColors=patientcolors)
null = dev.off()

print("done")
print("run parseYang16.py to get update gene names")

