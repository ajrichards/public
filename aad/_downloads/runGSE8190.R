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
library("GEOquery")
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
gseId <-  'GSE8190'
experimentId <- gseId
dataDir <- file.path(path.expand("~"),"data","samir","geo",gseId)
scriptDir <- file.path(path.expand("~"), "documents","ucdenver-docs","samir")

if (file.exists(file.path(path.expand("~"),"data","samir","raw",gseId)) == TRUE){
    rawDir <-  file.path(path.expand("~"),"data","samir","raw",gseId)
}else{
    rawDir <-  file.path("/","media",Sys.getenv("LOGNAME"),"mojo1","data","samir","raw",gseId)
}

## variables
if (outliers == TRUE){
    logFile = file.path(scriptDir,"results",paste(gseId,"-o-summary.log",sep=""))
}else{
    logFile = file.path(scriptDir,"results",paste(gseId,"-summary.log",sep=""))
}

cat(paste(format(Sys.time(), "%a %b %d %X %Y"),"\n",sep=""),file=logFile)
print(paste("saving...",logFile,sep=" "))

## load saved targets
targetsFile1 <- read.csv(file.path(scriptDir,"data","GSE8190-targets.csv"),header=TRUE)
targetsFile1$dx <- as.character(targetsFile1$dx)
targetsFile1$exposure <- as.character(targetsFile1$exposure)
targetsFile1$samples <- as.character(targetsFile1$samples)
targetsFile1$subject <- as.character(targetsFile1$subject)
targetsFile1$gender <- as.character(targetsFile1$gender)
targetsFile1$age <- as.character(targetsFile1$age)
targetsFile1$batch <- as.character(targetsFile1$batch)

control <- which(targetsFile1$dx==1)
allergy <- which(targetsFile1$dx==2)
asthma <- which(targetsFile1$dx==3)
asthmaAllergy <- which(targetsFile1$dx==4)
hdm <- which(targetsFile1$exposure=='HDM')
lps <- which(targetsFile1$exposure=='LPS')
saline <- which(targetsFile1$exposure=='saline')

tokeepControls <- intersect(control,saline)
tokeepAsthma <- intersect(c(asthma,asthmaAllergy),saline)
print(paste("used controls:",length(tokeepControls)))
print(paste("used asthma:",length(tokeepAsthma)))

## specify and remove nans
tokeepControls <- setdiff(tokeepControls,which(targetsFile1$age=='na'))
tokeepAsthma <- setdiff(tokeepAsthma,which(targetsFile1$age=='na')) 

tokeep <- c(tokeepControls,tokeepAsthma)
targets1 <- data.frame("sample"=targetsFile1$samples[tokeep],
                       "subject"=targetsFile1$subject[tokeep],
                       "phenotype"=targetsFile1$dx[tokeep],
                       "batch"=targetsFile1$batch[tokeep],
                       "age"=as.numeric(targetsFile1$age[tokeep]),
                       "gender"=targetsFile1$gender[tokeep],
                       "dx"=c(rep("control",length(tokeepControls)),rep("asthma",length(tokeepAsthma)))
                       )

targets1$FileName <- sapply(targets1$sample,function(x) paste(x,".txt",sep=""))

#print(targets$subject)
#print(targets$phenotype)

get_files <- function(subj,dx,column){
    indices <- which(targets1$subject == subj)
    dxIndices <- which(targets1$dx==dx)
    tokeep <- intersect(indices,dxIndices)
    
    if (column == 1){
        return(tokeep[c(TRUE, FALSE)])
    }else{
        return(tokeep[c(FALSE, TRUE)])
    }
}


print("...")

uniqueSubjects <- sort(as.character(unique(targets1$subject)))
print(uniqueSubjects)
asthma1 <- unlist(lapply(uniqueSubjects, get_files,dx='asthma',column=1))
asthma1 <- asthma1[!is.na(asthma1)]
asthma2 <- unlist(lapply(uniqueSubjects, get_files,dx='asthma',column=2))
asthma2 <- asthma2[!is.na(asthma2)]
control1 <- unlist(lapply(uniqueSubjects, get_files,dx='control',column=1))
control1 <- control1[!is.na(control1)]
control2 <- unlist(lapply(uniqueSubjects, get_files,dx='control',column=2))
control2 <- control2[!is.na(control2)]

tokeepAsthma <- c(asthma1,asthma2)
tokeepControls <- c(control1,control2)
cy3 <- c(rep("ref",length(control1)),rep("exp",length(control2)),rep("ref",length(asthma1)),rep("exp",length(asthma1)))
cy5 <- c(rep("exp",length(control1)),rep("ref",length(control2)),rep("exp",length(asthma1)),rep("ref",length(asthma1)))

tokeep <- c(tokeepControls,tokeepAsthma)
targets <- data.frame("subject"=targets1$subject[tokeep],
                      "phenotype"=targets1$dx[tokeep],
                      "batch"=targets1$batch[tokeep],
                      "age"=as.numeric(targets1$age[tokeep]),
                      "gender"=targets1$gender[tokeep],
                      "dx"=c(rep("control",length(tokeepControls)),rep("asthma",length(tokeepAsthma))),
                      "Cy3"=cy3,
                      "Cy5"=cy5,
                      "FileName"=targets1$FileName[tokeep]
                      )

targets$FileName <- sapply(targets$FileName,function(x) paste(x,".gz",sep=""))

## read raw data
print("...reading raw data")
normalization <- 'quantile'
drops <- c("Description")
threshold <- 4
exprFile <- file.path(rawDir,paste("expr-",normalization,".csv",sep=""))

## read in img data (two color)
print("loading data and saving expr file...")
x <- read.maimages(targets, source="agilent", path=rawDir, green.only=TRUE)

# background correct and normalize
y <- backgroundCorrect(x,method="normexp")
y <- normalizeBetweenArrays(y,method=normalization)
cat(paste("threshold: ",threshold,"\n",sep=""),file=logFile,append=TRUE)
cat(paste("original size:",dim(y)[1],'x', dim(y)[2],"\n",sep=" "),file=logFile,append=TRUE)
    
## filter out low expressed probes
neg95 <- apply(y$E[y$genes$ControlType==-1,],2,function(x) quantile(x,p=0.99))
cutoff <- matrix(1.1*neg95,nrow(y),ncol(y),byrow=TRUE)
isexpr <- rowSums(y$E > cutoff) >= threshold
y0 <- y[y$genes$ControlType==0 & isexpr,]
    
# write the expression values to file
write.csv(as.data.frame(y0),file=exprFile,quote=TRUE,row.names=FALSE)

## read exp file into eset
exprs <- read.csv(exprFile,header=TRUE, sep=",")

probeIds <- exprs[,"ProbeName"]
geneIds <- exprs[,"GeneName"]
systematicName <- exprs[,'SystematicName']
description <- exprs[,'Description']
probeUid <- exprs[,'ProbeUID']
exprs <- as.matrix(exprs[,9:dim(exprs)[2]])

print(paste("targets dimensions", dim(targets)[1],dim(targets)[2],sep=" "))
print(paste("features dimensions",dim(exprs)[1],dim(exprs)[2],sep=" "))

## average the intensities
controlExprs <- (exprs[,control1] + exprs[,control2]) / 2
asthmaExprs <- (exprs[,asthma1] + exprs[,asthma2]) / 2
exprs <- cbind(controlExprs,asthmaExprs)

## recreate the targets
covsTokeep = c(control1,asthma1)
targets <- data.frame("subject"=targets1$subject[covsTokeep],
                      "phenotype"=targets1$dx[covsTokeep],
                      "batch"=targets1$batch[covsTokeep],
                      "age"=as.numeric(targets1$age[covsTokeep]),
                      "gender"=targets1$gender[covsTokeep],
                      "dx"=c(rep("control",length(control1)),rep("asthma",length(asthma1)))
                      )

print(paste("features dimensions after avg",dim(exprs)[1],dim(exprs)[2],sep=" "))
print(paste("targets dimensions after avg", dim(targets)[1],dim(targets)[2],sep=" "))

## error checking
if (dim(targets)[1] != dim(exprs)[2]){
    print("ERROR: dimension mismatch")
    q()
}

targets$samples <- c(paste("control-",seq(1,length(control1)),sep=""),paste("asthma-",seq(1,length(asthma1)),sep=""))
rownames(targets) <- targets$samples
colnames(exprs) <- targets$samples

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
design <- model.matrix(~0+dx+age+gender+batch+subject)
cat(paste("model:","~0+dx+age+gender+batch+subject","\n",sep=" "),file=logFile,append=TRUE)
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
toRemove = c("Representative.Public.ID","Sequence Type","Sequence Source","Species.Scientific.Name","Annotation.Date","Sequence.Type","Gene.Title","Target.Description")
toRemove = c(toRemove,c("Gene.Ontology.Biological.Process","Gene.Ontology.Cellular.Component","Gene.Ontology.Molecular.Function","RefSeq.Transcript.ID"))
toRemove = c(toRemove,c("GO_ID","TIGR_ID","ACCESSION_STRING","CHROMOSOMAL_LOCATION","CYTOBAND","SEQUENCE","ENSEMBL_ID","DESCRIPTION","UNIGENE_ID"))
toRemove = c("foo")

results <- as.data.frame(results)[,!(names(results) %in% toRemove)]
results$systematicName <- systematicName
results$description <- description
results$probeUID <- probeUid

orderedInds <- order(results$adj.P.Val)

## write results
fileName <- file.path(scriptDir,'results',paste(experimentId,"-probes.tab",sep=""))
write.table(results[orderedInds,], file=fileName, sep="\t",row.names=FALSE)
   
## create a heatmap
print("...making heatmap")
results <- read.csv(fileName,header=TRUE,sep="\t")
print(names(results))
orderedInds <- order(results$adj.P.Val)
sigInds <- orderedInds[which(results$adj.P.Val[orderedInds] < 0.05)]

print(paste("There were ", length(sigInds), " significant probes",sep=""))
cat(paste("There were ", length(sigInds), " significant probes",sep=""),file=logFile,append=TRUE)

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
