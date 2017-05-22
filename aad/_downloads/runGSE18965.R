#!/usr/bin/Rscript
# source("https://bioconductor.org/biocLite.R")
# biocLite("hgu133plus2.db")

#library(affy)
library("limma")
library("GEOquery")
library("Biobase")

###############################################################################
args <-commandArgs(TRUE)
if (length(args) == 0){
    outliers = FALSE
}else{
    outliers = TRUE
}

## set dirs
gseId <-  'GSE18965'
dataDir <- file.path(path.expand("~"),"data","samir","geo",gseId)
scriptDir <- file.path(path.expand("~"),"documents","ucdenver-docs","samir")

## variables
if (outliers == TRUE){
    logFile = file.path(scriptDir,"results",paste(gseId,"-o-summary.log",sep=""))
}else{
    logFile = file.path(scriptDir,"results",paste(gseId,"-summary.log",sep=""))
}

cat(paste(format(Sys.time(), "%a %b %d %X %Y"),"\n",sep=""),file=logFile)
print(paste("saving...",logFile,sep=" "))

setwd(file.path(dataDir))
gseFiles <- Sys.glob("*.gz")
print("fetching...")
gseObjs <- getGEO(gseId,destdir=dataDir,GSEMatrix=TRUE)
cat(paste("total files: ",length(gseFiles),"\n",sep=""),file=logFile,append=TRUE)

## load saved targets
print("...reading targets")
targets = read.csv(file.path(scriptDir,"data","GSE18965-targets.csv"),header=TRUE)

## print the basic data
gseIndx = 1
#for (gseIndx in 1:length(gseFiles)){
experimentId <- paste(gseId,"-",gseIndx,sep="")
eset <- gseObjs[[gseIndx]]

cat(paste("file: ",gseIndx,"/",length(gseFiles),"\n",sep=""),file=logFile,append=TRUE)
cat(paste("file name: ",gseFiles[gseIndx],"\n",sep=""),file=logFile,append=TRUE)
cat(paste("status:",pData(eset)$status[1],"\n",sep=" "),file=logFile,append=TRUE)
cat(paste("data processing:",pData(eset)$data_processing[1],"\n",sep=" "),file=logFile,append=TRUE)
cat(paste("type:",pData(eset)$type[1],"\n",sep=" "),file=logFile,append=TRUE)

cat(paste("probes",dim(eset)[1],"\n"),file=logFile,append=TRUE)
cat(paste("samples",dim(eset)[2],"\n"),file=logFile,append=TRUE)
cat(paste("platform_id:",pData(eset)$platform_id[1],"\n",sep=" "),file=logFile,append=TRUE)

targetInds <- sapply(c(names(pData(eset))),function(x) grepl("characteristics",x))
covariates <- pData(eset)[,c(names(pData(eset))[targetInds])]
    
print(paste("targets dimensions", dim(targets)[1],dim(targets)[2],sep=" "))
print(paste("feature dimensions",dim(eset)[1],dim(eset)[2],sep=" "))

## match names
reorderedInds <- sapply(targets$sample,function(x) which(colnames(eset)==x))
eset <- eset[,reorderedInds]
        
## get variables    
aa = which(targets$phenotype=='AA')
hn = which(targets$phenotype=='HN')
rownames(targets) <- targets$sample

targets$phenotype <- as.character(targets$phenotype)
targets$phenotype[aa] <- 'asthma'
targets$phenotype[hn] <- 'control'
dx <- factor(targets$phenotype,levels=c("control","asthma"))
targets$dx <- dx
age <- targets$age
gender <- targets$gender 

outlierData <- read.csv(file.path(scriptDir,"sample-outliers.csv"),stringsAsFactors=FALSE)
outlierInds <- strsplit(outlierData$pca.outliers[which(outlierData$gseId==sub("-1","",gseId))],";")
outlierInds <- as.numeric(unlist(outlierInds))
    
if (outliers == TRUE){
    toKeep <- setdiff(1:length(age),outlierInds)
    experimentId <- paste(experimentId,"-o",sep="")

    dx <- dx[toKeep]
    age <- age[toKeep]
    gender <-  gender[toKeep]
    eset <- eset[,toKeep]
    targets <- targets[toKeep,]
}
    
## error checking
if (dim(targets)[1] != dim(eset)[2]){
    print("ERROR: dimension mismatch")
    q()
}

if (all(rownames(targets)==colnames(eset)) == FALSE){
    print("ERROR: ordering mismatch")
    q()
}
    
## write the data to a file
toWrite <- as.data.frame(exprs(eset))
probeIds <- rownames(toWrite)
print(names(toWrite))

write.csv(toWrite,file=file.path(scriptDir,"data",paste(experimentId,"-exp.csv",sep="")),quote=FALSE)
write.csv(as.data.frame(targets),file=file.path(scriptDir,"data",paste(experimentId,"-cov.csv",sep="")),quote=FALSE)
    
## plots
png(file=file.path(scriptDir,'figs',paste(experimentId,'-densities.png',sep="")),width=8,height=8,unit='in',res=300)
plotDensities(eset,legend=FALSE,main="normalized array densities")
dev.off()

png(file=file.path(scriptDir,'figs',paste(experimentId,'-limma.png',sep="")),width=8,height=4,unit='in',res=300)
par(mfrow = c(1,2))


############################################################################
## probe level fit
design <- model.matrix(~0+dx+age+gender)
cat(paste("model:","~0+dx+age+gender","\n",sep=" "),file=logFile,append=TRUE)
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
print("...filtering")
if ('genes.ENTREZ_GENE_ID' %in% names(as.data.frame(fitProbes2)) == FALSE){
    print("ERROR: need to specify correct gene id")
    print(names(as.data.frame(fitProbes2)))
    q()
}
geneIds = as.character(as.data.frame(fitProbes2)$genes.ENTREZ_GENE_ID)
a <- which(geneIds != "---")
b <- which(is.na(geneIds)==FALSE)
hasGene <- intersect(a,b)
cat(paste("gene filter ",length(hasGene),"/",length(geneIds),"\n",sep=""),file=logFile,append=TRUE)
fitProbes2 <- fitProbes2[hasGene,]
print(paste("...filtering on gene names keeping",length(hasGene),"/",length(geneIds),sep=""))

## filter out specific columns
results <- topTable(fitProbes2,n=nrow(eset),sort.by="none",coef="asthmaVsControl",adjust.method="BH")
toRemove = c("Representative.Public.ID","Sequence Type","Sequence Source","Species.Scientific.Name","Annotation.Date","Sequence.Type","Gene.Title","Target.Description")
toRemove = c(toRemove,c("Gene.Ontology.Biological.Process","Gene.Ontology.Cellular.Component","Gene.Ontology.Molecular.Function","RefSeq.Transcript.ID"))
toRemove = c(toRemove,c("GO_ID","TIGR_ID","ACCESSION_STRING","CHROMOSOMAL_LOCATION","CYTOBAND","SEQUENCE","ENSEMBL_ID","DESCRIPTION","UNIGENE_ID","GB_LIST","seqname","RANGE_GB","RANGE_START","RANGE_STRAND","RANGE_STOP","mrna_assignment","category","B","Sequence.Source"))
    
results <- as.data.frame(results)[,!(names(results) %in% toRemove)]
    
## write results
fileName <- file.path(scriptDir,'results',paste(experimentId,"-probes.tab",sep=""))
write.table(results, file=fileName, sep="\t",row.names=FALSE)
   
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
