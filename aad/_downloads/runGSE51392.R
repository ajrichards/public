#!/usr/bin/Rscript
# source("https://bioconductor.org/biocLite.R")
# biocLite("hgu133plus2.db")
# phenotype


#library(affy)
library("limma")
library("GEOquery")
library("Biobase")
library("affy")

###############################################################################
args <-commandArgs(TRUE)
if (length(args) == 0){
    outliers = FALSE
}else{
    outliers = TRUE
}

## set dirs
gseId <-  'GSE51392'
experimentId <-  gseId
dataDir <- file.path(path.expand("~"),"data","samir","geo",gseId)
rawDir <- file.path(path.expand("~"),"data","samir","raw",gseId)
scriptDir <- file.path(path.expand("~"),"documents","ucdenver-docs","samir")

## variables
if (outliers == TRUE){
    logFile = file.path(scriptDir,"results",paste(gseId,"-o-summary.log",sep=""))
}else{
    logFile = file.path(scriptDir,"results",paste(gseId,"-summary.log",sep=""))
}

cat(paste(format(Sys.time(), "%a %b %d %X %Y"),"\n",sep=""),file=logFile)
print(paste("saving...",logFile,sep=" "))

## load saved targets
targets <- read.csv(file.path(scriptDir,"data","GSE51392-targets.csv"),header=TRUE)
celFiles <- list.files(rawDir)

getFileName <- function(fileBase){
    fileBase <- paste("_",fileBase,sep="")
    fileName <- "nomatch"
    for (fname in celFiles){
        if (grepl(fileBase,fname,fixed=TRUE)){
            fileName <- fname
        }
    }
    return(fileName)
}


targets$fileName <- sapply(as.character(targets$celFile),getFileName)

## read raw data
print("...reading raw data")
ab <- ReadAffy(filenames=targets$fileName,celfile.path=rawDir)
eset <- rma(ab)

rownames(targets) = targets$fileName

## error checking
if (dim(targets)[1] != dim(eset)[2]){
    print("ERROR: dimension mismatch")
    q()
}

if (all(rownames(targets)==colnames(eset)) == FALSE){
    print("ERROR: ordering mismatch")
    q()
}


## remove rhinitis
healthy <- which(targets$phenotype=='healthy')
asthma <- which(targets$phenotype=='asthma')
tokeep <- union(asthma,healthy)

targets$tissue = as.character(targets$tissue)
targets$tissue[which(targets$tissue=='bronchial epithelium')] = 'bronchial'
targets$tissue[which(targets$tissue=='nasal epithelium')] = 'nasal'

oldTargets <- targets
targets <- data.frame("dx" = factor(as.character(oldTargets$phenotype)[tokeep],levels=c("healthy","asthma")),
                      "gender" = as.factor(as.character(oldTargets$gender)[tokeep]),
                      "tissue" = as.factor(as.character(oldTargets$tissue)[tokeep]),
                      "age" = as.numeric(oldTargets$age)[tokeep],
                      "subject" = as.factor(as.character(oldTargets$paired)[tokeep]),
                      "filename" = as.character(oldTargets$fileName)[tokeep])

rownames(targets) <- targets$filename
eset <- eset[,tokeep]
dx <- targets$dx
gender <- targets$gender
tissue <- targets$tissue
age <- as.numeric(targets$age)
subject <- targets$subject
probeIds <- rownames(eset)

outlierData <- read.csv("sample-outliers.csv",stringsAsFactors=FALSE)
outlierInds <- strsplit(outlierData$pca.outliers[which(outlierData$gseId==sub("-1","",gseId))],";")
outlierInds <- as.numeric(unlist(outlierInds))

if (outliers == TRUE){
                    
    toKeep <- setdiff(1:length(age),outlierInds)
    experimentId <- paste(experimentId,"-o",sep="")

    dx <- dx[toKeep]
    age <- age[toKeep]
    gender <-  gender[toKeep]
    subject <- subject[toKeep]
    tissue <-  tissue[toKeep]
    eset <- eset[,toKeep]
    targets <- targets[toKeep,]
}

print(paste("targets dimensions", dim(targets)[1],dim(targets)[2],sep=" "))
print(paste("features dimensions",dim(eset)[1],dim(eset)[2],sep=" "))

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
towriteEset <- as.data.frame(exprs(eset))

write.csv(towriteEset,file=file.path(scriptDir,"data",paste(experimentId,"-exp.csv",sep="")),quote=FALSE,row.names=TRUE)
write.csv(as.data.frame(targets),file=file.path(scriptDir,"data",paste(experimentId,"-cov.csv",sep="")),quote=FALSE)

## plots
png(file=file.path(scriptDir,'figs',paste(experimentId,'-densities.png',sep="")),width=8,height=8,unit='in',res=300)
plotDensities(eset,legend=FALSE,main="normalized array densities")
dev.off()

png(file=file.path(scriptDir,'figs',paste(experimentId,'-limma.png',sep="")),width=8,height=4,unit='in',res=300)
par(mfrow = c(1,2))

############################################################################
## probe level fit
print(summary(dx))
print(summary(age))
print(summary(gender))
print(summary(subject))


design <- model.matrix(~0+dx+age+gender+tissue)
cat(paste("model:","~0+dx+age+gender+tissue","\n",sep=" "),file=logFile,append=TRUE)
contMatrix <- makeContrasts(asthmaVsHealthy=dxasthma-dxhealthy,levels=design)
fitProbes <- lmFit(eset,design)
fitProbes2 <- contrasts.fit(fitProbes,contMatrix)
fitProbes2 <- eBayes(fitProbes2,trend=TRUE)

#design <- model.matrix(~0+dx+age+gender)
#corfit <- duplicateCorrelation(eset, design, block=targets$subject)
#cat(paste("model:","~0+dx+age+gender","\n",sep=" "),file=logFile,append=TRUE)
#cat(paste("model:","blocking factor = subject","\n",sep=" "),file=logFile,append=TRUE)
#fitProbes <- lmFit(eset, design, block=targets$subject, correlation=corfit$consensus.correlation)
#contMatrix <- makeContrasts(asthmaVsHealthy=dxasthma-dxhealthy,levels=design)
#fitProbes2 <- contrasts.fit(fitProbes,contMatrix)
#fitProbes2 <- eBayes(fitProbes2,trend=TRUE)
    
## volcano plot
print("...making volcano plot")
plotSA(fitProbes2, main="expression - probes")
volcanoplot(fitProbes2,highlight=8,names=probeIds,main="volcano - probes")
x <- dev.off()

#ids <- tbl[["ID"]]
#library("hgu95av2.db")
#print(names(fitProbes2))
## filter on gene names
#print("...filtering")
#if ('genes.ENTREZ_GENE_ID' %in% names(as.data.frame(fitProbes2)) == FALSE){
#    print(names(as.data.frame(fitProbes2)))
#    print("ERROR: need to specify correct gene id")
#    q()
#}
#geneIds = as.character(as.data.frame(fitProbes2)$genes.ENTREZ_GENE_ID)
#a <- which(geneIds != "")
#b <- which(is.na(geneIds)==FALSE)
#hasGene <- intersect(a,b)
#cat(paste("gene filter ",length(hasGene),"/",length(geneIds),"\n",sep=""),file=logFile,append=TRUE)
#fitProbes2 <- fitProbes2[hasGene,]
#print(paste("...filtering on gene names ",length(hasGene),"/",length(geneIds),sep=""))

## filter out specific columns
print("getting there...")
results <- topTable(fitProbes2,n=nrow(eset),sort.by="none",coef="asthmaVsHealthy",adjust.method="BH")
#print(names(results))
#print(names(as.data.frame(results)))
#q()

toRemove = c("Representative.Public.ID","Sequence Type","Sequence Source","Species.Scientific.Name","Annotation.Date","Sequence.Type","Gene.Title","Target.Description")
toRemove = c(toRemove,c("Gene.Ontology.Biological.Process","Gene.Ontology.Cellular.Component","Gene.Ontology.Molecular.Function","RefSeq.Transcript.ID"))
towrite = 
results <- as.data.frame(results)[,!(names(results) %in% toRemove)]
results$probes <- probeIds
orderedInds <- order(results$adj.P.Val)
    
## write results
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
