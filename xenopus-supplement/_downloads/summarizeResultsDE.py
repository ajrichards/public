#!/usr/bin/env python
"""
read the output of a DESeq2 analysis and summarize the results
create a new output file with associated gene names and go terms
"""

import os,sys,csv,getopt,re
import numpy as np
from htsint.database import db_connect,Gene,GoTerm,Taxon
from htsint.tools import read_matrix,read_de_results,print_rest_table_contents
from htsint.blast import BlastMapper
from htsint.database import db_connect,Taxon,Gene
from sqlalchemy.sql import select


session,engine = db_connect()

## read input 
if len(sys.argv) < 2:
    raise Exception(sys.argv[0] + " -a assembly [-a] is 'dn' or 'gg'")
try:
    optlist, args = getopt.getopt(sys.argv[1:], 'a:t:')
except getopt.GetoptError:
    raise Exception(sys.argv[0] + " -a assembly [-a] is 'dn' or 'gg'")

assembly = None
for o,a in optlist:
    if o == '-a':
        assembly = a

if assembly not in ['gg','dn','ref']:
    raise Exception('Invalid assembly argument')

homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus")
if assembly == 'ref':
    sourceDir = 'reference'
else:
    sourceDir = "%s-trinity"%(assembly)

## load ref2gene
reader = csv.reader(open("../gene2ref.tab","r"),delimiter="\t")
ref2gene = {}
for linja in reader:
    geneId = linja[1]
    proteinAccession = linja[5]
    ref2gene[proteinAccession] = geneId

## load gene data
session,engine = db_connect()
conn = engine.connect()
#taxId = session.query(Taxon).filter_by(ncbi_id='8364').first().id
s = select([Taxon.id,Taxon.ncbi_id,Taxon.name]).where(Taxon.ncbi_id.in_(['8364']))
_taxaQueries = conn.execute(s)
taxaQueries = _taxaQueries.fetchall()
gene2taxa,gene2desc,gene2sym = {},{},{}

for tquery in taxaQueries:
    s = select([Gene.taxa_id,Gene.ncbi_id,Gene.description,Gene.symbol],Gene.taxa_id==tquery['id'])
    _geneQueries = conn.execute(s)
    geneQueries = _geneQueries.fetchall()
    gene2taxa.update(dict([(str(r['ncbi_id']),str(r['taxa_id'])) for r in geneQueries]))
    gene2desc.update(dict([(str(r['ncbi_id']),str(r['description'])) for r in geneQueries]))
    gene2sym.update(dict([(str(r['ncbi_id']),str(r['symbol'])) for r in geneQueries]))

## load blast mapper        
summaryFile1 = os.path.join(homeDir,sourceDir,"blast-%s-parsed_summary.csv"%(assembly))
summaryFile2 = os.path.join(homeDir,sourceDir,"blast-xt-parsed_summary.csv")
bm = BlastMapper()
bmapSP = bm.load_summary(summaryFile1,trinityGene=False,best=True)
bmapXT = bm.load_summary(summaryFile2,trinityGene=False,best=True)

## load deseq results
featuresDir = os.path.join(homeDir,sourceDir,"features")
deseqResultsPath = os.path.join(featuresDir,"deseq.csv")
print deseqResultsPath
deseqIds, deseqColumns, deseqMat = read_de_results(deseqResultsPath,tool='DESeq')

## create a summary table and a csv file 
outPath = os.path.join(featuresDir,'de-summary-%s.csv'%(assembly))
print 'saving %s'%outPath
fidout = open(outPath,'w')
writer = csv.writer(fidout)
writer.writerow(['transcript ID','hitId','hitNcbiId','hitSpecies','e-value','DESeq-pval','DESeq-adj-pval'])

## prepare supplment output 
#columns = [21,22,65,17,17,17]
#row = "+"
#head = "+"
#for col in columns:
#    row += "-"*col+"-+"
#    head += "="*col+"=+"

#print("\nGene sets\n_____________________")
#print(row)
#items = ['transcript ID','hitId','species','e-value','DESeq-pval','DESeq-adj-pval']
#print_rest_table_contents(columns,items,withTrailing=False)
#print(head)

if assembly == 'dn':
    threshold = 0.05
else:
    threshold = 0.5

## filter out nans        
padjInd = np.where(deseqColumns == 'padj')[0]
size1 = deseqIds.shape[0]
filter1 = np.where(~np.isnan(deseqMat[:,padjInd]))[0]
deseqIds = deseqIds[filter1]
deseqMat = deseqMat[filter1,:]

## filter for only the most significant transcripts (max 50)
threshold = 0.5
size2 = deseqIds.shape[0]
filter2 = np.where(deseqMat[:,padjInd] <= threshold)[0][:50]
deseqIds = deseqIds[filter2]
deseqMat = deseqMat[filter2,:]

## prepare supplment output        
columns = [20,22,25,50,17,17,17]
row = "+"
head = "+"
for col in columns:
    row += "-"*col+"-+"
    head += "="*col+"=+"

print("\nGene sets\n_____________________")
print(row)
items = ['transcript ID','hitId','gene-symbol','species','e-value','DESeq-pval','DESeq-adj-pval']
print_rest_table_contents(columns,items,withTrailing=False)
print(head)

## loop through transcripts based on FDR                   
rankedInds = np.argsort(deseqMat[:,np.where(deseqColumns=='padj')[0][0]])
sphinxLinks = []
for rind in rankedInds:
    deseqFDR = deseqMat[rind,np.where(deseqColumns=='padj')[0][0]]
    deseqPval = deseqMat[rind,np.where(deseqColumns=='pvalue')[0][0]]
    transcriptId = deseqIds[rind]

    if assembly == 'ref':
        newTranscriptId = None
        for tid in  [transcriptId + "." + str(i) for i in range(1,9)]:
            if not bmapXT.has_key(tid):
                continue
            if not newTranscriptId:
                newTranscriptId = (tid,bmapXT[tid][4])
            if bmapXT[tid][4] < newTranscriptId[1]:
                newTranscriptId = (tid,bmapXT[tid][4])
        if newTranscriptId:
            transcriptId = newTranscriptId[0]

    bmap,geneId = None,None
    if bmapXT.has_key(transcriptId):
        bmap = bmapXT
    elif bmapSP.has_key(transcriptId):
        bmap = bmapSP

    if bmap:
        mId = bmap[transcriptId]
        if mId[1] != '-':
            geneId = mId[1]
        elif ref2gene.has_key(mId[0]):
            geneId = ref2gene[mId[0]]

        hitId,hitNcbiId,species,speciesNcbi,_evalue = bmap[transcriptId]
        _evalue = round(float(_evalue),7)

        if len(species) > 60:
            species = species[:61]+"..."

        uniLink =  "`%s`_"%hitId

        if hitNcbiId != '-':
            uniUrl = "http://www.uniprot.org/uniprot/%s"%(hitId)
        else:
            uniUrl = "http://www.ncbi.nlm.nih.gov/gene/?term=%s"%(hitId)

        species = re.sub("\s+$","",re.sub("\(.+\)","",species))
    else:
        hitId,geneId,species,speciesNcbi,_evalue,uniUrl,uniLink = '-','-','-','-','-','-','-'

    writer.writerow([transcriptId,hitId,geneId,species,_evalue,deseqPval,deseqFDR])

    ## print info for supplement
    if deseqFDR <= threshold and rind < 200:
        if hitId != '-':
            sphinxLinks.append(".. _%s: %s"%(hitId,uniUrl))

        geneSym = '-'
        if geneId == '-':
            geneSym = '-'
        elif gene2sym.has_key(geneId):
            geneSym = gene2sym[geneId]
        else:
            geneSym = 'unmapped'
        items = [transcriptId,uniLink,geneSym,species[:48],_evalue,round(deseqPval,7),round(deseqFDR,7)]
        print_rest_table_contents(columns,items)

for sl in sphinxLinks:
    print sl

print("done")
