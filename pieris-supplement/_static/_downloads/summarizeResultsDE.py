#!/usr/bin/env python
"""
read the output of a DESeq2 analysis and summarize the results
create a new output file with associated gene names and go terms
"""

import os,sys,csv,getopt,re
import numpy as np
from sqlalchemy.sql import select
from htsint.database import db_connect,Taxon,Gene
from htsint.tools import read_matrix,read_de_results,print_rest_table_contents
from htsint.blast import BlastMapper

#assembly = 'dn'
threshold = 0.05
homeDir = os.path.join(os.path.expanduser("~"),"sequencing","pieris")
transcript = 'isoforms'
featuresDir = os.path.join(homeDir,"features")

## load the differential expression results
evalue = 0.00001
deseqResultsPath = os.path.join(featuresDir,"deseq.csv")
deseqMatIds, deseqMatColumns, deseqMat = read_de_results(deseqResultsPath,tool='DESeq')

## create a summary table and a csv file
outPath = os.path.join(featuresDir,'de-summary-%s.csv'%(transcript))
fidout = open(outPath,'w')
writer = csv.writer(fidout)
writer.writerow(['transcript ID','hitId','hitNcbiId','hitSpecies','e-value','DESeq-pval','DESeq-adj-pval'])

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
s = select([Taxon.id,Taxon.ncbi_id,Taxon.name]).where(Taxon.ncbi_id.in_(['7227']))
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

## load the blast map
bm = BlastMapper()
summaryFile1 = os.path.join(homeDir,"dn-trinity",'blast-dn-parsed_summary.csv') 
summaryFile2 = os.path.join(homeDir,"dn-trinity",'blast-dm-parsed_summary.csv') 
summaryFile3 = os.path.join(homeDir,"dn-trinity","blast-mc-parsed_summary.csv")
summaryFile4 = os.path.join(homeDir,"dn-trinity",'blast-dp-parsed_summary.csv')

bmapSP = bm.load_summary(summaryFile1,trinityGene=False,best=True)
bmapDM = bm.load_summary(summaryFile2,trinityGene=False,best=True)
bmapMC = bm.load_summary(summaryFile3,trinityGene=False,best=True)
bmapDP = bm.load_summary(summaryFile4,trinityGene=False,best=True)

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
rankedInds = np.argsort(deseqMat[:,np.where(deseqMatColumns=='padj')[0][0]])
sphinxLinks = []
finalCount = 0
debug = 0
for rind in rankedInds:
    deseqFDR = deseqMat[rind,np.where(deseqMatColumns=='padj')[0][0]]
    deseqPval = deseqMat[rind,np.where(deseqMatColumns=='pvalue')[0][0]]
    transcriptId = deseqMatIds[rind]

    #if transcript == 'gene':
    #     transcriptId = re.sub("_i\d+","",transcriptId)

    bmap,geneId = None,None
    if bmapDM.has_key(transcriptId):
        bmap = bmapDM
    elif bmapSP.has_key(transcriptId):
        bmap = bmapSP
    elif bmapMC.has_key(transcriptId):
        bmap = bmapMC
    elif bmapDP.has_key(transcriptId):
        bmap = bmapDP
    
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
    if deseqFDR <= threshold:# and rind < 200:
        if bmap:
            finalCount += 1

        if hitId != '-':
            sphinxLinks.append(".. _%s: %s"%(hitId,uniUrl))
        debug += 1
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

print("final count:%s"%(finalCount))
print("debug:%s"%(debug))
print('done')
