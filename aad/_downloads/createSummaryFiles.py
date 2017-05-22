#!/usr/bin/env python
"""
make some plots
http://mbni.org/customcdf/20.0.0/entrezg.download/HGU133Plus2_Hs_ENTREZG_20.0.0.zip

"""

import os,sys,csv,re,time,cPickle,subprocess
import numpy as np
from htsint.database import Gene,Taxon,db_connect
sys.path.append(os.path.join("..","lib"))
from expressionLib import *

diagDict = {'healthy':1,'asthma':2}
sexDict = {'female':1,'male':2}
mappers = get_gene_mappers()
runAll = False
removeAll = True
outliers = True

if outliers:
    dataPickle = 'all-probes-o.pickle'
    gdataPickle = 'all-genes-o.pickle'
else:
    dataPickle = 'all-probes.pickle'
    gdataPickle = 'all-genes.pickle'

if removeAll:
    if os.path.exists(dataPickle):
        os.remove(dataPickle)
    if os.path.exists(gdataPickle):
        os.remove(gdataPickle)
    
def read_data(gseId):
    print("...reading %s"%gseId)
    
    if runAll:

        ## remove old
        if gseId in ["GSE51392","GSE8190","unpublished","yang16"]:
            experimentId = gseId
        else:
            experimentId = gseId + "-1"

        if outliers:
            experimentId += "-o"
                
        file1 = os.path.join("data",experimentId+"-exp.csv")
        file2 = os.path.join("data",experimentId+"-cov.csv")
        file3 = os.path.join("data",experimentId+"-probes.tab")
        file4 = os.path.join("data",experimentId+"-probes-parsed.tab")

        for f in [file1,file2,file3,file4]:
            if os.path.exists(f):
                print("...removing %s"%f)
                os.remove(f)
                
        ## run the R scripts
        if gseId == 'unpublished':
            cmd = 'Rscript runUnpublished.R'
        elif gseId == 'yang16':
            cmd = 'Rscript runYang16.R'
        else:
            cmd = 'Rscript run%s.R'%gseId

        if outliers:
            cmd = cmd + " -o"

        print(cmd)
        os.system(cmd)
        print("limma finished.")
        
        if gseId in ["GSE51392","GSE8190"]:
            cmd = 'python parse%s.py'%gseId

            if outliers:
                cmd = cmd + " -o"
                
            print(cmd)
            os.system(cmd)
        elif gseId == "unpublished":
            cmd = 'python parseUnpublished.py'

            if outliers:
                cmd = cmd + " -o"
        
            print(cmd)
            os.system(cmd)
        elif gseId == "yang16":
            cmd = 'python parseYang16.py'

            if outliers:
                cmd = cmd + " -o"
        
            print(cmd)
            os.system(cmd)
            
        else:
            print("...skipping parse")

            
    if gseId in ["GSE8190","GSE51392","unpublished","yang16"]:
        expId = gseId
        if outliers:
            resultsFile = os.path.join("results","%s-o-probes-parsed.tab"%expId)
        else:
            resultsFile = os.path.join("results","%s-probes-parsed.tab"%expId)
    else:
        expId = gseId + "-" + "1"
        if outliers:
            resultsFile = os.path.join("results","%s-o-probes.tab"%expId)
        else:
            resultsFile = os.path.join("results","%s-probes.tab"%expId)
            
    if outliers:
        expFile = os.path.join("data","%s-o-exp.csv"%expId)
        covsFile = os.path.join("data","%s-o-cov.csv"%expId)
    else:
        expFile = os.path.join("data","%s-exp.csv"%expId)
        covsFile = os.path.join("data","%s-cov.csv"%expId)
        
    if not os.path.exists(resultsFile):
        print("... %s"%resultsFile)
        raise Exception("cannot find results file \n %s"%resultsFile)

    probeIds, sampleIds, mat = read_exprs_data(expFile)
    mat = mat.T
    covs = read_generic(covsFile)
    results = read_limma_results(resultsFile)
    
    ## update out of date genes
    if gseId in ["GSE67472"]:
        results['gene-id'] = results["SPOT_ID"]
        results['probe-id'] = results["ID"]
    elif gseId in ['GSE43696']:
        results['gene-id'] = results["GENE"]
        results['probe-id'] = results["SPOT_ID"]
    elif gseId in ["GSE18965"]:
        results['gene-id'] = results["ENTREZ_GENE_ID"]
        results['probe-id'] = results["ID"]
        covs['samples'] = covs['sample']
    elif gseId in ["GSE19187"]:
        results['gene-id'] = [g.split("//")[-1] for g in results['gene_assignment']]
        results['gene-id'] = np.array([re.sub("\s+","",g) for g in results['gene-id']])
        results['probe-id'] = results["ID"]
    elif gseId in ["GSE8190","unpublished","yang16"]:
        results['gene-id'] = results["GeneName"]
        results['probe-id'] = results["probes"]
        if gseId == 'unpublished':
            covs['samples'] =  np.array(["X" + re.sub("\.txt","",s) for s in covs['FileName']])
            covs['samples'] =  np.array([re.sub("-",".",s) for s in covs['samples']])
        elif gseId == 'GSE8190':
            pass
        elif gseId == 'yang16':
            covs['samples'] =  np.array([re.sub("\.txt","",s) for s in covs['FileName']])
        else:
            covs['samples'] = np.array([s +".txt" for s in covs['sample']])
    elif gseId in ["GSE51392"]:
        results['gene-id'] = results["GeneName"]
        covs['samples'] = covs['filename']
        results['probe-id'] = results["probes"]
    else:
        for key,item in results.iteritems():
            print key, item[:5]
        raise Exception("need to specify genes")

    totalProbes = results['gene-id'].size

    ## take the first possibility
    results['gene-id'] = np.array([re.sub("\s+","",gene.split("///")[0]) for gene in results['gene-id']])

    ## update some old gene names
    genesToUpdate = {'83954':'387316','79614':'4883','646329':'378805','645460':'6426','645090':'286464',
                     '55410':'83869','149469':'23334','286367':'55335'}

    for key,item in genesToUpdate.iteritems():
        idx = np.where(results['gene-id'] == key)[0]
        results['gene-id'][idx] = item
    
    ## remove genes that are discontinued
    _discontinuedGenes = ['8475','84796',"100131298",'100507457','100505490','100505657','100506972','100507315',
                          '100507322','100507780']
    discontinuedGenes = [d for d in _discontinuedGenes]
    for g in results['gene-id']:
        if not mappers['symbol'].has_key(g):
            discontinuedGenes.append(g)
            
    print("...removing %s/%s discontinued genes"%(len(discontinuedGenes),results['gene-id'].size))
    dgIndices = set([])
    for dg in discontinuedGenes:
        inds = np.where(results["gene-id"] == dg)[0]
        if inds.size == 0:
            continue
        dgIndices.update(inds.tolist())
    dgIndices = np.array(list(dgIndices))
    sigIndices = np.where(results["adj.P.Val"]<0.1)[0]
    considerUpdating = np.intersect1d(dgIndices,sigIndices)
    if considerUpdating.size > 0:
        print("......there are %s you many want to consider updating"%considerUpdating.size)
        for indx in considerUpdating:            
            if results["gene-id"][indx] in _discontinuedGenes +['','-']:
                results["gene-id"][indx] = 'unmapped'
                continue

            print(".........%s (%s)"%(results["gene-id"][indx],results["adj.P.Val"][indx]))
            
    toKeep = np.setdiff1d(np.arange(0,totalProbes),dgIndices)
    unknownGenes = set([])
    for ug in ['','unknown','unmapped','NA','-']:
        inds = np.where(results['gene-id'] == ug)[0]
        if inds.size > 0:
            unknownGenes.update(inds.tolist())
    unknownGenes = np.array(list(unknownGenes))
    print("......removing %s unknown genes"%unknownGenes.size)
    toKeep = np.setdiff1d(toKeep,unknownGenes)
            
    for key,item in results.iteritems():
        results[key] = item[toKeep]
    
    ## reorder the matrix probes and samples to match results
    print "probeIds", probeIds[:5]
    for key,item in results.iteritems():
        print key, item[:4]
    probeInds = np.array([np.where(probeIds==p)[0][0] for p in results['probe-id']])
    print "reordering...", probeInds[1:5]
    probeIds = probeIds[probeInds]
    mat = mat[:,probeInds]

    if not covs.has_key('samples'):
        for key,item in covs.iteritems():
            print key, item[:5]
        raise Exception("set the samples column to 'samples' in covs")

    if not np.array_equal(covs['samples'],sampleIds):
        print("%s,%s"%(covs['samples'][0],sampleIds[0]))
        raise Exception("sample ids do not match")
    
    return {'p-value':results['P.Value'],
            'logFC':results['logFC'],
            'adj-p-value':results['adj.P.Val'],
            'probe':results['probe-id'],
            'gene':results['gene-id'],
            "mat":mat,
            "covs":covs}

## load data from all studies
dataSets = {"GSE67472":"Christenson",
            "GSE43696":"Voraphani",
            "GSE19187":"Giovannini-Chami",
            "GSE18965":"Kicic",
            "GSE8190":"Yang",
            "yang16":"yang16",
            "unpublished":"unpublished",
            "GSE51392":"Wagener"}

#dataSets = {"yang16":"yang16"}

## use a pickle file to speed things up
if not os.path.exists(dataPickle):
    print("...saving probes pickle")
    data = {}    
    for key in dataSets.iterkeys():
        data[key]=read_data(key)

    tmp = open(dataPickle,'w')
    cPickle.dump(data,tmp)
    tmp.close()
else:    
    print("...loading probes pickle")
    tmp = open(dataPickle,'r')
    data = cPickle.load(tmp)
    tmp.close()

## get master list of gene ids
allGenes = set([])
gData = {}
for study in dataSets.keys():
    allGenes.update(data[study]['gene'].tolist())
allGenes = np.array(list(allGenes))

## combine probes
print("combining probe-level data...")

if not os.path.exists(gdataPickle):
    print("...saving genes pickle")
    gdata = {}
    for study in dataSets.keys():
        data[study]['logFC'] = data[study]['logFC'].astype('float')
        samples = data[study]['covs']['samples'].size
        uniqueGenes = np.unique(data[study]['gene'])
        
        gdata[study] = {'p-value':[],'logFC':[],'adj-p-value':[],'gene':[],
                        'mat':np.zeros((samples,uniqueGenes.size),),
                        'covs':data[study]['covs']}
        print("%s\t %s/%s"%(study,uniqueGenes.size,data[study]['gene'].size))
        for g,gene in enumerate(uniqueGenes):
            geneIndices = np.where(data[study]['gene'] == gene)[0]
            gdata[study]['gene'].append(gene)
            gdata[study]['p-value'].append(np.median(data[study]['p-value'][geneIndices]))
            gdata[study]['adj-p-value'].append(np.median(data[study]['adj-p-value'][geneIndices]))
            gdata[study]['logFC'].append(np.median(data[study]['logFC'][geneIndices]))
            gdata[study]['mat'][:,g] = np.array([np.median(data[study]['mat'][s,geneIndices]) for s in range(data[study]['mat'].shape[0])])

        for key in ['gene','p-value','adj-p-value','logFC']:
            gdata[study][key] = np.array(gdata[study][key])        
    tmp = open(gdataPickle,'w')
    cPickle.dump(gdata,tmp)
    tmp.close()
else:    
    print("...loading genes pickle")
    tmp = open(gdataPickle,'r')
    gdata = cPickle.load(tmp)
    tmp.close()
    
## scan each study for results and assemble uniqueGene based lists
studyKeys = sorted(dataSets.keys())
studyNames = [dataSets[s] for s in studyKeys]
pvals,apvals,fcs,exprs = {},{},{},{}
dxHeader = []
studyHeader = []
genderHeader = []
ageHeader = []
sampleSizes = {}

for study in studyKeys:    
    if gdata[study]['covs'].has_key('dx'):
        dx = 'dx'
    elif gdata[study]['covs'].has_key('disease_state'):
        dx = 'disease_state'
    else:
        dx = 'phenotype'

    dxSamples = gdata[study]['covs'][dx]
    genderSamples = gdata[study]['covs']['gender']
    ageSamples = gdata[study]['covs']['age']
    if 'healthy' in dxSamples:
        dxSamples[np.where(dxSamples=='healthy')[0]] = 'control'

    print study, dx, np.unique(gdata[study]['covs'][dx]),dxSamples.size,'healthy'in dxSamples
    dxHeader.extend(dxSamples.tolist())
    studyHeader.extend([study]*dxSamples.size)
    genderHeader.extend(genderSamples)
    ageHeader.extend(ageSamples)
    sampleSizes[study] = dxSamples.size
    
    pvals[study] = []
    apvals[study] = []
    fcs[study] = []
    exprs[study] = []

for i in range(len(dxHeader)):
    genderHeader[i] = genderHeader[i].lower()
    if genderHeader[i] == 'f':
        genderHeader[i] = 'female'
    elif genderHeader[i] == 'm':
        genderHeader[i] = 'male'

for g,gene in enumerate(allGenes):
    if g%2000==0:
        print("%s/%s"%(g,allGenes.size))
    for study in studyKeys:
        indx = np.where(gdata[study]['gene']==gene)[0]
        if len(indx) > 0:
            indx = indx[0]
            pval = gdata[study]['p-value'][indx]
            apval = gdata[study]['adj-p-value'][indx]
            fc = gdata[study]['logFC'][indx]
            exp = gdata[study]['mat'][:,indx]
        else:
            pval,apval,fc,exp = 'NA','NA','NA',['NA']*sampleSizes[study]

        pvals[study].append(pval)
        apvals[study].append(apval)
        fcs[study].append(fc)
        exprs[study].append(exp)
            
## create the summary files
if outliers:
    fid1 = open("./data/merged-pvalues-o.csv","w")
else:
    fid1 = open("./data/merged-pvalues.csv","w")
writer1 = csv.writer(fid1)
writer1.writerow(['gene','symbol']+studyNames)

if outliers:
    fid2 = open("./data/merged-qvalues-o.csv","w")
else:
    fid2 = open("./data/merged-qvalues.csv","w")
writer2 = csv.writer(fid2)
writer2.writerow(['gene','symbol']+studyNames)

if outliers:
    fid3 = open("./data/merged-logfc-o.csv","w")
else:
    fid3 = open("./data/merged-logfc.csv","w")
writer3 = csv.writer(fid3)
writer3.writerow(['gene','symbol']+studyNames)

if outliers:
    fid4 = open("./data/merged-exp-o.csv","w")
else:
    fid4 = open("./data/merged-exp.csv","w")
writer4 = csv.writer(fid4)

writer4.writerow(['gene','symbol']+dxHeader)
writer4.writerow(['gene','symbol']+studyHeader)
writer4.writerow(['gene','symbol']+ageHeader)
writer4.writerow(['gene','symbol']+genderHeader)

for g,gene in enumerate(allGenes):
    symbol = mappers['symbol'][gene]
    writer1.writerow([gene,symbol] + [pvals[study][g] for study in studyKeys])
    writer2.writerow([gene,symbol] + [apvals[study][g] for study in studyKeys])
    writer3.writerow([gene,symbol] + [fcs[study][g] for study in studyKeys])

    towrite = []
    for study in studyKeys:
        towrite.extend(exprs[study][g])

    if len(towrite) != len(dxHeader):
        print towrite
        print len(towrite), len(dxHeader)
        raise Exception("invalid line to be written")
    writer4.writerow([gene,symbol] + towrite)
    
print("files written for %s genes"%allGenes.size)
