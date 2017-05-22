#!/usr/bin/env python
"""
make some plots

http://mbni.org/customcdf/20.0.0/entrezg.download/HGU133Plus2_Hs_ENTREZG_20.0.0.zip

"""
import os,sys,csv,re,time,cPickle
import numpy as np
import matplotlib.pyplot as plt
from sklearn import preprocessing
import matplotlib as mpl
sys.path.append(os.path.join("..","lib"))
from expressionLib import *
from scipy import stats as sstats
from matplotlib_venn import venn3, venn3_circles
mpl.rc('font', family='sans-serif')
mpl.rc('font', size=8)

## setup rpy
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
stats = importr('stats')

dataSets = {"GSE67472":"Christenson",
            "GSE43696":"Voraphani",
            "GSE19187":"Giovannini-Chami",
            "GSE18965":"Kicic",
            "GSE8190":"Yang",
            "unpublished":"unpublished",
            "GSE51392":"Wagener",
            "yang16":"yang16"}

#allStudies = np.array(['Giovannini-Chami', 'Kicic', 'Yang', 'Voraphani', 'unpublished','Wagener', 'Christenson','yang16'])
outfile = "significance-summary.csv"
header = ["study","samples","missing","outliers","invivo","0.05","0.1"]
writerSummary = csv.writer(open(outfile,'w'))
writerSummary.writerow(header)

for outliers in [True,False]:

    ## variables
    if outliers:
        probesPickle = 'all-probes-o.pickle'
        genesPickle = 'all-genes-o.pickle'
    else:
        probesPickle = 'all-probes.pickle'
        genesPickle = 'all-genes.pickle'

    if outliers:
        pvalsFile = os.path.join(".","data","merged-pvalues-o.csv")
        pvals = read_generic(pvalsFile)
        qvalsFile = os.path.join(".","data","merged-qvalues-o.csv")
        qvals = read_generic(qvalsFile)
        logFCsFile = os.path.join(".","data","merged-logfc-o.csv")
        lfcs = read_generic(logFCsFile)
    else:
        pvalsFile = os.path.join(".","data","merged-pvalues.csv")
        pvals = read_generic(pvalsFile)
        qvalsFile = os.path.join(".","data","merged-qvalues.csv")
        qvals = read_generic(qvalsFile)
        logFCsFile = os.path.join(".","data","merged-logfc.csv")
        lfcs = read_generic(logFCsFile)

    ## load gene level data
    print("...loading genes pickle")
    tmp = open(genesPickle,'r')
    gdata = cPickle.load(tmp)
    tmp.close()
    
    for invivo in [True,False]:
        if invivo == False:
            studies = ['Giovannini-Chami', 'Kicic', 'Yang', 'Voraphani', 'unpublished','Wagener', 'Christenson','yang16']
        else:
            studies = ['Giovannini-Chami', 'Yang', 'Voraphani', 'unpublished', 'Christenson','yang16']

        if not os.path.exists(genesPickle) or not os.path.exists(probesPickle):
            raise Exception("Run createSummaryFiles.py first")

        ## print summary info
        print("\n----------------------")
        significant,significant2 = {},{}
        sampleSizes = {}
        for dset,study in dataSets.iteritems():
            print(study)
            #print gdata[dset]['covs'].keys()
            if study == 'unpublished':
                dx = gdata[dset]['covs']['phenotype']
            elif study in ['Voraphani','Giovannini-Chami','Christenson']:
                dx = gdata[dset]['covs']['disease_state']
            else:
                dx = gdata[dset]['covs']['dx']

            sampleSizes[study] = dx.size
            print("\t"+";".join(["%s--%s"%(d,np.where(dx==d)[0].size) for d in np.unique(dx)]))
            print("\tmissing: %s"%np.where(qvals[study]=='NA')[0].size)
            sigIndices = np.where(qvals[study][np.where(qvals[study]!='NA')[0]].astype(float) < 0.05)[0]
            sigIndices2 = np.where(qvals[study][np.where(qvals[study]!='NA')[0]].astype(float) < 0.1)[0]
            sigGenes = qvals['symbol'][np.where(qvals[study]!='NA')[0]][sigIndices]
            sigGenes2 = qvals['symbol'][np.where(qvals[study]!='NA')[0]][sigIndices2]
        
            for sg in sigGenes:
                if not significant.has_key(sg):
                    significant[sg] = 0
                significant[sg] += 1

            for sg in sigGenes2:
                if not significant2.has_key(sg):
                    significant2[sg] = 0
                significant2[sg] += 1
   
            print("\tqval < 0.05: %s"%sigIndices.size)
            print("\tqval < 0.1: %s"%sigIndices2.size)    
            writerSummary.writerow([study,dx.size,np.where(qvals[study]=='NA')[0].size,str(outliers),str(invivo),sigIndices.size,sigIndices2.size])
            
        print("\n----------------------")
        sigGenes = np.array(significant.keys())
        sigCounts = np.array(significant.values())
        sigGenes2 = np.array(significant2.keys())
        sigCounts2 = np.array(significant2.values())

        if outliers and invivo:
            fid = open("./results/meta-analysis-invivo-outlier.csv","w")
        elif outliers and invivo == False:
            fid = open("./results/meta-analysis-noninvivo-outlier.csv","w")
        elif outliers == False and invivo == True:
            fid = open("./results/meta-analysis-invivo-nonoutlier.csv","w")
        elif outliers == False and invivo == False:
            fid = open("./results/meta-analysis-noninvivo-nonoutlier.csv","w")
        else:
            raise Exception("should not happen")

        writer = csv.writer(fid)
        writer.writerow(["gene","symbol","Stouffer-Z","Stouffer-pvalue","Stouffer-adjusted-pvalue"])

        ## calculate weighted Z-scores
        studyList = np.array(studies)
        stoufferZs,stoufferPvals = [],[]

        for g,gene in enumerate(qvals['gene']):
            symbol = qvals['symbol'][g]
            _pvalues = np.array([qvals[study][g] for study in studyList])
            nonMissing = np.where(_pvalues != 'NA')[0]
            pvalues = _pvalues[nonMissing].astype('float')
            studies = studyList[nonMissing]
            n = np.array([sampleSizes[s] for s in studies])
            Z,sP = sstats.combine_pvalues(pvalues,weights=n)
            stoufferZs.append(Z)
            stoufferPvals.append(sP)

        stoufferQvals = stats.p_adjust(FloatVector(stoufferPvals), method = 'BH')
        orderedInds = np.argsort(stoufferQvals)

        for i in orderedInds:
            gene = qvals['gene'][i]
            symbol = qvals['symbol'][i]
            sz = stoufferZs[i]
            spval = stoufferPvals[i]
            sqval = stoufferQvals[i]
            writer.writerow([gene,symbol,sz,spval,sqval])
        
        fid.close()
    
        print 'done.'

print("all done")

"""
## create venn diagrams
def create_venn(studies):
    print("...plotting %s"%";".join(studies))
    fig = plt.figure(figsize=(6,6))
    study1,study2,study3 = studies
    study1Inds = np.where(qvals[study1][np.where(qvals[study1]!='NA')[0]].astype(float) < 0.1)[0]
    study1Genes = set(qvals['symbol'][np.where(qvals[study1]!='NA')[0]][study1Inds].tolist())
    study2Inds = np.where(qvals[study2][np.where(qvals[study2]!='NA')[0]].astype(float) < 0.1)[0]
    study2Genes = set(qvals['symbol'][np.where(qvals[study2]!='NA')[0]][study2Inds].tolist())
    study3Inds = np.where(qvals[study3][np.where(qvals[study3]!='NA')[0]].astype(float) < 0.1)[0]
    study3Genes = set(qvals['symbol'][np.where(qvals[study3]!='NA')[0]][study3Inds].tolist())
    v = venn3([study1Genes,study2Genes,study3Genes],[re.sub("-","",study1),re.sub("-","",study2),re.sub("-","",study3)])
    figName = './figs/venn-%s-%s-%s-0.1.png'%(study1,study2,study3)
    print("...saving %s"%figName)
    plt.savefig(figName)    

create_venn(['Voraphani','Giovannini-Chami','Christenson'])
create_venn(['Voraphani','Yang','Christenson'])
create_venn(['Giovannini-Chami','Yang','Christenson'])
create_venn(['Voraphani','unpublished','Yang',])
create_venn(['Voraphani','Giovannini-Chami','unpublished',])
create_venn(['Voraphani','unpublished','Christenson'])

print("There are %s unique significant genes"%(sigGenes.size),sigGenes2.size)
for count in range(1,len(studies)+1):
    scount1 = np.where(sigCounts >= count)[0].size
    scount2 = np.where(sigCounts2 >= count)[0].size
    print("There are %s,%s significant genes present at least %s studies"%(scount1,scount2,count))

print("Top three: %s "%(";".join(sigGenes[np.where(sigCounts >= 4)[0]])))
print("Top 17: %s "%(";".join(sigGenes[np.where(sigCounts >= 3)[0]])))
    
print("\n----------------------")
print("...leaving out studies to see union")
def get_indices_in_common(studyList):
    inCommon = []
    for g,gene in enumerate(qvals['gene']):
        missing = np.where(np.array([qvals[study][g] for study in studyList]) == 'NA')[0]
        if len(missing) == 0:
            inCommon.append(g)
    
    return(inCommon)

for studyList in [['Giovannini-Chami', 'Kicic', 'Yang', 'Voraphani', 'unpublished','Wagener', 'Christenson'],
                  ['Kicic', 'Yang', 'Voraphani', 'unpublished','Wagener', 'Christenson'],
                  ['Giovannini-Chami','Yang', 'Voraphani', 'unpublished','Wagener', 'Christenson'],
                  ['Giovannini-Chami', 'Kicic', 'Voraphani', 'unpublished','Wagener', 'Christenson'],
                  ['Giovannini-Chami', 'Kicic', 'Yang', 'unpublished','Wagener', 'Christenson'],
                  ['Giovannini-Chami', 'Kicic', 'Yang', 'Voraphani', 'Wagener', 'Christenson'],
                  ['Giovannini-Chami', 'Kicic', 'Yang', 'Voraphani', 'unpublished', 'Christenson'],
                  ['Giovannini-Chami', 'Kicic', 'Yang', 'Voraphani', 'unpublished','Wagener'],
                  ['Giovannini-Chami', 'Kicic', 'Yang', 'Voraphani', 'Christenson']]:

    indices = get_indices_in_common(studyList)
    print ";".join(studyList),len(indices)
"""
