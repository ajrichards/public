#!/usr/bin/python
"""

"""

import os,csv,sys,re,getopt
import numpy as np
from htsint import run_subprocess

## specify the locations
homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus")
readsDir = os.path.join(homeDir,'reads')
    
## more variables
email = "adam.richards@ecoex-moulis.cnrs.fr"

## run edgeR and DESeq2 
def run_deseq(countsPath,outFile):
    cmd = "Rscript runDESeq.R %s %s"%(countsPath,outFile)
    print("running...\n%s"%cmd)
    run_subprocess(cmd)
    
def create_filtered(countsPath):
    if not os.path.exists(countsPath):
        raise Exception("Cannot find counts path")
    
    filteredCountsPath = re.sub("\.csv","_filtered.csv",countsPath)
    fid1 = open(countsPath,'r')
    fid2 = open(filteredCountsPath,'w')
    reader = csv.reader(fid1)
    writer = csv.writer(fid2)
    header = reader.next()
    writer.writerow(header)

    for linja in reader:
        if np.array([float(i) for i in linja[1:]]).sum() > 1:
            writer.writerow(linja)
    fid1.close()
    fid2.close()
    return filteredCountsPath

if __name__ == "__main__":

    print("all filtered matrices are saved.")
    for method in ['dn','gg']:
        featuresDir = os.path.join(homeDir,"%s-trinity"%method,"features")
        countsPath = os.path.join(featuresDir,"est_counts.csv")
        filteredCountsPath = create_filtered(countsPath)
        outFile = os.path.join(featuresDir,"deseq.csv")
        run_deseq(filteredCountsPath,outFile)

    ## the reference assembly
    featuresDir = os.path.join(homeDir,"reference","features")
    countsPath = os.path.join(featuresDir,'raw_counts_matrix.csv')
    filteredCountsPath = create_filtered(countsPath)
    outFile = os.path.join(featuresDir,"deseq.csv")
    run_deseq(filteredCountsPath,outFile)
