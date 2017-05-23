#!/usr/bin/python
"""
take the output from sailfis to get isoform count matrices
"""

__author__ = "Adam Richards"

import sys,os,getopt,csv,time,re,gc
import numpy as np

## declare variables
assembly = 'gg'
homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus","%s-trinity"%(assembly))

def load_file(sample):

    ## error check
    resultsFile = os.path.join(homeDir,'sailfish',sample,'quant.sf')
    if not os.path.exists(resultsFile):
        raise Exception("Cannot find results file %s"%resultsFile)

    ## infile
    fidin = open(resultsFile,'r')
    reader = csv.reader(fidin,delimiter="\t")

    debug = 0
    header = ['Transcript','Length','TPM','RPKM','KPKM','EstimatedNumKmers','EstimatedNumReads']
    results = {}

    for key in header:
        results[key] = []

    gc.disable()
    for linja in reader:
        if linja[0][0] == '#':
            continue
        
        results['Transcript'].append(linja[0])
        results['TPM'].append(linja[2])
        results['RPKM'].append(linja[3])
        results['EstimatedNumReads'].append(linja[6])

    gc.enable()
    return results

if __name__ == "__main__":

    sampleList = ["A", "B", "C", "D", "E", "F", "G", "H"]

    allResults = {}
    for sample in sampleList:
        allResults[sample] = load_file(sample)
    
    def write_matrix(column):
        numTranscripts = len(allResults['A']['Transcript']) 
        
        ## create a count matrix
        if column == 'EstimatedNumReads':
            outFile = os.path.join(homeDir,'features','est_counts.csv')
        else: 
            outFile = os.path.join(homeDir,'features','%s_counts.csv'%(column.lower()))
        fidout = open(outFile,'w')
        writer = csv.writer(fidout)
        writer.writerow(['transcript'] + sampleList)
        
        for row in range(numTranscripts):
            trans = allResults['A']['Transcript'][row]
            toWrite = [allResults[s][column][row] for s in sampleList]
            writer.writerow([trans] + toWrite)
    
    write_matrix('EstimatedNumReads')
    write_matrix('TPM')
    write_matrix('RPKM')

print('complete.')
