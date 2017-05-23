#!/usr/bin/python
# tutorials and docs
# http://www-huber.embl.de/users/anders/HTSeq/doc/tour.html



import os,sys,csv
import numpy as np
import HTSeq
from htsint import run_subprocess

## locations
seqDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus","reference")
starDir = os.path.join(seqDir,'star')
analysisDir = os.path.join(seqDir,"features")
if not os.path.isdir(analysisDir):
    os.mkdir(analysisDir)

def get_count_matrix(sample,overwrite=False):
    print("getting counts for sample %s"%sample)

    gtfFilePath  = os.path.realpath(os.path.join(seqDir,"Xentr7_2_Stable.gtf"))
    samFilePath = os.path.join(starDir,"%s_aligned_sorted.sam"%sample) 
    countsFilePath = os.path.join(analysisDir,'counts_%s.txt'%sample.lower())

    if not os.path.exists(samFilePath):
        raise Exception("cannot find sam file %s"%(samFilePath))

    if os.path.exists(countsFilePath) and overwrite:
        print("...removing old count file '%s'"%countsFilePath)
        os.remove(countsFilePath)
    elif os.path.exists(countsFilePath) and not overwrite:
        print("...'%s' exists use 'overwrite' to remove"%countsFilePath)
        return
    
    ## save the count matrix
    print('...running htseq-count')
    cmd = "/usr/bin/htseq-count %s %s > %s"%(samFilePath,gtfFilePath,countsFilePath)
    run_subprocess(cmd)


if __name__ == "__main__":
    sampleList = ["A","B","C","D","E","F","G","H"]

    ## write the individual count files
    overwrite = False
    for sample in sampleList:
        get_count_matrix(sample,overwrite=overwrite)

    ## create an empty matrix to hold counts
    fid = open(os.path.join(analysisDir,'counts_a.txt'),'r')
    reader = csv.reader(fid,delimiter="\t")
    lineCount = 0
    for linja in reader:
        lineCount += 1
    fid.close()
    mat = np.zeros((lineCount,len(sampleList)),dtype=int)
    rows = []

    ## load the samples and populate matrix
    for s,sample in enumerate(sampleList):
        fid = open(os.path.join(analysisDir,'counts_%s.txt'%sample.lower()),'r')
        reader = csv.reader(fid,delimiter="\t")
        print("\nadding sample %s"%sample)
        lineCount = 0 
        for linja in reader:
            if s==0:
                rows.append(linja[0])

            mat[lineCount,s] = int(linja[1])
            lineCount +=1
            
            if lineCount == 1:
                print 'debug',linja

        fid.close()

    rows = np.array(rows)

    ## save raw and filtered matrices to file
    print('writing matrix to file...')
    fidout1 = open(os.path.join(analysisDir,"ref_raw_counts_matrix.csv"),'w')
    writer1 = csv.writer(fidout1)
    writer1.writerow(['id']+sampleList)
    fidout2 = open(os.path.join(analysisDir,"ref_filtered_counts_matrix.csv"),'w')
    writer2 = csv.writer(fidout2)
    writer2.writerow(['id']+sampleList)

    for r in range(mat.shape[0]):
        writer1.writerow([rows[r]] + mat[r,:].tolist())

        if mat[r,:].sum() > 0:
            writer2.writerow([rows[r]] + mat[r,:].tolist())

    ## clean up
    fidout1.close()
    fidout2.close()
