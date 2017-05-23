#!/usr/bin/python
"""
tutorials and docs

Before we run the trinity assebmly we unzip the files an trim them
If the data are backed up we can remove the gz files to save space

from the trimmomatic documents
http://www.usadellab.org/cms/index.php?page=trimmomatic

java -jar <path to trimmomatic.jar> PE [-threads <threads] [-phred33 | -phred64] \
[-trimlog <logFile>] <input 1> <input 2> <paired output 1> <unpaired output 1> \ 
<paired output 2> <unpaired output 2> <step 1> ... 

PE [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...

"""

import os,sys,re,shutil
import HTSeq
from htsint import run_subprocess

## locations
seqDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus-data")
analysisDir = os.path.join(seqDir,"analyses")
if not os.path.isdir(analysisDir):
    os.mkdir(analysisDir)

## functions 
def unzip_file(source,target):
    print('unzipping...%s'%source)
    cmd = "gunzip -c %s > %s"%(source,target)
    print cmd
    run_subprocess(cmd)

def assemble_reads(sampleList):
    ## get the files
    def get_fq_files(sample):
        sampleDir = os.path.join(seqDir,"Sample_%s"%sample,"raw_illumina_reads")
        leftFile,rightFile = None,None
        for fileName in os.listdir(sampleDir):
            if re.search("^%s.*\.fastq$"%sample,fileName) and re.search("R1",fileName):
                leftFile = os.path.realpath(os.path.join(sampleDir,fileName))
            elif re.search("^%s.*.fastq$"%sample,fileName) and re.search("R2",fileName):
                rightFile = os.path.realpath(os.path.join(sampleDir,fileName))
        return leftFile,rightFile

    for sample in sampleList:
        print ("\nprocessing sample: %s"%sample)
        left_fq,right_fq = get_fq_files(sample)
        
        ## run trimmomatic (assuming FR)
        out1 = os.path.join(analysisDir,"%s_left_paired.fq"%sample)
        out2 = os.path.join(analysisDir,"%s_left_unpaired.fq"%sample)
        out3 = os.path.join(analysisDir,"%s_right_paired.fq"%sample)
        out4 = os.path.join(analysisDir,"%s_right_unpaired.fq"%sample)
        
        cmd = "java -jar /usr/share/java/trimmomatic-0.32.jar PE -threads 29 -phred33 " +\
              "%s %s "%(left_fq,right_fq) +\
              "%s %s "%(out1,out2)+\
              "%s %s "%(out3,out4)+\
              "ILLUMINACLIP:/usr/src/trinityrnaseq_r20140717/trinity-plugins/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:10 "+\
              "LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36"

        if not os.path.exists(out1):
            run_subprocess(cmd)

        ## remove 
        #for fileName in faFileList:
        #    print("removing %s"%fileName
        #    os.remove(fileName)

    print "complete"
 
if __name__ == "__main__":

    sampleList = ["A", "B", "C", "D", "E", "F", "G", "H"]
    assemble_reads(sampleList)
