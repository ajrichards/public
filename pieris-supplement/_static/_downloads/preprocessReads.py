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
from htsint import run_subprocess

## locations
seqDir = os.path.join(os.path.expanduser("~"),"sequencing","pieris")
analysisDir = os.path.join(seqDir,"reads")
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
    def get_gz_files(sample,laneDir):
        leftFile,rightFile = None,None
        for fileName in os.listdir(laneDir):
            if re.search("^%s.*\.fastq.gz$"%sample,fileName) and re.search("R1",fileName):
                leftFile = os.path.realpath(os.path.join(laneDir,fileName))
            elif re.search("^%s.*.fastq.gz$"%sample,fileName) and re.search("R2",fileName):
                rightFile = os.path.realpath(os.path.join(laneDir,fileName))
        return leftFile,rightFile

    lane1Dir = os.path.join(seqDir,"pieris-lane-1","RawData")
    lane2Dir = os.path.join(seqDir,"pieris-lane-2","RawData")

    for sample in sampleList:
        print ("\nprocessing sample: %s"%sample)
        left_gz1,right_gz1 = get_gz_files(sample,lane1Dir)
        left_gz2,right_gz2 = get_gz_files(sample,lane2Dir)

        ## unzip files into analysis dir
        for source in [left_gz1,left_gz2,right_gz1,right_gz2]:
            target = os.path.join(analysisDir,os.path.split(source)[-1][:-3])
            if not os.path.exists(target):
                unzip_file(source,target)

        ## concat files into analysis dir
        faFileList = [os.path.join(analysisDir,os.path.split(source)[-1][:-3]) for f in [left_gz1,left_gz2,right_gz1,right_gz2]]
        allLeftFile = os.path.join(analysisDir,"%s_left.fq"%sample)
        allRightFile = os.path.join(analysisDir,"%s_right.fq"%sample)
        catLeft = "cat %s %s > %s"%(faFileList[0],faFileList[1],allLeftFile)
        catRight = "cat %s %s > %s"%(faFileList[2],faFileList[3],allRightFile)
        
        print("concatenating lanes...")
        if not os.path.exists(allLeftFile):
            print catLeft
            run_subprocess(catLeft)
        if not os.path.exists(allRightFile):
            print catRight
            run_subprocess(catRight)

        ## run trimmomatic (assuming FR)
        out1 = os.path.join(analysisDir,"%s_left_paired.fq"%sample)
        out2 = os.path.join(analysisDir,"%s_left_unpaired.fq"%sample)
        out3 = os.path.join(analysisDir,"%s_right_paired.fq"%sample)
        out4 = os.path.join(analysisDir,"%s_right_unpaired.fq"%sample)
        
        cmd = "java -jar /usr/share/java/trimmomatic-0.32.jar PE -threads 29 -phred33 " +\
              "%s %s "%(allLeftFile,allRightFile) +\
              "%s %s "%(out1,out2)+\
              "%s %s "%(out3,out4)+\
              "ILLUMINACLIP:/usr/src/trinityrnaseq-2.0.4/trinity-plugins/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:10 "+\
              "LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:36"

        if not os.path.exists(out1):
            run_subprocess(cmd)

        ## remove 
        for fileName in faFileList:
            if os.path.exists(fileName):
                print("removing %s"%fileName)
                os.remove(fileName)

    print "complete"
 
if __name__ == "__main__":

    sampleList = ["17", "18", "33", "46", "56", "61", "DL47", "DL61", "D163", "D178", "D185", "D239"]
    assemble_reads(sampleList)
