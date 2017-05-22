#!/usr/bin/python
"""
paired-end data has two files ( R1.fastq and R2.fastq).
FR is the conventional order

http://trinityrnaseq.sourceforge.net/
"""

import os,sys,re,shutil
from htsint import run_subprocess

## locations
seqDir = os.path.join(os.path.expanduser("~"),"sequencing","pieris")
readsDir = os.path.join(seqDir,"reads")
assemblyDir = os.path.join(seqDir,"dn-trinity")
trinityHome = "/usr/src/trinityrnaseq-2.0.4"

def print_trinity_cmd(sampleList,clean=False):
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

    allLeft = []
    allRight = []
    for sample in sampleList:
        left_gz1,right_gz1 = get_gz_files(sample,lane1Dir)
        left_gz2,right_gz2 = get_gz_files(sample,lane2Dir)
        allLeft.extend([left_gz1,left_gz2])
        allRight.extend([right_gz1,right_gz2])


    ## construct the command
    cmd = "%s/Trinity --seqType fq --output %s --trimmomatic --full_cleanup "%(trinityHome,assemblyDir)+\
          "--SS_lib_type FR --max_memory 26G --CPU 29 --normalize_reads "+\
          "--left {},{},{},{},{},{},{},{} ".format(*allLeft)+\
          "--right {},{},{},{},{},{},{},{} 2>&1 | tee ./run-trinity-dn.log".format(*allRight)

    #for i,f in enumerate(allLeft):
    #    print f,allRight[i]

    print cmd
    #run_subprocess(cmd)
 
if __name__ == "__main__":

    #if os.path.isdir("trinity_out_dir"):
    #    shutil.rmtree("trinity_out_dir")
    #    os.mkdir("trinity_out_dir")

    sampleList = ["17", "18", "33", "46", "56", "61", "DL47", "DL61", "D163", "D178", "D185", "D239"]
    print_trinity_cmd(sampleList)
