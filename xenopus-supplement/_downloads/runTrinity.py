#!/usr/bin/python
"""
tutorials and docs

paired-end data has two files ( R1.fastq and R2.fastq).
FR is the conventional order
http://www.cureffi.org/2012/12/19/forward-and-reverse-reads-in-paired-end-sequencing/
http://trinityrnaseq.sourceforge.net/
http://www-huber.embl.de/users/anders/HTSeq/doc/tour.html
Trinity --genome genome.fasta --genome_guided_max_intron 10000 --genome_guided_sort_buffer 10G \
        --seqType fq --JM 4G --left reads_1.fq  --right reads_2.fq --CPU 10

## for a 32 core machine (32GB RAM)
--JM 4G  --genome_guided_sort_buffer 4G --CPU 8


## example from docs
../../Trinity --left top100k.Left.fq --right top100k.Right.fq --jaccard_clip --genome top100k.genome 
              --genome_guided_max_intron 1000 --JM 1G --seqType fq --SS_lib_type RF
              --output test_GG_Trinity_outdir

"""

import os,sys,re,shutil,getopt
from htsint import run_subprocess

## read input
if len(sys.argv) < 2:
    raise Exception(sys.argv[0] + " -m method [-m] is 'dn' or 'gg'")

try:
    optlist, args = getopt.getopt(sys.argv[1:], 'm:')
except getopt.GetoptError:
    raise Exception(sys.argv[0] + " -m method [-m] is 'dn' or 'gg'")

method = None
for o,a in optlist:
    if o == '-m':
        method = a

## locations
homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus")
readsDir = os.path.join(homeDir,'reads')
genomePath = os.path.join(homeDir,"star_all_reads_sorted.bam")
trinityHome = "/usr/src/trinityrnaseq-2.0.4"

if method == 'gg':
    assemblyDir = os.path.join(homeDir,"gg-trinity")
elif method == 'dn':
    assemblyDir = os.path.join(homeDir,"dn-trinity")
else:
    raise Exception("Invalid argument for assembly")

if method == 'gg' and not os.path.exists(genomePath):
    raise Exception("Invalid genomePath specified")

if not os.path.isdir(assemblyDir):
    os.mkdir(assemblyDir)

def assemble_reads(sampleList):

    ## get the files
    def get_fq_files(sample):
        sampleDir = os.path.join(homeDir,"Sample_%s"%sample,"raw_illumina_reads")
        leftFile,rightFile = None,None
        for fileName in os.listdir(sampleDir):
            if re.search("^%s.*\.fastq$"%sample,fileName) and re.search("R1",fileName):
                leftFile = os.path.realpath(os.path.join(sampleDir,fileName))
            elif re.search("^%s.*.fastq$"%sample,fileName) and re.search("R2",fileName):
                rightFile = os.path.realpath(os.path.join(sampleDir,fileName))
        return leftFile,rightFile

    allLeft = []
    allRight = []
    for sample in sampleList:
        left,right = get_fq_files(sample)
        allLeft.append(left)
        allRight.append(right)

    ## check
    if len(allLeft) != len(sampleList):
        raise Exception("Invalid number of left sequences")
    if len(allRight) != len(sampleList):
        raise Exception("Invalid number of right sequences")

    ## construct the command
    if method == 'gg':
        cmd = "%s/Trinity --seqType fq --output %s --trimmomatic --full_cleanup "%(trinityHome, assemblyDir)+\
              "--SS_lib_type FR --max_memory 26G --CPU 29 --normalize_reads "+\
              "--genome_guided_bam %s --genome_guided_max_intron 10000 "%(genomePath)+\
              "--left {},{},{},{},{},{},{},{} ".format(*allLeft)+\
              "--right {},{},{},{},{},{},{},{} 2>&1 | tee ./run-trinity-gg.log".format(*allRight)
    else:
        cmd = "%s/Trinity --seqType fq --output %s --trimmomatic --full_cleanup "%(trinityHome,assemblyDir)+\
              "--SS_lib_type FR --max_memory 26G --CPU 29 --normalize_reads "+\
              "--left {},{},{},{},{},{},{},{} ".format(*allLeft)+\
              "--right {},{},{},{},{},{},{},{} 2>&1 | tee ./run-trinity-dn.log".format(*allRight)

    print cmd
    #run_subprocess(cmd)
 
if __name__ == "__main__":

    sampleList = ["A", "B", "C", "D", "E", "F", "G", "H"]
    assemble_reads(sampleList)
