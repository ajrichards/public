#!/usr/bin/python
"""
Run STAR to align reads to transcriptome
"""

import os,sys,re,shutil,getopt,itertools
import HTSeq
from htsint import run_subprocess

## locations
homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus")
starDir =  os.path.join(homeDir,'star')
readsDir = os.path.join(homeDir,'reads')
genomePath = os.path.realpath(os.path.join(".","20100930","sequences","Xenopus_tropicalis.main_genome.scaffolds.fasta"))
gff3Path  = os.path.realpath(os.path.join(".","Xentr7_2_Stable.gff3"))
gtfPath  = os.path.realpath(os.path.join(starDir,"Xentr7_2_Stable.gtf"))
starPath = "/usr/src/STAR/source/STAR"

def get_reads(sampleList):

    def get_trimmed_files(sample):
        leftFile,rightFile = None,None
        for fileName in os.listdir(readsDir):
            if re.search("unpaired",fileName):
                continue

            if re.search("^%s.*\.fq$"%sample,fileName) and re.search("left",fileName):
                leftFile = os.path.realpath(os.path.join(readsDir,fileName))
            elif re.search("^%s.*.fq$"%sample,fileName) and re.search("right",fileName):
                rightFile = os.path.realpath(os.path.join(readsDir,fileName))
        return leftFile,rightFile

    allLeft = []
    allRight = []
    for sample in sampleList:
        left,right = get_trimmed_files(sample)
        allLeft.append(left)
        allRight.append(right)

    ## check
    if len(allLeft) != len(sampleList):
        raise Exception("Invalid number of left sequences")
    if len(allRight) != len(sampleList):
        raise Exception("Invalid number of right sequences")

    return allLeft,allRight

if __name__ == "__main__":

    sampleList = ["A", "B", "C", "D", "E", "F", "G", "H"]

    ## print version of STAR
    print("\n...get version")
    cmd1 = "%s --version"%(starPath)
    print(cmd1)

    ## convert gff3 to gtf
    cmd2 = "gffread %s -T -o %s"%(gff3Path,gtfPath)
    print("\n...convert gff3 into gtf file")
    print(cmd2)

    ## copy genome to star dir
    cmd3 = "cp %s %s"%(genomePath,os.path.join(starDir,'genome.fa'))
    print("\n...copy genome to STAR dir")
    print(cmd3)
        
    ## create the index
    cmd4 = "%s  --runMode genomeGenerate --runThreadN 24 --genomeDir %s "%(starPath,starDir)+\
           "--genomeFastaFiles %s --sjdbGTFfile %s --sjdbOverhang 100"%(os.path.join(starDir,'genome.fa'),gtfPath)
    print("\n...index the genome")
    print(cmd4)

    ## align reads
    print("\n...align reads")
    allLeft,allRight = get_reads(sampleList)
    cmd5 = ""
    for s,sample in enumerate(sampleList):
        if s != len(sampleList) -1:
            _cmd = "%s --genomeDir %s --runThreadN 24 "%(starPath,starDir)+\
                   "--readFilesIn %s %s --outFileNamePrefix %s_ && "%(allLeft[s],allRight[s],os.path.join(starDir,sample))
        else:
            _cmd = "%s --genomeDir %s --runThreadN 24 "%(starPath,starDir)+\
                   "--readFilesIn %s %s --outFileNamePrefix %s_"%(allLeft[s],allRight[s],os.path.join(starDir,sample))
        cmd5 += _cmd
    print cmd5

    overwrite = False

    ## convert the sam files to coordinate-sorted bam and files
    print("\n...sam/bam conversions and sort")
    for s,sample in enumerate(sampleList):
        samFile = os.path.join(starDir,"%s_Aligned.out.sam"%(sample))
        bamFile = os.path.join(readsDir,"%s_aligned.bam"%(sample))
        sbamFile = os.path.join(readsDir,"%s_aligned_sorted.bam"%(sample))
        ssamFile = os.path.join(readsDir,"%s_aligned_sorted.sam"%(sample))
        if not os.path.exists(samFile):
            raise Exception("cannot find sam file %s"%(samFile))

        if os.path.exists(ssamFile) and not overwrite:
            print("skipping sam to bam, align, bam to sam")
        else:
            cmd = "/usr/bin/samtools view -b -S %s > %s && "%(samFile,bamFile)
            cmd += "/usr/bin/samtools sort -n %s %s && "%(bamFile,sbamFile[:-4])
            cmd += "/usr/bin/samtools view -h %s > %s"%(sbamFile,ssamFile)
            print cmd
            run_subprocess(cmd)

    ## concat sam files and sort by coordinates
    print("\n...make single sorted bam file..")
    outBam = os.path.join(homeDir,"star_all_reads.bam")
    outSbam = os.path.join(homeDir,"star_all_reads_sorted.bam")
    cmdMerge = "/usr/bin/samtools merge %s"%(outBam)
    cmdSort = "/usr/bin/samtools sort %s %s"%(outBam,outSbam[:-4])
    for s,sample in enumerate(sampleList):
        sbamFile = os.path.join(readsDir,"%s_aligned_sorted.bam"%(sample))
        cmdMerge += " %s"%(sbamFile)
        
    cmdMergeSort = cmdMerge + " && %s"%(cmdSort)

    if os.path.exists(outSbam) and not overwrite:
        print("skipping concat bam files and sort")
    else:
        if os.path.exists(outSbam):
            os.remove(outSbam)
        if os.path.exists(outBam):
            os.remove(outBam)
        print cmdMergeSort
        run_subprocess(cmdMergeSort)
    print("\n")
