#!/usr/bin/python
"""
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Xenopus_Silurana_tropicalis/protein/protein.fa.gz
gunzip -c protein.fa.gz > Xenopus_Silurana_tropicalis_protein.fa
makeblastdb -in Xenopus_Silurana_tropicalis_protein.fa -dbtype prot

blastx - search protein db with translated nt query
blastp - search protein db with protein query

"""

import os,sys,getopt
from htsint.blast import ParallelBlast

## read input
if len(sys.argv) < 2:
    raise Exception(sys.argv[0] + " -a assembly [-c] \nwhere assembly is 'gg','dn' or 'ref' and '-c' specifies a cluster environment\n")
try:
    optlist, args = getopt.getopt(sys.argv[1:], 'a:c')
except getopt.GetoptError:
    raise Exception(sys.argv[0] + " -a assebbly [-c] \nwhere assembly is 'gg','dn' or 'ref' and '-c' specifies a cluster environment\n")
    sys.exit()

cluster = False
assembly= 'None'
for o,a in optlist:
    if o == '-a':
        assembly = a
    if o == '-c':
        cluster = True

minEvalue = 1e-04
outFmt = '5'
cores = 8
BLASTDB = "/usr/local/share/htsint"
homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus")

if assembly == 'gg':
    queryFilePath = os.path.realpath(os.path.join(homeDir,"gg-trinity","Trinity-GG.fasta"))
    targetDb = "uniprot_sprot.fasta.db"
    blastCmd = 'blastx'
elif assembly == 'dn':
    queryFilePath = os.path.realpath(os.path.join(homeDir,"dn-trinity","Trinity.fasta"))
    targetDb = "uniprot_sprot.fasta.db"
    blastCmd = 'blastx'
elif assembly == "ref":
    queryFilePath = os.path.realpath(os.path.join(".","Xentr7_2_Stable_Transcript.fa"))
    targetDb = "uniprot_sprot.fasta.db" 
    blastCmd = 'blastx'
else:
    raise Exception("Invalid assembly specified")

## error checking
if not os.path.exists(queryFilePath):
    raise Exception("Cannot find file: %s"%queryFilePath)
targetFilePath = os.path.join(BLASTDB,targetDb+".phr")
if not os.path.exists(targetFilePath):
    raise Exception("Cannot fine file: %s"%targetFilePath)

if not cluster:
    print("export BLASTDB='/usr/local/share/htsint'\n"+\
          "%s -query %s "%(blastCmd, queryFilePath)+\
          "-db %s -out %s.outfmt%s "%(targetDb,"blast-"+assembly,outFmt)+\
          "-evalue %s -num_threads %s "%(minEvalue,cores)+\
          "-max_target_seqs 1 -outfmt %s"%(outFmt))
else:
    parallelBlast = ParallelBlast(queryFilePath,targetDb,cmd=blastCmd,BLASTDB=BLASTDB)
    parallelBlast.evalue = minEvalue

    chunks = 29
    parallelBlast.create_scripts(chunks,"adam.richards@ecoex-moulis.cnrs.fr")
    parallelBlast.submit()
