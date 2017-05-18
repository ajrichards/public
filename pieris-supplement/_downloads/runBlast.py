#!/usr/bin/python
"""
blastx - search protein db with translated nt query
blastp - search protein db with protein query

"""

import os,sys,getopt
from htsint.blast import ParallelBlast

## read input
try:
    optlist, args = getopt.getopt(sys.argv[1:], 'c')
except getopt.GetoptError:
    raise Exception(sys.argv[0] + " '-c' optionally specifies a cluster environment\n")
    sys.exit()

cluster = False
for o,a in optlist:
    if o == '-c':
        cluster = True

assembly = 'dn'
minEvalue = 1e-04
outFmt = '5'
cores = 20
BLASTDB = "/usr/local/share/htsint"
homeDir = os.path.join(os.path.expanduser("~"),"sequencing","pieris")
queryFilePath = os.path.realpath(os.path.join(homeDir,"dn-trinity","Trinity.fasta"))
blastCmd = 'blastx'
targetDb = "uniprot_sprot.fasta.db"

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
