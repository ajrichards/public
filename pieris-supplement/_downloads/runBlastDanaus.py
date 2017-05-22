#!/usr/bin/python
"""
blastx - search protein db with translated nt query
blastp - search protein db with protein query
blastn - search nucl db with nucl query
tblastx - compares the six-frame translations of a nucleotide query sequence against the six-frame translations of a nucleotide sequence database.

### (tried this) wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-26/fasta/danaus_plexippus/dna/Danaus_plexippus.DanPle_1.0.26.dna.genome.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-28/fasta/danaus_plexippus/cdna/Danaus_plexippus.DanPle_1.0.28.cdna.all.fa.gz
mv Danaus_plexippus.DanPle_1.0.26.cdna.all.fa.gz /usr/local/share/htsint
gunzip -c Danaus_plexippus.DanPle_1.0.28.cdna.all.fa.gz > Danaus.fa
makeblastdb -in Danaus.fa -dbtype 'nucl' -out Danaus
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
cores = 8
BLASTDB = os.path.join(os.getenv("HOME"),"data","htsint")
os.system("export BLASTDB='%s/:$BLASTDB'"%BLASTDB)
homeDir = os.path.join(os.path.expanduser("~"),"sequencing","pieris")
queryFilePath = os.path.realpath(os.path.join(homeDir,"dn-trinity","Trinity.fasta"))
blastCmd = 'tblastx'
targetDb = "Danaus"

## error checking
if not os.path.exists(queryFilePath):
    raise Exception("Cannot find file: %s"%queryFilePath)
targetFilePath = os.path.join(BLASTDB,targetDb+".nhr")
if not os.path.exists(targetFilePath):
    raise Exception("Cannot fine file: %s"%targetFilePath)

if not cluster:
    print("export BLASTDB='%s'\n"%BLASTDB+\
          "%s -query %s "%(blastCmd, queryFilePath)+\
          "-db %s -out %s.outfmt%s "%(targetDb,"blast-dp",outFmt)+\
          "-evalue %s -num_threads %s "%(minEvalue,cores)+\
          "-max_target_seqs 1 -outfmt %s"%(outFmt))
else:
    parallelBlast = ParallelBlast(queryFilePath,targetDb,cmd=blastCmd,BLASTDB=BLASTDB)
    parallelBlast.evalue = minEvalue

    chunks = 29
    parallelBlast.create_scripts(chunks,"adam.richards@ecoex-moulis.cnrs.fr")
    parallelBlast.submit()
