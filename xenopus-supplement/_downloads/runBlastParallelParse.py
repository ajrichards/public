#!/usr/bin/python
"""
parse the results from xenopus BLAST search
"""

import os,sys,getopt
from htsint.blast import ParseParallelBlast


## read input
if len(sys.argv) < 2:
    raise Exception("\n\nInvalid syntax. use '%s'"%sys.argv[0] + " -d database\n")

try:
    optlist, args = getopt.getopt(sys.argv[1:], 'a:')
except getopt.GetoptError:
    raise Exception(sys.argv[0] + " -a assembly")
    sys.exit()

assembly = 'None'
for o,a in optlist:
    if o == '-a':
        assembly = a

homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus")
#targetDb = "uniprot_sprot.fasta.db"

if assembly == 'gg':
    queryFilePath = os.path.realpath(os.path.join(homeDir,"gg-trinity","Trinity-GG.fasta"))
    targetDb = "xtropicalis"
    blastCmd = 'blastx'
elif assembly == 'dn':
    queryFilePath = os.path.realpath(os.path.join(homeDir,"dn-trinity","Trinity.fasta"))
    targetDb = "xtropicalis"
    blastCmd = 'blastx'
elif assembly == "ref":
    queryFilePath = os.path.realpath(os.path.join(".","Xentr7_2_Stable_Transcript.fa"))
    targetDb = "xtropicalis" 
    blastCmd = 'blastx'
else:
    raise Exception("Invalid assembly specified")

if os.path.exists(queryFilePath) == False:
    raise Exception("Cannot find sequences")

## the chunks must be the same as ParallelBlast
chunks = 29
parser = ParseParallelBlast(queryFilePath)
parser.run(chunks)
