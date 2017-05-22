#!/usr/bin/python
"""
parse the results from xenopus BLAST search
"""

import os,sys,getopt
from htsint.blast import ParseParallelBlast


## read input
homeDir = os.path.join(os.path.expanduser("~"),"sequencing","pieris")
queryFilePath = os.path.realpath(os.path.join(homeDir,"dn-trinity","Trinity.fasta"))

if os.path.exists(queryFilePath) == False:
    raise Exception("Cannot find sequences")

## the chunks must be the same as ParallelBlast
chunks = 29
parser = ParseParallelBlast(queryFilePath)
parser.run(chunks)
