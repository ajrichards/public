#!/usr/bin/env python
"""
Take the parsed blast results and create a summary file to be read by BlastMapper

"""

import os,sys,csv,re,getopt,time

from htsint.blast import BlastMapper


## read input 
if len(sys.argv) < 2:
    raise Exception("\n\nInvalid syntax. use '%s'"%sys.argv[0] + " -m method\n")
try:
    optlist, args = getopt.getopt(sys.argv[1:], 'm:')
except getopt.GetoptError:
        raise Exception("\n\nInvalid syntax. use '%s'"%sys.argv[0] + " -m method\n")

method = 'None'
for o,a in optlist:
    if o == '-m':
        method = a

if method not in ['dn','gg','ref']:
    raise Exception("bad method specified")

xtropicalis = False

homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus")

if method == 'ref':
    sourceDir = "reference"
else:
    sourceDir = "%s-trinity"

if xtropicalis:
    parsedFilePath = os.path.realpath(os.path.join(homeDir,sourceDir,"blast-xt-parsed.csv"))
    refseq = True
    print("creating parsed summary file for %s BLASTED against xtropicalis"%(method))
else:
    parsedFilePath = os.path.realpath(os.path.join(homeDir,sourceDir,"blast-%s-parsed.csv"%(method)))
    refseq = False
    print("creating parsed summary file for %s BLASTED against SwissProt"%(method))

bm = BlastMapper()
bm.create_summarized(parsedFilePath,large=True,refseq=refseq)
