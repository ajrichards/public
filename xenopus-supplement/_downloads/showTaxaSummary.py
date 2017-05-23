#!/usr/bin/env python
"""
Take the parsed blast results and create a summary file to be read by BlastMapper

"""

import os,sys,csv,re,getopt,time

from htsint.blast import BlastMapper


## read input 
if len(sys.argv) < 2:
    raise Exception("\n\nInvalid syntax. use '%s'"%sys.argv[0] + " -s source\n")
try:
    optlist, args = getopt.getopt(sys.argv[1:], 's:')
except getopt.GetoptError:
        raise Exception("\n\nInvalid syntax. use '%s'"%sys.argv[0] + " -s source\n")

source = 'None'
for o,a in optlist:
    if o == '-s':
        source = a

homeDir = os.path.join(os.path.expanduser("~"),"sequencing","xenopus")
if source in ['dn','gg']:
    sourceDir = "%s-trinity"%(source)
elif source == 'ref':
    sourceDir = "reference"
else:
    raise Exception("Bad source")

summaryFile1 = os.path.join(homeDir,sourceDir,"blast-%s-parsed_summary.csv"%(source))
summaryFile2 = os.path.join(homeDir,sourceDir,"blast-xt-parsed_summary.csv")
bm = BlastMapper()

## load the gene and isoform maps
bmapIsoforms = bm.load_summary(summaryFile1,trinityGene=False,best=True)
bmapFrog = bm.load_summary(summaryFile1,trinityGene=False,best=True,taxaList=['8355','8364'])
bmapXT = bm.load_summary(summaryFile2,trinityGene=False,best=True)

print("-----------")
print("SwissProt - isoforms")
bm.print_summary(bmapIsoforms)
print("SwissProt [8355,8364] - isoforms")
bm.print_summary(bmapFrog)
print("X. tropicalis - isoforms")
bm.print_summary(bmapXT)


bm.make_taxa_pie_chart_and_table(bmapIsoforms,removeStrain=True,
                                 figName="%s-trinity-blast-pie-isoforms.png"%(source),
                                 csvName="%s-trinity-blast-species-isoforms.csv"%(source))

bm.make_taxa_pie_chart_and_table(bmapFrog,removeStrain=True,
                                 figName="%s-trinity-blast-pie-frog.png"%(source),
                                 csvName="%s-trinity-blast-species-frog.csv"%(source))
