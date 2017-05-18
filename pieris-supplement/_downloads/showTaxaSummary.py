#!/usr/bin/env python
"""
Take the parsed blast results and create a summary file to be read by BlastMapper

"""

import os,sys,csv,re,getopt,time

from htsint.blast import BlastMapper

homeDir = os.path.join(os.path.expanduser("~"),"sequencing","pieris")
parsedFilePath = os.path.realpath(os.path.join(homeDir,"dn-trinity","blast-dn-parsed.csv"))
summaryFile1 = os.path.join(homeDir,"dn-trinity","blast-dn-parsed_summary.csv")
summaryFile2 = os.path.join(homeDir,"dn-trinity",'blast-dm-parsed_summary.csv')
bm = BlastMapper()

## load the gene and isoform maps
bmapIsoforms = bm.load_summary(summaryFile1,trinityGene=False,best=True)
bmapInsects = bm.load_summary(summaryFile1,trinityGene=False,best=True,taxaList=['7227','7091'])
bmapDM = bm.load_summary(summaryFile2,trinityGene=False,best=True)

print("-----------")
print("SwissProt - isoforms")
bm.print_summary(bmapIsoforms)
print("SwissProt [7227,7091] - isoforms")
bm.print_summary(bmapInsects)
print("D. melanogaster - isoforms")
bm.print_summary(bmapDM)

bm.make_taxa_pie_chart_and_table(bmapIsoforms,removeStrain=True,
                                 figName="dn-trinity-blast-pie-isoforms.png",
                                 csvName="dn-trinity-blast-species-isoforms.csv")

bm.make_taxa_pie_chart_and_table(bmapInsects,removeStrain=True,
                                 figName="dn-trinity-blast-pie-insects.png",
                                 csvName="dn-trinity-blast-species-insects.csv")
