#!/usr/bin/env python
"""
Take the parsed blast results and create a summary file to be read by BlastMapper

"""

import os,sys,csv,re,getopt,time

from htsint.blast import BlastMapper

homeDir = os.path.join(os.path.expanduser("~"),"sequencing","pieris")
parsedFilePath = os.path.realpath(os.path.join(homeDir,"dn-trinity","blast-dn-parsed.csv"))
summaryFile1 = os.path.join(homeDir,"blast","blast-up-parsed_summary.csv")
summaryFile2 = os.path.join(homeDir,"blast",'blast-dm-parsed_summary.csv')
summaryFile3 = os.path.join(homeDir,"blast",'blast-dp-parsed_summary.csv')

bm = BlastMapper()

## load the gene and isoform maps
bmapSP = bm.load_summary(summaryFile1,trinityGene=False,best=True)
bmapDM = bm.load_summary(summaryFile2,trinityGene=False,best=True)
bmapDP = bm.load_summary(summaryFile3,trinityGene=False,best=True)

print("-----------")
print("SwissProt - isoforms")
bm.print_summary(bmapSP)
print("D. melanogaster - isoforms")
bm.print_summary(bmapDM)
print("Danaus plexippus - isoforms")
bm.print_summary(bmapDP)

bm.make_taxa_pie_chart_and_table(bmapSP,removeStrain=True,
                                 figName="dn-trinity-blast-pie-isoforms.png",
                                 csvName="dn-trinity-blast-species-isoforms.csv")

#bm.make_taxa_pie_chart_and_table(bmapInsects,removeStrain=True,
#                                 figName="dn-trinity-blast-pie-insects.png",
#                                 csvName="dn-trinity-blast-species-insects.csv")
