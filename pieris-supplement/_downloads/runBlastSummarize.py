#!/usr/bin/env python
"""
Take the parsed blast results and create a summary file to be read by BlastMapper

"""

import os,sys,csv,re,getopt,time

from htsint.blast import BlastMapper

homeDir = os.path.join(os.path.expanduser("~"),"sequencing","pieris")
parsedFilePath = os.path.realpath(os.path.join(homeDir,"dn-trinity","blast-dm-parsed.csv"))
bm = BlastMapper()
bm.create_summarized(parsedFilePath)
