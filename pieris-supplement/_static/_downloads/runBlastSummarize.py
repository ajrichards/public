#!/usr/bin/env python
"""
Take the parsed blast results and create a summary file to be read by BlastMapper

['#Format: tax_id GeneID Ensembl_gene_identifier RNA_nucleotide_accession.version Ensembl_rna_identifier protein_accession.version Ensembl_protein_identifier (tab is used as a separator', ' pound sign - start of a comment)']


"""

import os,sys,csv,re,getopt,time
from htsint.blast import BlastMapper
from sqlalchemy.sql import select
from htsint.database import db_connect,Taxon,Gene,Uniprot,uniprot_mapper

## read input
if len(sys.argv) < 2:
    raise Exception(sys.argv[0] + " -d database [-d] is 'dp','dm' or 'up'")
try:
    optlist, args = getopt.getopt(sys.argv[1:], 'd:')
except getopt.GetoptError:
    raise Exception(sys.argv[0] + " -d database [-d] is 'dp','dm' or 'up'")

db = None
for o,a in optlist:
    if o == '-d':
        db = a

if db not in ['dp','dm','up']:
    raise Exception("-d database [-d] is 'dp','dm' or 'up'")

if db == 'dp':
    species ='Danaus plexippus'
    taxaList=['13037']
    uniprot = False
elif db == 'dm':
    species = 'Drosophila melanogaster'
    taxaList = ['7227']
    uniprot = False
else:
    species = None
    taxaList = []
    uniprot = True
    
homeDir = os.path.join(os.path.expanduser("~"),"sequencing","pieris","blast")
parsedFilePath = os.path.realpath(os.path.join(homeDir,"blast-%s-parsed.csv"%(db)))
bm = BlastMapper()

## read in the gene2ensemble file
if db == 'dm':
    fid = open(os.path.join(homeDir,'gene2ensembl'),'r')
    reader = csv.reader(fid,delimiter="\t")
    header = reader.next()
    id2gene = {}
    for linja in reader:
        if linja[0] != taxaList[0]:
            continue
        id2gene[linja[4]] = linja[1]
elif db == 'dp':
    transcript2uniprot = {}
    fid = open(os.path.join(homeDir,'Danaus_plexippus.DanPle_1.0.28.uniprot.tsv'),'r')
    reader = csv.reader(fid,delimiter="\t")
    header = reader.next()
    for linja in reader:
        transcript2uniprot[linja[1]] = linja[3]

    print transcript2uniprot.keys()[:5]
    print len(transcript2uniprot.keys())
    print transcript2uniprot.values()[:5]
    #sys.exit()

    session,engine = db_connect()
    conn = engine.connect()
    
    #tqueryId = session.query(Taxon).filter(Taxon.ncbi_id=='13037').first().id
    #print tqueryId
    #uniprotQuery = session.query(Uniprot).filter_by(taxa_id=tqueryId).all()
    
    #for uq in uniprotQuery:
    #    print '...', uq
    #_codingGenes = list(set([u.gene_id for u in uniprotQuery if u]))
    #codingGenes = [g.ncbi_id for g in session.query(Gene).filter(Gene.id.in_(_codingGenes)).all() if g]
    #print 'coding genes', len(codingGenes)
    #sys.exit()
    
    upEntry2Gene = {}
    results = conn.execute(Uniprot.__table__.select(Uniprot.uniprot_ac.in_(transcript2uniprot.values())))
    blah = [r for r in results if r]

    hitEntries = transcript2uniprot.values()
    print 'hits', len(hitEntries)
    print 'matched', len(blah)
    
    #print 'results', len(results)
    #for row in results:
    #    upEntry2Gene[str(row.uniprot_entry)] = str(row.gene_id)
    #print upEntry2Gene
    #sys.exit()
    #    ## using htsint's mapper
    #upEntry2Gene,upEntry2Taxa = {},{}
    #uMapper = uniprot_mapper(session,uniprotIdList=hitEntries,gene=True,taxa=True)
    #for key,item in uMapper.iteritems():
    #    upEntry2Gene[str(key)] = str(item['gene_id'])
    #    upEntry2Taxa[str(key)] = str(item['taxa_id'])
    #               
    #print list(set(upEntry2Taxa.values()))
    #print len(upEntry2Gene.keys())
    #
    #
    sys.exit()
else:
    id2gene = None
        
## run
bm.create_summarized(parsedFilePath,species=species,taxaList=taxaList,hit2gene=id2gene,uniprot=uniprot)
