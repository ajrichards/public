.. xenopus supplement file, created by ARichards

==========
Database
==========

About the database
_____________________

The databased we used for this analysis used `PostgreSQL <http://www.postgresql.org>`_ and several python modules, including the object `relational mapper <http://en.wikipedia.org/wiki/Object-relational_mapping>`_  provided by `SQLAlchemy <http://www.sqlalchemy.org>`_.  The main purpose of the database is to be able to handle at a detailed level Gene Ontology [Ashburner00]_ information. Here is the relational schema produced with a tool called `sqlalchemy_schemadisplay <https://pypi.python.org/pypi/sqlalchemy_schemadisplay>`_. The database population scripts along with a number of classes to interface with the data are made possible through a project called `hts-integrate <https://github.com/ajrichards/hts-integrate>`_.

.. figure:: dbschema.png
   :scale: 70%
   :align: center
   :alt: database schema
   :figclass: align-center

Database contents
___________________

The following files were downloaded on **September 15, 2014**.

   * `go.obo <ftp://ftp.geneontology.org/pub/go/ontology/go.obo>`_
   * `taxdump.tar.gz <ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz>`_
   * `gene_info.gz <ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz>`_
   * `gene2go.gz <ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz>`_
   * `gene2refseq.gz <ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz>`_
   * `gene_association.goa_uniprot.gz <ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/UNIPROT/gene_association.goa_uniprot.gz>`_
   * `idmapping.dat.gz <ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping.dat.gz>`_
   * `uniprot_sprot.fasta.gz <ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz>`_
   * `uniprot_trembl.fasta.gz <ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz>`_

The uniprot and gene centric data from these files was parsed and used to populate the database.  The number of rows in each of the tables are shown below.

   .. code-block:: none
  
      There are 1262260 entries in the taxa table
      There are 681732 entries in the genes table
      There are 777608 entries in the uniprot table
      There are 42627 entries in the go_terms table
      There are 7463568 entries in the go_annotations table

Links
_________

   * `Gene Ontology <http://geneontology.org/>`_
   * `GO annotation file format <http://www.geneontology.org/GO.format.gaf-2_0.shtml>`_
   * `GO annotation README <http://www.geneontology.org/gene-associations/readme/goa.README>`_

   * `Refseq accession number and molecule types <http://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly>`_  
   * `NCBI FTP README <ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/README>`_

   * `Uniprot database <http://www.ebi.ac.uk/uniprot>`_
   * `Uniprot README <ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/README>`_
