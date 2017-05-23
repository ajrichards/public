.. xenopus supplement


.. _gene-set-methods:


Gene sets
==========

We partition the composite genomes of multiple species into functional gene sets.  This is because **Xenopus** is not among the most completely annotated taxa.  In this analysis we consider several ways to generate gene sets.  More specifically, the annotation information used may come from a list of different species or we may include/exclude specific evidence codes.  The result is always a list of functional gene sets. The majority of annotations in the GO are automatically inferred and the are non-curated.  They use specifically the ``IEA`` evidence code.  When ``IEA`` is excluded we use only curated annotations.

We are most interested in the following two species of frog.  However, there are other well-annotated species that can help create functional gene sets.

   * *Xenopus (Silurana) tropicalis* `[8364] <http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=8364>`_
   * *Xenopus laevis* `[8355] <http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=8355>`_

+-----------+-------------------------------------------------------------+-----------------------------------------+
| Taxon ID  | Scientific Name                                             | Common Name                             |
+===========+=============================================================+=========================================+
| 8355      | Xenopus laevis                                              | African clawed frog                     |
+-----------+-------------------------------------------------------------+-----------------------------------------+
| 8364      | Xenopus (Silurana) tropicalis                               | western clawed frog                     |
+-----------+-------------------------------------------------------------+-----------------------------------------+

Xenopus annotation Summary
_____________________________

In general annotations may be associated directly with an NCBI GeneID or they may be associated with a gene product that has an assigned Uniprot accession (UniprotKB-AC) ID.  All UniprotKB-AC accessions are mapped to UniprotKb-entries.  To clarify the following tables.

   * ``taxon`` - matches with the taxon summary table on the top of the page
   * ``genes (annots)`` - The number of unique GeneIDs from NCBI that are associated with a given taxon
   * ``uniprot (annots)``  - The identifiers from uniprot (non-unique) that are associated with a given taxon
   * ``coding-genes`` - GeneIDs that have at least one associated uniprot identifier
   * ``combined`` - uniprot and gene annotations with redundant annotations removed
   * ``coverage`` - ``combined`` divided by ``genes``

The ``combined`` number of annotations may be less than the total observed annotations because multiple uniprot ids may map to a single gene id.

Including IEA evidence
""""""""""""""""""""""""

+-----------+-----------------+-------------------+---------------+-----------+-----------+
| Taxon     | genes (annots)  | uniprot (annots)  | coding-genes  | combined  | coverage  |
+===========+=================+===================+===============+===========+===========+
| 8355      | 13531 (0)       | 3385 (3187)       | 2957          | 2773      | 20.4937   |
+-----------+-----------------+-------------------+---------------+-----------+-----------+
| 8364      | 24817 (0)       | 1701 (1565)       | 1615          | 1485      | 5.9838    |
+-----------+-----------------+-------------------+---------------+-----------+-----------+

Excluding IEA evidence
""""""""""""""""""""""""

+-----------+-----------------+-------------------+---------------+-----------+-----------+
| Taxon     | genes (annots)  | uniprot (annots)  | coding-genes  | combined  | coverage  |
+===========+=================+===================+===============+===========+===========+
| 8355      | 13531 (0)       | 3385 (1254)       | 2957          | 1145      | 8.4621    |
+-----------+-----------------+-------------------+---------------+-----------+-----------+
| 8364      | 24817 (0)       | 1701 (633)        | 1615          | 592       | 2.3855    |
+-----------+-----------------+-------------------+---------------+-----------+-----------+


Annotation of other model organisms
______________________________________

+-------+---------------------------------------------------+------------------------------+
| Taxon | Scientific Name                                   | Common Name                  |
+=======+===================================================+==============================+
| 10090 | Mus musculus                                      | house mouse                  |
+-------+---------------------------------------------------+------------------------------+
| 9606  | Homo sapiens                                      | human                        |
+-------+---------------------------------------------------+------------------------------+
| 7227  | Drosophila melanogaster                           | fruit fly                    |
+-------+---------------------------------------------------+------------------------------+
| 7955  | Danio rerio                                       | leopard danio                |
+-------+---------------------------------------------------+------------------------------+


Including IEA evidence
""""""""""""""""""""""""

+-----------+-----------------+-------------------+---------------+-----------+-----------+
| Taxon     | genes (annots)  | uniprot (annots)  | coding-genes  | combined  | coverage  |
+===========+=================+===================+===============+===========+===========+
| 10090     | 69347 (18938)   | 16671 (15735)     | 15619         | 19006     | 27.4071   |
+-----------+-----------------+-------------------+---------------+-----------+-----------+
| 9606      | 47736 (18114)   | 140289 (100435)   | 18921         | 18244     | 38.2185   |
+-----------+-----------------+-------------------+---------------+-----------+-----------+
| 7227      | 23387 (11563)   | 40464 (26100)     | 13543         | 11870     | 50.7547   |
+-----------+-----------------+-------------------+---------------+-----------+-----------+
| 7955      | 36579 (14990)   | 2926 (2697)       | 2669          | 15069     | 41.1958   |
+-----------+-----------------+-------------------+---------------+-----------+-----------+

Excluding IEA evidence
""""""""""""""""""""""""

+-----------+-----------------+-------------------+---------------+-----------+-----------+
| Taxon     | genes (annots)  | uniprot (annots)  | coding-genes  | combined  | coverage  |
+===========+=================+===================+===============+===========+===========+
| 10090     | 69347 (15763)   | 16671 (11428)     | 15619         | 15842     | 22.8445   |
+-----------+-----------------+-------------------+---------------+-----------+-----------+
| 9606      | 47736 (14077)   | 140289 (17803)    | 18921         | 14263     | 29.8789   |
+-----------+-----------------+-------------------+---------------+-----------+-----------+
| 7227      | 23387 (10338)   | 40464 (12648)     | 13543         | 10500     | 44.8967   |
+-----------+-----------------+-------------------+---------------+-----------+-----------+
| 7955      | 36579 (3962)    | 2926 (1261)       | 2669          | 4006      | 10.9516   |
+-----------+-----------------+-------------------+---------------+-----------+-----------+


Links
____________

   * `Statistics about the current GO annotations <http://www-test.geneontology.org/page/current-go-statistics>`_
   * :doc:`Information about the database used for these data </methods/database>`
   * `Guide to GO evidence codes <http://www.geneontology.org/page/guide-go-evidence-codes>`_
