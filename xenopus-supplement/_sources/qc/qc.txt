.. xenopus project

Quality control
=================

Summary
__________

The following samples were renamed to simplify the results presentation.

   * A_Aprime --> A
   * C_Cprime --> C
   * Dprime --> D
   * Eprime --> E
   * F_Fprime --> F
   * G_Gprime --> G
   * Hprime --> H

The following table summarizes the `quality scores <http://en.wikipedia.org/wiki/Phred_quality_score>`_ for the sequencing.

   +--------+----------------+--------------+---------+---------------------+
   | Sample | Yield (Mbases) | # reads      | Q30     | Mean quality score  |
   +========+================+==============+=========+=====================+
   | A      | 1792           | 17923310     |  94.02  | 36.42               |
   +--------+----------------+--------------+---------+---------------------+
   | B      | 3334           | 33344650     |  93.97  | 36.41               |
   +--------+----------------+--------------+---------+---------------------+
   | C      | 5702           | 57018924     |  93.71  | 36.31               |
   +--------+----------------+--------------+---------+---------------------+
   | D      | 2829           | 28291958     |  94.07  | 36.44               |
   +--------+----------------+--------------+---------+---------------------+
   | E      | 4133           | 41329780     |  94.19  | 36.49               |
   +--------+----------------+--------------+---------+---------------------+
   | F      | 3139           | 31394994     |  94.06  | 36.44               |
   +--------+----------------+--------------+---------+---------------------+
   | G      | 2329           | 23294654     |  94.37  | 36.54               |
   +--------+----------------+--------------+---------+---------------------+
   | H      | 1058           | 10576576     |  94.00  | 36.41               |
   +--------+----------------+--------------+---------+---------------------+


Q30 indicates the percentage of reads with >= Q30 (PF).


   * all lanes were run with spikein control data
   * the %PF for all lanes was 100
   * the percent of perfect read indices was 100 for all lanes

Phred quality scores can be intrepreted as follows.  A score of 10
means that there is a 1 in 10 chance of a bad call or a 0.90% chance
of accuracy. A score of 20 means a 1 in 100 chance and a score of 30
means a score of 1 in 1000 chance of being correct.  So quality score
greater than 30 means that there is > 99.9% chance that the call is
accurate.  Ideally a mean quality score should be above 30.


Duplicate check
___________________

The following section indicates the proportion of reads in each sample
which are duplicates. These could be biologically meaningful, but they
can indicate the presence of PCR duplicates which are not
representative of your biological sample. Only read 1 is checked.


   +--------+------------+----------------------+------------+-------------+--------------+-------------+----------+
   | sample | duplicates | collapsed duplicates | singletons | total reads | % singletons | % duplicate | % unique |
   +========+============+======================+============+=============+==============+=============+==========+
   | A      | 816924     | 6206014              | 2755641    | 8961655     | 0.308        | 0.693       | 0.132    |
   +--------+------------+----------------------+------------+-------------+--------------+-------------+----------+
   | B      | 1320499    | 11763702             | 4908623    | 16672325    | 0.294        | 0.706       | 0.112    |
   +--------+------------+----------------------+------------+-------------+--------------+-------------+----------+
   | C      | 2132432    | 21863379             | 6646083    | 28509462    | 0.233        | 0.767       | 0.097    |
   +--------+------------+----------------------+------------+-------------+--------------+-------------+----------+
   | D      | 1049441    | 9535735              | 4610244    | 14145979    | 0.326        | 0.674       | 0.110    |
   +--------+------------+----------------------+------------+-------------+--------------+-------------+----------+
   | E      | 1941564    | 15375232             | 5289658    | 20664890    | 0.256        | 0.744       | 0.126    |
   +--------+------------+----------------------+------------+-------------+--------------+-------------+----------+
   | F      | 1338812    | 10160529             | 5536968    | 15697497    | 0.353        | 0.647       | 0.132    |
   +--------+------------+----------------------+------------+-------------+--------------+-------------+----------+
   | G      | 900514     | 7177176              | 4470151    | 11647327    | 0.384        | 0.616       | 0.126    |
   +--------+------------+----------------------+------------+-------------+--------------+-------------+----------+
   | H      | 437911     | 3333163              | 1955125    | 5288288     | 0.370        | 0.630       | 0.131    |
   +--------+------------+----------------------+------------+-------------+--------------+-------------+----------+


* `duplicates` - the number of unique reads with duplicates
* `% unique`   - the unique reads with duplicates


Plots and additional information
____________________________________

.. toctree::
   :maxdepth: 1

   quality-scores
   base-distributions
   contaminant-check


Center specific details
___________________________

   * Generated on Illumina HiSeq 2500
   * CASAVA 1.8.2
   * RTA 1.13
   * Run ID: 130530_SN982_0225_AC20RLACXX
