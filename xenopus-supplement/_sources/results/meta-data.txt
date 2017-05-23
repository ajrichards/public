.. notes from the literature


Meta-data analysis
==========================

Target variable correlation
______________________________


There are four target variables 

   * ``vmax`` --- maximum swim speed
   * ``amax`` --- maximum swim acceleration
   * ``tmax-stam1`` --- stamina  (timed jumped till exhaustion)
   * ``dmax-stam1`` --- endurance (distance jumped till exhaustion)

We investigate the correlation among target variables, because classification results will be biased towards the correlated dimensions.

.. INCLUDE target_variable_correlation.png

.. image:: target_variable_correlation.png
   :width: 50%

There is a correlation between stamina and endurance.


Gaussian mixed model on target variables
__________________________________________

Because classes in the case of performance experiments are arbitrary, we explore the target variables with a Gaussian mixture model (GMM).  We are essentially, exploring several clustering configurations in order to understand the data in the 4 dimensions.

.. INCLUDE gmm_on_targets.png

.. image:: gmm_on_targets.png
   :width: 60%


Endurance Classes
____________________

We use two clusters to represent the data.  If we sample from the multi-variate Gaussian distributions with :math:`n=1000` we can produce the following plots that give us some insight into the structure of these data.

.. INCLUDE gmm_classes.png

.. image:: gmm_classes.png
   :width: 60%


Given the GMM shown in these dimensions, we can make some basic observations about the three classes.

+----------------------+-----------+-----------+------------+------------+
| Class                | amax      | vmax      | dmax-stam1 | tmax-stam1 |
+======================+===========+===========+============+============+
| blue (non-endurant)  | dispersed | high      | low        | low        |
+----------------------+-----------+-----------+------------+------------+
| black (endurant)     | dispersed | low       | high       | high       |
+----------------------+-----------+-----------+------------+------------+

We can distinguish the classes as **non-endurant** (blue) and  **endurant** (orange).  The probabilities of each classification are shown in parenthesis.

Classifications
^^^^^^^^^^^^^^^^^^

   * ``A`` = **['endurant']** (1.0)
   * ``B`` = **['non-endurant']** (1.0)
   * ``C`` = **['non-endurant']** (1.0)
   * ``D`` = **['endurant']** (0.9985)
   * ``E`` = **['endurant']** (1.0)
   * ``F`` = **['non-endurant']** (1.0)
   * ``G`` = **['endurant']** (1.0)
   * ``H`` = **['non-endurant']** (0.9995)

Feature selection
___________________

Because a large number of co-variates were measured throughout the course of the experiment we use feature selection to reduce the dimensionality of the problem and trim variables that are not descriptive of performance classificaiton.  Here we use a univariate feature selector and we show the weights from linear support vector machine (SVM) models [Chang08]_.

.. INCLUDE feature_selection.png

.. image:: feature_selection.png
   :width: 60%

Here are the features ranked based on a sum of the two scores.

   1. ``head-length`` 1.0,0.62,1.62
   2. ``radius`` 0.66,0.82,1.48
   3. ``head-width`` 0.05,1.0,1.05
   4. ``tibia`` 0.42,0.19,0.61
   5. ``finger`` 0.17,0.37,0.54
   6. ``hand`` 0.18,0.3,0.48
   7. ``ilium-length`` 0.22,0.22,0.44
   8. ``femur`` 0.3,0.09,0.39
   9. ``mass`` 0.35,0.02,0.38
   10. ``svl`` 0.1,0.24,0.34
   11. ``ilium-width`` 0.3,0.0,0.3
   12. ``head-height`` 0.08,0.22,0.3
   13. ``hindlimb-length`` 0.21,0.05,0.26
   14. ``foot`` 0.2,0.0,0.2
   15. ``humerus`` 0.01,0.12,0.12
   16. ``jaw-length`` 0.08,0.0,0.08
   17. ``fll`` 0.04,0.01,0.06
   18. ``toe`` 0.01,0.0,0.02
