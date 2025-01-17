.. _recipes_seaborn_diag:

Seaborn Diagnostics
===================

Overview
--------

These recipes showcase the use of the Seaborn diagnostic that provides a
high-level interface to `Seaborn <https://seaborn.pydata.org>`__ for ESMValTool
recipes used in the CMIP Rapid Evaluation Framework (REF).

CMIP Rapid Evaluation Framework
---------------------------------

`CMIP Rapid Evaluation Framework <https://wcrp-cmip.org/cmip7/rapid-evaluation-framework/>`

Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/ref/

* recipe_scatterplot.yml

Diagnostics are stored in diag_scripts/

* :ref:`seaborn_diag.py <api.esmvaltool.diag_scripts.seaborn_diag>`


Variables
---------

Clt, Swcre, Lwcre.


Observations and reformat scripts
---------------------------------

Arbitrary datasets are supported.


References
----------

* Waskom, M. L. (2021), seaborn: statistical data visualization, Journal of
  Open Source Software, 6(60), 3021, doi:10.21105/joss.03021.


Example plots
-------------

.. _fig_seaborn_1:
.. figure:: /recipes/figures/seaborn/clt_vs_swcre.png
   :align: center
   :width: 50%

.. _fig_seaborn_2:
.. figure:: /recipes/figures/seaborn/clt_vs_lwcre.png
   :align: center
   :width: 50%

.. _fig_seaborn_1:
.. figure:: /recipes/figures/seaborn/hist_clt_vs_swcre.png
   :align: center
   :width: 50%

.. _fig_seaborn_2:
.. figure:: /recipes/figures/seaborn/hist_clt_vs_lwcre.png
   :align: center
   :width: 50%

.. _fig_seaborn_1:
.. figure:: /recipes/figures/seaborn/hex_clt_vs_swcre.png
   :align: center
   :width: 50%

.. _fig_seaborn_2:
.. figure:: /recipes/figures/seaborn/hex_clt_vs_lwcre.png
   :align: center
   :width: 50%      

   