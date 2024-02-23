.. _recipes_psyplot_diag:

Psyplot Diagnostics
===================

Overview
--------

These recipes showcase the use of the Psyplot diagnostic that provides a
high-level interface to `Psyplot <https://psyplot.github.io/>`__ for ESMValTool
recipes.


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_psyplot.yml

Diagnostics are stored in diag_scripts/

   * :ref:`psyplot_diag.py <api.esmvaltool.diag_scripts.psyplot_diag>`


Variables
---------

Arbitrary variables are supported.


Observations and reformat scripts
---------------------------------

Arbitrary datasets are supported.


References
----------

* Sommer, (2017), The psyplot interactive visualization framework, Journal of
  Open Source Software, 2(16), 363, doi:10.21105/joss.00363


Example plots
-------------

.. _fig_psyplot_1:
.. figure:: /recipes/figures/psyplot/psyplot_CanESM5.jpg
   :align: center
   :width: 50%

   Historical near-surface air temperature climatology over Europe simulated by
   CanESM5 between 1995 and 2014. The plot visualizes the invidividual
   rectangular grid cells of the model's regular grid.

.. _fig_psyplot_2:
.. figure:: /recipes/figures/psyplot/psyplot_ICON-ESM-LR.jpg
   :align: center
   :width: 50%

   Historical near-surface air temperature climatology over Europe simulated by
   ICON-ESM-LR between 1995 and 2014. The plot visualizes the invidividual
   triangular grid cells of the model's unstructured grid.
