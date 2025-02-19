.. _recipes_tcre:

Transient Climate Response to Cumulative CO\ :sub:`2` Emissions (TCRE)
======================================================================

Overview
--------

The idea that global temperature rise is directly proportional to the total
amount of carbon dioxide (CO\ :sub:`2`) released into the atmosphere is
fundamental to climate policy.
The concept stems from research showing a clear linear relationship between
cumulative CO\ :sub:`2` emissions and global temperature change in climate
models (Allen et al. 2009; Matthews et al. 2009; Zickfeld et al. 2009).
This relationship is called the Transient Climate Response to Cumulative CO\
:sub:`2` Emissions (TCRE), which represents the amount of global warming caused
by each trillion tonnes of carbon emitted.
This simple yet powerful tool allows policymakers to directly link emission
budgets to specific temperature targets and compare the long-term effects of
different emissions scenarios.

For CMIP7, TCRE should be calculated from the ``esm-flat10`` experiment,
following the experimental protocol of Sanderson et al (2024).

Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

* recipe_tcre.yml
* ref/recipe_tcre.yml

Diagnostics are stored in diag_scripts/

* :ref:`climate_metrics/tcre.py <api.esmvaltool.diag_scripts.climate_metrics.tcre>`
* :ref:`climate_metrics/barplot.py <create_barplot.py>`


Variables
---------

* *tas* (atmos, monthly, longitude, latitude, time)
* *fco2antt* (atmos, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*


References
----------

* Sanderson, B. M., Brovkin, V., Fisher, R., Hohn, D., Ilyina, T., Jones, C.,
  Koenigk, T., Koven, C., Li, H., Lawrence, D., Lawrence, P., Liddicoat, S.,
  Macdougall, A., Mengis, N., Nicholls, Z., O'Rourke, E., Romanou, A.,
  Sandstad, M., Schwinger, J., Seferian, R., Sentman, L., Simpson, I., Smith,
  C., Steinert, N., Swann, A., Tjiputra, J., and Ziehn, T.: flat10MIP: An
  emissions-driven experiment to diagnose the climate response to positive,
  zero, and negative CO\ :sub:`2` emissions, EGUsphere [preprint],
  https://doi.org/10.5194/egusphere-2024-3356, 2024.
* Allen, M., Frame, D., Huntingford, C. et al.: Warming caused by cumulative
  carbon emissions towards the trillionth tonne, Nature 458, 1163–1166,
  https://doi.org/10.1038/nature08019, 2009.
* Matthews, H., Gillett, N., Stott, P. et al.: The proportionality of global
  warming to cumulative carbon emissions, Nature 459, 829–832,
  https://doi.org/10.1038/nature08047, 2009.
* K. Zickfeld, M. Eby, H. D. Matthews, & A. J. Weaver: Setting cumulative
  emissions targets to reduce the risk of dangerous climate change, Proc. Natl.
  Acad. Sci. U.S.A., 106 (38), 16129-16134,
  https://doi.org/10.1073/pnas.0805800106, 2009.


Example plots
-------------

.. _fig_tcre_1:
.. figure:: /recipes/figures/tcre/tcre.jpg
   :align: center
   :width: 50%

   Global mean surface air temperature anomaly versus cumulative CO\ :sub:`2`
   emissions for MPI-ESM1-2-LR using the emission-driven 1% CO\ :sub:`2`
   increase per year experiment.
