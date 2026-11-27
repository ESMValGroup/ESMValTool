.. _recipes_ocean_south:

Simple southern ocean diagnostics
=====

Overview
--------

A few simple recipes recreated in ESMValTool for running on bulk datasets e.g. CMIP6, from the `COSIMA cookbook <https://cosima-recipes.readthedocs.io/en/latest/index.html>`_.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

* recipe_ocean_south.yml

Diagnostics are stored in esmvaltool/diag_scripts/ocean_south/

* plot_barotropic.py: for plotting barotropic streamfunction
* plot_hovmoller.py: for plotting Hovmoller diagrams of ocean temperature and salinity
* plot_kinetic_energy.py: for plotting ocean kinetic energy


User settings in recipe
-----------------------

*None*


Variables
---------

* umo (ocean, monthly mean, longitude latitude time)
* so (ocean, monthly mean, longitude latitude time depth)
* thetao (ocean, monthly mean, longitude latitude time depth)
* uo (ocean, monthly mean, longitude latitude time depth)
* vo (ocean, monthly mean, longitude latitude time depth)



References
----------

* https://cosima-recipes.readthedocs.io/en/latest/02-Appetisers/Barotropic_Streamfunction.html
* https://cosima-recipes.readthedocs.io/en/latest/02-Appetisers/Hovmoller_Temperature_Depth.html
* https://cosima-recipes.readthedocs.io/en/latest/03-Mains/Eddy-Mean_Kinetic_Energy_Decomposition.html


Example plots
-------------

.. _fig_1:
.. figure::  /recipes/figures/ocean_south/hovmoller_ACCESS-OM2.png
   :align:   center

   Hovmoller diagram of ocean temperature and salinity for the ACCESS-OM2 model.

.. _fig_2:
.. figure::  /recipes/figures/ocean_south/kinetic_energy_ACCESS-OM2-025.png
   :align:   center

   Kinetic energy around Antarctica for the ACCESS-OM2-025 model.

.. _fig_3:
.. figure::  /recipes/figures/ocean_south/barotropic_streamfunction_ACCESS-OM2.png
   :align:   center

   Barotropic streamfunction around Antarctica for the ACCESS-OM2 model.
