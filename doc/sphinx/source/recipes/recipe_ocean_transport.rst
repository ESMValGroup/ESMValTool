.. _recipes_ocean_transport:

Ocean transport diagnostics
===========================

Overview
--------

This recipe computes three ocean transport diagnostics for a CMIP6 model:

* **Ocean heat transport:** the northward ocean heat transport in the
  Atlantic-Arctic basin, from the ``hfbasin`` variable. It is shown as a
  profile against latitude (in PW).
* **ACC transport:** the Antarctic Circumpolar Current transport through the
  Drake Passage, from the ``mfo`` variable. It is shown as a time series
  (in Tg s-1).
* **SSH gradient:** the magnitude of the horizontal sea surface height
  gradient, from the ``zos`` variable. It is shown as a map. Strong gradients
  mark ocean fronts such as the Gulf Stream, so this map helps to check the
  Gulf Stream separation.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

* recipe_ocean_transport.yml

Diagnostics are stored in esmvaltool/diag_scripts/ocean/

* diagnostic_ssh_gradient.py: computes and plots the SSH gradient magnitude.

The ocean heat transport and ACC transport diagnostics use standard
preprocessors together with the ``monitor/multi_datasets.py`` plot script.


User settings in recipe
-----------------------

#. Script diagnostic_ssh_gradient.py

   The script has no settings. It expects ``zos`` regridded to a regular
   latitude/longitude grid with the time axis kept, because the script computes
   the gradient magnitude for each time step and then averages over time.

The ocean heat transport and ACC transport diagnostics are controlled through
the preprocessors and the ``monitor/multi_datasets.py`` plot options in the
recipe.


Variables
---------

* hfbasin (ocean, monthly, latitude basin time)
* mfo (ocean, monthly, oline time)
* zos (ocean, monthly, longitude latitude time)


Observations and reformat scripts
---------------------------------

The recipe runs on a CMIP6 model alone. The SSH gradient can also be compared
with the ORAS5 reanalysis. ORAS5 is read on the fly by ESMValCore (project
``native6``) and needs the ``oras5_mesh_T.nc`` grid file, passed through the
``horizontal_grid`` facet. More details in the `ESMValCore documentation on
ORAS5`_. See the commented example in the recipe.

.. _ESMValCore documentation on ORAS5: https://docs.esmvaltool.org/projects/
   ESMValCore/en/latest/quickstart/find_data.html#oras5


References
----------

* de Mora, L., Yool, A., Palmieri, J., Sellar, A., Kuhlbrodt, T., Popova, E.,
  Jones, C., and Allen, J. I.: BGC-val: a model- and grid-independent Python
  toolkit to evaluate marine biogeochemical models, Geosci. Model Dev., 11,
  4215-4240, https://doi.org/10.5194/gmd-11-4215-2018, 2018.


Example plots
-------------

.. _fig_ocean_transport_1:
.. figure::  /recipes/figures/ocean_transport/ocean_heat_transport.png
   :align:   center

   Northward ocean heat transport in the Atlantic-Arctic basin (MPI-ESM1-2-LR).

.. _fig_ocean_transport_2:
.. figure::  /recipes/figures/ocean_transport/acc_drake_transport.png
   :align:   center

   ACC transport through the Drake Passage (MPI-ESM1-2-LR).

.. _fig_ocean_transport_3:
.. figure::  /recipes/figures/ocean_transport/ssh_gradient.png
   :align:   center

   Sea surface height gradient magnitude (MPI-ESM1-2-LR).
