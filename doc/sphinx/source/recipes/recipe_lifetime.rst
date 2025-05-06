.. _recipe_lifetime:

Lifetime
========

Overview
--------

The diagnostic calculates the lifetime of CH4 (or potentially other chemically
active trace gases) and is able to plot it in several different
representations.


Available recipes and diagnostics
---------------------------------

Recipes are stored in `recipes/lifetime`

* recipe_lifetime.yml

Diagnostics are stored in `diag_scripts/lifetime/`

* :ref:`lifetime.py <api.esmvaltool.diag_scripts.lifetime.lifetime>`


Recipe settings
~~~~~~~~~~~~~~~

A list of all possible configuration options that can be specified in the
recipe is given for each diagnostic individually (see previous section).


Example plots
-------------

.. _fig_1:
.. figure::  /recipes/figures/lifetime/timeseries_trop.png
   :align:   center
   :width:   14cm

Time series of CH4 lifetime in the troposphere.

.. _fig_2:
.. figure::  /recipes/figures/lifetime/1d_profile.png
   :align:   center
   :width:   14cm

1D profile of CH4 lifetime.

.. _fig_3:
.. figure::  /recipes/figures/lifetime/zonal_mean_profile.png
   :align:   center
   :width:   14cm

Zonal mean profile of CH4 lifetime.
