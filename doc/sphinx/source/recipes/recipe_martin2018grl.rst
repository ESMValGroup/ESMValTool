.. _recipes_martin2018grl:

Drough characteristics following Martin (2018)
==============================================

Overview
--------


Following `Martin (2018)`_ drough characteristica are calculated based on the standard precipitation index (SPI), see `Mckee et al. (1993)`_. These characteristics are number of drought event, average duration, SPI index and severity index of these events.

.. _`Martin (2018)`: https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2018GL079807
.. _`Mckee et al. (1993)`: https://www.nature.com/articles/nclimate3387


Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

   * recipe_martin2018grl.yml


Diagnostics are stored in diag_scripts/

   * droughtindex/diag_save_spi.R
   * droughtindex/collect_drought_obs_multi.py
   * droughtindex/collect_drought_model.py
   * droughtindex/collect_drought_func.py


User settings in recipe
-----------------------

The recipe can be run with different CMIP5 and CMIP6 models.



Variables
---------

* *pr* (atmos, monthly, longitude, latitude, time)
* *tas* (atmos, monthly, longitude, latitude, time)


Observations and reformat scripts
---------------------------------

*None*


References
----------

* Martin, E.R. (2018). Future Projections of Global Pluvial and Drought Event Characteristics. Geophysical Research Letters, 45, 11913-11920.

* McKee, T. B., Doesken, N. J., & Kleist, J. (1993). The relationship of drought frequency and duration to time scales. In Proceedings of the 8th Conference on Applied Climatology (Vol. 17, No. 22, pp. 179-183). Boston, MA: American Meteorological Society.

Example plots
-------------

.. _martin2018grl_fig1:
.. figure:: /recipes/figures/droughtindex/martin2018grl_fig1.png
   :align: center
   :width: 50%

   Obs_multi

.. _martin2018grl_fig1:
.. figure:: /recipes/figures/droughtindex/martin2018grl_fig2b.png
   :align: center
   :width: 50%

   model comparison.


