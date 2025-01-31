
.. _recipes_spei:

Droughts
========

Consecutive Dry Days (CDD)
--------------------------

Meteorological drought can in its simplest form be described by a lack of
precipitation. The Consecutive Dry Days (CDD) diagnostic calculates the longest
period frequency of dry days based on user defined thresholds.

More details and usage examples can be found in the
:ref:`CDD recipe documentation <recipe_consecdrydays>`.

Standardized Precipitation-Evapotranspiration Index (SPEI)
----------------------------------------------------------

Meteorological droughts are often described using the Standardized Precipitation
Index (SPI; McKee et al, 1993), which in a standardized way describes local
precipitation anomalies.

Because SPI does not account for evaporation from the ground, it lacks one it
lacks one component of the water fluxes at the surface and is therefore not
compatible with the concept of hydrological or agricultural drought. The
Standardized Precipitation-Evapotranspiration Index (SPEI; Vicente-Serrano et
al., 2010) has been developed to also account for temperature effects on the
surface water fluxes, by estimating the Potential Evapo-Transpiration (PET).

More details and usage examples can be found in the
:ref:`SPEI recipe documentation <recipe_spei>`.

Available recipes and diagnostics
---------------------------------


Recipes stored in ``recipes/droughts/``

    * recipe_cdd.yml
    * recipe_spi.yml
    * recipe_spei.yml


Diagnostics for index calculation:stored in ``diag_scripts/droughts/``

    * cdd.py: calculate Consecutive Dry Days
    * pet.R: calculate Potential Evapo-Transpiration
    * spi.R: calculate Standardized Precipitation Index
    * spei.R: calculate Standardized Evapo-Transpiration Index


User settings
-------------
The configuration that can be used in the recipe is documented in the
corresponding API documentation for python diagnostics (linked in diagnostic
list above). The R diagnostics are documented in
:ref:`SPEI recipe documentation <recipe_spei>`.



Variables
---------

* pr      (atmos, daily/monthly mean, time latitude longitude)
* tas     (atmos, monthly mean, time latitude longitude)
* tasmin     (atmos, monthly mean, time latitude longitude)
* tasmax     (atmos, monthly mean, time latitude longitude)
* sfcWind     (atmos, monthly mean, time latitude longitude)
* rsds     (atmos, monthly mean, time latitude longitude)
* clt    (atmos, monthly mean, time latitude longitude)
* hurs    (atmos, monthly mean, time latitude longitude)
* ps    (atmos, monthly mean, time latitude longitude)



References
----------
* McKee, T. B., Doesken, N. J., & Kleist, J. (1993). The relationship of drought frequency and duration to time scales. In Proceedings of the 8th Conference on Applied Climatology (Vol. 17, No. 22, pp. 179-183). Boston, MA: American Meteorological Society.

* Vicente-Serrano, S. M., Beguería, S., & López-Moreno, J. I. (2010). A multiscalar drought index sensitive to global warming: the standardized precipitation evapotranspiration index. Journal of climate, 23(7), 1696-1718.


Example plots
-------------

.. _fig_consecdrydays:
.. figure::  /recipes/figures/consecdrydays/consec_example_freq.png
   :align:   center
   :width:   14cm

   Example of the number of occurrences with consecutive dry days of more than five days in the period 2001 to 2002 for the CMIP5 model bcc-csm1-1-m.


.. _fig_spei:
.. figure::  /recipes/figures/spei/histogram_spei.png
   :align:   center
   :width:   14cm

   (top) Probability distribution of the standardized precipitation-evapotranspiration index of a sub-set of the CMIP5 models, and (bottom) bias relative to the CRU reference data set.

.. _fig_spi:
.. figure::  /recipes/figures/spei/histogram_spi.png
   :align:   center
   :width:   14cm

   (top) Probability distribution of the standardized precipitation index of a sub-set of the CMIP5 models, and (bottom) bias relative to the CRU reference data set.
