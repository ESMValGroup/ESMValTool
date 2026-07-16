
.. _recipes_droughts:

Droughts
========

Consecutive Dry Days (CDD)
--------------------------

Meteorological drought can in its simplest form be described by a lack of
precipitation. The Consecutive Dry Days (CDD) diagnostic calculates the longest
period of dry days based on user defined thresholds.

More details and usage examples can be found in the
:ref:`CDD recipe documentation <recipes_consecdrydays>`.


Potential or Reference Evapotranspiration (PET, ET0)
-----------------------------------------------------

The Potential Evapotranspiration (PET) is a measure of the evaporative demand
of the atmosphere. It represents the amount of water that would evaporate from
a reference surface, i.e. fully watered grass land. ``pet.R`` is able to
calculate PET based on a method of users choice (``pet_type``). The methods
require different input variables for example:

- Penman: tas, tasmin, tasmax, sfcWind, rsds, clt, hurs, ps
- Hargreaves: tasmin, tasmax, (rsds, pr)
- Thornthwaite: tas

A complete list and more details can be found in
:ref:`SPEI recipe documentation <recipes_spei>`.


Standardized Precipitation-Evapotranspiration Index (SPEI)
----------------------------------------------------------

Meteorological droughts are often described using the Standardized Precipitation
Index (SPI; McKee et al, 1993), which in a standardized way describes local
precipitation anomalies.

Because SPI does not account for evaporation from the ground, it lacks one
component of the water fluxes at the surface and is therefore not compatible
with the concept of hydrological or agricultural drought. The Standardized
Precipitation-Evapotranspiration Index
(SPEI; Vicente-Serrano et al., 2010) has been developed to also account for
temperature effects on the surface water fluxes, by estimating the
Potential Evapotranspiration (PET).

More details and usage examples can be found in the
:ref:`SPEI recipe documentation <recipes_spei>`.


Available recipes and diagnostics
---------------------------------

.. toctree::
   :maxdepth: 1

   droughts/recipe_consecdrydays
   droughts/recipe_spei
   droughts/recipe_martin18grl


Diagnostics are stored in ``diag_scripts/droughts/``. An incomplete list of
diagnostics that are used in different recipes is shown below. Some recipes
might use additional diagnostics. See the corresponding recipe documentation
for more details.:

- cdd.py: calculate Consecutive Dry Days (CDD)
- pet.R: calculate Potential Evaporanspiration (PET)
- spei.R: calculate Standardized Precipitation-Evapotranspiration Index (SPEI)


User settings
-------------
The configuration that can be used in the recipes is documented in the
corresponding recipe or diagnostics API documentation linked in the list above.
:ref:`SPEI recipe documentation <recipes_spei>`.


Variables
---------

- pr      (atmos, daily/monthly mean, time latitude longitude)
- tas     (atmos, monthly mean, time latitude longitude)
- tasmin     (atmos, monthly mean, time latitude longitude)
- tasmax     (atmos, monthly mean, time latitude longitude)
- sfcWind     (atmos, monthly mean, time latitude longitude)
- rsds     (atmos, monthly mean, time latitude longitude)
- rsdt     (atmos, monthly mean, time latitude longitude)
- clt    (atmos, monthly mean, time latitude longitude)
- hurs    (atmos, monthly mean, time latitude longitude)
- ps    (atmos, monthly mean, time latitude longitude)


References
----------

- Martin, E.R. (2018). Future Projections of Global Pluvial and Drought Event
  Characteristics. Geophysical Research Letters, 45, 11913-11920.
- McKee, T. B., Doesken, N. J., & Kleist, J. (1993). The relationship of drought
  frequency and duration to time scales. In Proceedings of the 8th Conference on
  Applied Climatology (Vol. 17, No. 22, pp. 179-183). Boston, MA: American
  Meteorological Society.
- Vicente-Serrano, S. M., Beguería, S., & López-Moreno, J. I. (2010). A
  multiscalar drought index sensitive to global warming: the standardized
  precipitation evapotranspiration index. Journal of climate, 23(7), 1696-1718.
