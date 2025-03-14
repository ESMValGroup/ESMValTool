.. _recipes_spei:

Standardized Precipitation-Evapotranspiration Index (SPEI)
==========================================================

Overview
--------
Droughts can be separated into three main types: meteorological, hydrological,
and agricultural drought.
Common for all types is that a drought needs to be put in context of local and
seasonal characteristics, i.e. a drought should not be defined with an absolute
threshold, but as an anomalous condition.

Meteorological droughts are often described using the
Standardized Precipitation Index (SPI; McKee et al, 1993), which in a
standardized way describes local precipitation anomalies.
It is calculated on monthly mean precipitation, and is therefore not accounting
for the intensity of precipitation and the runoff process.
Because SPI does not account for evaporation from the ground, it lacks one
component of the water fluxes at the surface and is therefore not compatible
with the concept of hydrological drought.

A hydrological drought occurs when low water supply becomes evident, especially
in streams, reservoirs, and groundwater levels, usually after extended periods
of meteorological drought.

GCMs normally do not simulate hydrological processes in sufficient detail to
give deeper insights into hydrological drought processes.
Neither do they properly describe agricultural droughts, when crops become
affected by the hydrological drought.
However, hydrological drought can be estimated by accounting for
evapotranspiration, and thereby estimate the surface retention of water.
The standardized precipitation-evapotranspiration index
(SPEI; Vicente-Serrano et al., 2010) has been developed to also account for
temperature effects on the surface water fluxes.
To consider evapotranspiration and plant water stress Potential
Evapotranspiration (PET) is calculated based on atmospheric variables.
Different methods to derive PET are described below.


This page documents a set of R diagnostics based on the
`SPEI.R library <https://CRAN.R-project.org/package=SPEI>`_.
``recipes/roughts/recipe_spei.yml`` is an example how to calculate and plot
SPEI using ``diag_scripts/droughts/pet.R`` and ``diag_scripts/droughts/spei.R``.



pet.R
-----

The Potential Evapo-Transpiration (PET) is a measure of the evaporative demand
of the atmosphere. It represents the amount of water, that would evaporate from
a reference surface, i.e. an fully watered alfalfa field. ``pet.R`` is able to
calculate PET based on a method of users choice (``pet_type``) using the
`SPEI.R package <https://doi.org/10.32614/CRAN.package.SPEI>`_. The
approximations require different input variables. To controll which variables
are available to which diagnostic script in a more complex recipe they can be
set explicitly as ancestors.

- Thornthwaite: tas
- Hargreaves: tasmin, tasmax, (rsdt, pr)
- Penman: tasmin, tasmax, sfcWind, ps, rsds, (rsdt, clt, hurs)

The Thornthwaite equation (Thornthwaite, 1948) is the simplest one based solely
on temperature. Hargreaves (1994) provides an equation based on daily minimum
(tasmin) and maximum temperature (tasmax) and external radiation (rsdt).
If precipitation data (pr) is provided and `use_pr: TRUE` it will be used as a
proxy for irradation to correct PET following Droogers and Allen (2002).
The Penman-Monteith formular additionally considers surface windspeed (sfcWind),
pressure (ps), and relative humidity (hurs). Some of these variables can be
approximated if not available (for example by providing clt instead of rsds).
There are further modifications to the Penman-Monteith equation, that can be
selected using the ``method`` key for ``pet_type: Penman``. Details about
the different method can be found in the SPEI.R package documentation
Beguería and Vicente-Serrano (2011).


User settings
~~~~~~~~~~~~~

pet_type: str
    The method used to calculate the potential evapotranspiration.
    Some settings and required variables depend on the calculation method.
    Options are: Penman, Thornthwaite, Hargreaves

use_pr: boolean, optional
    Use precipitation as proxy for irradation to correct PET. Only used for
    `pet_type: Hargreaves`.
    By default FALSE.

method: str, optional
    Method used for PET calculation. Only used for `pet_type: Penman`.
    Options are: ICID, FAO, ASCE.
    By default: ICID

crop: str, optional
    Crop type for PET calculation. Only used for `pet_type: Penman`.
    Options are: short, tall.
    By default: tall


spei.R
------

The Standardized Precipitation-Evapotranspiration Index (SPEI) is calculated by
fitting a probability distribution to the accumulated water budget
(pr-PET) for each grid cell and month of the year. The accumulation period is
configurable using the ``smooth_month`` setting.
The PDF is transformed to a normal distribution. Based on given percentiles,
index values from -2 (extreme droughts) to 2 (extreme wet spells) are assigned.
The distribution used to fit the water budget can be configured through the
``distribution`` setting. The required PET can be provided by setting the
pet.R diagnostic as an ancestor or adding an ``evspsblpot`` variable. To use
alternative variables (i.e. actual evapotranspiration) change the
``short_name`` setting accordingly.

    SPI is considered as a special case of the SPEI. Provide only precipitation
    as input variable or ancestor and set ``distribution="Gamma"``, to calculate
    SPI.

User settings
~~~~~~~~~~~~~

smooth_month: int
    The number of months (scale) to accumulate. Common choices are 3, 6, 9, 12.

write_coeffs: boolean, optional
    Save fitting coefficients.
    By default FALSE.

write_wb: boolean, optional
    Write water balance to netcdf file.
    By default FALSE.

short_name_pet: string, optional
    Short name of the variable to use as PET.
    By default "evspsblpot"

distributionn: string, optional
    Type of distribution used for SPEI calibration.
    Possible options are: "Gamma", "log-Logistic", "Pearson III".
    By default "log-Logistic".

refstart_year: int, optional
    First year of the reference period.
    By default first year of time series

refstart_month: int, optional
    First month of reference period.
    By default 1.

refend_year: int, optional
    Last year of the reference period.
    By default last year of time series.

refend_month: integer, optional
    Last month of reference period.
    By default 12.


References
----------
* McKee, T. B., Doesken, N. J., & Kleist, J. (1993). The relationship of drought frequency and duration to time scales. In Proceedings of the 8th Conference on Applied Climatology (Vol. 17, No. 22, pp. 179-183). Boston, MA: American Meteorological Society.

* Vicente-Serrano, S. M., Beguería, S., & López-Moreno, J. I. (2010). A multiscalar drought index sensitive to global warming: the standardized precipitation evapotranspiration index. Journal of climate, 23(7), 1696-1718.

* Beguería, S., & Vicente-Serrano, S. M. (2011). SPEI: Calculation of the Standardized Precipitation-Evapotranspiration Index (p. 1.8.1) [Dataset]. https://doi.org/10.32614/CRAN.package.SPEI

* Thornthwaite, C. W., (1948). An approach toward a rational classification of climate. Geogr. Rev., 38, 55-94. https://doi. org/10.1097/00010694-194807000-00007

* Hargreaves G.H., (1994). Defining and using reference evapotranspiration. Journal of Irrigation and Drainage Engineering 120: 1132-1139.

* Droogers P., Allen R. G., (2002). Estimating reference evapotranspiration under inaccurate data conditions. Irrigation and Drainage Systems 16: 33-45.

* Monteith, J.L., 1965. Evaporation and Environment. 19th Symposia of the Society for Experimental Biology, University Press, Cambridge, 19:205-234.


Example plots
-------------

.. _fig_spei_fig1:
.. figure:: /recipes/figures/droughts/spi_example.png
   :align: center
   :width: 80%

   Example plot of SPI averaged over the year 2005. The reference period for
   index calibration is 2000-2005.

.. _fig_spei_fig2:
.. figure:: /recipes/figures/droughts/spei_example.png
   :align: center
   :width: 80%

   Example plot of SPEI averaged over the year 2005. The reference period for
   index calibration is 2000-2005.
