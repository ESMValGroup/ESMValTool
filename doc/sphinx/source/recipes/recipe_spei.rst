.. _recipes_spei:

Standardized Precipitation-Evapotranspiration Index (SPEI)
==========================================================

Overview
--------
Droughts can be separated into three main types: meteorological, hydrological, and agricultural drought.

Common for all types is that a drought needs to be put in context of local and seasonal characteristics, i.e. a drought should not be defined with an absolute threshold, but as an anomalous condition.

Meteorological droughts are often described using the standardized precipitation index (SPI; McKee et al, 1993), which in a standardized way describes local precipitation anomalies. It is calculated on monthly mean precipitation, and is therefore not accounting for the intensity of precipitation and the runoff process. Because SPI does not account for evaporation from the ground, it lacks one component of the water fluxes at the surface and is therefore not compatible with the concept of hydrological drought.

A hydrological drought occurs when low water supply becomes evident, especially in streams, reservoirs, and groundwater levels, usually after extended periods of meteorological drought. GCMs normally do not simulate hydrological processes in sufficient detail to give deeper insights into hydrological drought processes. Neither do they properly describe agricultural droughts, when crops become affected by the hydrological drought. However, hydrological drought can be estimated by accounting for evapotranspiration, and thereby estimate the surface retention of water. The standardized precipitation-evapotranspiration index (SPEI; Vicente-Serrano et al., 2010) has been developed to also account for temperature effects on the surface water fluxes. Evapotranspiration is not normally calculated in GCMs, so SPEI often takes other inputs to estimate the evapotranspiration. Here, the Thornthwaite (Thornthwaite, 1948) method based on temperature is applied.



Available recipes and diagnostics
---------------------------------

Recipes are stored in ``recipes/droughts/``

* recipe_spei.yml


Diagnostics are stored in ``diag_scripts/droughts/``

* pet.R: calculate the PET as input variable for SPEI
* spei.R: calculate the SPI and SPEI index

User settings
-------------

pet.R
~~~~~

pet_type: str
    The method used to calculate the potential evapotranspiration.
    Some settings and required variables depend on the calculation method.
    Options are: Penman, Thornthwaite, Hargreaves


spei.R
~~~~~~

    SPI is considered as a special case of the SPEI. Provide only precipitation
    as input variable or ancestor and set ``distribution="Gamma"``, to calculate
    SPI.


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
    First of reference period.
    By default first year of time series

refstart_month: int, optional
    first month of reference period.
    By default 1.

refend_year: int, optional
    Last year of the reference period.
    By default last year of time series.

refend_month``: integer, optional
    Last month of reference period.
    By default 12.


Variables
---------

* pr      (atmos, monthly mean, time latitude longitude)
* tas     (atmos, monthly mean, time latitude longitude)

Which variables are required or used for index calculation depends on the index
and calculation method of the PET. The Thornthwaite method requires tas, while
the Hargreaves method requires tas and sfcWind. The Penman method requires tas,
sfcWind, rsds, clt, hurs, and ps, but it is possible to approximate some of the
variables if not available.


References
----------
* McKee, T. B., Doesken, N. J., & Kleist, J. (1993). The relationship of drought frequency and duration to time scales. In Proceedings of the 8th Conference on Applied Climatology (Vol. 17, No. 22, pp. 179-183). Boston, MA: American Meteorological Society.

* Vicente-Serrano, S. M., Beguería, S., & López-Moreno, J. I. (2010). A multiscalar drought index sensitive to global warming: the standardized precipitation evapotranspiration index. Journal of climate, 23(7), 1696-1718.
