TEMPLATE FOR NAMELIST
=====================

Overview
--------

The aerosol diagnostics currently implemented allow for three kinds of
comparisons:

* *station data*: concentrations of various aerosol species (aerosol sulfate,
  aerosol nitrate, aerosol ammonium, black carbon, organic carbon, PM2.5 and
  PM10) are compared with observational data from station networks (IMPROVE and
  CASTNET for North America, EMEP for Europe and EANET for Asia). Aerosol
  optical depth (AOD) at 550 nm is compared with AERONET data. The comparison
  between model and observations is performed considering all available
  observational data in the selected time period, on a monthly-mean basis. The
  model data is extracted in the grid boxes where the respective observational
  stations are located (co-located model and observational data). This
  diagnostic produces three types of plot: time series (model/observations vs.
  time), scatter plot (model versus observations) and a map plot (model as
  contour map, observations as dots using the same color coding).

* *satellite data*: AOD at 550 nm is compared with satellite data (MODIS and
  MISR). This diagnostic produces contour map plots of the AOD of model and
  observations, as well as difference plot between model and observations.

* *vertical profiles, size distributions*: vertical profiles of aerosol (mass
  and number) concentrations and of aerosol size distributions are compared to
  observational data from aircraft campaigns and ground based stations. The
  model data are extracted based on the campaign/station location (lat-lon box)
  and time period (on a climatological basis, i.e. selecting the same
  days/months, but regardless of the year). The basic statistics is calculated
  on the selected model data (mean, standard deviation, median,
  5-10-25-75-90-95 percentiles) and compared to the same quantities from the
  observations (where available). This diagnostics supports monthly-mean,
  daily-mean and instantaneous input data, with the latter two providing a more
  robust evaluation. For this diagnostic, rather specific variables are
  required (i.e., aerosol number concentration for particles with diameter
  larger than 14 nm) to match the properties of the instruments used during the
  campaign. New CMOR variables have been added and corresponding EMAC recipes
  have been defined.


Available Namelists and Diagnostics
-----------------------------------

Namelists are stored in nml/

* namelist_aerosol_CMIP5.xml
* namelist_aerosol_EMAC.xml

Diagnostics are stored in diag_scripts/

* aerosol_stations.ncl
* aerosol_satellite.ncl
* aerosol_profiles.ncl


User settings
-------------

TBD


Variables
---------

TBD


Observations and Reformat Scripts
---------------------------------

TBD



References
----------

TBD


Example plots
-------------

TBD
















