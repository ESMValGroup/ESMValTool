.. _recipe_ref_fire:

Climate drivers of fire
=======================

Overview
--------

This diagnostic includes:

* burnt fraction
* fire weather control
* fuel load/continuity control

The diagnostic relies on the processing of fire climate drivers through the
ConFire model (dedicated branch from the GitHub repository available at
https://github.com/douglask3/Bayesian_fire_models/tree/AR7_REF) and is based on
`Jones et al. (2024)`. The diagnostic computes the burnt fraction for each grid
cell based on a number of drivers. Additionally, the respective controls due to
fire weather and fuel load/continuity are computed. The stochastic control
corresponds to the unmodelled processed influencing to fire occurrence.
The ESMValTool diagnostic includes only the relevant part of the evaluation code
(see https://github.com/douglask3/Bayesian_fire_models/blob/AR7_REF/fire_model/ConFire.py
for now) as the model run is done offline beforehand. The corresponding result
files are made available through a Zenodo archive which can be retrieved inside
the diagnostic. The archive is published under the DOI 10.5281/zenodo.14917244,
but a specific version of the files can be specified using the corresponding DOI,
e.g. the default link to https://zenodo.org/records/14917245. The other option
is to provide a local directory which contains the necessary model files.

The ConFire model relies on a variety of observational datasets (see references):

* Global Fire Emissions Database version 5
* MODIS MOD44B
* ESA CCI Biomass
* ISIMIP3a GSWP3-W5E5 dataset

Note: If custom ConFire run files are to be used for the diagnostic,
they need to be present inside the auxiliary directory defined in the
user configuration file `config_user.yml`.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/:

* recipe_ref_fire.yml

Diagnostics are stored in esmvaltool/diag_scripts/fire/:

* fire_diagnostic.py: main diagnostic script calling a util function from diagnostic_run_ConFire.py.
* diagnostic_run_ConFire.py: script containing utils functions to run the ConFire model.


User settings in recipe
-----------------------

#. Script fire_diagnostic.py

   *Required settings for script*

   * var_order: list of climate drivers in the order corresponding to the one
     specified in the corresponding file from the confire_param directory.

   *Optional settings for script*

   * confire_param: path to the directory containing the required files to run
     the ConFire model or Zenodo URL to retrieve files from a Zenodo archive.
     If custom files are used, the corresponding directory needs to be present
     inside the auxiliary data directory defined inside the user configuration.
     This defaults to the original Zenodo archive otherwise.
   * remove_vpd_files: Removing or not the computed vapor pressure deficit files.
     It will only apply if the vapor pressure deficit is part of var_order.
     This defaults to False.
   * remove_confire_files: Removing or not the files produced during the ConFire
     model evaluation.
     This defaults to False.

   *Required settings for variables*

   *Optional settings for variables*

   *Required settings for preprocessor*

   *Optional settings for preprocessor*

#. Script diagnostic_run_ConFire.py

   *Required settings for script*

   * confire_param: path to the directory containing the required files to run
     the ConFire model (or downloaded files from a Zenodo archive).
   * files_input: list containing tuples of variable name and file path for each
     climate driver present in var_order.

   *Optional settings for script*

   * model_name: string containing the input data model for legending the plots.
   * timerange: string containing the time range of the input data for legending the plots.
   * project: string or list containing the input data project(s).
   * experiment: string or list containing the input data experiment(s).

   *Required settings for variables*

   *Optional settings for variables*

   *Required settings for preprocessor*

   *Optional settings for preprocessor*


Variables
---------

* pr (atmos, monthly mean, longitude latitude time)
* tasmax (atmos, monthly mean, longitude latitude time)
* treeFrac (land, monthly mean, longitude latitude time)
* vegFrac (land, monthly mean, longitude latitude time)
* cveg (land, monthly mean, longitude latitude time)
* tas (atmos, monthly mean, longitude latitude time): used to compute vpd.
* hurs (atmos, monthly mean, longitude latitude time): used to compue vpd.


References
----------

* Jones, M. W., Kelley, D. I., Burton, C. A., Di Giuseppe, F., Barbosa, M. L. F.,
  Brambleby, E., Hartley, A. J., Lombardi, A., Mataveli, G., McNorton, J. R.,
  Spuler, F. R., Wessel, J. B., Abatzoglou, J. T., Anderson, L. O., Andela, N.,
  Archibald, S., Armenteras, D., Burke, E., Carmenta, R., Chuvieco, E., Clarke, H.,
  Doerr, S. H., Fernandes, P. M., Giglio, L., Hamilton, D. S., Hantson, S.,
  Harris, S., Jain, P., Kolden, C. A., Kurvits, T., Lampe, S., Meier, S., New, S.,
  Parrington, M., Perron, M. M. G., Qu, Y., Ribeiro, N. S., Saharjo, B. H.,
  San-Miguel-Ayanz, J., Shuman, J. K., Tanpipat, V., van der Werf, G. R.,
  Veraverbeke, S., and Xanthopoulos, G.: State of Wildfires 2023–2024,
  Earth Syst. Sci. Data, 16, 3601–3685, https://doi.org/10.5194/essd-16-3601-2024, 2024.

* Yang Chen, Joanne Hall, Dave van Wees, Niels Andela, Stijn Hantson, Louis Giglio,
  Guido R. van der Werf, Douglas C. Morton, & James T. Randerson. (2023).
  Global Fire Emissions Database (GFED5) Burned Area (0.1) [Data set]. Zenodo.
  https://doi.org/10.5281/zenodo.7668424.

* DiMiceli, C., Sohlberg, R., Townshend, J. (2022). MODIS/Terra Vegetation Continuous
  Fields Yearly L3 Global 250m SIN Grid V061 [Data set]. NASA EOSDIS Land Processes
  Distributed Active Archive Center. Accessed 2025-04-01 from https://doi.org/10.5067/MODIS/MOD44B.061.

* Santoro, M.; Cartus, O. (2024): ESA Biomass Climate Change Initiative (Biomass_cci):
  Global datasets of forest above-ground biomass for the years 2010, 2015, 2016,
  2017, 2018, 2019, 2020 and 2021, v5.01. NERC EDS Centre for Environmental Data
  Analysis, 22 August 2024. https://dx.doi.org/10.5285/bf535053562141c6bb7ad831f5998d77.

* Stefan Lange, Matthias Mengel, Simon Treu, Matthias Büchner (2022): ISIMIP3a atmospheric
  climate input data (v1.0). ISIMIP Repository. https://doi.org/10.48364/ISIMIP.982724.


Example plots
-------------

.. _fig_ref_fire_burnt_area:
.. figure::  /recipes/figures/ref/burnt_fraction_MPI-ESM1-2-LR_historical_2013_2014.png
   :align:   center

   Burnt area fraction for the MPI-ESM1-2-LR model (CMIP-historical experiment)
   for the time period 2013-2014 as computed with the ConFire model `Jones et al. (2024)`.

.. _fig_ref_fire_fire_weather_control:
.. figure::  /recipes/figures/ref/fire_weather_control_MPI-ESM1-2-LR_historical_2013_2014.png
   :align:   center

   Fire weather control for the MPI-ESM1-2-LR model (CMIP-historical experiment)
   for the time period 2013-2014 as computed with the ConFire model `Jones et al. (2024)`.

.. _fig_ref_fire_fuel_load_continuity_control:
.. figure::  /recipes/figures/ref/fuel_load_continuity_control_MPI-ESM1-2-LR_historical_2013_2014.png
   :align:   center

   Fuel load continuity control for the MPI-ESM1-2-LR model (CMIP-historical experiment)
   for the time period 2013-2014 as computed with the ConFire model `Jones et al. (2024)`.
