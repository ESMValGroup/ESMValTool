.. _recipes_mjo:

MJO Wavenumber-Frequency Spectra
================================

Overview
--------

This recipe computes Wheeler-Kiladis space-time spectra to evaluate
Madden-Julian Oscillation (MJO)-related variability in the tropics.

The diagnostic calculates symmetric and anti-symmetric spectra in the
equatorial belt, derives a smoothed background spectrum, and produces
ratio spectra (signal/background). In addition, it computes seasonal
eastward/westward wavenumber-frequency spectra for boreal winter and summer.

The current recipe is configured for daily precipitation and zonal wind,
with zonal wind at 850 hPa and 200 hPa.


Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

* recipe_mjo.yml

Diagnostics are stored in esmvaltool/diag_scripts/mjo/

* mjo_diag.py: recipe entry point that loads each input cube and runs the
   spectra workflows.
* spectra_compute.py: Wheeler-Kiladis computations, smoothing,
   dispersion-curve overlays, plotting, and NetCDF output writing.


User settings in recipe
-----------------------

#. Script mjo_diag.py

    Required settings for script

    * script: ../diag_scripts/mjo/mjo_diag.py

    Optional (currently included in recipe_mjo.yml)

    * reference_datasets: list of dataset names used as references in metadata.
    * regrid_dataset: dataset used as target for regridding.
    * mip_name: project/mip label used in metadata.
    * base_range: baseline period.
    * analysis_range: analysis period.

    Required settings for variables

    * variable_group metadata is used internally to select processing branches.
    * eqw_spectra_pr requires variable pr with daily frequency.
    * eqw_spectra_winds requires variable ua with daily frequency.

    Optional settings for variables

    * Additional datasets/ensembles/experiments can be added in the datasets
       section of the recipe.

    Required settings for preprocessors

    * extract_region must provide a full cyclic longitude range and equatorial
       latitude coverage. The recipe uses:

       * start_longitude: 0.
       * end_longitude: 360.
       * start_latitude: -15.
       * end_latitude: 15.

    * pp_precip converts precipitation to mm day-1.
    * pp_levels_winds extracts 85000 Pa and 20000 Pa levels from ua.

    Notes and constraints

    * Input longitude must be periodic/cyclic.
    * Latitude and longitude coordinates must be 1D and monotonically ordered.
    * The algorithm requires at least 365 days and uses a 96-day analysis window.


Variables
---------

* pr (atmos, daily mean, longitude latitude time)
* ua (atmos, daily mean, longitude latitude time pressure)

Derived variable names in output are:

* Precipitation (from pr)
* x_wind_850hPa and x_wind_200hPa (from ua)


Observations and reformat scripts
---------------------------------

No dedicated observational reformat script is required by this recipe.
Any input dataset available through ESMValTool with the required variables,
daily frequency, and spatial/temporal coverage can be used.


References
----------

* Wheeler, M., and G. N. Kiladis (1999), Convectively Coupled Equatorial
   Waves: Analysis of Clouds and Temperature in the Wavenumber-Frequency
   Domain, Journal of the Atmospheric Sciences, 56, 374-399.
* Hayashi, Y. (1971), A Generalized Method of Resolving Disturbances into
   Progressive and Retrogressive Waves by Space and Fourier and Time Cross
   Spectral Analysis, Journal of the Meteorological Society of Japan, 49,
   125-128.


Example plots
-------------

The diagnostic produces several plot families per dataset and variable,
including:

* Raw anti-symmetric spectra
* Raw symmetric spectra
* Ratio anti-symmetric spectra (raw/background)
* Ratio symmetric spectra (raw/background)
* Seasonal (winter/summer) wavenumber-frequency spectra

At the moment, no static example image is shipped in
doc/sphinx/source/recipes/figures/mjo/.
