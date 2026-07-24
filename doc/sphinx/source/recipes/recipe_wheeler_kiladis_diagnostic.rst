.. _recipe_wheeler_kiladis_diagnostic:

Wheeler-Kiladis tropical wave spectra
=====================================

Overview
--------

This recipe computes normalized Wheeler-Kiladis wavenumber-frequency
spectra for tropical precipitation and outgoing longwave radiation.

The diagnostic is intended to evaluate tropical wave variability,
including Kelvin waves, equatorial Rossby waves, and broad intraseasonal
variability such as the Madden-Julian Oscillation.

The recipe uses standard ESMValTool preprocessors to convert units,
compute daily means, extract the tropical belt between 15S and 15N, and
regrid the data to a regular latitude-longitude grid. The diagnostic then
computes a cosine-latitude-weighted equatorial mean, removes leap days,
removes the time mean, linear trend, and annual harmonics, optionally
removes the zonal mean, computes two-dimensional Fourier spectra in time
and longitude, estimates a smoothed background spectrum, and plots
normalized spectra with optional shallow-water dispersion curves.

The diagnostic saves raw, background, and normalized Wheeler-Kiladis
spectra as NetCDF files and records provenance information for both data
outputs and figures.

Available recipes and diagnostics
---------------------------------

Recipes are stored in ``esmvaltool/recipes/``:

    * ``recipe_wheeler_kiladis_diagnostic.yml``

Diagnostics are stored in ``esmvaltool/diag_scripts/wheeler_kiladis/``:

    * ``wheeler_kiladis.py``: ESMValTool wrapper for running the
      Wheeler-Kiladis diagnostic.
    * ``spectra.py``: utilities for preparing equatorial fields and
      computing wavenumber-frequency spectra.
    * ``plot.py``: plotting utilities for normalized Wheeler-Kiladis
      spectra and theoretical dispersion curves.

User settings in recipe
-----------------------

#. Script ``wheeler_kiladis/wheeler_kiladis.py``

   *Optional settings for script*

   * ``annual_harmonics``: Number of annual harmonics removed during
     seasonal-cycle correction. Default: 3.

   * ``remove_zonal_mean``: If ``True``, remove the zonal mean before
     the wavenumber-frequency transform. This emphasizes propagating
     equatorial disturbances. Default: ``True``.

   * ``segment_length``: Length of each spectral segment in days.
     Default: 180.

   * ``segment_overlap``: Overlap between consecutive spectral segments
     in days. Default: 90.

   * ``sampling_frequency_per_day``: Sampling frequency of the input
     data in samples per day. For daily data, this should be 1.0.
     Default: 1.0.

   * ``sigma_freq``: Gaussian smoothing width in the frequency direction
     used to estimate the background spectrum. Default: 4.0.

   * ``sigma_wn``: Gaussian smoothing width in the zonal-wavenumber
     direction used to estimate the background spectrum. Default: 4.0.

   * ``max_wavenumber``: Maximum zonal wavenumber shown in the plot.
     Default: 15.

   * ``max_frequency``: Maximum frequency shown in the plot in cycles
     per day. Default: 0.5.

   * ``period_ticks``: Periods in days shown as labels on the frequency
     axis. Default: ``[2, 3, 5, 10, 20, 30, 60, 100]``.

   * ``equivalent_depths``: Equivalent depths in metres used to draw
     theoretical shallow-water dispersion curves. Default:
     ``[8, 12, 25, 50]``.

   * ``show_dispersion``: If ``True``, overlay theoretical Kelvin and
     equatorial Rossby wave dispersion curves. Default: ``True``.

   * ``mask_zero_wavenumber``: If ``True``, mask the
     zonal-wavenumber-zero column in the normalized spectrum plot. This
     is useful when the zonal mean has been removed. Default: ``True``.

Variables
---------

* pr (atmos, daily, longitude, latitude, time)
* rlut (atmos, daily, longitude, latitude, time)

Example plots
-------------

.. _wk_era5_pr_example:
.. figure:: /recipes/figures/wheeler_kiladis/ERA5_pr_wk_normalized_example.png
   :align: center

   Normalized Wheeler-Kiladis spectrum for ERA5 precipitation.

.. _wk_persiann_cdr_pr_example:
.. figure:: /recipes/figures/wheeler_kiladis/PERSIANN-CDR_pr_wk_normalized_example.png
   :align: center

   Normalized Wheeler-Kiladis spectrum for PERSIANN-CDR precipitation.

.. _wk_mpi_esm1_2_lr_rlut_example:
.. figure:: /recipes/figures/wheeler_kiladis/MPI-ESM1-2-LR_rlut_wk_normalized_example.png
   :align: center

   Normalized Wheeler-Kiladis spectrum for MPI-ESM1-2-LR outgoing
   longwave radiation.

References
----------

* Wheeler, M. C. and Kiladis, G. N. (1999): Convectively Coupled
  Equatorial Waves: Analysis of Clouds and Temperature in the
  Wavenumber-Frequency Domain.
* Kiladis, G. N. et al. (2009): Convectively Coupled Equatorial Waves.
* Hannah et al. (2020): Evaluation of tropical variability using
  Wheeler-Kiladis-type spectra.
