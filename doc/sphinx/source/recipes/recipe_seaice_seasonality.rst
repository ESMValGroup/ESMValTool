.. _recipe_seaice_seasonality:

Antarctic sea ice seasonality
=============================

``recipe_seaice_seasonality.yml`` calculates the day of sea ice advance,
retreat, and season duration from daily CMIP6 ``SIday siconc``. It produces
NetCDF fields for each model and a polar map with one row per model and
common colour scales. The default example compares ACCESS-CM2 and
BCC-CSM2-MR for the 2000/01 sea ice year. Both models publish daily
``SIday siconc`` for this period; ACCESS-ESM1-5 does not.

The calculation follows `Massom et al. (2013)
<https://doi.org/10.1371/journal.pone.0064756>`_. The ice year begins on
15 February. Advance is the first day in the first run of at least five
consecutive days with concentration at or above 15%. Retreat is the first
day below 15% after the final ice-covered day; it is capped at the last
day of the ice year. Duration is retreat minus advance. Perennial ice is
assigned the first and last days of the ice year. A cell with no
sustained advance or any missing daily values is masked.

This convention differs slightly from the `COSIMA advanced recipe
<https://github.com/COSIMA/cosima-recipes/blob/main/03-Advanced-Recipes/Sea_Ice_Seasonality_Statistics.ipynb>`_: that notebook reports the
fifth day of the first qualifying run as advance and the final
ice-covered day as retreat. The ESMValTool diagnostic follows the
published event definitions and states its day numbering explicitly.

The input must contain one uninterrupted daily sample per day from
15 February to 14 February. Monthly sea ice concentration cannot be
used. The example recipe extracts through 15 February of the next year
to retain 14 February samples timestamped at midday; any extra next-year
15 February sample is removed before calculating the ice season.
Nearest-neighbour regridding can shift coastal and marginal ice
cells. Model calendars can have different numbers of days; compare
the ordinal maps with that difference in mind. The example recipe does
not yet include observational sea ice concentration, so it supports
model-to-model comparison rather than an observational skill score.
