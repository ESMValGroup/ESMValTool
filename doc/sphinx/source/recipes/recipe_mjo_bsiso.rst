.. _recipe_mjo_bsiso:

Seasonal MJO and BSISO EEOF index
=================================

Overview
--------

``recipe_mjo_bsiso.yml`` ports the two ACCESS-NRI MJO/BSISO notebooks to
ESMValTool. The index diagnostic fits two area-weighted extended EOFs (EEOFs)
to daily outgoing longwave radiation (``rlut``) anomalies. A 139-tap
Hamming-window FIR filter retains periods of 25–90 days. Lagged fields at
days -10, -5, and 0 are assembled on uninterrupted daily data before DJF and
JJA training dates are selected. The two seasonal patterns are projected
onto each available day using the training PCA mean and normalization.

The mean-state diagnostic plots May–October precipitation and 850 hPa winds
from NCEP-NCAR-R1 for 2000–2020. Its two precipitation bias panels use
ACCESS-ESM1-5 historical 1985–2014 and ACCESS-CM2 piControl 1420–1449.
These periods and experiments differ, so the panels are descriptive
model-minus-reference comparisons.

The index and the Wheeler–Kiladis wavenumber-frequency spectra are
complementary diagnostics. The latter characterizes spectral power; it does
not provide the EEOF PCs and daily phase-space trajectory calculated here.

Available recipes and diagnostics
---------------------------------

* ``esmvaltool/recipes/recipe_mjo_bsiso.yml``
* ``esmvaltool/diag_scripts/mjo/bimodal_index.py``
* ``esmvaltool/diag_scripts/mjo/mean_state_bias.py``

The index diagnostic saves a NetCDF file containing seasonal EEOFs,
training and projected PCs, daily amplitudes, mode classifications, and
monthly occurrence estimates. It also saves EEOF maps, PC time series,
monthly occurrence bars, an amplitude comparison, and phase-space figures.
The bias diagnostic saves the plotted mean fields and biases as NetCDF and
a three-panel figure. All outputs include ESMValTool provenance.

User settings in recipe
-----------------------

The ``mjo_bsiso_index/index`` script accepts:

* ``low_period`` and ``high_period``: lower and upper filter periods in days
  (default 25 and 90).
* ``filter_window``: odd number of Hamming FIR taps (default 139).
* ``lags``: distinct signed lags in days (default ``[-10, -5, 0]``).

The recipe selects a 30-year ACCESS-ESM1-5 piControl sample to match the
source notebook. Set ``max_parallel_tasks`` to control memory use when
running the full example on Gadi.

Variables
---------

* ``rlut`` (atmosphere, daily, time × latitude × longitude): tropical OLR.
* ``pr`` (atmosphere, daily, time × latitude × longitude): precipitation.
* ``ua`` and ``va`` (atmosphere, daily, time × pressure × latitude ×
  longitude): winds selected at 850 hPa.

Scientific caveats
------------------

EOF signs are mathematically arbitrary. The BSISO EEOF2 and PC2 are flipped
together to follow the notebook convention, but absolute geographic phase
labels need independent comparison with a published reference. Accordingly,
the phase-space figure shows no geographic phase names. The daily amplitude
threshold of 1 and the larger-of-two-index mode classification reproduce the
source notebook; they do not uniquely identify the underlying physical
mechanism.

The FIR used here follows the source notebook's ``scipy.signal.firwin``
default Hamming window. The separate ESMValTool MJO Hovmöller diagnostic
uses a Lanczos filter, so its filtered anomalies will not be identical.

Example plots
-------------

.. figure:: /recipes/figures/mjo/mjo_bsiso_eeofs_access_esm1_5.png
   :align: center

   Seasonal lagged EEOF patterns for ACCESS-ESM1-5 piControl, 0296–0325.
   The BSISO modes are shown in EEOF1, EEOF2 order.

.. figure:: /recipes/figures/mjo/mjo_bsiso_mean_state_access.png
   :align: center

   NCEP-NCAR-R1 May–October precipitation and 850 hPa wind mean state,
   with the two model precipitation biases. The model periods and
   experiments differ and are labelled in the panels.

References
----------

* Kikuchi, K., Wang, B., and Kajikawa, Y. (2012), *Bimodal representation
  of the tropical intraseasonal oscillation*, Climate Dynamics 38,
  1989–2000, doi:10.1007/s00382-011-1159-1.
