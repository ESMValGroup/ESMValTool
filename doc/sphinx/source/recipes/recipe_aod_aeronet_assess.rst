.. _recipe_aod_aeronet_assess:

AOD AeroNET Assess
==================

Overview
--------

This diagnostic evaluates model aerosol optical depth (AOD) against ground
based observations from the AeroNET measurement network. Monthly mean AOD
data is downloaded from the AeroNET website and formatted (CMORized) using the
AERONET downloader and formatter within ESMValTool.

Multiannual seasonal means are calculated from the model output and compared
with a multiannual seasonal mean climatology generated from AeroNET
observational data. At each AeroNET station the data are screened for validity
according to the following default criteria:

  * 1. Monthly means must be generated from at least one AOD observation in that
    month.

  * 2. Seasonal means for DJF, MAM, JJA and SON must be calculated from three
    monthly means, i.e. a monthly mean from December January and Feburary.

  * 3. For a given year to be valid, there must be a seasonal mean for each climate
    season i.e. DJF, MAM, JJA and SON.

  * 4. For a multiannual seasonal means there must be at least five seasonaal means
    over the time range of interest.

NOTE: The code is designed to be flexible and the default criteria can be
changed according to the user's requirements (see the user settings below).

The evaluation is visualised by plotting model output as 2D filled contours and
overlaying AeroNET observations at model grid cells co-located with the AeroNET
measurement stations. Statistical data (root mean square error) is generated
using AeroNET observations at model grid cells co-located with the AeroNET
measurement stations.

Available recipes and diagnostics
---------------------------------

Recipes are stored in esmvaltool/recipes/

    * recipe_aod_aeronet_assess.yml

Diagnostics are stored in esmvaltool/diag_scripts/aerosols/

    * aod_aeronet_assess.py: Plot the AOD evaluation.
    * aero_utils.py: Utility functions commonly used by aerosol assessment routines.


User settings in recipe
-----------------------

#. Script aod_aeronet_assess.py

   *Required settings for script*

   * wavel: The wavelength of interest for the evaluation. Defaults to 440nm.
   * min_days_per_mon: The minimum number of days used to calculate the AOD monthly mean. Defaults to 1.
   * min_mon_per_seas: The minimum number of seasons used to calculate each
     seasonal mean. This must be between 1 and 3. Defaults to 3.
   * min_seas_per_year: The minimum number of seasonal means in each year. This
     must be between 1 and 4. Defaults to 4.
   * min_seas_per_clim: The minimum number of seasonal means used to calculate
     the multiannual seasonal mean. This must be btween 1 and the number of years
     of available AeroNET data. Defaults to 5.

   *Optional settings for script*

   * None

   *Required settings for variables*

   * None

   *Optional settings for variables*

   * None

   *Required settings for preprocessor*

   * None

   *Optional settings for preprocessor*

   * None

   *Color tables*

   * brewer_Spectral_11


Variables
---------

* od440aer (atmos, monthly mean, longitude latitude time)
* od550aer (atmos, monthly mean, longitude latitude time)
* od870aer (atmos, monthly mean, longitude latitude time)



Observations and reformat scripts
---------------------------------

* Note: (1) obs4MIPs data can be used directly without any preprocessing; (2)
  see headers of reformat scripts for non-obs4MIPs data for download
  instructions.

* The AeroNET data is downloaded from the AeroNET website using the downloader:

  .. code-block:: yaml

        $ esmvaltool data download AERONET.

* The AeroNET data is formatteed (CMORized) using the formatter:

  .. code-block:: yaml

        $ esmvaltool data format AERONET.



References
----------
* Holben B.N., T.F.Eck, I.Slutsker, D.Tanre, J.P.Buis, A.Setzer, E.Vermote, J.A.Reagan, Y.Kaufman, T.Nakajima, F.Lavenu, I.Jankowiak, and A.Smirnov, 1998: AERONET - A federated instrument network and data archive for aerosol characterization, Rem. Sens. Environ., 66, 1-16.

* Holben, B.N., D.Tanre, A.Smirnov, T.F.Eck, I.Slutsker, N.Abuhassan, W.W.Newcomb, J.Schafer, B.Chatenet, F.Lavenue, Y.J.Kaufman, J.Vande Castle, A.Setzer, B.Markham, D.Clark, R.Frouin, R.Halthore, A.Karnieli, N.T.O'Neill, C.Pietras, R.T.Pinker, K.Voss, and G.Zibordi, 2001: An emerging ground-based aerosol climatology: Aerosol Optical Depth from AERONET, J. Geophys. Res., 106, 12 067-12 097.

* Mulcahy, J. P., Johnson, C., Jones, C. G., Povey, A. C., Scott, C. E., Sellar, A., Turnock, S. T., Woodhouse, M. T., Abraham, N. L., Andrews, M. B., Bellouin, N., Browse, J., Carslaw, K. S., Dalvi, M., Folberth, G. A., Glover, M., Grosvenor, D. P., Hardacre, C., Hill, R., Johnson, B., Jones, A., Kipling, Z., Mann, G., Mollard, J., O’Connor, F. M., Palmiéri, J., Reddington, C., Rumbold, S. T., Richardson, M., Schitgens, N. A. J., Stier, P., Stringer, M., Tang, Y., Walton, J., Woodward, S., and Yool. A.: Description and evaluation of aerosol in UKESM1 and HadGEM3-GC3.1 CMIP6 historical simulations, Geosci. Model Dev., 13, 6383–6423, 2020

Example plots
-------------

.. _fig_aod_aeronet_assess_1:
.. figure::  /recipes/figures/aod_aeronet_assess/UKESM1-0-LL_CMIP_AERmon_historical_od440aer_gn_1994_2014_DJF.png
   :align:   center

   Evaluation of AOD at 440 nm from UKESM1 historical ensemble member r1i1p1f2 against the AeroNET climatology from ground-based observations for Dec-Jan-Feb. The multiannual seasonal mean is calculated for the model data for the period 1994-2014. The model output is overlaid with the observational climatology.

.. _fig_aod_aeronet_assess_2:
.. figure::  /recipes/figures/aod_aeronet_assess/UKESM1-0-LL_CMIP_AERmon_historical_od440aer_gn_1994_2014_MAM.png
   :align:   center

   Evaluation of AOD at 440 nm from UKESM1 historical ensemble member r1i1p1f2 against the AeroNET climatology from ground-based observations for Mar_Apr_May. The multiannual seasonal mean is calculated for the model data for the period 1994-2014. The model output is overlaid with the observational climatology.

.. _fig_aod_aeronet_assess_3:
.. figure::  /recipes/figures/aod_aeronet_assess/UKESM1-0-LL_CMIP_AERmon_historical_od440aer_gn_1994_2014_JJA.png
   :align:   center

   Evaluation of AOD at 440 nm from UKESM1 historical ensemble member r1i1p1f2 against the AeroNET climatology from ground-based observations for Jun-Jul-Aug. The multiannual seasonal mean is calculated for the model data for the period 1994-2014. The model output is overlaid with the observational climatology.

.. _fig_aod_aeronet_assess_4:
.. figure::  /recipes/figures/aod_aeronet_assess/UKESM1-0-LL_CMIP_AERmon_historical_od440aer_gn_1994_2014_SON.png
   :align:   center

   Evaluation of AOD at 440 nm from UKESM1 historical ensemble member r1i1p1f2 against the AeroNET climatology from ground-based observations for Sep-Oct-Nov. The multiannual seasonal mean is calculated for the model data for the period 1994-2014. The model output is overlaid with the observational climatology.

.. _fig_aod_aeronet_assess_5:
.. figure::  /recipes/figures/aod_aeronet_assess/UKESM1-0-LL_CMIP_AERmon_historical_od440aer_gn_1994_2014_scatter.png
   :align:   center

   Evaluation of AOD at 440 nm from UKESM1 historical ensemble member r1i1p1f2 against the AeroNET climatology from ground-based observations for Dec-Jan-Feb, Mar_Apr_May, Jun-Jul-Aug and Sep-Oct-Nov. The multiannual seasonal mean is calculated for the model data for the period 1994-2014.
