.. _recipes_aod_aeronet_assess:

AOD AeroNET Assess
==================

Overview
--------

This diagnostic evaluates model aerosol optical depth (AOD) against ground
based observations from the AeroNET measurement network. Multiannual seasonal
means are calculated from the model output and compared with a multiannual
seasonal mean climatology generated from AeroNET observational data.

The evaluation is visualised by plotting model output as 2D filled contours and
overlaying AeroNET observations at model grid cells co-located with the AeroNET measurement stations. Statistical data (root mean square error) is generated
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

   * wavel: The wavelength of interest for the evaluation, currently set up for 440nm.
   * aeronet_dir: String. Location of the observational AeroNET climatology.

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


Observations and reformat scripts
---------------------------------

*Note: (1) obs4MIPs data can be used directly without any preprocessing;
(2) see headers of reformat scripts for non-obs4MIPs data for download
instructions.*

* Pre-processed AOD climatologies are located in ''. Observations from AeroNET
stations are saved in individual NetCDF files with filenames of the format 'OBS_AERONET_ground_STATION_NAME_AERmon_od440aer_199301-201401.nc'.


References
----------
* Holben B.N., T.F.Eck, I.Slutsker, D.Tanre, J.P.Buis, A.Setzer, E.Vermote, J.A.Reagan, Y.Kaufman, T.Nakajima, F.Lavenu, I.Jankowiak, and A.Smirnov, 1998: AERONET - A federated instrument network and data archive for aerosol characterization, Rem. Sens. Environ., 66, 1-16.

* Holben, B.N., D.Tanre, A.Smirnov, T.F.Eck, I.Slutsker, N.Abuhassan, W.W.Newcomb, J.Schafer, B.Chatenet, F.Lavenue, Y.J.Kaufman, J.Vande Castle, A.Setzer, B.Markham, D.Clark, R.Frouin, R.Halthore, A.Karnieli, N.T.O'Neill, C.Pietras, R.T.Pinker, K.Voss, and G.Zibordi, 2001: An emerging ground-based aerosol climatology: Aerosol Optical Depth from AERONET, J. Geophys. Res., 106, 12 067-12 097.

* Mulcahy, J. P., Johnson, C., Jones, C. G., Povey, A. C., Scott, C. E., Sellar, A., Turnock, S. T., Woodhouse, M. T., Abraham, N. L., Andrews, M. B., Bellouin, N., Browse, J., Carslaw, K. S., Dalvi, M., Folberth, G. A., Glover, M., Grosvenor, D. P., Hardacre, C., Hill, R., Johnson, B., Jones, A., Kipling, Z., Mann, G., Mollard, J., O’Connor, F. M., Palmiéri, J., Reddington, C., Rumbold, S. T., Richardson, M., Schitgens, N. A. J., Stier, P., Stringer, M., Tang, Y., Walton, J., Woodward, S., and Yool. A.: Description and evaluation of aerosol in UKESM1 and HadGEM3-GC3.1 CMIP6 historical simulations, Geosci. Model Dev., 13, 6383–6423, 2020

Example plots
-------------

.. _fig_aod_aeronet_assess_1:
.. figure::  /recipes/figures/aod_aeronet_assess/UKESM1-0-LL_CMIP_AERmon_historical_od440aer_gn_1988_2008_DJF.png
   :align:   center

   Evaluation of AOD at 440 nm from UKESM1 historical ensemble member r1i1p1f2 against the AeroNET climatology from ground-based observations for Dec-Jan-Feb. The multiannual seasonal mean is calculated for the model data for the period 1998-2008. The model output is overlaid with the observational climatology.

.. _fig_aod_aeronet_assess_2:
.. figure::  /recipes/figures/aod_aeronet_assess/UKESM1-0-LL_CMIP_AERmon_historical_od440aer_gn_1988_2008_MAM.png
   :align:   center

   Evaluation of AOD at 440 nm from UKESM1 historical ensemble member r1i1p1f2 against the AeroNET climatology from ground-based observations for Mar_Apr_May. The multiannual seasonal mean is calculated for the model data for the period 1998-2008. The model output is overlaid with the observational climatology.

.. _fig_aod_aeronet_assess_3:
.. figure::  /recipes/figures/aod_aeronet_assess/UKESM1-0-LL_CMIP_AERmon_historical_od440aer_gn_1988_2008_JJA.png
   :align:   center

   Evaluation of AOD at 440 nm from UKESM1 historical ensemble member r1i1p1f2 against the AeroNET climatology from ground-based observations for Jun-Jul-Aug. The multiannual seasonal mean is calculated for the model data for the period 1998-2008. The model output is overlaid with the observational climatology.

.. _fig_aod_aeronet_assess_4:
.. figure::  /recipes/figures/aod_aeronet_assess/UKESM1-0-LL_CMIP_AERmon_historical_od440aer_gn_1988_2008_SON.png
   :align:   center

   Evaluation of AOD at 440 nm from UKESM1 historical ensemble member r1i1p1f2 against the AeroNET climatology from ground-based observations for Sep-Oct-Nov. The multiannual seasonal mean is calculated for the model data for the period 1998-2008. The model output is overlaid with the observational climatology.

.. _fig_aod_aeronet_assess_5:
.. figure::  /recipes/figures/aod_aeronet_assess/UKESM1-0-LL_CMIP_AERmon_historical_od440aer_gn_1988_2008_scatter.png
   :align:   center

   Evaluation of AOD at 440 nm from UKESM1 historical ensemble member r1i1p1f2 against the AeroNET climatology from ground-based observations for Dec-Jan-Feb, Mar_Apr_May, Jun-Jul-Aug and Sep-Oct-Nov. The multiannual seasonal mean is calculated for the model data for the period 1998-2008.
