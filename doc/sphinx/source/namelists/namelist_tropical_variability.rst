Tropical variability
====================

Overview
--------

The available diagnostics are motivated by the work of Li and Xie (2014). In particular, this diagnostics reproduces their Fig. 5 for models and observations/reanalyses, calculating equatorial mean (5°N-5°S), longitudinal sections of annual mean precipitation (pr), skin temperature (ts), horizontal winds (ua and va) and 925 hPa divergence (derived from the sum of the partial derivatives of the wind components extracted at the 925 hPa pressure level (that is du/dx + dv/dy). Latitude cross sections of the model variables are plotted for the equatorial Pacific, Indian and Atlantic Oceans with observational constraints provided by the TRMM-3B43-v7 for precipitation, the HadISST for SSTs, and ERA-interim reanalysis for temperature and winds. Latitudinal sections of absolute and normalized annual mean SST and precipitation are also calculated spatially averaged for the three ocean basins. Normalization follows the procedure outlined in Fig. 1 of Li and Xie (2014) whereby values at each latitude are normalized by the tropical mean (20°N-20°S) value of the corresponding parameter (e.g., annual mean precipitation at a given location is divided by the 20°N-20°S annual mean value). Finally, to assess how models capture observed relationships between SST and precipitation the co-variability of precipitation against SST is calculated for specific regions of the tropical Pacific. This analysis includes calculation of the Mean Square Error (MSE) between model SST/precipitation and observational equivalents.


Available Namelists and Diagnostics
-----------------------------------

Namelists are stored in nml/

* namelist_TropicalVariability.xml

Diagnostics are stored in diag_scripts/

* TropicalVariability.py
* TropicalVariability_EQ.py
* TropicalVariability_wind.py


User settings
-------------

TBD


Variables
---------

* ts: skin temperature (atmos, monthly mean, time latitude longitude)
* pr: precipitation (atmos, monthly mean, time latitude longitude)
* ua: u-wind (atmos, monthly mean, time plevel latitude longitude)
* va: v-wind (atmos, monthly mean, time plevel latitude longitude)


Observations and Reformat Scripts
---------------------------------

**Note: (1)** obs4mips data can be used directly without any preprocessing; (2) see headers of reformat scripts for non-obs4mips data for download instructions.

* HadISST: skin Temperature (ts) / sea surface temperature (SST)
  Reformat script: reformat_scripts/obs/reformat_obs_HadISST.ncl
* TRMM-L3 (pr, monthly means – obs4mips)
* ERA-Interim (u-wind, v-wind)
  Reformat script: reformat_scripts/obs/reformat_obs_ERA-Interim.ncl

References
----------

* Li, G., and , S.-P. Xie (2014), Tropical Biases in CMIP5 Multimodel Ensemble: The Excessive Equatorial Pacific Cold Tongue and Double ITCZ Problems. J. Climate, 27, 1765-1780. doi: http://dx.doi.org/10.1175/JCLI-D-13-00337.1.


Example plots
-------------


.. figure:: ../../source/namelists/figures/tropical_variability/fig1.png
   :scale: 50 %
   :alt: xxxx
   

.. figure:: ../../source/namelists/figures/tropical_variability/fig2.png
   :scale: 50 %
   :alt: xxxx
  
  
.. figure:: ../../source/namelists/figures/tropical_variability/fig3.png
   :scale: 50 %
   :alt: xxxx
   













