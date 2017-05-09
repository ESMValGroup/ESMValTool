South Asian Monsoon
=====================

Overview
--------

Three South Asian Summer Monsoon (SASM) namelists for the basic climatology, seasonal cycle, intra-seasonal and inter-annual variability and key teleconnections have been implemented into the ESMValTool focusing on SASM rainfall and horizontal winds in June-September (JJAS). The goal is to provide a comprehensive overview of the basic features of the South Asian Monsoon for a certain model compared to other models, observational and reanalysis data sets for a range of processes and variables relevant to the monsoon.


Available Namelists and Diagnostics
-----------------------------------

Namelists are stored in nml/

* namelist_SAMonsoon_AMIP.xml
* namelist_SAMonsoon_daily.xml
* namelist_SAMonsoon.xml

Diagnostics are stored in diag_scripts/

* SAMonsoon_precip_basic.ncl
* SAMonsoon_precip_daily.ncl
* SAMonsoon_precip_domain.ncl
* SAMonsoon_precip_IAV.ncl
* SAMonsoon_precip_propagation.ncl
* SAMonsoon_precip_seasonal.ncl
* SAMonsoon_teleconnections.ncl
* SAMonsoon_wind_basic.ncl
* SAMonsoon_wind_IAV.ncl
* SAMonsoon_wind_seasonal.ncl


User settings
-------------

TBD


Variables
---------

* pr (atmos, daily/monthly, longitude latitude time)
* ts (atmos, monthly, longitude latitude time)
* ua (atmos, monthly, longitude latitude lev time)
* va (atmos, monthly, longitude latitude lev time)


Observations and Reformat Scripts
---------------------------------

**Note**: (1) obs4mips data can be used directly without any preprocessing; (2) see headers of reformat scripts for non-obs4mips data for download instructions.

* ERA-Interim (ua, va, pr – reformat_scripts/obs/reformat_obs_ERA-Interim.ncl, reformat_obs_ERA-Interim_surffluxes.ncl)
* NCEP (ua, va – reformat_scripts/obs/reformat_obs_NCEP.ncl)
* TRMM-L3 (pr, monthly means – obs4mips)
* TRMM_3B42 (pr, daily means – reformat_scripts/obs/reformat_obs_TRMM-3B42-daily.ncl)
* HadISST (ts – reformat_scripts/obs/reformat_obs_HadISST.ncl)
* GPCP-1DD (pr, daily means – obs4mips)
* GPCP-SG (pr, monthly means – obs4mips)
* MERRA (pr – obs4mips)
* CMAP (pr – reformat_scripts/obs/reformat_obs_CMAP.ncl)



References
----------

* Sperber, K. R., et al., The Asian summer monsoon: an intercomparison of CMIP5 vs. CMIP3 simulations of the late 20th century, Clim Dyn (2013) 41:2711–2744, doi: 10.1007/s00382-012-1607-6, 2012.
* Lin, Jia-Lin, Klaus M. Weickman, George N. Kiladis, Brian E. Mapes, Siegfried D. Schubert, Max J. Suarez, Julio T. Bacmeister, Myong-In Lee, 2008: Subseasonal Variability Associated with Asian Summer Monsoon Simulated by 14 IPCC AR4 Coupled GCMs. J. Climate, 21, 4541-4567. doi: http://dx.doi.org/10.1175/2008JCLI1816.1.
* Webster, P. J., and S.Yang, 1992: Monsoon and ENSO: Selectively interactive systems. Quart. J. Roy. Meteor. Soc., 118, 877-926. (Webster-Yang dynamical monsoon index)
* Goswami, B. N., B. Krishnamurthy, and H. Annamalai, 1999: A broad-scale circulation index for interannual variability of the Indian summer monsoon. Quart. J. Roy. Meteor. Soc., 125, 611-633. (Goswami dynamical monsoon index)
* Wang, B., and Z. Fan, 1999: Choice of south Asian summer monsoon indices. Bull. Amer. Meteor. Soc., 80, 629-638. (Wang-Fan dynamical monsoon index)
* Wang B., J. Liu, H. J. Kim, P. J. Webster, and S. Y. Yim, Recent change of global monsoon precipitation (1979-2008), Climate Dynamics, doi: 10.1007/s00382-011-1266-z, 2011. (Intensity/Monsoon domain reference)


Example plots
-------------

.. figure:: ../../source/namelists/figures/south_asian_monsoon/fig1.png
   :scale: 50 %
   :alt: xxxx
   

.. figure:: ../../source/namelists/figures/south_asian_monsoon/fig2.png
   :scale: 50 %
   :alt: xxxx


.. figure:: ../../source/namelists/figures/south_asian_monsoon/fig3.png
   :scale: 50 %
   :alt: xxxx


.. figure:: ../../source/namelists/figures/south_asian_monsoon/fig4.png
   :scale: 50 %
   :alt: xxxx



.. figure:: ../../source/namelists/figures/south_asian_monsoon/fig5.png
   :scale: 50 %
   :alt: xxxx


.. figure:: ../../source/namelists/figures/south_asian_monsoon/fig6.png
   :scale: 50 %
   :alt: xxxx


.. figure:: ../../source/namelists/figures/south_asian_monsoon/fig7.png
   :scale: 50 %
   :alt: xxxx



.. figure:: ../../source/namelists/figures/south_asian_monsoon/fig8.png
   :scale: 50 %
   :alt: xxxx


.. figure:: ../../source/namelists/figures/south_asian_monsoon/fig9.png
   :scale: 50 %
   :alt: xxxx


.. figure:: ../../source/namelists/figures/south_asian_monsoon/fig10.png
   :scale: 50 %
   :alt: xxxx
   

.. figure:: ../../source/namelists/figures/south_asian_monsoon/fig11.png
   :scale: 50 %
   :alt: xxxx


