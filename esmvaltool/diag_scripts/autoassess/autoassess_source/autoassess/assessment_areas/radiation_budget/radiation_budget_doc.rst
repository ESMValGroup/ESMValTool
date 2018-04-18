===========================
Radiation Budget
===========================

======
Status
======
In development

===================
Assessment Overview
===================

The aim of monitoring for energy budget is understanding the (im)balance of energy flux between atmosphere and surface of model because it is linked with hydrological cycle and climate change. 
This radiation budget is defined which factors have more uncertainties compared with observation like CERES EBAF and statistics of the other models and figure out where this imbalance of fluxes comes from. 
It is available to use climate and numerial model, especially with forecast lead time such as T+24, T+72 and T+120. 
Also, it can be divided and calculated as seasonly mean like jja, djf, mam and son, and annualy mean.


=============================
Name of owners and developers
=============================

Owners
------
Prior to 2014 : Adapted from a routine written by M. Mizielinski and based on method developed by M.-E. Demory et al, 2014, Clim Dyn

Aug 2016 : B. Vanniere at Reading Uni.

Sep 2016 - Oct 2017 : Jiyoung Oh



Developers
----------

Prior to 2014 : M. Mizielinski and M.-E. Demory

Aug 2016 : B. Vanniere at Reading Uni.

Sep 2016 - Oct  2017 : Jiyoung Oh


=======================
Metrics and Diagnostics
=======================

Performance metrics:

* radiation_net_toa : difference between incoming shortwave radiation and the sum of outgoing short/longwave radiation at TOA
* toa_incoming_shortwave_flux : incoming shortwave radiation at TOA
* toa_outgoing_shortwave_flux : outgoing shortwave radiation from surface
* surface_downwelling_shortwave_flux_in_air : shortwave absorbed by surface plus the amount of reflected by surface 
* surface_net_downward_shortwave_flux : the value of incoming shortwave minus the reflection of cloud and absorbtion by atmosphere
* upward_sw_reflected_surface : total downward shortwave minus net downward shortwave in surface
* sw_reflected_clouds : difference between outgoing shortwave at TOA and the amount of reflected shortwave from surface
* sw_absorbed_atm : shortwave flux absorbed by atmosphere
* toa_outgoing_longwave_flux : outgoing longwave at TOA
* surface_downwelling_longwave_flux_in_air : back longwave radiation absorbed by surface
* surface_net_downward_longwave_flux : net longwave flux reached the surface 
* upward_lw_emitted_surface : surface downwelling longwave minus net downward longwave in surface
* net_surface_radiation : the summation of net downward short/longwave flux in surface
* surface_upward_sensible_heat_flux : sensible flux from surface
* surface_upward_latent_heat_flux : latent flux from surface
* radiation_adsorbed_surface : difference net surface radiation against sensible and latent heat flux
* sw_cloud_forcing : difference of outgoing shortwave between clear sky and all sky
* lw_cloud_forcing : same as sw cloud forcing but longwave instead of shortwave
* toa_outgoing_shortwave_flux_assuming_clear_sky : outgoing shortwave flux in clear sky at TOA
* toa_outgoing_longwave_flux_assuming_clear_sky : outgoing longwave flux in clear sky at TOA

Diagnostics:

* supermeans radiation flux, sensible/latent heat flux

==========
Model Data
==========

==========================================   ============= ============== ==============================================
          Variable/Field name                   realm        frequency                        Comment
==========================================   ============= ============== ==============================================
toa_incoming_shortwave_flux                    Atmosphere    supermean    - These nine variables having CF-name are
toa_outgoing_shortwave_flux                    Atmosphere    supermean         are necessary for calculation of the
surface_downwelling_shortwave_flux_in_air      Atmosphere    supermean         radiation budget. 
surface_net_downward_shortwave_flux            Atmosphere    supermean    - The others not having CF-name can get the
toa_outgoing_longwave_flux                     Atmosphere    supermean         combination of two/three variables.  
surface_downwelling_longwave_flux_in_air       Atmosphere    supermean    - CF name with "in_air" 
surface_net_downward_longwave_flux             Atmosphere    supermean
surface_upward_sensible_heat_flux              Atmosphere    supermean
surface_upward_latent_heat_flux                Atmosphere    supermean
sw_cloud_forcing                               Atmosphere    supermean
lw_cloud_forcing                               Atmosphere    supermean
toa_outgoing_sw_flux_assuming_clear_sky        Atmosphere    supermean
toa_outgoing_lw_flux_assuming_clear_sky        Atmosphere    supermean
==========================================   ============= ============== ==============================================


==========
References
==========

Graeme L. Stephens et al.(2012),An update on Earth's energy balance in light of the latest global observations,Nature Geoscience 5, 691-696(2012),doi:10.1038/ngeo1580

M.S.Demory et al. (2014), The role of horizontal resoluiton in simulating drivers of the global hydrological cycle,Clim Dyn(2014) 42:2201-2225, doi:10.1007/s00382-013-1924-4


======================
Observations Data sets
======================

How to obtain the data
----------------------
Demory_et_al_2014_obs_Energy_Budget.txt from Figure 2 in M.S.Demory et al. (2014)

Stephens_et_al_2012_obs_Energy_Budget.txt from Figure B1 in Graeme L. Stephens et al.(2012)

CERES ebaf data from /project/earthobs/CERES/EBAF_CMIP5/supermeansnc.new in Met Office.
(The data are from the code, which is in ~hadac/usr/scripts/ebaf_supermeans.sh)

========================
Sample Plots and metrics
========================

=============================================  ======================
                Metric name                      Value [units: W/m2]
=============================================  ======================
radiation_net_toa                                       0.6
toa_incoming_shortwave_flux                           340.2
toa_outgoing_shortwave_flux                           100.0
toa_outgoing_shortwave_flux_assuming_clear_sky         50.2
total_sw_cloud_forcing                                 47.5
surface_downwelling_shortwave_flux_in_air             188.0
surface_net_downward_shortwave_flux                   165.0
upward_sw_reflected_surface                            23.0
sw_reflected_clouds                                    74.7
sw_absorbed_atm                                        80.0
toa_outgoing_longwave_flux                            239.7
toa_outgoing_longwave_flux_assuming_clear_sky         266.4
total_lw_cloud_forcing                                 26.7
surface_downwelling_longwave_flux_in_air              345.6
net_downward_lw_surface                               -52.4
upward_lw_emitted_surface                             398.0
net_surface_radiation                                 112.6
surface_upward_sensible_heat_flux                      24.0
surface_upward_latent_heat_flux                        88.0
radiation_adsorbed_surface                              0.6
=============================================  ======================

.. figure:: images/Ebudget_{MODEL}_{RES}_{YEARs}_{YEARe}.png
      :scale: 100 %
   :alt: Ebudget_{MODEL}_{RES}_{YEARs}_{YEARe}.png

   Bias of radiation/heat fluxes in global.
