Southern Ocean
==============

Overview
--------

These diagnostics include polar stereographic (difference) maps to compare the monthly/annual mean sea surface temperature, salinity and wind stress from ESMs with ERA-Interim data. Furthermore, there are scripts to plot the differences in the area mean vertical profiles of temperature and salinity between models and data from the World Ocean Atlas (Antonov et al., 2010; Locarnini et al., 2010). The ocean mixed layer thickness from models can be compared with that obtained from the Argo floats (Dong et al., 2008), again using polar stereographic (difference) maps. Finally, the Antarctic Circumpolar Current strength, as measured by water mass transport through the Drake Passage, is calculated using the same method as in the CDFTOOLS package (CDFtools). This diagnostic can be used to calculate the transport through other section as well, but is only available for EC-Earth/NEMO output for which all grid information is available.


Available Namelists and Diagnostics
-----------------------------------

Namelists are stored in nml/

* namelist_SouthernOcean.xml

Diagnostics are stored in diag_scripts/

* SouthernOcean_polcon.ncl: create polar stereographic plots for ocean mixed layer thickness, sea surface salinity and sea surface temperature.
* SouthernOcean_polcon_diff.ncl: create polar stereographic plots of the difference between individual models and reference data for ocean mixed layer thickness, sea surface temperature and eastward and northward wind stress. All data are regridded to a common grid using the ESMF regridding software.
* SouthernOcean_vector_polcon_diff.ncl: create polar stereographic contour plots of the difference between individual model data and reference data similar to SouthernOcean_polcon_diff.ncl, but on top plots vectors (magnitude and direction) for both the individual models and the reference data. Currently it is used for wind stress, but it should be possible to use it for other variables with u and v components as well. All data are regridded to a common grid using the ESMF regridding software.
* SouthernOcean_areamean_vertconplot.ncl: calculate the average sea water salinity and temperature over a subdomain from model data and create a Hovmoller-like diagram with time and depth on the axes. All data are regridded to a common grid using the ESMF regridding software.
* SouthernOcean_transport.ncl: calculate the sea water volume transport across a section from the variables uo and vo using a similar approach as is done in the CDFTOOLS package. **Currently only available for EC-Earth/Nemo output** as the calculations are performed using uo and vo on a staggered grid and the grid dimensions of the u and v grids are required.


User settings
-------------

TBD


Variables
---------

TBD


Observations and Reformat Scripts
---------------------------------

TBD



References
----------

* Antonov, J. I., D. Seidov, T. P. Boyer, R. A. Locarnini, A. V. Mishonov, H. E. Garcia, O. K. Baranova, M. M. Zweng, and D. R. Johnson (2010). World Ocean Atlas 2009, Volume 2: Salinity. S. Levitus, Ed. NOAA Atlas NESDIS 69, U.S. Government Printing Office, Washington, D.C., 184 pp (available at ftp://ftp.nodc.noaa.gov/pub/WOA09/DOC/woa09_vol2_text.pdf).
* CDFtools: http://servforge.legi.grenoble-inp.fr/projects/CDFTOOLS
* Dong, S., J. Sprintall, S. T. Gille, and L. Talley (2008). Southern Ocean mixed-layer depth from Argo float profiles, J. Geophys. Res., 113, C06013, doi: 10.1029/2006JC004051.
* ERA-Interim: http://www.ecmwf.int/en/research/climate-reanalysis/era-interim 
* ESMF regridding with ncl: http://www.ncl.ucar.edu/Applications/ESMF.shtml 
* Locarnini, R. A., A. V. Mishonov, J. I. Antonov, T. P. Boyer, H. E. Garcia, O. K. Baranova, M. M. Zweng, and D. R. Johnson (2010). World Ocean Atlas 2009, Volume 1: Temperature. S. Levitus, Ed. NOAA Atlas NESDIS 68, U.S. Government Printing Office, * Washington, D.C., 184 pp. (available at ftp://ftp.nodc.noaa.gov/pub/WOA09/DOC/woa09_vol1_text.pdf).
* nco: http://nco.sourceforge.net 
* World Ocean Atlas 2009: https://www.nodc.noaa.gov/OC5/WOA09/pr_woa09.html


Example plots
-------------

.. figure:: ../../source/namelists/figures/southern_ocean/fig1.png
   :scale: 50 %
   :alt: xxxx
 
.. figure:: ../../source/namelists/figures/southern_ocean/fig2.png
   :scale: 50 %
   :alt: xxxx

.. figure:: ../../source/namelists/figures/southern_ocean/fig3.png
   :scale: 50 %
   :alt: xxxx

.. figure:: ../../source/namelists/figures/southern_ocean/fig4.png
   :scale: 50 %
   :alt: xxxx

.. figure:: ../../source/namelists/figures/southern_ocean/fig5.png
   :scale: 50 %
   :alt: xxxx

.. figure:: ../../source/namelists/figures/southern_ocean/fig6.png
   :scale: 50 %
   :alt: xxxx

.. figure:: ../../source/namelists/figures/southern_ocean/fig7.png
   :scale: 50 %
   :alt: xxxx




