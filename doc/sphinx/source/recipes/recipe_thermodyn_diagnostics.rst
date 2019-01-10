Runoff_ET
=========

Overview
--------

The tool allows to compute TOA, atmospheric and surface energy budgets, latent energy and water mass budgets, 
meridional heat transports, the Lorenz Energy Cycle (LEC), the material entropy production with the direct 
and indirect method.

The energy budgets are computed from monthly mean radiative and heat fluxes at the TOA and at the surface 
(cfr. Wild et al., 2013). The meridional heat transports are obtained from the latitudinal integration 
of the zonal mean energy budgets. When a land-sea mask is provided, results are also available for 
land and oceans, separately.

The water mass budget is obtained from monthly mean latent heat fluxes (for evaporation), total and snowfall 
precipitation (cfr. Liepert et al., 2012). Latent energy budget is obtained by multiplying each component of 
the water mass budget by the respective latent heat constant.  When a land-sea mask is provided, results are 
also available for land and oceans, separately.

The LEC is computed from 3D fields of daily mean velocity and temperature fields in the troposphere over 
pressure levels. The analysis is carried on in spectral fields, converting lonlat grids in Fourier coefficients. 
The components of the LEC are computed as in Ulbrich and Speth, 1991.

The material entropy production is computed by using the indirect or the direct method (or both), the former 
relying on the convergence of radiative heat in the atmosphere (cfr. Lucarini et al., 2011; Pascale et al., 2011), 
the latter on the computation of entropy production  of all viscous and non-viscous dissipative processes occurring 
in the atmosphere (namely the sensible heat fluxes, the hydrological cycle with its components and the kinetic energy 
dissipation).

A comprehensive report of this new method is currently in preparation for Geosciences Model Development.



Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_Thermodyn_diagtool.yml

Diagnostics are stored in diag_scripts/Thermodyn_diagtool/

    * thermodyn_diagnostics.py: the main script, where computations are performed and the multi-model ensemble plots
				created

    * mkthe.py: an auxiliary script for the computation of LCL height, boundary layer top height and temperature, potential
		temperature

    * fourier_coefficients.py: an auxiliary script for the computation of the Fourier coefficients from the lonlat grid

    * lorenz_cycle.py: an auxiliary script for the computation of the LEC components in Fourier coefficients

    * fluxogram.py: an auxiliary script for the computation of the block diagram displaying the reservoirs and conversion terms
		    of the LEC

    * diagram_module.py: an auxiliary script containing information that is needed by fluxogram.py

    * plot_script.py: an auxiliary script for the computation of maps, scatter plots, time series and meridional sections of some 
		      derived quantities for each model in the ensemble. The meridional heat and water mass transports are also
		      computed here, as well as the peak magnitudes and locations;


User settings
-------------

recipe_Thermodyn_diagtool.yml

   *Optional settings for variables*

   * eb (y or n): this flag is set to 'y' for computation of TOA, atmospheric and surface energy budgets and meridional energy transports.
		  In the current version, such computations are compulsory for the succesful completion of the program
   * wat (y or n): if set to 'y', this flag allows for computation of the water mass and latent energy budgets and transports
   * lsm (y or n): if set to 'y', this flag allows for computation of the energy budgets, meridional energy transports, 
		   water mass and latent energy budgets and transports separately over land and oceans
   * lec (y or n): if set to 'y', this flag allows for computation of the LEC
   * entr (y or n): if set to 'y', this flag allows for computation of the material entropy production
   * met (1, 2 or 3): this flas specifies if the computation of the material entropy production must be performed with the indirect method 
		      (1), the direct method (2), or both methods. If 2 or 3 options are chosen, the intensity of the LEC is needed for the 
		      entropy production related to the kinetic energy dissipation. If lec is set to 'n', a default value is provided.

   These options apply to all models provided for the multi-model ensemble computations


Variables
---------

* hfls    (atmos,  monthly mean, time latitude longitude)
* hfss    (atmos,  monthly mean, time latitude longitude)
* hus     (atmos,  monthly mean, time plev latitude longitude)
* pr      (atmos,  monthly mean, time latitude longitude)
* prsn    (atmos,  monthly mean, time latitude longitude)
* ps      (atmos,  monthly mean, time latitude longitude)
* rlds    (atmos,  monthly mean, time latitude longitude)
* rlus    (atmos,  monthly mean, time latitude longitude)
* rlut    (atmos,  monthly mean, time latitude longitude)
* rsds    (atmos,  monthly mean, time latitude longitude)
* rsdt    (atmos,  monthly mean, time latitude longitude)
* rsus    (atmos,  monthly mean, time latitude longitude)
* rsut    (atmos,  monthly mean, time latitude longitude)
* ta      (atmos,  daily   mean, time plev latitude longitude)
* tas     (atmos,  monthly mean, time latitude longitude)
* ts      (atmos,  monthly mean, time latitude longitude)
* ua      (atmos,  daily   mean, time plev latitude longitude)
* va      (atmos,  daily   mean, time plev latitude longitude)
* wap     (atmos,  daily   mean, time plev latitude longitude)


References
----------
* Lembo V, Lunkeit F, Lucarini V (2019) A new diagnostic tool for diagnosing water, energy and entropy budgets in climate models. Geophys Mod Dev, in prep.
* Liepert BG, Previdi M (2012) Inter-model variability and biases of the global water cycle in CMIP3 coupled climate models. Environ Res Lett 7:014006. doi: 10.1088/1748-9326/7/1/014006
* Lorenz EN (1955) Available Potential Energy and the Maintenance of the General Circulation. Tellus 7:157–167. doi: 10.1111/j.2153-3490.1955.tb01148.x
* Lucarini V, Fraedrich K, Ragone F (2010) New Results on the Thermodynamical Properties of the Climate System. J Atmo 68:. doi: 10.1175/2011JAS3713.1
* Lucarini V, Blender R, Herbert C, et al (2014) Reviews of Geophysics Mathematical and physical ideas for climate science. doi: 10.1002/2013RG000446
* Pascale S, Gregory JM, Ambaum M, Tailleux R (2011) Climate entropy budget of the HadCM3 atmosphere–ocean general circulation model and of FAMOUS, its low-resolution version. Clim Dyn 36:1189–1206. doi: 10.1007/s00382-009-0718-1
* Ulbrich U, Speth P (1991) The global energy cycle of stationary and transient atmospheric waves: Results from ECMWF analyses. Meteorol Atmos Phys 45:125–138. doi: 10.1007/BF01029650
* Wild M, Folini D, Schär C, et al (2013) The global energy balance from a surface perspective. Clim Dyn 40:3107–3134. doi: 10.1007/s00382-012-1569-8


Example plots
-------------

.. _fig_1:
.. figure::  
   :align:   center
   :width:   14cm

   Insert caption here
