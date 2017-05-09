NCAR's Climate Variability Diagnostics Package (CVDP)
=====================================================



Overview
--------

**The NCAR CVDP version currently implemented into the ESMValTool is v4.1.**

The Climate Variability Diagnostics Package (CVDP) developed by NCAR's Climate
Analysis Section (Phillips et al., 2014) has been implemented into the
ESMValTool in order to be able to run it within this framework and alongside
the ESGF on CMIP output. CVDP can be used to evaluate the major modes of
climate variability including ENSO, PDO, AMO, Northern and Southern Annular
Modes (NAM and SAM), North Atlantic Oscillation (NAO), Pacific North and South
American teleconnection patterns (PNA and PSA). In addition it calculates
global trend maps and index time series for the above modes and the North
Pacific Index, the Tropical North Atlantic SST, Tropical South Atlantic SST,
Tropical Indian Ocean SST, Niño1+2, Niño3, and Niño4 times series, and the
Indian Ocean Dipole (IOD).

CVDP is developed as a standalone tool outside the ESMValTool. Once a new
version of CVDP is released, the ESMValTool will be updated
accordingly. Therefore, the structure of CVDP was kept as is and a wrapper has
been written to be able to run CVDP directly within the ESMValTool.



Available Namelists and Diagnostics
-----------------------------------

Namelists are stored in nml/

	* namelist_CVDP.xml

Diagnostics are stored in diag_scripts/

*Wrapper scripts to run CVDP within the framework of the ESMValTool*

	* cvdp_obs.ncl: run for each variable separately with observational data available; renames the ESMValTool output (observations) following the filename conventions of the CVDP and creates the CVDP namelist “namelist_obs”.
	* cvdp_ocean.ncl: renames the ESMValTool output (ocean variables) following the filename conventions of the CVDP.
	* cvdp_atmos.ncl: renames the ESMValTool output (atmosphere variables) following the filename conventions of the CVDP and creates the CVDP namelist “namelist” containing the models. The script then runs the CVDP via a call to the wrapper script cvdp_driver.ncl.



User settings
-------------

User setting files (cfg files) are stored in nml/cfg_CVDP/

     (1)	cvdp_obs.ncl

     *Required diag_script_info attributes*

	* obs_ref: list of reference data sets (observations) (array)

     (2)	cvdp_driver.ncl (called by cvdp_atmos.ncl)

     *The wrapper script cvdp_driver.ncl sets the user options for the CVDP.*



Variables
---------

TBD


Observations and Reformat Scripts
---------------------------------

TBD



References
----------

TBD


Example plots
-------------

TBD
















