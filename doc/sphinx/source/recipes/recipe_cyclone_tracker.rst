.. _recipe_cyclone_tracker:

Cyclone tracker
===============


Overview
--------

This recipe prepares and formats the input variables to be read by the BSC cyclone tracker, 
a modified version of the GFDL vortex tracker V3.5b by the National Oceanic and Atmospheric Administration (NOAA).
Due to limitations in the tracker, each dataset has to contain at most six months of 6-hourly data.
If working with instantaneous 3-hourly data, the datasets can be resampled using the `resample_hours` preprocessor.
The recipe creates the namelist and tracker is run automatically if the path to the executable is provided. 
The executable can be obtained upon request. The tracker provides an estimate of the cyclone center position (latitude and longitude) 
along with metrics for intensity and structure at each time step.

Available recipes and diagnostics
---------------------------------

Recipes are stored in recipes/

    * recipe_cyclone_tracker.yml

Diagnostics are stored in diag_scripts/primavera/cyclone_tracker

    * cyclone_tracker.py: script to prepare the input and the namelist for the tracker.


User settings
-------------

script cyclone_tracker.py

    *Required settings for script*

    The following settings are required to write the namelist:

    * 'westbd': West boundary of the domain. Set to -1000 to include the whole domain.
    * 'eastbd': East boundary of the domain. Set to -1000 to include the whole domain.
    * 'northbd': North boundary of the domain. Set to -1000 to include the whole domain.
    * 'southbd': South boundary of the domain. Set to -1000 to include the whole domain.
    * 'type': Type of tracking to perform. Set to 'tcgen' for the BSC tracker. 
              The original GFDL tracker accepts other options.
    * 'mslpthresh': Minimum MSLP gradient threshold [hPa/km].
    * 'mslpthresh2': Maximum pressure threshold [hPa].
    * 'v850thresh': Minimum azimutally-average 850hPa windspeed threshold [m/s].
    * 'contint': Interval minimum pressure of a new low and the closed isobar[Pa].
    * 'wcore_depth': Contour interval [K].
    * 'verb': Verbosity of the tracker (0/1/2/3).
    * 'ikeflag': Compute Integrated Kinetic Energy (y/n).

    Other settings:

    * tracker_exe: Path to the tracker executable. If not provided the recipe will only save the namelist and the prepared input.

Variables
---------

* psl (atmos, 6hr or 3hr, time latitude longitude)
* ua (atmos, 6hr or 3hr, time plev latitude longitude)
* va (atmos, 6hr or 3hr, time plev latitude longitude)
* uas (atmos, 6hr or 3hr, time latitude longitude)
* vas (atmos, 6hr or 3hr, time latitude longitude)
* ta (atmos, 6hr or 3hr, time plev latitude longitude)
* zg (atmos, 6hr or 3hr, time plev latitude longitude)



References
----------

* Kreussler, P., Caron, L-P., Wild, S., Loosveldt Tomas, S., Chauvin, F., Moine, M-P., Roberts, M.J., 
  Ruprich‐Robert, Y., Seddon, J., Valcke, S., Vannière, B., Vidale, P.L., 2020.
  Tropical Cyclone Integrated Kinetic Energy in an Ensemble of HighResMIP Simulations. 
  Geophysical Research Letter, 48, 5.


