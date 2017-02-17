:mod:`xy_line`
==============
.. function:: profile_plev(wks_in[1], source, varname[1]: string)

   :param integer wks_in: workstations (graphic object or default will be used).
   :param integer source: data to be plotted or a NetCDF filename: Must have "plev" dimension & can have "models", "quantity". @ptop: controls pressure axis range (Default: full range). @pbot: controls pressure axis range (Default: full range). @zoom: controls x axis range (Default: "yes"). @font: font type to use @Refmodel: reference model (Default: first in models dimension) @long_name: long variable name (Default: var) @short_name: short variable name (Default: var) @units: variable units (Default: missing)
   :param  string varname: variable name, needed for netCDF files with multiple variables

   Source prototype
  
   Return value
      A graphic variable.
  
   Description
      Creates a plot of profile(s) with vertical pressure axis.
      Adds percentiles or stddev of first model (if available) as whiskers.
      Opens default wks, if not provided as argument of type "graphic".
      Defines default ressources, which are overridden by argument res.
      Creates plot, according to wks & res.
  
   Caveats
      Treatment of coordinate "quantity" not very general yet.
      Selection of defaults for res almost arbitrary.
      Please check results of all scripts that use this routine if
      modifying the defaults!
  
   Modification history
      * 20140214-A_gott_kl: written.
  
.. function:: aerosol_profile(wks_in[1], source_m, source_o, varname[1])

   :param integer wks_in: workstations (graphic object or default will be used).
   :param integer source_m: model data to be plotted or a NetCDF filename with data.
   :param integer source_o: observations data to be plotted or a NetCDF filename with data.
   :param integer varname: variable name in the file.

   Source prototype
      source[*][*][*]
      source!0 = model
      source!1 = statistic
      source!2 = plev or diam
  
   Return value
      A graphic variable.
  
   Description
      Creates a plot of vertical profile (level vs. data) or size distribution
      (data vs. diameter).
      Plots median, mean or both depending on availability of obervations.
  
   Caveats
  
   Modification history
      20140917-A_righ_ma: renamed to aerosol_profile and extended to for plotting
                       size distributions.
      20140705-A_righ_ma: written.
  
.. function::  xy_line(wks[1], source, source_x, source_stddev, res_in : logical, debuginfo[1] : logical)

   :param integer wks:  workstation, must be passed - no default used yet!
   :param integer source:        data to be plotted (no netCDF input possible yet)
   :param integer source_x:      x-axis of array to be plotted (e.g. source&time, ... )
   :param integer source_stddev: standard deviation of input, needed if diag_script_info@multi_model_mean is set to "y"
   :param  logical res_in:  diag_script-specific resources passed from diag_script
   :param  logical debuginfo:  description about diagnostic rendered onto plot

   Source prototype
  
   Description
      Defines default ressources, which are overridden by argument res.
      Creates an xy-plot, according to wks & res.
      Adds multi model mean and standard deviation if
      diag_script_info@multi_model_mean is set to "y".
  
   Caveats
  
   Modification history
      20140109-A_senf_da: written.
  
.. function:: timeseries_station(wks_in[1], source, varname[1]: string)

   :param integer wks_in: workstations (graphic object or default will be used).
   :param integer source: data to be plotted or a NetCDF filename with data. @stname: station name @stlat: station latitude @stlon: station longitude @stalt: station altitude
   :param  string varname: variable name, needed for netCDF files with multiple variables

   Source prototype
      source[*][*]
      source!0 = model
      source!1 = time or year
  
   Return value
      A graphic variable.
  
   Description
      Creates a time series plot for station data.
  
   Caveats:
      * selection of defaults for res almost arbitrary
      * Please check results of all scripts that use this routine if
        modifying the defaults!
  
   Modification history:
      * 20140325 written by Mattia Righi
  
.. function:: cycle_plot(wks_in[1], source, varname[1] : string)

   :param integer wks_in: workstations (graphic object or default will be used).
   :param integer source: data to be plotted or a NetCDF filename with data.
   :param  string varname: variable name in the file.

   Source prototype
      source[*,*]
      source!0 = model
      source!1 = month or season
  
   Return value
      A graphic variable.
  
   Description
      Draw an annual or seasonal cycle plot.
  
   Caveats
  
   Modification history
      20131206-A_fran_fr: written.
  
