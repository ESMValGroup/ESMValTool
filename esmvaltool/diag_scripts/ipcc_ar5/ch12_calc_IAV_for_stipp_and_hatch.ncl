; #############################################################################
;; INTERANNUAL VARIABILITY MULTI-MODEL MEAN FOR STIPPLING SIGNIFICANCE IPCCCH12
;; Author: Ruth Lorenz (ETH Zurich, Switzerland)
;  CRESCENDO project
; #############################################################################
; 
;  Description
;     Calculate Interannual variability from piControl runs for plotting
;     significance with stippling and hatching, regrid to common grid and save
;     multi-model mean.
;     Can either be calculated over full time period of piControl run
;     (no periodlength given)
;     or first averaged annually/seasonally over periodlength and then standard
;     deviation calculated based on averaged periods.
; 
;  Required diag_script_info attributes (diagnostics specific)
;     time_avg:  time period to calculate IAV, e.g. annualclim, seasonalclim
; 
;  Optional diag_script_info attributes (diagnostic specific)
;     periodlength: length of periods to calculate IAV accross, depends on
;                   period lengths used in
;                   map_diff_mmm_stippandhatch_ipcc12.ncl
;                   if not given whole time period calculated at once
;     grid:         grid for regridding (coarsest, finest, ref, 1x1, 2.5x2.5)
;     ref_model:    name of reference data set (only needed if grid = "ref")
; 
;  Required variable_info attributes (variable specific)
;     att1: short description
;     att2: short description
; 
;  Optional variable_info attributes (variable specific)
;     att1: short description
;     att2: short description
; 
;  Caveats
;     List possible caveats or limitations of this diagnostic
;     Features to-be-implemented shall also be mentioned here
; 
;  Modification history
;     20161219-A_lore_ru: remove seasonal cycle before std if seasonal
;     20161024-A_lore_ru: adapted to ESMValTool
;     20130501-A_sedl_ja: written for IPCC AR5 as get_natvar.ncl.
; 
; #############################################################################

; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
load "./interface_scripts/interface.ncl"

load "./diag_scripts/shared/statistics.ncl"

begin
  ; ##############################################################
  ; # Fetch general parameters, set in namelist_collins13ipcc.xml#
  ; # passed via environment variables by python code            #
  ; ##############################################################
  enter_msg(diag_script, "")

  ;  'diag_script' and 'variables' are fetched from the
  ;  above 'ncl.interface' file
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (var: " + variables(0) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ;  'models', 'variables' and 'field_types' are fetched from the
  ;  above 'ncl.interface' file
  dim_MOD = dimsizes(model_info@model)
  dim_VAR = dimsizes(variables)
  var0 = variables(0)
  field_type0 = field_types(0)

  ; ###########################################
  ; # References Information                  #
  ; ###########################################
  write_references(diag_script, \  ; script name
                   "A_lore_ru", \  ; authors
                   "", \           ; contributors
                   "D_collins13ipcc",    \  ; diag_references
                   "",          \  ; obs_references
                   "P_crescendo")   ; proj_references
  ; Check required diag_script_info attributes
  req_atts = (/"time_avg"/)
  exit_if_missing_atts(diag_script_info, req_atts)
end


begin
  ; Basename of diag_script
  diag_script_base = basename(diag_script)

  ; Output netcdf directory
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)  ; Needed?
end

begin
  ; #############################################
  ; # Get parameters from ./variable_defs/*.ncl #
  ; # passed via the 'info' attribute           #
  ; #############################################
  if (isvar("MyParam")) then
    delete(MyParam)
  end if
  if (isatt(variable_info, "long_name")) then
    MyParam = variable_info@long_name
    log_debug(" MyParam = " + MyParam)
  else
    error_msg("f", diag_script, "", "missing variable attribute " + \
              "'long_name'")
  end if
end

begin
  ; ###########################################
  ; # Get data and average time               #
  ; ###########################################
  ; get data from first model
  imod = 0  ; NCL array indicies start from zero
  log_debug("processing " + model_info@model(imod) + "_"\
            + model_info@exp(imod) + "_" \
            + model_info@ensemble(imod))
  outfile = "IAV_mmm_piControl_" + field_type0 + "_" \
            + var0 + "_" + diag_script_info@time_avg + \
            ".nc"

  ; See ./interface_scripts/data_handling.ncl
  A0 = read_data(imod, var0, field_type0)

  ; Check dimensions
  dims = getvardims(A0)
  rank = dimsizes(dims)
  if (dimsizes(dims) .lt. 3) then
    error_msg("f", diag_script, "", dimsizes(dims) + " dimensions, " + \
	      "need 3 or 4")
  end if
  idx = ind(dims .eq. "lat" .or. dims .eq. "rlat" .or. dims .eq. "j")
  if (ismissing(idx)) then
    error_msg("f", diag_script, "", "no lat dimension")
  end if
  lat_name = dims(idx)

  idx = ind(dims .eq. "lon" .or. dims .eq. "rlon" .or. dims .eq. "i")
  if (ismissing(idx)) then
    error_msg("f", diag_script, "", "no lon dimension")
  end if
  lon_name = dims(idx)

  index = ispan(0, dim_MOD - 1, 1)
  if (isatt(diag_script_info, "seasons")) then
    dim_seas = dimsizes(diag_script_info@seasons)
  else
    dim_seas = 1
    diag_script_info@seasons = (/0/)
  end if

  if (rank .eq. 4) then
    idx = ind(dims .eq. "lev" .or. dims .eq. "plev")
    if (ismissing(idx)) then
      error_msg("f", diag_script, "", "no level dimension")
    end if
    lev_name = dims(idx)
  end if

  ;  determine destination grid (regridding.ncl)
  if (.not.(isatt(diag_script_info, "grid"))) then
    error_msg("w", diag_script, "", \
	      "Regridding new grid not defined. " + "Using 2.5x2.5")
    selected_grid = "2.5x2.5"
  else
    selected_grid = diag_script_info@grid
  end if

  if (selected_grid .eq. "ref") then
    ; set reference model
    ref_ind = ind(model_info@model .eq. diag_script_info@ref_model(0))
    ref_model = diag_script_info@ref_model
    if (ismissing(ref_ind)) then
	if (diag_script_info@ref_model(0) .eq. "n.a.") then
	    ref_ind = 0
	else
	    error_msg("f", diag_script, "", \
		      "no adequate reference model provided")
	end if
    end if
    dest_grd = find_destination_grid(ref_ind, var0, \
				     field_type0, "coarsest")
  else if (selected_grid .eq. "1x1") then
    if (any(dims .eq. "plev")) then
	dest_grd     = new((/17, 180, 360/), float)
	dest_grd!0   = "plev"
	dest_grd!1   = "lat"
	dest_grd!2   = "lon"
	dest_grd&plev = (/100000., 92500., 85000., 70000., 60000., \
			 50000., 40000., 30000., 25000., 20000., 15000., \
			 10000., 7000., 5000., 3000., 2000., 1000./)
	dest_grd&lat = flt2dble(fspan(-89.5,  89.5, 180))  ; CMOR standard
	dest_grd&lon = flt2dble(fspan(0.5, 359.5, 360))    ; CMOR standard
	dest_grd&plev@standard_name = "air_pressure"
	dest_grd&lat@standard_name = "latitude"
	dest_grd&lon@standard_name = "longitude"
	dest_grd&plev@long_name     = "pressure"
	dest_grd&lat@long_name     = "latitude"
	dest_grd&lon@long_name     = "longitude"
	dest_grd&plev@units         = "Pa"
	dest_grd&lat@units         = "degrees_north"
	dest_grd&lon@units         = "degrees_east"
	dest_grd&lat@axis          = "Y"
	dest_grd&lon@axis          = "X"
	dest_grd&plev@axis          = "Z"
    else if (any(dims .eq. "lev")) then
	dest_grd     = new((/50, 180, 360/), float)
	dest_grd!0   = "lev"
	dest_grd!1   = "lat"
	dest_grd!2   = "lon"
	dest_grd&lev = (/5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115, \
			125, 135, 145, 155, 165, 175, 185, 195, 205, \
			216.846755981445, 241.349014282227, \
			280.780731201172, 343.250457763672, \
			427.315551757812, 536.715637207031, \
			665.414123535156, 812.781616210938, \
			969.065124511719, 1130.93493652344,  \
			1289.60461425781, 1455.77014160156, \
			1622.92565917969, 1801.55810546875, \
			1984.85461425781, 2182.90478515625, \
			2388.41748046875, 2610.93505859375, \
			2842.564453125, 3092.20483398438, \
			3351.29467773438, 3628.0576171875, \
			3913.26440429688, 4214.4951171875, \
			4521.91796875, 4842.56591796875, 5166.1298828125, \
			5499.2451171875, 5831.29443359375/)
	dest_grd&lat = flt2dble(fspan(-89.5,  89.5, 180))  ; CMOR standard
	dest_grd&lon = flt2dble(fspan(0.5, 359.5, 360))    ; CMOR standard
	dest_grd&lev@standard_name = "depth"
	dest_grd&lat@standard_name = "latitude"
	dest_grd&lon@standard_name = "longitude"
	dest_grd&lev@long_name     = "ocean depth coordinate"
	dest_grd&lat@long_name     = "latitude"
	dest_grd&lon@long_name     = "longitude"
	dest_grd&lev@units         = "m"
	dest_grd&lat@units         = "degrees_north"
	dest_grd&lon@units         = "degrees_east"
	dest_grd&lat@axis          = "Y"
	dest_grd&lon@axis          = "X"
	dest_grd&lev@axis          = "Z"
    else
	dest_grd     = new((/180, 360/), float)
	dest_grd!0   = "lat"
	dest_grd!1   = "lon"
	dest_grd&lat = flt2dble(fspan(-89.5,  89.5, 180))  ; CMOR standard
	dest_grd&lon = flt2dble(fspan(0.5, 359.5, 360))    ; CMOR standard
	dest_grd&lat@standard_name = "latitude"
	dest_grd&lon@standard_name = "longitude"
	dest_grd&lat@long_name     = "latitude"
	dest_grd&lon@long_name     = "longitude"
	dest_grd&lat@units         = "degrees_north"
	dest_grd&lon@units         = "degrees_east"
	dest_grd&lat@axis          = "Y"
	dest_grd&lon@axis          = "X"
    end if
    end if
  else if (selected_grid.eq."2.5x2.5") then
    if  (any(dims .eq. "plev")) then
	dest_grd     = new((/17, 72, 144/), float)
	dest_grd!0   = "plev"
	dest_grd!1   = "lat"
	dest_grd!2   = "lon"
	dest_grd&plev = (/100000., 92500., 85000., 70000., 60000., \
			 50000., 40000., 30000., 25000., 20000., 15000., \
			 10000., 7000., 5000., 3000., 2000., 1000./)
	dest_grd&lat = flt2dble(fspan(-88.75, 88.75, 72))  ; CMOR standard
	dest_grd&lon = flt2dble(fspan(1.25, 358.75, 144))  ; CMOR standard
	dest_grd&plev@standard_name = "air_pressure"
	dest_grd&lat@standard_name = "latitude"
	dest_grd&lon@standard_name = "longitude"
	dest_grd&plev@long_name     = "pressure"
	dest_grd&lat@long_name     = "latitude"
	dest_grd&lon@long_name     = "longitude"
	dest_grd&plev@units         = "Pa"
	dest_grd&lat@units         = "degrees_north"
	dest_grd&lon@units         = "degrees_east"
	dest_grd&lat@axis          = "Y"
	dest_grd&lon@axis          = "X"
	dest_grd&plev@axis          = "Z"
    else if (any(dims .eq. "lev")) then
	dest_grd     = new((/50, 72, 144/), float)
	dest_grd!0   = "lev"
	dest_grd!1   = "lat"
	dest_grd!2   = "lon"
	dest_grd&lev = (/5, 15, 25, 35, 45, 55, 65, 75, 85, 95, 105, 115, \
			125, 135, 145, 155, 165, 175, 185, 195, 205, \
			216.846755981445, 241.349014282227, \
			280.780731201172, 343.250457763672, \
			427.315551757812, 536.715637207031, \
			665.414123535156, 812.781616210938, \
			969.065124511719, 1130.93493652344, \
			1289.60461425781, 1455.77014160156, \
			1622.92565917969, 1801.55810546875, \
			1984.85461425781, 2182.90478515625, \
			2388.41748046875, 2610.93505859375, \
			2842.564453125, 3092.20483398438, \
			3351.29467773438, 3628.0576171875, \
			3913.26440429688, 4214.4951171875, 4521.91796875, \
			4842.56591796875, 5166.1298828125, \
			5499.2451171875, 5831.29443359375/)
	dest_grd&lat = flt2dble(fspan(-88.75, 88.75, 72))  ; CMOR standard
	dest_grd&lon = flt2dble(fspan(1.25, 358.75, 144))  ; CMOR standard
	dest_grd&lev@standard_name = "depth"
	dest_grd&lat@standard_name = "latitude"
	dest_grd&lon@standard_name = "longitude"
	dest_grd&lev@long_name     = "ocean depth coordinate"
	dest_grd&lat@long_name     = "latitude"
	dest_grd&lon@long_name     = "longitude"
	dest_grd&lev@units         = "m"
	dest_grd&lat@units         = "degrees_north"
	dest_grd&lon@units         = "degrees_east"
	dest_grd&lat@axis          = "Y"
	dest_grd&lon@axis          = "X"
	dest_grd&lev@axis          = "Z"
    else
	dest_grd     = new((/72, 144/), float)
	dest_grd!0   = "lat"
	dest_grd!1   = "lon"
	dest_grd&lat = flt2dble(fspan(-88.75, 88.75, 72))  ; CMOR standard
	dest_grd&lon = flt2dble(fspan(1.25, 358.75, 144))  ; CMOR standard
	dest_grd&lat@standard_name = "latitude"
	dest_grd&lon@standard_name = "longitude"
	dest_grd&lat@long_name     = "latitude"
	dest_grd&lon@long_name     = "longitude"
	dest_grd&lat@units         = "degrees_north"
	dest_grd&lon@units         = "degrees_east"
	dest_grd&lat@axis          = "Y"
	dest_grd&lon@axis          = "X"
    end if
    end if
  else
    dest_grd = find_destination_grid(index, var0, \
				     field_type0, selected_grid)
  end if
  end if
  end if
  ;annots = project_style(diag_script_info, "annots")
  Fill = default_fillvalue(typeof(A0))
  if (rank .eq. 3) then
    var_grd = new((/dim_MOD, dim_seas, dimsizes(dest_grd&lat), \
		   dimsizes(dest_grd&lon)/), typeof(A0), Fill)
    var_grd!0 = "models"
    var_grd!1 = "season"
    var_grd!2 = "lat"
    var_grd!3 = "lon"
    IAV_mmm = new((/dim_seas, dimsizes(dest_grd&lat), \
		   dimsizes(dest_grd&lon)/), typeof(A0), Fill)
  else if (rank .eq. 4) then
    var_grd = new((/dim_MOD, dim_seas, dimsizes(dest_grd&$lev_name$), \
		   dimsizes(dest_grd&lat), dimsizes(dest_grd&lon)/),\
		  typeof(A0), Fill)
    IAV_mmm = new((/dim_seas, dimsizes(dest_grd&$lev_name$), \
		   dimsizes(dest_grd&lat), dimsizes(dest_grd&lon)/),\
		  typeof(A0), Fill)
    var_grd!0 = "models"
    var_grd!1 = "season"
    var_grd!2 = lev_name
    var_grd&$lev_name$ = dest_grd&$lev_name$
    var_grd!3 = "lat"
    var_grd!4 = "lon"
  end if
  end if
  ;var_grd&models = annots
  if (isatt(diag_script_info, "seasons")) then
    var_grd&season = (/ispan(0, dim_seas - 1, 1)/)
  else
    var_grd&season = (/0/)
  end if
  var_grd&lat = dest_grd&lat
  var_grd&lon = dest_grd&lon

  do imod = 0, dim_MOD - 1
    if (imod .ne. 0) then
      A0 = read_data(imod, var0, field_type0)
      dims = getvardims(A0)
      idx = ind(dims .eq. "lat" .or. dims .eq. "rlat" .or. dims .eq. "j")
      if (ismissing(idx)) then
        error_msg("f", diag_script, "", "no lat dimension")
      end if
      lat_name = dims(idx)
      idx = ind(dims .eq. "lon" .or. dims .eq. "rlon" .or. dims .eq. "i")
      if (ismissing(idx)) then
        error_msg("f", diag_script, "", "no lon dimension")
      end if
      lon_name = dims(idx)
    end if
    ; determine if rectilinear or curvilinear original grid
    if (lat_name .ne. "lat") .or. (lon_name .ne. "lon") then
      regular_grid = False
      filename = interface_get_fullpath(var0, field_type0, imod)
      infile = addfile(filename, "r")
      lat2D = infile->lat
      lon2D = infile->lon
      if (selected_grid .eq. "1x1") then
        grid_res = "1degCMOR"
      else if (selected_grid .eq. "2.5x2.5") then
        grid_res = "2.5degCMOR"
      else
        error_msg("f", diag_script, "", "selected grid not possible"\
                  + " if original grid is irregular, select either "\
                  + " 1x1 or 2.5x2.5 as target grid")
      end if
      end if                
    else
      regular_grid = True
    end if
    ; Calculate IAV (function in ./diag_scripts/lib/ncl/statistics.ncl)
    ; cut the first 100 yr
    time = A0&time
    start_yr = cd_calendar(time(0), 0)
    ind_end = dimsizes(time) - 1
    end_yr = cd_calendar(time(ind_end), 0)
    if ((end_yr(0, 0) - start_yr(0, 0) + 1) .lt. 500) then
      print("warning: Control run is less than 500 years.")
    end if
    new_start_yr = toint((start_yr(0, 0) + 100))
    end_yr_int = toint(end_yr(0, 0))
    if (isatt(diag_script_info, "periodlength")) then
      length_of_period = toint(diag_script_info@periodlength)
      nr_periods = toint(floor((dimsizes(time) - 12 * 100) / \
                               (length_of_period * 12.)))
      if (nr_periods .lt. 1) then
        error_msg("w", diag_script, "", "time range too short, " + \
                  "less than 1 period covered, continue with " + \
                  "next model")
	continue
      end if
      yr_possible = (dimsizes(time) - 12 * 100) / (length_of_period * 12.)
      rest = yr_possible - nr_periods
      start_yrs = ispan(new_start_yr, end_yr_int, length_of_period)
      if (rank .eq. 3) then
        data_tmp = new((/nr_periods, dim_seas, dimsizes(A0&$lat_name$), \
                        dimsizes(A0&$lon_name$)/), typeof(A0), Fill)
      else if (rank .eq. 4) then
        data_tmp = new((/nr_periods, dim_seas, dimsizes(A0&$lev_name$), \
			dimsizes(A0&$lat_name$), dimsizes(A0&$lon_name$)/), \
		       typeof(A0), Fill)
      end if
      end if
      ;log_debug(data_tmp)
      if (diag_script_info@time_avg .eq. "seasonalclim") then
        ; remove seasonal cycle
        deseas_data = rmMonAnnCycTLL(A0)
        delete(A0)
        A0 = deseas_data
      end if
      do per = 0, nr_periods - 1
        if ((rest .gt. 0) .and. (per .eq. nr_periods - 1)) then
          new_end_yr = end_yr_int
        else
          new_end_yr = toint(start_yrs(per) + length_of_period - 1.0)
        end if
        ;  calculate seasonal/annual averages over periods
        if (rank .eq. 3) then
		if  ((dim_seas .eq. 1) .and. \
		    (diag_script_info@time_avg .eq. "annualclim")) then
		    data_tmp(per, 0, :, :) = \
			time_operations(A0, start_yrs(per), new_end_yr, \
					"average", \
					diag_script_info@time_avg, True)
		else if ((dim_seas .eq. 1) .and. \
			 (diag_script_info@time_avg .eq. "seasonalclim")) \
		    then
		    tmp = time_operations(A0, start_yrs(per), new_end_yr, \
					  "average",\
					  diag_script_info@time_avg, True)
		    data_tmp(per, 0, :, :) = \
			tmp(diag_script_info@seasons, :, :)
		    delete(tmp)
		else
		    data_tmp(per, :, :, :) = \
			time_operations(A0, start_yrs(per), new_end_yr, \
					"average", \
					diag_script_info@time_avg, True)
		end if
		end if
	else if (rank .eq. 4) then
		if ((dim_seas .eq. 1) .and. \
		    (diag_script_info@time_avg .eq. "annualclim")) then
		    data_tmp(per, 0, :, :, :) = \
			time_operations(A0, start_yrs(per), new_end_yr, \
					"average", \
					diag_script_info@time_avg, True)
		else if ((dim_seas .eq. 1) .and. \
			 (diag_script_info@time_avg .eq. "seasonalclim")) \
		    then
		    tmp = time_operations(A0, start_yrs(per), new_end_yr, \
					  "average",\
					  diag_script_info@time_avg, True)
		    data_tmp(per, 0, :, :, :) = \
			tmp(diag_script_info@seasons, :, :, :)
		    delete(tmp)
		else
		    data_tmp(per, :, :, :, :) = \
			time_operations(A0, start_yrs(per), new_end_yr,\
					"average", \
					diag_script_info@time_avg, True)
		end if
		end if
	end if
        end if
      end do
      data_dtr = dtrend_quadratic_msg_n(data_tmp, False, False, 0)
      data1 = dim_stddev_n(data_dtr, 0)
      if (rank .eq. 3) then
        data1!1 = lat_name
        data1!2 = lon_name
      else if (rank .eq. 4) then
        data1!1 = lev_name
        data1&$lev_name$ = data_tmp&$lev_name$
        data1!2 = lat_name
        data1!3 = lon_name
      end if
      end if
      data1&$lat_name$ = data_tmp&$lat_name$
      data1&$lon_name$ = data_tmp&$lon_name$
      delete(start_yrs)
      delete(data_tmp)
      delete(data_dtr)
      if (diag_script_info@time_avg .eq. "seasonalclim") then
        delete(deseas_data)
      end if
    else
      data1 = interannual_variability(A0, new_start_yr, end_yr_int, \
                                      diag_script_info@time_avg) ;, "quadratic")
    end if
    if (rank .eq. 3) then
	if (dim_seas .eq. 1) then
	    ; Regrid to common grid
	    ; (function in ./diag_scripts/lib/ncl/statistics.ncl)
	    if (regular_grid .eq. True) then
		tmp = linint2(data1&$lon_name$, data1&$lat_name$, data1,\
			      True, dest_grd&lon, dest_grd&lat, 0)
	    else
		tmp =  regrid_3D_to_rectilinear_grid(data1, lon2D, lat2D,\
						     grid_res, filename,\
						     False)
	    end if
	    var_grd(imod, 0, :, :) = (/rm_single_dims(tmp)/)
	    delete(tmp)
	else
	    ; Regrid to common grid
	    if (regular_grid .eq. True) then
		; (function in ./diag_scripts/lib/ncl/statistics.ncl)
		var_grd(imod, :, :, :) = (/linint2(data1&$lon_name$, \
						   data1&$lat_name$, \
						   data1, \
						   True, dest_grd&lon,\
						   dest_grd&lat, 0)/)
	    else
		; (function in ./diag_scripts/lib/ncl/regridding.ncl)
		var_grd(imod, :, :, :) = \
		    (/regrid_3D_to_rectilinear_grid(data1, lon2D, lat2D, \
						    grid_res, filename, \
						    False)/)
	    end if
	end if
    else if (rank .eq. 4) then
	dz = dimsizes(data1&$lev_name$)
	if (dim_seas .eq. 1) then
	    ; Regrid to common horizontal grid
	    if (regular_grid .eq. True) then
		; (function in ./diag_scripts/lib/ncl/statistics.ncl)
		regrd_tmp = rm_single_dims(linint2_Wrap(data1&$lon_name$, \
							data1&$lat_name$, \
							data1, True, \
							dest_grd&lon, \
							dest_grd&lat, 0))
	    else
		regrd_tmp = new((/dz, dimsizes(dest_grd&lat), \
				 dimsizes(dest_grd&lon)/), float)
		plev_tmp = new((/1, dimsizes(data1&$lat_name$), \
				dimsizes(data1&$lon_name$)/), float)
		do z = 0, dz - 1
		    plev_tmp(0, :, :) = data1(0, z, :, :)
		    regrd_tmp(z, :, :) = \
			(/regrid_3D_to_rectilinear_grid(plev_tmp, lon2D, \
							lat2D, grid_res, \
							filename, False)/)
		end do
	    regrd_tmp!0 = lev_name
	    regrd_tmp&$lev_name$ = data1&$lev_name$
	    end if
	    ; interpolate vertical levels to common levels
	    if (lev_name .eq. "plev") then
		var_grd(imod, 0, :, :, :) = \
		    (/rm_single_dims(int2p_n(regrd_tmp&$lev_name$, \
					     regrd_tmp, \
					     dest_grd&$lev_name$, 1, 0))/)
	    else
		var_grd(imod, 0, :, :, :) = \
		    (/rm_single_dims(linint1_n(regrd_tmp&$lev_name$, \
					       regrd_tmp, False, \
					       dest_grd&$lev_name$, 0, \
					       0))/)
	    end if
	else
	    ; Regrid to common grid
	    if (regular_grid .eq. True) then
	    ; (function in ./diag_scripts/lib/ncl/statistics.ncl)
		regrd_tmp = linint2(data1&$lon_name$, data1&$lat_name$, \
				    data1, True, dest_grd&lon, \
				    dest_grd&lat, 0)
	    else
		regrd_tmp = new((/dim_seas, dz, dimsizes(dest_grd&lat), \
				 dimsizes(dest_grd&lon)/), float)
		do z = 0, dz - 1
		    plev_tmp = data1(:, z, :, :)
		    regrd_tmp(:, z, :, :) = \
			regrid_3D_to_rectilinear_grid(plev_tmp, lon2D, \
						      lat2D, grid_res, \
						      filename, False)
		end do
	    end if
	    regrd_tmp!1 = lev_name
	    regrd_tmp&$lev_name$ = data1&$lev_name$
	    ; interpolate vertical levels to common levels
	    var_grd(imod, :, :, :, :) = (/int2p_n(regrd_tmp&$lev_name$, \
						  regrd_tmp, \
						  dest_grd&$lev_name$, 1, \
						  dz)/)
	end if
    end if
    end if
    if (rank .eq. 4) then
	delete(regrd_tmp)
	if (regular_grid .eq. False) then
	    delete(plev_tmp)
	end if
    end if
    delete(A0)
    if (regular_grid .eq. False) then
	delete(lat2D)
	delete(lon2D)
    end if
    delete(start_yr)
    delete(end_yr)
    delete(data1)
    delete(time)
  end do  ; imod
  ; ###########################################
  ; # Calculate multi-model mean of IAV       #
  ; ###########################################
  do s = 0, dim_seas - 1
    if (rank .eq. 3) then
      IAV_mmm(s, :, :) = \
      rm_single_dims(dim_avg_n(var_grd(:, s, :, :), 0) * sqrt(2.))
    else if (rank .eq. 4) then
      IAV_mmm(s, :, :, :) = \
      rm_single_dims(dim_avg_n(var_grd(:, s, :, :, :), 0) * sqrt(2.))
    end if
    end if
  end do
  if (rank .eq. 3) then
    IAV_mmm!0 = "season"
    IAV_mmm&season = diag_script_info@seasons
    IAV_mmm!1 = "lat"
    IAV_mmm!2 = "lon"
  else if (rank .eq. 4) then
    IAV_mmm!0 = "season"
    IAV_mmm&season = diag_script_info@seasons
    IAV_mmm!1 = lev_name
    IAV_mmm&$lev_name$ = var_grd&$lev_name$
    IAV_mmm!2 = "lat"
    IAV_mmm!3 = "lon"
  end if
  end if
  IAV_mmm&lat = dest_grd&lat
  IAV_mmm&lon = dest_grd&lon
  ; ###########################################
  ; # Output to netCDF                        #
  ; ###########################################
  IAV_mmm@diag_script = (/diag_script/)
  IAV_mmm@var = "iav"
  if (isatt(variable_info, "long_name")) then
    IAV_mmm@var_long_name = "iav of " + variable_info@long_name
  end if
  if (isatt(variable_info, "units")) then
    IAV_mmm@var_units     = variable_info@units
  end if
  IAV_mmm@comment = "models included in mmm are " + model_info@model
  ; Function in ~/interface_scripts/auxiliary.ncl
  ncdf_outfile = ncdf_write(IAV_mmm, work_dir + outfile)
end
