; #############################################################################
; Annual mean zonal mean multi model changes
; Author: Ruth Lorenz (ETH Zurich, Switzerland)
; CRESCENDO project
; #############################################################################
;
; Description
;    Multi-model annual or seasomal mean zonal mean changes
;    e.g. 2081-2100 relative to base period (1986-2005 in CMIP5)
;    for multiple RCPs
;
; Required diag_script_info attributes (diagnostics specific)
;    ncdf: filename for netcdf to write
;    scenarios: which scenarios to include
;    periods: list with start years of periods to be included
;    	       e.g. (/"2046", "2081", "2180"/)
;    time_avg: list with seasons or annual to be included e.g (/"annualclim"/)
;
; Optional diag_script_info attributes (diagnostic specific)
;    grid: grid for regridding (coarsest, finest, ref, 1x1, 2.5x2.5)
;    title: title of whole figure
;    label: labels of individual panels
;    max_hori: maximum number of panels in horizontal
;    max_vert: maximum number of panels in vertical
;    colormap: colortable different from default
;    explicit_cn_levels: range and levels for contour plot
;    base_cn: flag indicating if base period plotted in
;             contour lines (as in ipcc Fig. 12.19 for wind)
;    base_cnLevelSpacing: spacing for contour levels (e.g. 10)
;                         if base_cn is True
;    base_cnMinLevel: min value for base contours
;    base_cnMaxLevel: max level for base contours
;
; Required variable_info attributes (variable specific)
;    long_name: long name of variable
;    units: variable units
;
; Optional variable_info attributes (variable specific)
;    att1: short description
;    att2: short description
;
; Caveats
;    slow for irregular ocean grids (regridding on each depth level),
;    regridding only works for models with all relevant information provided
;    (therefore model name printed in loop, in case of crash know which model
;    caused crash)
;
; Modification history
;    20170517-A_lore_ru: changed to regrid to common grid first before
;                        calculating zonal average and added ability to use
;                        irregular grids (e.g. thetao)
;    20170127-A_lore_ru: changed routine for plotting stippling/hatching to
;                        zonalmean_profile to get dots and hatch
;    20161031-A_lore_ru: increased flexibility and add stippling and hatching
;    20160831-A_lore_ru: adapted for ESMValTool
;    20130503-A_sedl_ja: written for IPCC AR5
;
; #############################################################################
; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
load "./interface_scripts/interface.ncl"

load "./diag_scripts/shared/plot/style.ncl"
load "./diag_scripts/shared/statistics.ncl"
load "./diag_scripts/shared/ensemble.ncl"

begin
  ; ##############################################################
  ; # Fetch general parameters, set in namelist_collins13ipcc.xml#
  ; # passed via environment variables by python code            #
  ; ##############################################################
  enter_msg(diag_script, "")

  ; 'diag_script' and 'variables' are fetched from the
  ; above 'ncl.interface' file
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (var: " + variables(0) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; 'models', 'variables' and 'field_types' are fetched from the
  ; above 'ncl.interface' file
  dim_MOD = dimsizes(dataset_info@dataset)
  dim_VAR = dimsizes(variables)
  var0 = variables(0)
  field_type0 = field_types(0)

  ; ###########################################
  ; # References Information                  #
  ; ###########################################
;  write_references(diag_script, \  ; script name
;                   "A_lore_ru", \  ; authors
;                   "", \  ; contributors
;                   "D_collins13ipcc",    \  ; diag_references
;                   "",          \  ; obs_references
;                   "P_crescendo")   ; proj_references
  ; Check required diag_script_info attributes
  req_atts = (/"scenarios", "periods", "time_avg"/)
  exit_if_missing_atts(diag_script_info, req_atts)
end

begin
  ; Get environment variables
  if (isatt(diag_script_info, "plot_type")) then
    file_type = diag_script_info@plot_type
  else
    file_type = "pdf"
  end if

  ; Basename of diag_script
  diag_script_base = basename(diag_script)

  ; Output netcdf directory
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)
end

begin
  ; #############################################
  ; # Get parameters from ./variable_defs/*.ncl #
  ; # passed via the 'info' attribute           #
  ; #############################################
  if(isvar("MyParam")) then
      delete(MyParam)
  end if
  if(isatt(variable_info, "long_name")) then
    MyParam = variable_info@long_name
    log_info(" MyParam = " + MyParam)
  else
    error_msg("f", diag_script, "", "missing variable attribute " + \
              "'long_name'")
  end if

  ; ###########################################
  ; # Get data and average time               #
  ; ###########################################
  ; get data from first model
  imod = 0  ; NCL array indicies start from zero
  log_info("processing " + datset_info@dataset(imod) + "_"\
           + datset_info@exp(imod) + "_"\
           + datset_info@ensemble(imod))

  ; See ./interface_scripts/data_handling.ncl
  A0 = read_data(imod, var0, field__type0)

  ; Check dimensions
  dims = getvardims(A0)
  if (dimsizes(dims) .lt. 2) then
    error_msg("f", diag_script, "", dimsizes(dims) + " dimensions, " + \
              "need 2 or 3")
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

  idx = ind(dims .eq. "lev" .or. dims .eq. "plev")
  if (ismissing(idx)) then
    error_msg("f", diag_script, "", "no level dimension")
  end if
  lev_name = dims(idx)

    ; determine destination grid (regridding.ncl)
    index = ispan(0, dim_MOD - 1, 1)
    if (isatt(diag_script_info, "seasons")) then
        dim_seas = dimsizes(diag_script_info@seasons)
    else
        dim_seas = 1
    end if
    if (.not.(isatt(diag_script_info, "grid"))) then
        error_msg("w", diag_script, "", \
                  "Regridding method not defined. " + "Using 2.5x2.5")
        selected_grid = "2.5x2.5"
    else
        selected_grid = diag_script_info@grid
    end if
    log_info("Setting target grid: " + selected_grid)
    if (selected_grid.eq."ref") then
        ; set reference model
        ref_ind = ind(datset_info@dataset .eq. diag_script_info@ref_model(0))
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
            dest_grd!0   = lev_name
            dest_grd!1   = "lat"
            dest_grd!2   = "lon"
            dest_grd&plev = (/100000., 92500., 85000., 70000., 60000., \
                             50000., 40000., 30000., 25000., 20000., \
                             15000., 10000., 7000., 5000., 3000., 2000., \
                             1000./)
            dest_grd&lat = flt2dble(fspan(-89.5,  89.5, 180))  ; CMOR standard
            dest_grd&lon = flt2dble(fspan(0.5, 359.5, 360))    ; CMOR standard
            dest_grd&plev@standard_name = "air_pressure"
            dest_grd&plev@long_name     = "pressure"
            dest_grd&plev@units         = "Pa"
            dest_grd&plev@axis          = "Z"
        else if (any(dims .eq. "lev")) then
            dest_grd     = new((/50, 180, 360/), float)
            dest_grd!0   = lev_name
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
            dest_grd&lat = flt2dble(fspan(-89.5,  89.5, 180))  ; CMOR standard
            dest_grd&lon = flt2dble(fspan(0.5, 359.5, 360))    ; CMOR standard
            dest_grd&lev@standard_name = "depth"
            dest_grd&lev@long_name     = "ocean depth coordinate"
            dest_grd&lev@units         = "m"
            dest_grd&lev@axis          = "Z"
        else
            dest_grd     = new((/180, 360/), float)
            dest_grd!0   = "lat"
            dest_grd!1   = "lon"
            dest_grd&lat = flt2dble(fspan(-89.5,  89.5, 180))  ; CMOR standard
            dest_grd&lon = flt2dble(fspan(0.5, 359.5, 360))    ; CMOR standard
        end if
        end if
        dest_grd&lat@standard_name = "latitude"
        dest_grd&lon@standard_name = "longitude"
        dest_grd&lat@long_name     = "latitude"
        dest_grd&lon@long_name     = "longitude"
        dest_grd&lat@units         = "degrees_north"
        dest_grd&lon@units         = "degrees_east"
        dest_grd&lat@axis          = "Y"
        dest_grd&lon@axis          = "X"
    else if (selected_grid.eq."2.5x2.5") then
        if (any(dims .eq. "plev")) then
            dest_grd     = new((/17, 72, 144/), float)
            dest_grd!0   = lev_name
            dest_grd!1   = "lat"
            dest_grd!2   = "lon"
            dest_grd&plev = (/100000., 92500., 85000., 70000., 60000., \
                             50000., 40000., 30000., 25000., 20000., 15000., \
                             10000., 7000., 5000., 3000., 2000., 1000./)
            dest_grd&lat = flt2dble(fspan(-88.75, 88.75, 72))  ; CMOR standard
            dest_grd&lon = flt2dble(fspan(1.25, 358.75, 144))  ; CMOR standard
            dest_grd&plev@standard_name = "air_pressure"
            dest_grd&plev@long_name     = "pressure"
            dest_grd&plev@units         = "Pa"
            dest_grd&plev@axis          = "Z"
        else if (any(dims .eq. "lev")) then
            dest_grd     = new((/50, 72, 144/), float)
            dest_grd!0   = lev_name
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
                            3913.26440429688, 4214.4951171875, \
                            4521.91796875, 4842.56591796875, \
                            5166.1298828125, 5499.2451171875, \
                            5831.29443359375/)
            dest_grd&lat = flt2dble(fspan(-88.75, 88.75, 72))  ; CMOR standard
            dest_grd&lon = flt2dble(fspan(1.25, 358.75, 144))  ; CMOR standard
            dest_grd&lev@standard_name = "depth"
            dest_grd&lev@long_name     = "ocean depth coordinate"
            dest_grd&lev@units         = "m"
            dest_grd&lev@axis          = "Z"
        else
            dest_grd     = new((/72, 144/), float)
            dest_grd!0   = "lat"
            dest_grd!1   = "lon"
            dest_grd&lat = flt2dble(fspan(-88.75, 88.75, 72))  ; CMOR standard
            dest_grd&lon = flt2dble(fspan(1.25, 358.75, 144))  ; CMOR standard
        end if
        end if
        dest_grd&lat@standard_name = "latitude"
        dest_grd&lon@standard_name = "longitude"
        dest_grd&lat@long_name     = "latitude"
        dest_grd&lon@long_name     = "longitude"
        dest_grd&lat@units         = "degrees_north"
        dest_grd&lon@units         = "degrees_east"
        dest_grd&lat@axis          = "Y"
        dest_grd&lon@axis          = "X"
    else
        dest_grd = find_destination_grid(index, var0, \
                                         field_type0, selected_grid)
    end if
    end if
    end if
;    annots = project_style(diag_script_info, "annots")
  Fill = default_fillvalue(typeof(A0))

    var_grd = new((/dim_MOD, dim_seas, dimsizes(dest_grd&$lev_name$), \
                   dimsizes(dest_grd&lat), dimsizes(dest_grd&lon)/),\
                  typeof(A0), Fill)
    var_grd!0 = "models"
    var_grd!1 = "season"
    var_grd!2 = lev_name
    var_grd&$lev_name$ = dest_grd&$lev_name$
    var_grd!3 = "lat"
    var_grd!4 = "lon"

    var_reg = new((/dim_MOD, dim_seas, dimsizes(dest_grd&$lev_name$), \
                   dimsizes(dest_grd&lat)/), typeof(A0), Fill)

    var_reg!0 = "models"
    var_reg!1 = "season"
    var_reg!3 = "lat"
    var_reg&models = annots
    if (isatt(diag_script_info, "seasons")) then
        var_reg&season = (/ispan(0, dim_seas - 1, 1)/)
    else
        var_reg&season = (/0/)
    end if
    var_reg!2 = lev_name
    var_reg&$lev_name$ = dest_grd&$lev_name$
    var_reg&lat = dest_grd&lat
    if (isatt(diag_script_info, "iavmode")) then
        if (diag_script_info@iavmode .eq. "each") then
            iav = new((/dim_MOD, dim_seas, dimsizes(dest_grd&lat), \
                       dimsizes(dest_grd&lon)/), typeof(A0), Fill)
            copy_VarMeta(var_reg, iav)
        end if
    end if
  ; Average over time (function in ./diag_scripts/lib/ncl/statistics.ncl)
  data1 = time_operations(A0, -1, -1, "average", diag_script_info@time_avg, \
                          True)

  do imod = 0, dim_MOD - 1
    print(datset_info@dataset(imod))
    if (imod .ne. 0) then
      A0 = read_data(imod, var0, field_type0)
      data1 = time_operations(A0, -1, -1, "average", "annualclim", True)
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
	    error_msg("f", diag_script, "", "selected grid not possible "\
		      + "if original grid is irregular, select either "\
		      + "1x1 or 2.5x2.5 as target grid")
	end if
	end if
    else
	regular_grid = True
    end if
    dz = dimsizes(data1&$lev_name$)
    if (dim_seas .eq. 1) then
	; Regrid to common horizontal grid
	if (regular_grid .eq. True) then
	    ; (function in ./diag_scripts/lib/ncl/statistics.ncl)
	    regrd_tmp = linint2_Wrap(data1&$lon_name$, data1&$lat_name$, \
				     data1, True, dest_grd&lon, \
				     dest_grd&lat, 0)
	else
	    regrd_tmp = new((/dz, dimsizes(dest_grd&lat), \
			     dimsizes(dest_grd&lon)/), float)
	    plev_tmp = new((/1, dimsizes(data1&$lat_name$), \
			    dimsizes(data1&$lon_name$)/), float)
	    do z = 0, dz - 1
		plev_tmp(0, :, :) = data1(z, :, :)
		regrd_tmp(z, :, :) = \
		    (/regrid_3D_to_rectilinear_grid(plev_tmp, lon2D, \
						    lat2D, grid_res,\
						    filename, False)/)
	    end do
	regrd_tmp!0 = lev_name
	regrd_tmp&$lev_name$ = data1&$lev_name$
	end if
	; interpolate vertical levels to common levels
	if (lev_name .eq. "plev") then
	    var_grd(imod, 0, :, :, :) = \
		(/rm_single_dims(int2p_n(regrd_tmp&$lev_name$, \
					 regrd_tmp, dest_grd&$lev_name$,\
					 1, 0))/)
	else
	    var_grd(imod, 0, :, :, :) = \
		(/rm_single_dims(linint1_n(regrd_tmp&$lev_name$, \
					   regrd_tmp, False, \
					   dest_grd&$lev_name$, 0, 0))/)
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
		    regrid_3D_to_rectilinear_grid(plev_tmp, lon2D, lat2D, \
						  grid_res, filename, \
						  False)
	    end do
	end if
	regrd_tmp!1 = lev_name
	regrd_tmp&$lev_name$ = data1&$lev_name$
	; interpolate vertical levels to common levels
	var_grd(imod, :, :, :, :) = (/int2p_n(regrd_tmp&$lev_name$, \
					      regrd_tmp, \
					      dest_grd&$lev_name$, 1, dz)/)
    end if
    delete(regrd_tmp)
    if (regular_grid .eq. False) then
      delete(plev_tmp)
    end if
    delete([/A0, data1, dz/])
    if (regular_grid .eq. False) then
      delete([/lat2D, lon2D/])
    end if
    ; Average over longitudes for zonal mean
    var_reg = dim_avg_n_Wrap(var_grd, 4)

    dnames = getvardims(var_reg)

    if (isatt(diag_script_info, "iavmode")) then
      if (diag_script_info@iavmode .eq. "each") then
        ; #########################################################
        ; # Read natural variability for stippling for each model #
        ; #########################################################
        iav_file = addfile(work_dir + "/" + \
                           "IAV_for_stippandhatch_ipcc12/IAV_piControl_" + \
                           datset_info@dataset(imod) + "_" + \
                           datset_info@ensemble(imod) + "_" + field_type0 \
                           + "_" + var0 + "_" + \
                           diag_script_info@time_avg + ".nc", "r")
        iav_native = iav_file->iav
        ; Average over longitudes for zonal mean
        iav_zon = dim_avg_n_Wrap(iav_native, 2)
        dnames_iav = getvardims(iav_zon)
        tmp_regrd = linint1_n_Wrap(iav_zon&lat, iav_zon, False, \
                                   dest_grd&lat, 0, dx)
        iav_vreg = int2p_n_Wrap(tmp_regrd&lev_name, tmp_regrd, \
                                dest_grd&lev_name, 1, dz)
        if (dim_seas .eq. 1) then
          iav(imod, 0, :, :) = iav_vreg
        else
          iav(imod, :, :, :) = iav_vreg
        end if
        delete(iav_native)
      end if
    end if
  end do

  if (isatt(diag_script_info, "iavmode")) then
    if (diag_script_info@iavmode .eq. "each") then
      ; calculate multi-model mean of iav
      tmp_pow = dim_avg_n(iav ^ 2, 0)
      natvar_mmm = rm_single_dims(tmp_pow * sqrt(2.))
    else
      log_info("Error: This iavmode is not know")
    end if
  else
    ; read already calculated mmm iav
    natvar_file = addfile(work_dir + "/" + "IAV_mmm_for_stippandhatch" + \
                          "_ipcc12/IAV_mmm_piControl_" + field_type0 + \
                          "_" + var0 + "_" + diag_script_info@time_avg + \
                          ".nc", "r")
    natvar = natvar_file->iav
    ; Average over longitudes for zonal mean
    natvar_mmm = dim_avg_n_Wrap(natvar, 3)
  end if

  ; ##############################################
  ; # Calculate change from hist to periods      #
  ; ##############################################
  ; loop over rcps and periods
  dim_scen = dimsizes(diag_script_info@scenarios)
  dim_per = dimsizes(diag_script_info@periods)

  to_plot = new((/dim_scen * dim_per * dim_seas, \
                 dimsizes(dest_grd&$lev_name$), dimsizes(dest_grd&lat)/), \
                typeof(dest_grd), Fill)
  if (isatt(diag_script_info, "base_cn")) then
    to_plot_base = new((/dim_scen * dim_per * dim_seas, \
                       dimsizes(dest_grd&$lev_name$), \
                       dimsizes(dest_grd&lat)/), typeof(dest_grd), Fill)
  end if
  to_plot_signif = new((/dim_scen * dim_per * dim_seas, \
                        dimsizes(dest_grd&$lev_name$),  \
                        dimsizes(dest_grd&lat)/), \
                       typeof(dest_grd), Fill)
  to_plot_not_signif = new((/dim_scen * dim_per * dim_seas, \
                            dimsizes(dest_grd&$lev_name$),  \
                            dimsizes(dest_grd&lat)/), \
                           typeof(dest_grd), Fill)

  model_number = new((/dim_scen * dim_per * dim_seas/), integer, Fill)
  i = 0
  do seas = 0, dim_seas - 1  ; loop over seasons
    do rcp = 0, dim_scen - 1  ; loop over rcps
      do per = 0, dim_per - 1  ; loop over periods
        idx_rcp = ind((datset_info@exp .eq. \
                       diag_script_info@scenarios(rcp)) .and. \
                      (datset_info@start_year .eq. \
                       diag_script_info@periods(per)))
        proj = var_reg(idx_rcp, seas, :, :)
        proj_avg = dim_avg_n_Wrap(proj, 0)
        ; find historical runs from same models as in rcp
        idx_hist = find_model_index(datset_info@dataset(idx_rcp), \
                                    datset_info@dataset, "historical",\
                                    datset_info@exp, \
                                    datset_info@ensemble(idx_rcp), \
                                    datset_info@ensemble, 0)
        base = var_reg(idx_hist, seas, :, :)
        base_avg = dim_avg_n_Wrap(base, 0)
        var_diff = proj - base
        delete([/base, proj/])
        to_plot(i, :, :) = proj_avg - base_avg
        if (isatt(diag_script_info, "base_cn")) then
          to_plot_base(i, :, :) = base_avg
        end if
        ; determine sigma and sign for significance
        var_pos_signif = dim_num_n(where(var_diff .gt. 0., 1., \
                                         to_plot@_FillValue) .eq. 1., 0)
        var_neg_signif = dim_num_n(where(var_diff .lt. 0., 1., \
                                         to_plot@_FillValue) .eq. 1., 0)
        var_pos_signif_tmp = \
          var_pos_signif / where((var_pos_signif + var_neg_signif) \
                                 .ne. 0., var_pos_signif + var_neg_signif, \
                                 to_plot@_FillValue)
        var_neg_signif_tmp = \
          var_neg_signif / where((var_pos_signif + var_neg_signif) \
                                 .ne. 0., var_pos_signif + var_neg_signif, \
                                 to_plot@_FillValue)
        signif_and_sign = where(var_pos_signif_tmp .ge. 0.9 .or. \
                                var_neg_signif_tmp .ge. 0.9, 1., 0.)
        signif_and_sign = where(ismissing(signif_and_sign), 0., \
                                signif_and_sign)
        if (dim_seas .ne. 1) then
          sigma_and_sign = where(abs(to_plot(i, :, :)) .gt. \
                                 abs(2. * natvar_mmm(seas, :, :)), 1., 0.)
          to_plot_not_signif(i, :, :) = where(abs(to_plot(i, :, :)) .lt. \
                                              abs(natvar_mmm(seas, :, :)), \
                                              1., 0.)
        else
          sigma_and_sign = where(abs(to_plot(i, :, :)) .gt. \
                                 abs(2. * natvar_mmm(seas, :, :)), 1., 0.)
          to_plot_not_signif(i, :, :) = where(abs(to_plot(i, :, :)) .lt. \
                                              abs(natvar_mmm(seas, :, :)), \
                                              1., 0.)
        end if
        to_plot_signif(i, :, :) = where(signif_and_sign .eq. 1, \
                                        sigma_and_sign, 0.)
        model_number(i) = dimsizes(idx_rcp)
        i = i + 1
        delete([/idx_rcp, idx_hist, var_diff/])
      end do  ; per
    end do  ; rcp
  end do  ; seas
  to_plot!0 = "panel"
  to_plot&panel = diag_script_info@label

  to_plot!1 = lev_name
  to_plot&$lev_name$ = dest_grd&$lev_name$
  if (isatt(diag_script_info, "base_cn")) then
    to_plot_base!1 = lev_name
    to_plot_base&$lev_name$ = dest_grd&$lev_name$
  end if

  to_plot!2 = "lat"
  to_plot&lat = dest_grd&lat
  if (isatt(diag_script_info, "base_cn")) then
    copy_VarMeta(to_plot, to_plot_base)
  end if
  copy_VarMeta(to_plot, to_plot_signif)
  copy_VarMeta(to_plot, to_plot_not_signif)
  ; ###########################################
  ; # Other Metadata: diag_script, var        #
  ; ###########################################
  ; Add to to_plot, as attributes without prefix
  if(isatt(to_plot, "diag_script")) then  ; Add to existing entries
    temp = to_plot@diag_script
    delete(to_plot@diag_script)
    to_plot@diag_script = array_append_record(temp, (/diag_script/), 0)
    delete(temp)
  else  ; Add as new attribute
    to_plot@diag_script = (/diag_script/)
  end if
  to_plot@var = var0  ; Overwrite existing entry
  if(isatt(variable_info, "long_name")) then
    to_plot@var_long_name = "change in " + variable_info@long_name
  end if
  if(isatt(variable_info, "units")) then
    to_plot@var_units     = variable_info@units
  end if

  if (isatt(diag_script_info, "base_cn")) then
    if(isatt(to_plot_base, "diag_script")) then  ; Add to existing entries
      temp = to_plot_base@diag_script
      delete(to_plot_base@diag_script)
      to_plot_base@diag_script = array_append_record(temp, (/diag_script/), 0)
      delete(temp)
    else
      to_plot_base@diag_script = (/diag_script/)
    end if
    to_plot_base@var = var0 + "_base"
    to_plot_base@var_long_name = variable_info@long_name
    to_plot_base@var_units     = variable_info@units
  end if
  if(isatt(to_plot_signif, "diag_script")) then  ; Add to existing entries
    temp = to_plot_signif@diag_script
    delete(to_plot_signif@diag_script)
    to_plot_signif@diag_script = array_append_record(temp, (/diag_script/), 0)
    delete(temp)
  else  ; Add as new attribute
    to_plot_signif@diag_script = (/diag_script/)
  end if
  to_plot_signif@var = "signif"
  to_plot_signif@var_long_name = "significance: mean change larger than" +\
                                 " 2*natvar and 90% of models have the same"\
                                 + " sign"
  to_plot_signif@var_units     = 1

  if(isatt(to_plot_not_signif, "diag_script")) then  ; Add to existing entries
    temp = to_plot_not_signif@diag_script
    delete(to_plot_not_signif@diag_script)
    to_plot_not_signif@diag_script = array_append_record(temp, \
                                                         (/diag_script/), 0)
    delete(temp)
  else  ; Add as new attribute
    to_plot_not_signif@diag_script = (/diag_script/)
  end if
  to_plot_not_signif@var = "not_signif"
  to_plot_not_signif@var_long_name = "not significant: mean change smaller"\
                                     + " than 1*natvar"
  to_plot_not_signif@var_units     = 1

  model_number!0 = "panel"
  model_number&panel = diag_script_info@label
  if(isatt(model_number, "diag_script")) then  ; Add to existing entries
    temp = model_number@diag_script
    delete(model_number@diag_script)
    model_number@diag_script = array_append_record(temp, (/diag_script/), 0)
    delete(temp)
  else  ; Add as new attribute
    model_number@diag_script = (/diag_script/)
  end if
  model_number@var = "model_nr"
  model_number@var_long_name = "number of models in mmm"
  model_number@var_units     = 1

  ; ###########################################
  ; # Output to netCDF                        #
  ; ###########################################
  ; Optional output (controlled by diag_script_info)
  if(isatt(diag_script_info, "ncdf")) then
    to_plot@ncdf = variable_info@diagnostic + ".nc"
    ncdf_file = work_dir + "/" + to_plot@ncdf
    ncdf_file@existing = "overwrite"
    ; Function in ~/interface_scripts/auxiliary.ncl
    ncdf_outfile = ncdf_write(to_plot, ncdf_file)
    ncdf_file@existing = "append"
    ncdf_outfile = ncdf_write(to_plot_signif, ncdf_file)
    ncdf_outfile = ncdf_write(to_plot_not_signif, ncdf_file)
    ncdf_outfile = ncdf_write(model_number, ncdf_file)
    if ((isatt(diag_script_info, "base_cn")) .and. \
        (diag_script_info@base_cn .eq. True)) then
      ncdf_outfile = ncdf_write(to_plot_base, ncdf_file)
    end if
  end if

end
