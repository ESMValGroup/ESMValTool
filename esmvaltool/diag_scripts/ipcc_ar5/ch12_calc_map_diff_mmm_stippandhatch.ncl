; #############################################################################
; Multi-model mean change map with significance
; Author: Ruth Lorenz (ETH, Switzerland)
; CRESCENDO project
; #############################################################################
;
; Description
;    Annual mean or seasonal change in multi-model mean with signficance
;    in different RCPs in different periods (e.g. 2046-65, 2081-2100 and
;    2181-2200 with respect to 1986-2005)
;
; Required diag_script_info attributes (diagnostics specific)
;    scenarios: list with scenarios to be included in the
;               figure, e.g (/"rcp26","rcp45","rcp60","rcp85"/)
;    periods: list with start years of periods to be included
;    	       e.g. (/"2046","2081","2180"/)
;    time_avg: seasonal or annual means, either (/"annualclim"/) or
;              (/"seasonalclim"/)
;    label: labels on top of each panel, e.g. scenario + rcp,
;           loop to plot is seasons, scenarios, periods, needed to write
;           metadata into netcdf
;
; Optional diag_script_info attributes (diagnostic specific)
;    seasons: list with seasons index to be included if averagetime is
;             "seasonalclim" (required for seasonalclim),
;             DJF:0, MAM:1, JJA:2, SON:3
;    iavmode: use IAV data from multi-model mean calculated by
;             IAV_mmm_for_stippandhatch_ipcc12.ncl or include identical models
;             to the ones plotted calculated for each model on native grid
;             with IAV_for_stippandhatch_ipcc12.ncl default if not given is
;             mmm, otherwise set iavmode = "each"
;    grid:    grid for regridding (coarsest, finest, ref, 1x1, 2.5x2.5)
;    percent: 0 or 1, difference in percent = 1, default = 0
;;
; Required variable_info attributes (variable specific)
;    long_name: long name of variable
;    units: variable units
;
; Optional variable_info attributes (variable specific)
;    att1: short description
;    att2: short description
;
; Caveats
;    Only handles one variable at a time. Appends multiple output variables
;    (variable to_plot, signficance to_plot_signif and not significant
;    to_plot_not_signif and number of model in mmm to model_nr)
;    to the same file, if the file already exists overwrites
;;
; Modification history
;    20170120-A_lore_ru: separated plotting from calculations
;    20161027-A_lore_ru: added stippling and hatching for significance
;    20161024-A_lore_ru: adjusted to plot multiple figs with one diag script
;    20160621-A_lore_ru: adapted to ESMValTool
;    20130501-A_sedl_ja: written for IPCC AR5.
;;
; #############################################################################

; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
load "./interface_scripts/interface.ncl"

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

  ; 'datasets', 'variables' and 'field_types' are fetched from the
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
  ; Basename of diag_script
  diag_script_base = basename(diag_script)

  ; Output netcdf directory
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)  ; Needed?

  input_dir = diag_script_info@input_files
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
  ; get data from first dataset
  imod = 0  ; NCL array indicies start from zero
  log_debug("processing " + dataset_info@dataset(imod) + "_"\
            + dataset_info@exp(imod) + "_" \
            + dataset_info@ensemble(imod))

  ; See ./interface_scripts/data_handling.ncl
  A0 = read_data(imod, var0, field_type0)

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
  log_info("Setting target grid:")
    if (selected_grid.eq."ref") then
        ; set reference model
        ref_ind = ind(datasets@name.eq.diag_script_info@ref_dataset(0))
        ref_dataset = diag_script_info@ref_dataset
        if (ismissing(ref_ind)) then
            if (diag_script_info@ref_dataset(0).eq."n.a.") then
                ref_ind = 0
            else
                error_msg("f", diag_script, "", \
                          "no adequate reference dataset provided")
            end if
        end if
        dest_grd = find_destination_grid(ref_ind, var0, field_type0, \
                                         "coarsest")
    else if (selected_grid.eq."1x1") then
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
    else if (selected_grid.eq."2.5x2.5") then
        dest_grd     = new((/72, 144/), float)
        dest_grd!0   = "lat"
        dest_grd!1   = "lon"
        dest_grd&lat = flt2dble(fspan(-88.75, 88.75, 72))  ; CMOR standard
        dest_grd&lon = flt2dble(fspan(1.25, 358.75, 144))    ; CMOR standard
        dest_grd&lat@standard_name = "latitude"
        dest_grd&lon@standard_name = "longitude"
        dest_grd&lat@long_name     = "latitude"
        dest_grd&lon@long_name     = "longitude"
        dest_grd&lat@units         = "degrees_north"
        dest_grd&lon@units         = "degrees_east"
        dest_grd&lat@axis          = "Y"
        dest_grd&lon@axis          = "X"
    else
        dest_grd = find_destination_grid(index, var0, field_type0, \
                                         selected_grid)
    end if
    end if
    end if
    ;annots = project_style(diag_script_info, "annots")
    Fill = default_fillvalue(typeof(A0))
    var_grd = new((/dim_MOD, dim_seas, dimsizes(dest_grd&lat), \
                   dimsizes(dest_grd&lon)/), typeof(A0), Fill)
    var_grd!0 = "models"
    var_grd!1 = "season"
    var_grd!2 = "lat"
    var_grd!3 = "lon"
    var_grd&models = annots
    if (isatt(diag_script_info, "seasons")) then
        var_grd&season = (/ispan(0, dim_seas - 1, 1)/)
    else
        var_grd&season = (/0/)
    end if
    var_grd&lat = dest_grd&lat
    var_grd&lon = dest_grd&lon
    if (isatt(diag_script_info, "iavmode")) then
        if (diag_script_info@iavmode .eq. "each") then
            iav = new((/dim_MOD, dim_seas, dimsizes(dest_grd&lat), \
                       dimsizes(dest_grd&lon)/), typeof(A0), Fill)
            copy_VarMeta(var_grd, iav)
        end if
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
            error_msg("f", diag_script, "", "selected grid not possible if"\
                      + " original grid is irregular, select either 1x1 "\
                      + "or 2.5x2.5 as target grid")
        end if
        end if
    else
        regular_grid = True
    end if
      
  ; Average over time (function in ./diag_scripts/lib/ncl/statistics.ncl)
  if ((diag_script_info@time_avg .eq. "seasonalclim") .and. \
    (dim_seas .ne. 4)) then
    tmp_seas = time_operations(A0, -1, -1, "average", \
			       diag_script_info@time_avg, True)
    data1 = tmp_seas(diag_script_info@seasons, :, :)
    delete(tmp_seas)
  else
    data1 = time_operations(A0, -1, -1, "average", \
			    diag_script_info@time_avg, True)
  end if

  do imod = 0, dim_MOD - 1
    log_info(dataset_info@dataset(imod))
        if (imod .eq. 0) then
            if (dim_seas .eq. 1) then
                if (regular_grid .eq. True) then
                    var_grd(imod, 0, :, :) = linint2(data1&lon, data1&lat,\
                                                     data1, True, \
                                                     dest_grd&lon, \
                                                     dest_grd&lat, 0)
                else
                    tmp = new((/1, dimsizes(data1&$lat_name$), \
                               dimsizes(data1&$lon_name$)/), float)
                    tmp(0, :, :) = data1
                    tmp!0 = "season"
                    var_grd(imod, :, :, :) = \
                        regrid_3D_to_rectilinear_grid(tmp, lon2D, lat2D, \
                                                      grid_res, filename, \
                                                      False)
                    delete(tmp)
                end if
            else
                if (regular_grid .eq. True) then
                    var_grd(imod, :, :, :) = linint2(data1&lon, data1&lat, \
                                                     data1, True, \
                                                     dest_grd&lon,\
                                                     dest_grd&lat, 0)
                else
                    var_grd(imod, :, :, :) = \
                        regrid_3D_to_rectilinear_grid(data1, lon2D, lat2D, \
                                                      grid_res, filename, \
                                                      False)
                end if
            end if
            delete([/A0, data1/])
            if (regular_grid .eq. False) then
                delete([/lat2D, lon2D/])
            end if
        else
            A1 = extract_data(imod, data, -1, 0, 0)
            dims = getvardims(A1)
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
            ; Average over time,
            ; function in ./diag_scripts/lib/ncl/statistics.ncl)
            if ((diag_script_info@time_avg .eq. "seasonalclim") .and. \
                (dim_seas .ne. 4)) then
                tmp_seas = time_operations(A1, -1, -1, "average", \
                                           diag_script_info@time_avg, True)
                data2 = tmp_seas(diag_script_info@seasons, :, :)
                delete(tmp_seas)
            else
                data2 = time_operations(A1, -1, -1, "average", \
                                        diag_script_info@time_avg, True)
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
                    error_msg("f", diag_script, "", "selected grid not "\
                              + "possible if original grid is irregular, "\
                              + "select either 1x1 or 2.5x2.5 as target grid")
                end if
                end if
            else
                regular_grid = True
            end if
            ; Regrid to common grid
            if (dim_seas .eq. 1) then
                if (regular_grid .eq. True) then
                    var_grd(imod, 0, :, :) = linint2(data2&$lon_name$, \
                                                     data2&$lat_name$, \
                                                     data2, True, \
                                                     dest_grd&lon, \
                                                     dest_grd&lat, 0)
                else
                    tmp = new((/1, dimsizes(data2&$lat_name$), \
                               dimsizes(data2&$lon_name$)/), float)
                    tmp(0, :, :) = data2
                    var_grd(imod, :, :, :) = \
                        regrid_3D_to_rectilinear_grid(tmp, lon2D, lat2D, \
                                                      grid_res, filename, \
                                                      False)
                    delete(tmp)
                end if
            else
                if (regular_grid .eq. True) then
                    var_grd(imod, :, :, :) = linint2(data2&$lon_name$, \
                                                     data2&$lat_name$, \
                                                     data2, True, \
                                                     dest_grd&lon, \
                                                     dest_grd&lat, 0)
                else
                    var_grd(imod, :, :, :) = \
                        regrid_3D_to_rectilinear_grid(data2, lon2D, lat2D, \
                                                      grid_res, filename, \
                                                      False)
                end if
            end if
            delete(A1)
            delete(data2)
            if (regular_grid .eq. False) then
                delete(lat2D)
                delete(lon2D)
            end if
        end if
    if (isatt(diag_script_info, "iavmode")) then
      if (diag_script_info@iavmode .eq. "each") then
	; #########################################################
	; # Read natural variability for stippling for each model #
	; #########################################################
        log_info("Read natural variability for stippling for each model")
        iav_file = addfile(input_dir + "/" + "IAV_for_stippandhatch_" +\
                           "ipcc12/IAV_piControl_" + dataset_info@dataset(imod)\
                           + "_" + dataset_info@ensemble(imod) + "_" + \
                           field_type0 + "_" + var0 + "_" + \
                           diag_script_info@time_avg + ".nc", "r")
        iav_native = iav_file->iav
        tmp_regrd = linint2_Wrap(iav_native&$lon_name$, \
                                 iav_native&$lat_name$, iav_native, \
                                 True, dest_grd&lon, dest_grd&lat, 0)
	if (dim_seas .eq. 1) then
          iav(imod, 0, :, :) = rm_single_dims(tmp_regrd)
	else if (dim_seas .ne. 4) then
          iav(imod, :, :, :) = tmp_regrd(diag_script_info@seasons, :, :)
	else
          iav(imod, :, :, :) = tmp_regrd
        end if
        end if
        delete(iav_native)
      end if
    end if
  end do

  if (isatt(diag_script_info, "iavmode")) then
    if (diag_script_info@iavmode .eq. "each") then
      ; calculate multi-model mean of iav
      log_info("Calculate multi-model mean of natural variability")
      natvar_mmm = dim_avg_n(iav, 0) * sqrt(2.)
    else
      error_msg("f", diag_script, "", "Error: This iavmode is not know")
    end if
  else
    ; read already calculated mmm iav
    log_info("Read precalculated natural variability for multi-model")
    natvar_file = addfile(input_dir + "/" + "IAV_mmm_piControl_" + \
                          field_type0 + \
			  "_" + var0 + "_" + diag_script_info@time_avg + \
			  ".nc", "r")
    natvar_mmm = natvar_file->iav
    if ((dim_seas .ne. 4) .and. \
      (diag_script_info@time_avg .ne. "annualclim")) then
      tmp = natvar_mmm
      delete(natvar_mmm)
      natvar_mmm = new((/dimsizes(diag_script_info@seasons), \
			dimsizes(dest_grd&lat), \
			dimsizes(dest_grd&lon)/), float)
      natvar_mmm(:, :, :) = tmp(diag_script_info@seasons, :, :)
    end if
  end if

  ; ##############################################
  ; # Calculate change from hist to periods      #
  ; ##############################################
  ; loop over rcps and periods
  dim_scen = dimsizes(diag_script_info@scenarios)
  dim_per = dimsizes(diag_script_info@periods)

  to_plot = new((/dim_scen * dim_per * dim_seas, dimsizes(dest_grd&lat), \
	         dimsizes(dest_grd&lon)/), typeof(var_grd), Fill)
  if (isatt(diag_script_info, "percent")) then
    percent = diag_script_info@percent
    ref = new((/dim_scen * dim_per * dim_seas, dimsizes(dest_grd&lat), \
	       dimsizes(dest_grd&lon)/), typeof(var_grd), Fill)
  else
    percent = 0
  end if
  to_plot_signif = new((/dim_scen * dim_per * dim_seas,  \
		        dimsizes(dest_grd&lat), dimsizes(dest_grd&lon)/), \
		       typeof(var_grd), Fill)
  to_plot_not_signif = new((/dim_scen * dim_per * dim_seas, \
			    dimsizes(dest_grd&lat), \
			    dimsizes(dest_grd&lon)/), typeof(var_grd), Fill)
  model_number = new((/dim_scen * dim_per * dim_seas/), integer)
  if (var0 .eq. "mrsos" .or. var0 .eq. "mrro") then
    ; read mask to mask data where not defined (ocean and ice)
    in = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
    msk_data = in->LSMASK
    msk = landsea_mask(msk_data, dest_grd&lat, dest_grd&lon)
  end if
  i = 0
  do seas = 0, dim_seas - 1  ; loop over seasons
    do rcp = 0, dim_scen - 1  ; loop over rcps
      do per = 0, dim_per - 1  ; loop over periods
	idx_rcp = ind((dataset_info@exp .eq. \
		       diag_script_info@scenarios(rcp)) .and. \
		      (dataset_info@start_year .eq. \
		       diag_script_info@periods(per)))
	proj = var_grd(idx_rcp, seas, :, :)
	; **cut low values ***
	proj = where(abs(proj) .gt. 1.e-14, proj, proj@_FillValue)
	proj_avg = dim_avg_n_Wrap(proj, 0)
	; find historical runs from same datasets as in rcp
	idx_hist = find_model_index(dataset_info@dataset(idx_rcp), \
				    dataset_info@dataset, "historical",\
				    dataset_info@exp, \
				    dataset_info@ensemble(idx_rcp), \
				    dataset_info@ensemble, 0)
	base = var_grd(idx_hist, seas, :, :)
	; **cut low values ***
	base = where(abs(base) .gt. 1.e-14, base, base@_FillValue)
	base_avg = dim_avg_n_Wrap(base, 0)
	var_diff = proj - base
	delete([/base, proj/])
	to_plot(i, :, :) = proj_avg - base_avg
	if (isatt(diag_script_info, "percent")) then
	  ref(i, :, :) = base_avg
	end if
	; determine sigma and sign for significance
	var_pos_signif = dim_num_n(where(var_diff .gt. 0., 1., \
					 to_plot@_FillValue) .eq. 1., 0)
	var_neg_signif = dim_num_n(where(var_diff .lt. 0., 1., \
					 to_plot@_FillValue) .eq. 1., 0)
	var_pos_signif_tmp = \
	    var_pos_signif / where((var_pos_signif + var_neg_signif) \
				   .ne. 0., var_pos_signif + \
				   var_neg_signif, to_plot@_FillValue)
	var_neg_signif_tmp = \
	    var_neg_signif / where((var_pos_signif + var_neg_signif) \
				   .ne. 0., var_pos_signif + \
				   var_neg_signif, to_plot@_FillValue)
	signif_and_sign = where(var_pos_signif_tmp .ge. 0.9 .or. \
				var_neg_signif_tmp .ge. 0.9, 1., 0.)
	signif_and_sign = where(ismissing(signif_and_sign), 0.,\
				signif_and_sign)
	if (dim_seas .ne. 1) then
	  sigma_and_sign = where(abs(to_plot(i, :, :)) .gt. \
				 abs(2. * natvar_mmm(seas, :, :)), 1., 0.)
	  to_plot_not_signif(i, :, :) = where(abs(to_plot(i, :, :)) .lt. \
                                              abs(natvar_mmm(seas, :, :)),\
                                              1., 0.)
	else
	  sigma_and_sign = where(abs(to_plot(i, :, :)) .gt. \
				 abs(2. * natvar_mmm(0, :, :)), 1., 0.)
	  to_plot_not_signif(i, :, :) = where(abs(to_plot(i, :, :)) .lt. \
                                              abs(natvar_mmm(0, :, :)),	1., 0.)
	end if
	to_plot_signif(i, :, :) = where(signif_and_sign .eq. 1, \
					sigma_and_sign, 0.)
	if (var0 .eq. "mrsos" .or. var0 .eq. "mrro") then
	  to_plot(i, :, :) = mask(to_plot(i, :, :), msk .eq. 0, False)
	  to_plot(i, :, :) = mask(to_plot(i, :, :), msk .eq. 4, False)
	  to_plot_signif(i, :, :) = mask(to_plot_signif(i, :, :), \
					 msk .eq. 0, False)
	  to_plot_signif(i, :, :) = mask(to_plot_signif(i, :, :), \
					 msk .eq. 4, False)
          to_plot_not_signif(i, :, :) = mask(to_plot_not_signif(i, :, :), \
                                             msk .eq. 0, False)
          to_plot_not_signif(i, :, :) = mask(to_plot_not_signif(i, :, :), \
                                             msk .eq. 4, False)
	end if
	model_number(i) = dimsizes(idx_rcp)
	i = i + 1
	delete([/idx_hist, idx_rcp, var_diff/])
      end do  ; per
    end do  ; rcp
  end do  ; seas
  to_plot!0 = "panel"
  to_plot&panel = diag_script_info@label
  to_plot!1 = "lat"
  to_plot&lat = dest_grd&lat
  to_plot!2 = "lon"
  to_plot&lon = dest_grd&lon
  copy_VarMeta(to_plot, to_plot_signif)
  copy_VarMeta(to_plot, to_plot_not_signif)
  if (isatt(diag_script_info, "percent")) then
    percent = diag_script_info@percent
    if (percent .eq. 1) then
      to_plot = 100 * to_plot / where(ref .ne. 0., ref, ref@_FillValue)
      to_plot@units     = "%"
      delete(ref)
    end if
  end if
  ; ###########################################
  ; # Other Metadata: diag_script, var        #
  ; ###########################################
  ; Add to to_plot, as attributes without prefix
  if (isatt(to_plot, "diag_script")) then  ; Add to existing entries
    temp = to_plot@diag_script
    delete(to_plot@diag_script)
    to_plot@diag_script = array_append_record(temp, (/diag_script/), 0)
    delete(temp)
  else  ; Add as new attribute
    to_plot@diag_script = (/diag_script/)
  end if
  to_plot@var = var0  ; Overwrite existing entry
  if (isatt(variable_info, "long_name")) then
    to_plot@var_long_name = variable_info@long_name
  end if
  ; check unit and adjust for plotting if necessary
  if (isatt(variable_info, "units")) then
    if ((var0 .eq. "mrro" .or. var0 .eq. "evspsbl") \
	.and. variable_info@units .eq. "kg m-2 s-1") then
      ; convert [kg m-2 s-1] to [mm day-1]
      to_plot = to_plot * 60 * 60 * 24
      variable_info@units = "mm day-1"
      to_plot@var_units   = variable_info@units
    else if (var0 .eq. "mrsos" .and. \
	     variable_info@units .eq. "kg m-2") then
      ; convert [kg m-2] to [mm]
      variable_info@units = "mm"
      to_plot@var_units   = variable_info@units
    else if (var0 .eq. "hfls" .and. variable_info@units .eq. "W m-2") then
      to_plot = to_plot / (2.5 * 10 ^ 6) * 24 * 60 * 60
      variable_info@units = "mm day-1"
      to_plot@var_units   = variable_info@units
    else if (variable_info@units .eq. "Pa") then
      ; convert to hPa
      to_plot = to_plot / 100
      to_plot@var_units     = "hPa"
    else
      to_plot@var_units     = variable_info@units
    end if
    end if
    end if
    end if
  end if
  ; swap downward/upward radiation fluxes according to IPCC Chapter12
  if ((var0 .eq. "rsut") .or. (var0 .eq. "rlut")) then
    to_plot = - to_plot
  end if

  if(isatt(to_plot_signif, "diag_script")) then  ; Add to existing entries
    temp = to_plot_signif@diag_script
    delete(to_plot_signif@diag_script)
    to_plot_signif@diag_script = array_append_record(temp, \
						     (/diag_script/), 0)
    delete(temp)
  else  ; Add as new attribute
    to_plot_signif@diag_script = (/diag_script/)
  end if

  to_plot_signif@var = "signif"
  to_plot_signif@var_long_name = "significance: mean change larger than "\
			         + "2*natvar and 90% of models have the " \
			         + "same sign"
  to_plot_signif@var_units     = 1

  if(isatt(to_plot_not_signif, "diag_script")) then
    temp = to_plot_not_signif@diag_script
    delete(to_plot_not_signif@diag_script)
    to_plot_not_signif@diag_script = array_append_record(temp, \
							 (/diag_script/), 0)
    delete(temp)
  else  ; Add as new attribute
    to_plot_not_signif@diag_script = (/diag_script/)
  end if
  to_plot_not_signif@var = "not_signif"
  to_plot_not_signif@var_long_name = "not significant: mean change " \
    				     + "smaller than 1*natvar"
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
  ; # Optional output to netCDF               #
  ; ###########################################
  ; Optional output (controlled by diag_script_info)
  if (isatt(diag_script_info, "ncdf")) then
    to_plot@ncdf = diag_script_info@ncdf
    ncdf_file = netcdf_dir + "/" + diag_script_info@ncdf
    ncdf_file@existing = "overwrite"
    ; Function in ~/interface_scripts/auxiliary.ncl
    ncdf_outfile = ncdf_write(to_plot, ncdf_file)
    ncdf_file@existing = "append"
    ncdf_outfile = ncdf_write(to_plot_signif, ncdf_file)
    ncdf_outfile = ncdf_write(to_plot_not_signif, ncdf_file)
    ncdf_outfile = ncdf_write(model_number, ncdf_file)
  end if
end
