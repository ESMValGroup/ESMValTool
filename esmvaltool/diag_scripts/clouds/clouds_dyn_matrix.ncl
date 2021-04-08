; CLOUDS_DYN_MATRIX
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; ############################################################################
; Description
;   Calculates mean values of variable z per bin of variable x and y.
;   The results are displayed as a matrix.
;
; Required diag_script_info attributes (diagnostic specific)
;   var_x: short name of variable on x-axis
;   var_y: short name of variable on y-axis
;   var_z: short name of variable to be binned
;   xmin:  min x value for generating x bins
;   xmax:  max x value for generating x bins
;   ymin:  min y value for generating y bins
;   ymax:  max y value for generating y bins
;
; Optional diag_script_info attributes (diagnostic specific)
;   filename_add: optionally add this string to plot filesnames
;   nbins: number of equally spaced bins (var_x), default = 20
;   zdmin:  min z value for labelbar (difference plots)
;   zdmax:  max z value for labelbar (difference plots)
;   zmin:  min z value for labelbar
;   zmax:  max z value for labelbar
;
; Required variable attributes (variable specific)
;   none
;
; Optional variable_info attributes (variable specific)
;   reference_dataset: reference dataset
;
; Caveats
;   none
;
; Modification history
;   20210408-lauer_axel: written.
;
; ############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/style.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  diag = "clouds_dyn_matrix.ncl"
  variables = metadata_att_as_array(variable_info, "short_name")

  ; Check required diag_script_info attributes
  exit_if_missing_atts(diag_script_info, (/"var_x", "var_y", "var_z", \
                                           "xmin", "xmax", "ymin", "ymax"/))

  file_type = output_type()

  ; make sure required variables are available
  var_x = diag_script_info@var_x
  var_y = diag_script_info@var_y
  var_z = diag_script_info@var_z

  ; special case: columnicefrac = clivi / (clivi + lwp)
  ; note: clwvi is not used since it contains lwp only for some models
  ;       (by error)

  calcicefrac = False
  calctcwp = False
  if (var_z .eq. "columnicefrac") then
    calcicefrac = True
    varlist = (/var_x, var_y, "clivi", "lwp"/)
  else if (var_z .eq. "totalcwp") then
    calctcwp = True
    varlist = (/var_x, "var_y", "clivi", "lwp"/)
  else
    varlist = (/var_x, var_y, var_z/)
  end if
  end if
  idx = new(dimsizes(varlist), integer)

  nVAR = dimsizes(varlist)
  refname = new(nVAR, string)

  do i = 0, nVAR - 1
    idx(i) = ind(variables .eq. varlist(i))
  end do

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + variables(idx) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  if (any(ismissing(idx))) then
    errstr = "diagnostic " + diag + " requires the following variable(s): " \
             + str_join(varlist, ", ")
    error_msg("f", DIAG_SCRIPT, "", errstr)
  end if

  ; save input files for writing provenance

  infiles = metadata_att_as_array(input_file_info, "filename")

  ; get reference datasets (if present) and check that number of datasets
  ; is equal for each variable

  do i = 0, nVAR - 1
    if (isatt(variable_info[idx(i)], "reference_dataset")) then
      refname(i) = variable_info[idx(i)]@reference_dataset
    end if
    info = select_metadata_by_name(input_file_info, variables(idx(i)))
    if (i .eq. 0) then
      dim_MOD = ListCount(info)
    else
      dim_test = ListCount(info)
      if (dim_test .ne. dim_MOD) then
        error_msg("f", DIAG_SCRIPT, "", "number of datasets for variable " \
                  + variables(i) + " does not match number of datasets for " \
                  + variables(0))
      end if
    end if
    delete(info)
  end do

  ; Set default values for non-required diag_script_info attributes

  set_default_att(diag_script_info, "filename_add", "")
  set_default_att(diag_script_info, "nbins", 20)

  if (diag_script_info@filename_add .ne. "") then
    filename_add = "_" + diag_script_info@filename_add
  else
    filename_add = ""
  end if

  nbins = toint(diag_script_info@nbins)

  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

end

begin
  ; ############
  ; # get data #
  ; ############

  info_x = select_metadata_by_name(input_file_info, varlist(0))
  names_x = metadata_att_as_array(info_x, "dataset")
  projects_x = metadata_att_as_array(info_x, "project")
  info_y = select_metadata_by_name(input_file_info, varlist(1))
  names_y = metadata_att_as_array(info_y, "dataset")
  projects_y = metadata_att_as_array(info_y, "project")
  info_z = select_metadata_by_name(input_file_info, varlist(2))
  names_z = metadata_att_as_array(info_z, "dataset")
  projects_z = metadata_att_as_array(info_z, "project")

  refidx_x = ind(names_x .eq. refname(0))
  refidx_y = ind(names_y .eq. refname(1))
  refidx_z = ind(names_z .eq. refname(2))

  if (ismissing(refidx_x) .or. ismissing(refidx_y) .or. ismissing(refidx_z)) \
    then
    refidx_x = -1
    refidx_y = -1
    refidx_z = -1
  end if

  result = new((/dim_MOD, nbins, nbins/), float)
  bincenter_x = new((/nbins/), float)
  bincenter_y = new((/nbins/), float)

  xmax = diag_script_info@xmax
  xmin = diag_script_info@xmin
  ymax = diag_script_info@ymax
  ymin = diag_script_info@ymin
  binsize_x = tofloat(xmax - xmin) / nbins
  binsize_y = tofloat(ymax - ymin) / nbins

  do n = 0, nbins - 1
    x0 = n * binsize_x
    x1 = x0 + binsize_x
    bincenter_x(n) = xmin + 0.5 * (x0 + x1)
    y0 = n * binsize_y
    y1 = y0 + binsize_y
    bincenter_y(n) = ymin + 0.5 * (y0 + y1)
  end do

  do ii = 0, dim_MOD - 1
    atts_x = True
    atts_x@short_name = varlist(0)

    atts_y = True
    atts_y@short_name = varlist(1)

    atts_z = True
    atts_z@short_name = varlist(2)

    ; reference datasets may have different names
    if (ii .eq. refidx_x) then
      atts_x@dataset = refname(0)
      atts_y@dataset = refname(1)
      atts_z@dataset = refname(2)
    else  ; all other datasets: force same dataset name for var_x, var_y, var_z
      atts_x@dataset = names_x(ii)
      atts_y@dataset = names_x(ii)
      atts_z@dataset = names_x(ii)
    end if

    ; read var_x

    info = select_metadata_by_atts(input_file_info, atts_x)
    x = read_data(info[0])
    delete(info)

    ; read var_y

    info = select_metadata_by_atts(input_file_info, atts_y)
    y = read_data(info[0])
    delete(info)

    ; read var_z

    info = select_metadata_by_atts(input_file_info, atts_z)
    z = read_data(info[0])
    delete(info)

    if (calcicefrac) then
      atts_z@short_name = varlist(3)
      info = select_metadata_by_atts(input_file_info, atts_z)
      z2 = read_data(info[0])
      delete(info)
      min_mass = 1.0e-6
      ; filter valid values (needed for some models)
      z = where(z .lt. 0.0, z@_FillValue, z)
      z = where(isnan_ieee(z), y@_FillValue, z)
      z2 = where(z2 .lt. 0.0, z2@_FillValue, z2)
      z2 = where(isnan_ieee(z2), z@_FillValue, z2)
      mass = z + z2
      delete(z)
      delete(z2)
      mass = where(mass .lt. min_mass, mass@_FillValue, mass)

      ; ice fraction = ice / (ice + lwp) * 100%
      z = 100.0 * z / mass

      z@units = "%"
      z@long_name = "cloud ice fraction"
      z@var = "columnicefrac"
    end if

    ; calculate total cloud water path as sum of liquid water path (lwp)
    ; and ice water path (clivi);
    ; we do not use the CMOR variable clwvi directly as this variable
    ; erroneously contains only cloud liquid water for some models

    if (calctcwp) then
      atts_z@short_name = varlist(3)
      info = select_metadata_by_atts(input_file_info, atts_z)
      z2 = read_data(info[0])
      delete(info)
      z = z + z2
      delete(z2)
      z@long_name = "Condensed Water Path"
      z@var = "totalcwp"
    end if

    ; check dimensions

    dims_x = dimsizes(x)
    dims_y = dimsizes(y)
    dims_z = dimsizes(z)

    dimerror = False

    if (dimsizes(dims_x) .eq. dimsizes(dims_y)) then
      if (any(dims_x - dims_y .ne. 0)) then
        dimerror = True
      end if
    else
      dimerror = True
    end if

    if (dimsizes(dims_x) .eq. dimsizes(dims_z)) then
      if (any(dims_x - dims_z .ne. 0)) then
        dimerror = True
      end if
    else
      dimerror = True
    end if

    if (dimerror) then
      error_msg("f", DIAG_SCRIPT, "", "dimensions of datasets " \
                + atts_x@dataset + " (variable " + var_x + ") and " \
                + atts_y@dataset + " (variable " + var_y + ") and " \
                + atts_z@dataset + " (variable " + var_z + ") do not match.")
    end if

    ; check dimensions

    if ((dimsizes(dims_x) .ne. dimsizes(dims_y)) .or. \
        (dimsizes(dims_x) .lt. 3) .or. (dimsizes(dims_x) .gt. 4)) then
      dimerror = True
    end if

    if ((dimsizes(dims_y) .ne. dimsizes(dims_z)) .or. \
        (dimsizes(dims_y) .lt. 3) .or. (dimsizes(dims_y) .gt. 4)) then
      dimerror = True
    end if

    if (dimerror) then
      error_msg("f", DIAG_SCRIPT, "", "all variables need to have the " + \
                "same number of dimensions (time, [optional: level], " + \
                "latitude, longitude)")
    end if

    do i = 0, nVAR - 1
      var = variables(idx(i))
      if (var .eq. varlist(0)) then
        dims = getvardims(x)
      else if (var .eq. varlist(1)) then
        dims = getvardims(y)
      else
        dims = getvardims(z)
      end if
      end if
      testidx = ind(dims .eq. "lon")
      if (ismissing(testidx)) then
        error_msg("f", DIAG_SCRIPT, "", var + ": no lon dimension")
      end if
      testidx = ind(dims .eq. "lat")
      if (ismissing(testidx)) then
        error_msg("f", DIAG_SCRIPT, "", var + ": no lat dimension")
      end if
      testidx = ind(dims .eq. "time")
      if (ismissing(testidx)) then
        error_msg("f", DIAG_SCRIPT, "", var + ": no time dimension")
      end if
      delete(dims)
    end do

    delete(dims_x)
    delete(dims_y)
    delete(dims_z)
    delete(testidx)

    ref_ind = refidx_x
    if (ismissing(ref_ind)) then
      ref_ind = -1
    end if
    names = names_x
    projects = projects_x

    if (refidx_x .ge. 0) then
      ; if reference datasets for var_x and var_y are from different sources
      if ((refname(1) .ne. refname(0)) .or. (refname(1) .ne. refname(2))) then
        names(refidx_x) = names(refidx_x) + "/" + refname(1)
      end if
      if ((refname(2) .ne. refname(0)) .and. (refname(2) .ne. refname(1))) then
        names(refidx_x) = names(refidx_x) + "/" + refname(2)
      end if
    end if

    ; save attributes long_name and units
    long_name = z@long_name
    xunits = x@units
    yunits = y@units
    zunits = z@units

    x1d = ndtooned(x)
    delete(x)
    y1d = ndtooned(y)
    delete(y)
    z1d = ndtooned(z)
    delete(z)

    r = bin_avg(x1d, y1d, z1d, bincenter_x, bincenter_y, False)
    result(ii, :, :) = r(0, :, :)

    delete(x1d)
    delete(y1d)
    delete(z1d)
  end do  ; ii-loop (models)

  delete(r)

  ; ###########################################
  ; # netCDF output                           #
  ; ###########################################

  nc_filename = work_dir + "clouds_scatter_" + var_x + "_" + var_y + "_" + \
    var_z + filename_add + ".nc"

  result!0 = "model"
  result!1 = "bin_y"
  result!2 = "bin_x"

  result&model = str_sub_str(names, "/", "-")
  result&bin_x = bincenter_x
  result&bin_y = bincenter_y

  result@diag_script = (/DIAG_SCRIPT/)
  result@var = var_z
  result@var_long_name = long_name
  result@var_units = zunits
  result@_FillValue = 1.0e20

  nc_outfile = ncdf_write(result, nc_filename)

  ; ###########################################
  ; # create the plots                        #
  ; ###########################################

  plots = new(dim_MOD, graphic)
  dplots = new(dim_MOD, graphic)

  res = True
  res@cnFillOn           = True               ; color Fill
  res@cnFillMode         = "RasterFill"       ; Raster Mode
  res@cnLinesOn          =  False             ; Turn off contour lines
  res@tiXAxisString      = var_x + " (" + xunits + ")"
  res@tiYAxisString      = var_y + " (" + yunits + ")"
  res@lbOrientation      = "Vertical"

  dres = res

  if (isatt(diag_script_info, "zmin") .and. isatt(diag_script_info, "zmax")) \
    then
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = fspan(diag_script_info@zmin, diag_script_info@zmax, 19)
  end if

  if (isatt(diag_script_info, "zdmin") .and. \
      isatt(diag_script_info, "zdmax")) then
    dres@cnLevelSelectionMode = "ExplicitLevels"
    dres@cnLevels = fspan(diag_script_info@zdmin, diag_script_info@zdmax, 19)
  end if

  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_dyn_matrix_" + \
                var_x + "_" + var_y + "_" + var_z + filename_add)
  dwks = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_dyn_matrix_bias_" + \
                 var_x + "_" + var_y + "_" + var_z + filename_add)

  do ii = 0, dim_MOD - 1
    res@tiMainString =  names(ii)
    plots(ii) = gsn_csm_contour(wks, result(ii, :, :), res)
    if (isdefined("dresplot")) then
      delete(dresplot)
    end if
    if (ii .eq. ref_ind) then
      dresplot = res
      dresplot@tiMainString =  "REF"
      diff = result(ii, :, :)
    else
      dresplot = dres
      dresplot@tiMainString =  names(ii) + " - REF"
      diff = result(ii, :, :) - result(ref_ind, :, :)
    end if
    dplots(ii) = gsn_csm_contour(dwks, diff, dresplot)
  end do

  pres                      = True    ; needed to override
                                      ; panelling defaults
  pres@gsnPanelCenter                   = False

  pres@gsnPanelFigureStrings = names
  pres@gsnPanelFigureStringsFontHeightF = min((/0.008, 0.008 * 6.0 \
                                              / tofloat((dim_MOD + 1) / 2)/))
  pres@lbLabelFontHeightF               = min((/0.01, 0.01 * 6.0 \
                                              / tofloat((dim_MOD + 1) / 2)/))
  outfile = panelling(wks, plots, (dim_MOD + 3) / 4, 4, pres)
  doutfile = panelling(dwks, dplots, (dim_MOD + 3) / 4, 4, pres)

  log_info("Wrote " + outfile)
  log_info("Wrote " + doutfile)

  ; ==========================================================================

  ; ----------------------------------------------------------------------
  ; write provenance to netcdf output (and plot file)
  ; ----------------------------------------------------------------------

  statistics = (/"clim", "mean"/)
  domain = "reg"
  plottype = "scatter"
  caption = "Scatterplot of " + var_x + " (x) vs. " + var_y + " (y)."
  log_provenance(nc_outfile, outfile, caption, statistics, \
                 domain, plottype, "", "", infiles)

  leave_msg(DIAG_SCRIPT, "")

end
