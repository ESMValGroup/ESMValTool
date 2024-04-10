; CLOUDS_SCATTER
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; ############################################################################
; Description
;   Calculates mean values of variable y per bin of variable x and displays
;   the results as scatter plot.
;
; Required diag_script_info attributes (diagnostic specific)
;   var_x: short name of variable on x-axis
;   var_y: short name of variable on y-axis
;   xmin:  min x value for generating bins
;   xmax:  max x value for generating bins
;
; Optional diag_script_info attributes (diagnostic specific)
;   filename_add: optionally add this string to plot filesnames
;   nbins: number of equally spaced bins (var_x), default = 20
;   ymin_mm:  min y value for plotting MultiModelMean
;   ymax_mm:  max y value for plotting MultiModelMean
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
;   20230117-lauer_axel: added support for ICON (code from Manuel)
;   20210210-lauer_axel: written.
;
; ############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/dataset_selection.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/style.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  diag = "clouds_scatter.ncl"
  variables = get_unique_values(metadata_att_as_array(variable_info, \
                                                      "short_name"))

  ; Check required diag_script_info attributes
  exit_if_missing_atts(diag_script_info, (/"var_x", "var_y", "xmin", "xmax"/))

  file_type = output_type()

  ; make sure required variables are available
  var_x = diag_script_info@var_x
  var_y = diag_script_info@var_y

  ; special case: columnicefrac = clivi / (clivi + lwp)
  ; note: clwvi is not used since it contains lwp only for some models
  ;       (by error)

  calcicefrac = False
  calctcwp = False
  if (var_y .eq. "columnicefrac") then
    calcicefrac = True
    varlist = (/var_x, "clivi", "lwp"/)
  else if (var_y .eq. "totalcwp") then
    calctcwp = True
    varlist = (/var_x, "clivi", "lwp"/)
  else
    varlist = (/var_x, var_y/)
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
    var = variables(idx(i))
    var_info = select_metadata_by_name(variable_info, var)
    var_info := var_info[0]
    if (isatt(var_info, "reference_dataset")) then
      refname(i) = var_info@reference_dataset
    end if
    info = select_metadata_by_name(input_file_info, var)
    if (i .eq. 0) then
      dim_MOD = ListCount(info)
    else
      dim_test = ListCount(info)
      if (dim_test .ne. dim_MOD) then
        error_msg("f", DIAG_SCRIPT, "", "number of datasets for variable " \
                  + var + " does not match number of datasets for " \
                  + variables(idx(0)))
      end if
    end if
    delete(info)
    delete(var)
    delete(var_info)
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

  refidx_x = ind(names_x .eq. refname(0))
  refidx_y = ind(names_y .eq. refname(1))

  if (ismissing(refidx_x) .or. ismissing(refidx_y)) then
    refidx_x = -1
    refidx_y = -1
  end if

  result_avg = new((/dim_MOD, nbins/), float)
  result_std = new((/dim_MOD, nbins/), float)
  bincenter = new((/nbins/), float)

  xmax = diag_script_info@xmax
  xmin = diag_script_info@xmin
  binsize = tofloat(xmax - xmin) / nbins

  do n = 0, nbins - 1
    x0 = n * binsize
    x1 = x0 + binsize
    bincenter(n) = xmin + 0.5 * (x0 + x1)
  end do

  do ii = 0, dim_MOD - 1
    atts_x = True
    atts_x@short_name = varlist(0)

    atts_y = True
    atts_y@short_name = varlist(1)

    ; reference datasets may have different names
    if (ii .eq. refidx_x) then
      atts_y@dataset = refname(1)
      atts_x@dataset = refname(0)
    ; all other datasets: force same dataset name for var_x and var_y
    else
      atts_y@dataset = names_x(ii)
      atts_x@dataset = names_x(ii)
    end if

    ; read var_x

    info = select_metadata_by_atts(input_file_info, atts_x)
    x = read_data(info[0])
    delete(info)

    ; read var_y

    info = select_metadata_by_atts(input_file_info, atts_y)
    y = read_data(info[0])
    delete(info)

    if (calcicefrac) then
      atts_y@short_name = varlist(2)
      info = select_metadata_by_atts(input_file_info, atts_y)
      z = read_data(info[0])
      delete(info)
      min_mass = 1.0e-6
      ; filter valid values (needed for some models)
      y = where(y .lt. 0.0, y@_FillValue, y)
      y = where(isnan_ieee(y), y@_FillValue, y)
      z = where(z .lt. 0.0, z@_FillValue, z)
      z = where(isnan_ieee(z), z@_FillValue, z)
      mass = y + z
      delete(z)
      mass = where(mass .lt. min_mass, mass@_FillValue, mass)

      ; ice fraction = ice / (ice + lwp) * 100%
      y = 100.0 * y / mass
      delete(mass)

      y@units = "%"
      y@long_name = "cloud ice fraction"
      y@var = "columnicefrac"
    end if

    ; calculate total cloud water path as sum of liquid water path (lwp)
    ; and ice water path (clivi);
    ; we do not use the CMOR variable clwvi directly as this variable
    ; erroneously contains only cloud liquid water for some models

    if (calctcwp) then
      atts_y@short_name = varlist(2)
      info = select_metadata_by_atts(input_file_info, atts_y)
      z = read_data(info[0])
      delete(info)
      y = y + z
      delete(z)
      y@long_name = "Condensed Water Path"
      y@var = "totalcwp"
    end if

    ; check dimensions

    dims_x = dimsizes(x)
    dims_y = dimsizes(y)

    dimerror = False

    if (dimsizes(dims_x) .eq. dimsizes(dims_y)) then
      if (any(dims_x - dims_y .ne. 0)) then
        dimerror = True
      end if
    else
      dimerror = True
    end if

    if (dimerror) then
      error_msg("f", DIAG_SCRIPT, "", "dimensions of datasets " \
                + atts_x@dataset + " (variable " + var_x + ") and " \
                + atts_y@dataset + " (variable " + var_y + ") do not match.")
    end if

    ; check dimensions

    if ((dimsizes(dims_x) .ne. dimsizes(dims_y)) .or. \
        (dimsizes(dims_x) .lt. 3) .or. (dimsizes(dims_x) .gt. 4)) then
      error_msg("f", DIAG_SCRIPT, "", "all variables need to have the " + \
                "same number of dimensions (time, [optional: level], " + \
                "latitude, longitude)")
    end if

    do i = 0, nVAR - 1
      var = variables(idx(i))
      if (var .eq. varlist(0)) then
        dims = getvardims(x)
      else
        dims = getvardims(y)
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
    delete(testidx)

    ref_ind = refidx_x
    if (ismissing(ref_ind)) then
      ref_ind = -1
    end if
    names = names_x
    projects = projects_x

    if (refidx_x .ge. 0) then
      ; if reference datasets for var_x and var_y are from different sources
      if (refname(0) .ne. refname(1)) then
        names(refidx_x) = refname(0) + "/" + refname(1)
      end if
    end if

    ; save attributes long_name and units
    long_name = y@long_name
    units = y@units
    xunits = x@units

    x1d = ndtooned(x)
    delete(x)
    y1d = ndtooned(y)
    delete(y)

    do n = 0, nbins - 1
      x0 = xmin + n * binsize
      x1 = x0 + binsize
      idx0 = ind((x1d .gt. x0) .and. (x1d .le. x1))
      if (.not.all(ismissing(idx0))) then
        result_avg(ii, n) = avg(y1d(idx0))
        result_std(ii, n) = stddev(y1d(idx0))

      else
        result_avg(ii, n) = result_avg@_FillValue
        result_std(ii, n) = result_std@_FillValue
      end if
      delete(idx0)
    end do

    delete(x1d)
    delete(y1d)
  end do  ; ii-loop (models)

  ; if multiple models are present, calculate standard deviation of all models

  ; find all indices of models w/o MultiModelMean/MultiModelMedian (if present)

  idxmod = get_mod(names, projects)

  if (idxmod(0) .eq. -1) then
    flag_multimod = False
    mm_ind = -1
  elseif (dimsizes(idxmod) .eq. 1) then
    flag_multimod = False
    mm_ind = -1
  else
    flag_multimod = True
    mmavg = new((/1, nbins/), float)
    mmstd = new((/1, nbins/), float)
    mmp10 = new((/1, nbins/), float)
    mmp90 = new((/1, nbins/), float)
    do n = 0, nbins - 1
      mmavg(0, n) = avg(result_avg(idxmod, n))
      mmstd(0, n) = stddev(result_avg(idxmod, n))
      selection = result_avg(idxmod, n)
      itmp = ind(.not.ismissing(selection))
      if (.not. ismissing(itmp(0))) then
        sorted = selection(itmp)
        qsort(sorted)
        i10 = toint(dimsizes(sorted) * 0.1 + 0.5)
        i90 = toint(dimsizes(sorted) * 0.9 - 0.5)
        mmp10(0, n) = sorted(i10)
        mmp90(0, n) = sorted(i90)
        delete(sorted)
      else
        mmp10(0, n) = mmp10@_FillValue
        mmp90(0, n) = mmp90@_FillValue
      end if
      delete(selection)
      delete(itmp)
    end do
    mm_ind = dim_MOD
    dim_MOD = dim_MOD + 1
    result_avg := array_append_record(result_avg, mmavg, 0)
    result_std := array_append_record(result_std, mmstd, 0)
    result_std(dim_MOD - 1, :) = 0.0
    names := array_append_record(names, (/"Multi-model average"/), 0)
  end if

  ; ###########################################
  ; # netCDF output                           #
  ; ###########################################

  nc_filename = work_dir + "clouds_scatter_" + var_x + "_" + var_y + \
    filename_add + ".nc"

  result_avg!0 = "model"
  result_avg!1 = "bin"

  result_avg&model = str_sub_str(names, "/", "-")
  result_avg&bin = bincenter

  result_avg@diag_script = (/DIAG_SCRIPT/)
  result_avg@var = var_y
  result_avg@var_long_name = long_name
  result_avg@var_units = units

  nc_outfile = ncdf_write(result_avg, nc_filename)

  ; ###########################################
  ; # create the plots                        #
  ; ###########################################

  plots = new(dim_MOD, graphic)
  stdbar = new((/nbins, dim_MOD/), graphic)
  centers = new((/nbins, dim_MOD/), graphic)
  centersout = new((/nbins, dim_MOD/), graphic)
  stdbarR = new((/nbins, dim_MOD/), graphic)
  centersR = new((/nbins, dim_MOD/), graphic)
  centersRout = new((/nbins, dim_MOD/), graphic)
  res = True

  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_scatter_" + \
                var_x + "_" + var_y + filename_add)

  data = new((/2, nbins/), float)
  if (ref_ind .gt. 0) then
    data(0, :) = result_avg(ref_ind, :)
  else
    data(0, :) = data@_FillValue
  end if

  do ii = 0, dim_MOD - 1
    if (ii .eq. refidx_x) then
      continue
    end if

    res@gsnDraw            = False  ; do not draw yet
    res@gsnFrame           = False  ; don't advance frame
    res@xyMarkLineMode     = "MarkLines"
    res@xyDashPatterns     = (/0., 0./)
    res@tmLabelAutoStride  = True
    res@xyLineThicknesses  = (/2.0, 2.0/)
    res@xyLineColors       = (/"black", "red"/)
    res@tiMainFontHeightF  = 0.025
    res@tiYAxisFontHeightF = 0.025
    res@tiXAxisFontHeightF = 0.025
    res@tiXAxisString      = var_x + " (" + xunits + ")"
    res@tiYAxisString      = var_y + " (" + units + ")"

    if (ii .eq. mm_ind) then
      if (isatt(diag_script_info, "ymin_mm")) then
        res@trYMinF = diag_script_info@ymin_mm
      end if
      if (isatt(diag_script_info, "ymax_mm")) then
        res@trYMaxF = diag_script_info@ymax_mm
      end if
    else
      if (isatt(res, "trYMinF")) then
        delete(res@trYMinF)
      end if
      if (isatt(res, "trYMaxF")) then
        delete(res@trYMaxF)
      end if
    end if

    polyres                  = True
    polyres@gsMarkerSizeF    = 0.01
    polyres@gsLineColor      = "red"
    polyres@gsLineThicknessF = 1.0

    polyresRef = True
    polyresRef@gsMarkerSizeF    = 0.01
    polyresRef@gsLineColor      = "black"
    polyresRef@gsLineThicknessF = 1.0

    data(1, :) = result_avg(ii, :)

    res@tiMainString = names(ii)

    plots(ii) = gsn_csm_xy(wks, result_avg&bin, data, res)

    if (ii .eq. mm_ind) then
      res_std = True
      res_std@gsnDraw         = False   ; do not draw yet
      res_std@gsnFrame        = False  ; don't advance frame
      res_std@gsnXYFillColors = (/1.0, 0.9, 0.9/)  ; "lightpink"
      res_std@xyLineColor     = -1  ; Make lines transparent
      mmstddev = new((/2, nbins/), float)
;      mmstddev(0, :) = mmavg(0, :) - mmstd(0, :)
;      mmstddev(1, :) = mmavg(0, :) + mmstd(0, :)
      mmstddev(0, :) = mmp10(0, :)
      mmstddev(1, :) = mmp90(0, :)
      plotstd = gsn_csm_xy(wks, result_avg&bin, mmstddev, res_std)
      delete(mmstddev)
      overlay(plots(ii), plotstd)
    end if

    do i = 0, nbins - 1
      y0 = result_avg(ii, i)
      if (.not.ismissing(y0)) then
        x0 = result_avg&bin(i)
        stdbar(i, ii) = gsn_add_polyline(wks, plots(ii), (/x0, x0/), \
                                         (/y0 + result_std(ii, i), y0 - \
                                         result_std(ii, i)/), polyres)
        polyres@gsMarkerIndex = 16
        polyres@gsMarkerColor = "red"
        centers(i, ii) = gsn_add_polymarker(wks, plots(ii), (/x0, x0/), \
                                            (/y0, y0/), polyres)
        polyres@gsMarkerIndex = 4
        polyres@gsMarkerColor = "black"
        centersout(i, ii) = gsn_add_polymarker(wks, plots(ii), (/x0, x0/), \
                                               (/y0, y0/), polyres)
      end if
      y0 = result_avg(ref_ind, i)
      if (.not.ismissing(y0)) then
        x0 = result_avg&bin(i)
        stdbarR(i, ii) = gsn_add_polyline(wks, plots(ii), (/x0, x0/), \
                                          (/y0 + result_std(ref_ind, i), \
                                          y0 - result_std(ref_ind, i)/), \
                                          polyresRef)
        polyresRef@gsMarkerIndex = 16
        polyresRef@gsMarkerColor = "white"
        centersR(i, ii) = gsn_add_polymarker(wks, plots(ii), (/x0, x0/), \
                                             (/y0, y0/), polyresRef)
        polyresRef@gsMarkerIndex = 4
        polyresRef@gsMarkerColor = "black"
        centersRout(i, ii) = gsn_add_polymarker(wks, plots(ii), (/x0, x0/), \
                                                (/y0, y0/), polyresRef)
      end if
    end do

    draw(plots(ii))
    frame(wks)

  end do

  pres                      = True    ; needed to override
                                      ; panelling defaults
  pres@gsnPanelCenter                   = False

  idx0 = ind(.not.ismissing(plots))
  n = dimsizes(idx0)
  pres@gsnPanelFigureStrings = names(idx0)
  pres@gsnPanelFigureStringsFontHeightF = min((/0.008, 0.008 * 6.0 \
                                              / tofloat((dim_MOD + 1) / 2)/))
  pres@lbLabelFontHeightF               = min((/0.01, 0.01 * 6.0 \
                                              / tofloat((dim_MOD + 1) / 2)/))
  outfile = panelling(wks, plots(idx0), (n + 3) / 4, 4, pres)
  delete(idx0)

  log_info("Wrote " + outfile)

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

  ; ----------------------------------------------------------------------
  ; write mmm and ref to additional netcdf
  ; ----------------------------------------------------------------------

  if ((mm_ind .ge. 0) .and. (ref_ind .ge. 0)) then
    mmm = result_avg(mm_ind, :)
    ref = result_avg(ref_ind, :)

    mmm@var = var_y + "_mmm"
    ref@var = var_y + "_ref"

    ratio = mmm
    ratio = ratio / ref
    ratio@average = avg(ratio)
    ratio@var = var_y + "_ratio"

    nc_filename2 = work_dir + "clouds_scatter_" + var_x + "_" + var_y + \
      filename_add + "_ref_mmm_ratio.nc"

    nc_outfile2 = ncdf_write(mmm, nc_filename2)
    nc_filename2@existing = "append"
    nc_outfile2 = ncdf_write(ref, nc_filename2)
    nc_outfile2 = ncdf_write(ratio, nc_filename2)
  end if

  ; ----------------------------------------------------------------------

  leave_msg(DIAG_SCRIPT, "")

end
