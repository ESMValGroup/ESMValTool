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
;   clevels: explicit values for probability labelbar (array)
;   filename_add: optionally add this string to plot filesnames
;   nbins: number of equally spaced bins (var_x), default = 100
;   sidepanels: show/hide side panels
;   xlabel: label overriding variable name for x-axis (e.g. SST)
;   ylabel: label overriding variable name for y-axis (e.g. omega500)
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
;   20230117-lauer_axel: added support for ICON (code from Manuel)
;   20220126-lauer_axel: added optional variable labels for x- and y-axes
;   20211118-lauer_axel: added output of frequency distributions
;   20210607-lauer_axel: added multi-model-average (= average over all models)
;                        as an individual plot
;   20210408-lauer_axel: written
;
; ############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/dataset_selection.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  diag = "clouds_dyn_matrix.ncl"
  variables = get_unique_values(metadata_att_as_array(variable_info, \
                                                      "short_name"))

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

  calcicefrac = new(3, logical)
  calctcwp = new(3, logical)
  calcswclt = new(3, logical)
  calclwclt = new(3, logical)

  calcicefrac = False
  calctcwp = False
  calcswclt = False
  calclwclt = False

  ; if present, replace special variables "columnicefrac" and "totalcwp" in
  ; list with variables used for variable derivation

  tmplist = (/var_x, var_y, var_z/)
  checklist = (/"columnicefrac", "totalcwp", "swcreclt", "lwcreclt"/)
  substlist = (/(/"clivi", "lwp"/), (/"clivi", "lwp"/), (/"swcre", "clt"/), \
                (/"lwcre", "clt"/)/)
  varstart = new(3, integer)
  varstart(0) = 0

  do i = 0, dimsizes(tmplist) - 1
    if (i .eq. 0) then
      if (any(checklist .eq. tmplist(i))) then
        j = ind(checklist .eq. tmplist(i))
        varlist = substlist(j, :)
        varstart(1) = 2
      else
        varlist = tmplist(i)
        varstart(1) = 1
      end if
    else
      if (any(checklist .eq. tmplist(i))) then
        j = ind(checklist .eq. tmplist(i))
        varlist := array_append_record(varlist, substlist(j, :), 0)
        if (i .lt. dimsizes(tmplist) - 1) then
          varstart(i + 1) = varstart(i) + 2
        end if
      else
        varlist := array_append_record(varlist, tmplist(i), 0)
        if (i .lt. dimsizes(tmplist) - 1) then
          varstart(i + 1) = varstart(i) + 1
        end if
      end if
    end if
  end do

  do i = 0, dimsizes(tmplist) - 1
    if (tmplist(i) .eq. checklist(0)) then
      calcicefrac(i) = True
    else if (tmplist(i) .eq. checklist(1)) then
      calctcwp(i) = True
    else if (tmplist(i) .eq. checklist(2)) then
      calcswclt(i) = True
    else if (tmplist(i) .eq. checklist(3)) then
      calclwclt(i) = True
    end if
    end if
    end if
    end if
  end do

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
  set_default_att(diag_script_info, "nbins", 100)
  set_default_att(diag_script_info, "sidepanels", False)
  set_default_att(diag_script_info, "xlabel", var_x)
  set_default_att(diag_script_info, "ylabel", var_y)

  if (diag_script_info@filename_add .ne. "") then
    filename_add = "_" + diag_script_info@filename_add
  else
    filename_add = ""
  end if

  if (diag_script_info@sidepanels) then
    flag_sidepanels = True
  else
    flag_sidepanels = False
  end if

  nbins = toint(diag_script_info@nbins)

  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir
  ; Create work dir
  system("mkdir -p " + work_dir)

end

begin
  ; ############
  ; # get data #
  ; ############

  info_x = select_metadata_by_name(input_file_info, varlist(varstart(0)))
  names_x = metadata_att_as_array(info_x, "dataset")
  projects_x = metadata_att_as_array(info_x, "project")
  info_y = select_metadata_by_name(input_file_info, varlist(varstart(1)))
  names_y = metadata_att_as_array(info_y, "dataset")
  projects_y = metadata_att_as_array(info_y, "project")
  info_z = select_metadata_by_name(input_file_info, varlist(varstart(2)))
  names_z = metadata_att_as_array(info_z, "dataset")
  projects_z = metadata_att_as_array(info_z, "project")

  refidx_x = ind(names_x .eq. refname(varstart(0)))
  refidx_y = ind(names_y .eq. refname(varstart(1)))
  refidx_z = ind(names_z .eq. refname(varstart(2)))

  if (ismissing(refidx_x) .or. ismissing(refidx_y) .or. ismissing(refidx_z)) \
    then
    refidx_x = -1
    refidx_y = -1
    refidx_z = -1
  end if

  ref_ind = refidx_x
  names = names_x
  projects = projects_x

  if (ref_ind .ge. 0) then
    ; if reference datasets for var_x, var_y, var_z are from different
    ; sources
    uninames = get_unique_values(refname)
    names(ref_ind) = str_join(uninames, "/")
    delete(uninames)
  end if

  ; find all indices of models w/o MultiModelMean/MultiModelMedian (if present)

  idxmod = get_mod(names_x, projects_x)

  if (idxmod(0) .eq. -1) then  ; no model found
    flag_multimod = False
  elseif (dimsizes(idxmod) .eq. 1) then  ; one model found
    flag_multimod = False
  else  ; more than one model found
    flag_multimod = True
  end if

  result = new((/dim_MOD, nbins, nbins/), float)
  count = new((/dim_MOD, nbins, nbins/), float)
  bincenter_x = new((/nbins/), float)
  bincenter_y = new((/nbins/), float)
  bin_x0 = new((/nbins/), float)
  bin_x1 = new((/nbins/), float)
  bin_y0 = new((/nbins/), float)
  bin_y1 = new((/nbins/), float)

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
    bin_x0(n) = bincenter_x(n) - 0.5 * binsize_x
    bin_x1(n) = bincenter_x(n) + 0.5 * binsize_x
    y0 = n * binsize_y
    y1 = y0 + binsize_y
    bincenter_y(n) = ymin + 0.5 * (y0 + y1)
    bin_y0(n) = bincenter_y(n) - 0.5 * binsize_y
    bin_y1(n) = bincenter_y(n) + 0.5 * binsize_y
  end do

  atts_x = True
  atts_y = True
  atts_z = True

  do ii = 0, dim_MOD - 1
    atts_x@short_name = varlist(varstart(0))
    atts_y@short_name = varlist(varstart(1))
    atts_z@short_name = varlist(varstart(2))

    ; reference datasets may have different names
    if (ii .eq. refidx_x) then
      atts_x@dataset = refname(varstart(0))
      atts_y@dataset = refname(varstart(1))
      atts_z@dataset = refname(varstart(2))
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

    atts_list = NewList("fifo")

    ListAppend(atts_list, atts_x)
    ListAppend(atts_list, atts_y)
    ListAppend(atts_list, atts_z)

    vars_list = NewList("fifo")

    ListAppend(vars_list, x)
    ListAppend(vars_list, y)
    ListAppend(vars_list, z)

    do i = 0, 2
      ; read second variable needed to derive icefrac/tcwp/swcreclt/lwcreclt
      if (calcicefrac(i) .or. calctcwp(i) .or. calcswclt(i) .or. \
          calclwclt(i)) then
        atts_list[i]@short_name = varlist(varstart(i) + 1)
        if (ii .eq. refidx_x) then
          atts_list[i]@dataset = refname(varstart(i) + 1)
        end if
        info = select_metadata_by_atts(input_file_info, atts_list[i])
        var2 = read_data(info[0])
        delete(info)
        var2 = where(isnan_ieee(var2), var2@_FillValue, var2)
      end if

      ; calculate column ice fraction

      if (calcicefrac(i)) then
        min_mass = 1.0e-6
        ; filter invalid values (needed for some models)
        vars_list[i] = where(vars_list[i] .lt. 0.0, vars_list[i]@_FillValue, \
                             vars_list[i])
        vars_list[i] = where(isnan_ieee(vars_list[i]), \
                             vars_list[i]@_FillValue, vars_list[i])
        var2 = where(var2 .lt. 0.0, var2@_FillValue, var2)
        mass = vars_list[i] + var2
        mass = where(mass .lt. min_mass, mass@_FillValue, mass)

        ; ice fraction = ice / (ice + lwp) * 100%
        vars_list[i] = 100.0 * vars_list[i] / mass
        delete(mass)

        vars_list[i]@units = "%"
        vars_list[i]@long_name = "cloud ice fraction"
        vars_list[i]@var = "columnicefrac"
      end if

      ; calculate total cloud water path as sum of liquid water path (lwp)
      ; and ice water path (clivi);
      ; we do not use the CMOR variable clwvi directly as this variable
      ; erroneously contains only cloud liquid water for some models

      if (calctcwp(i)) then
        vars_list[i] = vars_list[i] + var2
        vars_list[i]@long_name = "Condensed Water Path"
        vars_list[i]@var = "totalcwp"
      end if

      ; calculate swcre divided by cloud fraction

      if (calcswclt(i)) then
        var2 = where(var2 .le. 1.0e-6, var2@_FillValue, var2)
        vars_list[i] = vars_list[i] / var2 * 100.0
        vars_list[i]@long_name = \
          "Shortwave cloud radiative effect / cloud cover"
        vars_list[i]@var = "swcreclt"
      end if

      ; calculate swcre divided by cloud fraction

      if (calclwclt(i)) then
        var2 = where(var2 .le. 1.0e-6, var2@_FillValue, var2)
        vars_list[i] = vars_list[i] / var2 * 100.0
        vars_list[i]@long_name = \
          "Longwave cloud radiative effect / cloud cover"
        vars_list[i]@var = "lwcreclt"
      end if

      if (isvar("var2")) then
        delete(var2)
      end if
    end do

    delete(atts_list)

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

    do i = 0, 2
      dims = getvardims(vars_list[i])
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

    delete(vars_list)

    delete(dims_x)
    delete(dims_y)
    delete(dims_z)
    delete(testidx)

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

    ; This approach to grid the data is significantly slower than
    ; using "bin_avg" but still fine for moderate grids (e.g. 100x100 grid
    ; cells). In NCL 6.6.2, the function bin_avg crashes after being
    ; called several times (segmentation violation). This is the reason
    ; for choosing this "manual" approach.

    do n = 0, nbins - 1
      selidx = ind((x1d .ge. bin_x0(n)) .and. (x1d .lt. bin_x1(n)))
      if (all(ismissing(selidx))) then
        delete(selidx)
        result(ii, :, n) = result@_FillValue
        count(ii, :, n) = count@_FillValue
        continue
      end if
      xsel = x1d(selidx)
      ysel = y1d(selidx)
      zsel = z1d(selidx)
      delete(selidx)
      do m = 0, nbins - 1
        selidx = ind((ysel .ge. bin_y0(m)) .and. (ysel .lt. bin_y1(m)))
        if (.not.all(ismissing(selidx))) then
          result(ii, m, n) = avg(zsel(selidx))
          count(ii, m, n) = num(.not.ismissing(zsel(selidx)))
        else
          result(ii, m, n) = result@_FillValue
          count(ii, m, n) = count@_FillValue
        end if
        delete(selidx)
      end do
      delete(xsel)
      delete(ysel)
      delete(zsel)
    end do

;    r = bin_avg(x1d, y1d, z1d, bincenter_x, bincenter_y, False)
;    result(ii, :, :) = r(0, :, :)

    delete(x1d)
    delete(y1d)
    delete(z1d)
;    delete(r)

    count(ii, :, :) = count(ii, :, :) / sum(count(ii, :, :)) * 1.0e2

  end do  ; ii-loop (models)

  ; ###########################################
  ; # netCDF output                           #
  ; ###########################################

  nc_filename = work_dir + "clouds_scatter_" + var_x + "_" + var_y + "_" + \
    var_z + filename_add + ".nc"
  nc_filename2 = work_dir + "clouds_scatter_prob_" + var_x + "_" + var_y + \
    filename_add + ".nc"

  result!0 = "model"
  result!1 = "bin_y"
  result!2 = "bin_x"

  result&model = str_sub_str(names, "/", "-")
  result&bin_y = bincenter_y
  result&bin_x = bincenter_x

  result@diag_script = (/DIAG_SCRIPT/)
  result@var = var_z
  result@var_long_name = long_name
  result@var_units = zunits
  result@_FillValue = 1.0e20

  copy_VarCoords(result, count)
  count@diag_script = (/DIAG_SCRIPT/)
  count@var = "count"
  count@var_units = "1e-3 %"

  nc_outfile = ncdf_write(result, nc_filename)
  nc_outfile2 = ncdf_write(count, nc_filename2)

  ; ###########################################
  ; # create the plots                        #
  ; ###########################################

  nplots = dim_MOD
  if (flag_multimod) then
    nplots = nplots + 1
  end if

  plots = new((/nplots, 2/), graphic)
  baseplots = new(nplots, graphic)
  xaddplots = new(nplots, graphic)
  yaddplots = new(nplots, graphic)

  dplots = new(nplots, graphic)

  plots_c = new((/nplots, 2/), graphic)
  baseplots_c = new(nplots, graphic)
  xaddplots_c = new(nplots, graphic)
  yaddplots_c = new(nplots, graphic)

  res = True
  res@cnFillOn       = True               ; color Fill
  res@cnFillMode     = "RasterFill"       ; Raster Mode
  res@cnLinesOn      =  False             ; Turn off contour lines
  res@tiXAxisString  = diag_script_info@xlabel + " (" + xunits + ")"
  res@tiYAxisString  = diag_script_info@ylabel + " (" + yunits + ")"
  res@lbOrientation  = "Vertical"
  res@lbTitleString      = zunits
  res@lbTitlePosition    = "Right"
  res@lbTitleFontHeightF = 0.0275
  res@lbLabelFontHeightF = 0.0275
  res@tmXMajorGrid = True
  res@tmXMinorGrid = True
  res@tmYMajorGrid = True
  res@tmYMinorGrid = True

  dres = res

  if (flag_sidepanels) then
    res@pmLabelBarOrthogonalPosF = 0.35
  end if

  cres = res
  cres@lbTitleString = count@var_units

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

  if (isatt(diag_script_info, "clevels")) then
    cres@cnLevelSelectionMode = "ExplicitLevels"
    cres@cnLevels = diag_script_info@clevels
  end if

  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_dyn_matrix_" + \
                var_x + "_" + var_y + "_" + var_z + filename_add)
  dwks = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_dyn_matrix_bias_" + \
                 var_x + "_" + var_y + "_" + var_z + filename_add)
  cwks = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_dyn_matrix_prob_" + \
                 var_x + "_" + var_y + filename_add)

  res@gsnDraw  = False   ; don't draw yet
  res@gsnFrame = False   ; don't advance frame yet

  cres@gsnDraw  = False   ; don't draw yet
  cres@gsnFrame = False   ; don't advance frame yet

  xyres1                 = True         ; xy plot mods desired
  xyres1@tmXBMinorOn     = False        ; no minor tickmarks
  xyres1@tmXBLabelStride = 2            ; label stride
  xyres1@gsnDraw         = False        ; don't draw yet
  xyres1@gsnFrame        = False        ; don't advance frame yet
  xyres1@xyLineThicknessF = 2.0
  xyres1@tmEqualizeXYSizes = True
  xyres1@tmXBLabelFontHeightF = 0.0275

  xyres2 = xyres1

  xyres1@vpHeightF = .20          ; set width of second plot
  xyres1@trXMinF = xmin
  xyres1@trXMaxF = xmax
  xyres1@tiXAxisString = diag_script_info@xlabel + " (" + xunits + ")"

  xyres2@vpWidthF = .20          ; set width of second plot
  xyres2@trYMinF = ymin
  xyres2@trYMaxF = ymax
  xyres2@tiXAxisSide = "Top"
  xyres2@tmXTLabelsOn = True
  xyres2@tmXBLabelsOn = False

  do ii = 0, nplots - 1
    if (ii .lt. dim_MOD) then
      plotdata = result(ii, :, :)
      countdata = count(ii, :, :)
      plotname = names(ii)
    else
      plotdata = dim_avg_n_Wrap(result(idxmod, :, :), 0)
      countdata = dim_avg_n_Wrap(count(idxmod, :, :), 0)
      plotname = "Multi-model average"
    end if

    countdata = countdata * 1.0e3

    res@tiMainString = plotname
    res@tiMainFontHeightF = 0.03

    if (ii .eq. ref_ind) then
      if (refname(varstart(0)) .eq. refname(varstart(1))) then
        plotname = refname(varstart(0))
      else
        plotname = refname(varstart(0)) + "/" + refname(varstart(1))
      end if
    end if

    cres@tiMainString = plotname

    ; z-values

    baseplots(ii) = gsn_csm_contour(wks, plotdata, res)

    if (flag_sidepanels) then
      xyres1@tiYAxisString = zunits
      xyres2@tiXAxisString = zunits

      xaddplots(ii) = gsn_csm_xy(wks, result&bin_x, \
                                 dim_avg_n(plotdata, 0), xyres1)
      yaddplots(ii) = gsn_csm_xy(wks, dim_avg_n(plotdata, 1), \
                                 result&bin_y, xyres2)
      xyres1@gsnAttachPlotsXAxis = True
      plots(ii, 0) = gsn_attach_plots(baseplots(ii), xaddplots(ii), \
                                      res, xyres1)
      delete(xyres1@gsnAttachPlotsXAxis)
      plots(ii, 1) = gsn_attach_plots(baseplots(ii), yaddplots(ii), \
                                      res, xyres2)
    end if

    draw(baseplots(ii))
    frame(wks)

    ; differences z-values

    if (isdefined("dresplot")) then
      delete(dresplot)
    end if
    if (ii .eq. ref_ind) then
      dresplot = res
      dresplot@tiMainString =  "REF"
      diff = plotdata
    else
      dresplot = dres
      dresplot@tiMainString =  plotname + " - REF"
      diff = plotdata - result(ref_ind, :, :)
    end if
    copy_VarCoords(plotdata, diff)
    dplots(ii) = gsn_csm_contour(dwks, diff, dresplot)

    ; probability values (%)

    xyres1@tiYAxisString = count@var_units
    xyres2@tiXAxisString = count@var_units

    baseplots_c(ii) = gsn_csm_contour(cwks, countdata, cres)

    if (flag_sidepanels) then
      xaddplots_c(ii) = gsn_csm_xy(cwks, count&bin_x, \
                                   dim_avg_n(countdata, 0), xyres1)
      yaddplots_c(ii) = gsn_csm_xy(cwks, dim_avg_n(countdata, 1), \
                                   count&bin_y, xyres2)

      xyres1@gsnAttachPlotsXAxis = True
      plots_c(ii, 0) = gsn_attach_plots(baseplots_c(ii), \
                                        xaddplots_c(ii), cres, xyres1)
      delete(xyres1@gsnAttachPlotsXAxis)
      plots_c(ii, 1) = gsn_attach_plots(baseplots_c(ii), \
                                        yaddplots_c(ii), cres, xyres2)
    end if

    draw(baseplots_c(ii))
    frame(cwks)

    delete(plotdata)
    delete(countdata)
  end do

;  pres                      = True    ; needed to override
;                                      ; panelling defaults
;  pres@gsnPanelCenter                   = False
;
;  pres@gsnPanelFigureStrings = names
;  pres@gsnPanelFigureStringsFontHeightF = min((/0.008, 0.008 * 6.0 \
;                                              / tofloat((dim_MOD + 1) / 2)/))
;  pres@lbLabelFontHeightF               = min((/0.01, 0.01 * 6.0 \
;                                              / tofloat((dim_MOD + 1) / 2)/))
;  outfile = panelling(wks, baseplots, (dim_MOD + 3) / 4, 4, pres)
;  doutfile = panelling(dwks, dplots, (dim_MOD + 3) / 4, 4, pres)

  outfile = wks@fullname
  doutfile = dwks@fullname
  coutfile = cwks@fullname

  log_info("Wrote " + outfile)
  log_info("Wrote " + doutfile)
  log_info("Wrote " + coutfile)

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
  ; write relative differences of the probabilities (count) to netcdf
  ; ----------------------------------------------------------------------

  ref = where(count(ref_ind, :, :) .gt. 1.0e-3, count(ref_ind, :, :), \
              count@_FillValue)
  if (dimsizes(idxmod) .gt. 1) then
    modavg = dim_avg_n_Wrap(count(idxmod, :, :), 0)
  else
    modavg = count(idxmod, :, :)
  end if
  modavg = where(modavg .gt. 1.0e-3, modavg, modavg@_FillValue)

  diff = modavg
  diff = (diff - ref) / ref * 100.0

  nc_filename3 = work_dir + "clouds_scatter_prob_" + var_x + "_" + var_y + \
    "_mmm-ref-reldiff" + filename_add + ".nc"

  diff@var_long_name = "(count(MMM) - count(REF)) / count(REF) * 100%"
  diff@var_units = "%"

  nc_outfile3 = ncdf_write(diff, nc_filename3)

  leave_msg(DIAG_SCRIPT, "")

end
