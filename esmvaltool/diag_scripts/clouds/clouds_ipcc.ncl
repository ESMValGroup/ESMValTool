; CLOUDS_IPCC
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; PROJECT-NAME EMBRACE
; ############################################################################
; Description
;   Calculates the multi-model mean bias of annual mean 2-d cloud variables
;   compared with a reference data set (observations). In addition,
;   zonal averages of the individual models, the multi-model mean and the
;   reference data set (observations) are calculated.
;   If more than one variable is specified in the namelist, variables 2 to n
;   are assumed to be observational error estimates. These error estimates are
;   added to the zonal mean plot as light red shading.
;
; Required diag_script_info attributes (diagnostic specific)
;   none
;
; Optional diag_script_info attributes (diagnostic specific)
;   explicit_cn_levels: contour levels
;   mask_ts_sea_ice:    - True = mask T < 272 K as sea ice (only for
;                         variable "ts")
;                       - False = no additional grid cells masked for
;                         variable "ts"
;   projection:         map projection, e.g., Mollweide, Mercator
;   styleset:           style set for zonal mean plot ("CMIP5", "DEFAULT")
;   timemean:           time averaging, i.e. "seasonalclim" (DJF, MAM, JJA,
;                       SON), "annualclim" (annual mean)
;   valid_fraction:     used for creating sea ice mask
;                       (mask_ts_sea_ice = true): fraction of valid time steps
;                       required to mask grid cell as valid data
;
; Required variable_info attributes (variable specific)
;   reference_dataset:  name of reference data set
;
; Optional variable_info attributes (variable specific)
;   long_name:          description of variable
;   units:              variable units
;
; Caveats
;   KNOWN ISSUES
;     1) specifying more than one data set for the observational uncertainties
;        may lead to unexpected or undefined results
;     2) Bias and zonal means cannot be written to the same netCDF because
;        function ncdf_write does not support writing variables with different
;        dimensions to the same output file yet. As the Python function for
;        writing the provenance information called by function log_provenance
;        does not support different captions for netCDF and (associated)
;        plotfile, the caption does not represent the contents of the netCDF
;        exactly but rather those of the plot. The data used for creating the
;        plot is written to 2 separate netCDFs.
;
; Modification history
;   20190222-lauer_axel: added output of provenance (v2.0)
;   20181119-lauer_axel: adapted code to multi-variable capable framework
;   20180923-lauer_axel: added writing of results to netcdf
;   20180529-lauer_axel: code rewritten for ESMValTool v2.0
;   20170620-lauer_axel: added tags for reporting
;   20160920-lauer_axel: added optional shading of observational uncertainties
;                        to the zonal mean plot
;   20160901-lauer_axel: added regridding option 1 deg x 1 deg
;   20151027-lauer_axel: moved call to 'write_references' to the beginning
;                        of the code
;   20150428-lauer_axel: written.
;
; ############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/contour_maps.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  tmp = metadata_att_as_array(input_file_info, "short_name")
  variables = get_unique_values(tmp)
  delete(tmp)
  n = dimsizes(variables)

  if (n.gt.2) then
    log_info(DIAG_SCRIPT + ": WARNING - max 2 variables supported. " \
             + "Using only first two variables.")
  end if

  mainvarind = 0
  errvarind = -1

  if (n .gt. 1) then
    do i = 0, 1
      if (isStrSubset(variables(i), "err")) then
        errvarind = i
        exit
      end if
    end do
    if (errvarind .eq. 0) then
      mainvarind = 1
    end if
  end if

  var0 = variable_info[mainvarind]@short_name
  info0 = select_metadata_by_name(input_file_info, var0)
  dim_MOD0 = ListCount(info0)
  if (isatt(variable_info[mainvarind], "reference_dataset")) then
    refname = variable_info[mainvarind]@reference_dataset
  end if
  if (isatt(variable_info[mainvarind], "alternative_dataset")) then
    refname2 = variable_info[mainvarind]@alternative_dataset
  end if
  names = metadata_att_as_array(info0, "dataset")
  infiles = metadata_att_as_array(info0, "filename")

  if (errvarind .ge. 0) then
    var1 = variable_info[errvarind]@short_name
    info1 = select_metadata_by_name(input_file_info, var1)
    dim_MOD1 = ListCount(info1)
    names_err = metadata_att_as_array(info1, "dataset")
    infiles_err = metadata_att_as_array(info1, "filename")
  else
    dim_MOD1 = 0
  end if

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Set default values for non-required diag_script_info attributes
  set_default_att(diag_script_info, "mask_ts_sea_ice", False)
  set_default_att(diag_script_info, "projection", "CylindricalEquidistant")
  set_default_att(diag_script_info, "timemean", "annualclim")
  set_default_att(diag_script_info, "valid_fraction", 0.5)

  ; time averaging: at the moment, only "annualclim" and "seasonalclim"
  ; are supported

  numseas = 1          ; default
  season = (/"annual"/)

  timemean = diag_script_info@timemean

  if (timemean.eq."seasonalclim") then
    delete(season)
    numseas = 4
    season = (/"DJF", "MAM", "JJA", "SON"/)
  end if

  ; create string for caption (netcdf provenance)

  allseas = season(0)
  do is = 1, numseas - 1
    allseas = allseas + "/" + season(i)
  end do

  ; make sure path for (mandatory) netcdf output exists

  write_nc = True
  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

  if (diag_script_info@write_plots.eq."True") then
    write_plots = True
  else
    write_plots = False
  end if

end

begin
  ; ========================================================================
  ; ========================== initialization ==============================
  ; ========================================================================

  diag_script_base = basename(DIAG_SCRIPT)

  ; check for reference dataset definition
  if (.not.(isvar("refname"))) then
    error_msg("f", DIAG_SCRIPT, "", "no reference dataset defined in recipe")
  end if

  ; set reference dataset

  ref_ind = ind(names .eq. refname)
  if (ismissing(ref_ind)) then
    error_msg("f", DIAG_SCRIPT, "", "reference dataset (" \
              + refname + ") is missing")
  end if

  if (isvar("refname2")) then
    ref_ind2 = ind(names .eq. refname2)
    if (ismissing(ref_ind2)) then
      ref_ind2 = -1
    end if
  else
    ref_ind2 = -1
  end if

  ; get multi-model mean index

  mm_ind = ind(names .eq. "MultiModelMean")

  if (ismissing(mm_ind)) then
    error_msg("f", DIAG_SCRIPT, "", "multi-model mean is missing (required)")
  end if

  mask_ts_sea_ice = diag_script_info@mask_ts_sea_ice

  if ((var0 .eq. "ts") .and. (mask_ts_sea_ice))
    ; Create a missing value mask (from reference data)
    ; The dim_avg_wgt function with option limit is used on the time
    ; coordinate: in this way each grid box containing more than
    ; valid_fraction of missing values along the time series is
    ; set to missing.

    data_temp = read_data(info0[ref_ind])

    ww = data_temp&time
    ww = 1.
    limit  = toint(diag_script_info@valid_fraction * \
                   dimsizes(data_temp&time))

    if (limit.lt.1) then
      limit = 1
    end if

    ; temperatures below 272 K are sea ice ---> mask
    data_temp = where(data_temp.lt.272.0, data_temp@_FillValue, data_temp)

    ; dim_avg_wgt_n_Wrap produces an invalid result for rank = 1
    if (dimsizes(dimsizes(data_temp)).eq.1) then
      tmp = dim_avg_wgt_Wrap(data_temp, ww, limit)
    else
      tmp = dim_avg_wgt_n_Wrap(data_temp, ww, limit, 0)
    end if

    global_mask = where(ismissing(tmp), tmp, 0.)

    delete(tmp)
    delete(ww)
    delete(data_temp)
  end if

  ; process uncertainty estimate (if present)

  if (errvarind .ge. 0) then
    do imod = 0, dim_MOD1 - 1
      log_info("processing error estimate from " + names_err(imod))

      A0 = read_data(info1[imod])

      log_info(" ++++++++++++++ Treating variable " + \
               var1 + " as error variable when averaging ")
      A0 = A0 * A0

      ; time mean

      data1 = time_operations(A0, -1, -1, "average", timemean, True)
      delete(A0)

      ; zonal mean

      err = dim_avg_n_Wrap(data1, dimsizes(dimsizes(data1)) - 1)
      ; err is still the squared error (squared for averaging)!
      err = sqrt(err)
      delete(data1)

      if (.not.isdefined("err_zm")) then
        tmp = dimsizes(err)
        dims = array_append_record(dim_MOD0, tmp, 0)
        delete(tmp)
        err_zm = new(dims, float)
        err_zm!0 = "model"
        err_zm&model = names
        if (numseas.eq.1) then
          err_zm!1 = "lat"
          err_zm&lat = err&lat
        else
          err_zm!1 = "season"
          err_zm&season = err&season
          err_zm!2 = "lat"
          err_zm&lat = err&lat
        end if
        delete(dims)
      end if

      idx = ind(names .eq. names_err(imod))
      if (all(ismissing(idx)) .or. (dimsizes(idx) .gt. 1)) then
        error_msg("f", DIAG_SCRIPT, "", "error estimates for dataset " + \
                  names_err(imod) + " do not have a (unique) corresponding " \
                  + "dataset for variable " + var0)
      end if
      ; err_zm&model(imod) = names_err(imod)
      if (numseas.gt.1) then
        err_zm(idx, :, :) = err(:, :)
      else
        err_zm(idx, :) = err(:)
      end if

      delete(err)

    end do  ; imod-loop
  end if  ; if error estimate is present

  ; ========================================================================
  ; =========================== calculations ===============================
  ; ========================================================================

  ; debug output

;  system ("rm debug.nc")
;  debugfile = addfile("debug.nc","c")
;  debugfile->mask = global_mask

  ; calculate zonal means (main variable)

  do imod = 0, dim_MOD0 - 1
    ; note: 1) masking is handled by the backend
    ;       2) multi-model mean is calculated by the backend

    A0 = read_data(info0[imod])
    data = time_operations(A0, -1, -1, "average", timemean, True)
    delete(A0)

    if (isdefined("global_mask")) then
      if (numseas.eq.1) then
        data = data + global_mask
      else
        do is = 0, numseas - 1
          data(is, :, :) = data(is, :, :) + global_mask
        end do
      end if
    end if

    ; zonal mean

    tmp = dim_avg_n_Wrap(data, dimsizes(dimsizes(data)) - 1)

    if (.not.isdefined("zm")) then
      tmpdim = dimsizes(tmp)
      dims = array_append_record(dim_MOD0, tmpdim, 0)
      delete(tmpdim)
      zm = new(dims, float)
      zm!0 = "model"
      zm&model = names
      if (numseas.eq.1) then
        zm!1 = "lat"
        zm&lat = tmp&lat
      else
        zm!1 = "season"
        zm&season = tmp&season
        zm!2 = "lat"
        zm&lat = tmp&lat
      end if
      delete(dims)
    end if

    if (numseas.gt.1) then
      zm(imod, :, :) = tmp(:, :)
    else
      zm(imod, :) = tmp(:)
    end if

    delete(tmp)

    ; save maps of multi-model mean and reference data

    if (imod.eq.mm_ind) then
      mmdata = data
      copy_VarMeta(data, mmdata)
    end if
    if (imod.eq.ref_ind) then
      refdata = data
      copy_VarMeta(data, refdata)
    end if

    ; debug output

;    debugfile->$names(imod)$ = data1

    delete(data)

  end do  ; imod

  ; differences between multi-model mean and reference data set

  diff = mmdata - refdata
  copy_VarMeta(refdata, diff)

;  debugfile->diff = diff

  ; we order the zonal mean array in a way so that
  ; the lines for the multi-model mean and reference model will
  ; be drawn on top of the lines for the individual models, i.e.:
  ;  (1) individual model(s)
  ;  (2) reference model(s) (observations)
  ;  (3) multi-model mean

  dims = dimsizes(zm)
  zonalmean = new(dims, float)
  copy_VarMeta(zm, zonalmean)

  ; model indices with no reference model(s) and no multi-model mean

  model_ind = ispan(0, dim_MOD0 - 1, 1)
  model_ind(ref_ind) = -1
  if (ref_ind2 .ge. 0) then
    model_ind(ref_ind2) = -1
  end if
  model_ind(mm_ind) = -1

  modelsonly_ind = ind(model_ind.ge.0)
  delete(model_ind)

  n = dimsizes(modelsonly_ind) - 1

  ; first entries in "zonalmean" = individual models

  if (numseas.gt.1) then
    zonalmean(0:n, :, :) = zm(modelsonly_ind, :, :)
  else
    zonalmean(0:n, :) = zm(modelsonly_ind, :)
  end if

  zonalmean&model(0:n) = zm&model(modelsonly_ind)

  ; observation(s)

  n = n + 1
  zonalmean&model(n) = zm&model(ref_ind)

  if (numseas.gt.1) then
    zonalmean(n, :, :) = zm(ref_ind, :, :)
    if (ref_ind2 .ge. 0) then
      n = n + 1
      zonalmean(n, :, :) = zm(ref_ind2, :, :)
    end if
  else
    zonalmean(n, :) = zm(ref_ind, :)
    if (ref_ind2 .ge. 0) then
      n = n + 1
      zonalmean(n, :) = zm(ref_ind2, :)
    end if
  end if

  if (ref_ind2 .ge. 0) then
    zonalmean&model(n) = zm&model(ref_ind2)
  end if

  ; last entry in "zonalmean" = multi-model mean

  n = n + 1

  if (numseas.gt.1) then
    zonalmean(n, :, :) = zm(mm_ind, :, :)  ; multi-model mean
  else
    zonalmean(n, :) = zm(mm_ind, :)  ; multi-model mean
  end if

  zonalmean&model(n) = zm&model(mm_ind)

  ; update indices for array "zonalmean"

  mm_ind = n

  delete(zm)

  ; ========================================================================
  ; ============================= plotting =================================
  ; ========================================================================

  climofiles = new(dim_MOD0 + dim_MOD1, string)
  climofiles(0:dim_MOD0 - 1) = infiles
  if (dim_MOD1 .gt. 0) then
    climofiles(dim_MOD0:dim_MOD0 + dim_MOD1 - 1) = infiles_err
  end if

  if (.not. isatt(diag_script_info, "explicit_cn_levels")) then
    log_info(DIAG_SCRIPT + " (var: " + var0 + "):")
    log_info("info: using default contour levels")
    cnLevels = fspan(min(diff), max(diff), 20)
  end if

  diff@diag_script = DIAG_SCRIPT
  diff@res = True

  diff@res_gsnMaximize     = True      ; use full page for the plot
  diff@res_cnFillOn        = True      ; color plot desired
  diff@res_cnLineLabelsOn  = False     ; contour lines
  diff@res_cnLinesOn       = False
  diff@res_tiMainOn        = False
  diff@res_mpPerimOn       = False
  diff@res_mpGridLineColor = -1
  diff@res_mpGridAndLimbOn = True

  diff@res_gsnLeftStringFontHeightF = 0.015
  diff@res_cnLevelSelectionMode     = "ExplicitLevels"
  diff@res_mpOutlineOn     = True
  diff@res_cnLevels        = cnLevels
  diff@res_mpFillOn        = False

  diff@res_lbLabelBarOn          = True
  diff@res_gsnRightString        = ""
  diff@res_mpFillDrawOrder       = "PostDraw"    ; draw map fill last
  diff@res_cnMissingValFillColor = "Gray"
  diff@res_tmYLLabelsOn          = False
  diff@res_tmYLOn                = False
  diff@res_tmYRLabelsOn          = False
  diff@res_tmYROn                = False
  diff@res_tmXBLabelsOn          = False
  diff@res_tmXBOn                = False
  diff@res_tmXTLabelsOn          = False
  diff@res_tmXTOn                = False
  diff@res_cnInfoLabelOn         = False    ; turn off cn info label
  diff@res_mpProjection          = diag_script_info@projection
  if (isatt(diag_script_info, "explicit_cn_levels")) then
    if (isatt(diff, "res_cnLevels")) then
      delete(diff@res_cnLevels)
    end if
    diff@res_cnLevelSelectionMode = "ExplicitLevels"
    diff@res_cnLevels = diag_script_info@explicit_cn_levels
  end if

  diff@var = var0  ; overwrite existing entry
  if (.not. isatt(variable_info[mainvarind], "long_name")) then
    variable_info[mainvarind]@long_name = ""
  end if
  if (.not. isatt(variable_info[mainvarind], "units")) then
    variable_info[mainvarind]@units = ""
  end if

  diff@var_long_name = variable_info[mainvarind]@long_name
  diff@var_units     = variable_info[mainvarind]@units

  if (var0.eq."clt") then
    if (isatt(diff, "res_cnLevels")) then
      delete(diff@res_cnLevels)
    end if
    diff@res_cnLevels = fspan(-25, 25, 11)
  end if

  if ((var0.eq."lwp").or.(var0.eq."clivi")) then
    if (isatt(diff, "res_cnLevels")) then
      delete(diff@res_cnLevels)
    end if
    diff@res_cnLevels = fspan(-100, 100, 11)
    diff = diff * 1000
    zonalmean = zonalmean * 1000
    diff@var_units = "g m-2"
    variable_info[mainvarind]@units = diff@var_units
    zonalmean@units = diff@var_units
  end if

  if ((var0.eq."pr").or.(var0.eq."pr-mmday")) then
    diff@res_mpLandFillColor      = "Black"
    pal = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                             + "ipcc-precip-delta.rgb")
    diff@res_cnFillColors         = pal
  end if

  plots = new((/2, numseas/), graphic)
  plotfile = new(numseas, string)
  plotfile(:) = ""

  do is = 0, numseas - 1

    ; --------------------------------------------------------------------
    ; create workspace

    if (isvar("wks")) then
      delete(wks)
    end if

    wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_ipcc_" + var0 \
                  + "_" + season(is))
;    drawNDCGrid(wks) ; debugging option

    ; --------------------------------------------------------------------
    ; plot contour map

    diff@res_gsnDraw       = False  ; do not draw yet
    diff@res_gsnFrame      = False  ; don't advance frame

    diff@res_lbTitleString = "~F8~D~F21~" + diff@var_long_name + \
                             " (" + diff@var_units + ")"
    diff@res_lbTitlePosition = "Bottom"
    diff@res_lbLabelFontHeightF = 0.015
    diff@res_lbAutoManage    = False
;    diff@res_lbTopMarginF    = 0.1
    diff@res_lbTitleFontHeightF = 0.015

    ; plot contour map

    if (numseas.gt.1) then
      plots(0, is) = contour_map(wks, diff(is, :, :), var0)
    else
      plots(0, 0) = contour_map(wks, diff, var0)
    end if

    ; --------------------------------------------------------------------
    ; plot zonal means

    lat = (/zonalmean&lat/)

    ; the attribute "long_name" will be used as title string for the x-axis

    lat@long_name = "Latitude"
    lat@units = "degrees_north"

    ; the attribute "long_name" will be used as title string for the y-axis

    zonalmean@long_name = zonalmean@long_name + " (" + \
      zonalmean@units + ")"

    linethickness       = new(dim_MOD0, float)
    linecolor           = new(dim_MOD0, string)
    linedash            = new(dim_MOD0, integer)

    n = dimsizes(modelsonly_ind) - 1

    ; settings for all models that have been used to calculate the
    ; multi-model mean (= first entries in "zonalmean")

    linethickness(0:n)   = 1.0
    linecolor(0:n)       = "(/0.5, 0.5, 0.5/)"
    linedash             = 0

    ; for all other models, styleset defined in config file is
    ; used (if defined), otherwise use default values

    if (isatt(diag_script_info, "styleset")) then
      colortab     = project_style(info0, diag_script_info, "colors")
      dashtab      = project_style(info0, diag_script_info, "dashes")
      thicknesstab = project_style(info0, diag_script_info, "thicks")
      do i = n + 1, dim_MOD0 - 1
        if (isdefined("idx")) then
          delete(idx)
        end if
        idx = ind(names .eq. zonalmean&model(i))
        if (all(ismissing(idx))) then
        else
          linecolor(i) = colortab(idx(0))
          linethickness(i) = thicknesstab(idx(0))
          linedash(i) = dashtab(idx(0))
        end if
      end do
    else
      linethickness(n+1:dim_MOD0-1) = 4.0      ; reference dataset(s)
      linethickness(mm_ind)         = 4.0      ; multi-model mean
      linecolor(n+1:dim_MOD0-1)     = "Black"  ; reference data set
      linecolor(mm_ind)             = "Red"    ; multi-model mean
    end if

    res                   = True
    res@xyDashPattern     = linedash
    res@xyMonoLineColor   = False
    res@xyLineThicknesses = linethickness
    res@xyLineColors      = linecolor
    res@tmXBLabelFontHeightF = 0.023
    res@tmYLLabelFontHeightF = 0.023
    res@gsnDraw           = False  ; donot draw yet
    res@gsnFrame          = False  ; don't advance frame

    ; this controls the size and location of the second plot

    res@vpWidthF        = 0.6
    res@vpHeightF       = 0.475

    if (numseas.gt.1) then
      plots(1, is) = gsn_csm_xy(wks, lat, zonalmean(:, is, :), res)
      if (isdefined("err_zm")) then
        y = new((/2, dimsizes(err_zm&lat)/), float)
        res@gsnXYFillColors = "(/1.0, 0.8, 0.8/)"
        res@xyLineColors(:) = "transparent"
        do imod = 0, dim_MOD0 - 1
          erridx = ind(err_zm&model.eq.zonalmean&model(imod))
          if (.not.all(ismissing(erridx))) then
            y(0, :) = zonalmean(imod, is, :) - err_zm(erridx, is, :)
            y(1, :) = zonalmean(imod, is, :) + err_zm(erridx, is, :)
            errorshade = gsn_csm_xy(wks, err_zm&lat, y, res)
            overlay(plots(1, is), errorshade)
          end if
        end do
      end if
    else
      plots(1, 0) = gsn_csm_xy(wks, lat, zonalmean, res)

      if (isdefined("err_zm")) then
        y = new((/2, dimsizes(err_zm&lat)/), float)
        res@gsnXYFillColors = "(/1.0, 0.8, 0.8/)"
        res@xyLineColors(:) = "transparent"
        do imod = 0, dim_MOD - 1
          erridx = ind(err_zm&model.eq.zonalmean&model(imod))
          if (.not.all(ismissing(erridx))) then
            y(0, :) = zonalmean(imod, :) - err_zm(erridx, :)
            y(1, :) = zonalmean(imod, :) + err_zm(erridx, :)
            errorshade = gsn_csm_xy(wks, err_zm&lat, y, res)
            overlay(plots(1, 0), errorshade)
          end if
        end do
      end if
    end if

    ; --------------------------------------------------------------------

    if (write_plots) then
      pres                = True
      pres@gsnPanelCenter = False
      pres@gsnPanelXF     = (/0.075, 0.625/)  ; hor. pos. of sub-plots
      pres@txString       = season(is)

      outfile = panelling(wks, plots(:, is), 1, 2, pres)
      log_info("Wrote " + wks@fullname)

      plotfile(is) = wks@fullname
    end if  ; if write_plots
  end do  ; is-loop (seasons)

  ; ###########################################
  ; # Output to netCDF                        #
  ; ###########################################

  ; note: function ncdf_write currently does not support writing variables
  ;       diff and zonalmean to the same netCDF.

  nc_filename = work_dir + "clouds_ipcc_" + var0 + "_bias.nc"
  diff@var = var0 + "_bias"
  diff@diag_script = DIAG_SCRIPT
  nc_outfile_bias = ncdf_write(diff, nc_filename)

  nc_filename = work_dir + "clouds_ipcc_" + var0 + "_zonal.nc"
  nc_filename@existing = "append"
  zonalmean@var = var0 + "_zonal"
  zonalmean@diag_script = DIAG_SCRIPT
  nc_outfile_zonal = ncdf_write(zonalmean, nc_filename)
  if (isvar("err_zm")) then
    err_zm@var = var1 + "_zonal"
    nc_outfile_zonal = ncdf_write(err_zm, nc_filename)
  end if

  ; ------------------------------------------------------------------------
  ; write provenance to netcdf output and plot file(s)
  ; ------------------------------------------------------------------------

  statistics = "clim"
  domain = "global"
  plottype = (/"geo", "zonal"/)
  caption = "Multi model mean bias (left) and zonal averages (right) " \
            + "for variable " + var0 + " (" + allseas \
            + "), reference = " + names(ref_ind) + "."

  do is = 0, numseas - 1
    log_provenance(nc_outfile_bias, plotfile(is), caption, statistics, \
                   domain, plottype, "", "", climofiles)
    log_provenance(nc_outfile_zonal, plotfile(is), caption, statistics, \
                   domain, plottype, "", "", climofiles)
  end do

  leave_msg(DIAG_SCRIPT, "")

end
