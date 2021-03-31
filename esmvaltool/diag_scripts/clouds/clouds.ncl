; CLOUDS
; ############################################################################
; Author: Axel Lauer (DLR, Germany)
; PROJECT-NAME EMBRACE
; ############################################################################
; Description
;   Calculates annual/seasonal means of 2-d (cloud) parameters for comparison
;   with a reference data set. Optionally, differences to the reference data
;   set are also plotted.
;
; Required diag_script_info attributes (diagnostic specific)
;   none
;
; Optional diag_script_info attributes (diagnostic specific)
;   embracesetup:       True = 2 plots per line, False = 4 plots per line
;                       (default)
;   explicit_cn_levels: explicit contour levels (array)
;   extralegend:        plot legend(s) to extra file(s)
;   filename_add:       optionally add this string to plot filesnames
;   multiobs_uncertainty: calculate uncertainty estimates from multiple
;                       observational datasets (true, false)
;   panel_labels:       label individual panels (true, false)
;   PanelTop:           manual override for "@gnsPanelTop" used by panel
;                       plot(s)
;   projection:         map projection for plotting (default =
;                       "CylindricalEquidistant")
;   showdiff            calculate and plot differences (default = False)
;   rel_diff:           if showdiff = True, then plot relative differences (%)
;                       (default = False)
;   ref_diff_min:       lower cutoff value in case of calculating relative
;                       differences
;                       (in units of input variable)
;   region:             show only selected geographic region given as latmin,
;                       latmax, lonmin, lonmax
;   timemean:           time averaging - "seasonal" = DJF, MAM, JJA, SON),
;                                        "annualclim" = annual mean
;   treat_var_as_error: treat variable as error when averaging (true, false)
;                       true:  avg = sqrt(mean(var*var))
;                       false: avg = mean(var)
;   var:                short_name of variable to process (default = "" - use
;                       first variable in variable list)
;
; Required variable attributes (variable specific)
;   none
;
; Optional variable_info attributes (variable specific)
;   long_name:         variable description
;   reference_dataset: reference dataset; REQUIRED when calculating
;                      differences (showdiff = True)
;   units:             variable units (for labeling plot only)
;
; Caveats
;   none
;
; Modification history
;   20210325-lauer_axel: added option to estimate observational uncertainty
;                        from multiple observational datasets
;   20210318-lauer_axel: added option to speficfy variable if more than one
;                        variable is present
;   20190220-lauer_axel: added output of provenance (v2.0)
;   20181119-lauer_axel: adapted code to multi-variable capable framework
;   20180923-lauer_axel: added writing of results to netcdf
;   20180518-lauer_axel: code rewritten for ESMValTool v2.0
;   20170621-lauer_axel: reworked code to add tags for reporting
;   20160901-lauer_axel: added regridding option 1 deg x 1 deg
;   20151027-lauer_axel: moved call to 'write_references' to the beginning
;                        of the code
;   20150415-lauer_axel: written.
;
; ############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/contour_maps.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  ; Set default values for non-required diag_script_info attributes

  set_default_att(diag_script_info, "embrace_setup", False)
  set_default_att(diag_script_info, "extralegend", False)
  set_default_att(diag_script_info, "filename_add", "")
  set_default_att(diag_script_info, "multiobs_uncertainty", False)
  set_default_att(diag_script_info, "panel_labels", True)
  set_default_att(diag_script_info, "rel_diff", False)
  set_default_att(diag_script_info, "rel_diff_min", -1.0e19)
  set_default_att(diag_script_info, "showdiff", False)
  set_default_att(diag_script_info, "timemean", "annualclim")
  set_default_att(diag_script_info, "treat_var_as_error", False)
  set_default_att(diag_script_info, "var", "")

  if (diag_script_info@var .eq. "") then
    var0 = variable_info[0]@short_name
  else
    var0 = diag_script_info@var
  end if

  variables = metadata_att_as_array(variable_info, "short_name")
  varidx = ind(variables .eq. var0)
  if (ismissing(varidx)) then
    errstr = "diagnostic " + diag + " requires the following variable: var0"
    error_msg("f", DIAG_SCRIPT, "", errstr)
  end if

  info0 = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info0)

  if (isatt(variable_info[varidx], "reference_dataset")) then
    refname = variable_info[varidx]@reference_dataset
  end if
  names = metadata_att_as_array(info0, "dataset")
  projects = metadata_att_as_array(info0, "project")

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  flag_diff = diag_script_info@showdiff
  flag_rel_diff = diag_script_info@rel_diff
  flag_rel_diff_min = diag_script_info@rel_diff_min
  flag_multiobs_unc = diag_script_info@multiobs_uncertainty

  if (.not.flag_diff .and. flag_rel_diff) then
    log_info("flag_rel_diff = True has no effect until flag_diff is also " \
             + "set to True")
  end if

  if (diag_script_info@filename_add .ne. "") then
    filename_add = "_" + diag_script_info@filename_add
  else
    filename_add = ""
  end if

  embracesetup = diag_script_info@embrace_setup

  if (isatt(diag_script_info, "projection")) then
    projection = diag_script_info@projection
    perim = False
  else
    projection = "CylindricalEquidistant"
    perim = True
  end if

  ; time averaging: at the moment, only "annualclim" and "seasonalclim"
  ; are supported

  timemean = diag_script_info@timemean
  numseas = 1          ; default
  season = (/"annual"/)

  if (timemean.eq."seasonalclim") then
    numseas = 4
    delete(season)
    season = (/"DJF", "MAM", "JJA", "SON"/)
  end if

  if (flag_multiobs_unc .and. timemean .ne. "annualclim") then
    log_info("multiobs_uncertainty = True is currently supported for annual" \
             + " means only (timemean = annualclim). Setting " \
             + " multiobs_uncertainty to False.")
    flag_multiobs_unc = False
  end if

  ; find indices of all OBS and obs4mips datasets

  if (flag_multiobs_unc) then
    idxtmp = new(dimsizes(projects), integer)
    do i = 0, dimsizes(projects) - 1
      if (isStrSubset(str_lower(projects(i)), "obs"))
        idxtmp(i) = i
      else
        idxtmp(i) = -1
      end if
    end do
    idxobs = ind(idxtmp .ge. 0)
    delete(idxtmp)
    if (all(ismissing(idxobs))) then
      log_info("No OBS or obs4mips datasets found. Setting " \
               + " multiobs_uncertainty to False.")
      flag_multiobs_unc = False
    else
      refname = "REF"
      ref_ind = dimsizes(names)
      names := array_append_record(names, (/refname/), 0)
      log_info("Using the following datasets to calculate multi-" \
               + "observational average (REF): " \
               + str_join(names(idxobs), ", "))
    end if
  else
    ref_ind = -1  ; set to invalid value
    ; if attribute is present, use it so correlations can be calculated
    if (isvar("refname")) then
      ; set reference model
      ref_ind = ind(names .eq. refname)
      if (ismissing(ref_ind)) then
        log_info("warning: reference dataset (" + refname + ") not found.")
        ref_ind = -1
      end if
    end if
  end if

  ; create string for caption (netcdf provenance)

  allseas = season(0)
  do is = 1, numseas - 1
    allseas = allseas + "/" + season(i)
  end do

  panel_labels = diag_script_info@panel_labels

  treat_var_as_error = diag_script_info@treat_var_as_error

  extralegend = diag_script_info@extralegend

  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

  if (config_user_info@write_plots.eq."True") then
    write_plots = True
  else
    write_plots = False
  end if

  climofiles = metadata_att_as_array(info0, "filename")

  outfile = new(numseas, string)
  outfile(:) = ""

  if (flag_diff) then
    outfile_d = new(numseas, string)
    outfile_d(:) = ""

    ; check for reference model definition
    if (.not.isvar("refname")) then
      error_msg("f", DIAG_SCRIPT, "", \
                "no reference dataset defined in recipe")
    end if

    ; set reference model
    if (ref_ind .lt. 0) then
      error_msg("f", DIAG_SCRIPT, "", "cannot calculate differences as " \
                + "reference dataset (" + refname + ") is missing")
    end if
  end if

end

begin
  ; ###########################################
  ; # get data and average time               #
  ; ###########################################

  ; ---------------------------------------------------------
  ; if requested, calculate multi-observational mean and standard deviation

  if (flag_multiobs_unc) then
    nobs = dimsizes(idxobs)

    ; step 1: calculate multi-obs mean

    do i = 0, nobs - 1
      A0 = read_data(info0[idxobs(i)])
      ; calculate yearly averages
      mean = time_operations(A0, -1, -1, "average", "annualclim", True)
      delete(A0)
      ; if requested, extract geographical region
      if (isatt(diag_script_info, "region")) then
        region = diag_script_info@region
        mean := area_operations(mean, region(0), region(1), region(2), \
                                region(3), "extract", False)
      end if
      if (i .eq. 0) then
        ref_avg = mean
      else
        ref_avg = ref_avg + mean
      end if
    end do
    delete(mean)
    ref_avg = ref_avg / tofloat(nobs)

    ; step 2: calculate standard deviation of all obs datasets using
    ;         the multi-obs mean

    sigma2 = new(dimsizes(ref_avg), float)

    do i = 0, nobs - 1
      A0 = read_data(info0[idxobs(i)])
      ; calculate yearly averages
      ymean = time_operations(A0, -1, -1, "average", "yearly", True)
      delete(A0)
      ; if requested, extract geographical region
      if (isatt(diag_script_info, "region")) then
        region = diag_script_info@region
        ymean := area_operations(ymean, region(0), region(1), region(2), \
                                 region(3), "extract", False)
      end if
      dims = dimsizes(ymean)
      nyears = dims(0)
      delete(dims)
      sigma2 = 0.0
      do j = 0, nyears - 1  ; dimension 0 = time
        sigma2 = sigma2 + (ref_avg - ymean(j, :, :)) ^ 2
      end do
      if (nyears .gt. 1) then
        sigma2 = sigma2 / tofloat(nyears - 1)
      end if
      if (i .eq. 0) then
        ref_std = sigma2
      else
        ref_std = ref_std + sigma2
      end if
    end do
    copy_VarMeta(ref_avg, ref_std)
    delete(sigma2)
    ref_std = sqrt(ref_std / tofloat(nobs))
  end if

  ; ---------------------------------------------------------

  nplots = dim_MOD
  if (flag_multiobs_unc) then
    nplots = nplots + 1
  end if

  maps = new((/nplots, 4/), graphic)
  maps_d = new((/nplots, 4/), graphic)

  ind_all_sorted = ispan(0, nplots, 1)  ; create array

  if (ref_ind .ge. 0) then
    ind_wo_ref = ind(names .ne. refname)
    ind_all_sorted(0) = ref_ind
    n = dimsizes(names)
    ind_all_sorted(1:n - 1) = ind_wo_ref
  end if

  corr = new((/numseas/), float)
  gavg = new((/numseas/), float)
  rmsd = new((/numseas/), float)
  bias = new((/numseas/), float)

  ; filenames for netcdf output

  nc_filename_bias = work_dir + "clouds_" + var0 + "_bias.nc"
  nc_filename_bias@existing = "append"
  nc_filename_mean = work_dir + "clouds_" + var0 + "_mean.nc"
  nc_filename_mean@existing = "append"

  res = True

  do ii = 0, nplots - 1

    imod = ind_all_sorted(ii)
    log_info("processing " + names(imod) + " ***")

    if (isvar("data1")) then
      delete(data1)
    end if

    if (isvar("A0")) then
      delete(A0)
    end if

    if (imod .ne. ref_ind .or. .not.flag_multiobs_unc) then
      A0 = read_data(info0[imod])

      ; check dimensions

      dims = getvardims(A0)
      if (dimsizes(dims) .lt. 2) then
        error_msg("f", DIAG_SCRIPT, "", dimsizes(dims) + \
                  " dimensions, need 2 or 3")
      end if
      idx = ind(dims .eq. "lat")
      if (ismissing(idx)) then
        error_msg("f", DIAG_SCRIPT, "", "no lat dimension")
      end if
      idx = ind(dims .eq. "lon")
      if (ismissing(idx)) then
        error_msg("f", DIAG_SCRIPT, "", "no lon dimension")
      end if

      ; average over time

      ; if variable is an error variable, we have to square it before
      ; averaging and then calculate the square-root afterwards

      if (treat_var_as_error) then
        log_info(" ++++++++++++++ Treating variable as error " + \
                 "variable when averaging ")
        A0 = A0 * A0
      end if

      data1 = time_operations(A0, -1, -1, "average", timemean, True)

      if (treat_var_as_error) then
        data1 = sqrt(data1)
      end if

      delete(A0)
    else
      data1 = ref_avg
    end if

    ; if requested, select geographical region

    if (isatt(diag_script_info, "region")) then
      region = diag_script_info@region
      data1 := area_operations(data1, region(0), region(1), region(2), \
                               region(3), "extract", False)
      if (region(2).eq.0. .and. region(3).eq.360.) then
      else
        res@gsnAddCyclic = False
      end if
      res@mpMinLatF = region(0)    ; range to zoom in on
      res@mpMaxLatF = region(1)
      res@mpMinLonF = region(2)
      res@mpMaxLonF = region(3)
      res@mpCenterLonF = 0.5 * (region(2) + region(3))
      delete(region)
    end if

    ; ###########################################
    ; # Style dependent annotation              #
    ; ###########################################
    ; retrieve unique strings describing the data
    ; function in ./diag_scripts/shared/plot/style.ncl

    ; ###########################################
    ; # plot ressources                         #
    ; ###########################################

    res@cnFillOn       = True      ; color plot desired
    res@cnLineLabelsOn = False     ; contour lines

    ; colors
    ; http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml

    ; annotation

;    res@tiMainOn             = False
    res@tiMainString = names(imod)

    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLinesOn            = False

    res@mpOutlineOn          = True
    res@mpFillOn             = False

    ; variable specific plotting settings

    if (var0.eq."pr") then
      res@cnLevels = fspan(0.5, 10, 20)
      ; convert from kg m-2 s-1 to mm day-1
      data1 = data1 * 86400.0
      data1@units = "mm day-1"
    end if

    if (var0.eq."lwp") then
      res@cnLevels            = ispan(10, 200, 10) * 0.001
;      res@mpOutlineOn         = False
;      res@mpFillOn            = True
;      res@mpLandFillColor     = "Black"
      pal = read_colormap_file("$diag_scripts/shared/plot/rgb/qcm3.rgb")
      res@cnFillColors        = pal
    end if

    if (var0.eq."tas") then
      res@cnLevels            = ispan(-30, 30, 3)
      pal = read_colormap_file("$diag_scripts/shared/plot/rgb/ipcc-tas.rgb")
      res@cnFillColors        = pal
      ; convert from K to degC
      data1 = data1 - 273.15
      data1@units = "degC"
    end if

    if ((var0.eq."clt") .or. (var0.eq."cltisccp")) then
      res@cnLevels            = fspan(5, 100, 20)
    end if

    if (var0.eq."clivi") then
      res@cnLevels            = ispan(10, 200, 10) * 0.001
    end if

    if (var0.eq."clwvi") then
      res@cnLevels            = ispan(10, 300, 10) * 0.001
    end if

    if (var0.eq."swcre") then
      res@cnLevels            = ispan(-100, 0, 10)
    end if

    if (var0.eq."lwcre") then
      res@cnLevels            = ispan(0, 100, 10)
    end if

    if (var0.eq."netcre") then
      res@cnLevels            = ispan(-70, 70, 10)
    end if

;    res@lbLabelBarOn         = False
    res@gsnRightString       = ""

    res@mpFillDrawOrder       = "PostDraw"    ; draw map last
    res@cnMissingValFillColor = "Gray"

    ; no tickmarks and no labels

    res@tmYLLabelsOn       = False
    res@tmYLOn             = False
    res@tmYRLabelsOn       = False
    res@tmYROn             = False
    res@tmXBLabelsOn       = False
    res@tmXBOn             = False
    res@tmXTLabelsOn       = False
    res@tmXTOn             = False
    res@cnInfoLabelOn      = False    ; turn off cn info label
    res@mpPerimOn          = perim    ; draw line around map

    ; specified in namelist

    res@mpProjection       = projection

    ; set explicit contour levels

    if (isatt(diag_script_info, "explicit_cn_levels")) then
      res@cnLevelSelectionMode = "ExplicitLevels"
      res@cnLevels = diag_script_info@explicit_cn_levels
    end if

    if (.not. isatt(res, "cnLevels")) then
      log_info(DIAG_SCRIPT + " (var: " + var0 + "):")
      log_info("info: using default contour levels")
      res@cnLevels = fspan(min(data1), max(data1), 20)
    end if

    ; ###########################################
    ; # other Metadata: diag_script, var        #
    ; ###########################################
    ; add to data1 as attributes without prefix

    if (isatt(data1, "diag_script")) then  ; add to existing entries
      temp = data1@diag_script
      delete(data1@diag_script)
      data1@diag_script = array_append_record(temp, (/DIAG_SCRIPT/), 0)
      delete(temp)
    else  ; add as new attribute
      data1@diag_script = (/DIAG_SCRIPT/)
    end if

    if (isatt(variable_info[varidx], "long_name")) then
      data1@var_long_name = variable_info[varidx]@long_name
    end if

    data1@var = var0

    if (isatt(variable_info[varidx], "units")) then
      data1@var_units = variable_info[varidx]@units
    else
      data1@var_units = ""
    end if

    if (.not. isvar("ref_data")) then
      ref_data = data1
    end if

    ; check if data are on same grid (for calculating difference, RMSD,
    ; correlation)

    same_grid = False

    if (all(dimsizes(ref_data) .eq. dimsizes(data1))) then
      if (max(abs(ref_data&lat - data1&lat)) .le. 1.0e-6) then
        if (max(abs(ref_data&lon - data1&lon)) .le. 1.0e-6) then
          same_grid = True
        end if
      end if
    end if

    if (flag_diff .and. .not.same_grid) then
      flag_diff = False
      error_msg("f", DIAG_SCRIPT, "", \
                "Data are not on same grid, cannot calculate differences. " \
                + "Set showdiff to False in namelist or regrid data to " \
                + "common grid (check/adjust " \
                + "preprocessor settings in recipe).")
    end if

    corr = corr@_FillValue
    gavg = gavg@_FillValue

    if (.not.all(ismissing(data1))) then
      if (numseas.gt.1) then
        do is = 0, numseas - 1
          if (same_grid .and. (ref_ind .ge. 0)) then
            corr(is) = calculate_metric(ref_data(is, :, :), data1(is, :, :), \
                                        "correlation")
          end if
          gavg(is) = area_operations(data1(is, :, :), -90., 90., 0., 360., \
                                     "average", True)
        end do
      else
        if (same_grid .and. (ref_ind .ge. 0)) then
          corr(0) = calculate_metric(ref_data, data1, "correlation")
        end if
        gavg(0) = area_operations(data1, -90., 90., 0., 360., "average", True)
      end if
    end if

;    res@gsnLeftStringFontHeightF  = min((/0.025, 0.015 * 6.0 \
;                                              / tofloat((dim_MOD + 1) / 2)/))
;    res@gsnRightStringFontHeightF = min((/0.025, 0.015 * 6.0 \
;                                              / tofloat((dim_MOD + 1) / 2)/))

    ; ###########################################
    ; # create the plot                         #
    ; ###########################################

;    res@gsnDraw        = False  ; do not draw yet
;    res@gsnFrame       = False  ; don't advance frame

    ; function in aux_plotting.ncl

    if (ii.eq.0) then
      ; note: an array of workspaces (i.e. wks(numseas)) does not work as
      ;       attributes cannot be assigned to each array element
      ;       individually
      wks0 = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 + \
                     "_" + season(0) + filename_add)
      ; difference plots will be saved to a different file
      if (flag_diff) then
        wks0d = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 + \
                        "_bias_" + season(0) + filename_add)
      end if
      if (numseas.gt.1) then
        wks1 = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 + \
                       "_" + season(1) + filename_add)
        wks2 = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 + \
                       "_" + season(2) + filename_add)
        wks3 = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 + \
                       "_" + season(3) + filename_add)
        ; difference plots will be saved to a different files
        if (flag_diff) then
          wks1d = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 + \
                          "_bias_" + season(1) + filename_add)
          wks2d = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 + \
                          "_bias_" + season(2) + filename_add)
          wks3d = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 + \
                          "_bias_" + season(3) + filename_add)
        end if
      end if
    end if

    if (numseas.gt.1) then
      do is = 0, numseas - 1
        if (.not.ismissing(corr(is))) then
          res@gsnRightString = "corr = " + sprintf("%6.3f", corr(is))
        else
          res@gsnRightString = ""
        end if
        if (.not.ismissing(gavg(is))) then
          res@gsnLeftString = "mean = " + sprintf("%6.3f", gavg(is))
        else
          res@gsnLeftString = ""
        end if

        if (imod.eq.ref_ind) then  ; remove corr. string for reference dataset
          res@gsnRightString = ""
        end if

        if (is.eq.0) then
          maps(imod, is) = gsn_csm_contour_map(wks0, data1(is, :, :), res)
        end if
        if (is.eq.1) then
          maps(imod, is) = gsn_csm_contour_map(wks1, data1(is, :, :), res)
        end if
        if (is.eq.2) then
          maps(imod, is) = gsn_csm_contour_map(wks2, data1(is, :, :), res)
        end if
        if (is.eq.3) then
          maps(imod, is) = gsn_csm_contour_map(wks3, data1(is, :, :), res)
        end if
      end do
    else
      if (.not.ismissing(corr(0))) then
        res@gsnRightString = "corr = " + sprintf("%6.3f", corr(0))
      else
        res@gsnRightString = ""
      end if
      if (.not.ismissing(gavg(0))) then
        res@gsnLeftString = "mean = " + sprintf("%6.3f", gavg(0))
      else
        res@gsnLeftString = ""
      end if

      if (imod.eq.ref_ind) then  ; remove corr. string for reference dataset
        res@gsnRightString = ""
      end if

      maps(imod, 0) = gsn_csm_contour_map(wks0, data1, res)
    end if

    ; mandatory netcdf output

    data1@var = var0 + "_mean_" + names(imod)
    nc_outfile_mean = ncdf_write(data1, nc_filename_mean)

    ; =======================================================================
    ; Create difference plots (if requested)
    ; =======================================================================

    if (flag_diff) then

      dres = True

      if (imod .ne. ref_ind) then
        diff = data1
        if (flag_rel_diff) then
          diff = (diff - ref_data) / ref_data * 100.0
          diff = where(ref_data .le. rel_diff_min, diff@_FillValue, diff)
        else
          diff = diff - ref_data
        end if

        dres@gsnLeftString  = ""
        dres@gsnRightString = ""

        dres@tiMainString =  names(imod) + " - " + refname

        rmsd = rmsd@_FillValue
        bias = bias@_FillValue

        if (numseas.gt.1) then
          do is = 0, numseas - 1
            if (.not. flag_rel_diff) then
              if (same_grid) then
                rmsd(is) = calculate_metric(ref_data(is, :, :), \
                                            data1(is, :, :), "RMSD")
              end if
              bias(is) = area_operations(diff(is, :, :), -90., 90., 0., 360., \
                                         "average", True)
            end if
          end do
        else
          if (.not. flag_rel_diff) then
            if (same_grid) then
              rmsd(0) = calculate_metric(ref_data, data1, "RMSD")
            end if
            bias(0) = area_operations(diff, -90., 90., 0., 360., "average", \
                                      True)
          end if
        end if
      else if (flag_multiobs_unc) then
        diff = ref_std
        if (.not.isatt(diff, "diag_script")) then
          diff@diag_script = (/DIAG_SCRIPT/)
        end if
        dres@gsnLeftString  = ""
        dres@gsnRightString = ""
        dres@tiMainString =  refname + " uncertainty"
        rmsd = rmsd@_FillValue
        bias = bias@_FillValue
      else
        continue
      end if
      end if

      ; ----------------------------------------------------------------------

      ; ###########################################
      ; # plot ressources                         #
      ; ###########################################

;      dres@gsnLeftStringFontHeightF  = min((/0.025, 0.015 * 6.0 \
;                                           / tofloat((dim_MOD + 1) / 2)/))
;      dres@gsnRightStringFontHeightF = min((/0.025, 0.015 * 6.0 \
;                                           / tofloat((dim_MOD + 1) / 2)/))

;      dres@tiMainOn       = False

      dres@cnFillOn       = True      ; color plot desired
      dres@cnLineLabelsOn = False     ; contour lines
      dres@cnLinesOn      = False

      ; colors
      ; http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml

      ; annotation

      dres@cnLevelSelectionMode = "ExplicitLevels"
      dres@mpOutlineOn          = True
      dres@mpFillOn             = False

      ; variable specific plotting settings

      ; set contour levels / colors

      if (.not.isvar("cnLevels")) then

        if (isatt(dres, "cnLevels")) then
          delete(dres@cnLevels)
        end if
        if (isatt(dres, "cnFillColors")) then
          delete(dres@cnFillColors)
        end if
        if (isvar("pal")) then
          delete(pal)
        end if

        if (var0.eq."pr") then
          dres@cnLevels = ispan(-30, 30, 5) * 0.1
          pal = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                   + "ipcc-precip-delta.rgb")
          dres@cnFillColors = pal
          dres@lbOrientation = "horizontal"
        end if

        if ((var0.eq."tas") .or. (var0.eq."ts")) then
          pal = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                   + "ipcc-tas-delta.rgb")
          dres@cnFillPalette = pal
          if (var0.eq."ts") then
            dres@cnLevels = ispan(-5, 5, 1) * 0.5
          end if
        end if

        if (var0.eq."lwp") then
          dres@cnLevels           = ispan(-50, 50, 10) * 0.001
          pal = read_colormap_file("$diag_scripts/shared/plot/rgb/qcm3.rgb")
          dres@cnFillColors       = pal
        end if

        if (var0.eq."clt") then
          dres@cnLevels           = fspan(-25, 25, 11)
        end if

        if (var0.eq."clivi") then
          dres@cnLevels           = ispan(-70, 70, 10) * 0.001
        end if

        if (var0.eq."clwvi") then
          dres@cnLevels           = ispan(-50, 50, 10) * 0.001
        end if

        if (var0.eq."swcre") then
          dres@cnLevels          = ispan(-30, 30, 5)
        end if

        if (var0.eq."lwcre") then
          dres@cnLevels          = ispan(-30, 30, 5)
        end if

        if (var0.eq."netcre") then
          dres@cnLevels          = ispan(-30, 30, 5)
        end if

        ; ******************************************************
        ; *** relative differences: use specific color table ***
        ; ******************************************************

        if (flag_rel_diff) then
          if (isatt(dres, "cnLevels")) then
            delete(dres@cnLevels)
          end if
          if (isatt(dres, "cnFillColors")) then
            delete(dres@cnFillColors)
          end if
          dres@cnLevels = fspan(-100, 100, 21)
          if (isvar("pal")) then
            delete(pal)
          end if
          pal = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                   + "percent100.rgb")
          dres@cnFillColors = pal
        end if

        ; ******************************************************

        if (.not. isatt(dres, "cnLevels")) then
          log_info(DIAG_SCRIPT + " (var: " + var0 + "):")
          log_info("info: using default contour levels")
          dres@cnLevels = fspan(min(diff), max(diff), 20)
        end if

        cnLevels = dres@cnLevels
        if (isatt(dres, "cnFillColors")) then
          cnFillColors = dres@cnFillColors
        end if

      else  ; use previously defined colors and contour intervals

        if (isatt(dres, "cnLevels")) then
          delete(dres@cnLevels)
        end if
        if (isatt(dres, "cnFillColors")) then
          delete(dres@cnFillColors)
        end if

        dres@cnLevels = cnLevels

        if (isvar("cnFillColors")) then
          dres@cnFillColors = cnFillColors
        end if

      end if  ; if .not.isvar("cnLevels")

;      if (imod.eq.ref_ind) then
;        dres@lbLabelBarOn = True
;      else
;        dres@lbLabelBarOn = False
;      end if

      ; map attributes

      dres@mpFillDrawOrder       = "PostDraw"    ; draw map last
      dres@cnMissingValFillColor = "Gray"

      ; no tickmarks and no labels

      dres@tmYLLabelsOn       = False
      dres@tmYLOn             = False
      dres@tmYRLabelsOn       = False
      dres@tmYROn             = False
      dres@tmXBLabelsOn       = False
      dres@tmXBOn             = False
      dres@tmXTLabelsOn       = False
      dres@tmXTOn             = False
      dres@cnInfoLabelOn      = False    ; turn off cn info label

      ; specified in namelist

      dres@mpProjection       = projection

      ; set explicit contour levels

      if (isatt(diag_script_info, "explicit_cn_levels")) then
        dres@cnLevelSelectionMode = "ExplicitLevels"
        if (isatt(dres, "cnLevels")) then
          delete(dres@cnLevels)
        end if
        dres@cnLevels = diag_script_info@explicit_cn_levels
      end if

      ; ###########################################
      ; # other Metadata: diag_script, var        #
      ; ###########################################
      ; add to diff as attributes without prefix

      if (isatt(variable_info, "long_name")) then
        diff@var_long_name = variable_info@long_name
      end if
      if (isatt(variable_info, "units")) then
        diff@var_units = variable_info@units
      else
        diff@var_units = ""
      end if

      ; ###########################################
      ; # create the plot                         #
      ; ###########################################

      if (flag_multiobs_unc) then
        dres@gsnDraw        = False  ; do not draw yet
        dres@gsnFrame       = False  ; don't advance frame
      end if

      ; ----------------------------------------------------------------------

      if (numseas.gt.1) then
        do is = 0, numseas - 1
          if (.not.ismissing(rmsd(is))) then
            dres@gsnRightString = "rmsd = " + sprintf("%6.3f", rmsd(is))
          else
            dres@gsnRightString = ""
          end if
          if (.not.ismissing(bias(is))) then
            dres@gsnLeftString = "bias = " + sprintf("%6.3f", bias(is))
          else
            dres@gsnLeftString = ""
          end if

          if (is.eq.0) then
            maps_d(imod, is) = gsn_csm_contour_map(wks0d, diff(is, :, :), dres)
          end if
          if (is.eq.1) then
            maps_d(imod, is) = gsn_csm_contour_map(wks1d, diff(is, :, :), dres)
          end if
          if (is.eq.2) then
            maps_d(imod, is) = gsn_csm_contour_map(wks2d, diff(is, :, :), dres)
          end if
          if (is.eq.3) then
            maps_d(imod, is) = gsn_csm_contour_map(wks3d, diff(is, :, :), dres)
          end if
        end do
      else
        if (.not.ismissing(rmsd(0))) then
          dres@gsnRightString = "rmsd = " + sprintf("%6.3f", rmsd(0))
        else
          dres@gsnRightString = ""
        end if
        if (.not.ismissing(bias(0))) then
          dres@gsnLeftString = "bias = " + sprintf("%6.3f", bias(0))
        else
          dres@gsnLeftString = ""
        end if
        maps_d(imod, 0) = gsn_csm_contour_map(wks0d, diff, dres)

        ; obs uncertainty

        if (flag_multiobs_unc) then
          xmask = where(abs(diff) .gt. ref_std, 0.0, 1.0)
          copy_VarMeta(diff, xmask)

;          if (ii .eq. 0) then
;            system("rm debug.nc")
;            debugfile = addfile("debug.nc", "c")
;          end if
;          xm = "xmask" + names(imod)
;          debugfile->$xm$ = xmask
;          debugfile->obs_std = ref_std
;          dn = "diff" + names(imod)
;          debugfile->$dn$ = diff

          xres = True
          xres@gsnDraw        = False  ; do not draw yet
          xres@gsnFrame       = False  ; don't advance frame
          xres@cnMissingValFillColor = -1
          xres@cnLevelSelectionMode = "ExplicitLevels"
          xres@cnLevels = 0.5
          xres@cnFillColors = (/"black", "transparent"/)
          xres@cnFillPattern = 17
          xres@cnFillOn       = True      ; color plot desired
          xres@cnInfoLabelOn = False
          xres@cnLinesOn = False
          xres@cnLineLabelsOn = False
          xres@lbLabelBarOn = False
          xres@gsnRightString   = ""
          xres@gsnLeftString    = ""
          xres@gsnCenterString  = ""
          xres@tiYAxisOn = False
          xres@tmXBBorderOn = False
          xres@tmXTBorderOn = False
          xres@tmYLBorderOn = False
          xres@tmYRBorderOn = False
          xres@tmXBLabelsOn = False
          xres@tmYLLabelsOn = False
          xres@tmXBOn = False
          xres@tmXTOn = False
          xres@tmYLOn = False
          xres@tmYROn = False
          xres@cnConstFEnableFill = True

          if (imod .ne. ref_ind) then
            plotmask = gsn_csm_contour(wks0d, xmask, xres)
            overlay(maps_d(imod, 0), plotmask)
            delete(plotmask)
          end if
          delete(xmask)
        end if
      end if  ; if flag_multiobs_unc

      draw(maps_d(imod, 0))
      frame(wks0d)

      ; mandatory netcdf output

      diff@var = var0 + "_bias_" + names(imod)
      nc_outfile_bias = ncdf_write(diff, nc_filename_bias)

    end if  ; if flag_diff

    ; =======================================================================

  end do  ; ii-loop (models)

  if (write_plots) then
    ; save default color map in case it is needed later for optionally
    ; plotting color bar to a separate file

    tmp_colors = gsn_retrieve_colormap(wks0)
    cdims = dimsizes(tmp_colors)
    nboxes = dimsizes(res@cnLevels)
    clen = cdims(0)
    stride = max((/1, ((clen(0)-1) - 2) / nboxes /))
    fill_colors = ispan(2, clen(0) - 1, stride)
    mean_colors = tmp_colors(fill_colors, :)
    delete(tmp_colors)
    delete(fill_colors)
    delete(cdims)

    ; sort plots if needed (observations go first)

    plottmp = ispan(0, dim_MOD - 1, 1)
    plotind = plottmp

    ; move plots of observational datasets (if present) into the first line(s)
    ; of the panel plot

    j = 0
    do i = 0, dimsizes(plottmp) - 1
      if (i.eq.ref_ind) then
        plotind(j) = plottmp(i)
        j = j + 1
      else if (plottmp(i) .lt. dimsizes(projects)) then
        if (isStrSubset(str_lower(projects(plottmp(i))), \
                        "obs")) then
          plotind(j) = plottmp(i)
          j = j + 1
        end if
      end if
      end if
    end do

    do i = 0, dimsizes(plottmp) - 1
      if ((isStrSubset(str_lower(projects(plottmp(i))), \
                       "obs")).or.(i.eq.ref_ind)) then
      else
        plotind(j) = plottmp(i)
        j = j + 1
      end if
    end do

    pres                      = True    ; needed to override
                                        ; panelling defaults
;    pres@gsnPanelLabelBar     = True    ; add common colorbar
    if (panel_labels) then
      ; print dataset name on each panel
      pres@gsnPanelFigureStrings = names(plotind)
    end if
    pres@gsnPanelFigureStringsFontHeightF = min((/0.01, 0.01 * 6.0 \
                                                / tofloat((dim_MOD + 1) / 2)/))
    pres@lbLabelFontHeightF               = min((/0.015, 0.01 * 6.0 \
                                                / tofloat((dim_MOD + 1) / 2)/))
    pres@lbAutoManage                     = False
    pres@lbTopMarginF                     = 0.1
    pres@lbTitleOn                        = True
    pres@lbTitleFontHeightF               = min((/0.015, 0.01 * 6.0 \
                                                / tofloat((dim_MOD + 1) / 2)/))
    pres@lbTitlePosition                  = "Bottom"
    pres@lbTitleString                    = data1@long_name + " (" \
      + data1@units + ")"
    pres@lbPerimOn                        = False   ; draw line around label
                                                    ; bar area
    pres@gsnPanelCenter                   = False
    if (dim_MOD.le.8) then
      pres@pmLabelBarOrthogonalPosF       = -0.03
    else
      pres@pmLabelBarOrthogonalPosF       = -0.01   ; shift label bar a bit to
                                                    ; the bottom
    end if

    if (embracesetup) then
      if (numseas.gt.1) then
        pres@txString = season(0)
        outfile(0) = panelling(wks0, maps(plotind, 0), (dim_MOD + 3) / 4, \
                               4, pres)

        pres@txString = season(1)
        outfile(1) = panelling(wks1, maps(plotind, 1), (dim_MOD + 3) / 4, \
                               4, pres)

        pres@txString = season(2)
        outfile(2) = panelling(wks2, maps(plotind, 2), (dim_MOD + 3) / 4, \
                               4, pres)

        pres@txString = season(3)
        outfile(3) = panelling(wks3, maps(plotind, 3), (dim_MOD + 3) / 4, \
                               4, pres)
        log_info(" Wrote " + outfile)
      else
        pres@gsnPanelRowSpec = True             ; tell panel what order to plt
        pres@gsnPanelYWhiteSpacePercent = 5
        pres@gsnPanelXWhiteSpacePercent = 5
        if (isatt(diag_script_info, "PanelTop")) then
          top = tofloat(diag_script_info@PanelTop)
        else
          top = 0.99  ; default
        end if
        pres@gsnPanelTop = top

        if (isvar("plotsperline")) then
          delete(plotsperline)
        end if

        plotsperline = new((dim_MOD + 1) / 2, integer)
        plotsperline = 2

        if ((isStrSubset(str_lower(projects(plotind(0))), \
                         "obs")).and. \
           .not.(isStrSubset(str_lower(projects(plotind(1))), \
                             "obs"))) then
          plotsperline(0) = 1
        end if

        if (sum(plotsperline).gt.dimsizes(plotind)) then
          plotsperline(dimsizes(plotsperline) - 1) = 1
        end if

        if (sum(plotsperline).lt.dimsizes(plotind)) then
          xadd = 1
          xtmp = array_append_record(plotsperline, xadd, 0)
          delete(plotsperline)
          plotsperline = xtmp
          delete(xtmp)
        end if

        gsn_panel(wks0, maps(plotind, 0), plotsperline, pres)
        outfile(0) = wks0@fullname
      end if
    else  ; if embracesetup
      if (numseas.gt.1) then
        pres@txString = season(0)
        outfile(0) = panelling(wks0, maps(plotind, 0), (dim_MOD + 3) / 4, \
                               4, pres)

        pres@txString = season(1)
        outfile(1) = panelling(wks1, maps(plotind, 1), (dim_MOD + 3) / 4, \
                               4, pres)

        pres@txString = season(2)
        outfile(2) = panelling(wks2, maps(plotind, 2), (dim_MOD + 3) / 4, \
                               4, pres)

        pres@txString = season(3)
        outfile(3) = panelling(wks3, maps(plotind, 3), (dim_MOD + 3) / 4, \
                               4, pres)
      else
        outfile(0) = panelling(wks0, maps(plotind, 0), (dim_MOD + 3) / 4, \
                               4, pres)
      end if
    end if  ; if embracesetup

    do is = 0, numseas - 1
      log_info("Wrote " + outfile(is))
    end do

    ; ------------------------------------------------------------------------
    ; write provenance to netcdf output and plot file(s) (mean)
    ; ------------------------------------------------------------------------

    statistics = (/"clim", "mean"/)
    if (isatt(diag_script_info, "region")) then
      domain = "reg"
    else
      domain = "global"
    end if
    plottype = "geo"

    do is = 0, numseas - 1
      caption = "Mean values for variable " + var0 \
                + " (" + allseas + ")."
      log_provenance(nc_outfile_mean, outfile(is), caption, statistics, \
                     domain, plottype, "", "", climofiles)
    end do

    ; ========================================================================

    if (flag_diff) then
      pres@lbTitleString = "~F33~D~F21~" + diff@long_name + " (" + \
                           diff@units + ")"

      ; save default color map in case it is needed later for optionally
      ; plotting color bar to a separate file

      if (isvar("nboxes")) then
        delete(nboxes)
      end if

      tmp_colors = gsn_retrieve_colormap(wks0d)
      cdims = dimsizes(tmp_colors)
      nboxes = dimsizes(dres@cnLevels)
      clen = cdims(0)
      stride = max((/1, ((clen(0)-1) - 2) / nboxes /))
      fill_colors = ispan(2, clen(0) - 1, stride)
      diff_colors = tmp_colors(fill_colors, :)
      delete(tmp_colors)
      delete(fill_colors)
      delete(cdims)

      if (isvar("plottmp")) then
        delete(plottmp)
      end if

      if (isvar("plotind")) then
        delete(plotind)
      end if

      plottmp = ind(ispan(0, dim_MOD - 1, 1).ne.ref_ind)
      plotind = plottmp

      ; if there is a second observational dataset, move the corresponding
      ; plot to the first line of the panel plot

      j = 0
      do i = 0, dimsizes(plottmp) - 1
        if (isStrSubset(str_lower(projects(plottmp(i))), "obs")) then
          plotind(j) = plottmp(i)
          j = j + 1
        end if
      end do
      do i = 0, dimsizes(plottmp) - 1
        if (isStrSubset(str_lower(projects(plottmp(i))), "obs")) then
        else
          plotind(j) = plottmp(i)
          j = j + 1
        end if
      end do

      if (isatt(pres, "gsnPanelFigureStrings")) then
        delete(pres@gsnPanelFigureStrings)
      end if
      if (panel_labels) then
        pres@gsnPanelFigureStrings = names(plotind)
      end if

      if (dimsizes(plotind).eq.1) then
        pres@gsnPanelRight = 0.5
      end if

      if (embracesetup) then
        if (numseas.gt.1) then
          pres@txString = season(0)
          outfile_d(0) = panelling(wks0d, maps_d(plotind, 0), \
                                   (dim_MOD + 3) / 4, 4, pres)

          pres@txString = season(1)
          outfile_d(1) = panelling(wks1d, maps_d(plotind, 1), \
                                   (dim_MOD + 3) / 4, 4, pres)

          pres@txString = season(2)
          outfile_d(2) = panelling(wks2d, maps_d(plotind, 2), \
                                   (dim_MOD + 3) / 4, 4, pres)

          pres@txString = season(3)
          outfile_d(3) = panelling(wks3d, maps_d(plotind, 3), \
                                   (dim_MOD + 3) / 4, 4, pres)
        else
          pres@gsnPanelRowSpec = True           ; tell panel what order to plt
          pres@gsnPanelYWhiteSpacePercent = 5
          pres@gsnPanelXWhiteSpacePercent = 5
          pres@gsnPanelTop = tofloat(diag_script_info@PanelTop)

          if (isvar("plotsperline")) then
            delete(plotsperline)
          end if

          plotsperline = new(max((/1, dim_MOD / 2/)), integer)
          plotsperline = 2

          if (dimsizes(plotind).gt.1) then
            if ((isStrSubset(str_lower(projects(plotind(0))), "obs")).and. \
               .not. \
                (isStrSubset(str_lower(projects(plotind(1))), "obs"))) then
              plotsperline(0) = 1
            end if
          end if

          if (sum(plotsperline).gt.dimsizes(plotind)) then
            plotsperline(dimsizes(plotsperline) - 1) = 1
          end if

          if (sum(plotsperline).lt.dimsizes(plotind)) then
            xadd = 1
            xtmp = array_append_record(plotsperline, xadd, 0)
            delete(plotsperline)
            plotsperline = xtmp
            delete(xtmp)
          end if

          gsn_panel(wks0d, maps_d(plotind, 0), plotsperline, pres)
          outfile_d(0) = wks0d@fullname
        end if
      else  ; embracesetup = False
        if (numseas.gt.1) then
          pres@txString  = season(0)
          outfile_d(0) = panelling(wks0d, maps_d(plotind, 0), \
                                   (dim_MOD + 3) / 4, 4, pres)

          pres@txString  = season(1)
          outfile_d(1) = panelling(wks1d, maps_d(plotind, 1), \
                                   (dim_MOD + 3) / 4, 4, pres)

          pres@txString  = season(2)
          outfile_d(2) = panelling(wks2d, maps_d(plotind, 2), \
                                   (dim_MOD + 3) / 4, 4, pres)

          pres@txString  = season(3)
          outfile_d(3) = panelling(wks3d, maps_d(plotind, 3), \
                                   (dim_MOD + 3) / 4, 4, pres)
        else
          outfile_d(0) = panelling(wks0d, maps_d(plotind, 0), \
                                   (dim_MOD + 3) / 4, 4, pres)
        end if
      end if  ; end if embracesetup

      do is = 0, numseas - 1
        log_info(" Wrote " + outfile(is))

        ; --------------------------------------------------------------------
        ; write provenance to netcdf output and plot file(s) (bias)
        ; --------------------------------------------------------------------

        statistics = (/"clim", "diff"/)
        if (isatt(diag_script_info, "region")) then
          domain = "reg"
        else
          domain = "global"
        end if
        plottype = "geo"

        ; note: because function log_provenance does not yet support to attach
        ;       different captions to netcdf (contains all seasons) and plots
        ;       (contain one season each), the caption cannot specifiy the
        ;       season plotted; using "annual" or "DJF/MAM/JJA/SON" instead.

        caption = "Differences for variable " + var0 \
                  + " (" + allseas + "), reference = " + refname + "."
        log_provenance(nc_outfile_bias, outfile_d(is), caption, statistics, \
                       domain, plottype, "", "", climofiles)
      end do

    end if  ; if flag_diff

    ; optionally save legend(s) to extra file(s)

    if (extralegend) then
      nboxes = dimsizes(res@cnLevels) + 1
      wksleg = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 \
                       + "_legend")
      pres@lbMonoFillPattern = True
      pres@lbOrientation = "Horizontal"
      pres@vpWidthF = 0.7
      pres@vpHeightF = 0.1
      pres@lbLabelFontHeightF = 0.015
      pres@lbLabelAlignment = "InteriorEdges"
      pres@lbTitleFontHeightF = 0.015
      pres@lbTitleString = data1@long_name + " (" + data1@units + ")"

      labels = tostring(res@cnLevels)

      ; remove trailing zeros from strings

      do i = 0, dimsizes(labels) - 1
        i1  = str_index_of_substr(labels(i), ".", -1)
        if (.not.ismissing(i1)) then
          tmp = stringtochar(labels(i))
          do j = dimsizes(tmp) - 2, i1, 1
            if ((tmp(j).ne.".").and.(tmp(j).ne."0")) then
              break
            end if
          end do
          labels(i) = chartostring(tmp(0:j))
          delete(tmp)
        end if
      end do

      if (isatt(res, "cnFillColors")) then
        pres@lbFillColors = res@cnFillColors
      else if (isatt(res, "cnFillPalette")) then
        pres@lbFillColors = res@cnFillPalette
      else
        pres@lbFillColors = mean_colors  ; default colors
      end if
      end if

      gsn_labelbar_ndc(wksleg, nboxes, labels, 0.1, 0.9, pres)

      delete(wksleg)
      delete(labels)
      delete(pres@lbFillColors)

      if (flag_diff) then
        nboxes = dimsizes(dres@cnLevels) + 1
        wksleg = get_wks("dummy_for_wks", DIAG_SCRIPT, "clouds_" + var0 \
                         + "_diff_legend")

        labels = tostring(dres@cnLevels)

        ; remove trailing zeros from strings

        do i = 0, dimsizes(labels) - 1
          i1  = str_index_of_substr(labels(i), ".", -1)
          if (.not.ismissing(i1)) then
            tmp = stringtochar(labels(i))
            do j = dimsizes(tmp) - 2, i1, 1
              if ((tmp(j).ne.".").and.(tmp(j).ne."0")) then
                break
              end if
            end do
            labels(i) = chartostring(tmp(0:j))
            delete(tmp)
          end if
        end do

        if (flag_rel_diff) then
          pres@lbTitleString = "~F33~D~F21~" + data1@long_name + " (%)"
        else
          pres@lbTitleString = "~F33~D~F21~" + data1@long_name + " (" + \
                        data1@units + ")"
        end if

        if (isatt(dres, "cnFillColors")) then
          pres@lbFillColors = dres@cnFillColors
        else if (isatt(dres, "cnFillPalette")) then
          pres@lbFillColors = dres@cnFillPalette
        else
          pres@lbFillColors = diff_colors  ; default colors
        end if
        end if

        gsn_labelbar_ndc(wksleg, nboxes, labels, 0.1, 0.9, pres)
      end if  ; if (flag_diff)
    end if  ; if (extralegend)
  end if  ; if write_plots

  ; ==========================================================================

  leave_msg(DIAG_SCRIPT, "")

end
