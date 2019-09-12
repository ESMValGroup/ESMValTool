; #############################################################################
; eyring13jgr_fig01.ncl
; #############################################################################
; Author: Irene Cionni (ENEA, Italy)
; ESMVal project PI: Veronika Eyring (DLR, Germany)
; #############################################################################
;
; Description:
;   Calculates seasonal cycles of zonally averaged total ozone columns.
;
; Required diag_script_info attributes:
;   none
;
; Optional diag_script_info attributes:
;   e13fig01_diff_ref: subtract reference dataset (True, False);
;     true requires definition of "reference_dataset" (see below)
;   e13fig01_list_chem_mod: list of models in the group “chem” (array of
;     strings, default = (/"All"/))
;   e13fig01_list_chem_mod_string: plotting label for group “chem”,
;     e.g., "CMIP5"
;   e13fig01_list_nochem_mod: list of models in the group “nochem”
;     (array of strings, default = (/""/))
;   e13fig01_list_nochem_mod_string: plotting label for group “nochem”,
;     e.g., "NoChem"
;   eyring13_cnlvl: contour levels for plots (e.g. from Eyring et al. (2013))
;   eyring13_difflvl: contour levels for difference plots
;   eyring13_reflvl: contour levels for plots (reference dataset)
;   font: overrides default font (e.g., 21, see
;     www.ncl.ucar.edu/Document/Graphics/Resources/tx.shtml#txFont)
;   lbLabelBarOn: plot a label bar (True, False)
;   range_option: 0 = as in nml, 1 = overlapping time period
;   rgb_file: path + filename of color table
;     (e.g., "diag_scripts/shared/plot/rgb/eyring_toz.rgb")
;   show_mask: plot mask (True, False)
;
; Required variable_info attributes:
;
; Optional variable_info attributes:
;   reference_dataset: name of reference dataset (observations)
;
; Caveats:
;
; Modification history:
;   20190503-A_laue_ax: code rewritten for v2, regridding now done by backend
;   20170810-A_bock_ls: added tags for reporting
;   20160926-A_laue_ax: bugfix for global mean values: missing values are
;                       now taken into account for all 12 month;
;                       added option to hide/show mask in all plots
;                       (diag_script_info@show_mask)
;   20160825-A_laue_ax: replaced "area_hi2lores_Wrap" with
;                       "linint2_Wrap" to prevent problems with grids
;                       that have a lower resolution than the destination
;                       grid in only one coordinate (lat OR lon), which
;                       makes area_hi2lores fail.
;   20151208-A_laue_ax: code clean-up, moved call to write_references to
;                       beginning of code, using explicit color scale
;                       levels from Eyring et al. (2013), using color map
;                       specified by diag_script_info@rgb_file
;   20140528 A_GoKl written based on ancyc_lat.ncl
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/style.ncl"


; fetch general parameters
begin
  enter_msg(DIAG_SCRIPT, "")

  var0 = variable_info[0]@short_name
  info0 = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info0)
  if (isatt(variable_info[0], "reference_dataset")) then
    refname = variable_info[0]@reference_dataset
  end if
  names = metadata_att_as_array(info0, "dataset")
  projects = metadata_att_as_array(info0, "project")
  exps = metadata_att_as_array(info0, "exp")
  ensembles = metadata_att_as_array(info0, "ensemble")

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Set default values for optional diag_script_info attributes

  set_default_att(diag_script_info, "e13fig01_diff_ref", False)
  set_default_att(diag_script_info, "e13fig01_list_chem_mod", "All")
  set_default_att(diag_script_info, "font", 21)
  set_default_att(diag_script_info, "lbLabelBarOn", True)
  set_default_att(diag_script_info, "range_option", 0)

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
end

; main code

begin
  ; ********************diag_script_info*****************

  if (isvar("refname")) then
    OBS_ref = refname
    ind_OBS_ref = ind(names.eq.OBS_ref)
    if (ismissing(ind_OBS_ref) .and. diag_script_info@e13fig01_diff_ref)
      error_msg("f", DIAG_SCRIPT, "", "reference dataset not found")
    end if
  else
    OBS_ref = ""
    ind_OBS_ref = -1
    if (diag_script_info@e13fig01_diff_ref)
      error_msg("f", DIAG_SCRIPT, "", "no reference dataset defined")
    end if
  end if

  if (diag_script_info@e13fig01_list_chem_mod.eq."All") then
    tmp = ind(projects.ne."OBS" .and. projects.ne."obs4mips")
    ind_chem_mod = ind(tmp.ne.ind_OBS_ref)
    list_chem_mod = names(ind_chem_mod)
    delete(tmp)
    attrib_chem = " "
  else
    list_chem_mod = diag_script_info@e13fig01_list_chem_mod
    ind_chem_mod = ind(names.eq.list_chem_mod)
    attrib_chem = "_CHEM"
  end if

  if (all(ismissing(ind_chem_mod))) then
    error_msg("f", DIAG_SCRIPT, "", "no chem_mod datasets found")
  else
    num_chem_mod = dimsizes(ind_chem_mod)
  end if

  if ((isatt(diag_script_info, "e13fig01_list_nochem_mod"))) then
    list_nochem_mod = diag_script_info@e13fig01_list_nochem_mod
    ind_nochem_mod = ind(names.eq.list_nochem_mod)
    if (all(ismissing(ind_nochem_mod))) then
      log_info("Warning: none of the datasets specified in " + \
               "e13fig01_list_nochem_mod found.")
      num_nochem_mod = 0
    else
      num_nochem_mod = dimsizes(ind_nochem_mod)
    end if
  else
    num_nochem_mod = 0
  end if

  ; *************************************************

  ; default time range & annots
  y1 = metadata_att_as_array(info0, "start_year")
  y2 = metadata_att_as_array(info0, "end_year")
  annots = project_style(info0, diag_script_info, "annots")

  range_option = diag_script_info@range_option
  if ((range_option .ne. 0) .and. (range_option .ne. 1)) then
    error_msg("f", DIAG_SCRIPT, "", \
              "invalid range_option; valid settings: 0 or 1")
  end if

  ; Use overlapping time period or as given in nml?
  if (range_option .eq. 1) then
    y1 = max(toint(y1))
    y2 = min(toint(y2))
  end if

  ; Loop over models and create mask for missing values
  do idx_mod = 0, dim_MOD - 1
    ; Extract data
    A0 = read_data(info0[idx_mod])
    temp = time_operations(A0, y1(idx_mod), y2(idx_mod), "extract", "", 0)

    ; check for supported variable dimensions
    dims = getvardimnames(A0)
    if (any(dims .eq. "lon") .and. any(dims .eq. "lat") .and. \
        any(dims .eq. "time")) then
    else
      error_msg("f", DIAG_SCRIPT, "", \
                "only variables with dimensions lon, lat, time " \
                + "supported by this diagnostic script.")
    end if

    delete(A0)

    ; seasonal cycle (note: returns time dimension "month")
    A1 = time_operations(temp, -1, -1, "average", "monthlyclim", True)
    delete(temp)

    ; initialize mask...
    if (idx_mod.eq.0) then
      if (.not.isatt(A1, "_FillValue")) then
        Fill = default_fillvalue(typeof(A1))
      else
        Fill = A1@_FillValue
      end if

      var_grd = new((/dim_MOD, 12, dimsizes(A1&lat), \
                    dimsizes(A1&lon)/), typeof(A1), Fill)
      var_grd!0 = "models"
      var_grd!1 = "month"
      var_grd!2 = "lat"
      var_grd!3 = "lon"
      var_grd&models = annots
      var_grd&month = A1&month
      var_grd&lat = A1&lat
      var_grd&lon = A1&lon

      mask0 = A1
      mask0(:, :, :) = 1

      var_grd@_FillValue = Fill  ; ensure consistent _FillValue
      delete(dims)
    end if

    var_grd(idx_mod, :, :, :) = where(ismissing(A1), Fill, A1)

    ; Build the mask: missing in at least one model = missing in the mask
    mask0 = where(ismissing(A1), mask0@_FillValue, mask0)

    ; clean up
    delete(A1)
  end do  ; idx_MOD

  if ((isatt(diag_script_info, "show_mask"))) then
    if (diag_script_info@show_mask) then
      do imod = 0, dim_MOD - 1
        var_grd(imod, :, :, :) =  var_grd(imod, :, :, :) * mask0
      end do
    end if
  end if

  log_info(DIAG_SCRIPT + ": Mean value calculated after weighted " + \
           "averaging over time, regridding (lon, lat) to common grid, " + \
           "masking all points that are missing in at least one model, " + \
           "considering area weighting")

  ; General settings for this plot
  res = True
  res@cnFillOn = True
  res@cnLinesOn = False

  if (isatt(diag_script_info, "rgb_file")) then
    pal = read_colormap_file(diag_script_info@rgb_file)
    res@cnFillPalette = pal
  end if

  res@cnLevelSelectionMode = "ExplicitLevels"
  if (isatt(diag_script_info, "eyring13_cnlvl")) then
    ; scale bar from Eyring et al. (2013)
    res@cnLevels = diag_script_info@eyring13_cnlvl
  else
    res@cnLevels = fspan(min(var_grd), max(var_grd), 20)
  end if

  res@sfXArray = ispan(0, 12, 1)
  res@txFont = diag_script_info@font
  res@tmXBLabelFont = diag_script_info@font
  res@tmYLLabelFont = diag_script_info@font
  res@tmXBMode = "Explicit"
  res@tmXBValues = ispan(0, 12, 1)
  res@tmXBLabels = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", "O", \
                     "N", "D", "J"/)
  res@lbLabelBarOn = diag_script_info@lbLabelBarOn
  res@lbLabelFont = diag_script_info@font
  res@lbTitleFont = diag_script_info@font
  res@lbTitleFontHeightF = 0.01  ; default: 0.025
  res@lbTitleOffsetF = -0.25  ; default: 0.03 (negative values move it down)
  res@lbTopMarginF = 0.1  ; default: 0.05 (negative values move it up)
  res@lbTitleJust = "BottomCenter"  ; default: CenterCenter
  res@lbLabelFontHeightF = 0.025

  ; Loop over datasets for plotting
  do idx_mod = 0, dim_MOD - 1
    ; calculate mean
    modeldata = var_grd(idx_mod, :, :, :) * mask0
    copy_VarMeta(mask0, modeldata)
    area_mean = area_operations(modeldata, -90, 90, 0, 360, "average", True)
    time_mean = avg(area_mean)
    delete(modeldata)
    delete(area_mean)

    if (time_mean.ge.1000) then
      fmt = tostring(toint(time_mean))
    else if (time_mean.lt.1000 .and. time_mean.ge.100) then
      fmt = sprintf("%5.1f", time_mean)
    else if (time_mean.lt.100 .and. time_mean.ge.10) then
      fmt = sprintf("%5.2f", time_mean)
    else if (time_mean.lt.10 .and. time_mean.ge.1) then
      fmt = sprintf("%5.3f", time_mean)
    else
      fmt = tostring(time_mean)
    end if
    end if
    end if
    end if

    ; Define output workstation
    wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring13f01_" \
                  + names(idx_mod) + "_" + var0)
    output_filename = wks@fullname

    ; Data to plot
    A3 = dim_avg_Wrap(var_grd(idx_mod, :, :, :))

    ; model specific resources
    res@tiMainFontHeightF = 0.03
    res@tiMainString = annots(idx_mod)
    res@gsnLeftString = var0 + " [" + variable_info[0]@units + "]"
    res@gsnRightString = fmt + " " + variable_info[0]@units
    res@gsnRightStringFontHeightF = 0.023
    res@gsnLeftStringFontHeightF = 0.023

    plotfield = A3(lat|:, month|:)

    ; Plotting
    if (write_plots) then
      map = gsn_csm_lat_time(wks, plotfield, res)
      log_info(" Wrote " + output_filename)
    end if

    nc_filename = work_dir + "eyring13f01_" + names(idx_mod) + "_" + \
      var0 + ".nc"

    plotfield@var = var0
    plotfield@diag_script = (/DIAG_SCRIPT/)

    nc_outfile = ncdf_write(plotfield, nc_filename)
    delete(plotfield)

    ; ------------------------------------------------------------------------
    ; write provenance to netcdf output and plot file
    ; ------------------------------------------------------------------------

    atts = True
    atts@dataset = names(idx_mod)
    info_dataset = select_metadata_by_atts(info0, atts)
    infiles = metadata_att_as_array(info_dataset, "filename")

    statistics = (/"clim", "mean"/)
    domain = "global"
    plottype = (/"seas", "zonal"/)
    contrib_authors = (/"cionni_irene", "lauer_axel"/)
    diag_ref = "eyring13jgr"

    caption = "Mean seasonal cycle of " + var0 \
              + " averaged between " + y1(idx_mod) + " and " + y2(idx_mod) \
              + ", similar to Eyring et al., JGR, 2013, fig. 1."

    log_provenance(nc_outfile, output_filename, caption, statistics, \
                   domain, plottype, contrib_authors, diag_ref, infiles)

    delete(caption)
    delete(info_dataset)
    delete(atts)
    delete(infiles)

    ; clean up
    delete([/wks, A3, time_mean/])
  end do  ; idx_MOD

  ; *********plot MMM plot list_chem plot list_nochem***********

  var_grd_mod = var_grd(ind_chem_mod, :, :, :)
  if (dimsizes(ind_chem_mod) .gt. 1) then
    MMM_chem = dim_avg_n_Wrap(var_grd_mod, 0)
  else
    MMM_chem = var_grd_mod
  end if
  delete(var_grd_mod)

  if (diag_script_info@e13fig01_diff_ref) then
    var_grd_ref = var_grd(ind_OBS_ref, :, :, :)
    diff_var_grd = MMM_chem - var_grd_ref
    copy_VarCoords(MMM_chem, diff_var_grd)
  end if

  if (num_nochem_mod .ge. 1) then
    var_grd_nochem_mod = var_grd(ind_nochem_mod, :, :, :)
    if (num_nochem_mod .gt. 1) then
      MMM_nochem = dim_avg_n_Wrap(var_grd_nochem_mod, 0)
    else
      MMM_nochem = var_grd_nochem_mod
    end if
    delete(var_grd_nochem_mod)
    if (diag_script_info@e13fig01_diff_ref) then
      diff_var_grd_nochem = MMM_nochem - var_grd_ref
      copy_VarCoords(MMM_nochem, diff_var_grd_nochem)
    end if
  end if

  ; **************plot panel**********

  ; ---------------------------------------------------------------------------
  ; filename for netCDF output (panel plot)

  nc_filename = work_dir + "eyring13f01_panel_" + var0 + ".nc"
  nc_filename@existing = "append"

  ; create indices for all datasets processed/plotted in the panel plot
  ; (0 = dataset not used, 1 = dataset used)

  indmask = new(dim_MOD, integer)
  indmask(:) = 0

  indmask(ind_chem_mod) = 1

  if (num_nochem_mod .ge. 1) then
    indmask(ind_nochem_mod) = 1
  end if

  if (ind_OBS_ref .ge. 0) then
    indmask(ind_OBS_ref) = 1
  end if

  ; create list of climofiles used for panel plot (for provenance)

  allind = ispan(0, dim_MOD - 1, 1)
  panelind = allind(ind(indmask .eq. 1))
  infiles = new(dimsizes(panelind), string)

  atts = True
  do ii = 0, dimsizes(panelind) - 1
    atts@dataset = names(panelind(ii))
    if (.not.ismissing(ensembles(panelind(ii)))) then
      atts@ensemble = ensembles(panelind(ii))
    else
      if (isatt(atts, "ensemble")) then
        delete(atts@ensemble)
      end if
    end if
    info_dataset = select_metadata_by_atts(info0, atts)
    infiles(ii) = metadata_att_as_array(info_dataset, "filename")
  end do

  delete(panelind)
  delete(indmask)
  delete(allind)
  ; ---------------------------------------------------------------------------

  ; calculate mean MMM

  modeldata = MMM_chem(:, :, :) * mask0
  copy_VarMeta(mask0, modeldata)
  area_mean = area_operations(modeldata, -90, 90, 0, 360, "average", True)
  time_mean = avg(area_mean)
  delete(modeldata)
  delete(area_mean)

  if (time_mean.ge.1000) then
    fmt = tostring(toint(time_mean))
  else if (time_mean.lt.1000 .and. time_mean.ge.100) then
    fmt = sprintf("%5.1f", time_mean)
  else if(time_mean.lt.100 .and. time_mean.ge.10) then
    fmt = sprintf("%5.2f", time_mean)
  else if (time_mean.lt.10 .and. time_mean.ge.1) then
    fmt = sprintf("%5.3f", time_mean)
  else
    fmt = tostring(time_mean)
  end if
  end if
  end if
  end if

  ; Define output workstation

  aux_title_info1 = "MMM_chem_" + OBS_ref + "_panel" + "_" + exps(0) + "_" + \
                    sprinti("%0.4i", y1(0)) + "-" + sprinti("%0.4i", y2(0))
  wks1 = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring13f01_" \
                 + aux_title_info1 + "_" + var0)
  output_filename1 = wks1@fullname

  ip = 0
  map_p = new(5, "graphic")

  if ((isatt(diag_script_info, "e13fig01_list_chem_mod_string"))) then
    res@tiMainString = diag_script_info@e13fig01_list_chem_mod_string
  else
    res@tiMainString = "MMM_chem"
  end if

  ; just one panel plot at the end
  res@gsnDraw  = False                          ; don't draw
  res@gsnFrame = False                          ; don't advance frame

  res@tiMainFontHeightF = 0.03
  res@gsnLeftString = var0 + " [" + variable_info[0]@units + "]"
  res@gsnRightString = fmt + " " + variable_info[0]@units
  res@gsnRightStringFontHeightF = 0.023
  res@gsnLeftStringFontHeightF = 0.023

  ; Plotting
  plotfield = dim_avg_Wrap(MMM_chem(lat|:, month|:, lon|:))
  map_p(ip) = gsn_csm_lat_time(wks1, plotfield, res)
  ip = ip + 1

  plotfield@var = var0 + "_MMMchem"
  plotfield@diag_script = (/DIAG_SCRIPT/)
  nc_outfile = ncdf_write(plotfield, nc_filename)

  ; clean up
  delete(time_mean)
  delete(plotfield)

  ; *****************NO CHEM ***************************

  if (num_nochem_mod .ge. 1) then
    modeldata = MMM_nochem(:, :, :) * mask0
    copy_VarMeta(mask0, modeldata)
    area_mean = area_operations(modeldata, -90, 90, 0, 360, "average", True)
    time_mean = avg(area_mean)
    delete(modeldata)
    delete(area_mean)

    if (time_mean.ge.1000) then
      fmt = tostring(toint(time_mean))
    else if (time_mean.lt.1000 .and. time_mean.ge.100) then
      fmt = sprintf("%5.1f", time_mean)
    else if (time_mean.lt.100 .and. time_mean.ge.10) then
      fmt = sprintf("%5.2f", time_mean)
    else if (time_mean.lt.10 .and. time_mean.ge.1) then
      fmt = sprintf("%5.3f", time_mean)
    else
      fmt = tostring(time_mean)
    end if
    end if
    end if
    end if

    if ((isatt(diag_script_info, "e13fig01_list_nochem_mod_string"))) then
      res@tiMainString = e13fig01_list_nochem_mod_string
    else
      res@tiMainString = "MMM_nochem"
    end if

    res@tiMainFontHeightF = 0.03
    res@gsnLeftString = var0 + " [" + variable_info[0]@units + "]"
    res@gsnRightString = fmt + " " + variable_info[0]@units
    res@gsnLeftStringFontHeightF = 0.023
    res@gsnRightStringFontHeightF = 0.023

    ; Plotting
    plotfield = dim_avg_Wrap(MMM_nochem(lat|:, month|:, lon|:))
    map_p(ip) = gsn_csm_lat_time(wks1, plotfield, res)
    ip = ip + 1

    plotfield@var = var0 + "_MMMnochem"
    plotfield@diag_script = (/DIAG_SCRIPT/)
    nc_outfile = ncdf_write(plotfield, nc_filename)

    ; clean up
    delete(time_mean)
    delete(plotfield)
  end if

  ; ****************************************************
  ; ******************REF*******************************

  if (OBS_ref .ne. "") then
    modeldata = var_grd_ref * mask0
    copy_VarMeta(mask0, modeldata)
    area_mean = area_operations(modeldata, -90, 90, 0, 360, "average", True)
    time_mean = avg(area_mean)
    delete(modeldata)
    delete(area_mean)

    if (time_mean.ge.1000) then
      fmt = tostring(toint(time_mean))
    else if (time_mean.lt.1000 .and. time_mean.ge.100) then
      fmt = sprintf("%5.1f", time_mean)
    else if (time_mean.lt.100 .and. time_mean.ge.10) then
      fmt = sprintf("%5.2f", time_mean)
    else if (time_mean.lt.10 .and. time_mean.ge.1) then
      fmt = sprintf("%5.3f", time_mean)
    else
      fmt = tostring(time_mean)
    end if
    end if
    end if
    end if

    res@tiMainString = OBS_ref
    res@tiMainFontHeightF = 0.03
    res@gsnLeftString = var0 + " [" + variable_info[0]@units + "]"
    res@gsnRightString = fmt + " " + variable_info[0]@units
    res@gsnRightStringFontHeightF = 0.023
    res@gsnLeftStringFontHeightF = 0.023
    res@cnLevelSelectionMode = "ExplicitLevels"
    delete(res@cnLevels)
    if (isatt(diag_script_info, "eyring13_reflvl")) then
      res@cnLevels = diag_script_info@eyring13_reflvl
    else
      res@cnLevels = fspan(min(var_grd_ref), max(var_grd_ref), 20)
    end if

    ; Plotting
    plotfield = dim_avg_Wrap(var_grd_ref(lat|:, month|:, lon|:))
    map_p(ip) = gsn_csm_lat_time(wks1, plotfield, res)
    ip = ip + 1

    plotfield@var = var0 + "_ref"
    plotfield@diag_script = (/DIAG_SCRIPT/)
    nc_outfile = ncdf_write(plotfield, nc_filename)

    ; clean up
    delete(time_mean)
    delete(plotfield)
  end if

  ; *****************diff**************************************

  if (diag_script_info@e13fig01_diff_ref) then
    modeldata = diff_var_grd * mask0
    copy_VarMeta(mask0, modeldata)
    area_mean = area_operations(modeldata, -90, 90, 0, 360, "average", True)
    time_mean = avg(area_mean)
    delete(modeldata)
    delete(area_mean)

    if (time_mean.ge.1000) then
      fmt = tostring(toint(time_mean))
    else if (time_mean.lt.1000 .and. time_mean.ge.100) then
      fmt = sprintf("%5.1f", time_mean)
    else if (time_mean.lt.100 .and. time_mean.ge.10) then
      fmt = sprintf("%5.2f", time_mean)
    else if (time_mean.lt.10 .and. time_mean.ge.1) then
      fmt = sprintf("%5.3f", time_mean)
    else
      fmt = tostring(time_mean)
    end if
    end if
    end if
    end if

    if ((isatt(diag_script_info, "e13fig01_list_chem_mod_string"))) then
      res@tiMainString = diag_script_info@e13fig01_list_chem_mod_string + \
                  " - " + OBS_ref
    else
      res@tiMainString = "MMM_chem - " + OBS_ref
    end if

    res@tiMainFontHeightF = 0.03
    res@gsnLeftString = var0 + " [" + variable_info[0]@units + "]"
    res@gsnRightString = fmt + " " + variable_info[0]@units
    res@gsnRightStringFontHeightF = 0.023
    res@gsnLeftStringFontHeightF = 0.023

    ; Plotting
    delete(res@cnLevelSelectionMode)
    delete(res@cnLevels)

    if (isatt(diag_script_info, "eyring13_difflvl")) then
      res@cnLevelSelectionMode = "ExplicitLevels"
      res@cnLevels = diag_script_info@eyring13_difflvl
    end if

    plotfield = dim_avg_Wrap(diff_var_grd(lat|:, month|:, lon|:))
    map_p(ip) = gsn_csm_lat_time(wks1, plotfield, res)
    ip = ip + 1

    plotfield@var = var0 + "MMMchem-ref"
    plotfield@diag_script = (/DIAG_SCRIPT/)
    nc_outfile = ncdf_write(plotfield, nc_filename)

    ; clean up
    delete(time_mean)
    delete(plotfield)

    if (num_nochem_mod .ge. 1) then
      modeldata = diff_var_grd_nochem * mask0
      copy_VarMeta(mask0, modeldata)
      area_mean = area_operations(modeldata, -90, 90, 0, 360, \
                                  "average", True)
      time_mean = avg(area_mean)
      delete(modeldata)
      delete(area_mean)

      if (time_mean.ge.1000) then
        fmt = tostring(toint(time_mean))
      else if (time_mean.lt.1000 .and. time_mean.ge.100) then
        fmt = sprintf("%5.1f", time_mean)
      else if (time_mean.lt.100 .and. time_mean.ge.10) then
        fmt = sprintf("%5.2f", time_mean)
      else if (time_mean.lt.10 .and. time_mean.ge.1) then
        fmt = sprintf("%5.3f", time_mean)
      else
        fmt = tostring(time_mean)
      end if
      end if
      end if
      end if

      if ((isatt(diag_script_info, "e13fig01_list_nochem_mod_string"))) then
        res@tiMainString = e13fig01_list_nochem_mod_string + " - " \
                           + OBS_ref
      else
        res@tiMainString = "MMM_nochem - " + OBS_ref
      end if

      res@tiMainFontHeightF = 0.03
      res@gsnLeftString = var0 + " [" + variable_info[0]@units + "]"
      res@gsnRightString = fmt + " " + variable_info[0]@units
      res@gsnRightStringFontHeightF = 0.023
      res@gsnLeftStringFontHeightF = 0.023

      ; Plotting
      delete(res@cnLevelSelectionMode)
      delete(res@cnLevels)

      if (isatt(diag_script_info, "eyring13_difflvl")) then
        res@cnLevelSelectionMode = "ExplicitLevels"
        res@cnLevels = diag_script_info@eyring13_difflvl
      end if

      plotfield = diff_var_grd_nochem(lat|:, month|:, lon|:)
      map_p(ip) = gsn_csm_lat_time(wks1, plotfield, res)
      ip = ip + 1

      plotfield@var = var0 + "MMMnochem-ref"
      plotfield@diag_script = (/DIAG_SCRIPT/)
      nc_outfile = ncdf_write(plotfield, nc_filename)

      ; clean up
      delete(time_mean)
      delete(plotfield)
    end if  ; if num_nochem_mod .ge. 1
  end if  ; if diag_script_info@e13fig01_diff_ref

  if (write_plots) then
    resP = True
    gsn_panel(wks1, map_p(0:ip - 1), (/1, ip/), resP)
    log_info(" Wrote " + output_filename1)
  end if

  ; ------------------------------------------------------------------------
  ; write provenance to netcdf output and plot file
  ; ------------------------------------------------------------------------

  caption = "Climatological mean seasonal cycle of " + var0 \
            + " averaged between " + y1(0) + " and " + y2(0) \
            + ", similar to Eyring et al., JGR, 2013, fig. 1."

  log_provenance(nc_outfile, output_filename, caption, statistics, \
                 domain, plottype, contrib_authors, diag_ref, infiles)

  delete(infiles)

  ; ********plot panel *****************************************************

  leave_msg(DIAG_SCRIPT, "")
end

