;##############################################################################
;# eyring13jgr_fig04.ncl                                                      #
;##############################################################################
;# Author: Irene Cionni (ENEA, Italy)                                         #
;# ESMVal project PI: Veronika Eyring (DLR, Germany)                          #
;##############################################################################
;#                                                                            #
;# Description:                                                               #
;#    * Calculates climatological annual mean tropospheric ozone columns      #
;#      and plots the geographical distribution.                              #
;#                                                                            #
;# Required diag_script_info attributes:                                      #
;#     none                                                                   #
;#                                                                            #
;# Optional diag_script_info attributes:                                      #
;#     acp_cnlvl: contour levels for plotting all datasets except ref dataset #
;#     anomalylevels: user-defined levels for plotting anomalies              #
;#     e13fig04_list_chem_mod: list of datasets in the group “chem” (array of #
;#         strings, default = all datasets except OBS, obs4mips, ref)         #
;#     e13fig04_list_chem_mod_string: plotting label for group “chem”,        #
;#         e.g., "Chem", default = "MMM_chem"                                 #
;#     e13fig04_list_nochem_mod: list of models in the group “nochem”         #
;#         (array of strings, default = (/""/))                               #
;#     e13fig04_list_nochem_mod_string: plotting label for group “nochem”,    #
;#         e.g., "NoChem", default = ""MMM_nochem"                            #
;#     font: overrides default font (e.g., 21, see                            #
;#         www.ncl.ucar.edu/Document/Graphics/Resources/tx.shtml#txFont)      #
;#     lbLabelBarOn: plot a label bar (True, False)                           #
;#     map_diff_Levels: contour levels for plotting differences               #
;#     map_ref_Levels: contour levels for plotting reference dataset          #
;#     mpProjection: map projection                                           #
;#     plot_units: units of variable used for plotting                        #
;#     range_option: 0 = as in recipe, 1 = overlapping time period            #
;#     subtract_global_mean: subtract global mean before plotting             #
;#                                                                            #
;# Required variable_info attributes:                                         #
;#     reference_dataset: name of reference dataset (observations)            #
;#                                                                            #
;# Caveats:                                                                   #
;#                                                                            #
;# Modification history:                                                      #
;#    20180409-A_laue_ax: code rewritten for v2, regridding now done by       #
;#                        backend                                             #
;#    20170811-A_bock_ls: added tags for reporting                            #
;#    20160916-A-laue_ax: added option to subtract global mean before         #
;#                        plotting (--> anomalies)                            #
;#    20151208-A-laue_ax: code clean-up, moved call to write_references to    #
;#                        beginning of code                                   #
;#    20140528-A_GoKl: written based on ancyc_lat.ncl                         #
;#                                                                            #
;##############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/contour_maps.ncl"


;; fetch general parameters
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

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Set default values for optional diag_script_info attributes

  set_default_att(diag_script_info, "font", 21)
  set_default_att(diag_script_info, "lbLabelBarOn", True)
  set_default_att(diag_script_info, "mpProjection", "CylindricalEquidistant")
  set_default_att(diag_script_info, "range_option", 0)
  set_default_att(diag_script_info, "subtract_global_mean", False)

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

;; main code
begin

  subtract_global_mean = diag_script_info@subtract_global_mean

  if (isvar("refname")) then
    OBS_ref = refname
    ind_OBS_ref = ind(names.eq.OBS_ref)
  else
    error_msg("f", DIAG_SCRIPT, "", "no reference dataset defined")
  end if

  if (ismissing(ind_OBS_ref)) then
    error_msg("f", DIAG_SCRIPT, "", "reference dataset not found")
  end if

  if ((isatt(diag_script_info, "anomalylevels"))) then
    anomalylevels = diag_script_info@anomalylevels
  else
    anomalylevels = fspan(-10., 10., 21)
  end if

  if ((isatt(diag_script_info,"e13fig04_list_chem_mod"))) then
    list_chem_mod = diag_script_info@e13fig04_list_chem_mod
    ind_chem_mod = ind(models@name.eq.list_chem_mod)
  else
    tmp = ind(projects.ne."OBS" .and. projects.ne."obs4mips")
    ind_chem_mod = ind(tmp .ne. ind_OBS_ref)
    list_chem_mod = names(ind_chem_mod)
    delete(tmp)
  end if

  if ((isatt(diag_script_info, "e13fig04_list_nochem_mod"))) then
    list_nochem_mod = diag_script_info@e13fig04_list_nochem_mod
    ind_nochem_mod = ind(names.eq.list_nochem_mod)
  end if

  if ((isatt(diag_script_info, "e13fig04_diff_ref"))) then
    log_info("e13fig04_diff_ref is no longer used and will be irgnored. " \
             + "The reference dataset is now defined via the variable " \
             + "attribute 'reference_dataset' (ESMValTool v2 standard). ")
  end if

  ;*************************************************

  ;; default time range & annots
  y1 = metadata_att_as_array(info0, "start_year")
  y2 = metadata_att_as_array(info0, "end_year")
  annots = project_style(info0, diag_script_info, "annots")

  range_option = diag_script_info@range_option
  if ((range_option .ne. 0) .and. (range_option .ne. 1)) then
    error_msg("f", DIAG_SCRIPT, "", \
              "invalid range_option; valid settings: 0 or 1")
  end if

  ;; Use overlapping time period or as given in nml?
  if (range_option .eq. 1) then
    y1 = max(toint(y1))
    y2 = min(toint(y2))
  end if

  ;; Loop over models and create mask for missing values
  do idx_mod = 0, dim_MOD - 1
    ;; Extract data
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

    ;; annual mean
    A1 = time_operations(temp, -1, -1, "average", "annualclim", True)
    delete(temp)

    ;; initialize mask...
    if (idx_mod.eq.0) then
      if (.not.isatt(A1, "_FillValue")) then
        Fill = default_fillvalue(typeof(A1))
      else
        Fill = A1@_FillValue
      end if

      var_grd = new((/dim_MOD, dimsizes(A1&lat), dimsizes(A1&lon)/), \
                    typeof(A1), Fill)
      var_grd!0 = "models"
      var_grd!1 = "lat"
      var_grd!2 = "lon"
      var_grd&models = annots
      var_grd&lat = A1&lat
      var_grd&lon = A1&lon

      mask0 = A1
      mask0(:, :) = 1

      var_grd@_FillValue = Fill  ; ensure consistent _FillValue
      delete(dims)
    end if

    var_grd(idx_mod, :, :) = where(ismissing(A1), Fill, A1)

    ;; Build the mask: missing in at least one model = missing in the mask
    mask0 = where(ismissing(A1), mask0@_FillValue, mask0)

    ;; clean up
    delete(A1)
  end do  ; idx_MOD

  ;; Convert units if requested
  if (isatt(diag_script_info, "plot_units")) then
    var_grd = convert_units(var_grd, diag_script_info@plot_units)
  end if

  ;; General settings for this plot
  res = True
  res@cnFillOn = True
  res@cnLinesOn = False
  res@cnLevelSelectionMode = "ExplicitLevels"
  if (isatt(diag_script_info, "acp_cnlvl")) then
    res@cnLevels = diag_script_info@acp_cnlvl
  else
    res@cnLevels = fspan(min(var_grd), max(var_grd), 20)
  end if

  res@mpProjection = diag_script_info@mpProjection
  if (res@mpProjection.eq."Robinson") then
    res@mpPerimOn = False  ; turn off perimeter around map
    res@mpGridLineColor = -1
    res@mpGridAndLimbOn = True
  end if

  res@txFont = diag_script_info@font
  res@tmXBLabelFont = diag_script_info@font
  res@tmYLLabelFont = diag_script_info@font
  res@lbLabelBarOn = diag_script_info@lbLabelBarOn
  res@lbLabelFont = diag_script_info@font
  res@lbTitleFont = diag_script_info@font
  res@lbTitleFontHeightF = 0.01  ; default: 0.025
  res@lbTitleOffsetF = -0.25  ; default: 0.03 (negative values move it down)
  res@lbTopMarginF = 0.1  ; default: 0.05 (negative values move it up)
  res@lbTitleJust = "BottomCenter"  ; default: CenterCenter
  res@lbLabelFontHeightF = 0.025

  ;; Loop over datasets for plotting
  do idx_mod = 0, dim_MOD - 1
    ;; calculate mean
    A2 = where(ismissing(mask0), mask0, var_grd(idx_mod, :, :))
    copy_VarMeta(mask0, A2)
    time_mean = area_operations(A2, -90, 90, 0, 360, "average", True)
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

    ;; Define output workstation

    wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring13f04_" \
                  + names(idx_mod) + "_" + var0)
    output_filename = wks@fullname

    ;; Data to plot
    A3 = var_grd(idx_mod, :, :)

    ;; model specific ressources
    res@tiMainFontHeightF = 0.03
    res@tiMainString = annots(idx_mod)
    res@gsnLeftString = var0 + " [" + variable_info[0]@units + "]"
    res@gsnRightString = fmt + " " + variable_info[0]@units
    res@gsnRightStringFontHeightF = 0.023
    res@gsnLeftStringFontHeightF = 0.023

    ;; Plotting
    if (subtract_global_mean) then
      plotfield = A3 - area_operations(A3, -90, 90, 0, 360, "average", True)
      delete(res@cnLevels)
      res@cnLevels = anomalylevels
    else
       plotfield = A3
    end if
    copy_VarMeta(A3, plotfield)

    if (write_plots) then
      map = gsn_csm_contour_map(wks, plotfield, res)
      log_info(" Wrote " + output_filename)
    end if

    ; write mandatory netCDF

    plotfield@var = var0
    plotfield@diag_script = (/DIAG_SCRIPT/)

    nc_filename = work_dir + "eyring13f04_" + names(idx_mod) + "_" \
                  + var0 + ".nc"
    nc_outfile = ncdf_write(plotfield, nc_filename)
    delete(plotfield)

    ; ------------------------------------------------------------------------
    ; write provenance to netcdf output and plot file
    ; ------------------------------------------------------------------------

    statistics = (/"clim", "mean"/)
    domain = "glob"
    plottype = "geo"
    if (subtract_global_mean) then
      caption = "Deviation of " + var0 + " from global annual mean" \
                + " (averaged between " + y1(idx_mod) + " and " + y2(idx_mod) + ")"
    else
      caption = "Climatological annual mean of " + var0 \
                + " averaged between " + y1(idx_mod) + " and " + y2(idx_mod) \
                + ", similar to Eyring et al., JGR, 2013, fig. 4."
    end if

    atts = True
    atts@dataset = names(idx_mod)
    info_dataset = select_metadata_by_atts(info0, atts)
    infiles = metadata_att_as_array(info_dataset, "filename")

    contrib_authors = (/"cion_ir", "laue_ax"/)

print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
print(nc_outfile)
print(output_filename)
print(infiles)
print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

    log_provenance(nc_outfile, output_filename, caption, statistics, \
                   domain, plottype, contrib_authors, "", infiles)

    delete(infiles)

    ;; clean up
    delete([/wks, A2, A3, time_mean/])
  end do  ; idx_MOD, for plotting

  ; ---------------------------------------------------------------------------
  ; filename for netCDF output (panel plot)

  nc_filename = work_dir + "eyring13f04_panel_" + var0 + ".nc"
  nc_filename@existing = "append"

  ; create indices for all datasets processed/plotted in the panel plot
  ; (0 = dataset not used, 1 = dataset used)

  indmask = new(dim_MOD, integer)
  indmask(:) = 0

  indmask(ind_chem_mod) = 1
  if ((isatt(diag_script_info, "e13fig04_list_nochem_mod"))) then
    if (dimsizes(ind_nochem_mod) .ge. 1) then
      indmask(ind_nochem_mod) = 1
    end if
  end if
  indmask(ind_OBS_ref) = 1

  ; create list of climofiles used for panel plot (for provenance)

  allind = ispan(0, dim_MOD - 1, 1)
  panelind = allind(ind(indmask .eq. 1))
  infiles = new(dimsizes(panelind), string)

  atts = True
  do ii = 0, dimsizes(panelind) - 1
    atts@dataset = names(panelind(ii))
    info_dataset = select_metadata_by_atts(info0, atts)
    infiles(ii) = metadata_att_as_array(info_dataset, "filename")
  end do

  delete(panelind)
  delete(indmask)
  delete(allind)

  ; ---------------------------------------------------------------------------

  ;*********plot MMM plot lista_chem plot lista_nochem***********

  var_grd_mod = var_grd(ind_chem_mod, :, :)
  if (dimsizes(ind_chem_mod) .gt. 1) then
    MMM_chem = dim_avg_n_Wrap(var_grd_mod, 0)
  else
    MMM_chem = var_grd_mod
  end if
  delete(var_grd_mod)

  var_grd_ref = var_grd(ind_OBS_ref, :, :)

  if ((isatt(diag_script_info, "e13fig04_list_nochem_mod"))) then
    if (dimsizes(ind_nochem_mod) .ge. 1) then
      var_grd_nochem_mod = var_grd(ind_nochem_mod, :, :)
      if (dimsizes(ind_nochem_mod) .gt. 1) then
        MMM_nochem = dim_avg_n_Wrap(var_grd_nochem_mod, 0)
      else
        MMM_nochem = var_grd_nochem_mod
      end if
      delete(var_grd_nochem_mod)
      diff_var_grd_nochem = MMM_nochem - var_grd_ref
    end if
  end if

  diff_var_grd = MMM_chem - var_grd_ref

  ;**************plot add**********

  ;; calculate mean MMM
  MMM_chem = where(ismissing(mask0), mask0, MMM_chem)
  copy_VarMeta(mask0, MMM_chem)
  time_mean = area_operations(MMM_chem, -90, 90, 0, 360, "average", True)

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

  ;; Define output workstation

  aux_title_info1 = "MMM_panel" + "_" + exps(0) + "_" + \
                    sprinti("%0.4i", y1(0)) + "-" + sprinti("%0.4i", y2(0))

  wks1 = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring13f04_" \
                 + aux_title_info1 + "_" + var0)
  output_filename1 = wks1@fullname

  ip = 0
  map_p = new(5, "graphic")

  if ((isatt(diag_script_info, "e13fig04_list_chem_mod_string"))) then
    res@tiMainString = diag_script_info@e13fig04_list_chem_mod_string
  else
    res@tiMainString = "MMM_chem"
  end if

  ;just one panel plot at the end
  res@gsnDraw  = False                          ; don't draw
  res@gsnFrame = False                          ; don't advance frame

  res@tiMainFontHeightF = 0.03
  res@gsnLeftString = var0 + " [" + variable_info[0]@units + "]"
  res@gsnRightString = fmt + " " + variable_info[0]@units
  res@gsnRightStringFontHeightF = 0.023
  res@gsnLeftStringFontHeightF = 0.023
  delete(res@cnLevels)
  if (isatt(diag_script_info, "acp_cnlvl")) then
    res@cnLevels = diag_script_info@acp_cnlvl
  else
    res@cnLevels = fspan(min(MMM_chem), max(MMM_chem), 20)
  end if

  ;; Plotting
  if (subtract_global_mean) then
    plotfield = MMM_chem - area_operations(MMM_chem, -90, 90, 0, 360, \
                "average", True)
    delete(res@cnLevels)
    res@cnLevels = anomalylevels
  else
    plotfield = MMM_chem
  end if
  copy_VarMeta(MMM_chem, plotfield)
  map_p(ip) = gsn_csm_contour_map(wks1, plotfield, res)

  ip = ip + 1

  plotfield@var = var0 + "_MMMchem"
  plotfield@diag_script = (/DIAG_SCRIPT/)
  nc_outfile = ncdf_write(plotfield, nc_filename)

  ;; clean up
  delete(time_mean)
  delete(plotfield)

  ;*****************NO CHEM ***************************

  if ((isatt(diag_script_info, "e13fig04_list_nochem_mod"))) then
    if (dimsizes(ind_nochem_mod) .ge. 1) then
      MMM_nochem = where(ismissing(mask0), mask0, MMM_nochem)
      copy_VarMeta(mask0, MMM_nochem)
      time_mean = area_operations(MMM_nochem, -90, 90, 0, 360, "average", True)
      if (time_mean.ge.1000) then
        fmt = tostring(toint(time_mean))
      else if (time_mean.lt.1000 .and. time_mean.ge.100) then
        fmt = sprintf("%5.1f", time_mean)
      else if(time_mean.lt.100 .and. time_mean.ge.10) then
        fmt = sprintf("%5.2f", time_mean)
      else if(time_mean.lt.10 .and. time_mean.ge.1) then
        fmt = sprintf("%5.3f", time_mean)
      else
        fmt = tostring(time_mean)
      end if
      end if
      end if
      end if

      if ((isatt(diag_script_info, "e13fig04_list_nochem_mod_string"))) then
        res@tiMainString = e13fig04_list_nochem_mod_string
      else
        res@tiMainString = "MMM_nochem"
      end if

      res@tiMainFontHeightF = 0.03
      res@gsnLeftString = var0 + " [" + variable_info[0]@units + "]"
      res@gsnRightString = fmt + " " + variable_info[0]@units
      res@gsnLeftStringFontHeightF = 0.023
      res@gsnRightStringFontHeightF = 0.023
      delete(res@cnLevels)
      if (isatt(diag_script_info, "acp_cnlvl")) then
        res@cnLevels = diag_script_info@acp_cnlvl
      else
        res@cnLevels = fspan(min(MMM_nochem), max(MMM_nochem), 20)
      end if

      ;; Plotting
      if (subtract_global_mean) then
        plotfield = MMM_nochem - area_operations(MMM_nochem, -90, 90, \
                                                 0, 360, "average", True)
        delete(res@cnLevels)
        res@cnLevels = anomalylevels
      else
        plotfield = MMM_nochem
      end if
      copy_VarMeta(MMM_nochem, plotfield)
      map_p(ip) = gsn_csm_contour_map(wks1, plotfield, res)

      ip = ip + 1

      plotfield@var = var0 + "_MMMnochem"
      plotfield@diag_script = (/DIAG_SCRIPT/)
      nc_outfile = ncdf_write(plotfield, nc_filename)

      ;; clean up
      delete(time_mean)
      delete(plotfield)
    end if
  end if

  ;****************************************************
  ;******************REF*******************************

  var_grd_ref = where(ismissing(mask0), mask0, var_grd_ref)
  copy_VarMeta(mask0, var_grd_ref)
  time_mean = area_operations(var_grd_ref, -90, 90, 0, 360, "average", True)
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

  aux_title_info = refname
  res@tiMainString = refname
  res@tiMainFontHeightF = 0.03
  res@gsnLeftString = var0 + " [" + variable_info[0]@units + "]"
  res@gsnRightString = fmt + " " + variable_info[0]@units
  res@gsnRightStringFontHeightF = 0.023
  res@gsnLeftStringFontHeightF = 0.023
  delete(res@cnLevels)
  if (isatt(diag_script_info, "map_ref_Levels")) then
    res@cnLevels = diag_script_info@map_ref_Levels
  else
    res@cnLevels = fspan(min(var_grd_ref), max(var_grd_ref), 20)
  end if

  ;; Plotting
  if (subtract_global_mean) then
    plotfield = var_grd_ref - area_operations(var_grd_ref, -90, 90, \
                                              0, 360, "average", True)
    delete(res@cnLevels)
    res@cnLevels = anomalylevels
  else
    plotfield = var_grd_ref
  end if
  copy_VarMeta(var_grd_ref, plotfield)
  map_p(ip) = gsn_csm_contour_map(wks1, plotfield, res)

  ip = ip + 1

  plotfield@var = var0 + "_ref"
  plotfield@diag_script = (/DIAG_SCRIPT/)
  nc_outfile = ncdf_write(plotfield, nc_filename)

  ;; clean up
  delete(time_mean)
  delete(plotfield)

  ;*****************diff**************************************

  diff_var_grd = where(ismissing(mask0), mask0, diff_var_grd)
  copy_VarMeta(mask0, diff_var_grd)
  time_mean = area_operations(diff_var_grd, -90, 90, 0, 360, "average", True)

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

  if ((isatt(diag_script_info, "e13fig04_list_chem_mod_string"))) then
    res@tiMainString = diag_script_info@e13fig04_list_chem_mod_string \
                       + " - " + refname
  else
    res@tiMainString = "MMM_chem"
  end if

  res@tiMainFontHeightF = 0.03
  res@gsnLeftString = var0 + " [" + variable_info[0]@units + "]"
  res@gsnRightString = fmt + " " + variable_info[0]@units
  res@gsnRightStringFontHeightF = 0.023
  res@gsnLeftStringFontHeightF = 0.023

  ;; Plotting
  delete(res@cnLevels)
  if (isatt(diag_script_info, "map_diff_Levels")) then
    res@cnLevels = diag_script_info@map_diff_Levels
  else
    res@cnLevels = fspan(min(diff_var_grd), max(diff_var_grd), 20)
  end if
  map_p(ip) = gsn_csm_contour_map(wks1, diff_var_grd, res)

  ip = ip + 1

  diff_var_grd@var = var0 + "MMMchem-ref"
  diff_var_grd@diag_script = (/DIAG_SCRIPT/)
  nc_outfile = ncdf_write(diff_var_grd, nc_filename)

  ;; clean up
  delete(time_mean)
  delete(diff_var_grd)

  if ((isatt(diag_script_info, "e13fig04_list_nochem_mod"))) then
    if (dimsizes(ind_nochem_mod) .ge. 1) then
      diff_var_grd_nochem = where(ismissing(mask0), mask0, diff_var_grd_nochem)
      copy_VarMeta(mask0, diff_var_grd_nochem)
      time_mean = area_operations(diff_var_grd_nochem, -90, 90, 0, 360, \
                                  "average", True)
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

      if ((isatt(diag_script_info, "e13fig04_list_nochem_mod_string"))) then
        res@tiMainString = e13fig04_list_nochem_mod_string + " - " + refname
      else
        res@tiMainString = "MMM_nochem"
      end if

      res@tiMainFontHeightF = 0.03
      res@gsnLeftString = var0 + " [" + variable_info[0]@units + "]"
      res@gsnRightString = fmt + " " + variable_info[0]@units
      res@gsnRightStringFontHeightF = 0.023
      res@gsnLeftStringFontHeightF = 0.023

      ;; Plotting
      delete(res@cnLevels)
      if (isatt(diag_script_info, "map_diff_Levels")) then
        res@cnLevels = diag_script_info@map_diff_Levels
      else
        res@cnLevels = fspan(min(diff_var_grd_nochem), \
                             max(diff_var_grd_nochem), 20)
      end if

      map_p(ip) = gsn_csm_contour_map(wks1, diff_var_grd_nochem, res)
      ip = ip + 1

      diff_var_grd_nochem@var = var0 + "MMMnochem-ref"
      diff_var_grd_nochem@diag_script = (/DIAG_SCRIPT/)
      nc_outfile = ncdf_write(diff_var_grd_nochem, nc_filename)

      ;; clean up
      delete(time_mean)
      delete(diff_var_grd_nochem)
    end if
  end if

  if (write_plots) then
    resP = True
    gsn_panel(wks1, map_p(0:ip - 1), (/1, ip/), resP)
    log_info(" Wrote " + output_filename1)
  end if

  ; ------------------------------------------------------------------------
  ; write provenance to netcdf output and plot file
  ; ------------------------------------------------------------------------

  statistics = (/"clim", "mean"/)
  domain = "glob"
  plottype = "geo"
  if (subtract_global_mean) then
    caption = "Multi-model deviation of " + var0 + " from global annual mean" \
              + " (averaged between " + y1(0) + " and " + y2(0) + ")"
  else
    caption = "Multi-model climatological annual mean of " + var0 \
              + " averaged between " + y1(0) + " and " + y2(0) \
              + ", similar to Eyring et al., JGR, 2013, fig. 4."
  end if

  contrib_authors = (/"cion_ir", "laue_ax"/)

  log_provenance(nc_outfile, output_filename, caption, statistics, \
                 domain, plottype, contrib_authors, "", infiles)

  delete(infiles)

  ;***********************************************************************
  ;********plot diff *****************************************************

  leave_msg(DIAG_SCRIPT, "")
end
