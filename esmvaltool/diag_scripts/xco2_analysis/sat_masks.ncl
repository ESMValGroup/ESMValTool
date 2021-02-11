; #############################################################################
;  xco2_analysis/sat_masks.ncl
; #############################################################################
; Description
;    Plotting maps for missing values in satellite data.
;
; Required diag_script_info attributes:
;
; Optional diag_script_info attributes:
;    output_file_type: output file type for plots. Default: png
;    var_plotname: String formatting how variable should be named in plots
;                  defaults to short_name if not assigned
;    c3s_plots: Missing value plots seperated by timeseries of c3s satellites
;
; Caveats
;
; Modification history
;    20201119-gier_bettina: Added provenance and cleaned up code.
;    20200226-gier_bettina: Adapted to version 2
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/scaling.ncl"
load "$diag_scripts/shared/set_operators.ncl"
load "$diag_scripts/shared/statistics.ncl"

load "$diag_scripts/shared/plot/scatterplot.ncl"
load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/xy_line.ncl"

load "$diag_scripts/xco2_analysis/stat.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")
  AUTHORS = (/"gier_bettina"/)
  REFERENCES = (/"gier20bg"/)
  DOMAIN = (/"global"/)

  ; Variable
  var0 = variable_info[0]@short_name

  ; Input data
  INFO0 = select_metadata_by_name(input_file_info, var0)
  DATASETS = metadata_att_as_array(INFO0, "dataset")
  DATASETS := array_append_record(DATASETS, "multi-model mean", 0)
  experiments = metadata_att_as_array(INFO0, "exp")
  dim_MOD = dimsizes(DATASETS)
  ALL_FILES = metadata_att_as_array(INFO0, "filename")

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")
end

begin

  ; Plot file type
  if (isatt(diag_script_info, "output_file_type")) then
    file_type = diag_script_info@output_file_type
  elseif (isatt(config_user_info, "output_file_type")) then
    file_type = config_user_info@output_file_type
  else
    file_type = "png"
  end if

  ; Output plot directory
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)

  ; Formatted varname for plots
  if (isatt(diag_script_info, "var_plotname")) then
      var0_plotname = diag_script_info@var_plotname
  else
      var0_plotname = var0
  end if

  ; Formatted varname for plots
  if (isatt(diag_script_info, "c3s_plots")) then
      c3s_mask = True
  else
      c3s_mask = False
  end if

end

begin
  ; Prepare arrays

  ; Make array to collect timeseries data
  start_years = metadata_att_as_array(INFO0, "start_year")
  start_year = min(start_years)
  end_years = metadata_att_as_array(INFO0, "end_year")
  end_year = max(end_years)
  nyear = (end_year - start_year + 1)
  time = new(nyear*12, float)
  do yy = start_year, end_year
    do mm = 0, 11
      time(12 * (yy - start_year) + mm) = yy + (mm + 0.5)/12.
      end do
  end do

end

begin
  ; First read observations for observational masking
  if (isatt(INFO0[0], "reference_dataset")) then
    ref_ind = ind(DATASETS.eq.INFO0[0]@reference_dataset)
  else
    error_msg("w", DIAG_SCRIPT, "", \
              "this diagnostic requires a reference dataset")
  end if

  obs_data =  read_data(INFO0[ref_ind])

  ; ---------------------------------------------
  ; Determine missing value masks
  ; ----------------- OVERALL MASKS
  time_mnth = cd_calendar(obs_data&time, 0)

  global_mask = obs_data(0, :, :)
  global_mask = 0.

  monthly_mask = conform_dims((/12, dimsizes(global_mask(:, 0)), \
                               dimsizes(global_mask(0, :))/), global_mask, \
                              (/1, 2/))
  monthly_mask!0 = "month"
  monthly_mask&month = (/"January", "February", "March", "April", "May", \
                         "June", "July", "August", "September", "October",\
                         "November", "December"/)
  monthly_mask!1 = "lat"
  monthly_mask&lat = global_mask&lat
  monthly_mask!2 = "lon"
  monthly_mask&lon = global_mask&lon
  monthly_mask@var = "data_coverage_percent"
  monthly_mask@diag_script = DIAG_SCRIPT

  copy_VarAtts(obs_data, monthly_mask)

  do i_tim = 0, dimsizes(obs_data&time) - 1
      imon = toint(time_mnth(i_tim, 1))
      mnth_mv = where(ismissing(obs_data(i_tim, :, :)), 1, 0)
      global_mask = global_mask + mnth_mv
      monthly_mask(imon-1, :, :) = monthly_mask(imon-1, :, :) + mnth_mv
  end do
  global_mask = 1 - (global_mask / tofloat(dimsizes(obs_data&time)))
  monthly_mask = 1 - (monthly_mask / (dimsizes(obs_data&time)/12.))

  ; --------------------------------------------------------------
  ; ------ PER SATELLITE MASKS -> C3S special masks
  ; --------------------------------------------------------------

  if c3s_mask then

    ; Calculate spatial mask approximations

    ; Sciamachy Masks! 2003-2008
    scia_ind = max(ind(time_mnth(:, 0).eq.2008))
    scia_l_mask = global_mask
    scia_l_mask = 0.

    do i_tim = 0, dimsizes(obs_data&time(:scia_ind)) - 1
        mnth_mv = where(ismissing(obs_data(i_tim, :, :)), 0, 1)
        scia_l_mask = scia_l_mask + mnth_mv
    end do
    scia_l_mask = scia_l_mask / tofloat(dimsizes(obs_data&time(:scia_ind)))

    ; Gosat Masks! 2013 - .. (overlap 2009 to 2012)
    gosat_ind = min(ind(time_mnth(:, 0).eq.2013))
    gosat_l_mask = global_mask
    gosat_l_mask = 0.

    do i_tim = gosat_ind, dimsizes(obs_data&time) - 1
        mnth_mv = where(ismissing(obs_data(i_tim, :, :)), 0, 1)
        gosat_l_mask = gosat_l_mask + mnth_mv
    end do
    gosat_l_mask = gosat_l_mask / tofloat(dimsizes(obs_data&time(gosat_ind:)))

    ; SCIA + GOSAT joint mask
    sgj_i_ind = min(ind(time_mnth(:, 0).eq.2009))
    sgj_f_ind = max(ind(time_mnth(:, 0).eq.2012))
    sgj_l_mask = global_mask
    sgj_l_mask = 0.

    do i_tim = sgj_i_ind, sgj_f_ind
        mnth_mv = where(ismissing(obs_data(i_tim, :, :)), 0, 1)
        sgj_l_mask = sgj_l_mask + mnth_mv
    end do
    sgj_l_mask = sgj_l_mask / \
      tofloat(dimsizes(obs_data&time(sgj_i_ind:sgj_f_ind)))

    ; Add attributes for netcdf files
  scia_l_mask@var = "data_coverage_percent_scia"
  scia_l_mask@diag_script = DIAG_SCRIPT
  gosat_l_mask@var = "data_coverage_percent_gosat"
  gosat_l_mask@diag_script = DIAG_SCRIPT
  sgj_l_mask@var = "data_coverage_percent_overlap"
  sgj_l_mask@diag_script = DIAG_SCRIPT

  end if  ; end C3S mask special

  ; ------------------------ PLOTS --------------------------------------------
  ; -------------- Missing value Plots ----------------------------------------
  mv_cmap = read_colormap_file("WhiteGreen")

  ; Plot missing values per month
  filename = var0 + "_missing_values_monthly_" + start_year + "-" + end_year
  outfile = plot_dir + filename
  outfile_netcdf = work_dir + "/" + filename + ".nc"
  ncdf_outfile = ncdf_write(monthly_mask, outfile_netcdf)

  wks = gsn_open_wks(file_type, outfile)
  mres = True
  mres@cnFillOn = True
  mres@cnLinesOn = False
  mres@gsnDraw = False
  mres@gsnFrame = False
  mres@lbLabelBarOn = False
  mres@gsnRightString = ""
  mres@gsnStringFontHeightF = 0.03
  mres@tmXBLabelFontHeightF = 0.025
  mres@tmYLLabelFontHeightF = 0.025
  mres@tmLabelAutoStride = True
  mres@cnConstFLabelOn = False
  mres@cnFillPalette    = mv_cmap   ; set color map

  mon_plots = new(12, graphic)

  do i = 0, 11
      mres@gsnLeftString = monthly_mask&month(i)
      mon_plots(i) = gsn_csm_contour_map(wks, monthly_mask(i, :, :), mres)
  end do

  pres = True
  pres@gsnMaximize            = True    ; Maximize in frame
  pres@gsnFrame               = False   ; Don't advance frame
  pres@gsnPanelLabelBar = True
  pres@lbTitleOn = True
  pres@lbTitlePosition = "Bottom"
  pres@lbTitleString = var0_plotname + " Data Coverage"
  pres@lbLabelFontHeightF = 0.012
  pres@lbTitleFontHeightF = 0.015
  pres@lbBottomMarginF = -1
  pres@lbTopMarginF = 0.7
  pres@lbTitleOffsetF = 0.2
  pres@gsnPanelYWhiteSpacePercent = 2.

  gsn_panel(wks, mon_plots, (/4, 3/), pres)

  delete(pres)
  delete(mres)
  frame(wks)

  ; Provenance
  log_provenance(ncdf_outfile, \
                 outfile, \
                 var0 + "fractional data coverage " + start_year + "-" \
                 + end_year + ". Similar to Gier et al 2020, Fig 1.", \
                 (/"mean"/), \
                 DOMAIN, \
                 "histogram", \
                 AUTHORS, \
                 REFERENCES, \
                 ALL_FILES)

  ; -------------------------------------------------------------------------
  ; --------------------Panel Plot of spatial masks

  if c3s_mask then
    ; Plot spatial mask scia
    filename_c3s = var0 + "_spatial_mask_panels"
    outfile_c3s = plot_dir + "/" + filename_c3s
    outfile_netcdf_c3s = work_dir + "/" + filename_c3s + ".nc"

    ; Write data to netcdf
    outfile_netcdf_c3s@existing = "overwrite"
    ncdf_outfile = ncdf_write(scia_l_mask, outfile_netcdf_c3s)
    outfile_netcdf_c3s@existing = "append"
    ncdf_outfile = ncdf_write(gosat_l_mask, outfile_netcdf_c3s)
    ncdf_outfile = ncdf_write(sgj_l_mask, outfile_netcdf_c3s)

    wks = gsn_open_wks(file_type, outfile_c3s)
    res = True
    res@lbLabelBarOn        = False
    res@cnFillOn = True
    res@cnLinesOn = False
    res@gsnDraw = False
    res@gsnFrame = False
    res@gsnLeftString = "SCIAMACHY"
    res@gsnRightString = "2003 - 2008"
    res@tiMainString = ""
    res@gsnStringFontHeightF = 0.02
    res@tmXBLabelFontHeightF = 0.02
    res@tmYLLabelFontHeightF = 0.02
    res@cnConstFLabelOn = False
    res@cnLineLabelsOn      = False
    res@tmLabelAutoStride = True
    mv_cmap = read_colormap_file("WhiteGreen")
    res@cnFillPalette    = mv_cmap
    plot_scia = gsn_csm_contour_map(wks, scia_l_mask, res)
    ; Add polylines for 30 + 60N
    lnres = True
    lnres@gsLineThicknessF = 2.0
    lnid30 = gsn_add_polyline(wks, plot_scia, (/0, 120, 240, 360/), \
                              (/30, 30, 30, 30/), lnres)
    lnid60 = gsn_add_polyline(wks, plot_scia, (/0, 120, 240, 360/), \
                              (/60, 60, 60, 60/), lnres)

    ; Plot for stippling

    res2 = True                            ; res2 probability plots
    res2@lbLabelBarOn        = False
    res2@gsnDraw             = False       ; Do not draw plot
    res2@gsnFrame            = False       ; Do not advance frame
    res2@cnInfoLabelOn       = False
    res2@cnLinesOn           = False       ; do not draw contour lines
    res2@cnLineLabelsOn      = False       ; do not draw contour labels
    res2@cnConstFLabelOn = False
    res2@cnFillColor = "black"

    res2@gsnLeftString       = ""
    res2@gsnRightString      = ""

    plot_scias = gsn_csm_contour(wks, scia_l_mask, res2)
    opt     = True
    opt@gsnShadeFillType = "pattern"
    opt@gsnShadeHigh = 12
    opt@gsnShadeFillScaleF = 0.6
    ; shade all areas with less than 50% missing values
    plot_scias = gsn_contour_shade(plot_scias, 0.00, 0.50, opt)
    overlay(plot_scia, plot_scias)

    res@gsnLeftString = "SCIAMACHY + GOSAT"
    res@gsnRightString = "2009 - 2012"
    plot_sgj = gsn_csm_contour_map(wks, sgj_l_mask, res)
    lnid30sgj = gsn_add_polyline(wks, plot_sgj, (/0, 120, 240, 360/), \
                                 (/30, 30, 30, 30/), lnres)
    lnid60sgj = gsn_add_polyline(wks, plot_sgj, (/0, 120, 240, 360/), \
                                 (/60, 60, 60, 60/), lnres)
    plot_sgjs = gsn_csm_contour(wks, sgj_l_mask, res2)
    plot_sgjs  = gsn_contour_shade(plot_sgjs, 0.00, 0.50, opt)
    overlay(plot_sgj, plot_sgjs)

    res@gsnLeftString = "GOSAT"
    res@gsnRightString = "2013 - " + end_year
    plot_gosat = gsn_csm_contour_map(wks, gosat_l_mask, res)
    lnid30gosat = gsn_add_polyline(wks, plot_gosat, (/0, 120, 240, 360/), \
                                   (/30, 30, 30, 30/), lnres)
    lnid60gosat = gsn_add_polyline(wks, plot_gosat, (/0, 120, 240, 360/), \
                                   (/60, 60, 60, 60/), lnres)
    plot_gosats = gsn_csm_contour(wks, gosat_l_mask, res2)
    plot_gosats   = gsn_contour_shade(plot_gosats, 0.00, 0.50, opt)
    overlay(plot_gosat, plot_gosats)

    resP                  = True
    resP@gsnFrame         = False
    resP@gsnPanelLabelBar = True
    resP@lbTitleOn = True
    resP@lbTitlePosition = "Bottom"
    resP@lbTitleString = "XCO~B~2~N~ Data Coverage"
    resP@lbTitleFontHeightF = 0.015
    resP@lbLabelFontHeightF = 0.012
    resP@lbBottomMarginF = -0.2
    resP@lbTopMarginF = 0.3
    resP@lbTitleOffsetF = 0.2
    resP@gsnPanelBottom   = 0.05
    resP@gsnPanelYWhiteSpacePercent = 2.
    gsn_panel(wks, (/plot_scia, plot_sgj, plot_gosat/), (/3, 1/), resP)
    frame(wks)
    delete(outfile)
    delete(wks)
    delete(res)

    log_provenance(ncdf_outfile, \
                   outfile_c3s, \
                   var0 + "fractional data coverage for different obs " \
                   + "periods. Similar to Gier et al 2020, Fig 8.", \
                   (/"mean"/), \
                   DOMAIN, \
                   "histogram", \
                   AUTHORS, \
                   REFERENCES, \
                   ALL_FILES)

  end if  ; end special c3s mask plot

  leave_msg(DIAG_SCRIPT, "")
end
