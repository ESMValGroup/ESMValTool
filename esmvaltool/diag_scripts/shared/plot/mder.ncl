; #############################################################################
; shared/plot/mder.ncl
; #############################################################################
; GENERAL PLOTTING ROUTINES FOR MULTIPLE DIAGNOSTIC ENSEMBLE REGRESSION
; #############################################################################
; Check the header of each routine for documentation.
;
; Contents:
;     function error_bar_plot: plots an error bar for every single value.
;     function time_raster_plot: plots corellation coeficients in a raster.
;     function weight_lin_plot: plots a time series together with weighted
;                               model ensemble and observations.
;     function mder_scatter_plot: creates a scatterplot with regression error
;                                 and xy-stddev of each model value.
;     function squared_error_plot: plots the squared error for various model
;                                  ensembles.
;     function diag_scatter: plots scatterplot for each diagnostic of
;                            historical simulation vs projection.
;
; #############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/scatterplot.ncl"

load "$diag_scripts/shared/mder.ncl"


; #############################################################################
undef("error_bar_plot")
function error_bar_plot(Ydata[*] : numeric, \
                        Ystd[*][*] : numeric, \
                        XBname[*] : string, \
                        color[*] : string, \
                        MainStg[1] : string, \
                        YStg[1] : string, \
                        wks[1], \
                        opt[1] : string)
;
; Arguments:
;     Ydata: array for plotting on Y-axis.
;     Ystd: array with + and - 95% confidence intervals of Ydata.
;     XBname: array with X-axis label names.
;     color: array indicating positive or negative values.
;     MainStg: header string for plot.
;     YStg: Y-axis string.
;     wks: workstation for were to draw the plot.
;     opt: type of operations:
;          "draw": draws the plot directly and returnes logical.
;          "panel": returnes plotvalue for pannel plot.
;
; Modification history:
;     20180726-schlund_manuel: ported to v2.0.
;     201402??-wenzel_sabrina: written.
;     201411??-wenzel_sabrina: added collor code for pos and ne values.
;
local plot, ti, error_bar, error_barend, DIM_VAR
begin

  funcname = "error_bar_plot"
  scriptname = "diag_scripts/shared/plot/mder.ncl"
  enter_msg(scriptname, funcname)

  ; Check for valid data
  if (all(ismissing(Ydata(:)))) then
    log_info("No data available, skipping")
    return(0)
  end if

  ; Define names for x-ticks and dimsizes
  DIM_VAR  = dimsizes(Ydata(:))
  error_bar = new(DIM_VAR, graphic)
  error_barend = new(DIM_VAR, graphic)
  marker = new(DIM_VAR, graphic)

  res = True
  res@gsnDraw = False                   ; Do not draw plot
  res@gsnFrame = False                  ; Do not advance frame
  res@gsnMaximize = True

  res@tiXAxisFontHeightF = 0.020
  res@tiYAxisFontHeightF = 0.020
  res@tiMainFontHeightF  = 0.025
  res@tiMainString = MainStg
  res@tiXAxisString = YStg

  res@vpHeightF = 0.60                  ; change aspect ratio of plot
  res@vpWidthF = 0.30

  res@xyMarker = 16
  res@xyMarkerSizeF = 0.0012
  res@xyMarkLineMode = "Markers"
  res@xyMonoMarkerColor = False
  res@xyMonoMarker = True
  res@xyMarkerColors = color(0)         ; black

  ti = ispan(0, DIM_VAR - 1, 1)
  res@trYMinF = -0.5
  res@trYMaxF = DIM_VAR-1+0.5
  res@trYReverse = True                 ; reverse y-axis
  res@tmYLMode = "Explicit"
  res@tmYLValues = ispan(0, DIM_VAR - 1, 1)
  res@tmYLLabels = XBname
  res@trXMinF = 0.0
  res@trXMaxF = 1.0

  ; Create plot
  plot = gsn_csm_xy(wks, \
                    (/Ydata(0), Ydata(0)/), \
                    (/ti(0), ti(0)/), \
                    res)

  polyres = True
  polyres@gsMarkerIndex = 16
  polyres@gsMarkerSizeF = 0.012
  pollyres = True
  pollyres@gsLineThicknessF = 2

  do t = 0, DIM_VAR - 1
    polyres@gsMarkerColor = color(t)
    marker(t) = gsn_add_polymarker(wks, plot, Ydata(t), t, polyres)
    error_bar(t) = gsn_add_polyline(wks, plot, (/t, t/), \
                                    (/Ystd(0, t), Ystd(1, t)/), polyres)
    error_barend(t) = gsn_add_polyline(wks, plot, \
                                       (/Ystd(0, t), Ystd(0, t), Ystd(0, t), \
                                        Ystd(1, t), Ystd(1, t), Ystd(1, t)/), \
                                       (/t - 0.25, t + 0.25, t, t, t - 0.25, \
                                        t + 0.25/), pollyres)
  end do

  leave_msg(scriptname, funcname)

  if (opt .eq. "draw") then
    draw(plot)
    frame(wks)
    return(0)
  end if

  if (opt .eq. "panel") then
    return(plot)
  end if

end


; #############################################################################
undef("time_raster_plot")
function time_raster_plot(Ydata[*][*] : numeric, \
                          XBname[*] : string, \
                          YLname[*] : string, \
                          MainStg[1] : string, \
                          wks[1], \
                          opt[1] : string)
;
; Arguments:
;     Ydata: array for plotting.
;     XBname: array with X-axis label names.
;     YBname: array with Y-axis label names.
;     MainStg: header string for plot.
;     wks: workstation for were to draw the plot.
;     opt: type of operations:
;          "draw": draws the plot directly and returnes logical.
;          "panel": returnes plotvalue for pannel plot.
;
; Modification history
;     20180726-schlund_manuel: ported to v2.0.
;     201504??-wenzel_sabrina: written.
;
local plot, ti, error_bar, error_barend, DIM_VAR, pxi
begin

  funcname = "time_raster_plot"
  scriptname = "diag_scripts/shared/plot/mder.ncl"
  enter_msg(scriptname, funcname)

  ; Check required attributes of diag_script_info
  req_atts = (/"p_step"/)
  exit_if_missing_atts(diag_script_info, req_atts)
  p_step = diag_script_info@p_step

  ; Check for valid data
  if (all(ismissing(Ydata(:, :)))) then
    log_info("No data available, skipping")
    return(0)
  end if

  ; Define names for x-ticks and dimsizes
  dim_X = dimsizes(Ydata(0, :))
  dim_Y = dimsizes(Ydata(:, 0))

  ; Color map
  gsn_define_colormap(wks, "temp_19lev")

  res = True
  res@gsnDraw = False                   ; Do not draw plot
  res@gsnFrame = False                  ; Do not advance frame
  res@gsnMaximize = True
  res@gsnSpreadColors = True

  res@cnFillOn = True                   ; Turn on contour fill
  res@cnFillMode = "RasterFill"         ; Turn on raster fill
  res@cnLinesOn = False                 ; Turn off contour lines
  res@tiMainFontHeightF = 0.025
  res@tiMainString = MainStg

  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF = -0.8
  res@cnMaxLevelValF = 0.8
  res@cnLevelSpacingF = 0.01
  res@lbTitleFontHeightF = 0.020

  res@lbLabelBarOn = True
  res@lbOrientation = "vertical"        ; vertical label bar
  res@lbLabelStride = 10
  res@lbBoxLinesOn = False

  res@trXMinF = -0.5
  res@trXMaxF = dim_X - 0.5
  res@tmXBMode = "Explicit"
  res@tmXBValues = ispan(0, dim_X - 1, 1)
  res@tmXBLabels = XBname
  res@tmXBLabelAngleF = 35
  res@tmXBLabelJust   = "TopRight"
  res@tmXBLabelFontHeightF   = 0.020
  res@tmXBMinorValues = ispan(0, dim_X - 1, 1) - 0.5
  res@tmXMinorGrid = True               ; implement y grid
  res@tmXMinorGridThicknessF = 2.0      ; 2.0 is default
  res@sfXArray = res@tmXBValues

  res@trYMinF = -0.5
  res@trYMaxF = dim_Y - 0.5
  res@trYReverse = True                 ; reverse y-axis
  res@tmYLMode = "Explicit"
  res@tmYLValues = ispan(0, dim_Y - 1, 1)
  res@tmYLLabels = YLname
  res@tmYLLabelFontHeightF   = 0.022
  res@tmYLMinorValues = ispan(0, dim_Y - 1, 1) - 0.5
  res@tmYMinorGrid = True               ; implement y grid
  res@tmYMinorGridThicknessF = 2.0      ; 2.0 is default
  res@sfYArray = res@tmYLValues

  ; Create plot
  plot = gsn_csm_contour(wks, Ydata, res)

  ; Plot crosses
  climo_dir = getenv("ESMValTool_climo_dir")
  diag_file = config_user_info@work_dir + variables(0) + "_" + \
    input_file_info@exp(dimsizes(input_file_info@exp) - 1) + "_" + \
    (p_step + 1) + "ystep_diagnostics.txt"
  if (isfilepresent(diag_file)) then
    diag_terms = asciiread(diag_file, -1, "string")
    pdx = dimsizes(diag_terms) / dimsizes(res@sfXArray)
    respoly = True
    respoly@gsMarkerColor = "black"
    respoly@gsMarkerIndex = 5
    respoly@gsMarkerThicknessF = 2
    respoly@gsMarkerSizeF = 30
    pxi = new((/dimsizes(res@sfXArray), dimsizes(res@sfYArray)/), graphic)
    do xi = 0, dimsizes(res@sfXArray) - 1
      idx = diag_terms(xi * pdx: xi * pdx + pdx - 1)
      do yi = 0, dimsizes(res@sfYArray) - 1
        do di = 0, dimsizes(idx) - 1
          if (.not. ismissing(idx(di))) then
            if (YLname(yi) .eq. idx(di)) then
              pxi(xi, yi) = gsn_add_polymarker(wks, \
                                               plot, \
                                               res@tmXBValues(xi), \
                                               res@tmYLValues(yi), \
                                               respoly)
            end if
          end if
        end do
      end do
      delete(idx)
    end do
  end if

  leave_msg(scriptname, funcname)

  if (opt .eq. "draw") then
    draw(plot)
    frame(wks)
    return(0)
  end if

  if (opt .eq. "panel") then
    return(plot)
  end if

end


; #############################################################################
undef("weight_lin_plot")
function weight_lin_plot(Ydata[*][*] : numeric, \
                         Ywterror[*][*] : numeric, \
                         Ywerror[*][*] : numeric, \
                         Yobs[*][*] : numeric, \
                         MainStg[1] : string, \
                         YStg[1] : string, \
                         wks[1], \
                         opt[1] : string, \
                         items: list)
;
; Arguments:
;     Ydata: array with Y data for each model.
;     Ywterror: error of MDER results.
;     Ywerror: error of time-dependent MDER results.
;     Yobs: array with Y data for observations.
;     MainStg: plot header string.
;     YStg: Y-axis string.
;     wks: workstation for were to draw the plot.
;     opt: type of operations:
;          "draw": draws the plot directly and returnes logical.
;          "panel": returnes plotvalue for pannel plot.
;    items: list of input_file_info metadata
;
; Modification history
;     20180726-schlund_manuel: ported to v2.0.
;     201410??-wenzel_sabrina: written.
;
local plot_array, plot, ploto, lbid, annoid, wtmmm, ummm, colors, dashes, \
  thicks, plot_num, res, bres, tline, zline, wmmm, atts, obs, dim_obs, \
  lg_labels, lg_colors, item_type
begin

  funcname = "weight_lin_plot"
  scriptname = "diag_scripts/shared/plot/mder.ncl"
  enter_msg(scriptname, funcname)

  ; Check required attributes of diag_script_info
  req_atts = (/"calc_type"/)
  exit_if_missing_atts(diag_script_info, req_atts)
  calc_type = diag_script_info@calc_type

  ; Check for valid data
  if (all(ismissing(Ydata(:, :)))) then
    log_info("No data available, skipping")
    return(0)
  end if

  ; Check for plotting of whole MDER time series
  plot_mder_time_series = False
  if (isatt(diag_script_info, "plot_mder_time_series")) then
    plot_mder_time_series = diag_script_info@plot_mder_time_series
  end if

  ; Define line colors (first color is used for the observations)
  dashes = project_style(items, diag_script_info, "dashes")
  obs = get_obs_list(items)
  dim_obs = ListCount(obs)

  if (any(calc_type .eq. (/"pos", "int"/))) then
    lg_labels = new((/3 + dim_obs/), string)
    lg_colors = new((/3 + dim_obs/), string)
    item_type = new((/3 + dim_obs/), string)

    OBSdashes = (/4, 2, 0, 9/)
    lg_labels(dim_obs + 1) = "CMIP5 models"
    lg_labels(dim_obs + 2) = Ydata&dataset(1)
    lg_colors(dim_obs + 1) = "grey50"         ; "blue2";"slateblue4";colors
    lg_colors(dim_obs + 2) = "blue2"          ; "slateblue4"
    item_type(dim_obs + 1) = "Lines"
    item_type(dim_obs + 2) = "Markers"
  else if (calc_type .eq. "trend") then
    lg_labels = new((/2 + dim_obs/), string)
    lg_colors = new((/2 + dim_obs/), string)
    item_type = new((/2 + dim_obs/), string)

    OBSdashes = (/4, 2, 9/)
    lg_labels(dim_obs + 1) = Ydata&dataset(1)
    lg_colors(dim_obs + 1) = "blue2"          ; "slateblue4"
    item_type(dim_obs + 1) = "Markers"
  else
    error_msg("f", DIAG_SCRIPT, "", "calc_type " + calc_type + \
              " not supported")
  end if
  end if

  ; Define array with legend names
  lg_labels(:dim_obs - 1) = Yobs&name
  lg_labels(dim_obs) = Ydata&dataset(0)

  ; Define array with colors
  lg_colors(:dim_obs - 1) = "orange"
  lg_colors(dim_obs) = "red"

  ; Define item types
  item_type(:dim_obs) = "Lines"

  error_bar = new(dimsizes(Ydata&time), graphic)
  error_barend = new(dimsizes(Ydata&time), graphic)

  ; Resources
  res = True
  res@gsnDraw = False                 ; Do not draw plot
  res@gsnFrame = False                ; Do not advance frame
  res@gsnMaximize = True

  res@txFontHeightF = 0.025
  res@tiXAxisString = "Year"
  res@tiXAxisFontHeightF = 0.025
  res@tiYAxisString = YStg
  res@tiYAxisFontHeightF = 0.025
  res@tiMainFontHeightF = 0.025
  res@tiMainString  = MainStg

  res@vpHeightF = 0.4                 ; change aspect ratio of plot
  res@vpWidthF = 0.6

  res@trXMinF = min(Ydata&time)
  res@trXMaxF = max(Ydata&time) - 1
  res@trYMinF = min(Ydata(4:, :)) - 0.02 * abs(avg(Ydata(4:, :)))
  res@trYMaxF = max(Ydata(4:, :)) + 0.02 * abs(avg(Ydata(4:, :)))

  ; Plot error bar of time-dependent MDER results
  if (plot_mder_time_series .and. (calc_type .eq. "pos")) then
    res@xyLineColor = "skyblue1"
    res@gsnXYAboveFillColors = "skyblue1"
  else
    res@xyLineColor = "white"
    res@gsnXYAboveFillColors = "white"
  end if
  plot_array = gsn_csm_xy(wks, Ydata&time, Ywerror(:, :), res)

  ; Plot time-dependent MDER results
  if (plot_mder_time_series .and. (calc_type .eq. "pos")) then
    mder = True
    mder@gsLineDashPattern = 4
    mder@xyMonoLineColor = False
    mder@gsLineColor = "slateblue4"
    mder@gsLineThicknessF = 3
    wmmm = gsn_add_polyline(wks, plot_array, Ydata&time, Ydata(2, :), mder)
  end if

  ; Plot models
  lres = True
  lres@xyMonoDashPattern = False
  lres@gsLineColor = "grey50"
  lres@gsLineThicknessF = 2
  plotm = new((/dimsizes(Ydata(3:, 0))/), graphic)
  do idat = 3, dimsizes(Ydata(:, 0)) - 1
    lres@gsLineDashPattern = dashes(idat)
    plotm(idat - 3) = gsn_add_polyline(wks, plot_array, Ydata&time, \
                                       Ydata(idat, :), lres)
  end do

  ; Plot observations
  mres = True
  mres@gsLineColor = "orange"
  mres@gsLineThicknessF = 3
  ploto = new((/dimsizes(Yobs(:, 0))/), graphic)
  do iobs = 0, dimsizes(Yobs(:, 0)) - 1
    mres@gsLineDashPattern = OBSdashes(iobs)
    ploto(iobs) = gsn_add_polyline(wks, plot_array, Yobs&time, \
                                   Yobs(iobs, :), mres)
  end do

  ; Plot uMMM
  mres@gsLineColor = "red"
  mres@gsLineThicknessF = 3
  mres@gsLineDashPattern = 2
  ummm = gsn_add_polyline(wks, plot_array, Ydata&time, Ydata(0, :), mres)

  ; Plot error bars of MDER results
  erres = True
  erres@gsLineColor = "blue2"         ; "slateblue4"
  erres@gsLineThicknessF = 3
  do tt = 0, dimsizes(Ydata&time) - 1
    if (.not. ismissing(Ywterror(0, tt))) then
      error_bar(tt) = gsn_add_polyline(wks, plot_array, \
                                       (/Ydata&time(tt), Ydata&time(tt)/), \
                                       (/Ywterror(0, tt), \
                                        Ywterror(1, tt)/), \
                                       erres)
      error_barend(tt) = gsn_add_polyline(wks, plot_array, \
                                          (/Ydata&time(tt) - 2.0, \
                                           Ydata&time(tt) + 2.0, \
                                           Ydata&time(tt), \
                                           Ydata&time(tt), \
                                           Ydata&time(tt) - 2.0, \
                                           Ydata&time(tt) + 2.0/), \
                                          (/Ywterror(0, tt), \
                                           Ywterror(0, tt), \
                                           Ywterror(0, tt), \
                                           Ywterror(1, tt), \
                                           Ywterror(1, tt), \
                                           Ywterror(1, tt)/), \
                                          erres)
    end if
  end do

  ; Plot means of MDER results
  mder = True
  mder@gsMarkerColor = "blue2"        ; "slateblue4"
  mder@gsMarkerIndex = 9
  mder@gsMarkerThicknessF = 3
  mder@gsMarkerSizeF = 0.01
  wtmmm = gsn_add_polymarker(wks, plot_array, Ydata&time, Ydata(1, :), mder)

  ; Plot line at Y = 0
  lres = True
  lres@gsLineColor = "black"
  lres@gsLineThicknessF = 1
  zline = gsn_add_polyline(wks, plot_array, (/1900, 2100/), (/0, 0/), lres)

  leave_msg(scriptname, funcname)

  if (opt.eq."draw") then
    lgres = True
    lgres@lgPerimOn = False
    ; lgres@lgMonoitem_type = False
    ; lgres@lgitem_types = item_type
    lgres@vpWidthF = 0.30
    lgres@vpHeightF = 0.025 * dimsizes(OBSdashes)
    lgres@lgLabelFontHeightF = 0.20
    lgres@lgLineColors = lg_colors
    lgres@lgDashIndexes = OBSdashes
    lgres@lgLineThicknessF = 2
    OBSdashes(:dimsizes(OBSdashes) - 2) = -1
    lgres@lgMarkerIndexes = OBSdashes
    lgres@lgMarkerSizeF = 0.009
    lgres@lgMarkerThicknessF = 3
    lbid = gsn_create_legend(wks, dimsizes(OBSdashes), lg_labels, lgres)

    amres = True
    amres@amParallelPosF = 0.25
    amres@amOrthogonalPosF = -0.37
    annoid = gsn_add_annotation(plot_array, lbid, amres)

    str1 = unique_string("annoid")
    plot_array@$str1$ = annoid

    draw(plot_array)
    frame(wks)
    return(0)
  end if

  if (opt .eq. "panel") then
    return(plot_array)
  end if

end


; #############################################################################
undef("mder_scatter_plot")
function mder_scatter_plot(Xdata[*] : numeric, \
                           Xerr[*] : numeric, \
                           Ydata[*] : numeric, \
                           Yerr[*] : numeric, \
                           Yfit[*] : numeric, \
                           Yfit_err[*][*] : numeric, \
                           MainStg[1] : string, \
                           XStg[1] : string, \
                           YStg[1] : string, \
                           wks[1], \
                           opt[1] : string, \
                           items: list)
;
; Arguments:
;     Xdata: array with X data for each model.
;     Xerr: array with error values for X.
;     Ydata: array with Y data for each model.
;     Yfit[*]: array with regression values.
;     Yfit_err: array with error regression values.
;     Xobs: OBS value.
;     Xobs_std: OBS error value.
;     MainStg: plot header string.
;     XStg: X-axis string.
;     YStg: Y-axis string.
;     wks: workstation for were to draw the plot.
;     opt: type of operations:
;          "draw": draws the plot directly and returnes logical.
;          "panel": returnes plotvalue for pannel plot.
;    items: list of input_file_info metadata
;
local plot_array, plot, colors, dashes, thicks, res, bres
begin

  funcname = "mder_scatter_plot"
  scriptname = "diag_scripts/shared/plot/mder.ncl"
  enter_msg(scriptname, funcname)

  ; Check for valid data
  if (all(ismissing(Xdata))) then
      log_info("No data available, skipping")
      return(0)
  end if

  ; Define line colors (first color is used for the observations)
  colors = project_style(items, diag_script_info, "colors")
  dashes = project_style(items, diag_script_info, "dashes")
  thicks = project_style(items, diag_script_info, "thicks")
  lg_labels = unique_labels_min(items, (/"dataset", "exp", "ensemble"/))

  plot = new((/dimsizes(Xdata)/), graphic)
  plot_m = new((/dimsizes(Xdata)/), graphic)
  Xerror = new((/dimsizes(Xdata)/), graphic)
  Xerror_end = new((/dimsizes(Xdata)/), graphic)
  Yerror = new((/dimsizes(Xdata)/), graphic)
  Yerror_end = new((/dimsizes(Xdata)/), graphic)

  ; Basic resources
  res = True
  res@gsnDraw = False                 ; Do not draw plot
  res@gsnFrame = False                ; Do not advance frame
  res@gsnMaximize = True

  res@txFontHeightF = 0.025
  res@tiXAxisString = XStg
  res@tiXAxisFontHeightF = 0.022
  res@tiYAxisString = YStg
  res@tiYAxisFontHeightF = 0.025
  res@tiMainFontHeightF = 0.025
  res@tiMainString = MainStg

  res@vpHeightF = 0.6                 ; change aspect ratio of plot
  res@vpWidthF = 0.6

  res@trYMinF = min(Ydata) - 0.2 * (abs(min(Ydata)) / 8)
  res@trYMaxF = max(Ydata) + 0.2 * (abs(max(Ydata)) / 8)
  res@trXMinF = max((/res@trYMinF, min(Yfit)/))
  res@trXMaxF = min((/res@trYMaxF, max(Yfit)/))

  ; Plot regression errors (filled area between curves)
  res@xyLineColor = "grey85"
  res@gsnXYFillColors = "grey85"
  plot_array = gsn_csm_xy(wks, Yfit, Yfit_err, res)

  ; Plot observations
  fres = True
  fres@gsLineColor = "orange"
  fres@gsFillColor = "orange"
  fres@gsLineDashPattern = 0
  plot_obs  = gsn_add_polygon(wks, plot_array, \
                              (/Xdata(1) + Xerr(1), Xdata(1) - Xerr(1), \
                               Xdata(1) - Xerr(1), Xdata(1) + Xerr(1)/), \
                              (/100.0, 100.0, -100.0, -100.0/), \
                              fres)

  ; Plot regression line
  fres@gsLineColor = "slateblue4"
  fres@gsLineThicknessF  = 3
  fres@gsLineDashPattern = 0
  plot_yfit = gsn_add_polyline(wks, plot_array, Yfit, Yfit, fres)

  ; Plot observation fit (MDER result)
  fres@gsLineDashPattern = 1
  plot_mder = gsn_add_polyline(wks, plot_array, \
                               (/-100.0, Xdata(1), Xdata(1)/), \
                               (/Ydata(1), Ydata(1), -100.0/), \
                               fres)
  ; Plot uMMM line
  fres@gsLineColor = "red"
  fres@gsLineDashPattern = 3
  plot_ummm = gsn_add_polyline(wks, plot_array, (/100.0, -100.0/), \
                               (/Ydata(0), Ydata(0)/), fres)

  ; Plot models and errors
  bres = True
  bres@gsLineColor = "black"
  bres@gsLineThicknessF = 1.5
  bres@gsMarkerThicknessF = 1.5
  bres@gsMarkerSizeF = 0.01
  bres@gsMarkerColor = "black"
  bres@gsMarkerIndex = 1
  Xbend = 0.01 * (res@trYMaxF - res@trYMinF)
  Ybend = 0.01 * (res@trXMaxF - res@trXMinF)
  txres = True
  txres@txFontColor = "black"
  txres@txFontHeightF = 0.016
  txres@txJust = "TopLeft"
  do idat = 1, dimsizes(Xdata)-1
    ; bres@gsMarkerColor = colors(idat)
    ; bres@gsMarkerIndex = thicks(idat)
    if (.not. ismissing(Yerr(idat)) .and. .not. ismissing(Xerr(idat))) then

      ; Model error bars
      Xerror(idat) = gsn_add_polyline(wks, plot_array, \
                                      (/Xdata(idat) + Xerr(idat), \
                                       Xdata(idat) - Xerr(idat)/), \
                                      (/Ydata(idat), Ydata(idat)/), bres)
      Xerror_end(idat) = gsn_add_polyline(wks, plot_array, \
                                          (/Xdata(idat) + Xerr(idat), \
                                           Xdata(idat) + Xerr(idat), \
                                           Xdata(idat) + Xerr(idat), \
                                           Xdata(idat) - Xerr(idat), \
                                           Xdata(idat) - Xerr(idat), \
                                           Xdata(idat) - Xerr(idat)/), \
                                          (/Ydata(idat) - Xbend, \
                                           Ydata(idat) + Xbend, \
                                           Ydata(idat), Ydata(idat), \
                                           Ydata(idat) - Xbend, \
                                           Ydata(idat) + Xbend/), bres)
      Yerror(idat) = gsn_add_polyline(wks, plot_array, \
                                      (/Xdata(idat), Xdata(idat)/), \
                                      (/Ydata(idat) + Yerr(idat), \
                                       Ydata(idat) - Yerr(idat)/), bres)
      Yerror_end(idat) = gsn_add_polyline(wks, plot_array, \
                                          (/Xdata(idat) - Ybend, \
                                           Xdata(idat) + Ybend, \
                                           Xdata(idat), Xdata(idat), \
                                           Xdata(idat) - Ybend, \
                                           Xdata(idat) + Ybend/), \
                                          (/Ydata(idat) + Yerr(idat), \
                                           Ydata(idat) + Yerr(idat), \
                                           Ydata(idat) + Yerr(idat), \
                                           Ydata(idat) - Yerr(idat), \
                                           Ydata(idat) - Yerr(idat), \
                                           Ydata(idat) - Yerr(idat)/), bres)
    end if

    ; Model means
    plot_m(idat) = gsn_add_polymarker(wks, plot_array, Xdata(idat), \
                                      Ydata(idat), bres)

    ; Model numbers
    if (isatt(Xdata, "model_numbers")) then
      plot(idat) = gsn_add_text(wks, plot_array, \
                                sprinti("%-0.2hi", \
                                        Xdata@model_numbers(idat)), \
                                Xdata(idat), Ydata(idat), txres)
    else
      plot(idat) = gsn_add_text(wks, plot_array, sprinti("%-0.2hi", idat), \
                                Xdata(idat), Ydata(idat), txres)
    end if
  end do

  leave_msg(scriptname, funcname)

  if (opt .eq. "draw") then
    lgres = True
    lgres@lgPerimOn = False
    ; lgres@lgitem_type = "Markers"
    lgres@lgMonoMarkerColor = False
    lgres@lgMonoMarkerIndex = False
    lgres@vpWidthF = 0.18
    lgres@lgLabelFontHeightF = 0.15
    lgres@lgMarkerSizeF = 0.01

    amres = True
    amres@amOrthogonalPosF = 0.
    amres@amParallelPosF = 0.67
    lgres@vpHeightF = 0.6

    lgres@lgMarkerColors = colors(2:)
    lgres@lgMarkerIndexes = thicks(2:)
    lbid = gsn_create_legend(wks, dimsizes(Xdata) - 2, lg_labels(2:), lgres)
    ; annoid1 = gsn_add_annotation(plot_array, lbid, amres)

    draw(plot_array)
    frame(wks)
    return(0)
  end if

  if (opt .eq. "panel") then
    return(plot_array)
  end if

end


; #############################################################################
undef("squared_error_plot")
function squared_error_plot(Ydata[*][*] : numeric, \
                            Ymean[*] : numeric, \
                            Ysqrd[*][*] : numeric, \
                            MainStg[1] : string, \
                            XStg[*] : string, \
                            YStg[1] : string, \
                            wks[1], \
                            opt[1] : string)
;
; Arguments:
;     Ydata: array, where the rightmost dimension contains the box plot
;            reference pt. data:
;            y(n, 0) = bottom_value,
;            y(n, 1) = bottom_value_of_box,
;            y(n, 2) = mid-value_of_box,
;            y(n, 3) = top_value_of_box,
;            y(n, 4) = top_value.
;     Ymean: array containing the multi-model-average of the squared errors.
;     Ysqrd: array containing the models' squared errors.
;     MainStg: plot header strings.
;     XStg: X-axis string.
;     YStg: Y-axis string.
;     wks: workstation for were to draw the plot.
;     opt: type of operations:
;          "draw": draws the plot directly and returnes logical.
;          "panel": returnes plotvalue for pannel plot.
;
local plot, res, llres, mres, sres
begin

  funcname = "squared_error_plot"
  scriptname = "diag_scripts/shared/plot/mder.ncl"
  enter_msg(scriptname, funcname)

  ; Check for valid data
  if (all(ismissing(Ydata))) then
    log_info("No data available, skipping")
    return(0)
  end if

  ; Define names for x-ticks and dimsizes
  DIM_VAR = dimsizes(Ymean)
  Xdata = ispan(0, DIM_VAR - 1, 1)
  dumAVG = new(DIM_VAR, graphic)
  dumSQRT = new((/DIM_VAR, dimsizes(Ysqrd(0, :))/), graphic)

  res = True
  res@gsnMaximize = True

  res@tiYAxisFontHeightF = 0.020
  res@tiMainFontHeightF = 0.025
  res@tiMainString = MainStg
  res@tiYAxisString = YStg

  res@vpHeightF = 0.4                   ; change aspect ratio of plot
  res@vpWidthF = 0.6

  res@trXMinF = -0.5
  res@trXMaxF = DIM_VAR - 1.0 + 0.5
  res@tmXBLabelFontHeightF = 0.018
  res@tmXBMajorOutwardLengthF = 0.015
  res@tmXBMode = "Explicit"
  res@tmXBLabels = XStg
  if (dimsizes(XStg) .gt. 4) then
    res@tmXBLabelJust = "BottomRight"
    res@tmXBLabelAngleF = 28
  end if
  res@tmXTBorderOn = True
  res@tmYRBorderOn = True
  res@trYMinF = -0.2
  res@trYMaxF = max(Ysqrd) + 0.05 * (max(Ysqrd))

  ; Plot boxes
  llres = True
  llres@gsLineThicknessF = 3            ; line thickness
  opti = True
  opti@boxWidth = 0.8                   ; Width of box (x units)
  boxColors = new(dimsizes(Xdata), string)
  boxColors(0::2) = (/"red"/)
  boxColors(1::2) = (/"slateblue4"/)
  opti@boxColors = boxColors            ; Color of box(es)
  plot = boxplot(wks, Xdata, Ydata, opti, res, llres)

  ; Plot models (weighted and unweighted)
  mres = True                           ; marker mods desired
  mres@gsMarkerIndex = 5                ; polymarker style
  mres@gsMarkerThicknessF = 2
  mres@gsMarkerColor = "black"          ; polymarker color

  sres = True                           ; marker mods desired
  sres@gsMarkerIndex = 16               ; polymarker style
  sres@gsMarkerSizeF = 0.004
  sres@gsMarkerColor = "gray55"         ; polymarker color

  do ivar = 0, DIM_VAR - 1
    do ip = 0, dimsizes(Ysqrd(ivar, :)) - 1

      ; Plot models
      dumSQRT(ivar, ip) = gsn_add_polymarker(wks, plot, Xdata(ivar), \
                                             Ysqrd(ivar, ip), sres)
    end do

    ; Plot multi-model means
    dumAVG(ivar) = gsn_add_polymarker(wks, plot, Xdata(ivar), Ymean(ivar), \
                                      mres)
  end do

  leave_msg(scriptname, funcname)

  if (opt .eq. "draw") then
    draw(plot)
    frame(wks)
    return(0)
  end if

  if (opt .eq. "panel") then
    return(plot)
  end if

end


; #############################################################################
undef("diag_scatter")
function diag_scatter(res_in[1] : logical, \
                      diag_all : numeric, \
                      future_diag : numeric, \
                      target_var : string, \
                      items: list)
;
; Description:
;     Scatterplot for each diagnostic of historical simulation vs projection.
;
; Return value:
;     plotnames: array of strings containing plotnames.
;
; Arguments:
;     res_in: resources, must be given as attributes with "res_" prefix.
;     diag_all: [dim_diags x dim_datasets], coordinates must contain names of
;               diagnostics and datasets, respectively.
;     future_diag: [dim_datasets], target diagnostic of projections, can have
;                  attribute "future_exp" for plot title.
;     target_var: name of target variable, needs attribute 'units'.
;     items: list of input_file_info metadata
;
; Modification history:
;     20191121-schlund_manuel: added provenance tracking.
;     20180806-schlund_manuel: ported to v2.0.
;     201806??-senftleben_daniel: written.
;
local dim_diag, dim_model, diags, datasets, x_var, wks, plot
begin

  funcname = "diag_scatter"
  scriptname = "diag_scripts/shared/plot/mder.ncl"
  enter_msg(scriptname, funcname)

  all_diagnostics = diag_all
  all_diagnostics!0 = "diags"
  all_diagnostics!1 = "datasets"
  diags = all_diagnostics&diags
  datasets = all_diagnostics&datasets
  plotnames = new(dimsizes(diags), string)

  ; Name of future projection
  if (isatt(future_diag, "future_exp")) then
    future_exp = " (" + future_diag@future_exp + ")"
  else
    future_exp = ""
  end if

  do id = 0, dimsizes(diags) - 1
    wks = get_wks(0, DIAG_SCRIPT, target_var + "_" + diags(id))
    plotnames(id) = wks@fullname
    x_var = (/all_diagnostics(id, :), future_diag/)
    x_var!1 = "datasets"
    x_var&datasets = datasets
    copy_VarAtts(res_in, x_var)  ; pass resources
    x_var@res_tiMainString  = "Correlation between " + diags(id) + \
      " and " + target_var + future_exp
    x_var@res_tiXAxisString = diags(id)
    x_var@res_tiYAxisString = target_var + " [" + target_var@units + "]"
    x_var@res_trXMinF = min(x_var(0, :))
    x_var@res_trXMaxF = max(x_var(0, :))
    x_var@res_trYMinF = min(x_var(1, :))
    x_var@res_trYMaxF = max(x_var(1, :))
    plot = scatterplot(wks, x_var, target_var, False, False, items)
    draw(plot)
    frame(wks)

    ; Write ncdf file
    work_dir = config_user_info@work_dir
    system("mkdir -p " + work_dir)
    new_path = work_dir + target_var + "_" + diags(id) + ".nc"
    x_var!0 = "variable"
    x_var&variable = (/diags(id), target_var/)
    x_var@var = "value"
    x_var@diag_script = DIAG_SCRIPT
    x_var@ncdf = new_path
    outfile = ncdf_write(x_var, new_path)

    ; Provenance tracking
    if (isStrSubset(diags(id), "_t")) then
      statistics := (/"trend"/)
    else if(isStrSubset(diags(id), "_c")) then
      statistics := (/"mean", "clim"/)
    else
      statistics := (/"mean"/)
    end if
    end if
    authors = (/"schlund_manuel"/)
    domain = diag_script_info@domain
    plot_type = "scatter"
    references = (/"wenzel16jclim"/)
    caption = "Scatterplot between " + target_var + " and " + diags(id) + "."
    ancestors = metadata_att_as_array(items, "filename")
    log_provenance(outfile, plotnames(id), caption, statistics, domain, \
                   plot_type, authors, references, ancestors)
  end do

  leave_msg(scriptname, funcname)
  return(plotnames)

end
