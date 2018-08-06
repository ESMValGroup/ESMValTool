; #############################################################################
; GENERAL PLOTTING ROUTINES FOR MULTIPLE DIAGNOSTIC ENSEMBLE REGRESSION
; #############################################################################
; Check the header of each routine for documentation.
;
; Contents:
;     function error_bar_plot : plots an error bar for every single value.
;     function time_raster_plot : plots corellation coeficients in a raster.
;     function weight_lin_plot : plots a time series together with with
;                                weighted model enseble and observations.
;     function mder_scatter_plot : creates a scatterplot with regression error
;                                  and xy-stddev of each model value.
;     function sqrt_error_plot   : plots the sqared error for various model
;                                  ensembles.
;
; #############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

load "interface_scripts/interface.ncl"

load "./diag_scripts/shared/plot/style.ncl"
load "./diag_scripts/shared/plot/scatterplot.ncl"


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
;     Ydata : array for plotting on Y-axis
;     Ystd : array with + and - 95% confidence intervals of Ydata
;     XBname : array with X-axis label names
;     color : array indicating pos or neg values
;     MainStg : hedar string for plot
;     YStg : Y-axis string
;     wks : workstation for were to draw the plot
;     opt : type of operations:
;       "draw" : draws the plot directly and returnes logical
;       "panel" : returnes plotvalue for pannel plot
;
; Modification history:
;     20180726_A-schl_ma: ported to v2.0
;     201402??_A-wenz_sa: written
;     201411??_A-wenz_sa: added collor code for pos and ne values
;
local plot, ti, error_bar, error_barend, DIM_VAR, lg_labels
begin

  ; Check for valid data
  if (all(ismissing(Ydata(:)))) then
    log_info("No data available, skipping")
    return(0)
  end if

  ; Define names for x-ticks and dimsizes
  lg_labels = dataset_info@dataset
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

  if (opt .eq. "draw") then
    draw(plot)
    frame(wks)
    delete(plot)
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
;     Ydata : array for plotting
;     XBname: array with X-axis label names
;     YBname: array with Y-axis label names
;     MainStg : hedar string for plot
;     wks : workstation for were to draw the plot
;     opt : type of operations:
;       "draw": draws the plot directly and returnes logical
;       "panel" : returnes plotvalue for pannel plot
;
; Modification history
;     20180726_A-schl_ma: ported to v2.0
;     201504??_A-wenz_sa: written
;
local plot, ti, error_bar, error_barend, DIM_VAR, lg_labels
begin

  ; Check for valid data
  if (all(ismissing(Ydata(:, :)))) then
    log_info("No data available, skipping")
    return(0)
  end if

  ; Define names for x-ticks and dimsizes
  lg_labels = dataset_info@dataset
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

  ; Add Crosses (FIXME for v2.0 --> diag_file?)
  ; if (isatt(diag_script_info, "pt_time")) then
  ;   p_time = diag_script_info@pt_time
  ;   p_step = (diag_script_info@pstep)
  ; end if
  ; climo_dir = getenv("ESMValTool_climo_dir")
  ; diag_file = climo_dir + "/" + variables(0) + "_" + \
  ;   input_file_info@exp(dimsizes(input_file_info@exp)-1)+"_"+\
  ;                           (p_step+1)+"ystep_diagnostics.txt"
  ;
  ; if (isfilepresent(diag_file)) then
  ;   diag_terms = asciiread(diag_file, -1, "string")
  ;   pdx = dimsizes(diag_terms) / dimsizes(res@sfXArray)
  ;   respoly = True
  ;   respoly@gsMarkerColor = "black"
  ;   respoly@gsMarkerIndex = 5
  ;   respoly@gsMarkerThicknessF = 2
  ;   respoly@gsMarkerSizeF = 30
  ;   pxi = new((/dimsizes(res@sfXArray), dimsizes(res@sfYArray)/), graphic)
  ;   do xi = 0, dimsizes(res@sfXArray) - 1
  ;     idx = diag_terms(xi * pdx: xi * pdx + pdx - 1)
  ;     do yi = 0, dimsizes(res@sfYArray) - 1
  ;       do di = 0, dimsizes(idx) - 1
  ;         if (.not. ismissing(idx(di))) then
  ;           if (YLname(yi) .eq. idx(di)) then
  ;             pxi(xi, yi) = gsn_add_polymarker(wks, \
  ;                                             plot, \
  ;                                             res@tmXBValues(xi), \
  ;                                             res@tmYLValues(yi), \
  ;                                             respoly)
  ;           end if
  ;         end if
  ;       end do
  ;     end do
  ;     delete(idx)
  ;   end do
  ; end if

  if (opt .eq. "draw") then
    draw(plot)
    frame(wks)
    delete(pxi)
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
                         opt[1] : string)
;
; Arguments:
;     Ydata : array with Y data for each model
;     Yweight : array with weighted Y data for each model
;     MainStg : plot hedar string
;     YStg : Y-axis string
;     wks : workstation for were to draw the plot
;     opt : type of operations:
;       "draw" : draws the plot directly and returnes logical
;       "panel" : returnes plotvalue for pannel plot
;
; Modification history
;     20180726_A-schl_ma: ported to v2.0
;     201410??_A-wenz_sa: written
;
local plot_array, plot, colors, dashes, thicks, plot_num, res, bres, tline, \
  zline
begin

  ; Check for valid data
  if (all(ismissing(Ydata(:, :)))) then
    log_info("No data available, skipping")
    return(0)
  end if

  ; Define line colors (first color is used for the observations)
  dashes = project_style(diag_script_info, "dashes")

  i_OBS = ind(dataset_info@project .eq. "OBS")
  DIM_OBS = dimsizes(i_OBS)

  if (diag_script_info@typ .eq. "pos") then
    lg_labels = new((/3 + DIM_OBS/), string)
    lg_colors = new((/3 + DIM_OBS/), string)
    item_type = new((/3 + DIM_OBS/), string)

    OBSdashes = (/4, 2, 0, 9/)
    lg_labels(DIM_OBS + 1) = "CMIP5 models"   ; dataset_info@dataset(2)
    lg_labels(DIM_OBS + 2) = dataset_info@dataset(1)
    lg_colors(DIM_OBS + 1) = "grey50"         ; "blue2";"slateblue4";colors
    lg_colors(DIM_OBS + 2) = "blue2"          ; "slateblue4"
    item_type(DIM_OBS + 1) = "Lines"
    item_type(DIM_OBS + 2) = "Markers"
  else if (diag_script_info@typ .eq. "trend") then
    lg_labels = new((/2 + DIM_OBS/), string)
    lg_colors = new((/2 + DIM_OBS/), string)
    item_type = new((/2 + DIM_OBS/), string)

    OBSdashes = (/4, 2, 9/)
    lg_labels(DIM_OBS + 1) = dataset_info@dataset(1)
    lg_colors(DIM_OBS + 1) = "blue2"          ; "slateblue4"
    item_type(DIM_OBS + 1) = "Markers"
  end if
  end if

  ; Define array with legend names
  lg_labels(:DIM_OBS - 1) = Yobs@name
  lg_labels(DIM_OBS) = dataset_info@dataset(0)

  ; Define array with colors
  lg_colors(:DIM_OBS - 1) = "orange"
  lg_colors(DIM_OBS) = "red"

  ; Define item types
  item_type(:DIM_OBS) = "Lines"

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
  res@trXMaxF = 2060                  ; max(Ydata&time) - 1
  res@trYMinF = min(Ydata) - 0.02 * (abs(avg(Ydata)))
  res@trYMaxF = max(Ydata) + 0.02 * (abs(avg(Ydata)))

  ; Create a plot with filled area between both curves
  ; if (diag_script_info@typ .eq. "pos") then
  ;   res@xyLineColor = "skyblue1"
  ;   res@gsnXYAboveFillColors = "skyblue1"
  ; else
  res@xyLineColor          = "white"
  res@gsnXYAboveFillColors = "white"
  ; end if

  ; Create plot
  plot_array = gsn_csm_xy(wks, Ydata&time, Ywerror(:, :), res)

  ; Add MDER to plot
  if (diag_script_info@typ .eq. "pos") then
    mder = True
    mder@gsLineDashPattern = 4
    mder@xyMonoLineColor = False
    mder@gsLineColor = "slateblue4"
    mder@gsLineThicknessF = 3
    ; wmmm = gsn_add_polyline(wks,plot_array,Ydata&time,Ydata(2,:),mder)
  end if

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

  ; Add OBS to plot
  mres = True
  mres@gsLineColor = "orange"
  mres@gsLineThicknessF = 3
  ploto = new((/dimsizes(Yobs(:, 0))/), graphic)
  do iobs = 0, dimsizes(Yobs(:, 0)) - 1
    mres@gsLineDashPattern = OBSdashes(iobs)
    ploto(iobs) = gsn_add_polyline(wks, plot_array, Yobs&time, \
                                   Yobs(iobs, :), mres)
  end do

  ; Add uMMM to plot
  mres@gsLineColor = "red"
  mres@gsLineThicknessF = 3
  mres@gsLineDashPattern = 2
  ummm = gsn_add_polyline(wks, plot_array, Ydata&time, Ydata(0, :), mres)

  ; Add error bars to dependent MDER-markers
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

  ; Add time dependent MDER to plot
  mder = True
  mder@gsMarkerColor = "blue2"        ; "slateblue4"
  mder@gsMarkerIndex = 9
  mder@gsMarkerThicknessF = 3
  mder@gsMarkerSizeF = 0.01
  wtmmm = gsn_add_polymarker(wks, plot_array, Ydata&time, Ydata(1, :), mder)

  lres = True
  lres@gsLineColor = "black"
  lres@gsLineThicknessF = 1
  zline = gsn_add_polyline(wks, plot_array, (/1900, 2100/), (/0, 0/), lres)

  if (opt.eq."draw") then
    lgres = True
    lgres@lgPerimOn = False
    lgres@lgMonoitem_typee = False
    lgres@lgitem_typees = item_type
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
    delete([/plot_array, ploto/])
    delete([/lbid, annoid/])
    delete([/wtmmm, ummm, zline/])

    ; if (diag_script_info@typ.eq."pos") then
    ;   delete(wmmm)
    ; end if
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
                           opt[1] : string)
;
; Arguments:
;     Xdata : array with Y data for each model
;     Yfit[*] : arrey with regression values
;     Yfit_err : arrey with error regression values
;     Xobs : OBS value
;     Xobs_std : OBS error value
;     MainStg : plot hedar string
;     XStg : X-axis string
;     YStg : Y-axis string
;     wks : workstation for were to draw the plot
;     opt : type of operations:
;         "draw" : draws the plot directly and returnes logical
;         "panel" : returnes plotvalue for pannel plot
;
local plot_array, plot, colors, dashes, thicks, res, bres
begin

  ; Check for valid data
  if (all(ismissing(Xdata))) then
      log_info("No data available, skipping")
      return(0)
  end if

  ; Define line colors (first color is used for the observations)
  colors = project_style(diag_script_info, "colors")
  dashes = project_style(diag_script_info, "dashes")
  thicks = project_style(diag_script_info, "thicks")
  lg_labels = unique_labels_min((/"name", "ensemble"/))

  plot = new((/dimsizes(Xdata)/), graphic)
  plot_m = new((/dimsizes(Xdata)/), graphic)
  Xerror = new((/dimsizes(Xdata)/), graphic)
  Xerror_end = new((/dimsizes(Xdata)/), graphic)
  Yerror = new((/dimsizes(Xdata)/), graphic)
  Yerror_end = new((/dimsizes(Xdata)/), graphic)

  ; Resources
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

  ; Create a plot with filled area between both curves
  res@xyLineColor = "grey85"
  res@gsnXYFillColors = "grey85"

  ; Create plot
  plot_array = gsn_csm_xy(wks, Yfit, Yfit_err, res)

  fres = True
  fres@gsLineColor = "orange"
  fres@gsFillColor = "orange"
  fres@gsLineDashPattern = 0
  plot_obs  = gsn_add_polygon(wks, plot_array, \
                              (/Xdata(1) + Xerr(1), Xdata(1) - Xerr(1), \
                               Xdata(1) - Xerr(1), Xdata(1) + Xerr(1)/), \
                              (/100.0, 100.0, -100.0, -100.0/), \
                              fres)
  fres@gsLineColor = "slateblue4"
  fres@gsLineThicknessF  = 3
  fres@gsLineDashPattern = 0
  plot_yfit = gsn_add_polyline(wks, plot_array, Yfit, Yfit, fres)
  if (isatt(diag_script_info, "loglin")) then
    if (diag_script_info@loglin) then
      plot_yfit := gsn_add_polyline(wks, plot_array, Xdata(2:), Yfit, fres)
    end if
  end if
  fres@gsLineDashPattern = 1
  plot_mder = gsn_add_polyline(wks, plot_array, \
                               (/-100.0, Xdata(1), Xdata(1)/), \
                               (/Ydata(1), Ydata(1), -100.0/), \
                               fres)
  fres@gsLineColor = "red"
  fres@gsLineDashPattern = 3
  plot_ummm = gsn_add_polyline(wks, plot_array, (/100.0, -100.0/), \
                               (/Ydata(0), Ydata(0)/), fres)

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
    plot_m(idat) = gsn_add_polymarker(wks, plot_array, Xdata(idat), \
                                      Ydata(idat), bres)
    if (isatt(Xdata, "model_numbers")) then
      plot(idat) = gsn_add_text(wks, plot_array, \
                                sprinti("%-0.2hi", \
                                        Xdata@model_numbers(idat)), \
                                Xdata(idat), Ydata(idat), txres)
    else
      plot(idat) = gsn_add_text(wks, plot_array, sprinti("%-0.2hi", idat), \
                                Xdata(idat), Ydata(idat), txres)
    end if
    ; gsn_add_polymarker(wks, plot_array, Xdata(idat), Ydata(idat), bres)
  end do

  if (opt .eq. "draw") then
    lgres = True
    lgres@lgPerimOn = False
    lgres@lgitem_typee = "Markers"
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
    delete(plot_array)
    return(0)
  end if

  if (opt .eq. "panel") then
    return(plot_array)
  end if

end


; #############################################################################
undef("sqrt_error_plot")
function sqrt_error_plot(Ydata[*][*] : numeric, \
                         Ymean[*] : numeric, \
                         Ysqrd[*][*] : numeric, \
                         MainStg[1] : string, \
                         XStg[*] : string, \
                         Xtitle[1] : string, \
                         YStg[1] : string, \
                         wks[1], \
                         opt[1] : string)
;
; Arguments:
;   Ydata : array, where the rightmost dimension contains the box plot
;           reference pt. data:
;           y(n,0)=bottom_value,
;           y(n,1)=bottom_value_of_box,
;           y(n,2)=mid-value_of_box,
;           y(n,3)=top_value_of_box,
;           y(n,4)=top_value.
;   Ymean : array containing the average of the pt. data
;   Ysqrd : arrey containing the sqared error
;   MainStg : plot hedar string
;   XStg : X-axis string
;   YStg : Y-axis string
;   wks : workstation for were to draw the plot
;   opt : type of operations:
;     "draw": draws the plot directly and returnes logical
;     "panel" : returnes plotvalue for pannel plot
;
local plot, res, llres, mres, sres
begin

  ; Check for valid data
  if (all(ismissing(Ydata))) then
    log_info("No data available, skipping")
    return(0)
  end if

  ; Define names for x-ticks and dimsizes
  lg_labels = dataset_info@dataset
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

  llres = True
  llres@gsLineThicknessF = 3            ; line thickness

  opti = True
  opti@boxWidth = 0.8                   ; Width of box (x units)
  boxColors = new(dimsizes(Xdata), string)
  boxColors(0::2) = (/"red"/)
  boxColors(1::2) = (/"slateblue4"/)
  opti@boxColors = boxColors            ; Color of box(es)

  ; Create plot
  plot = boxplot(wks, Xdata, Ydata, opti, res, llres)

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
      dumSQRT(ivar, ip) = gsn_add_polymarker(wks, plot, Xdata(ivar), \
                                             Ysqrd(ivar, ip), sres)
    end do
    dumAVG(ivar) = gsn_add_polymarker(wks, plot, Xdata(ivar), Ymean(ivar), \
                                      mres)
  end do

  if (opt .eq. "draw") then
    draw(plot)
    frame(wks)
    delete(plot)
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
                      target_var : string)
;
; Description:
;   Scatterplot for each diagnostic of historical simulation vs projection
;
; Return value:
;   plotnames : array of strings containing plotnames
;
; Arguments:
;   res_in: resources, must be given as attributes with "res_" prefix
;   diag_all: [dim_diags x dim_datasets], coordinates must contain
;   names of diagnostics and datasets, respectively
;   future_diag: [dim_diags], target diagnostic of projections
;   target_var: name of target variable, needs attribute 'units'
;
; Modification history:
;   * 20180806-A_schl_ma: ported to v2.0
;   * 201806??-A_senf_da: written.
;
local dim_diag, dim_model, diags, datasets, wks, plot
begin

  diag_all!0 = "diags"
  diag_all!1 = "datasets"
  diags = diag_all&diags
  datasets = diag_all&datasets
  plotnames = new(dimsizes(diags), string)

  do id = 0, dimsizes(diags) - 1
    wks = get_wks(0, diag_script, target_var + "_" + diags(id))
    plotnames(id) = wks@fullname
    x_var = (/diag_all(id, :), future_diag/)
    x_var!1 = "datasets"
    x_var&datasets = datasets
    copy_VarAtts(res_in, x_var)  ; pass resources
    x_var@res_tiMainString  = "Correlation between " + diags(id) + \
      " and " + target_var
    x_var@res_tiXAxisString = diags(id)
    x_var@res_tiYAxisString = target_var + " [" + target_var@units + "]"
    x_var@res_trXMinF = min(x_var(0, :))
    x_var@res_trXMaxF = max(x_var(0, :))
    x_var@res_trYMinF = min(x_var(1, :))
    x_var@res_trYMaxF = max(x_var(1, :))
    plot = scatterplot(wks, x_var, target_var, True, False)
    delete([/wks, x_var, plot/])
  end do

  return(plotnames)

end
; #############################################################################
