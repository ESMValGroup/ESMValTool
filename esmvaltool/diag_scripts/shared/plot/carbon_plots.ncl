;#######################################################################
; carbon_plots.ncl
;#######################################################################
; Plots:
; * single_barchart         : one value per model
; * carbon_Line_plot        : time series as line plot
; * carbon_PDF_plot         : line plot for PDF
; * carbon_MarkLine_plot    : time series as mark-line plot
; * carbon_scatter_plot     : scatter plot of X and Y variable
; * carbon_XYscatter_plot   : scatter plot with X and Y error
; * carbon_scatter_num_plot : scatter plot with years as markers
; * carbon_2Line_plot       : line plot comparing two vars
; * carbon_dec_mean_season_plot : line plot of decadel mean seasonal cycle
; * carbon_scatter_circle_plot : plots ellips instead of x and y error bars
;
; Written by Sabrina Wenzel (sabrina.wenzel@dlr.de)
;#######################################################################
load "$diag_scripts/../interface_scripts/auxiliary.ncl"
load "$diag_scripts/../interface_scripts/logging.ncl"

;#######################################################################
undef("barchart")
function barchart(Ydata[*] : numeric,\
                MainStg[1] : string,\
                   XStg[*] : string,\
                   YStg[1] : string,\
                    wks[1], opt[1] : string)
; Arguments:
;  * Ydata   : array for plotting on Y-axis
;  * MainStg : hedar string for plot
;  * YStg    : Y-axis string
;  * wks     : workstation for were to draw the plot
;  * opt     : type of operations:
;            "draw"  : draws the plot directly and returnes logical
;            "panel" : returnes plotvalue for pannel plot

local plot, idx_obs, idx_mod, sres, dim_MOD, colors, newYdata, newXStg, newcolors
begin

  ; Check for valid data
  funcname = "barchart"
  scriptname = "diag_scripts/shared/plot/carbon_plots.ncl"
  enter_msg(scriptname, funcname)

  ; Define names for x-ticks and dimsizes 
  dim_MOD  = dimsizes(Ydata(:))
  colors   = project_style(Ydata@info, diag_script_info, "colors")
  newYdata = Ydata
  newXStg  = XStg
  newcolors  = colors
  newPattern = new(dim_MOD,integer)
  newPattern(0) = 0

  do xi = 1,dimsizes(XStg)-1
    newPattern(xi) = where(XStg(xi-1).eq.XStg(xi), 16, 0)
  end do

  ; Resources
  sres = True
  sres@gsnDraw       = False;True
  sres@gsnFrame      = False
  sres@gsnMaximize   = True
  sres@gsnXYBarChart = True
  sres@gsnXYBarChartBarWidth = 0.30  ; change bar widths

  sres@trXMinF   = 1-0.5
  sres@trXMaxF   = dim_MOD+0.5

  sres@vpWidthF  = 0.5
  sres@vpHeightF = 0.4

  sres@tmXBMode             = "Explicit"  ; explicit labels
  sres@tmXBValues           = ispan(1,dim_MOD,1)
  sres@tmXBLabels           = newXStg
  sres@tmXBLabelFontHeightF = 0.025
  sres@tmXBLabelAngleF      = 30
  sres@tmXBLabelJust        = "TopRight"
  sres@tmXTLabelFontHeightF = 0.020
  sres@tmYLLabelFontHeightF = 0.025

  sres@tiMainFontHeightF  = 0.025
  sres@tiXAxisFontHeightF = 0.025
  sres@tiMainString       = MainStg
  sres@tiYAxisString      = YStg
  ; if (.not.any(models@name.eq."OBS")) then 
  ;   sres@trYMinF = 1.
  ; end if
  sres@gsnYRefLine   = 0.
  sres@gsnXYBarChartColors2 = newcolors 
  sres@gsnXYBarChartPatterns = newPattern

  plot = gsn_csm_xy(wks,fspan(1,dim_MOD,dim_MOD),newYdata,sres)	

  if (opt.eq."draw") then
    draw(plot)
    frame (wks)
    delete(plot)
    return(0)
  end if

  if (opt.eq."panel") then
    return(plot)
  end if

end
;#######################################################################
undef("error_bar_plot")
function error_bar_plot(Ydata[*] : numeric,\
                         Ystd[*] : numeric,\
                      MainStg[1] : string,\
                         YStg[1] : string,\
                          wks[1], opt[1] : string)
;; Arguments:
;;  * Ydata   : array for plotting on Y-axis
;;  * Yerror  : array with + and - 95% confidence intervals of Ydata
;;  * XBname  : array with X-axis label names
;;  * MainStg : hedar string for plot
;;  * YStg    : Y-axis string
;;  * wks     : workstation for were to draw the plot
;;  * opt     : type of operations:
;;      "draw"  : draws the plot directly and returnes logical
;;      "panel" : returnes plotvalue for pannel plot
;;
;; Modification history
;;    201402_A_We_Sa: written 
;;    201411_A_We_Sa: added
local plot,ti,error_bar,error_barend,dim_VAR,lgLabels,dashes
begin

    ;; Check for valid data

  ;; Define names for x-ticks and dimsizes
  dim_MOD  = dimsizes(Ydata(:))
  colors   = project_style(info, diag_script_info, "colors")
  lgLabels = models@name  
  error_bar    = new(dim_MOD,graphic)
  error_barend = new(dim_MOD,graphic)
  marker = new(dim_MOD,graphic)

  res = True
  res@gsnDraw            = False   ; Do not draw plot 
  res@gsnFrame           = False   ; Do not advance frame 
  res@gsnMaximize        = True 
  
  res@tiXAxisFontHeightF = 0.025 
  res@tiYAxisFontHeightF = 0.025 
  res@tiMainFontHeightF  = 0.025 
  res@tiYAxisFontHeightF = 0.028
  res@tiMainString       = MainStg
  res@tiYAxisString      = YStg

  res@vpHeightF	= 0.5             ; change aspect ratio of plot 
  res@vpWidthF 	= 0.5 
 
  res@xyMarker       = 16
  res@xyMarkerSizeF  = 0.01
  res@xyMarkLineMode = "Markers"
  res@xyMonoMarkerColor = False
  res@xyMonoMarker   = False
  res@xyMarkerColors = "white"

  ti = fspan(1,dim_MOD,dim_MOD)
  res@trXMinF   = 1-0.5
  res@trXMaxF   = dim_MOD+0.5
  res@tmXBMode        = "Explicit"
  res@tmXBValues      = ispan(1,dim_MOD,1)
  res@tmXBLabels      = lgLabels
  res@tmXBLabelAngleF = 60
  res@tmXBLabelJust   = "TopRight"

  plot = gsn_csm_xy(wks,fspan(1,dim_MOD,dim_MOD),Ydata,res)	

  polyres = True     
  polyres@xyMonoMarkerColor = False
  polyres@xyMonoMarker   = False
  polyres@gsMarkerIndex  = 16
  polyres@xyMarkerSizeF  = 0.01
  do t=0,dim_MOD-1
    if (t.eq.dim_MOD-1) then
      polyres@gsMarkerColor = "red"
      polyres@gsLineColor = "red"
    else
      polyres@gsMarkerColor = "black"
      polyres@gsLineColor = "black"
    end if
    marker(t) = gsn_add_polymarker(wks,plot,ti(t),Ydata(t),polyres)
    error_bar(t) = gsn_add_polyline(wks,plot,(/ti(t),ti(t)/),(/Ydata(t)+Ystd(t),Ydata(t)-Ystd(t)/),polyres)
    error_barend(t) = gsn_add_polyline(wks,plot,(/ti(t)-0.25,ti(t)+0.25,ti(t),ti(t),ti(t)-0.25,ti(t)+0.25/),\
                                                (/Ydata(t)+Ystd(t),Ydata(t)+Ystd(t),Ydata(t)+Ystd(t),\
                                                  Ydata(t)-Ystd(t),Ydata(t)-Ystd(t),Ydata(t)-Ystd(t)/),polyres)
  end do

  if (opt.eq."draw") then
    draw(plot)
    frame (wks)
    delete(plot)
    return(0)
  end if

  if (opt.eq."panel") then
      str1 = unique_string("marker")
      str2 = unique_string("error_bar")
      str3 = unique_string("error_barend")
      plot@$str1$ = marker
      plot@$str2$ = error_bar 
      plot@$str3$ = error_barend
    return(plot)
  end if
end
;#######################################################################
undef("carbon_Line_plot")
function carbon_Line_plot(Xdata : numeric,\
                          Ydata : numeric,\
                        MainStg[1] : string,\
                           YStg[1] : string,\
                           xMin[1] : string,\
                           xMax[1] : string,\
                            wks[1], opt[1] : string)
;; Arguments:
;;  * Xdata   : array (containing time) for plotting on X-axis
;;  * Ydata   : array for plotting
;;  * MainStg : plot hedar string 
;;  * YStg    : Y-axis string
;;  * xMin    : time for X-axis lower limit
;;  * xMax    : time for X-axis upper limit
;;  * wks     : workstation for were to draw the plot
;;  * opt     : type of operations:
;;            "draw"  : draws the plot directly and returnes logical
;;            "panel" : returnes plotvalue for pannel plot

local plot_array,colors,dashes,thicks,plot_num,res 
begin

    ;; Check for valid data

    ;; Define line colors (first color is used for the observations)
    ;if dimsizes()
    if (dimsizes(models@name).eq.2) then
      colors = (/"black","red"/);
    else
      colors = project_style(info, diag_script_info, "colors");(/"black","red"/);
    end if
    dashes = 0;project_style(info, diag_script_info, "dashes")
    ;thicks = project_style(info, diag_script_info, "avgstd")

    ;; Resources
    res          = True
    res@gsnDraw  = False           ; Do not draw plot
    res@gsnFrame = False           ; Do not advance frame

    res@txFontHeightF      = 0.025
    res@tiXAxisString      = "Year"
    res@tiXAxisFontHeightF = 0.025
    res@tiYAxisString      = YStg
    res@tiYAxisFontHeightF = 0.025
    res@tiMainFontHeightF  = 0.025
    res@tiMainString       = MainStg
    ;res@tiMainJust         = "topLeft"
    ;res@tiMainOffsetXF     = -0.10
 
    res@vpHeightF       = 0.5               ; change aspect ratio of plot
    res@vpWidthF        = 0.8

    res@xyDashPattern   = 0;dashes
    res@xyMonoLineColor = False
    res@xyMonoLineThickness = False
    res@xyLineColors    = colors
    res@xyLineThicknesses = where(models@name .eq."OBS",2,1)

    res@trXMinF = xMin;toint(min(Xdata));diag_script_info@cl_xmin_year
    res@trXMaxF = xMax;toint(max(Xdata));diag_script_info@cl_xmax_year
    res@trYMinF = min(Ydata)
    res@trYMaxF = max(Ydata)

    plot_array = gsn_csm_xy(wks,Xdata,Ydata,res)

    if (opt.eq."draw") then
      if (dimsizes(dimsizes(Ydata)).eq.2) then
      ;;add legend:
      lgres                 = True
      lgres@lgPerimOn       = False
      lgres@lgMonoDashIndex = True
      lgres@lgLineDashSegLenF  = 0.4
      lgres@vpWidthF           = 0.25
      lgres@vpHeightF          = 0.02+0.03*dimsizes(Ydata(:,0))
      lgres@lgLabelFontHeightF = 0.1
      lgres@lgLineColors       = colors
      lgres@lgDashIndex        = 0;dashes
      lgres@lgLineThicknessF   = 2
      lbid = gsn_create_legend(wks,dimsizes(Ydata(:,0)),models@name,lgres)

      amres                    = True
      amres@amParallelPosF     = -0.30
      amres@amOrthogonalPosF   = -0.20
      ;annoid = gsn_add_annotation(plot_array,lbid,amres)
      delete(lbid)
      ;delete(annoid)
    end if

      draw(plot_array)
      frame(wks)
      ;delete(plot_array)
      return(0)
    end if

    if (opt.eq."panel") then
      return(plot_array)
    end if

end

;#######################################################################
undef("carbon_PDF_plot")
function carbon_PDF_plot(Xdata[*] : numeric,\
                          Ydata[*] : numeric,\
                        MainStg[1] : string,\
                           YStg[1] : string,\
                           XStg[1] : string,\
                           xMin[1] : string,\
                           xMax[1] : string,\
                          color[*] : string,\
                            wks[1], opt[1] : string)
;; Arguments:
;;  * Xdata   : array (containing time) for plotting on X-axis
;;  * Ydata   : array for plotting
;;  * MainStg : plot hedar string 
;;  * YStg    : Y-axis string
;;  * XStg    : X-axis string
;;  * xMin    : time for X-axis lower limit
;;  * xMax    : time for X-axis upper limit
;;  * wks     : workstation for were to draw the plot
;;  * opt     : type of operations:
;;            "draw"  : draws the plot directly and returnes logical
;;            "panel" : returnes plotvalue for pannel plot

local plot_array,colors,dashes,thicks,plot_num 
begin

    ;; Check for valid data

    ;; Define line colors (first color is used for the observations)
    colors = project_style(info, diag_script_info, "colors")
    dashes = 0;project_style(info, diag_script_info, "dashes")
    thicks = project_style(info, diag_script_info, "thicks")

    ;; Resources
    res          = True
    res@gsnDraw  = False           ; Do not draw plot
    res@gsnFrame = False           ; Do not advance frame

    res@tiXAxisString      = XStg
    res@tiXAxisFontHeightF = 0.025
    res@tiYAxisString      = YStg
    res@tiYAxisFontHeightF = 0.025

    res@tiMainFontHeightF  = 0.025
    res@tiMainString       = MainStg
    ;res@tiMainJust         = "topLeft"
    ;res@tiMainOffsetXF     = -0.10
 
    res@vpHeightF       = 0.6      ; change aspect ratio of plot
    res@vpWidthF        = 0.6

    res@txFontHeightF     = 0.025
    res@xyDashPattern     = 0;2
    res@xyLineColor       = color;"red";"black"
    res@xyLineThicknessF  = 2

    res@trXMinF = 0.6;xMin;diag_script_info@cl_xmin_year
    res@trXMaxF = 2.0;xMax;diag_script_info@cl_xmax_year
    res@trYMinF = 0;min(Ydata)
    res@trYMaxF = 10.0;max(Ydata)+max(Ydata)/20

    plot_array = gsn_csm_xy(wks,Xdata,Ydata,res)

    if (opt.eq."draw") then
      draw(plot_array)
      frame (wks)
      delete(plot_array)
      return(0)
    end if

    if (opt.eq."panel") then
      return(plot_array)
    end if

end
;#######################################################################
undef("carbon_MarkLine_plot")
function carbon_MarkLine_plot(Xdata[*][*] : numeric,\
                              Ydata[*][*] : numeric,\
                            MainStg[1] : string,\
                               YStg[1] : string,\
                               xMin[1] : string,\
                               xMax[1] : string,\
                                wks[1], opt[1] : string)
;; Arguments:
;;  * Xdata   : array (containing time) for plotting on X-axis
;;  * Ydata   : array for plotting on Y-axis
;;  * MainStg : plot hedar string 
;;  * YStg    : Y-axis string
;;  * xMin    : time for X-axis lower limit
;;  * xMax    : time for X-axis upper limit
;;  * wks     : workstation for were to draw the plot
;;  * opt     : type of operations:
;;            "draw"   : draws the plot directly and returnes logical
;;            "panel" : returnes plotvalue for pannel plot

local plot_array,colors,dashes,thicks,lgLabels,plot_num
begin

    ;; Check for valid data

    ;; Define line colors (first color is used for the observations)
    colors = project_style(info, diag_script_info, "colors")
    dashes = 0;project_style(info, diag_script_info, "dashes")
    thicks = project_style(info, diag_script_info, "thicks")
    lgLabels = models@name  
    plot_num = integertochar(ispan(97,122,1))

    ;; Resources
    res                    = True
    res@gsnDraw            = False           ; Do not draw plot
    res@gsnFrame           = False           ; Do not advance frame

    res@tiXAxisString      = "Year"
    res@tiXAxisFontHeightF = 0.025
    res@tiYAxisString      = YStg
    res@tiYAxisFontHeightF = 0.025

    res@tiMainFontHeightF  = 0.025
    res@tiMainString       = MainStg
    ;res@tiMainJust         = "topLeft"
    ;res@tiMainOffsetXF     = -0.10
 
    res@vpHeightF       = 0.5                    ; change aspect ratio of plot
    res@vpWidthF        = 0.8

    res@txFontHeightF   = 0.025

    res@xyMonoLineColor   = False
    res@xyLineColors      = colors
    res@xyDashPattern     = dashes
    res@xyLineThicknesses = thicks          ; LineThick
    res@xyMarkLineMode    = "MarkLines"
    res@xyMarker          = 1 
    res@xyMarkerSizeF     = 0.015
    res@xyMonoMarkerColor = False
    res@xyMarkerColors    = colors

    res@trXMinF = toint(xMin);-0.5*toint(xMin);diag_script_info@cl_xmin_year
    res@trXMaxF = toint(xMax);+0.5*toint(xMax);diag_script_info@cl_xmax_year
    res@trYMinF = min(Ydata)
    res@trYMaxF = max(Ydata)

    plot_array = gsn_csm_xy(wks, Xdata, Ydata, res)

    if (opt.eq."draw") then
      draw(plot_array)
      frame (wks)
      delete(plot_array)
      return(0)
    end if

    if (opt.eq."panel") then
      return(plot_array)
    end if

end

;#######################################################################
undef("carbon_scatter_plot")
function carbon_scatter_plot(Xdata[*] : numeric,\
                           ; Xerror[*] : numeric,\
                             Ydata[*] : numeric,\
                           MainStg[1] : string,\
                              YStg[1] : string,\
                              XStg[1] : string,\
                               wks[1], opt[1] : string)

;; Arguments:
;;  * Xdata   : array for plotting on X-axis
;;  * Xerror  : array containing errors of Xdata for error bars
;;  * Ydata   : array for plotting on Y-axis
;;  * MainStg : plot hedar string 
;;  * YStg    : Y-axis string
;;  * XStg    : X-axis string
;;  * wks     : workstation for were to draw the plot
;;  * opt     : type of operations:
;;            "draw"   : draws the plot directly and returnes logical
;;            "pannel" : returnes plotvalue for pannel plot

local plot_array, marker, error_bar, error_bar_end, bend,\
      str1, str2, str3, xmin, xmax, res, pmres, dim_MOD, plot_num
begin

    ;; Check for valid data

    ;; Plot array
    dim_MOD       = dimsizes(models@name)
    plot_array    = new(1, graphic)
    marker        = new((/1,dim_MOD/),graphic) 
    error_bar     = new((/1,dim_MOD/),graphic)
    error_bar_end = new((/1,dim_MOD/),graphic)

    ;; Define line colors (first color is used for the observations)
    colors = project_style(info, diag_script_info, "colors")
    thicks = project_style(info, diag_script_info, "avgstd");("thicks")
    lgLabels = models@name  

    plot_num = integertochar(ispan(97,122,1))

    ;; Resources
    res          = True 
    res@gsnDraw  = False               ; Do not draw plot 
    res@gsnFrame = False               ; Do not advance frame 
    ;res@gsnMaximize        = True 
  
    res@tiXAxisString      = XStg
    res@tiXAxisFontHeightF = 0.025 
    res@tiYAxisString      = YStg
    res@tiYAxisFontHeightF = 0.025 
  
    res@tiMainFontHeightF = 0.025 
    res@tiMainString      = MainStg
    ;res@tiMainJust        = "TopLeft"
    ;res@tiMainOffsetXF    = -0.10
    res@tiMainFontColor   = "black"

    res@vpHeightF	= 0.5          ; change aspect ratio of plot 
    res@vpWidthF 	= 0.5
  
    res@xyMarkers         = 4
    res@xyMarkerSizeF     = 0.015
    res@xyMarkerColors    = "black";"white"
    res@xyMarkLineModes   = "Markers"
    res@xyMonoMarkerColor = False 

    xmin = ind(min(Xdata))
    xmax = ind(max(Xdata))
    ;res@trXMinF  = min(Xdata)-(2*Xerror(xmin))
    ;res@trXMaxF  = max(Xdata)+(2*Xerror(xmax))

    plot_array = gsn_csm_xy(wks,Xdata,Ydata,res)  

    pmres = True
    pmres@xyMonoMarkerColor  = False
    pmres@gsMarkerSizeF      = 0.01
    pmres@gsMarkerThicknessF = 2
    bend = max(Ydata(:))/20

    ;do imod = 0,dim_MOD-1
    ;  pmres@gsMarkerIndex = thicks(imod)
    ;  pmres@gsMarkerColor = colors(imod) 
    ;  if (.not.ismissing(Ydata(imod)).and..not.ismissing(Xdata(imod))) then
    ;    error_bar(0,imod) = gsn_add_polyline(wks,plot_array(0),\
    ;                                        (/Xdata(imod)+Xerror(imod),Xdata(imod)-Xerror(imod)/),\
    ;                                        (/Ydata(imod),Ydata(imod)/),pmres)      
    ;    error_bar_end(0,imod) = gsn_add_polyline(wks,plot_array(0),\
    ;                                            (/Xdata(imod)+Xerror(imod),Xdata(imod)+Xerror(imod),\
    ;                                              Xdata(imod)+Xerror(imod),Xdata(imod)-Xerror(imod),\
    ;                                              Xdata(imod)-Xerror(imod),Xdata(imod)-Xerror(imod)/),\
    ;                                            (/Ydata(imod)-bend,Ydata(imod)+bend,Ydata(imod),\
    ;                                              Ydata(imod),Ydata(imod)-bend,Ydata(imod)+bend/),pmres)

        ;marker(0,imod) = gsn_add_polymarker(wks,plot_array(0),Xdata(imod),Ydata(imod),pmres)
     ; end if
    ;end do

    if (opt.eq."draw") then
      lgres                    = True
      lgres@lgPerimOn          = False
      lgres@lgItemType         = "Markers"
      lgres@lgMonoMarkerColor  = False
      lgres@lgMonoMarkerIndex  = False
      lgres@vpWidthF           = 0.20
      lgres@vpHeightF          = 0.21
      lgres@lgLabelFontHeightF = 0.1
      lgres@lgMarkerColors     = colors(1:)
      lgres@lgMarkerIndexes    = thicks(1:)
      lgres@lgMarkerSizeF      = 0.010
      lbid = gsn_create_legend(wks,dim_MOD-1,lgLabels(1:),lgres)

      amres                    = True
      amres@amParallelPosF     = -0.30
      amres@amOrthogonalPosF   = -0.30
      ;annoid = gsn_add_annotation(plot_array,lbid,amres)

      draw(plot_array)
      frame (wks)
      delete(plot_array)
      return(0)
    end if

    if (opt.eq."panel") then
      str1 = unique_string("marker")
      str2 = unique_string("error_bar")
      str3 = unique_string("error_bar_end")
      plot_array@$str1$ = marker
      plot_array@$str2$ = error_bar 
      plot_array@$str3$ = error_bar_end
      return(plot_array)
    end if
    
end
;#######################################################################
undef("carbon_XYscatter_plot")
function carbon_XYscatter_plot(Xdata[*] : numeric,\
                            Xerror[*] : numeric,\
                             Ydata[*] : numeric,\
                            Yerror[*] : numeric,\
                           MainStg[1] : string,\
                              YStg[1] : string,\
                              XStg[1] : string,\
                               wks[1], opt[1] : string)

;; Arguments:
;;  * Xdata   : array for plotting on X-axis
;;  * Xerror  : array containing errors of Xdata for error bars
;;  * Ydata   : array for plotting on Y-axis
;;  * Yerror  : array containing errors of Ydata for error bars
;;  * MainStg : plot hedar string 
;;  * YStg    : Y-axis string
;;  * XStg    : X-axis string
;;  * wks     : workstation for were to draw the plot
;;  * opt     : type of operations:
;;            "draw"   : draws the plot directly and returnes logical
;;            "pannel" : returnes plotvalue for pannel plot

local plot_array, marker, errorX_bar, error_Xbar_end, Xbend, error_Ybar,\
      error_Ybar_end, Ybend, str1, str2, str3, str4, str5, xmin, xmax, \
      res, pmres, dim_MOD, plot_num
begin

    ;; Check for valid data

    ;; Plot array
    dim_MOD       = dimsizes(Xdata(:));dimsizes(models@name)

    plot_array    = new(1, graphic)
    marker        = new((/1,dim_MOD/),graphic) 
    error_Xbar     = new((/1,dim_MOD/),graphic)
    error_Xbar_end = new((/1,dim_MOD/),graphic)
    error_Ybar     = new((/1,dim_MOD/),graphic)
    error_Ybar_end = new((/1,dim_MOD/),graphic)

    ;; Define line colors (first color is used for the observations)
    colors = project_style(info, diag_script_info, "colors")
    thicks = project_style(info, diag_script_info, "avgstd");("thicks")
    lgLabels = models@name  

    ;; Resources
    res          = True 
    res@gsnDraw  = False               ; Do not draw plot 
    res@gsnFrame = False               ; Do not advance frame 
  
    res@tiXAxisString      = XStg
    res@tiXAxisFontHeightF = 0.025 
    res@tiYAxisString      = YStg
    res@tiYAxisFontHeightF = 0.025 

    res@tiMainFontHeightF = 0.025 
    res@tiMainString      = MainStg
    ;res@tiMainJust        = "TopLeft"
    ;res@tiMainOffsetXF    = -0.10
    ;res@tiMainFontColor   = "black"
  
    res@vpHeightF	= 0.5                    ; change aspect ratio of plot 
    res@vpWidthF 	= 0.5 
  
    res@xyMarkers         = 1;thicks;(:)
    res@xyMarkerSizeF     = 0.015
    res@xyMarkerColors    = "white";(:) 
    res@xyMarkLineModes   = "Markers"
    res@xyMonoMarkerColor = False 

    ;xmin = ind(min(Xdata))
    ;xmax = ind(max(Xdata))
    ;res@trXMinF  = min(Xdata)-(2*Xerror(xmin))
    ;res@trXMaxF  = max(Xdata)+(2*Xerror(xmax))

    plot_array = gsn_csm_xy(wks,Xdata,Ydata,res)  

    pmres = True
    pmres@xyMonoMarkerColor  = False
    pmres@gsMarkerSizeF      = 0.01
    pmres@gsMarkerThicknessF = 2
    Xbend = max(Ydata(:))/50
    Ybend = max(Xdata(:))/50

    do imod = 0,dim_MOD-1
      pmres@gsMarkerIndex = thicks(imod)
      pmres@gsMarkerColor = colors(imod) 

      if (.not.ismissing(Ydata(imod)).and..not.ismissing(Xdata(imod))) then
        ;if (models@name(imod).ne."OBS") then
          error_Xbar(0,imod) = gsn_add_polyline(wks,plot_array(0),\
                                             (/Xdata(imod)+Xerror(imod),\
                                               Xdata(imod)-Xerror(imod)/),\
                                             (/Ydata(imod),Ydata(imod)/),pmres)      

          error_Xbar_end(0,imod) = gsn_add_polyline(wks,plot_array(0),\
                                (/Xdata(imod)+Xerror(imod),Xdata(imod)+Xerror(imod),\
                                  Xdata(imod)+Xerror(imod),Xdata(imod)-Xerror(imod),\
                                  Xdata(imod)-Xerror(imod),Xdata(imod)-Xerror(imod)/),\
                                (/Ydata(imod)-Xbend,Ydata(imod)+Xbend,Ydata(imod),\
                                  Ydata(imod),Ydata(imod)-Xbend,Ydata(imod)+Xbend/),pmres)

          error_Ybar(0,imod) = gsn_add_polyline(wks,plot_array(0),\
                                             (/Xdata(imod),Xdata(imod)/),\
                                             (/Ydata(imod)+Yerror(imod),\
                                               Ydata(imod)-Yerror(imod)/),pmres)      

          error_Ybar_end(0,imod) = gsn_add_polyline(wks,plot_array(0),\
                                (/Xdata(imod)-Ybend,Xdata(imod)+Ybend,Xdata(imod),\
                                  Xdata(imod),Xdata(imod)-Ybend,Xdata(imod)+Ybend/),\
                                (/Ydata(imod)+Yerror(imod),Ydata(imod)+Yerror(imod),\
                                  Ydata(imod)+Yerror(imod),Ydata(imod)-Yerror(imod),\
                                  Ydata(imod)-Yerror(imod),Ydata(imod)-Yerror(imod)/),pmres)

          marker(0,imod) = gsn_add_polymarker(wks,plot_array(0),Xdata(imod),Ydata(imod),pmres)
        end if
      ;end if
    end do

    if (opt.eq."draw") then
      lgres                    = True
      lgres@lgPerimOn          = False
      lgres@lgItemType         = "Markers"
      lgres@lgMonoMarkerColor  = False
      lgres@lgMonoMarkerIndex  = False
      lgres@vpWidthF           = 0.20
      lgres@vpHeightF          = 0.21
      lgres@lgLabelFontHeightF = 0.1
      lgres@lgMarkerColors     = colors
      lgres@lgMarkerIndexes    = thicks
      lgres@lgMarkerSizeF      = 0.010
      lbid = gsn_create_legend(wks,dim_MOD,lgLabels,lgres)

      amres                    = True
      amres@amParallelPosF     = 0.30
      amres@amOrthogonalPosF   = -0.30
      ;annoid = gsn_add_annotation(plot_array,lbid,amres)

      draw(plot_array)
      frame (wks)
      delete(plot_array)
      return(0)
    end if

    if (opt.eq."panel") then
      str1 = unique_string("marker")
      str2 = unique_string("error_Xbar")
      str3 = unique_string("error_Xbar_end")
      str4 = unique_string("error_Ybar")
      str5 = unique_string("error_Ybar_end")
      plot_array@$str1$ = marker
      plot_array@$str2$ = error_Xbar 
      plot_array@$str3$ = error_Xbar_end
      plot_array@$str4$ = error_Ybar 
      plot_array@$str5$ = error_Ybar_end
      return(plot_array)
    end if
    
end
;#######################################################################
undef("carbon_scatter_circle_plot")
function carbon_scatter_circle_plot(Xdata[*] : numeric,\
                            Xerror[*] : numeric,\
                             Ydata[*] : numeric,\
                            Yerror[*] : numeric,\
                           MainStg[1] : string,\
                              YStg[1] : string,\
                              XStg[1] : string,\
                               wks[1], opt[1] : string)

;; Arguments:
;;  * Xdata   : array for plotting on X-axis
;;  * Xerror  : array containing errors of Xdata for error bars
;;  * Ydata   : array for plotting on Y-axis
;;  * Yerror  : array containing errors of Ydata for error bars
;;  * MainStg : plot hedar string 
;;  * YStg    : Y-axis string
;;  * XStg    : X-axis string
;;  * wks     : workstation for were to draw the plot
;;  * opt     : type of operations:
;;            "draw"   : draws the plot directly and returnes logical
;;            "pannel" : returnes plotvalue for pannel plot

local plot_array, marker, errorX_bar, error_Xbar_end, Xbend, error_Ybar,\
      error_Ybar_end, Ybend, str1, str2, str3, str4, str5, xmin, xmax, \
      res, pmres, dim_MOD, plot_num
begin

    ;; Check for valid data

    ;; Plot array
    dim_MOD       = dimsizes(Xdata(:));dimsizes(models@name)

    plot_array   = new(1, graphic)
    marker       = new((/1,dim_MOD/),graphic) 
    error_fillcircle = new((/1,dim_MOD/),graphic)
    error_linecircle = new((/1,dim_MOD/),graphic)

    ;; Define line colors (first color is used for the observations)
    colors = project_style(info, diag_script_info, "colors")
    thicks = project_style(info, diag_script_info, "avgstd");("thicks")
    lgLabels = models@name  

    ;; Resources
    res          = True 
    res@gsnDraw  = False               ; Do not draw plot 
    res@gsnFrame = False               ; Do not advance frame 
  
    res@tiXAxisString      = XStg
    res@tiXAxisFontHeightF = 0.022 
    res@tiYAxisString      = YStg
    res@tiYAxisFontHeightF = 0.022 

    res@tiMainFontHeightF = 0.025 
    res@tiMainString      = MainStg
    ;res@tiMainJust        = "TopLeft"
    ;res@tiMainOffsetXF    = -0.10
    ;res@tiMainFontColor   = "black"
  
    res@vpHeightF	= 0.5                    ; change aspect ratio of plot 
    res@vpWidthF 	= 0.5 
  
    res@xyMarkers         = 1;thicks;(:)
    res@xyMarkerSizeF     = 0.015
    res@xyMarkerColors    = "white";(:) 
    res@xyMarkLineModes   = "Markers"
    res@xyMonoMarkerColor = False 

    ;xmin = ind(min(Xdata))
    ;xmax = ind(max(Xdata))
    ;res@trXMinF  = min(Xdata)-(2*Xerror(xmin))
    ;res@trXMaxF  = max(Xdata)+(2*Xerror(xmax))

    plot_array = gsn_csm_xy(wks,Xdata,Ydata,res)
  
    pmres = True
    pmres@xyMonoMarkerColor  = False
    pmres@gsMarkerSizeF      = 0.015
    pmres@gsMarkerThicknessF = 2

    getvalues plot_array
    "trXMinF":xmin
    "trXMaxF":xmax
    "trYMinF":ymin
    "trYMaxF":ymax
    "vpWidthF":vpw
    "vpHeightF":vph
    end getvalues

    th = ispan(0,360,1)
    pi = acos(-1.0)

    do imod = 0,dim_MOD-1
      pmres@gsMarkerColor = colors(imod)
      pmres@gsMarkerIndex = thicks(imod)
      gsres = True
      gsres@gsLineColor   = colors(imod) 
      gsres@gsLineThicknessF = 2.0

      if (.not.ismissing(Ydata(imod)).and..not.ismissing(Xdata(imod))) then
        rx = Xerror(imod)/(xmax-xmin)
        ry = Yerror(imod)/(ymax-ymin)

        xpts1 = cos(th*pi/180)
        ypts1 = sin(th*pi/180)

        xpts = rx*xpts1*(xmax-xmin) + Xdata(imod)
        ypts = ry*ypts1*(ymax-ymin) + Ydata(imod)

        if (models@name(imod).eq."OBS") then
            gsres@gsFillIndex   = 6    
            gsres@gsFillScaleF = .7
            gsres@gsLineColor   = colors(imod) 
            gsres@gsFillColor   = "grey80";colors(imod)	
            error_fillcircle(0,imod) = gsn_add_polygon(wks,plot_array(0),xpts,ypts,gsres)	
        end if

        error_linecircle(0,imod) = gsn_add_polyline(wks,plot_array(0),xpts,ypts,gsres)
        delete([/rx,ry,xpts1,ypts1,xpts,ypts,gsres/])

        marker(0,imod) = gsn_add_polymarker(wks,plot_array(0),Xdata(imod),Ydata(imod),pmres)
      end if

    end do

    if (opt.eq."draw") then
      lgres                    = True
      lgres@lgPerimOn          = False
      lgres@lgItemType         = "Markers"
      lgres@lgMonoMarkerColor  = False
      lgres@lgMonoMarkerIndex  = False
      lgres@vpWidthF           = 0.20
      lgres@vpHeightF          = 0.21
      lgres@lgLabelFontHeightF = 0.05
      lgres@lgMarkerColors     = colors
      lgres@lgMarkerIndexes    = thicks
      lgres@lgMarkerSizeF      = 0.010
      lbid = gsn_create_legend(wks,dim_MOD,lgLabels,lgres)

      amres                    = True
      amres@amParallelPosF     = 0.30
      amres@amOrthogonalPosF   = -0.30
      ;annoid = gsn_add_annotation(plot_array,lbid,amres)

      draw(plot_array)
      frame (wks)
      delete(plot_array)
      return(0)
    end if

    if (opt.eq."panel") then
      str1 = unique_string("marker")
      str2 = unique_string("error_Xbar")
      str3 = unique_string("error_Xbar_end")
      str4 = unique_string("error_Ybar")
      str5 = unique_string("error_Ybar_end")
      plot_array@$str1$ = marker
      plot_array@$str2$ = error_Xbar 
      plot_array@$str3$ = error_Xbar_end
      plot_array@$str4$ = error_Ybar 
      plot_array@$str5$ = error_Ybar_end
      return(plot_array)
    end if
    
end

