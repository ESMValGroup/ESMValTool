; **************************************************************************
; This library contains functions in support of regional_downscaling
; the functions included are:
; * function:
;   select_region1(region:string)
; * function:
;    box_plot_cr(wks:graphic, x[* ]:numeric, y[* ][* ]:numeric, \
;                    boxOpts:logical, plotres:logical, lineres:logical)
; *  function:
; *    box_plot_lg(wks:graphic, x[* ]:numeric, y[* ][* ]:numeric, \
;                    boxOpts:logical, plotres:logical, lineres:logical)
; *  function:
;       select_region_938(region:string, lsdata:numeric)
; *  function
;       select_region_srex_poly(region:string)
; *  function:
;       select_region_srex(region:string, lsdata:numeric)
; *  function:
;       get_unique_intersection( a, b )
; *  function:
;       get_unique_difference( a, b )
; *  procedure:
;       add_legenda_page(wks, MODEL:string, colors:string, dashs:numeric,\
;                         scatters:string, ticknesses:numeric, place:string)
; # Caveats:                                                    #
; #                                                             #
; # Modification history:
; 20220314-weigel_katja: get_unique_intersection and get_unique_difference
;                       (ESMValTool v1 functions to allow labelling of
;                        multiple projects)
; 2021????-weigel_katja: replaces select_region_939 and select_region_940
;                        with select_region_srex, added select_region_srex_poly
; **************************************************************************
undef("select_region1")
function select_region1(region:string)
;
; Arguments
;    region: a string specifying the region to be selected.
;
; Return value
;    An array with the region boundary as (latmin, latmax, lonmin, lonmax)
;    with the name of the region as a string attribute @name.
;
; Description
;    Translates a region specification into lat/lon boundaries and a region
;    name as an attribute.
;
; Modification history
;
;    20141205 -A_gott_kl: adjusted names to Righi et al. (2015).
;    20140410 -A_fran_fr: extended to midlat, equatorial and polar regions.
;    20140129 -A_fran_fr: written.
;
local funcname, scriptname, verbosity, region
begin

  funcname = "select_region1"
  scriptname = "diag_scripts/lib/ncl/latlon.ncl"
  verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
  ; enter_msg(scriptname, funcname, 10)

  if (region.eq."Global") then
    region_array = (/ -90., 90., 0., 360./)
    region_array@name = "Glob"
    return(region_array)
  end if
  if (region.eq."Tropics") then
    region_array = (/ -20., 20., 0., 360./)
    region_array@name = "Trop"
    return(region_array)
  end if
  if (region.eq."NH extratropics") then
    region_array = (/20., 90., 0., 360./)
    region_array@name = "NHext"
    return(region_array)
  end if
  if (region.eq."SH extratropics") then
    region_array = (/ -90.,  -20., 0., 360./)
    region_array@name = "SHext"
    return(region_array)
  end if
  if (region.eq."NH equatorial") then
    region_array = (/0., 30., 0., 360./)
    region_array@name = "NHtrop"
    return(region_array)
  end if
  if (region.eq."SH equatorial") then
    region_array = (/ -30., 0., 0., 360./)
    region_array@name = "SHtrop"
    return(region_array)
  end if
  if (region.eq."Northern Hemisphere") then
    region_array = (/20., 90., 0., 360./)
    region_array@name = "NH"
    return(region_array)
  end if
  if (region.eq."Southern Hemisphere") then
    region_array = (/ -90.,  -20., 0., 360./)
    region_array@name = "SH"
    return(region_array)
  end if
  if (region.eq."NH midlatitudes") then
    region_array = (/35., 60., 0., 360./)
    region_array@name = "NHmidlat"
    return(region_array)
  end if
  if (region.eq."SH midlatitudes") then
    region_array = (/ -60.,  -35., 0., 360./)
    region_array@name = "SHmidlat"
    return(region_array)
  end if
  if (region.eq."Arctic") then
    region_array = (/60., 90., 0., 360./)
    region_array@name = "NHpolar"
    return(region_array)
  end if
  if (region.eq."Antarctic") then
    region_array = (/ -90.,  -60., 0., 360./)
    region_array@name = "SHpolar"
    return(region_array)
  end if
  if (region.eq."Equatorial") then
    region_array = (/ -10., 10., 0., 360./)
    region_array@name = "EQ"
    return(region_array)
  end if
  ; Region latitudes and longitudes following Figure 9.38 Chapter 9 AR5
  ; (as Comunication from Chou Sin Chan)
  ; 1 - ENA Eastern North America
  ; Lon = 250 to 350
  ; Lat = 30 to 67

  if (region.eq."ENA") then
    region_array = (/30., 67., 250., 350./)
    region_array@name = "ENA"
    return(region_array)
  end if
  ; 2 - WNA Western North America
  ; Lon = 190 to 250
  ; Lat = 30 to 67
  if (region.eq."WNA") then
    region_array = (/30., 67., 190., 250./)
    region_array@name = "WNA"
    return(region_array)
  end if
  ; 3  - CAM Central America
  ; Lon = 240 to 300
  ; Lat = 10 to 30
  if (region.eq."CAM") then
    region_array = (/10., 30., 240., 300./)
    region_array@name = "CAM"
    return(region_array)
  end if
  ; 4 - AMZ Amazon
  ; Lon = 277 to 329
  ; Lat =  -20 to 10
  if (region.eq."TSA") then
    region_array = (/ -20., 10., 277., 329./)
    region_array@name = "TSA"
    return(region_array)
  end if
  ; 5 - SSA Southern South America
  ; Lon = 277 to 329
  ; Lat =  -56 to  -20
  if (region.eq."SSA") then
    region_array = (/ -56.,  -20., 277., 329./)
    region_array@name = "SSA"
    return(region_array)
  end if
  ; 6 - EUM Europe and Mediterranean
  ; Lon = 350 to 40
  ; Lat = 30 to 67
  if (region.eq."EUM") then
    region_array = (/30., 67., 350., 40./)
    region_array@name = "ENA"
    return(region_array)
  end if
  ; 7 - NAF Northern Africa
  ; Lon = 342 to 65
  ; Lat = 16 to 30
  if (region.eq."NAF") then
    region_array = (/16., 30., 342., 65./)
    region_array@name = "NAF"
    return(region_array)
  end if
  ; 8 - SAF South Africa
  ; Lon = 10 to 51
  ; Lat =  -35 to  -10
  if (region.eq."SAF") then
    region_array = (/ -35.,  -10., 10., 51./)
    region_array@name = "SAF"
    return(region_array)
  end if
  ; 9 - CAF central Africa
  ; Lon = 342 to 60
  ; Lat =  -10 to 16
  if (region.eq."CAF") then
    region_array = (/ -10., 16., 342., 60./)
    region_array@name = "CAF"
    return(region_array)
  end if
  ; 10 - NAS Nothern Asia
  ; Lon = 40 to 167
  ; Lat = 48 to 67
  if (region.eq."NAS") then
    region_array = (/48., 67., 40., 167./)
    region_array@name = "NAS"
    return(region_array)
  end if
  ; 11 - EAS Eastern Asia
  ; Lon = 100 to 167
  ; Lat = 20 to 48
  if (region.eq."EAS") then
    region_array = (/20., 48., 100., 167./)
    region_array@name = "EAS"
    return(region_array)
  end if
  ; 12 - CAS Central Asia
  ; Lon = 40 to 100
  ; Lat = 30 to 48
  if (region.eq."CAS") then
    region_array = (/30., 48., 40., 100./)
    region_array@name = "CAS"
    return(region_array)
  end if
  ; 13 - SEA Southeast Asia
  ; Lon = 95 to 151
  ; Lat =  -11 to 20
  if (region.eq."SEA") then
    region_array = (/ - 11., 20., 95., 151./)
    region_array@name = "SEA"
    return(region_array)
  end if
  ; 14 - SAS South Asia
  ; Lon = 65 to 95
  ; Lat = 5 to 30
  if (region.eq."SAS") then
    region_array = (/5., 30., 65., 95./)
    region_array@name = "SAS"
    return(region_array)
  end if
  ; 15 - AUS Oceania
  ; Lon = 112 to 180
  ; Lat =  -48 to  -11
  if (region.eq."AUS") then
    region_array = (/ - 48.,  - 11., 112., 180./)
    region_array@name = "AUS"
    return(region_array)
  end if
  if (region.eq."Arctic") then
    region_array = (/67.5, 90., 0., 360./)
    region_array@name = "Arctic"
    return(region_array)
  end if
end
;  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

undef("box_plot_cr")
function box_plot_cr(wks:graphic, x[*]:numeric, y[*][*]:numeric, \
                     boxOpts:logical, plotres:logical, lineres:logical)

begin
  dimquery = dimsizes(y)
  numbox = dimquery(0)
  boxWidths = new((/numbox/), float)
  if (numbox.ne.dimsizes(x)) then
    print("boxplot: Fatal: X must be one - dimensional and both X and Y" + \
          "have the same rightmost dimension")
    exit
  end if
  if (any(ismissing(x))) then
    print("boxplot: Fatal: X array cannot contain missing data, exiting")
    exit
  end if
; Developing x -axis
  xAxis = new(numbox + 2, typeof(x))
  xAxis(1:numbox) = x
  if (numbox.ne.1) then
    dx = x(1) - x(0)
    xAxis(0) = x(0) - dx
    xAxis(numbox + 1) = x(numbox - 1) + dx
  else
    dx = 1
    xAxis(0) = x - dx
    xAxis(2) = x + dx
  end if
  if (boxOpts) then
    if (isatt(boxOpts, "boxWidth")) then
      if (dimsizes(boxOpts@boxWidth).ne.1.and.dimsizes(boxOpts@boxWidth)\
                                    .ne.numbox) then
        print("boxplot: Number of input box widths must either equal 1 or" +\
              "the number of boxes (" + numbox + "). Using first specified" +\
              "box width only.")
        boxWidths(:) = boxOpts@boxWidth(0)
      else
        boxWidths = boxOpts@boxWidth
      end if
    else
      boxWidths(:) = dx * .4
    end if
  else
    boxWidths(:) = dx * .4
  end if
  labarr = new(numbox + 2, "string")  ; Prepare actual X -axis labels...
  labarr(0) = ""
  labarr(numbox + 1) = ""
  labarr(1:numbox) = xAxis(1:numbox)

  ; Whether to maximize plot in frame.
  maximize = get_res_value(plotres, "gsnMaximize", False)

  if (plotres) then
    ;  print("Plot resources detected, accepting")
    fsatts = getvaratts(plotres)
    do ty = 0, dimsizes(fsatts) - 1
      if (fsatts(ty).eq."tmXBLabels") then  ; XBLabels
        if (dimsizes(plotres@tmXBLabels).ne.numbox) then
          print("boxplot: Fatal:Number of XB Labels does not match number" + \
                "of boxes, exiting")
          exit
        else
          labarr(1:numbox) = plotres@$fsatts(ty)$
        end if
        xblab = plotres@tmXBLabels
        delete(plotres@tmXBLabels)  ; Delete so tmXBLabels is not used
      end if                       ; atts are assigned below...
    end do
    delete(fsatts)
  end if
  ind_blue = plotres@ind_blue
  ind_red = plotres@ind_red
  ind_blue = ind_blue + 1
  ind_red = ind_red + 1
  plot = create "plot" logLinPlotClass wks
    "trYMinF" : min(y) - 2
    "trYMaxF" : max(y) + 2
    "trXMinF" : min(xAxis)
    "trXMaxF" : max(xAxis)
    "pmTickMarkDisplayMode" : "Always"
    "tmXBMode"   : "Explicit"
    "tmXBValues" : xAxis(ind_blue)
    "tmXBLabels" : labarr(ind_blue)
    "tmXBLabelFontHeightF" : 0.014
    "tmXBLabelAngleF"      : 90.
    "tmXBLabelJust"        : "TopRight"
    "tmXBLabelFontColor" : "blue"
    "tmYROn"  : False
    "tmXTOn"  : False
    "tmYRBorderOn" : True
    "tmXTBorderOn" : True
    "pmTitleDisplayMode": "Always"
    "tiMainOn"     : True
    "tiMainString"   : ""
  end create
  if (.not.all(ismissing(plotres@ind_red))) then
    plot1 = create "plot" logLinPlotClass wks
      "trYMinF" : min(y) - 2
      "trYMaxF" : max(y) + 2
      "trXMinF" : min(xAxis)
      "trXMaxF" : max(xAxis)
      "pmTickMarkDisplayMode" : "Always"
      "tmXBMode"   : "Explicit"
      "tmXBValues" : xAxis(ind_red)
      "tmXBLabels" : labarr(ind_red)
      "tmXBLabelFontHeightF" : 0.0105
      "tmXBLabelAngleF"      : 90.
      "tmXBLabelFontColor" : "red"
      "tmXBLabelJust"        : "TopRight"
      "tmYROn"  : False
      "tmXTOn"  : False
      "tmYRBorderOn" : True
      "tmXTBorderOn" : True
      "pmTitleDisplayMode": "Always"  ; allow titles
      "tiMainOn"     : True
      "tiMainString"   : ""
    end create
    overlay(plot, plot1)
  end if

  if (plotres) then
    attsetvalues(plot, plotres)
  end if

  polyres = True                 ; Set up defaults
  polyres@gsLineColor = "black"        ; color of lines
  polyres@gsLineThicknessF = 1.5                          ; thickness of lines
  polyres@gsLineDashPattern = 0
  if (lineres) then
    fsatts = getvaratts(lineres)
    do ty = 0, dimsizes(fsatts) - 1
      polyres@$fsatts(ty)$ = lineres@$fsatts(ty)$
    end do
  end if

  if (boxOpts) then
    if (isatt(boxOpts, "boxColors")) then
      boxcolor = boxOpts@boxColors
      if (dimsizes(boxcolor).eq.1.or.dimsizes(boxcolor).ne.numbox) then
        if (dimsizes(boxcolor).ne.numbox) then
          print("boxplot: warning: Number of input colors must either " + \
                "equal 1 or the number of boxes (" + numbox + ")." + \
                "Using first specified color only.")
        end if
        polyres@gsLineColor = boxcolor(0)
        cflag = 1
      else
        cflag = 2
      end if
    else
      cflag = 1
    end if
  else
    cflag = 1
  end if

  dum = new((/numbox, 9/), graphic)
  do gg = 0, numbox - 1
    ff = xAxis(gg + 1)
    if (cflag.eq.2) then
      polyres@gsLineColor = boxcolor(gg)
    end if
    yy = (/y(gg, 4), y(gg, 4)/)
    xx = (/(ff - (boxWidths(gg) / 8.)), (ff + (boxWidths(gg) / 8.))/)
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
      dum(gg, 0) = gsn_add_polyline(wks, plot, xx, yy, polyres)
    end if

    yy = (/y(gg, 3), y(gg, 4)/)
    xx = (/ff, ff/)
    polyres@gsLineDashPattern = 0
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
      dum(gg, 1) = gsn_add_polyline(wks, plot, xx, yy, polyres)
    end if
    polyres@gsLineDashPattern = 0

    yy = (/y(gg, 3), y(gg, 3)/)
    xx = (/(ff - (boxWidths(gg) / 2.)), (ff + (boxWidths(gg) / 2.))/)
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
      dum(gg, 2) = gsn_add_polyline(wks, plot, xx, yy, polyres)
    end if

    yy = (/y(gg, 1), y(gg, 3)/)
    xx = (/(ff - (boxWidths(gg) / 2.)), (ff - (boxWidths(gg) / 2.))/)
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
      dum(gg, 3) = gsn_add_polyline(wks, plot, xx, yy, polyres)
    end if

    yy = (/y(gg, 2), y(gg, 2)/)
    xx = (/(ff - (boxWidths(gg) / 2.)), (ff + (boxWidths(gg) / 2.))/)
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
      dum(gg, 4) = gsn_add_polyline(wks, plot, xx, yy, polyres)
    end if

    yy = (/y(gg, 1), y(gg, 3)/)
    xx = (/(ff + (boxWidths(gg) / 2.)), (ff + (boxWidths(gg) / 2.))/)
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
      dum(gg, 5) = gsn_add_polyline(wks, plot, xx, yy, polyres)
    end if

    yy = (/y(gg, 1), y(gg, 1)/)
    xx = (/(ff - (boxWidths(gg) / 2.)), (ff + (boxWidths(gg) / 2.))/)
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
    dum(gg, 6) = gsn_add_polyline(wks, plot, xx, yy, polyres)
    end if

    yy = (/y(gg, 0), y(gg, 1)/)
    xx = (/ff, ff/)
    polyres@gsLineDashPattern = 0
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
      dum(gg, 7) = gsn_add_polyline(wks, plot, xx, yy, polyres)
    end if
    polyres@gsLineDashPattern = 0

    yy = (/y(gg, 0), y(gg, 0)/)
    xx = (/(ff - (boxWidths(gg) / 8.)), (ff + (boxWidths(gg) / 8.))/)
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
      dum(gg, 8) = gsn_add_polyline(wks, plot, xx, yy, polyres)
    end if
  end do

  dumname = unique_string("dum")
  plot@$dumname$ = dum
  if (isvar("xblab").eq."True") then
    plotres@tmXBLabels = xblab   ; reassign XBLabels for possible future plots
  end if
  if (maximize) then
    mres = True
    mres@gsnDraw = False
    mres@gsnFrame = False
    maximize_output(wks, mres)
  end if
  return(plot)
end
; -------------------------------------
; Adam Phillips
;
undef("box_plot_lg")
function box_plot_lg(wks:graphic, x[*]:numeric, y[*][*]:numeric, \
                     boxOpts:logical, plotres:logical, lineres:logical)
begin
  dimquery = dimsizes(y)
  numbox = dimquery(0)
  boxWidths = new((/numbox/), float)
  if (numbox.ne.dimsizes(x)) then
    print("boxplot: Fatal: X must be one -dimensional and both X and Y " + \
          "must have the same rightmost dimension")
    exit
  end if
  if (any(ismissing(x))) then
    print("boxplot: Fatal: X array cannot contain missing data, exiting")
    exit
  end if

  ;    Developing x -axis
  xAxis = new(numbox + 2, typeof(x))
  xAxis(1:numbox) = x
  if (numbox.ne.1) then
    dx = x(1) - x(0)
    xAxis(0) = x(0) - dx
    xAxis(numbox + 1) = x(numbox - 1) + dx
  else
    dx = 1
    xAxis(0) = x - dx
    xAxis(2) = x + dx
  end if

  if (boxOpts) then
    if (isatt(boxOpts, "boxWidth")) then
      if (dimsizes(boxOpts@boxWidth).ne.1.and.(dimsizes(boxOpts@boxWidth).ne.\
                                               numbox)) then
        print("boxplot: Number of input box widths must either equal 1 " + \
              "or the number of boxes (" + numbox + "). Using first " + \
              "specified box width only.")
        boxWidths(:) = boxOpts@boxWidth(0)
      else
        boxWidths = boxOpts@boxWidth
      end if
    else
      boxWidths(:) = dx * .4
    end if
  else
    boxWidths(:) = dx * .4
  end if

  labarr = new(numbox + 2, "string")
  labarr(0) = ""
  labarr(numbox + 1) = ""
  labarr(1:numbox) = xAxis(1:numbox)

; Whether to maximize plot in frame.
  maximize = get_res_value(plotres, "gsnMaximize", False)

  if (plotres) then
    print("Plot resources detected, accepting")
    fsatts = getvaratts(plotres)
    do ty = 0, dimsizes(fsatts) - 1
      if (fsatts(ty).eq."tmXTLabels") then
        print(plotres@tmXTLabels)
        if (dimsizes(plotres@tmXTLabels).ne.numbox) then
          delete(labarr)
          labarr = new(dimsizes(plotres@tmXTLabels) + 2, "string")
          labarr(0) = ""
          labarr(dimsizes(plotres@tmXTLabels) + 1) = ""
          labarr(1:dimsizes(plotres@tmXTLabels)) = plotres@$fsatts(ty)$
        else
          labarr(1:numbox) = plotres@$fsatts(ty)$
        end if
        xblab = plotres@tmXTLabels
        delete(plotres@tmXTLabels)
      end if
    end do
    delete(fsatts)
  end if
  plot = create "plot" logLinPlotClass wks
    "trYMinF" : - 25
    "trYMaxF" : 125
    "trXMinF" : 0.5
    "trXMaxF" : dimsizes(labarr) + 0.5
    "pmTickMarkDisplayMode" : "Always"
    "tmXTMode"   : "Explicit"
    "tmXTValues" : ispan(0, dimsizes(labarr) - 2, 1)
    "tmXTLabels" : labarr(1:)
    "tmYLMode"   : "Explicit"
    "tmYLValues" : (/5, 25, 50, 75, 95/)
    "tmYLLabels" : (/"5%", "25%", "50%", "75%", "95%"/)
    "tmXTLabelFontColor" : "black"
    "tmYLOn"  : True
    "tmXBOn"  : False
    "tmXTOn"  : True
    "tmXBLabelsOn": False
    "tmXTLabelsOn": True
    "tmXTMinorOn" :False
    "tmXBMinorOn" :False
    "tmYLMinorOn" :False
    "tmXUseBottom": False
    "tmYRBorderOn" : False
    "tmXTBorderOn" : False
    "tmYLBorderOn" : False
    "tmXBBorderOn" : False
    "tmXTMajorLineColor"  : "Transparent"
    "tmYLMajorLineColor"  : "Transparent"
    "tmYLLabelFontHeightF" : 0.02
    "tmXTLabelFontHeightF" : 0.02
    "pmTitleDisplayMode": "Always"
    "tiMainOn"     : True
    "tiMainString"   : ""
  end create

  if (plotres) then
    attsetvalues(plot, plotres)
  end if
  polyres = True
  polyres@gsLineColor  = "black"
  polyres@gsLineThicknessF = 3.
  polyres@gsLineDashPattern = 0
  if (lineres) then
    fsatts = getvaratts(lineres)
    do ty = 0, dimsizes(fsatts) - 1
      polyres@$fsatts(ty)$ = lineres@$fsatts(ty)$
    end do
  end if

  if (boxOpts) then
    if (isatt(boxOpts, "boxColors")) then
      boxcolor = boxOpts@boxColors
      if (dimsizes(boxcolor).eq.1.or.dimsizes(boxcolor).ne.numbox) then
        if (dimsizes(boxcolor).ne.numbox) then
          print("boxplot: warning: Number of input colors must either" + \
                "equal 1 or the number of boxes (" + numbox + ")." + \
                "Using first specified color only.")
        end if
        polyres@gsLineColor = boxcolor(0)
        cflag = 1
      else
        cflag = 2
      end if
    else
      cflag = 1
    end if
  else
    cflag = 1
  end if
  dum = new((/numbox, 9/), graphic)
  do gg = 0, numbox - 1
    ff = xAxis(gg + 1)
    if (cflag.eq.2) then
      polyres@gsLineColor = boxcolor(gg)
    end if
    yy = (/y(gg, 4), y(gg, 4)/)
    xx = (/(ff - (boxWidths(gg)/8.)), (ff + (boxWidths(gg)/8.))/)
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
      dum(gg, 0) = gsn_add_polyline(wks, plot, xx, yy, polyres)
    end if

    yy = (/y(gg, 3), y(gg, 4)/)
    xx = (/ff, ff/)
    polyres@gsLineDashPattern = 0
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
      dum(gg, 1) = gsn_add_polyline(wks, plot, xx, yy, polyres)
    end if
    polyres@gsLineDashPattern = 0

    yy = (/y(gg, 3), y(gg, 3)/)
    xx = (/(ff - (boxWidths(gg)/2.)), (ff + (boxWidths(gg)/2.))/)
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
      dum(gg, 2) = gsn_add_polyline(wks, plot, xx, yy, polyres)
    end if

    yy = (/y(gg, 1), y(gg, 3)/)
    xx = (/(ff - (boxWidths(gg)/2.)), (ff - (boxWidths(gg)/2.))/)
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
      dum(gg, 3) = gsn_add_polyline(wks, plot, xx, yy, polyres)
    end if

    yy = (/y(gg, 2), y(gg, 2)/)
    xx = (/(ff - (boxWidths(gg)/2.)), (ff + (boxWidths(gg)/2.))/)
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
      dum(gg, 4) = gsn_add_polyline(wks, plot, xx, yy, polyres)
    end if

    yy = (/y(gg, 1), y(gg, 3)/)
    xx = (/(ff + (boxWidths(gg)/2.)), (ff + (boxWidths(gg)/2.))/)
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
      dum(gg, 5) = gsn_add_polyline(wks, plot, xx, yy, polyres)
    end if

    yy = (/y(gg, 1), y(gg, 1)/)
    xx = (/(ff - (boxWidths(gg)/2.)), (ff + (boxWidths(gg)/2.))/)
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
      dum(gg, 6) = gsn_add_polyline(wks, plot, xx, yy, polyres)
    end if

    yy = (/y(gg, 0), y(gg, 1)/)
    xx = (/ff, ff/)
    polyres@gsLineDashPattern = 0
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
      dum(gg, 7) = gsn_add_polyline(wks, plot, xx, yy, polyres)
    end if
    polyres@gsLineDashPattern = 0

    yy = (/y(gg, 0), y(gg, 0)/)
    xx = (/(ff - (boxWidths(gg)/8.)), (ff + (boxWidths(gg)/8.))/)
    if (.not.(any(ismissing(xx).or.ismissing(yy)))) then
      dum(gg, 8) = gsn_add_polyline(wks, plot, xx, yy, polyres)
    end if
  end do

  dumname = unique_string("dum")
  plot@$dumname$ = dum
  if (isvar("xblab").eq."True") then
    ; plotres@tmXBLabels =
    plotres@tmXBLabels = xblab
  end if
  if (maximize) then
    mres = True
    mres@gsnDraw = False
    mres@gsnFrame = False
    maximize_output(wks, mres)
  end if
  return(plot)
end
;  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
; #####################################################################
undef("select_region_938")
function select_region_938(region:string, lsdata:numeric)
; Arguments
;    region: a string specifying the region to be selected.
;
; Return value
;    An array with the region boundary as (lon1, lat1, lon2,\
;    lat2, lon3, lat3, lon4, lat4 ect) based on the borders
;    given through select_region1
;    with the name of the region as a string attribute @name.
;
; Description
;    Translates a region specification into lat/lon boundaries and a region
;    name as an attribute for non SREX regions
;
; Modification history
begin

  funcname = "select_region_938"
  scriptname = "diag_scripts/lib/ncl/latlon.ncl"
  verbosity = stringtointeger(getenv("ESMValTool_verbosity"))

  region_array_hlp = select_region1(region)
  region_array = (/(/region_array_hlp(2),   region_array_hlp(0)/), \
                   (/region_array_hlp(2),   region_array_hlp(1)/), \
                   (/region_array_hlp(3),   region_array_hlp(1)/), \
                   (/region_array_hlp(3),   region_array_hlp(0)/)/)

  lat = lsdata&lat
  lon = lsdata&lon
  inout = new((/dimsizes(lsdata)/), "float")
  in = 0
  inn = 0
  do ilat = 0, dimsizes(lat) - 1
    do ilon = 0, dimsizes(lon) - 1
      if (dimsizes(region_array(:, 0)).le.2) then
        if ((lat(ilat).ge.region_array(0, 1)).and.(lat(ilat)\
            .lt.region_array(1, 1)).and.(lon(ilon).ge.region_array(0, 0))\
            .and.(lon(ilon).lt.region_array(1, 0))) then
          inou = "True"
          if (inou.eq."True") then
            inout(ilat, ilon) = 1.
            in = in + 1
          end if
          delete(inou)
          inn = inn + 1
        end if
      else
        inou = gc_inout(lat(ilat), lon(ilon), region_array(:, 1),\
                        region_array(:, 0))

        if(inou.eq."True") then
          inout(ilat, ilon) = 1.
          in = in + 1
          delete(inou)
        end if
        inn = inn + 1
      end if
    end do
  end do

  return(inout)
end
; °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
;  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
; #####################################################################
undef("select_region_srex_poly")
function select_region_srex_poly(region:string)
;
; Arguments
;    region: a string specifying the region to be selected.
;
; Return value
;    An array with the position for the string in the region boundary as
;    (lon, lat)
;
; Description
;    Position for the name of the region on a map.
;
; Modification history
begin

  funcname = "select_region_srex_poly"
  scriptname = "diag_scripts/lib/ncl/latlon.ncl"
  verbosity = stringtointeger(getenv("ESMValTool_verbosity"))

  if (region.eq."ALAs") then
    region_array = (/200.000,   68.000/)
  end if
  if (region.eq."AMZs") then
    region_array = (/285.,   0.000/)
  end if
  if (region.eq."CAMs") then
    region_array = (/265.,  15./)
  end if
  if (region.eq."CASs") then
    region_array = (/62.000,   40.000/)
  end if
  if (region.eq."CEUs") then
    region_array = (/ 15.000,    52.000/)
  end if
  if (region.eq."CGIs") then
    region_array = (/265.000,  65.000/)
  end if
  if (region.eq."CNAs") then
    region_array = (/255.000,  40.000/)
  end if
  if (region.eq."EAFs") then
    region_array = (/30.000,    -5.0/)
  end if
  if (region.eq."EASs") then
    region_array = (/110.000,   30.000/)
  end if
  if (region.eq."ENAs") then
    region_array = (/280.000,  30.000/)
  end if
  if (region.eq."MEDs") then
    region_array = (/355.000,  38.000/)
  end if
  if (region.eq."NASs") then
    region_array = (/50.000,    72.000/)
  end if
  if (region.eq."NAUs") then
    region_array = (/120.000,    -20.000/)
  end if
  if (region.eq."NEBs") then
    region_array = (/315.000,    -6.000/)
  end if
  if (region.eq."NEUs") then
    region_array = (/355.000,   65.000/)
  end if
  if (region.eq."SAFs") then
    region_array = (/320.000,    -20.000/)
  end if
  if (region.eq."SAHs") then
    region_array = (/350.000,  22.000/)
  end if
  if (region.eq."SASs") then
    region_array = (/70.000,   10.000/)
  end if
  if (region.eq."SAUs") then
    region_array = (/120.000,    -42.000/)
  end if
  if (region.eq."SSAs") then
    region_array = (/295.0,    -40.000/)
  end if
  if (region.eq."SEAs") then
    region_array = (/105.000,     -0.000/)
  end if
  if (region.eq."TIBs") then
    region_array = (/85.000,  40.000/)
  end if
  if (region.eq."WAFs") then
    region_array = (/350.000,   0.0/)
  end if
  if (region.eq."WASs") then
    region_array = (/45.000,  25.000/)
  end if
  if (region.eq."WSAs") then
    region_array = (/280.0,   -30.0/)
  end if
  if (region.eq."WNAs") then
    region_array = (/233.000,  35.0/)
  end if
    if (region.eq."Arctic") then
    region_array = (/300.000,  75.000/)
  end if
  if (region.eq."Antarctic") then
    region_array = (/300.000,   -75.000/)
  end if
  if (region.eq."Caribbean") then
    region_array = (/290.0,   15.0/)
  end if
  if (region.eq."WesternIndianOcean") then
    region_array = (/60.000,   -20.000/)
  end if
  if (region.eq."NorthernIndianOcean") then
    region_array = (/70.000,  10.000/)
  end if
  if (region.eq."NorthernTropicalPacific") then
    region_array = (/175.000,  10.000/)
  end if
  if (region.eq."EquatorialTropicalPacific") then
    region_array = (/185.000,   0.000/)
  end if
  if (region.eq."SouthernTropicalPacific") then
    region_array = (/195.000,   -10.000/)
  end if
  if (region.eq."World") then
    region_array = (/300.000,   0.000/)
  end if

  return(region_array)

end

;  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
undef("select_region_srex")
function select_region_srex(region:string, lsdata:numeric)
; Arguments
;    region: a string specifying the region to be selected.
;
; Return value
;    An array with the region boundary as (lon1, lat1, lon2,\
;    lat2, lon3, lat3, lon4, lat4 ect)
;    with the name of the region as a string attribute @name.
;
; Description
;    Translates a region specification into lat/lon boundaries and a region
;    name as an attribute.
;
; Modification history
begin

  funcname = "select_region_srex"
  scriptname = "diag_scripts/lib/ncl/latlon.ncl"
  verbosity = stringtointeger(getenv("ESMValTool_verbosity"))

    if (region.eq."ALAs") then
    region_array = (/(/255.000,   60.000/), \
                     (/191.978,   60.000/), \
                     (/191.978,   72.554/), \
                     (/255.000,   72.554/)/)
  end if
  if (region.eq."AMZs") then
    region_array = (/(/293.623,   - 20.000/), \
                     (/280.271,   - 1.239/), \
                     (/291.2,    11.439/), \
                     (/310.00,   11.439/), \
                     (/310.00,    - 20.00/)/)
  end if
  if (region.eq."CAMs") then
    region_array = (/(/291.2,  11.439/), \
                     (/280.271,     -1.239/), \
                     (/241.677,    28.566/), \
                     (/269.685,    28.566/)/)
  end if
  if (region.eq."CASs") then
    region_array = (/(/60.000,   30.000/), \
                     (/60.000,   50.000/), \
                     (/75.000,   50.000/), \
                     (/75.000,   30.000/)/)
  end if
  if (region.eq."CEUs") then
    region_array = (/(/ -10.000,    45.000/), \
                     (/ -10.000,    48.000/), \
                     (/ 40.000,    61.320/), \
                     (/ 40.000,    45.000/)/)
  end if
  if (region.eq."CGIs") then
    region_array = (/(/350.000,  50.000/), \
                     (/255.000,  50.000/), \
                     (/255.000,  85.000/), \
                     (/350.000,  85.000/)/)
  end if
  if (region.eq."CNAs") then
    region_array = (/(/275.000,  50.000/), \
                     (/275.000,  28.566/), \
                     (/255.000,  28.566/), \
                     (/255.000,  50.000/)/)
  end if
  if (region.eq."EAFs") then
    region_array = (/(/25.000,    -11.365/), \
                     (/25.000,   15.000/), \
                     (/51.990,   15.000/), \
                     (/51.990,    -11.365/)/)
  end if
  if (region.eq."EASs") then
    region_array = (/(/100.000,   20.000/), \
                     (/100.000,   50.000/), \
                     (/145.000,   50.000/), \
                     (/145.000,   20.000/)/)
  end if
  if (region.eq."ENAs") then
    region_array = (/(/300.000,  25.000/), \
                     (/275.000,  25.000/), \
                     (/275.000,  50.000/), \
                     (/300.000,  50.000/)/)
  end if
  if (region.eq."MEDs") then
    region_array = (/(/350.000,  30.000/), \
                     (/350.000,  45.000/), \
                     (/40.000,  45.000/), \
                     (/40.000,  30.000/)/)
  end if
  if (region.eq."NASs") then
    region_array = (/(/40.000,    50.000/), \
                     (/40.000,    70.000/), \
                     (/180.000,   70.000/), \
                     (/180.000,   50.000/)/)
  end if
  if (region.eq."NAUs") then
    region_array = (/(/110.000,    -30.000/), \
                     (/110.000,    -10.000/), \
                     (/155.000,    -10.000/), \
                     (/155.000,    -30.000/)/)
  end if
  if (region.eq."NEBs") then
    region_array = (/(/326.000,    -20.000/), \
                     (/310.000,    -20.000/), \
                     (/310.000,    0.000/), \
                     (/326.000,    0.000/)/)
  end if
  if (region.eq."NEUs") then
    region_array = (/(/350.000,   48.000/), \
                     (/350.000,   75.000/), \
                     (/40.000,   75.000/), \
                     (/40.000,   61.320/)/)
  end if
  if (region.eq."SAFs") then
    region_array = (/(/350.000,    -35.000/), \
                     (/350.000,    -11.365/), \
                     (/308.01,     -11.365/), \
                     (/308.01,     -35.000/)/)
  end if
  if (region.eq."SAHs") then
    region_array = (/(/340.000,  15.000/), \
                     (/340.000,  30.000/), \
                     (/ 40.000,  30.000/), \
                     (/ 40.000,  15.000/)/)
  end if
  if (region.eq."SASs") then
    region_array = (/(/60.000,   5.000/), \
                     (/60.000,   30.000/), \
                     (/100.000,  30.000/), \
                     (/100.000,  20.000/), \
                     (/95.000,   20.000/), \
                     (/95.000,   5.000/)/)
  end if
  if (region.eq."SAUs") then
    region_array = (/(/110.000,    -50.000/), \
                     (/110.000,    -30.000/), \
                     (/180.000,    -30.000/), \
                     (/180.000,    -50.000/)/)
  end if
  if (region.eq."SSAs") then
    region_array = (/(/320.624,    -20.000/), \
                     (/320.624,    -56.704/), \
                     (/292.652,    -56.704/), \
                     (/287.859,    -50.000/), \
                     (/293.623,    -20.000/)/)
  end if
  if (region.eq."SEAs") then
    region_array = (/(/95.000,     -10.000/), \
                     (/95.000,    20.000/), \
                     (/155.000,   20.000/), \
                     (/155.000,    -10.000/)/)
  end if
  if (region.eq."TIBs") then
    region_array = (/(/75.000,  30.000/), \
                     (/75.000,  50.000/), \
                     (/100.000,  50.000/), \
                     (/100.000,  30.000/)/)
  end if
  if (region.eq."WAFs") then
    region_array = (/(/340.000,   -11.365/), \
                     (/340.000,  15.000/), \
                     (/25.000,   15.000/), \
                     (/25.000,    -11.365/)/)
  end if
  if (region.eq."WASs") then
    region_array = (/(/40.000,  15.000/), \
                     (/40.000,  50.000/), \
                     (/60.000,  50.000/), \
                     (/60.000,  15.000/)/)
  end if
  if (region.eq."WSAs") then
    region_array = (/(/280.271,   -1.239/), \
                     (/293.623,   -20.000/), \
                     (/287.859,   -50.000/), \
                     (/292.652,   -56.704/), \
                     (/277.978,   -56.704/), \
                     (/277.978,  0.530/)/)
  end if

  if (region.eq."WNAs") then
    region_array = (/(/255.000,  28.566/), \
                     (/230.000,  28.566/), \
                     (/230.000,  60.000/), \
                     (/255.000,  60.000/)/)
  end if

  if (region.eq."Arctic") then
    region_array = (/(/0.000,   67.500/), \
                     (/360.000,  90.000/)/)

  end if
  if (region.eq."Antarctic") then
    region_array = (/(/0.000,   -90.000/), \
                     (/360.000,   -50.000/)/)
  end if
  if (region.eq."Caribbean") then
    region_array = (/(/291.200,   11.400/), \
                     (/274.200,   25.000/), \
                     (/300.000,   25.000/), \
                     (/300.000,   11.400/)/)
  end if
  if (region.eq."WesternIndianOcean") then
    region_array = (/(/52.000,   -25.000/), \
                     (/75.000,  5.000/)/)
  end if
  if (region.eq."NorthernIndianOcean") then
    region_array = (/(/60.000,  5.000/), \
                     (/95.000,  30.000/)/)
  end if
  if (region.eq."NorthernTropicalPacific") then
    region_array = (/(/155.000,  5.000/), \
                     (/210.000,  25.000/)/)
  end if
  if (region.eq."EquatorialTropicalPacific") then
    region_array = (/(/155.000,   -5.000/), \
                     (/230.000,  5.000/)/)
  end if
  if (region.eq."SouthernTropicalPacific") then
    region_array = (/(/155.000,   -25.000/), \
                     (/230.000,   -5.000/)/)
  end if
  if (region.eq."World") then
    region_array = (/(/0.000,   -90.000/), \
                     (/360.000,  90.000/)/)
  end if

  lat = lsdata&lat
  lon = lsdata&lon
  inout = new((/dimsizes(lsdata)/), "float")
  in = 0
  inn = 0
  do ilat = 0, dimsizes(lat) - 1
    do ilon = 0, dimsizes(lon) - 1
      if (dimsizes(region_array(:, 0)).le.2) then
        if ((lat(ilat).ge.region_array(0, 1)).and.(lat(ilat)\
            .lt.region_array(1, 1)).and.(lon(ilon).ge.region_array(0, 0))\
            .and.(lon(ilon).lt.region_array(1, 0))) then
          inou = "True"
          if (inou.eq."True") then
            inout(ilat, ilon) = 1.
            in = in + 1
          end if
          delete(inou)
          inn = inn + 1
        end if
      else
        inou = gc_inout(lat(ilat), lon(ilon), region_array(:, 1),\
                        region_array(:, 0))

        if(inou.eq."True") then
          inout(ilat, ilon) = 1.
          in = in + 1
          delete(inou)
        end if
        inn = inn + 1
      end if
    end do
  end do

  return(inout)
end
; °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
undef("add_legenda_page")
procedure add_legenda_page(wks, MODEL:string, colors:string, dashs:numeric,\
                           scatters:string, ticknesses:numeric, place:string)
local MODEL, dim_mod, colors, dashs, edge_x, edge_y, edge_res, tx1res, \
      dim_lb, xleg, xl, yleg, yl, xtxt, ytxt, place, dim_lines,\
      space_lines, lin_sp, scatters, ticknesses
begin
;  *  * * Create Legend * * * * * * * * * * * * * * * * * * * * * * * *
  dim_mod = dimsizes(MODEL)
  dim_lines = floattoint(dim_mod/3) + 1
  if ((mod(dim_mod, 3)).eq.0) then
    dim_lines = floattoint(dim_mod/3) + 1
  else
    dim_lines = floattoint(dim_mod/3) + 2
  end if
  lin_sp = 0.025
  space_lines = (dim_lines + 1) * lin_sp

  if (place.eq."bottom")
    top_line = 0.20
  end if
  if (place.eq."middle")
    top_line = 0.42
  end if
  if (place.eq."top")
    top_line = 0.99
  end if
  if ((place.ne."bottom").and.(place.ne."middle").and.(place.ne."top")) then
    top_line = tofloat(place)
  end if
  bot_line = top_line - space_lines
  ; edge_x = (/ 0.12, 0.92, 0.92, 0.12, 0.12/)

  edge_x = (/ 0.22, 0.76, 0.76, 0.22, 0.22/)
  edge_y = (/space_lines + bot_line, space_lines + bot_line,\
            bot_line, bot_line, space_lines + bot_line/)
  edge_res = True             ; Indicate you want to set some resources.
  edge_res@gsLineColor = "white"
  gsn_polyline_ndc(wks, edge_x, edge_y, edge_res)

; * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  dim_lb = dimsizes(MODEL)

  gs1res = True
  tx1res = True
  tx1res@txFontHeightF = lin_sp / 5 * 2

  xleg = new((/3 * dim_lines, 4/), "float")
  xl = new((/3, 4/), "float")

  do i = 0, 2
    xl(0, :) = fspan((edge_x(0) + 0.02), (edge_x(0) + 0.07), 4)
    xl(1, :) = fspan(((edge_x(2) + edge_x(0))/2 - 0.09),\
                     ((edge_x(2) + edge_x(0))/2 - 0.04), 4)
    xl(2, :) = fspan((edge_x(2) - 0.20), (edge_x(2) - 0.15), 4)
  end do

  nt = 0
  do j = 0, dim_lines - 1
    do i = 0, 2
      xleg(nt, :) = (/xl(i, :)/)
      nt = nt + 1
    end do
  end do

  yleg = new((/3 * dim_lines, 4/), "float")
  yl = new(dim_lines, "float")

  do i = 0, dim_lines - 1
    yl(dim_lines - 1 - i) = edge_y(3) + (((edge_y(0) - \
                                         (edge_y(3)))/dim_lines) * \
                                         (i + 1)) - (lin_sp)
  end do

  nt = 0
  do j = 0, dim_lines - 1
    do i = 0, 2
      yleg(nt, :) = (/yl(j)/)
      nt = nt + 1
    end do
  end do

  xtxt = new((/dim_lines * 3/), "float")
  ytxt = new((/dim_lines * 3/), "float")

  nt = 0
  do j = 0, dim_lines - 1
    do i = 0, 2
      xtxt(nt) = xl(i, 3) + lin_sp/5.
      nt = nt + 1
    end do
  end do
  nt = 0
  do i = 0, dim_lines - 1
    do j = 0, 2
      ytxt(nt) = (/yl(i)/)
      nt = nt + 1
    end do
  end do

  do i = 0, (dimsizes(MODEL) - 1)
    gs1res@gsLineColor = colors(i)
    gs1res@gsLineDashPattern = dashs(i)
    gs1res@gsLineThicknessF = ticknesses(i)
    tx1res@txFontColor = colors(i)
    tx1res@txJust = "CenterLeft"             ; change justification
    if (scatters(i).eq."Markers")
      gs1res@gsMarkerColor = colors(i)
      gs1res@gsMarkerIndex = dashs(i)
      gs1res@gsMarkerSizeF = 0.01
      gs1res@gsMarkerThicknessF = 1.5
      gsn_polymarker_ndc(wks, xleg(i, 2), yleg(i, 2), gs1res)
    else

      gsn_polyline_ndc(wks, xleg(i, :), yleg(i, :), gs1res)
    end if

    gsn_text_ndc(wks, MODEL(i), xtxt(i), ytxt(i), tx1res)
  end do

end

; ---------------------------------------------------------------------
; The INTERSECTION of two sets is the set of elements which are in both sets.
; For example: let A = (1,2,3) and B = (3,4,5).
; The INTERSECTION of A and B, written A intersection B = (3).
; ---------------------------------------------------------------------
undef("get_unique_intersection")
function get_unique_intersection(a, b)
local atype, btype, au, bu, nau, nbu, abi, k, n, m
begin
  atype = typeof(a)
  btype = typeof(b)
  if (atype.ne.btype) then
      print("get_unique_intersection: arguments a and b must be the same type")
      print("                         typeof(a)="+atype)
      print("                         typeof(b)="+btype)
      exit
  end if
  ; 'get_unique_values' sorts the elements
  ;  and returns a one-dimensional arry
  au  = get_unique_values(a)   ; avoid duplicate 'a' entries
  bu  = get_unique_values(b)   ; avoid duplicate 'b' entries

  nau = dimsizes(au)
  nbu = dimsizes(bu)
  ; reserve space
  abi = new(min((/nau, nbu/)), typeof(a), default_fillvalue(atype))
  k   = -1

  na  = 0                                 ; MH suggested loop
  nb  = 0
  do while(na.lt.nau.and.nb.lt.nbu)
    if (bu(nb).gt.au(na)) then
        na = na + 1
    else
        if (au(na).gt.bu(nb)) then
            nb = nb + 1
        else                              ; they are equal
            k = k+1
            abi(k) = au(na)
            nb = nb + 1
            na = na + 1
        end if
    end if
  end do

  if (k.eq.-1) then
      abi    := new(1, typeof(a), default_fillvalue(atype))
      abi@tag = "get_unique_intersection; NULL set; no intersection pairs"
      return(abi)
  end if

  abi@tag = "NCL: get_unique_intersection has created the result"
  return(abi(0:k))
end

; ---------------------------------------------------------------------
; The DIFFERENCE of two sets is the set of elements,
; which are unique to each set.
; For example: let A = (1,2,3) and B = (3,4,5).
; The DIFFERENCE of A and B, written A difference B = (1,2,4,5).
; ---------------------------------------------------------------------
undef("get_unique_difference")
function get_unique_difference(a, b)
local atype, btype, abi, nabi, au, bu, abd, AB, nAB, n
begin
  atype = typeof(a)
  btype = typeof(b)
  if (atype.ne.btype) then
      print("get_unique_difference: arguments a and b must be the same type")
      print("                       typeof(a)=" + atype)
      print("                       typeof(b)=" + btype)
      exit
  end if
  ; 'get_unique_values' sorts the elements
  ;  and returns a one-dimensional arry
  abi  = get_unique_intersection(a, b)    ; intersect (common) values
  if (ismissing(abi(0))) then             ; NO intersect
      ; 'get_unique_values' sorts the elements
      au  = get_unique_values(a)          ; unique 'a' [au]
      bu  = get_unique_values(b)          ; unique 'b' [bu]
      ; all unique valuesl [d] difference
      abd = array_append_record(au, bu, 0)
      delete([/au, bu /])
      ; must sort again cuz 'ab' may not be in order
      if (atype.eq."string") then
          sqsort(abd)
      else
          qsort(abd)
      end if
      abd@tag = "NCL: get_unique_difference; arrays are DISJOINT; " + \
                "no intersection"
      return(abd)
  end if
  nabi = dimsizes(abi)

  au   = get_unique_values(a)   ; 'a' [u]nique entries
  bu   = get_unique_values(b)   ; 'b' [u]nique entries
  nau  = dimsizes(au)
  nbu  = dimsizes(bu)
  ; reserve space
  abd  = new(nau + nbu,  typeof(a), default_fillvalue(atype))
  ABu  = array_append_record(au, bu, 0)    ; unique AB
  delete([/au, bu /])
  nABu = dimsizes(ABu)

  abd  = new(nABu, typeof(a), default_fillvalue(atype))
  k    = -1
  do n = 0, nABu - 1
    if (.not.any(abi.eq.ABu(n))) then
      k = k+1
      abd(k) = ABu(n)
    end if
  end do
  ; must sort again cuz 'ab' may not be in order
  if (atype.eq."string") then
      sqsort(abd(0:k))
  else
      qsort(abd(0:k))
  end if
  abd@tag = "NCL: get_unique_difference has created the result"
  return(abd(0:k))
end

; ######################################################################
