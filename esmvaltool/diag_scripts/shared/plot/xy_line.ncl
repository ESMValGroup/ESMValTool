; #############################################################################
; PLOT SCRIPTS FOR XY LINE PLOTS
; #############################################################################
; Please consider using or extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
; Contents:
;
;    function profile_plev
;    function aerosol_profile
;    function aerosol_sizedist
;    procedure xy_line
;    procedure xy_line_anom
;    function timeseries_station
;    function cycle_plot
;    function errorbar_plot
;
; #############################################################################

load "$diag_scripts/shared/set_operators.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"

; #############################################################################
undef("profile_plev")
function profile_plev(wks_in[1],
                      source,
                      varname[1]: string)
;
; Arguments
;    wks_in: workstations (graphic object or default will be used).
;    source: data to be plotted or a NetCDF filename:
;            Must have "plev" dimension & can have "models", "quantity".
;              @ptop: controls pressure axis range (Default: full range).
;              @pbot: controls pressure axis range (Default: full range).
;              @zoom: controls x axis range (Default: "yes").
;              @font: font type to use
;              @Refmodel: reference model (Default: first in models dimension)
;              @long_name: long variable name (Default: var)
;              @short_name: short variable name (Default: var)
;              @units: variable units (Default: missing)
;    varname: variable name, needed for netCDF files with multiple variables
;
; Source prototype
;
; Return value
;    A graphic variable.
;
; Description
;    Creates a plot of profile(s) with vertical pressure axis.
;    Adds percentiles or stddev of first model (if available) as whiskers.
;    Opens default wks, if not provided as argument of type "graphic".
;    Defines default ressources, which are overridden by argument res.
;    Creates plot, according to wks & res.
;
; Caveats
;    Treatment of coordinate "quantity" not very general yet.
;    Selection of defaults for res almost arbitrary.
;    Please check results of all scripts that use this routine if
;    modifying the defaults!
;
; Modification history
;    20140214-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, diag_script, data, var, font, font_size_fac, \
  i_ref, FuncCode, i_N, i_mean, i_stddev, i_min, i_max, i_5, i_25, \
  i_median, i_75, i_95, target, i_quantity, data2plot, i_plev, ymin, ymax, \
  xmin, xmax, colors, dashes, thicks, annots, avgstd, res, res_whiskers, \
  plot, lgLabels, defaults, zoom, left, right, ixmin, ixmax, res_markers, \
  xmk, ymk, border, info, short_name, long_name, temp
begin

  funcname = "profile_plev"
  scriptname = "diag_scripts/shared/plot/xy_line.ncl"
  enter_msg(scriptname, funcname)

  ; Set hardcoded parameters
  border = 0.04  ; factor to extent y- & x-range (except @ top)
  font_size_fac = 1.5  ; factor to change all font sizes (w.r.t. to defaults)
  i_ref = 0  ; hardcode first model as reference -> generalize, if needed

  ; Get data, either directly or via netCDF file
  if (typeof(source) .eq. "string") then
    data = ncdf_read(source, varname)  ; can deal with defaults
  else
    data = source
    copy_VarMeta(source, data)
  end if

  ; Retrieve some metadata from data
  defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
  if (any(varname .eq. defaults)) then  ; Use default
    var = att2var(data, "var")
  else                                  ; Use explicitely given name
    var = varname
  end if
  diag_script = att2var(data, "diag_script")
  ptop = att2var_default(data, "ptop", "full", 0)
  pbot = att2var_default(data, "pbot", "full", 0)
  font = att2var_default(data, "font", 25, 0)  ; default: Times New Roman
  long_name = att2var_default(data, "long_name", var, 0)
  short_name = att2var_default(data, "short_name", var, 0)
  zoom = att2var_default(data, "zoom", "yes", 0)
  if (dimsizes(zoom).eq.2) then
    temp = zoom
    delete(zoom)
    zoom = "yes"
    zoom@range = tofloat(temp)
    delete(temp)
  end if

  ; Position of reference model
  if (isdim(data, "models")) then
    Refmodel = att2var_default(data, "Refmodel", "missing", 0)
    if (Refmodel.ne."missing") then
      i_ref = ind(data&models.eq.Refmodel)
      if (all(ismissing(i_ref))) then
        error_msg("w", scriptname, funcname, "Refmodel " + Refmodel + \
                  "not found, using first")
        i_ref = 0
      else
        if (dimsizes(i_ref).ne.1) then
          error_msg("f", scriptname, funcname, "Refmodel " + \
                    Refmodel + "not unique")
        end if
      end if
    end if
  end if

  ; determine handling of legend
  legend = att2var(data, "legend")  ; mode for legend: inline/separate/none
  valid_legends = (/"inline", "separate"/)
  if (any(legend .eq. valid_legends)) then
    if (.not. iscoord(data, "models")) then
      log_debug("warning: No 'models' coordinate found.")
      log_debug("warning:  -> No legend in profile_plev")
    end if
  else
    legend = "none"
  end if

  ; Check if a valid wks has been provided, otherwise invoke default
  wks = get_wks(wks_in, diag_script, var)  ; Function in aux_plotting.ncl

  if (isdim(data, "quantity")) then
    i_N      = ind(data&quantity.eq."N")
    i_mean   = ind(data&quantity.eq."mean")
    i_stddev = ind(data&quantity.eq."stddev")
    i_min    = ind(data&quantity.eq."min")
    i_max    = ind(data&quantity.eq."max")
    i_5      = ind(data&quantity.eq."5")
    i_25     = ind(data&quantity.eq."25")
    i_median = ind(data&quantity.eq."median")
    i_75     = ind(data&quantity.eq."75")
    i_95     = ind(data&quantity.eq."95")
    if (ismissing(i_mean)) then
      ; "median" is evaluated and plotted, if no "mean" is there
      target = "median"  ; 2nd priority
      i_quantity = (/i_5, i_25, i_median, i_75, i_95/)
    else
      target = "mean"  ; 1st priority
      i_quantity = (/i_min, i_stddev, i_mean, i_stddev, i_max/)
    end if
  end if

  ; Extract data to plot & check for missing values
  ; in preparation to determine vertical range
  if (isdim(data, "quantity") .and. isdim(data, "models")) then
    data2plot = data(models|:, plev|:, quantity|i_quantity)
    if (target.eq."mean") then  ; mean +/- stddev
      data2plot(:, :, 1) = data2plot(:, :, 2) - data2plot(:, :, 1)
      data2plot(:, :, 3) = data2plot(:, :, 2) + data2plot(:, :, 3)
    end if
    i_plev = ind(.not. ismissing(data2plot(i_ref, :, 2)))
  else if (isdim(data, "quantity") .and. .not.isdim(data, "models")) then
    data2plot = data(plev|:, quantity|i_quantity)
    if (target.eq."mean") then  ; mean +/- stddev
      data2plot(:, 1) = data2plot(:, 2) - data2plot(:, 1)
      data2plot(:, 3) = data2plot(:, 2) + data2plot(:, 3)
    end if
    i_plev = ind(.not. ismissing(data2plot(:, 2)))
  else if (isdim(data, "models") .and. .not.isdim(data, "quantity")) then
    data2plot = data(models|:, plev|:)
    i_plev = ind(.not. ismissing(data2plot(i_ref, :)))
  else
    data2plot = data(plev|:)
    i_plev = ind(.not. ismissing(data2plot))
  end if
  end if
  end if

  ; y-axis range
  if (ptop.ne."full") then
    if (ptop .eq. "auto") then
      ymin = min(data2plot&plev(i_plev))
    else
      ymin = tofloat(ptop) * 100.  ; convert to Pa
    end if
  else
    ymin = min(data2plot&plev)
  end if
  if (pbot.ne."full") then
    if (pbot .eq. "auto") then
      ymax = max(data2plot&plev(i_plev))
    else
      ymax = tofloat(pbot) * 100.  ; convert to Pa
    end if
  else
    ymax = max(data2plot&plev)
  end if
  delete(i_plev)  ; update i_plev
  i_plev = ind(data2plot&plev.ge.ymin .and. data2plot&plev.le.ymax)

  ; Determine min & max values of x-axis, and arrow positions
  ; to indicate cutoff.
  if (isdim(data, "quantity") .and. isdim(data, "models")) then
    if (zoom.eq."yes") then
      if (isatt(zoom, "range")) then
        xmin = zoom@range(0)
        xmax = zoom@range(1)
      else
        xmin = min(data2plot(:, i_plev, 2))
        xmax = max(data2plot(:, i_plev, 2))
      end if
      ixmin = ind(data2plot(i_ref, :, 0).lt.xmin .or.\
                  data2plot(i_ref, :, 1).lt.xmin)  ; for arrows
      ixmax = ind(data2plot(i_ref, :, 4).gt.xmax .or.\
                  data2plot(i_ref, :, 3).gt.xmax)  ; to indicate cutoff
    else
      ixmin = new(1, "integer")  ; missing: no arrows
      ixmax = new(1, "integer")  ; to indicate cutoff
      xmin = min(data2plot(:, i_plev, :))
      xmax = max(data2plot(:, i_plev, :))
    end if
  else if (isdim(data, "quantity") .and. .not.isdim(data, "models")) then
    if (zoom.eq."yes") then
      if (isatt(zoom, "range")) then
        xmin = zoom@range(0)
        xmax = zoom@range(1)
      else
        xmin = min(data2plot(i_plev, 2))
        xmax = max(data2plot(i_plev, 2))
      end if
      ixmin = ind(data2plot.lt.xmin)  ; for arrows
      ixmax = ind(data2plot.gt.xmax)  ; to indicate cutoff
    else
      ixmin = new(1, "integer")  ; missing: no arrows
      ixmax = new(1, "integer")  ; to indicate cutoff
      xmin = min(data2plot(i_plev, :))
      xmax = max(data2plot(i_plev, :))
    end if
  else if (isdim(data, "models") .and. .not.isdim(data, "quantity")) then
    if (zoom.eq."yes" .and. isatt(zoom, "range")) then
      xmin = zoom@range(0)
      xmax = zoom@range(1)
      ixmin = ind(data2plot.lt.xmin)  ; for arrows
      ixmax = ind(data2plot.gt.xmax)  ; to indicate cutoff
    else
      ixmin = new(1, "integer")  ; missing: no arrows
      ixmax = new(1, "integer")  ; to indicate cutoff
      xmin = min(data2plot(:, i_plev))
      xmax = max(data2plot(:, i_plev))
    end if
  else
    if (zoom.eq."yes" .and. isatt(zoom, "range")) then
      xmin = zoom@range(0)
      xmax = zoom@range(1)
      ixmin = ind(data2plot.lt.xmin)  ; for arrows
      ixmax = ind(data2plot.gt.xmax)  ; to indicate cutoff
    else
      ixmin = new(1, "integer")  ; missing: no arrows
      ixmax = new(1, "integer")  ; to indicate cutoff
      xmin = min(data2plot(i_plev))
      xmax = max(data2plot(i_plev))
    end if
  end if
  end if
  end if

  ; Get default style parameters (diag_scripts/shared/plot/style.ncl)
  ; NOTE: These may or may not contain the reference model
  if (.not.isatt(diag_script_info, "styleset")) then
    diag_script_info@styleset = "DEFAULT"
  end if
  ; FIX-ME: project_style function has changed
  ; colors = project_style(diag_script_info, "colors")
  ; dashes = project_style(diag_script_info, "dashes")
  ; thicks = tofloat(project_style(diag_script_info, "thicks"))
  ; annots = project_style(diag_script_info, "annots")
  ; avgstd = project_style(diag_script_info, "avgstd")

  ; Set default ressources
  res = True
  res@gsnDraw           = True  ; Draw
  res@gsnFrame          = False  ; Advance frame
  res@tiMainString      = "Output of " + diag_script  ; Main title
  res@gsnLeftString     = ""  ; Annotation
  res@trYReverse        = True  ; reverse Y-axis
  res@trYAxisType       = "LogAxis"
  res@tmYLMode          = "Explicit"  ; explicit labels
  res@tmYLValues        = (/1000, 700, 500, 400, 300, 200, \
                           100, 50, 30, 10, 5, 1/) * 100  ; [Pa]
  res@tmYLLabels        = "" + res@tmYLValues / 100.  ; make strings [hPa]
  res@tiYAxisString     = "Pressure / hPa"
  res@xyMonoDashPattern = False
  res@tmEqualizeXYSizes = True
  res@txFontQuality     = "High"
  res@txFont            = font
  res@txFontHeightF     = 0.02  ; for subheader
  res@tiMainFont        = font
  res@tiXAxisFont       = font
  res@tiYAxisFont       = font
  res@tmXBLabelFont     = font
  res@tmYLLabelFont     = font
  res@tiMainFontHeightF    = 0.025 * font_size_fac  ; set font_size_fac above
  res@tiXAxisFontHeightF   = 0.02 * font_size_fac
  res@tiYAxisFontHeightF   = 0.02 * font_size_fac
  res@tmYLLabelFontHeightF = 0.02 * font_size_fac
  res@tmXBLabelFontHeightF = 0.02 * font_size_fac
  res@tmYUseLeft        = True
  res@trXMinF           = xmin - border * (xmax - xmin)
  res@trXMaxF           = xmax + border * (xmax - xmin)
  res@trYMaxF           = ymax + border * (ymax - ymin)
  ; Do not limit top to show models reaching higher
  ; One or more model lines are cut at res@trYMinF
  ;     for some mysterious reason -> do not consider border here.
  res@trYMinF           = ymin  ; - border * (ymax - ymin)

  ; Styles
  res@xyLineColors      = colors
  res@xyLineThicknesses = thicks
  res@xyDashPatterns    = dashes

  ; tiXAxisString
  res@tiXAxisString = short_name
  if (isatt(data, "units")) then
    res@tiXAxisString = \
      res@tiXAxisString + " [" + format_units(data@units) + "]"
  end if

  ; inline legend
  if (legend.eq."inline") then
    lgLabels = data&models
    res@pmLegendDisplayMode    = "Always"  ; turn on legend
    res@pmLegendSide           = "Right"
    res@lgJustification = "CenterLeft"
    res@lgLabelFontHeightF     = .02
    res@pmLegendHeightF        = \
      res@lgLabelFontHeightF * (1 + dimsizes(lgLabels))
    res@pmLegendWidthF         = 0.05  ; length of the lines
    res@pmLegendOrthogonalPosF = \
      -0.027 * max(strlen(lgLabels)) - res@pmLegendWidthF * font_size_fac
    res@lgLabelFont            = font
    res@lgLabelFontHeightF     = .02  ; change font height
    res@lgPerimOn              = False  ; no box around
    res@xyExplicitLegendLabels = " " + lgLabels
    res@trXMaxF = res@trXMaxF + (res@trXMaxF - res@trXMinF) * \
      abs(res@pmLegendOrthogonalPosF) * font_size_fac
  end if

  ; Override defaults with "res_" attributes of "data"
  ; Function in ~/interface_scripts/auxiliary.ncl
  res_new = att2var(data, "res_")
  copy_VarMeta(res_new, res)

  ; resources for whiskers
  res_whiskers                   = True
  res_whiskers@gsLineColor       = res@xyLineColors(i_ref)
  res_whiskers@gsLineDashPattern = res@xyDashPatterns(i_ref)
  res_whiskers@gsLineThicknessF  = res@xyLineThicknesses(i_ref) / 2.

  ; resources for arrow heads & refresh markers
  res_markers               = True
  res_markers@gsMarkerSizeF = 0.001 * res@xyLineThicknesses(i_ref)
  res_markers@gsMarkerColor = res@xyLineColors(i_ref)

  ; Plotting
  if (isdim(data, "quantity") .and. isdim(data, "models")) then

    ; mean / median
    plot = gsn_csm_xy(wks, data2plot(:, :, 2), data2plot&plev, res)

    ; sdtddev / 25%-75% (limit for blank border)
    left = where(data2plot(i_ref, :, 1).lt.xmin, xmin, \
                 data2plot(i_ref, :, 1))
    right = where(data2plot(i_ref, :, 3).gt.xmax, xmax, \
                  data2plot(i_ref, :, 3))
    horizontal_whiskers(wks, plot, res_whiskers, left, right, \
                        data2plot&plev)

    ; min-max / 5%-95% (limit for blank border)
    res_whiskers@gsLineDashPattern(:) = 2
    res_whiskers@gsLineThicknessF  = res@xyLineThicknesses(i_ref) / 4.
    left = where(data2plot(i_ref, :, 0).lt.xmin, xmin, \
                 data2plot(i_ref, :, 0))
    right = where(data2plot(i_ref, :, 4).gt.xmax, xmax, \
                  data2plot(i_ref, :, 4))
    horizontal_whiskers(wks, plot, res_whiskers, left, right, \
                        data2plot&plev)
    delete([/left, right/])

    ; arrow heads for cut off whiskers due to zooming
    if (.not. all(ismissing(ixmin))) then
      res_markers@gsMarkerIndex = 10  ; left arrow head
      add_markers(wks, plot, res_markers, \
                  xmin, data2plot&plev(ixmin))
    end if
    if (.not. all(ismissing(ixmax))) then
      res_markers@gsMarkerIndex = 11  ; right arrow head
      add_markers(wks, plot, res_markers, \
                  xmax, data2plot&plev(ixmax))
    end if

    ; refresh mean / median with markers
    res_markers@gsMarkerIndex = 16  ; solid circle
    add_markers(wks, plot, res_markers, \
                data2plot(i_ref, :, 2), data2plot&plev)

    else if (isdim(data, "quantity") .and. .not.isdim(data, "models")) then

      ; mean / median
      plot = gsn_csm_xy(wks, data2plot(:, 2), data2plot&plev, res)

      ; sdtddev / 25%-75% (limit for blank border)
      left = where(data2plot(:, 1).lt.xmin, xmin, data2plot(:, 1))
      right = where(data2plot(:, 3).gt.xmax, xmax, data2plot(:, 3))
      horizontal_whiskers(wks, plot, res_whiskers, left, right, \
                          data2plot&plev)

      ; min-max / 5%-95% (limit for blank border)
      res_whiskers@gsLineDashPattern(:) = 2
      res_whiskers@gsLineThicknessF  = res@xyLineThicknesses(i_ref) / 4.
      left = where(data2plot(:, 0).lt.xmin, xmin, data2plot(:, 0))
      right = where(data2plot(:, 4).gt.xmax, xmax, data2plot(:, 4))
      horizontal_whiskers(wks, plot, res_whiskers, left, right, \
                          data2plot&plev)
      delete([/left, right/])

      ; arrow heads for cut off whiskers due to zooming
      if (.not. all(ismissing(ixmin))) then
        res_markers@gsMarkerIndex = 10  ; left arrow head
        add_markers(wks, plot, res_markers, xmin, data2plot&plev(ixmin))
      end if
      if (.not. all(ismissing(ixmax))) then
        res_markers@gsMarkerIndex = 11  ; right arrow head
        add_markers(wks, plot, res_markers, xmax, data2plot&plev(ixmax))
      end if

      ; refresh mean / median with markers
      res_markers@gsMarkerIndex = 16  ; solid circle
      add_markers(wks, plot, res_markers, data2plot(:, 2), data2plot&plev)

    else

      ; plot everything, if there is no "quantity" coordinate
      plot = gsn_csm_xy(wks, data2plot, data2plot&plev, res)

      ; refresh mean / median with markers
      res_markers@gsMarkerIndex = 16  ; solid circle
      add_markers(wks, plot, res_markers, data2plot, data2plot&plev)

    end if
  end if

  ; outfile name
  if (isatt(wks, "fullname")) then
    plot@outfile = wks@fullname
  else
    plot@outfile = wks@name
    error_msg("w", scriptname, funcname, "wks@fullname missing, " + \
              "consider to use get_wks to open wks")
  end if

  ; confirmation of operation
  if (isatt(res, "gsnDraw")) then
    if (res@gsnDraw .eq. True) then
      log_info("Wrote " + plot@outfile)
    end if
  else
    log_info("Wrote " + plot@outfile)
  end if

  leave_msg(scriptname, funcname)
  return(plot)

end

; #############################################################################
undef("aerosol_profile")
function aerosol_profile(wks_in[1],
                         source_m,
                         source_o,
                         varname[1])
;
; Arguments
;    wks_in: workstations (graphic object or default will be used).
;    source_m: model data to be plotted or a NetCDF filename with data.
;    source_o: observations data to be plotted or a NetCDF filename with data.
;    varname: variable name in the file.
;
; Source_m prototype
;    source[*][*][*]
;    source!0 = model
;    source!1 = statistic
;    source!2 = plev
;
; Source_o prototype
;    source[*][*][*]
;    source!0 = case
;    source!1 = statistic
;    source!2 = plev
;
; Return value
;    A graphic variable.
;
; Description
;    Creates a plot of vertical profile (level vs. data)
;    Plots median, mean or both depending on availability of obervations.
;
; Caveats
;
; Modification history
;    20161013-righi_mattia: moved size distributions to a separate routine.
;    20140917-righi_mattia: renamed to aerosol_profile and extended for
;                           plotting size distributions.
;    20140705-righi_mattia: written.
;
local funcname, scriptname, tmp, colors, annots, res, res_new
begin

  funcname = "aerosol_profile"
  scriptname = "diag_scripts/shared/plot/xy_line.ncl"
  enter_msg(scriptname, funcname)

  ; Get model data, either directly or via netCDF file
  if (typeof(source_m) .eq. "string") then
    model = ncdf_read(source_m, varname)
  else
    model = source_m
    copy_VarMeta(source_m, model)
  end if

  ; Get observational data, either directly or via netCDF file
  if (typeof(source_o) .eq. "string") then
    obser = ncdf_read(source_o, varname)
  else
    obser = source_o
    copy_VarMeta(source_o, obser)
  end if

  ; Retrieve basic metadata from data
  defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
  if (any(varname .eq. defaults)) then
    var = att2var(model, "var")
  else
    var = varname
  end if

  ; Check if a valid wks has been provided, otherwise invoke default
  wks = get_wks(wks_in, diag_script, var)

  ; Define line colors (first color is used for the observations)
  ; FIX-ME: project_style function has changed
  ; colors = project_style(diag_script_info, "colors")
  ; annots = project_style(diag_script_info, "annots")
  gsn_define_colormap(wks, array_append_record((/"white", "black"/), \
                      colors, 0))

  ; Define markers and dashes for median and mean
  marker_mean = 16
  dash_median = 0
  dash_mean = 1
  dash_stddev = 2

  ; Resources
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@trXAxisType = "LogAxis"
  res@trYAxisType = "LogAxis"
  res@tmXBMinorPerMajor = 8
  res@tmXTMinorPerMajor = 8
  res@tmYLMinorPerMajor = 8
  res@tmYRMinorPerMajor = 8
  res@tiXAxisString = varname + " [" + format_units(model@units) + "]"
  res@tiYAxisString = "Pressure [" + format_units(model&plev@units) + "]"
  res@gsnCenterStringFontHeightF = 0.015

  ; Lines
  resL = True
  resL@gsLineThicknessF = 2.0
  resL@gsLineDashSegLenF = 0.25

  ; Markers
  resM = True

  ; Polygons
  resP = True

  ; Text
  resT = True
  resT@txFont = 21

  ; Markers
  resM = True

  ; Override defaults with "res_" attributes of "data"
  res_new = att2var(model, "res_")
  copy_VarMeta(res_new, res)

  ; Set number of models
  nmodels = dimsizes(annots)

  ; Set number of cases
  ncases = dimsizes(obser&case)

  ; Set the index for the various statistics
  i_N      = ind(model&statistic.eq."N")
  i_mean   = ind(model&statistic.eq."mean")
  i_stddev = ind(model&statistic.eq."stddev")
  i_min    = ind(model&statistic.eq."min")
  i_max    = ind(model&statistic.eq."max")
  i_5      = ind(model&statistic.eq."5")
  i_10     = ind(model&statistic.eq."10")
  i_25     = ind(model&statistic.eq."25")
  i_median = ind(model&statistic.eq."median")
  i_75     = ind(model&statistic.eq."75")
  i_90     = ind(model&statistic.eq."90")
  i_95     = ind(model&statistic.eq."95")

  ; Set plot ranges
  modelmin = min((/min(model(:, i_mean, :)), min(model(:, i_median, :))/))
  modelmax = max((/max(model(:, i_mean, :)), max(model(:, i_median, :))/))
  obsermin = min((/min(obser(:, i_mean, :)), min(obser(:, i_median, :))/))
  obsermax = max((/max(obser(:, i_mean, :)), max(obser(:, i_median, :))/))
  res@tmXBMode = "Explicit"
  res@tmYLMode = "Explicit"
  res@trYReverse = True

  ; Set X-axis labels
  res@trXMinF = min((/modelmin, obsermin/)) / 2.
  res@trXMaxF = 2 * max((/modelmax, obsermax/))
  res@tmXBValues = set_log_ticks(res@trXMinF, res@trXMaxF, "major")
  res@tmXBLabels = res@tmXBValues
  res@tmXBMinorValues = set_log_ticks(res@trXMinF, res@trXMaxF, "minor")

  ; Set Y-axis labels
  res@trYMinF = 100
  res@trYMaxF = 1050
  res@tmYLValues = ispan(100, 1000, 100)
  res@tmYLLabels = res@tmYLValues

  ; Check available statistics in the observations
  l_median = .not.(all(ismissing(obser(:, i_median, :))))
  l_median25 = l_median .and. \
    .not.(all(ismissing(obser(:, i_25, :)))) .and. \
    .not.(all(ismissing(obser(:, i_75, :))))
  l_median10 = l_median .and. \
    .not.(all(ismissing(obser(:, i_10, :)))) .and. \
    .not.(all(ismissing(obser(:, i_90, :))))
  l_median5 = l_median .and. \
    .not.(all(ismissing(obser(:, i_5, :)))) .and. \
    .not.(all(ismissing(obser(:, i_95, :))))
  l_mean = .not.(all(ismissing(obser(:, i_mean, :))))
  l_stddev = .not.(all(ismissing(obser(:, i_stddev, :))))
  l_both = l_median.and.l_mean

  ; If both 10-90 and 5-95 are available, priority to 10-90
  if (l_median10.and.l_median5) then
    l_median5 = False
  end if

  ; Check if this is a summary plot
  l_summary = False
  if (ncases.gt.1) then
    l_summary = True
    l_stddev = False
    ; Only one stat and priority to median
    l_both = False
    if (l_median) then
      l_mean = False
    end if
  end if

  ; Set title for median only
  if (l_median.and..not.l_both) then
    res@gsnCenterString = "Median (lines), 25-75% percentiles (shades)"
    if (l_median10) then
      res@gsnCenterString = "Median (lines), " + \
        "25-75% and 10-90% percentiles (shades)"
    end if
    if (l_median5) then
      res@gsnCenterString = "Median (lines), " + \
        "25-75% and 5-95% percentiles (shades)"
    end if
    if (l_summary) then
      res@gsnCenterString = \
        str_sub_str(res@gsnCenterString, "(lines)", "(lines and crosses)")
    end if
  end if

  ; Set titles for mean only
  if (l_mean.and..not.l_both) then
    if (l_stddev) then
      res@gsnCenterString = "            " + \
        "Mean (dots and long-dashed lines) and " + \
        "~C~ standard deviation (whiskers and short-dashed lines)"
    else
      res@gsnCenterString = "Mean (dots and long-dashed lines) and " + \
        "standard deviation (short-dashed)"
    end if
    if (l_summary) then
      res@gsnCenterString = "Mean (crosses and long_dashed line) " + \
        "standard deviation (short-dashed)"
    end if
  end if

  ; Set titles for both median and mean
  if (l_both) then
    res@gsnCenterString =  "           " + \
      "Median and 25-75% percentiles " + \
      "(solid lines and shades), ~C~ mean and std. dev. "
    if (l_median10) then
      res@gsnCenterString =  "   " + \
        "Median, 25-75% and 10-90% percentiles " + \
        "(solid lines and shades), ~C~ mean and std. dev. "
    end if
    if (l_median5) then
      res@gsnCenterString =  "    " + \
        "Median, 25-75% and 5-95% percentiles " + \
        "(solid lines and shades), ~C~ mean and std."
    end if
    if (l_stddev) then
      res@gsnCenterString = res@gsnCenterString + \
        "(long- and short-dashed lines, dots and whiskers)"
    else
      res@gsnCenterString = res@gsnCenterString + \
        "(long- and short-dashed lines)"
    end if
  end if

  ; Draw a blank plot function, lines are dots are added below for
  ; better controlling the draw order
  plot = gsn_csm_blank_plot(wks, res)

  ; Draw percentiles as shaded areas
  if (l_median) then
    shade_flag_m = (/l_median5, l_median10, l_median/)
    shade_flag_o = (/l_median5, l_median10, l_median25/)
    shade_idx1 = (/i_5, i_10, i_25/)
    shade_idx2 = (/i_95, i_90, i_75/)
    shade_opac = (/0.1, 0.1, 0.15/)
    shade_name = (/"5_95", "10_90", "25_75"/)
    do shID = 0, dimsizes(shade_flag_m) - 1

      ; Model
      if (shade_flag_m(shID)) then
        resP@gsFillOpacityF = shade_opac(shID)
        g_obj = new(nmodels, graphic)
        do ii = 0, nmodels - 1
          resP@gsFillColor = colors(ii)
          x1 = model(ii, shade_idx1(shID), :)
          x2 = model(ii, shade_idx2(shID), :)
          y1 = model&plev
          y2 = y1
          xx = array_append_record( \
            array_append_record(x1, x2(::-1), 0), x1(0), 0)
          yy = array_append_record( \
            array_append_record(y1, y2(::-1), 0), y1(0), 0)
          q = ind(.not.ismissing(xx))
          if (dimsizes(q).ne.dimsizes(xx)) then
            log_debug("Some missing values in model percentiles")
          end if
          if (.not.all(ismissing(q))) then
            g_obj(ii) = gsn_add_polygon(wks, plot, xx(q), yy(q), resP)
          else
            log_info("All values in model percentiles are missing (skipping)")
          end if
          delete([/x1, x2, y1, y2, xx, yy, q/])
        end do
        str = "model_" + shade_name(shID)
        plot@$str$ = g_obj
        delete(g_obj)
        delete(str)
      end if

      ; Observations
      if (shade_flag_o(shID).and..not.l_summary) then
        resP@gsFillColor = "black"
        g_obj = new(ncases, graphic)
        do ii = 0, ncases - 1
          x1 = obser(ii, shade_idx1(shID), :)
          x2 = obser(ii, shade_idx2(shID), :)
          y1 = obser&plev
          y2 = y1
          xx = array_append_record( \
            array_append_record(x1, x2(::-1), 0), x1(0), 0)
          yy = array_append_record( \
            array_append_record(y1, y2(::-1), 0), y1(0), 0)
          q = ind(.not.ismissing(xx))
          if (dimsizes(q).ne.dimsizes(xx)) then
            log_debug("Some missing values in obs percentiles")
          end if
          if (.not.all(ismissing(q))) then
            g_obj(ii) = gsn_add_polygon(wks, plot, xx(q), yy(q), resP)
          else
            log_info("All values in obs percentiles are missing (skipping)")
          end if
          delete([/x1, x2, y1, y2, xx, yy, q/])
        end do
        str = "obser_" + shade_name(shID)
        plot@$str$ = g_obj
        delete(g_obj)
        delete(str)
      end if
    end do
  end if

  ; Draw median as solid lines
  if (l_median) then

    ; Model
    g_obj = new(nmodels, graphic)
    do ii = 0, nmodels - 1
      resL@gsLineColor = colors(ii)
      resL@gsLineDashPattern = dash_median
      xx = model(ii, i_median, :)
      yy = model&plev
      q = ind(.not.ismissing(xx))
      if (dimsizes(q).ne.dimsizes(xx)) then
        log_debug("Some missing values in model median")
      end if
      if (.not.all(ismissing(q))) then
        g_obj(ii) = gsn_add_polyline(wks, plot, xx(q), yy(q), resL)
      else
        log_info("All values in model median are missing (skipping)")
      end if
      delete([/xx, yy, q/])
    end do
    plot@model_median = g_obj
    delete(g_obj)

    ; Observations
    if (l_summary) then
      resM@gsMarkerColor = "black"
      resM@gsMarkerIndex = 5
      resM@gsMarkerSizeF = 0.015
    else
      resL@gsLineColor = "black"
    end if
    g_obj = new(ncases, graphic)
    do ii = 0, ncases - 1
      xx = obser(ii, i_median, :)
      yy = obser&plev
      q = ind(.not.ismissing(xx))
      if (dimsizes(q).ne.dimsizes(xx)) then
        log_debug("Some missing values in obs. median")
      end if
      if (.not.all(ismissing(q))) then
        if (l_summary) then
          g_obj(ii) = gsn_add_polymarker(wks, plot, xx(q), yy(q), resM)
        else
          g_obj(ii) = gsn_add_polyline(wks, plot, xx(q), yy(q), resL)
        end if
      else
        log_info("All values in obs median are missing (skipping)")
      end if
      delete([/xx, yy, q/])
    end do
    plot@obser_median = g_obj
    delete(g_obj)

  end if

  ; Draw mean and standard deviation as dashed/dotted (model) and dots/bars
  ; (observations)
  if (l_mean) then

    ; Model mean and std. dev.
    g_obj = new(nmodels, graphic)
    do ii = 0, nmodels - 1
      resL@gsLineColor = colors(ii)
      resL@gsLineDashPattern = dash_mean
      xx = model(ii, i_mean, :)
      yy = model&plev
      q = ind(.not.ismissing(xx))
      if (dimsizes(q).ne.dimsizes(xx)) then
        log_debug("Some missing values in model mean")
      end if
      if (.not.all(ismissing(q))) then
        g_obj(ii) = gsn_add_polyline(wks, plot, xx(q), yy(q), resL)
      else
        log_info("All values in model mean are missing (skipping)")
      end if
      delete([/xx, yy, q/])
    end do
    plot@model_mean = g_obj
    delete(g_obj)

    g_obj = new(nmodels, graphic)
    do ii = 0, nmodels - 1
      resL@gsLineColor = colors(ii)
      resL@gsLineDashPattern = dash_stddev
      xx = model(ii, i_mean, :) + model(ii, i_stddev, :)
      yy = model&plev
      q = ind(.not.ismissing(xx))
      if (dimsizes(q).ne.dimsizes(xx)) then
        log_debug("Some missing values in model std. dev.")
      end if
      if (.not.all(ismissing(q))) then
        g_obj(ii) = gsn_add_polyline(wks, plot, xx(q), yy(q), resL)
      else
        log_info("All values in model std. dev. are missing (skipping)")
      end if
      delete([/xx, yy, q/])
    end do
    plot@model_stddev = g_obj
    delete(g_obj)

    ; Observations mean (dots or lines, depending on standard deviation)
    if (l_stddev) then
      resM@gsMarkerColor = "black"
      resM@gsMarkerIndex = marker_mean
      resM@gsMarkerSizeF = 0.01
    else
      resL@gsLineColor = "black"
      resL@gsLineDashPattern = dash_mean
    end if
    g_obj = new(ncases, graphic)
    do ii = 0, ncases - 1
      xx = obser(ii, i_mean, :)
      yy = obser&plev
      q = ind(.not.ismissing(xx))
      if (dimsizes(q).ne.dimsizes(xx)) then
        log_debug("Some missing values in obs. mean")
      end if
      if (.not.all(ismissing(q))) then
        if (l_stddev) then
          g_obj(ii) = gsn_add_polymarker(wks, plot, xx(q), yy(q), resM)
        else
          g_obj(ii) = gsn_add_polyline(wks, plot, xx(q), yy(q), resL)
        end if
      else
        log_info("All values in obs. mean are missing (skipping)")
      end if
      delete([/xx, yy, q/])
    end do
    plot@obser_mean = g_obj
    delete(g_obj)

    ; Observations std. dev.
    if (l_stddev) then
      resL@gsLineColor = "black"
      resL@gsLineDashPattern = 0
      resL@gsLineThicknessF = 1
      do ii = 0, ncases - 1
        g_obj = new(dimsizes(obser&plev), graphic)
        do kk = 0, dimsizes(obser&plev) - 1
          xx = (/obser(ii, i_mean, kk), \
                obser(ii, i_mean, kk) + obser(ii, i_stddev, kk)/)
          yy = (/obser&plev(kk), obser&plev(kk)/)
          if (.not.ismissing(xx(0)) .and. .not.ismissing(xx(1))) then
            g_obj(kk) = gsn_add_polyline(wks, plot, xx, yy, resL)
          end if
        end do
        str = "obser" + ii + "_stddev"
        plot@$str$ = g_obj
        delete(str)
        delete([/xx, yy/])
        delete(g_obj)
      end do
    end if
  end if
  draw(plot)

  ; Append legend
  xpos = 0.21
  resT@txJust = "CenterLeft"
  ypos = 0.78
  resT@txFontHeightF = 0.020
  resT@txFontColor   = "black"
  gsn_text_ndc(wks, diag_script_info@campaign, xpos, ypos, resT)
  do mID = 0, dimsizes(annots) - 1
    ypos = ypos - 0.03
    resT@txFontColor   = colors(mID)
    gsn_text_ndc(wks, annots(mID), xpos, ypos, resT)
  end do

  frame(wks)
  leave_msg(scriptname, funcname)
  return(plot)

end

; #############################################################################
undef("aerosol_sizedist")
function aerosol_sizedist(wks_in[1],
                          source_m,
                          source_o,
                          varname[1])
;
; Arguments
;    wks_in: workstations (graphic object or default will be used).
;    source_m: model data to be plotted or a NetCDF filename with data.
;    source_o: observations data to be plotted or a NetCDF filename with data.
;    varname: variable name in the file.
;
; Source_m prototype
;    source[*][*][*]
;    source!0 = model
;    source!1 = statistic
;    source!2 = plev
;
; Source_o prototype
;    source[*][*][*]
;    source!0 = case
;    source!1 = statistic
;    source!2 = plev
;
; Return value
;    A graphic variable.
;
; Description
;    Creates a size distribution plot (data vs. diameter).
;    Plots median, mean or both depending on availability of obervations.
;
; Caveats
;
; Modification history
;    20161013-righi_mattia: written based on aerosol_profile.
;
local funcname, scriptname, tmp, colors, annots, res, res_new
begin

  funcname = "aerosol_sizedist"
  scriptname = "diag_scripts/shared/plot/xy_line.ncl"
  enter_msg(scriptname, funcname)

  ; Get model data, either directly or via netCDF file
  if (typeof(source_m) .eq. "string") then
    model = ncdf_read(source_m, varname)
  else
    model = source_m
    copy_VarMeta(source_m, model)
  end if

  ; Get observational data, either directly or via netCDF file
  if (typeof(source_o) .eq. "string") then
    obser = ncdf_read(source_o, varname)
  else
    obser = source_o
    copy_VarMeta(source_o, obser)
  end if

  ; Retrieve basic metadata from data
  defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
  if (any(varname .eq. defaults)) then
    var = att2var(model, "var")
  else
    var = varname
  end if

  ; Check if a valid wks has been provided, otherwise invoke default
  wks = get_wks(wks_in, diag_script, var)

  ; Define line colors (first color is used for the observations)
  ; FIX-ME: project_style function has changed
  ; colors = project_style(diag_script_info, "colors")
  ; annots = project_style(diag_script_info, "annots")
  gsn_define_colormap(wks, array_append_record((/"white", "black"/), \
                      colors, 0))

  ; Define markers and dashes for median and mean
  marker_mean = 16
  dash_median = 0
  dash_mean = 1
  dash_stddev = 2

  ; Resources
  res = True
  res@gsnDraw = False
  res@gsnFrame = False
  res@trXAxisType = "LogAxis"
  res@trYAxisType = "LogAxis"
  res@tmXBMinorPerMajor = 8
  res@tmXTMinorPerMajor = 8
  res@tmYLMinorPerMajor = 8
  res@tmYRMinorPerMajor = 8
  res@tiYAxisString = varname + " [" + format_units(model@units) + "]"
  res@tiXAxisString = "Diameter [" + format_units(model&diam@units) + "]"
  res@gsnCenterStringFontHeightF = 0.015

  ; Lines
  resL = True
  resL@gsLineThicknessF = 2.0
  resL@gsLineDashSegLenF = 0.25

  ; Markers
  resM = True
  resM@gsMarkerSizeF = 0.01

  ; Polygons
  resP = True

  ; Text
  resT = True
  resT@txFont = 21

  ; Markers
  resM = True

  ; Override defaults with "res_" attributes of "data"
  res_new = att2var(model, "res_")
  copy_VarMeta(res_new, res)

  ; Set number of models
  nmodels = dimsizes(annots)

  ; Set number of cases
  ncases = dimsizes(obser&case)

  ; Set the index for the various statistics
  i_N      = ind(model&statistic.eq."N")
  i_mean   = ind(model&statistic.eq."mean")
  i_stddev = ind(model&statistic.eq."stddev")
  i_min    = ind(model&statistic.eq."min")
  i_max    = ind(model&statistic.eq."max")
  i_5      = ind(model&statistic.eq."5")
  i_10     = ind(model&statistic.eq."10")
  i_25     = ind(model&statistic.eq."25")
  i_median = ind(model&statistic.eq."median")
  i_75     = ind(model&statistic.eq."75")
  i_90     = ind(model&statistic.eq."90")
  i_95     = ind(model&statistic.eq."95")

  ; Set plot ranges
  modelmin = min((/min(model(:, i_mean, :)), min(model(:, i_median, :))/))
  modelmax = max((/max(model(:, i_mean, :)), max(model(:, i_median, :))/))
  obsermin = min((/min(obser(:, i_mean, :)), min(obser(:, i_median, :))/))
  obsermax = max((/max(obser(:, i_mean, :)), max(obser(:, i_median, :))/))
  res@tmXBMode = "Explicit"
  res@tmYLMode = "Explicit"

  ; Set X-axis labels
  res@trXMinF = 0.001
  res@trXMaxF = 10.
  res@tmXBValues = set_log_ticks(res@trXMinF, res@trXMaxF, "major")
  res@tmXBLabels = res@tmXBValues
  res@tmXBMinorValues = set_log_ticks(res@trXMinF, res@trXMaxF, "minor")

  ; Set Y-axis labels
  res@trYMaxF = 2 * max((/modelmax, obsermax/))
  res@trYMinF = res@trYMaxF / 1.e5
  res@tmYLValues = set_log_ticks(res@trYMinF, res@trYMaxF, "major")
  res@tmYLLabels = res@tmYLValues
  res@tmYLMinorValues = set_log_ticks(res@trYMinF, res@trYMaxF, "minor")

  ; Check available statistics in the observations
  l_median = .not.(all(ismissing(obser(:, i_median, :))))
  l_median25 = l_median .and. \
    .not.(all(ismissing(obser(:, i_25, :)))) .and. \
    .not.(all(ismissing(obser(:, i_75, :))))
  l_median10 = l_median .and. \
    .not.(all(ismissing(obser(:, i_10, :)))) .and. \
    .not.(all(ismissing(obser(:, i_90, :))))
  l_median5 = l_median .and. \
    .not.(all(ismissing(obser(:, i_5, :)))) .and. \
    .not.(all(ismissing(obser(:, i_95, :))))
  l_mean = .not.(all(ismissing(obser(:, i_mean, :))))
  l_stddev = .not.(all(ismissing(obser(:, i_stddev, :))))
  l_both = l_median.and.l_mean

  ; If both 10-90 and 5-95 are available, priority to 10-90
  if (l_median10.and.l_median5) then
    l_median5 = False
  end if

  ; Set title for median only
  if (l_median.and..not.l_both) then
    res@gsnCenterString = "Median (lines), 25-75% percentiles (shades)"
    if (l_median10) then
      res@gsnCenterString = "Median (lines), " + \
        "25-75% and 10-90% percentiles (shades)"
    end if
    if (l_median5) then
      res@gsnCenterString = "Median (lines), " + \
        "25-75% and 5-95% percentiles (shades)"
    end if
  end if

  ; Set titles for mean only
  if (l_mean.and..not.l_both) then
    if (l_stddev) then
      res@gsnCenterString = "            " + \
        "Mean (dots and long-dashed lines) and " + \
        "~C~ standard deviation (whiskers and short-dashed lines)"
    else
      res@gsnCenterString = "Mean (long-dashed) and " + \
        "standard deviation (short-dashed)"
    end if
  end if

  ; Set titles for both median and mean
  if (l_both) then
    res@gsnCenterString =  "           " + \
      "Median and 25-75% percentiles " + \
      "(solid lines and shades), ~C~ mean and std. dev. "
    if (l_median10) then
      res@gsnCenterString =  "   " + \
        "Median, 25-75% and 10-90% percentiles " + \
        "(solid lines and shades), ~C~ mean and std. dev. "
    end if
    if (l_median5) then
      res@gsnCenterString =  "    " + \
        "Median, 25-75% and 5-95% percentiles " + \
        "(solid lines and shades), ~C~ mean and std."
    end if
    if (l_stddev) then
      res@gsnCenterString = res@gsnCenterString + \
        "(long- and short-dashed lines, dots and whiskers)"
    else
      res@gsnCenterString = res@gsnCenterString + \
        "(long- and short-dashed lines)"
    end if
  end if

  ; Draw a blank plot function, lines are dots are added below for
  ; better controlling the draw order
  plot = gsn_csm_blank_plot(wks, res)

  ; Draw percentiles as shaded areas
  if (l_median) then
    shade_flag_m = (/l_median5, l_median10, l_median/)
    shade_flag_o = (/l_median5, l_median10, l_median25/)
    shade_idx1 = (/i_5, i_10, i_25/)
    shade_idx2 = (/i_95, i_90, i_75/)
    shade_opac = (/0.1, 0.1, 0.15/)
    shade_name = (/"5_95", "10_90", "25_75"/)
    do shID = 0, dimsizes(shade_flag_m) - 1

      ; Model
      if (shade_flag_m(shID)) then
        resP@gsFillOpacityF = shade_opac(shID)
        g_obj = new(nmodels, graphic)
        do ii = 0, nmodels - 1
          resP@gsFillColor = colors(ii)
          y1 = model(ii, shade_idx1(shID), :)
          y2 = model(ii, shade_idx2(shID), :)
          x1 = model&diam
          x2 = x1
          xx = array_append_record( \
            array_append_record(x1, x2(::-1), 0), x1(0), 0)
          yy = array_append_record( \
            array_append_record(y1, y2(::-1), 0), y1(0), 0)
          q = ind(.not.ismissing(yy))
          if (dimsizes(q).ne.dimsizes(yy)) then
            log_debug("Some missing values in model percentiles")
          end if
          if (.not.all(ismissing(q))) then
            g_obj(ii) = gsn_add_polygon(wks, plot, xx(q), yy(q), resP)
          else
            log_info("All values in model percentiles are missing (skipping)")
          end if
          delete([/x1, x2, y1, y2, xx, yy, q/])
        end do
        str = "model_" + shade_name(shID)
        plot@$str$ = g_obj
        delete(g_obj)
        delete(str)
      end if

      ; Observations
      if (shade_flag_o(shID)) then
        resP@gsFillColor = "black"
        g_obj = new(ncases, graphic)
        do ii = 0, ncases - 1
          y1 = obser(ii, shade_idx1(shID), :)
          y2 = obser(ii, shade_idx2(shID), :)
          x1 = obser&diam
          x2 = x1
          xx = array_append_record( \
            array_append_record(x1, x2(::-1), 0), x1(0), 0)
          yy = array_append_record( \
            array_append_record(y1, y2(::-1), 0), y1(0), 0)
          q = ind(.not.ismissing(yy))
          if (dimsizes(q).ne.dimsizes(yy)) then
            log_debug("Some missing values in obs percentiles")
          end if
          if (.not.all(ismissing(q))) then
            g_obj(ii) = gsn_add_polygon(wks, plot, xx(q), yy(q), resP)
          else
            log_info("All values in obs percentiles are missing (skipping)")
          end if
          delete([/x1, x2, y1, y2, xx, yy, q/])
          str = "obser_" + shade_name(shID)
          plot@$str$ = g_obj
          delete(g_obj)
          delete(str)
        end do
      end if

    end do
  end if

  ; Draw median as solid lines
  if (l_median) then

    ; Model
    g_obj = new(nmodels, graphic)
    do ii = 0, nmodels - 1
      resL@gsLineColor = colors(ii)
      resL@gsLineDashPattern = dash_median
      xx = model&diam
      yy = model(ii, i_median, :)
      q = ind(.not.ismissing(yy))
      if (dimsizes(q).ne.dimsizes(yy)) then
        log_debug("Some missing values in model median")
      end if
      if (.not.all(ismissing(q))) then
        g_obj(ii) = gsn_add_polyline(wks, plot, xx(q), yy(q), resL)
      else
        log_info("All values in model median are missing (skipping)")
      end if
      delete([/xx, yy, q/])
    end do
    plot@model_median = g_obj
    delete(g_obj)

    ; Observations
    resL@gsLineColor = "black"
    g_obj = new(ncases, graphic)
    do ii = 0, ncases - 1
      xx = obser&diam
      yy = obser(ii, i_median, :)
      q = ind(.not.ismissing(yy))
      if (dimsizes(q).ne.dimsizes(yy)) then
        log_debug("Some missing values in obs. median")
      end if
      if (.not.all(ismissing(q))) then
        g_obj(ii) = gsn_add_polyline(wks, plot, xx(q), yy(q), resL)
      else
        log_info("All values in obs median are missing (skipping)")
      end if
      delete([/xx, yy, q/])
    end do
    plot@obser_median = g_obj
    delete(g_obj)

  end if

  ; Draw mean and standard deviation as dashed/dotted (model) and dots/bars
  ; (observations)
  if (l_mean) then

    ; Model
    g_obj = new(nmodels, graphic)
    do ii = 0, nmodels - 1
      resL@gsLineColor = colors(ii)
      resL@gsLineDashPattern = dash_mean
      xx = model&diam
      yy = model(ii, i_mean, :)
      q = ind(.not.ismissing(yy))
      if (dimsizes(q).ne.dimsizes(yy)) then
        log_debug("Some missing values in model mean")
      end if
      if (.not.all(ismissing(q))) then
        g_obj(ii) = gsn_add_polyline(wks, plot, xx(q), yy(q), resL)
      else
        log_info("All values in model mean are missing (skipping)")
      end if
      delete([/xx, yy, q/])
    end do
    plot@model_mean = g_obj
    delete(g_obj)

    g_obj = new(nmodels, graphic)
    do ii = 0, nmodels - 1
      resL@gsLineColor = colors(ii)
      resL@gsLineDashPattern = dash_stddev
      xx = model&diam
      yy = model(ii, i_mean, :) + model(ii, i_stddev, :)
      q = ind(.not.ismissing(yy))
      if (dimsizes(q).ne.dimsizes(yy)) then
        log_debug("Some missing values in model std. dev.")
      end if
      if (.not.all(ismissing(q))) then
        g_obj(ii) = gsn_add_polyline(wks, plot, xx(q), yy(q), resL)
      else
        log_info("All values in model std. dev. are missing (skipping)")
      end if
      delete([/xx, yy, q/])
    end do
    plot@model_stddev = g_obj
    delete(g_obj)

        ; Observations (dots or lines, depending on standard deviation)
    if (l_stddev) then
      resM@gsMarkerColor = "black"
      resM@gsMarkerIndex = marker_mean
    else
      resL@gsLineColor = "black"
      resL@gsLineDashPattern = dash_mean
    end if
    g_obj = new(ncases, graphic)
    do ii = 0, ncases - 1
      xx = obser&diam
      yy = obser(ii, i_mean, :)
      q = ind(.not.ismissing(yy))
      if (dimsizes(q).ne.dimsizes(yy)) then
        log_debug("Some missing values in obs. mean")
      end if
      if (.not.all(ismissing(q))) then
        if (l_stddev) then
          g_obj(ii) = gsn_add_polymarker(wks, plot, xx(q), yy(q), resM)
        else
          g_obj(ii) = gsn_add_polyline(wks, plot, xx(q), yy(q), resL)
        end if
      else
        log_info("All values in obs. mean are missing (skipping)")
      end if
      delete([/xx, yy, q/])
    end do
    plot@obser_mean = g_obj
    delete(g_obj)

    ; Observations
    if (l_stddev) then
      resL@gsLineColor = "black"
      resL@gsLineDashPattern = 0
      resL@gsLineThicknessF = 1
      do ii = 0, ncases - 1
        g_obj = new(dimsizes(obser&diam), graphic)
        do kk = 0, dimsizes(obser&diam) - 1
          xx = (/obser&diam(kk), obser&diam(kk)/)
          yy = (/obser(ii, i_mean, kk), \
                obser(ii, i_mean, kk) + obser(ii, i_stddev, kk)/)
          if (.not.ismissing(xx(0)) .and. .not.ismissing(xx(1))) then
            g_obj(kk) = gsn_add_polyline(wks, plot, xx, yy, resL)
          end if
        end do
        str = "obser" + ii + "_stddev"
        plot@$str$ = g_obj
        delete(str)
        delete([/xx, yy/])
        delete(g_obj)
      end do
    end if
  end if
  draw(plot)

  ; Append legend
  xpos = 0.79
  resT@txJust = "CenterRight"
  ypos = 0.78
  resT@txFontHeightF = 0.020
  resT@txFontColor   = "black"
  gsn_text_ndc(wks, diag_script_info@campaign, xpos, ypos, resT)
  do mID = 0, dimsizes(annots) - 1
    ypos = ypos - 0.03
    resT@txFontColor   = colors(mID)
    gsn_text_ndc(wks, annots(mID), xpos, ypos, resT)
  end do

  frame(wks)

  leave_msg(scriptname, funcname)
  return(plot)

end

; #############################################################################
undef("xy_line")
procedure xy_line(wks[1],
                  source,
                  source_x,
                  source_stddev,
                  res_in:logical,
                  items:list)
;
; Arguments:
;    wks: workstation, must be passed - no default used yet!
;    source: data to be plotted (no netCDF input possible yet)
;    source_x: x-axis of array to be plotted (e.g. source&time, ... )
;    source_stddev: standard deviation of input, needed for multi_model_mean
;    res_in: diag_script-specific resources passed from diag_script
;    items: list of input_file_info metadata
;
; Source prototype
;
; Description
;    Defines default ressources, which are overridden by argument res.
;    Creates an xy-plot, according to wks & res.
;    Adds multi model mean and standard deviation if
;    diag_script_info@multi_model_mean is set to True.
;
; Caveats
;
; Modification history
;    20150511-senftleben_daniel: modified legend.
;    20140109-senftleben_daniel: written.
;
local funcname, scriptname, res, res_in, res_stddev, source, source_x, \
  source_stddev, wks, wks_in, colors, colors_mm, dashes, dashes_mm, thicks, \
  thicks_mm, annots, annots_mm, avgstd, avgstd_mm, temp, plot, \
  shading_plot, mm, lgres, nitems, lbid, amres, annoid, labels, lg_outfile, \
  psres, vpx, vph, vpy, vpw, bpres, tmborder, items_mmm, mmm_var
begin

  funcname = "xy_line"
  scriptname = "diag_scripts/shared/plot/xy_line.ncl"
  enter_msg(scriptname, funcname)

  ; Select colors and other plotting attributes
  ; (see ./diag_scripts/shared/plot/style.ncl)
  colors = project_style(items, diag_script_info, "colors")
  dashes = project_style(items, diag_script_info, "dashes")
  thicks = project_style(items, diag_script_info, "thicks")
  annots = project_style(items, diag_script_info, "annots")
  avgstd = project_style(items, diag_script_info, "avgstd")

  ; Select colors and other plotting attributes for multi-model mean
  if (diag_script_info@multi_model_mean) then
    items_mmm = NewList("fifo")
    mmm_var = True
    mmm_var@dataset = "MultiModelMean"
    ListPush(items_mmm, mmm_var)
    colors_mmm = project_style(items_mmm, diag_script_info, "colors")
    dashes_mmm = project_style(items_mmm, diag_script_info, "dashes")
    thicks_mmm = project_style(items_mmm, diag_script_info, "thicks")
    annots_mmm = project_style(items_mmm, diag_script_info, "annots")
    avgstd_mmm = project_style(items_mmm, diag_script_info, "avgstd")
  end if

  ; Set/copy resources
  res                   = True
  res@xyLineColors      = colors  ; change line color
  res@xyLineThicknesses = thicks  ; make 2nd lines thicker
  res@xyMonoDashPattern = False
  res@xyDashPatterns    = dashes
  res@gsnDraw           = False   ; don't draw yet
  res@gsnFrame          = False   ; don't advance frame yet
  res@tiMainFontHeightF = .025    ; font height
  res@gsnMaximize       = True
  res@vpHeightF = 0.4
  res@vpWidthF = 0.8
  res@trYMinF  = min(source) - 0.05 * (max(source) - min(source))
  res@trYMaxF  = max(source) + 0.05 * (max(source) - min(source))
  copy_VarMeta(res_in, res)  ; copy passed resources

  ; Plot
  plot = gsn_csm_xy(wks, source_x, source, res)

  ; Add multi model mean and stddev
  if (diag_script_info@multi_model_mean) then
    ; Stddev
    res_stddev = True
    copy_VarMeta(res, res_stddev)
    res_stddev@gsnXYFillColors = "LightGrey"
    delete(res_stddev@xyLineColors)
    res_stddev@xyLineColor = -1  ; Make lines transparent
    shading_plot = gsn_csm_xy(wks, source_x, source_stddev(2:3, :), \
                              res_stddev)
    overlay(plot, shading_plot)
    ; MMM
    delete([/res@xyLineThicknesses, res@xyLineColors, res@xyDashPatterns/])
    res@xyLineThicknesses = thicks_mmm
    res@xyLineColors = colors_mmm
    res@xyDashPatterns = dashes_mmm
    mmm = gsn_csm_xy(wks, source_x, source_stddev(0, :), res)
    overlay(plot, mmm)
  end if

  ; Resources for a customized legend.
  lgres                     = True
  lgres@lgMonoLineThickness = False
  lgres@lgLabelFontHeightF  = .08      ; legend label font thickness
  lgres@vpWidthF            = 0.15     ; width of legend (NDC)
  lgres@vpHeightF           = 0.4      ; height of legend (NDC)
  lgres@lgPerimOn           = True
  lgres@lgPerimColor        = 0
  lgres@lgMonoDashIndex     = False
  lgres@lgBoxBackground     = 0
  lgres@lgPerimFill         = 0
  lgres@lgPerimFillColor    = 0
  if (.not.isatt(diag_script_info, "EMs_in_lg")) then
    diag_script_info@EMs_in_lg = True  ; Set default
  end if
  if (.not.diag_script_info@EMs_in_lg) then
    datasets = metadata_att_as_array(items, "dataset")
    idcs_modelnames = UNIQ(datasets)
    colors := colors(idcs_modelnames)
    dashes := dashes(idcs_modelnames)
    thicks := thicks(idcs_modelnames)
    annots := datasets(idcs_modelnames)
  end if
  if (diag_script_info@multi_model_mean) then
    lgres@lgLineColors      = array_append_record(colors, colors_mmm, 0)
    lgres@lgDashIndexes     = array_append_record(dashes, dashes_mmm, 0)
    lgres@lgLineThicknesses = \
      array_append_record(thicks, thicks_mmm, 0) + 0.5
    labels = array_append_record(annots, annots_mmm, 0)
  else
    lgres@lgLineColors      = colors
    lgres@lgDashIndexes     = dashes
    lgres@lgLineThicknesses = thicks + 0.5
    labels = annots
  end if
  nitems = dimsizes(labels)
  lgres@lgItemOrder = ispan(nitems - 1, 0, 1)

  ; Create legend
  lbid = gsn_create_legend(wks, nitems, labels, lgres)

  ; Add legend at the correct position
  ; Point (0, 0) is the dead center of the plot. Point (0, .5) is center,
  ; flush bottom. Point (0.5, 0.5) is flush bottom, flush right.
  amres                  = True
  amres@amJust           = "TopRight"  ; reference corner of box
  amres@amParallelPosF   =  0.5        ; Move legend to +right, -left
  amres@amOrthogonalPosF = -0.5        ; +down, -up

  ; Draw and frame
  psres = True
  psres@gsnDraw = False
  psres@gsnFrame = False
  maximize_output(wks, psres)

  ; Get plot size
  getvalues plot
    "vpXF"      : vpx
    "vpYF"      : vpy
    "vpWidthF"  : vpw
    "vpHeightF" : vph
  end getvalues
  bres = True
  bres@gsnDraw      = False
  bres@gsnFrame     = False
  bres@tmXBOn       = False
  bres@tmYLOn       = False
  bres@tmXTOn       = False
  bres@tmYROn       = False
  bres@tmXBBorderOn = False
  bres@tmXTBorderOn = False
  bres@tmYLBorderOn = False
  bres@tmYRBorderOn = False

  ; Create a blank plot with the same size as plot, attach legend
  bres@vpXF      = vpx
  bres@vpYF      = vpy
  bres@vpWidthF  = vpw
  bres@vpHeightF = vph
  blank_plot = gsn_csm_blank_plot(wks, bres)

  ; Add legend to plot
  if (isatt(diag_script_info, "xy_line_legend")) then
    if (diag_script_info@xy_line_legend .ne. False) then
      annoid = gsn_add_annotation(blank_plot, lbid, amres)
    end if
  else
    annoid = gsn_add_annotation(blank_plot, lbid, amres)
  end if

  ; Put legend into an extra file
  if (isatt(source, "legend_outside")) then
    if (source@legend_outside) then
      lg_outfile = wks@legendfile
      styles = True
      styles@colors = colors
      styles@dashes = dashes
      styles@thicks = thicks
      create_legend_lines(annots, styles, lg_outfile, "lines")
    end if
  end if

  ; Create another blank plot to make sure plot border thickness is even
  bres@tmXBBorderOn = True
  bres@tmXTBorderOn = True
  bres@tmYLBorderOn = True
  bres@tmYRBorderOn = True
  bres@tmBorderThicknessF = 3
  blank_plot2 = gsn_csm_blank_plot(wks, bres)

  ; Draw first plot with the actual values (+ grid lines if
  ; tmXMajorGrid/gsnYRefLine are set)
  draw(plot)

  ; Draw second plot with legend on top of previous plot. This
  ; is, as far as I know, the only way to draw the legend on top
  ; of the grid lines
  draw(blank_plot)

  ; Redraw plot borders since the legend may (partially) cover some
  ; of the borders drawn in the first 'plot'
  draw(blank_plot2)

  frame(wks)

  ; outfile name
  if (isatt(wks, "fullname")) then
    outfile = wks@fullname
  else
    outfile = wks@name
  end if
  log_info(" Wrote " + outfile)
  leave_msg(scriptname, funcname)

end

; #############################################################################
undef("xy_line_anom")
procedure xy_line_anom(wks[1],
                       source_mean,
                       source,
                       source_x,
                       source_stddev,
                       ref_start,
                       ref_end,
                       res_in: logical,
                       res0_in: logical,
                       items: list)
;
; Arguments:
;    wks:  workstation, must be passed - no default used yet!
;    source_mean:   source_mean
;    source:        data to be plotted (no netCDF input possible yet)
;    source_x:      x-axis of array to be plotted (e.g. source&time, ... )
;    source_stddev: standard deviation of input, needed multi_model_mean
;    ref_start:     start year of the reference dataset
;    ref_end:       end year of the reference dataset
;    res_in:        diag_script-specific resources passed from diag_script
;    res0_in:       res0_in
;    items: list of input_file_info metadata
;
; Source prototype
;
; Description
;    Defines default ressources, which are overridden by argument res.
;    Creates an xy-plot, according to wks & res.
;    Adds multi model mean and standard deviation if
;    diag_script_info@multi_model_mean is set to True.
;
; Caveats
;
; Modification history
;    20160822-bock_lisa: written
;

local funcname, scriptname, verbosity, res, res_in, res_stddev, source, \
    source_x, source_stddev, wks, wks_in, colors, colors_mm, dashes, \
    dashes_mm, thicks, thicks_mm, annots, annots_mm, avgstd, avgstd_mm, temp, \
    plot, shading_plot, mm, lgres, nitems, lbid, amres, annoid, labels, \
    psres, vpx, vph, vpy, vpw, bpres, tmborder

begin

  funcname = "xy_line_anom"
  scriptname = "plot_scripts/ncl/xy_line.ncl"
  enter_msg(scriptname, funcname)

  ; Select colors and other plotting attributes
  ; (see ./diag_scripts/lib/ncl/style.ncl)
  colors = project_style(items, diag_script_info, "colors")
  dashes = project_style(items, diag_script_info, "dashes")
  thicks = project_style(items, diag_script_info, "thicks")
  annots = project_style(items, diag_script_info, "annots")
  avgstd = project_style(items, diag_script_info, "avgstd")

  ; individual case for HadCRUT4 observations
  ; FIX-ME: mean value comes from climatology file (absolute.nc).
  ; There are no missing values as in the anomaly data.
  datasetnames = metadata_att_as_array(items, "dataset")
  if (any(datasetnames.eq."HadCRUT4-clim")) then
    ind_wo_clim = ind(datasetnames .ne. "HadCRUT4-clim")
    tmp1 = colors(ind_wo_clim)
    delete(colors)
    colors = tmp1
    tmp2 = dashes(ind_wo_clim)
    delete(dashes)
    dashes = tmp2
    tmp3 = thicks(ind_wo_clim)
    delete(thicks)
    thicks = tmp3
    tmp4 = annots(ind_wo_clim)
    delete(annots)
    annots = tmp4
    tmp5 = avgstd(ind_wo_clim)
    delete(avgstd)
    avgstd = tmp5
    delete([/tmp1, tmp2, tmp3, tmp4, tmp5/])
  end if

  ; ************************************************
  ; plotting parameters
  ; ************************************************

  plot = new(1, graphic)

  res0 = True
  res0@gsnDraw   = False
  res0@gsnFrame  = False
  res0@vpHeightF = 0.5
  res0@vpWidthF  = 0.05
  res0@trYMinF = \
    min(source_mean) - 0.05 * (max(source_mean) - min(source_mean))
  res0@trYMaxF = \
    max(source_mean) + 0.05 * (max(source_mean) - min(source_mean))
  res0@trXMinF = 0.
  res0@trXMaxF = 2.
  res0@tmXTOn = False
  res0@tmXBOn = False
  res0@tmYLLabelsOn = False
  res0@tmYRLabelsOn = True
  res0@tmYRLabelFontHeightF = 0.016
  res0@tmYLLabelFontHeightF = 0.016
  res0@tiXAxisFontHeightF = 0.016
  res0@tiYAxisFontHeightF = 0.016
  res0@tiYAxisSide   = "Right"
  res0@tiYAxisAngleF = 90.
  res0@tiYAxisOn = True
  res0@pmLegendDisplayMode = "Never"
  res0@tmXBMajorOutwardLengthF = 0.006
  res0@tmYLMajorOutwardLengthF = 0.006
  res0@tmXBMinorOutwardLengthF = 0.003
  res0@tmYLMinorOutwardLengthF = 0.003
  res0@tmXBMajorLengthF = 0.006
  res0@tmYLMajorLengthF = 0.006
  res0@tmXBMinorLengthF = 0.003
  res0@tmYLMinorLengthF = 0.003

  res0@xyDashPatterns    = dashes
  res0@xyLineThicknesses = thicks
  res0@xyLineColors      = colors

  copy_VarMeta(res0_in, res0)  ; copy passed resources

  res           = True
  res@gsnDraw   = False
  res@gsnFrame  = False

  res@vpXF      = 0.05
  res@vpYF      = 0.7
  res@vpHeightF = 0.4
  res@vpWidthF  = 0.7
  res@pmLegendDisplayMode = "Never"
  res@tmYRLabelFontHeightF = 0.016
  res@tmYLLabelFontHeightF = 0.016
  res@tiXAxisFontHeightF = 0.016
  res@tiYAxisFontHeightF = 0.016
  res@tmXBMajorOutwardLengthF = 0.006
  res@tmYLMajorOutwardLengthF = 0.006
  res@tmXBMinorOutwardLengthF = 0.003
  res@tmYLMinorOutwardLengthF = 0.003
  res@tmXBMajorLengthF = 0.006
  res@tmYLMajorLengthF = 0.006
  res@tmXBMinorLengthF = 0.003
  res@tmYLMinorLengthF = 0.003

  res@trYMinF  = min(source) - 0.05 * (max(source) - min(source))
  res@trYMaxF  = max(source) + 0.05 * (max(source) - min(source))
  res@tiYAxisOn = True
  res@tiXAxisString = "Year"
  res@gsnStringFontHeightF = 0.016

  res@xyDashPatterns    = dashes
  res@xyLineThicknesses = thicks  ; make 2nd lines thicker
  res@xyLineColors      = colors  ; change line color

  copy_VarMeta(res_in, res)  ; copy passed resources

  ; ***************************************
  ; panel first two plots
  ; ***************************************

  plot(0) = gsn_csm_xy(wks, source_x, source, res)  ; create plot

  if (diag_script_info@ts_anomaly .eq. "anom") then

    if (diag_script_info@ref_value) then
      ; Add right panel with mean values in anomaly plot

      var = fspan(0., 2., 3)
      mean = new((/dimsizes(source_mean), 3/), double)
      mean(:, 0) = source_mean(:)
      mean(:, 1) = source_mean(:)
      mean(:, 2) = source_mean(:)

      plot2 = gsn_csm_xy(wks, var, mean, res0)  ; create plot
    end if

    ; ---------------------------------------------------------------------
    ; yellow shading of reference period
    ; ---------------------------------------------------------------------
    gsres                = True
    gsres@gsFillColor    = "yellow"
    gsres@gsFillOpacityF = 0.1
    xbox = (/ref_start, ref_end, ref_end, ref_start, ref_start/)
    ybox = (/res@trYMinF, res@trYMinF, res@trYMaxF, res@trYMaxF, res@trYMinF/)
    newplot00 = gsn_add_polygon(wks, plot(0), xbox, ybox, gsres)

    txtres               = True
    txtres@txFont        = "helvetica-bold"
    txtres@txFontColor   = "yellow3"
    txtres@txFontHeightF = 0.013
    ref_txt = gsn_add_text(wks, plot, "reference period", \
                           (ref_start + ref_end) / 0.5, -1., txtres)

  end if

  ; ---------------------------------------------------------------------
  ; add lines and names for volcanic eruptions
  ; ---------------------------------------------------------------------
  if (diag_script_info@volcanoes .eq. True) then
    vres = True
    vres@gsLineColor       = "grey20"
    vres@gsLineThicknessF  = 1.
    vres@gsLineDashPattern = 2.
    txres               = True
    txres@txAngleF      = 90
    txres@txFont        = "helvetica-bold"
    txres@txFontColor   = "grey20"
    txres@txFontHeightF = 0.013

    yy = (/res@trYMinF, res@trYMaxF/)

    vol_name1 = gsn_add_text(wks, plot, "Krakatoa", 1885, 1.1, txres)
    xx = (/1883.0, 1883.0/)
    vol1 = gsn_add_polyline(wks, plot(0), xx, yy, vres)
    vol_name2 = gsn_add_text(wks, plot, "Santa Maria", 1904, 1.1, txres)
    xx = (/1902.0, 1902.0/)
    vol2 = gsn_add_polyline(wks, plot(0), xx, yy, vres)
    vol_name3 = gsn_add_text(wks, plot, "Agung", 1965, 1.1, txres)
    xx = (/1963.0, 1963.0/)
    vol3 = gsn_add_polyline(wks, plot(0), xx, yy, vres)
    vol_name4 = gsn_add_text(wks, plot, "El Chichon", 1984, 1.1, txres)
    xx = (/1982.0, 1982.0/)
    vol4 = gsn_add_polyline(wks, plot(0), xx, yy, vres)
    vol_name5 = gsn_add_text(wks, plot, "Pinatubo", 1993, 1.1, txres)
    xx = (/1991.0, 1991.0/)
    vol5 = gsn_add_polyline(wks, plot(0), xx, yy, vres)
  end if

  ; ---------------------------------------------------------------------
  ; add line at 0.0
  ; ---------------------------------------------------------------------
  res_lines                   = True      ; polyline mods desired
  res_lines@tfPolyDrawOrder   = "Predraw"
  res_lines@gsLineColor       = "grey20"  ; line color
  res_lines@gsLineThicknessF  = 1.        ; line thicker
  res_lines@gsLineDashPattern = 1.        ; dash pattern

  xx = (/res@trXMinF, res@trXMaxF/)
  yy = (/0.0, 0.0/)
  dum0 = gsn_add_polyline(wks, plot(0), xx, yy, res_lines)
  newplot = gsn_attach_plots(plot(0), plot2, res, res0)

  ; ***********************************************
  ; legend resources
  ; ***********************************************

  lgres                    = True
  lgres@lgItemType         = "Lines"   ; show lines only (default)
  lgres@lgLabelFontHeightF = .05       ; set the legend label font thickness
  lgres@vpWidthF           = 0.15      ; width of legend (NDC)
  lgres@vpHeightF          = 0.5       ; height of legend (NDC)
  lgres@lgPerimColor       = "gray"    ; draw the box perimeter in orange
  lgres@lgPerimThicknessF  = 1.0       ; thicken the box perimeter

  lgres@lgLineColors      = colors
  lgres@lgDashIndexes     = dashes
  lgres@lgLineThicknesses = thicks + 0.5
  labels = annots
  nitems = dimsizes(labels)
  lgres@lgItemOrder = ispan(nitems - 1, 0, 1)

  ; Create legend
  lbid = gsn_create_legend(wks, nitems, labels, lgres)

  amres = True
  amres@amParallelPosF   = 0.81
  amres@amOrthogonalPosF = 0.0
  annoid1 = gsn_add_annotation(plot(0), lbid, amres)

  resP = True
  resP@gsnMaximize = True
  resP@gsnPaperOrientation = "portrait"
  resP@gsnPaperMargin =  0.8

  gsn_panel(wks, plot, (/1, 1/), resP)

  leave_msg(scriptname, funcname)

end

; #############################################################################
undef("timeseries_station")
function timeseries_station(wks_in[1],
                            source,
                            varname[1]: string)
;
; Arguments
;    wks_in: workstations (graphic object or default will be used).
;    source: data to be plotted or a NetCDF filename with data.
;              @stname: station name
;              @stlat: station latitude
;              @stlon: station longitude
;              @stalt: station altitude
;    varname: variable name, needed for netCDF files with multiple variables
;
; Source prototype
;    source[*][*]
;    source!0 = model
;    source!1 = time or year
;
; Return value
;    A graphic variable.
;
; Description
;    Creates a time series plot for station data.
;
; Caveats:
;    Selection of defaults for res almost arbitrary.
;
; Modification history:
;    20140325-righi_mattia: written.
;
local funcname, scriptname
begin

  funcname = "timeseries_station"
  scriptname = "diag_scripts/shared/plot/xy_line.ncl"
  enter_msg(scriptname, funcname)

  ; Get data, either directly or via netCDF file
  if (typeof(source) .eq. "string") then
    data = ncdf_read(source, varname)
  else
    data = source
    copy_VarMeta(source, data)
  end if

  ; Retrieve basic metadata from data
  defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
  if (any(varname .eq. defaults)) then
    var = att2var(data, "var")
  else
    var = varname
  end if

  ; Check if a valid wks has been provided, otherwise invoke default
  wks_out = get_wks(wks_in, diag_script, var)

  ; Default resources
  res = True
  res@gsnDraw             = False
  res@gsnFrame            = False
  res@tiMainFont          = 22
  res@tiMainString        = var
  res@tiYAxisString       = "[" + format_units(data@units) + "]"
  res@trYMinF             = 0.
  res@trYMaxF             = max(1.3 * data)
  res@tmXTOn              = False
  res@tmXBLabelAngleF     = 45.
  res@tmXBLabelJust       = "TopRight"
  res@tmXBMode            = "Explicit"
  res@tmLabelAutoStride   = True
  if (diag_script_info@time_avg.eq."monthly") then
    res@trXMinF             = min(data&time) - 50
    res@trXMaxF             = max(data&time) + 50
    res@tmXBValues          = data&time(::12)
    res@tmXBMinorValues     = data&time
    years                   = cd_calendar(data&time, 0)
    res@tmXBLabels          = years(::12, 0)
  else
    res@trXMinF             = min(data&year) - 0.2
    res@trXMaxF             = max(data&year) + 0.2
    res@tmXBValues          = data&year
    res@tmXBLabels          = data&year
  end if

  if (isatt(data, "long_name")) then
    res@gsnCenterString = data@long_name
  end if

  ; Override defaults with "res_" attributes of "data"
  res_new = att2var(data, "res_")
  copy_VarMeta(res_new, res)

  ; Plot
  dimnames = getvardims(data)
  plot = gsn_csm_xy(wks_out, data&$dimnames(1)$, data, res)

  ; Append legend
  resT = True
  resT@txFont = 21
  resT@txFontColor = "black"
  resT@txFontHeightF = 0.016
  if (.not.diag_script_info@legend_outside) then
    resT@txJust = "CenterLeft"
    resL = True
    xpos = 0.215
    ypos = 0.785
    do mID = 0, dimsizes(data&model) - 1
      resL@gsLineColor = res@xyLineColors(mID)
      resL@gsLineThicknessF = res@xyLineThicknesses(mID)
      resL@gsLineDashPattern = res@xyDashPatterns(mID)
      gsn_polyline_ndc(wks_out, (/xpos, xpos + 0.045/), \
                       (/ypos, ypos/), resL)
      gsn_text_ndc(wks_out, data&model(mID), xpos + 0.055, ypos, resT)
      ypos = ypos - 0.022
    end do
  end if

  ; Append station info
  resT@txJust = "CenterRight"
  xpos = 0.785
  ypos = 0.785

  if (isatt(data, "stname")) then
    gsn_text_ndc(wks_out, data@stname, xpos, ypos, resT)
    ypos = ypos - 0.022
  end if

  if (isatt(data, "stlat")) then
    latinfo = sprintf("%4.2f", abs(data@stlat)) + \
      where(data@stlat .lt. 0, " S", " N")
    gsn_text_ndc(wks_out, latinfo, xpos, ypos, resT)
    ypos = ypos - 0.022
  end if

  if (isatt(data, "stlon")) then
    if (data@stlon .gt. 180) then
      loninfo = sprintf("%5.2f", abs(data@stlon - 360.)) + " W"
    else if (data@stlon .lt. 0) then
      loninfo = sprintf("%5.2f", abs(data@stlon)) + " W"
    else
      loninfo = sprintf("%5.2f", data@stlon) + " E"
    end if
    end if
    gsn_text_ndc(wks_out, loninfo, xpos, ypos, resT)
    ypos = ypos - 0.022
  end if

  if (isatt(data, "stalt")) then
    if (.not.ismissing(data@stalt)) then
      gsn_text_ndc(wks_out, sprintf("%4.0f", data@stalt) + " m", \
                   xpos, ypos, resT)
    end if
  end if

  leave_msg(scriptname, funcname)
  return(plot)

end

; #############################################################################
undef("cycle_plot")
function cycle_plot(wks_in[1],
                    source,
                    varname[1]: string,
                    items: list)
;
; Arguments
;    wks_in: workstations (graphic object or default will be used).
;    source: data to be plotted or a NetCDF filename with data.
;    varname: variable name in the file.
;    items: list of input_file_info metadata
;
; Source prototype
;    source(*, *, 2)
;    source!0 = model
;    source!1 = month or season
;    source!2 = statistic
;    source@legend_outside = draw a legend withing the plot or in a separate
;                            file
;
; Return value
;    A graphic variable.
;
; Description
;    Draw an annual or seasonal cycle plot.
;
; Caveats
;
; Modification history
;    20131206-frank_franziska: written.
;
local funcname, scriptname, wks_out, wks_in, data, source, res, atts, base, \
  varname
begin

  funcname = "cycle_plot"
  scriptname = "diag_scripts/shared/plot/xy_line.ncl"
  enter_msg(scriptname, funcname)

  ; Get data, either directly or via netCDF file
  if (typeof(source) .eq. "string") then
    data = ncdf_read(source, varname)
  else
    data = source
    copy_VarMeta(source, data)
  end if

  ; Retrieve basic metadata from data
  defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
  if (any(varname .eq. defaults)) then
    var = att2var(data, "var")
    DIAG_SCRIPT = att2var(data, "diag_script")
  else
    var = varname
  end if

  if (.not.isatt(diag_script_info, "styleset")) then
    diag_script_info@styleset = "DEFAULT"
  end if

  ; Select colors and other plotting attributes from the project style files
  ; See ./diag_scripts/shared/plot/style.ncl
  colors = project_style(items, diag_script_info, "colors")
  dashes = project_style(items, diag_script_info, "dashes")
  thicks = project_style(items, diag_script_info, "thicks")
  annots = project_style(items, diag_script_info, "annots")
  avgstd = project_style(items, diag_script_info, "avgstd")

  ; Check if a valid wks has been provided, otherwise invoke default
  wks_out = get_wks(wks_in, DIAG_SCRIPT, var)

  ; Define default plot resources
  res                  = True
  res@gsnDraw          = False
  res@gsnFrame         = False
  res@gsnMaximize      = True
  res@gsnPaperOrientation = "landscape"
  res@xyLineThicknesses = thicks + 0.5
  res@xyLineColors      = colors
  res@xyMonoDashPattern = False
  res@xyDashPatterns    = dashes

  coord_names = getvardims(data)
  x_axis = data&$coord_names(1)$
  x_axis_num = ispan(1, dimsizes(data&$coord_names(1)$), 1)

  res@trXMinF  = min(x_axis_num) - 0.05 * (max(x_axis_num) - min(x_axis_num))
  res@trXMaxF  = max(x_axis_num) + 0.25 * (max(x_axis_num) - min(x_axis_num))
  if (data@legend_outside) then
    res@trXMinF  = min(x_axis_num) - 0.3
    res@trXMaxF  = max(x_axis_num) + 0.3
  end if
  ymax = max(data(:, :, 0) + data(:, :, 1))
  ymin = min(data(:, :, 0) - data(:, :, 1))
  res@trYMinF = ymin - 0.05 * (ymax - ymin)
  res@trYMaxF = ymax + 0.05 * (ymax - ymin)
  res@tiMainFontHeightF    = 0.025
  res@vpHeightF = 0.4
  res@vpWidthF  = 0.8
  res@tmXBMode   = "Explicit"
  res@tmXTOn = False
  res@tmXBValues = x_axis_num
  res@tmXBLabels = x_axis

  res@tiMainString = "Annual cycle"
  if (isatt(data, "units"))
    res@tiYAxisString = "[" + format_units(data@units) + "]"
  else
    res@tiYAxisString = varname
  end if

  ; Override defaults with "res_" attributes of "data"
  res_new = att2var(data, "res_")
  copy_VarMeta(res_new, res)

  ; Create plot
  plot = gsn_csm_xy(wks_out, x_axis_num, data(:, :, 0), res)
  plot@outfile = wks_out@name

  ; Add standard deviation
  if (isatt(data, "plot_stddev")) then
    if (.not.(data@plot_stddev.eq."none")) then

      ; Stddev
      res_stddev = True
      res_stddev@tfPolyDrawOrder = "Predraw"
      res_stddev@gsFillOpacityF = 0.15

      ; Assign edges of polygons, second round reverted
      data_new = new((/dimsizes(data&model), 2 * dimsizes(x_axis)/), float)
      data_new(:, 0:dimsizes(x_axis) - 1) = data(:, :, 0)
      data_new(:, 0:dimsizes(x_axis) - 1) = data(:, :, 0) + data(:, :, 1)
      data_new(:, dimsizes(x_axis):) = data(:, ::-1, 0)
      data_new(:, dimsizes(x_axis):) = data(:, ::-1, 0) - data(:, ::-1, 1)
      data_new&$coord_names(1)$(0:dimsizes(x_axis) - 1) = x_axis_num
      data_new&$coord_names(1)$(dimsizes(x_axis):) = x_axis_num(::-1)

      ;: Loop over models
      shading = new((/dimsizes(data(:, 0, 0))/), "graphic")
      do imod = 0, dimsizes(shading) - 1
        if ((any(data@plot_stddev.eq.annots(imod))).or.\
            (data@plot_stddev.eq."all").or.\
            ((data@plot_stddev.eq."ref_model").and.\
            (any(data@ref_model.eq.annots(imod))))) then
          if (dimsizes(dimsizes(colors)).eq.2) then
            ; Assign model corresponding color (RGB arrays)
            res_stddev@gsFillColor = colors(imod, :)
          else
            ; Assign model corresponding color (String arrays)
            res_stddev@gsFillColor = colors(imod)
          end if
          shading(imod) = \
            gsn_add_polygon(wks_out, plot, \
                            tofloat(data_new&$coord_names(1)$), \
                            data_new(imod, :), res_stddev)
          str = "shade" + imod
          plot@$str$ = shading(imod)
        end if
      end do
      delete(data_new)
    end if
  end if

  if (data@legend_outside) then

    lg_outfile = wks_out@legendfile

    ; Attach line styles
    styles = True
    styles@colors = colors
    styles@dashes = dashes
    styles@thicks = thicks
    create_legend_lines(annots, styles, lg_outfile, "lines")

  else

    nitems = dimsizes(annots)

    ; Resources for a customized legend
    lgres                    = True
    lgres@lgMonoLineThickness = False
    lgres@lgLabelFontHeightF = 0.08
    lgres@vpWidthF           = 0.15
    lgres@vpHeightF          = 0.4
    lgres@lgPerimOn          = False
    lgres@lgMonoDashIndex    = False
    lgres@lgLineColors       = colors(::-1)
    lgres@lgDashIndexes      = dashes(::-1)
    lgres@lgLineThicknesses  = thicks(::-1) + 0.5

    ; Create legends within plot
    lbid_ta = gsn_create_legend(wks_out, nitems, annots(::-1), lgres)

    ; Add legend at the correct position
    ; Point (0, 0) is the dead center of the plot. Point (0, .5)
    ; is center, flush bottom. Point (0.5, 0.5) is flush bottom,
    ; flush right.
    amres                  = True
    amres@amJust           = "TopRight"  ; reference corner of box
    amres@amParallelPosF   = 0.5         ; Move legend to +right, -left
    amres@amOrthogonalPosF = -0.5        ; +down, -up

    ; Add legend to plot
    annoid = gsn_add_annotation(plot, lbid_ta, amres)
    plot@annoid = annoid

  end if

  leave_msg(scriptname, funcname)
  return(plot)

end

; #############################################################################
undef("errorbar_plot")
function errorbar_plot(wks_in[1],
                       source,
                       varname[1]:string)
;
; Arguments
;    wks_in: workstations (graphic object or default will be used).
;    source: data to be plotted or a NetCDF filename with data.
;    varname: variable name in the file.
;
; Source prototype
;
; Return value
;    A graphic variable.
;
; Description
;    source = (2, npoints)
;    source(0, :) = mean
;    source(1, :) = standard deviation
;    source!0 = statistic
;    source!1 = model
;
; Caveats
;
; Modification history
;    20151105-righi_mattia: written.
;
local funcname, scriptname, data, defaults, var, wks, res, ren_new, resL, \
  error_bar
begin

  funcname = "errorbar_plot"
  scriptname = "diag_scripts/shared/plot/xy_line.ncl"
  enter_msg(scriptname, funcname)

  ; Get data, either directly or via netCDF file
  if (typeof(source) .eq. "string") then
    data = ncdf_read(source, varname)
  else
    data = source
    copy_VarMeta(source, data)
  end if

  ; Retrieve basic metadata from data
  defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
  if (any(varname .eq. defaults)) then
    var = att2var(data, "var")
    diag_script = att2var(data, "diag_script")
  else
    var = varname
  end if

  ; Check if a valid wks has been provided, otherwise invoke default
  wks = get_wks(wks_in, diag_script, var)

  ; Default resources
  res = True
  res@gsnFrame = False
  res@gsnDraw = False
  res@xyMarkLineMode = "Markers"
  res@xyMarkers = 16
  res@xyMarkerSizeF = 0.01
  res@tiYAxisString = varname
  if (isatt(data, "units")) then
    res@tiYAxisString = \
      res@tiYAxisString + " [" + format_units(data@units) + "]"
  end if
  res@vpWidthF = 0.7
  res@vpHeightF = 0.25
  res@txFontHeightF = 0.01
  res@tiMainFontHeightF = 0.025
  res@tmXBMode = "Explicit"
  res@trXMinF = -1
  res@trXMaxF = dimsizes(data&model)
  res@tmXBValues = ispan(0, dimsizes(data&model) - 1, 1)
  res@trYMinF = min(data(0, :) - data(1, :) / 2.)
  res@trYMinF = where(res@trXMinF.lt.0, 1.1 * res@trYMinF, 0.9 * res@trYMinF)
  res@trYMaxF = max(data(0, :) + data(1, :) / 2.)
  res@trYMaxF = where(res@trXMaxF.gt.0, 1.1 * res@trYMaxF, 0.9 * res@trYMaxF)
  res@tiMainString = varname + " - " + diag_script_info@region
  res@tmXBLabels = data&model
  res@tmXBLabelAngleF = 45.
  res@tmXBLabelJust = "CenterRight"

  ; Override defaults with "res_" attributes of "data"
  res_new = att2var(data, "res_")
  copy_VarMeta(res_new, res)

  ; Draw plot
  plot = gsn_csm_y(wks, data(0, :), res)

  ; Draw errorbars
  resL = True
  resL@gsLineThicknessF = 1.
  error_bar = new(dimsizes(data&model), graphic)
  do imod = 0, dimsizes(data&model) - 1
    error_bar(imod) = \
      gsn_add_polyline(wks, plot, (/imod, imod/), \
                       (/data(0, imod) + data(1, imod) / 2, \
                         data(0, imod) - data(1, imod) / 2/), resL)
    str = "bar" + imod
    plot@$str$ = error_bar(imod)
  end do

  leave_msg(scriptname, funcname)
  return(plot)

end
