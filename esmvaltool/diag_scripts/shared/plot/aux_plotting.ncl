; #############################################################################
; AUXILIARY PLOTTING ROUTINES
; #############################################################################
; Please consider using or extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
; Contents:
;
;    procedure create_legend_lines
;    function output_type
;    procedure copy_VarAtt_sel
;    function panelling
;    function get_plot_dir
;    function get_outfile_name
;    function get_wks
;    procedure add_markers
;    procedure add_num_markers
;    procedure add_errorbar
;    procedure horizontal_whiskers
;    procedure add_prediction_error
;    function month_sel
;    function lat_names
;    procedure add_line
;    procedure add_scatt
;    procedure add_legend
;    function calcRegCoeffs
;    function genZonalMeans
;    function calcMeanAnnCycleMonthly
;    function calcMeanAnnCycleAnnual
;    procedure rmMeanAnnCycle
;    function apfiltersmooth
;    procedure smoothAnomalies
;    function clmMon2clmDayn
;
; #############################################################################

load "$diag_scripts/../interface_scripts/auxiliary.ncl"
load "$diag_scripts/../interface_scripts/logging.ncl"

; #############################################################################
undef("create_legend_lines")
procedure create_legend_lines(labels:string, \
                              styles, \
                              outfile:string,\
                              opt)
;
; Arguments
;    labels:  labels of the lines
;    styles: style of the lines in the plot
;              @colors
;              @dashes  -> optional (line dashes
;              @thicks  -> optional (line/marker thickness)
;              @markers -> optional (marker index)
;              @sizes   -> optional (marker size)
;    outfile: outfile directory
;    opt: "lines" or "markers" legend
;
; Description
;    Draws an extra plot with a legend, specified by labels and styles.
;
; Caveats
;
; Modification history
;    20150511-lauer_axel: added safe lower limits for panelling plot
;    20150508-righi_mattia: added lines/markers option.
;    20150120-gottschaldt_klaus-dirk: remove pre-existing file type suffix
;    20140305-righi_mattia: modified to plot always as epsi format.
;    20140219-frank_franziska: written.
;
local funcname, scriptname, region, temp, outfile, n_suff, wks_legend
begin

  funcname = "create_legend_lines"
  scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
  enter_msg(scriptname, funcname)

  ; Open workstation
  wks_legend = gsn_open_wks("epsi", outfile)

  ; General resources
  res = True
  res@tmBorderThicknessF = 0.
  res@tmBorderLineColor = "white"
  res@tmXTOn = False
  res@tmXBOn = False
  res@tmYLOn = False
  res@tmYROn = False
  res@trXMinF = 0
  res@trXMaxF = 1
  res@trYMinF = 0
  res@trYMaxF = 1
  res@vpXF = 0.0
  res@vpYF = 1.0

  ; Resources for polyline and text
  resL = True
  resT = True
  resT@txFontHeightF = 0.013
  resT@txJust = "CenterLeft"

  ; Set number of lines and labels
  dim_LAB = dimsizes(labels)

  ; Rows and columns
  deltax = 0.06
  deltay = 0.03
  nrowmax = 32
  ncols = 3
  res@vpWidthF = 0.8
  if (dim_LAB .gt. nrowmax * ncols) then
    ncols = 4
    res@vpWidthF = 1.0
    if (dim_LAB .gt. nrowmax * ncols) then
      error_msg("f", "xy_line.ncl", funcname, "too many datasets")
    end if
  end if
  nrows = toint(dim_LAB/tofloat(ncols))

  line_xpos = new((/ncols, 2/), float)
  do jj = 0, ncols - 1
    line_xpos(jj, 0) = 0.01 + jj / tofloat(ncols)
  end do
  line_xpos(:, 1) = line_xpos(:, 0) + deltax
  text_pos = line_xpos(:, 1) + 0.01
  res@vpHeightF = 0.04 + deltay * (nrows - 1)

  ; Draw blank plot
  blank = gsn_blank_plot(wks_legend, res)

  ; Loop over labels, draw lines and text
  ypos = 0.98
  ii = 0
  do while (ii.lt.dim_LAB)
    ; Set color
    resL@gsLineColor = styles@colors(ii)
    resL@gsMarkerColor = styles@colors(ii)
    ; Set dash pattern
    if (isatt(styles, "dashes")) then
      resL@gsLineDashPattern    = styles@dashes(ii)
    end if
    ; Set thickness
    if (isatt(styles, "thicks")) then
      resL@gsLineThicknessF    = styles@thicks(ii)
      resL@gsMarkerThicknessF  = styles@thicks(ii)
    end if
    ; Set marker index
    if (isatt(styles, "markers")) then
      resL@gsMarkerIndex       = styles@markers(ii)
    end if
    ; Set marker size
    if (isatt(styles, "sizes")) then
      resL@gsMarkerSizeF = styles@sizes(ii)
    end if
    ; Draw legend
    if (opt.eq."lines") then
      gsn_polyline_ndc(wks_legend, line_xpos(mod(ii, ncols), :), \
                       (/ypos, ypos/), resL)
      gsn_text_ndc(wks_legend, labels(ii), text_pos(mod(ii, ncols)), \
                   ypos, resT)
    else if (opt.eq."markers") then
      gsn_polymarker_ndc(wks_legend, avg(line_xpos(mod(ii, ncols), :)), \
                         ypos, resL)
      gsn_text_ndc(wks_legend, labels(ii), text_pos(mod(ii, ncols)) - 0.02, \
                   ypos, resT)
    else
      error_msg("f", "xy_line.ncl", funcname, "unrecognized option " + opt)
    end if
    end if
    if (mod(ii, ncols).eq.ncols - 1) then
      ypos = ypos - deltay  ; go to the next row
    end if
    ii = ii + 1
  end do
  draw(blank)
  frame(wks_legend)

  leave_msg(scriptname, funcname)

end

; #############################################################################
undef("output_type")
function output_type()
;
; Arguments
;
; Return value
;    A string with the output file type
;
; Description
;    Provides a default, if file type is not explicitly specified
;
; Caveats
;
; Modification history
;    20131028-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, file_type
begin

  funcname = "output_type"
  scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
  enter_msg(scriptname, funcname)

  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if
  file_type = str_lower(file_type)

  leave_msg(scriptname, funcname)
  return(file_type)

end

; #############################################################################
undef("copy_VarAtt_sel")
procedure copy_VarAtt_sel(var1, var2, sel: string)
;
; Arguments
;    var1: variable of any type and dimension with source attributes
;    var2: variable of any type and dimension, shall receive selected
;          attributes from var1
;    sel: string (or list of strings) that specify the BEGINNING letters
;         of the attributes to copy
;
; Description
;    Intended to copy selected plot ressources for use in a different
;    plot routine that may not allow all ressources from var1.
;    Written for function legend_lines (legends.ncl).
;
; Caveats:
;
; Modification history
;    20141227-gottschaldt_klaus-dirk written.
;
local funcname, scriptname, atts, i1, i2, n
begin

  funcname = "copy_VarAtt_sel"
  scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
  enter_msg(scriptname, funcname)

  ; Check if var1 has attributes at all
  atts = getvaratts(var1)
  if (all(atts.eq."missing")) then
    error_msg("w", scriptname, funcname, "no attributes to copy")
    return(var2)  ; return var2 unchanged
  end if

  ; Copy selected attributes from var1
  ; Existing attributes of the same name will be overwritten in var2.
  do i2 = 0, dimsizes(sel) - 1
    n = strlen(sel(i2))
    do i1 = 0, dimsizes(atts) - 1
      if (str_get_cols(atts(i1), 0, n - 1) .eq. sel(i2)) then
        if (isatt(var2, atts(i1))) then
          delete(var2@$atts(i1)$)
        end if
        var2@$atts(i1)$ = var1@$atts(i1)$
      end if
    end do
  end do

  leave_msg(scriptname, funcname)

end

; #############################################################################
undef("panelling")
function panelling(wks,
                   plots:graphic,
                   nvert[1]:integer,
                   nhori[1]:integer,
                   pres_in[1]:logical)
;
; Arguments
;    wks: graphics object associated with a file.
;    plots: graphics variable containing one or more plots.
;    nvert: Maximum allowed number of plots per page (vertical).
;    nhori: Maximum allowed number of plots per page (horizontal).
;    pres_in: if it is a logical variable, attributes are used to override
;             panelling defaults
;
; Return value
;    A string with the output file name
;
; Description
;    Writes a graphics file
;
; Caveats
;    Might not be fit to be used with non-NCL routines.
;
; Modification history
;    20131112-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, plots, info, wks, outfile
begin

  funcname = "panelling"
  scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
  enter_msg(scriptname, funcname)

  ; Check if valid wks has been provided
  if (typeof(wks).ne."graphic") then
    error_msg("f", scriptname, funcname, "panelling of " + \
              "multiple plots needs valid wks")
  end if

  ; Get outfile name
  if (isatt(wks, "fullname")) then
    outfile = wks@fullname
  else
    outfile = wks@name
    error_msg("w", scriptname, funcname, "wks@fullname " + \
              "missing, consider to use get_wks to open wks")
  end if

  ; Panel parameters
  npages = toint(ceil(tofloat(dimsizes(plots)) / tofloat(nvert * nhori)))

  ; Draw the panel
  pres = True
  pres@gsnMaximize       = True
  pres@gsnPanelYWhiteSpacePercent = 5
  pres@gsnPanelXWhiteSpacePercent = 5
  if (isatt(diag_script_info, "PanelTop")) then
    top = tofloat(diag_script_info@PanelTop)
  else
    top = 0.99  ; default
  end if
  pres@gsnPanelTop = top
  if (npages .gt. 1) then
    pres@gsnPaperOrientation   = "portrait"
    if (typeof(pres_in).eq."logical") then  ; Override defaults
      copy_VarMeta(pres_in, pres)
    end if
    ia = 0
    ie = nvert * nhori - 1
    ie = min((/dimsizes(plots) - 1, ie/))
    do ipage = 0, npages - 1
      gsn_panel(wks, plots(ia:ie), (/nvert, nhori/), pres)
      ia = ie + 1
      ie = ia + (nvert * nhori) - 1
      ie = min((/dimsizes(plots) - 1, ie/))
    end do
  else
    pres@gsnPaperOrientation   = "auto"
    if (typeof(pres_in).eq."logical") then  ; Override defaults
      copy_VarMeta(pres_in, pres)
    end if
    nrow = toint(ceil(tofloat(dimsizes(plots)) / tofloat(nhori)))
    ncol = min((/nhori, dimsizes(plots)/))
    gsn_panel(wks, plots, (/nrow, ncol/), pres)
  end if

  leave_msg(scriptname, funcname)
  return(outfile)

end

; #############################################################################
undef("get_plot_dir")
function get_plot_dir()
;
; Arguments
;
; Return value
;    A string root path for plots
;
; Description
;    Provides a default, if plot_dir is not explicitely specified
;
; Caveats
;
; Modification history
;    20131104-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, plot_dir
begin

  funcname = "get_plot_dir"
  scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
  enter_msg(scriptname, funcname)

  plot_dir = config_user_info@plot_dir
  if (ismissing(plot_dir)) then
    plot_dir = "./plots/"
  end if
  if (.not. fileexists(plot_dir)) then
    system("mkdir -p " + plot_dir)
  end if

  leave_msg(scriptname, funcname)
  return(plot_dir)

end

; #############################################################################
undef("get_outfile_name")
function get_outfile_name(add_specs[1]:string)
;
; Arguments
;    add_specs: string containing specific elements to be added to the file
;               name if not of type string, only variable name is used.
;
; Return value
;    Complete outfile name incl. path, additional specifications and file type
;    plot_dir + add_specs + file_type
;
; Description
;    Fetches file_type, plot_dir
;    Fetches string with additional elements specified within diag script
;
; Modification history
;    20131204-senftleben_daniel: generalized naming.
;    20131104-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, outfile, file_type, plot_dir, diag_script_base, \
  add_specs, output_dir
begin

  funcname = "get_outfile_name"
  scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
  enter_msg(scriptname, funcname)

  ; Collect parts of outfile name
  file_type = output_type()
  plot_dir = get_plot_dir()

  ; Assemble final name
  if (typeof(add_specs).ne."string") then
    outfile = plot_dir + var + "." + file_type
  else
    outfile = plot_dir + add_specs + "." + file_type
  end if

  leave_msg(scriptname, funcname)
  return(outfile)

end

; #############################################################################
undef("get_wks")
function get_wks(wks_in,
                 diag_script[1]: string,
                 add_specs[1]: string)
;
; Arguments
;    wks_in: dummy or graphic object
;    diag_script: name of diagnostic script
;    add_specs: own specificactions to be added to file name, i.e. variable
;               name, etc. - needed for function get_outfile_name
;
; Return value
;    wks: graphic object
;
; Description
;    Provides a default wks, if wks_in is not of type "graphic".
;    Attribute wks@fullname is used to transfer the output file name, since
;    wks@name cuts off the path to the file name.
;
; Caveats
;
; Modification history
;    20131113-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, wks_in, diag_script, add_specs, file_type
begin

  funcname = "get_wks"
  scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
  enter_msg(scriptname, funcname)

  ; Check if a valid wks has been provided, otherwise invoke default
  if (typeof(wks_in).ne."graphic") then  ; Use defaults
    file_type = output_type()
    copy_VarMeta(wks_in, file_type)
    outfile = get_outfile_name(add_specs)
    wks = gsn_open_wks(file_type, outfile)
    wks@fullname = outfile
  else
    wks = wks_in
    copy_VarMeta(wks_in, wks)
  end if

  leave_msg(scriptname, funcname)
  return(wks)

end

; #############################################################################
undef("add_markers")
procedure add_markers(wks[1]:graphic,
                      plot[1]:graphic,
                      res_in[1]:logical,
                      xpos_in:numeric,
                      ypos_in:numeric)
;
; Arguments
;    wks: valid workstation, e.g. created by get_wks
;    plot: plot identifier, e.g. created by gsn_*
;    res_in: plot ressources that may override local function defaults
;    xpos_in: horizontal marker position(s)
;    ypos_in: vertical marker position(s)
;
; Return value
;    Attaches polyline IDs as attributes to plot.
;
; Description:
;    Adds markers to an existing plot.
;
; Caveats:
;
; Modification history
;    20140224-gottschaldt_klaus-dirk: written for use with profile_plev.ncl
;                                     in Emmons.ncl
;
local funcname, scriptname, res, str, atts, iatt, color
begin

  funcname = "add_markers"
  scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
  enter_msg(scriptname, funcname)

  ; Set default plot ressources
  res               = True
  res@gsMarkerSizeF = 0.006
  res@gsMarkerColor = "black"
  res@gsMarkerIndex = 16  ; solid circle

  ; Override defaults (after cleaning up for possible vector input)
  if (res_in) then
    atts = getvaratts(res_in)
    if (.not.all(ismissing(atts))) then
      do iatt = 0, dimsizes(atts) - 1
        if (isatt(res, atts(iatt))) then
          delete(res@$atts(iatt)$)
        end if
      end do
    end if
    copy_VarMeta(res_in, res)
  end if

  ; Plot symbol rim (slightly bigger, white polymarker background)
  color = res@gsMarkerColor  ; remember original color
  res@gsMarkerSizeF = res@gsMarkerSizeF + res@gsMarkerSizeF * 0.2
  res@gsMarkerColor(:) = "white"
  str = unique_string("polymarker")
  plot@$str$ = gsn_add_polymarker(wks, plot, xpos_in, ypos_in, res)

  ; Plot symbols
  res@gsMarkerColor(:) = color  ; restore original color
  res@gsMarkerSizeF = res@gsMarkerSizeF - res@gsMarkerSizeF * 0.2
  str = unique_string("polymarker")
  plot@$str$ = gsn_add_polymarker(wks, plot, xpos_in, ypos_in, res)

  leave_msg(scriptname, funcname)

end

; #############################################################################
undef("add_num_markers")
procedure add_num_markers(wks[1]:graphic,
                          plot[1]:graphic,
                          res_in[1]:logical,
                          xpos_in:numeric,
                          ypos_in:numeric)
;
; Arguments
;    wks: valid workstation, e.g. created by get_wks
;    plot: plot identifier, e.g. created by gsn_*
;    res_in: plot ressources that may override local function defaults
;    xpos_in: horizontal marker position(s), needs year as a coordinate
;    ypos_in: vertical marker position(s)
;
; Return value
;    Attaches polytext IDs as attributes to plot.
;
; Description:
;    Adds markers to an existing plot, being the last two intagers of a year.
;    If a horizontal (vertical) coordinate
;    has only one element, then this position is used for all markers.
;
; Caveats:
;
; Modification history
;    20150914-wenzel_sabrina: written.
;
local funcname, scriptname, xpos, ypos, res, str, atts, iatt, color
begin

  funcname = "add_num_markers"
  scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
  enter_msg(scriptname, funcname)

  ; set default plot ressources
  res = True
  res@txFontHeightF = 0.025
  res@txFontColor   = "black"

  ; override defaults (after cleaning up for possible vector input)
  if (res_in) then
    atts = getvaratts(res_in)
    if (.not.all(ismissing(atts))) then
      do iatt = 0, dimsizes(atts) - 1
        if (isatt(res, atts(iatt))) then
          delete(res@$atts(iatt)$)
        end if
      end do
    end if
    copy_VarMeta(res_in, res)
  end if

  atext = new(dimsizes(xpos_in), integer)

  do y = 0, dimsizes(xpos_in) - 1
    strs     = str_split_by_length(tostring(xpos_in&year(y)), 1)
    strs_y   = toint(strs(0) + strs(1) + "00")
    atext(y) = xpos_in&year(y) - strs_y
    delete([/strs, strs_y/])
  end do

  ; create vector, if horizontal or vertical position is given as scalar
  if (dimsizes(xpos_in).eq.1) then
    xpos = new(dimsizes(ypos_in), typeof(xpos_in))
    xpos = (/xpos_in/)
  else
    xpos = (/xpos_in/)
  end if
  if (dimsizes(ypos_in).eq.1) then
    ypos = new(dimsizes(xpos_in), typeof(ypos_in))
    ypos = (/ypos_in/)
  else
    ypos = (/ypos_in/)
  end if

  ; plot symbol
  str = unique_string("polymarker")
  plot@$str$ = \
    gsn_add_text(wks, plot, sprinti("%-0.2hi", atext), xpos, ypos, res)

  leave_msg(scriptname, funcname)

end

; #############################################################################
undef("add_errorbar")
procedure add_errorbar(wks[1]:graphic,
                       plot[1]:graphic,
                       res_in[1]:logical,
                       xpos_in:numeric,
                       ypos_in:numeric)
;
; Arguments
;    wks: valid workstation, e.g. created by get_wks
;    plot: plot identifier, e.g. created by gsn_*
;    res_in: plot ressources that may override local function defaults
;    xpos_in: 2D array with horizontal marker position(s) and stddev's,
;    ypos_in: 2D array with vertical marker position(s) and stddev's
;
; Return value
;    Attaches polytext IDs as attributes to plot.
;
; Description:
;    Adds error bars to an existing plot. If a horizontal (vertical)
;    coordinate has only one element, then this position is used for all
;    markers.
;
; Caveats:
;
; Modification history
;    20150914-wenzel_sabrina: written.
;
local funcname, scriptname, xpos, ypos, res, str, atts, iatt, color
begin

  funcname = "add_errorbar"
  scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
  enter_msg(scriptname, funcname)

  ; Set default plot ressources
  res = True
  res@gsLineThicknessF = 2

  colors = new(dimsizes(xpos_in(0, :)), string)
  colors = "black"

  ; Override defaults (after cleaning up for possible vector input)
  if (res_in) then
    atts = getvaratts(res_in)
    if (.not.all(ismissing(atts))) then
      do iatt = 0, dimsizes(atts) - 1
        if (isatt(res, atts(iatt))) then
          delete(res@$atts(iatt)$)
        end if
      end do
    end if
    copy_VarMeta(res_in, res)

    if (isatt(res, "gsLineColor")) then
      colors = res@gsLineColor
      delete(res@gsLineColor)
    end if
  end if

  ; Create vector, if horizontal or vertical position is given as scalar
  if (dimsizes(xpos_in(0, :)).eq.1) then
    xpos = new(dimsizes(ypos_in), typeof(xpos_in))
    xpos(:, :) = xpos_in
  else
    xpos = xpos_in
  end if

  if (dimsizes(ypos_in(0, :)).eq.1) then
    ypos = new(dimsizes(xpos_in), typeof(ypos_in))
    ypos(:, :) = ypos_in
  else
    ypos = ypos_in
  end if

  Xend = 0.01 * ((max(ypos(0, :)) + 0.2 * (abs(max(ypos(0, :))))) - \
                 (min(ypos(0, :)) - 0.2 * (abs(min(ypos(0, :))))))
  Yend = 0.01 * ((max(xpos(0, :)) + 0.2 * (abs(avg(xpos(0, :))))) - \
                 (min(xpos(0, :)) - 0.2 * (abs(avg(xpos(0, :))))))

  do i = 0, dimsizes(xpos(0, :)) - 1
    res@gsLineColor = colors(i)

    if (.not.any(ismissing((/xpos(0, i), ypos(0, i)/))))
      Xbarstr = unique_string("Xbarerror")
      Ybarstr = unique_string("Ybarerror")
      Xendstr = unique_string("Xenderror")
      Yendstr = unique_string("Yenderror")

      plot@$Xbarstr$ = \
        gsn_add_polyline(wks, plot, \
                         (/xpos(0, i) + xpos(1, i), \
                           xpos(0, i) - xpos(1, i)/), \
                         (/ypos(0, i), ypos(0, i)/), res)
      plot@$Ybarstr$ = \
        gsn_add_polyline(wks, plot, \
                         (/xpos(0, i), xpos(0, i)/), \
                         (/ypos(0, i) + ypos(1, i), \
                           ypos(0, i) - ypos(1, i)/), res)
      plot@$Xendstr$ = \
        gsn_add_polyline(wks, plot, \
                         (/xpos(0, i) + xpos(1, i), \
                           xpos(0, i) + xpos(1, i), \
                           xpos(0, i) + xpos(1, i), \
                           xpos(0, i) - xpos(1, i), \
                           xpos(0, i) - xpos(1, i), \
                           xpos(0, i) - xpos(1, i)/), \
                         (/ypos(0, i) - Xend, \
                           ypos(0, i) + Xend, \
                           ypos(0, i), ypos(0, i), \
                           ypos(0, i) - Xend, \
                           ypos(0, i) + Xend/), res)
      plot@$Yendstr$ = \
        gsn_add_polyline(wks, plot, \
                         (/xpos(0, i) - Yend, \
                           xpos(0, i) + Yend, \
                           xpos(0, i), \
                           xpos(0, i), \
                           xpos(0, i) - Yend, \
                           xpos(0, i) + Yend/), \
                         (/ypos(0, i) + ypos(1, i), \
                           ypos(0, i) + ypos(1, i), \
                           ypos(0, i) + ypos(1, i), \
                           ypos(0, i) - ypos(1, i), \
                           ypos(0, i) - ypos(1, i), \
                           ypos(0, i) - ypos(1, i) /), res)
    end if
  end do

  leave_msg(scriptname, funcname)

end

; #############################################################################
undef("horizontal_whiskers")
procedure horizontal_whiskers(wks[1]:graphic,
                              plot[1]:graphic,
                              res_in[1]:logical,
                              xmin_in:numeric,
                              xmax_in:numeric,
                              ypos_in: numeric)
;
; Arguments
;    wks: valid workstation, e.g. created by get_wks.
;    plot: plot identifier, e.g. created by gsn_*.
;    res_in: plot ressources that may override local function defaults.
;    xmin_in: vector of whiskers' left ends (same size as xmax & y).
;    xmax_in: vector of whiskers' right ends (same size as xmin & y).
;    ypos_in: vector of whiskers' vertical positions (must have same size as
;             xmax & xmin).
;
; Return value
;    Attaches polyline IDs as attributes to plot.
;
; Description
;    Creates vectors suitable as input for gsn_add_polyline:
;      x = (/xmin1,xmax1,_FillValue,xmin2,xmax2,_FillValue, .../)
;      y = (/ypos1,ypos1,_FillValue,ypos2,ypos2,_FillValue, .../)
;    The separation by _FillValue results in individual whiskers.
;    No whisker is created where xmin, xmax or ypos is missing.
;
; Caveats
;
; References
;    www.ncl.ucar.edu/Document/Graphics/Interfaces/gsn_add_polyline.shtml
;
; Modification history
;    20140224-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, n, i, j, xmin, xmax, ypos, imiss, ivalid, result, \
  x, y, res, str
begin

  funcname = "horizontal_whiskers"
  scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
  enter_msg(scriptname, funcname)

  ; Check sizes of input vectors
  n = dimsizes(ypos_in)
  if (dimsizes(xmin_in).ne.n .or. dimsizes(xmax_in).ne.n) then
    error_msg("f", scriptname, funcname, "all input vectors " + \
              "must have same size")
  end if

  ; Unify _FillValue and type
  ypos = tofloat(ypos_in)  ; #####################
  copy_VarMeta(ypos_in, ypos)
  ypos@_FillValue = default_fillvalue("float")
  ypos = where(ismissing(ypos_in), ypos@_FillValue, ypos)
  xmin = tofloat(xmin_in)  ; #####################
  copy_VarMeta(xmin_in, xmin)
  xmin@_FillValue = default_fillvalue("float")
  xmin = where(ismissing(xmin_in), xmin@_FillValue, xmin)
  xmax = tofloat(xmax_in)  ; #####################
  copy_VarMeta(xmax_in, xmax)
  xmax@_FillValue = default_fillvalue("float")
  xmax = where(ismissing(xmax_in), xmax@_FillValue, xmax)

  ; Find positions where no input element is missing
  imiss = ispan(0, n - 1, 1)
  imiss@_FillValue = default_fillvalue("integer")
  imiss = where(ismissing(ypos), imiss@_FillValue, imiss)
  imiss = where(ismissing(xmin), imiss@_FillValue, imiss)
  imiss = where(ismissing(xmax), imiss@_FillValue, imiss)
  ivalid = ind(.not. ismissing(imiss))

  if (all(ismissing(ivalid))) then
    log_debug(" No valid data to plot whiskers: Skipping")
  else

    ; Create vectors suited for gsn_add_polyline
    x = new(3 * dimsizes(ivalid) - 1, typeof(xmin))
    x@_FillValue = xmin@_FillValue
    x(:) = x@_FillValue
    y = new(3 * dimsizes(ivalid) - 1, typeof(ypos))
    y@_FillValue = xmax@_FillValue
    y(:) = y@_FillValue
    do i = 0, dimsizes(ivalid) - 1
      j = ivalid(i)
      x(i * 3)   = xmin(j)
      x(i * 3 + 1) = xmax(j)
      y(i * 3)   = ypos(j)
      y(i * 3 + 1) = ypos(j)
    end do

    ; Resources
    res = True
    if (res_in) then
      copy_VarMeta(res_in, res)
    end if

    ; plot whiskers
    ; (see Example 2 of gsn_add_polyline help page):
    ; Attaching the return value of gsn_add_polyline
    ; to the plot variable ensures that this
    ; reference can live outside this script
    str = unique_string("polyline")
    plot@$str$ = gsn_add_polyline(wks, plot, x, y, res)
  end if

  leave_msg(scriptname, funcname)

end

; #############################################################################
undef("add_prediction_error")
procedure add_prediction_error(wks[1]:graphic,
                               plot[1]:graphic,
                               res_in[1]:logical,
                               xpos_in:numeric,
                               ypos_in:numeric)
;
; Arguments
;    wks: valid workstation, e.g. created by get_wks
;    plot: plot identifier, e.g. created by gsn_*
;    res_in: plot ressources that may override local function defaults
;    xpos_in: array defining the x-position of the four corners
;    ypos_in: vertical marker position(s)
;
; Return value
;    Attaches polytext IDs as attributes to plot.
;
; Description:
;    Adds markers to an existing plot, being the last two intagers of a year.
;    If a horizontal (vertical) coordinate
;    has only one element, then this position is used for all markers.
;
; Caveats:
;
; Modification history
;    20150914-wenzel_sabrina: written.
;
local funcname, scriptname, xpos, ypos, res, str, atts, iatt, color
begin

  funcname = "add_prediction_error"
  scriptname = "diag_scripts/shared/plot/aux_plotting.ncl"
  enter_msg(scriptname, funcname)

  ; Set default plot ressources
  res = True

  ; Create a plot with filled area between both curves
  res@xyLineColor       = "grey85"
  res@gsFillColor       = "grey85"
  res@gsLineDashPattern = 0

  ; Override defaults (after cleaning up for possible vector input)
  if (res_in) then
    atts = getvaratts(res_in)
    if (.not.all(ismissing(atts))) then
      do iatt = 0, dimsizes(atts) - 1
        if (isatt(res, atts(iatt))) then
          delete(res@$atts(iatt)$)
        end if
      end do
    end if
    copy_VarMeta(res_in, res)
  end if

  ; Plot symbol
  str = unique_string("polygon")
  plot@$str$ = gsn_add_polygon(wks, plot, xpos_in, ypos_in, res)

  leave_msg(scriptname, funcname)

end

; #############################################################################
undef("month_sel")
function month_sel(month_names:string)
;
; Arguments
;
;
; Return value
;
;
; Description
;
;
; Caveats
;
; Modification history
;    20190605-righi_mattia: ported from v1.
;
local month_names, month_n, i
begin

  if (.not.any(ismissing(ind(month_names.eq."ANN")))) then
    month_n = new(12, "integer")
    month_n = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)
  else
    month_n = new(dimsizes(month_names), "integer")
    do i = 0, dimsizes(month_names)-1
      if (.not.ismissing(ind(month_names(i).eq."JAN"))) then
        month_n(i) = 1
      end if
      if (.not.ismissing(ind(month_names(i).eq."FEB"))) then
        month_n(i) = 2
      end if
      if (.not.ismissing(ind(month_names(i).eq."MAR"))) then
        month_n(i) = 3
      end if
      if (.not.ismissing(ind(month_names(i).eq."APR"))) then
        month_n(i) = 4
      end if
      if (.not.ismissing(ind(month_names(i).eq."MAY"))) then
        month_n(i) = 5
      end if
      if (.not.ismissing(ind(month_names(i).eq."JUN"))) then
        month_n(i) = 6
      end if
      if (.not.ismissing(ind(month_names(i).eq."JUL"))) then
        month_n(i) = 7
      end if
      if (.not.ismissing(ind(month_names(i).eq."AUG"))) then
        month_n(i) = 8
      end if
      if (.not.ismissing(ind(month_names(i).eq."SEP"))) then
        month_n(i) = 9
      end if
      if (.not.ismissing(ind(month_names(i).eq."OCT"))) then
        month_n(i) = 10
      end if
      if (.not.ismissing(ind(month_names(i).eq."NOV"))) then
        month_n(i) = 11
      end if
      if (.not.ismissing(ind(month_names(i).eq."DEC"))) then
        month_n(i) = 12
      end if
    end do
  end if

  return(month_n)

end

; #############################################################################
undef("lat_names")
function lat_names(lat_val:numeric)
;
; Arguments
;
;
; Return value
;
;
; Description
;
;
; Caveats
;
; Modification history
;    20190605-righi_mattia: ported from v1.
;
local lat_name, lat_val, i
begin

  lat_name = new((/dimsizes(lat_val)/), "string")
  do i = 0, dimsizes(lat_val) - 1
    if (lat_val(i).gt.0.) then
      lat_name(i) = lat_val(i) + "N"
    end if
    if (lat_val(i).eq.0.) then
      lat_name(i) = "Equator"
    end if
    if (lat_val(i).lt.0.) then
      lat_name(i) = -lat_val(i) + "S"
    end if
  end do

  return (lat_name)

end

; #############################################################################
undef("add_line")
procedure add_line(wks,
                   plot,
                   x,
                   yval,
                   line_color,
                   line_dash)
;
; Arguments
;
;
; Description
;
;
; Caveats
;
; Modification history
;    20190605-righi_mattia: ported from v1.
;
local plres, str
begin

  ; gsn_add_polyline crashes if the arrays x and/or yval contain only missing
  if (all(ismissing(x)).or.all(ismissing(yval))) then
    return
  end if

  plres = True
  plres@gsLineColor = line_color
  plres@gsLineDashPattern = line_dash
  plres@gsLineThicknessF = 2.5
  str = unique_string("polyline")
  plot@$str$ = gsn_add_polyline(wks, plot, x, yval, plres)

end

; #############################################################################
undef("add_scatt")
procedure add_scatt(wks,
                    plot,
                    x,
                    yval,
                    mark_color,
                    mark_ind)
;
; Arguments
;
;
; Description
;
;
; Caveats
;
; Modification history
;    20190605-righi_mattia: ported from v1.
;
local plres, str
begin

  plres = True
  plres@gsMarkerColor = mark_color
  plres@gsMarkerIndex = mark_ind
  str = unique_string("polymarker")
  plot@$str$ = gsn_add_polymarker(wks, plot, x, yval, plres)

end

; #############################################################################
undef("add_legend")
procedure add_legend(wks,
                     model:string,
                     colors:string,
                     dashs:numeric,
                     scatters:string,
                     ticknesses:numeric,
                     place:string)
;
; Arguments
;
;
; Return value
;
;
; Description
;
;
; Caveats
;
; Modification history
;    20190605-righi_mattia: ported from v1.
;
local model, dim_mod, colors, dashs, edge_x, edge_y, edge_res, tx1res, \
  dim_lb, xleg, xl, yleg, yl, xtxt, ytxt, place, dim_lines, space_lines, \
  lin_sp, scatters, ticknesses
begin

  dim_mod = dimsizes(model)
  dim_lines = floattoint(dim_mod / 3)+1
  if ((mod(dim_mod, 3)).eq.0) then
    dim_lines = floattoint(dim_mod /3) + 1
  else
    dim_lines = floattoint(dim_mod / 3) + 2
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
  bot_line = top_line-space_lines
  edge_x = (/ 0.12, 0.92, 0.92, 0.12, 0.12/)
  edge_y = (/space_lines+bot_line, space_lines+bot_line, bot_line, \
            bot_line, space_lines+bot_line/)
  edge_res = True
  edge_res@gsLineColor = "white"
  gsn_polyline_ndc(wks, edge_x, edge_y, edge_res)

  dim_lb = dimsizes(model)

  gs1res               = True
  tx1res               = True
  tx1res@txFontHeightF = lin_sp / 5 * 2

  xleg = new((/3 * dim_lines, 4/), "float")
  xl = new((/3, 4/), "float")

  do i = 0, 2
    xl(0, :) = fspan((edge_x(0) + 0.02), (edge_x(0) + 0.07), 4)
    xl(1, :) = fspan(((edge_x(2) + edge_x(0)) / 2 - 0.09), \
                     ((edge_x(2) + edge_x(0)) / 2 - 0.04), 4)
    xl(2, :) = fspan((edge_x(2) - 0.20), (edge_x(2) - 0.15), 4)
  end do

  nt = 0
  do  j = 0, dim_lines - 1
    do i = 0, 2
      xleg(nt, :) = (/xl(i, :)/)
      nt = nt + 1
    end do
  end do

  yleg = new((/3 * dim_lines, 4/), "float")
  yl = new(dim_lines, "float")

  do i = 0, dim_lines - 1
    yl(dim_lines - 1 - i) = edge_y(3) + \
      (((edge_y(0) - (edge_y(3))) / dim_lines) * (i + 1)) - (lin_sp)
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
      xtxt(nt) = xl(i, 3) + lin_sp / 5.
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

  do i = 0, (dimsizes(model) - 1)
    gs1res@gsLineColor = colors(i)
    gs1res@gsLineDashPattern = dashs(i)
    gs1res@gsLineThicknessF = ticknesses(i)
    tx1res@txFontColor = colors(i)
    tx1res@txJust = "CenterLeft"

    if (scatters(i).eq."Markers")
      gs1res@gsMarkerColor = colors(i)
      gs1res@gsMarkerIndex = dashs(i)
      gs1res@gsMarkerSizeF = 0.01
      gs1res@gsMarkerThicknessF = 1.5
      gsn_polymarker_ndc(wks,  xleg(i, 2), yleg(i, 2), gs1res)
    else
      gsn_polyline_ndc(wks,  xleg(i, :), yleg(i, :), gs1res)
    end if
    gsn_text_ndc(wks, model(i), xtxt(i), ytxt(i), tx1res)
  end do

end

; #############################################################################
undef("calcRegCoeffs")
function calcRegCoeffs(y:numeric)
;
; Arguments
;
;
; Return value
;
;
; Description
;
;
; Caveats
;
; Modification history
;    20190605-righi_mattia: ported from v1.
;
local y, d_y, d_t, ma, x, funcx, pi, t, a
begin
  d_y = dimsizes(y)
  d_t = d_y(0)

  ma = 14
  x = fspan(1, d_t, d_t)
  a = new((/ma/), typeof(y))

  funcx = new((/ma, d_t/), "float")
  pi = 4.0 * atan(1.0)
  t = 2 * pi * (x - 0.5) / 12

  funcx(0, :) = 1
  funcx(1, :) = sin(t)
  funcx(2, :) = cos(t)
  funcx(3, :) = sin(2 * t)
  funcx(4, :) = cos(2 * t)
  funcx(5, :) = sin(3 * t)
  funcx(6, :) = cos(3 * t)
  funcx(7, :) = sin(4 * t)
  funcx(8, :) = cos(4 * t)
  funcx(9, :) = x
  funcx(10, :) = sin(t) * x
  funcx(11, :) = cos(t) * x
  funcx(12, :) = sin(2 * t) * x
  funcx(13, :) = cos(2 * t) * x

  a(:) = reg_multlin(y(:), funcx, False)
  a@_FillValue = 1e20
  delete_VarAtts(a,  "constant")

  a@long_name = "Coefficients A and B"
  a!0 = "number of coefficients"

  return(a)

end

; ###########################################################################
undef("genZonalMeans")
function genZonalMeans(x:numeric,
                       lat:numeric,
                       startLat:numeric,
                       endLat:numeric)
;
; Arguments
;
;
; Return value
;
;
; Description
;
;
; Caveats
;
; Modification history
;    20190605-righi_mattia: ported from v1.
;
local dimx, ntim, nlat, nzone, z, t, i, j, cnt, idx_st, idx_ed, xZonMean
begin
  dimx = dimsizes(x)
  ntim = dimx(0)
  nlat = dimx(1)
  nlon = dimx(2)
  nzone = dimsizes(startLat)

  nmos = 12
  modCheck("genZonalMeans", ntim, nmos)

  xZonMean = new((/ntim, nzone/), typeof(x))

  pi  = 4. * atan(1.0)
  rad = (pi / 180.)

  ; loop through lat_band
  do t = 0, (ntim-1)
    do z = 0, (nzone-1)
      idx_st = 0
      do while (lat(idx_st).lt.startLat(z))
        idx_st = idx_st + 1
      end do
      idx_ed = nlat - 1
      do while (lat(idx_ed).gt.endLat(z))
        idx_ed = idx_ed - 1
      end do

      zoneTotalX = 0.0
      ZoneTotalArea = 0.0

      do i = idx_st, idx_ed
        if (i.eq.idx_st) then
          Lat1 = startLat(z)
        else
          Lat1 = (lat(i) + lat(i - 1))/2
        end if
        if (i.eq.idx_ed) then
          Lat2 = endLat(z)
        else
          Lat2 = (lat(i) + lat(i + 1))/2
        end if
        lat_slice = x(t, i, :)
        idx_lon = ind(.not.ismissing(lat_slice))
        if (.not.all(ismissing(idx_lon))) then
          CellArea = abs(sin(Lat1 * rad) - sin(Lat2 * rad))
          zoneTotalX = zoneTotalX + sum(lat_slice) * CellArea
          ZoneTotalArea = ZoneTotalArea + dimsizes(idx_lon) * CellArea
        end if
        delete(idx_lon)
      end do

      if (ZoneTotalArea.ne.0) then
        xZonMean(t, z) = zoneTotalX / ZoneTotalArea
      end if
    end do  ; loop(z)
  end do  ; loop(t)

  ; Create an informational attribute:
  xZonMean@info = "zonal mean: CCMOzoneModelData.ncl"

  return (xZonMean)

end

; #############################################################################
undef("calcMeanAnnCycleMonthly")
function calcMeanAnnCycleMonthly(RegCoeffAs:numeric)
;
; Arguments
;
;
; Return value
;
;
; Description
;
;
; Caveats
;
; Modification history
;    20190605-righi_mattia: ported from v1.
;
local dimx, z, month, t, xMeanAnnCycleMon
begin
  dimx = dimsizes(RegCoeffAs)

  pi  = 4. * atan(1.0)

  month = fspan(1, 12, 12)
  t = 2 * pi * (month - 0.5) / 12.0

  xMeanAnnCycleMon = new((/12/), "float")
  xMeanAnnCycleMon(:) = RegCoeffAs(0) + \
    RegCoeffAs(1) * sin(t) + RegCoeffAs(2) * cos(t) + \
    RegCoeffAs(3) * sin(2 * t) + RegCoeffAs(4) * cos(2 * t) + \
    RegCoeffAs(5) * sin(3 * t) + RegCoeffAs(6) * cos(3 * t) + \
    RegCoeffAs(7) * sin(4 * t) + RegCoeffAs(8) * cos(4 * t)

  ; Create an informational attribute:
  xMeanAnnCycleMon@info = "mean annual cycle monthly: CCMOzoneModelData.ncl"

  return(xMeanAnnCycleMon)

end

; #############################################################################
undef("calcMeanAnnCycleAnnual")
function calcMeanAnnCycleAnnual(RegCoeffAs:numeric)
;
; Arguments
;
;
; Return value
;
;
; Description
;
;
; Caveats
;
; Modification history
;    20190605-righi_mattia: ported from v1.
;
local dimx, day, t, xMeanAnnCycleDay
begin

  dimx = dimsizes(RegCoeffAs)

  pi = 4. * atan(1.0)

  day = fspan(1, 366, 366)
  t = 2 * pi * (day - 0.5) / 366.0

  xMeanAnnCycleDay = new((/366/), "float")

  xMeanAnnCycleDay(:) = RegCoeffAs(0) + \
    RegCoeffAs(1) * sin(t) + RegCoeffAs(2) * cos(t) + \
    RegCoeffAs(3) * sin(2 * t) + RegCoeffAs(4) * cos(2 * t) + \
    RegCoeffAs(5) * sin(3 * t) + RegCoeffAs(6) * cos(3 * t) + \
    RegCoeffAs(7) * sin(4 * t) + RegCoeffAs(8) * cos(4 * t)

  ; Create an informational attribute:
  xMeanAnnCycleDay@info = "mean annual cycle daily: CCMOzoneModelData.ncl"

  return (xMeanAnnCycleDay)

end

; #############################################################################
undef("rmMeanAnnCycle")
procedure rmMeanAnnCycle(x:numeric,
                         xMeanAnnCycle:numeric,
                         RegCoeffA1:numeric,
                         xMonAnom:numeric,
                         xAnnAnom:numeric)
;
; Arguments
;
;
; Return value
;
;
; Description
;
;
; Caveats
;
; Modification history
;    20190605-righi_mattia: ported from v1.
;
local dimx, ntim, t, i, idx_st, idx_ed, nyr, xdata
begin
  dimx = dimsizes(x)
  ntim = dimx(0)

  nmos  = 12
  modCheck("rmMeanAnnCycle", ntim, nmos)  ; error check

  nyr = ntim/nmos

  ; loop through lat_band
  do t = 0, (ntim - 1)
    mon = mod(t, 12)
    xMonAnom(t) = (x(t) - xMeanAnnCycle(mon)) * 100.0 / RegCoeffA1
  end do  ; loop(t)

  do t = 0, (nyr - 1)
    idx_st = t * 12
    idx_ed = (t + 1) * 12 - 1
    xdata = xMonAnom(idx_st:idx_ed)
    if(.not.any(ismissing(xdata))) then
      xAnnAnom(t) = avg(xdata)
    end if
  end do  ; loop(t)

  ; Create an informational attribute:
  xMonAnom@info = "monthly anomalies: CCMOzoneModelData.ncl"
  xAnnAnom@info = "annual anomalies: CCMOzoneModelData.ncl"

end

; #############################################################################
undef("apfiltersmooth")
function apfiltersmooth(xMonAnom:numeric,
                        filter:numeric,
                        iterations:integer)
;
; Arguments
;
;
; Return value
;
;
; Description
;
;
; Caveats
;
; Modification history
;    20190605-righi_mattia: ported from v1.
;
local dimx, dimf, ntim, nwin, z, t, i, j, hfw, idx_ed, avgwidth, \
  totalwgt, total, num, minval, finished, data, hold, xSmthMonAnom
begin
  dimx = dimsizes(xMonAnom)
  ntim = dimx(0)

  nmos = 12
  modCheck("apfiltersmooth", ntim, nmos)  ; error check

  dimf = dimsizes(filter)
  nwin = dimf(0)
  hfw = nwin / 2

  data = new(nwin, typeof(xMonAnom))
  hold = new(ntim, typeof(xMonAnom))
  xSmthMonAnom = new((/ntim/), typeof(xMonAnom))

  do i = 1, iterations

    ; Transfer the data to the storage array.
    if (i.eq.1) then
      hold(:) = xMonAnom(:, z)
    else
      hold(:) = xSmthMonAnom(:, z)
    end if

    ; Fill the data array for the first time.
    do j = 1, hfw
      data(j - 1) = -999.99
    end do
    do j = hfw, (nwin - 1)
      data(j) = hold(j - hfw)
    end do
    idx_ed = nwin-hfw

    ; Process the data array.
    do t = 0, (ntim - 1)
      if (data(hfw).lt.-999.0) then
        xSmthMonAnom(t, z) = -999.99
      else
        avgwidth = 0
        totalwgt = filter(hfw)
        total = filter(hfw) * data(hfw)
        num = 1
        finished = False
        do while (.not.finished)
          avgwidth = avgwidth + 1
          if ((data(hfw - avgwidth).gt.-999.0).and. \
              (data(hfw + avgwidth).gt.-999.0)) then
            totalwgt = totalwgt + filter(hfw - avgwidth) + \
              filter(hfw + avgwidth)
            total = total + (filter(hfw - avgwidth) * data(hfw - avgwidth)) + \
              (filter(hfw + avgwidth) * data(hfw + avgwidth))
            num = num+2
          else
            finished = True
          end if
          if (avgwidth.eq.hfw) then
            finished = True
          end if
        end do
        if (i.eq.iterations) then
          minval = hfw
        else
          minval = 0
        end if
        if (num.gt.minval) then
          xSmthMonAnom(t) = total / totalwgt
        else
          xSmthMonAnom(t) = -999.99
        end if
      end if

      do j = 1, (nwin - 1)
        data(j - 1) = data(j)
      end do
      idx_ed = idx_ed + 1
      if (idx_ed.gt.ntim) then
        data(nwin - 1) = -999.99
      else
        data(nwin - 1) = hold(idx_ed - 1)
      end if
    end do  ; loop(t)
  end do  ; loop(i)

  ; Create an informational attribute:
  xSmthMonAnom@info = "smoothed monthly anomalies: CCMOzoneModelData.ncl"
  xSmthMonAnom@_FillValue = -999.99

  return (xSmthMonAnom)

end

; #############################################################################
undef("smoothAnomalies")
procedure smoothAnomalies(xAnnAnom:numeric,
                          iterations:integer)
;
; Arguments
;
;
; Return value
;
;
; Description
;
;
; Caveats
;
; Modification history
;    20190605-righi_mattia: ported from v1.
;
local dimx, dimf, ntim, nset, nwin, s, t, i, j, hfw, idx_ed, avgwidth, \
  totalwgt, total, num, minval, finished, smthdata, hold, buffer, filter, \
  method
begin

  filter = (/1., 2., 1./)
  method = "EdgeP"

  dimx = dimsizes(xAnnAnom)
  ntim = dimx(0)

  dimf = dimsizes(filter)
  nwin = dimf(0)
  hfw = nwin/2

  buffer = new(nwin, typeof(xAnnAnom))
  smthdata = new(ntim, typeof(xAnnAnom))
  hold = new(ntim, typeof(xAnnAnom))

  do i = 1, iterations
    ; Transfer the data to the storage array.
    if (i.eq.1) then
      hold = (/xAnnAnom/)
    else
      hold = (/smthdata/)
    end if

    ; Fill the data array for the first time.
    do j = 1, hfw
      buffer(j - 1) = -999.99
    end do
    do j = hfw, (nwin - 1)
      buffer(j) = (/hold(j - hfw)/)
    end do
    idx_ed = nwin - hfw

    ; Process the data array.
    do t = 0, (ntim - 1)
      if (buffer(hfw).lt.-999.0) then
        smthdata(t) = -999.99
      else
        avgwidth = 0
        totalwgt = filter(hfw)
        total = filter(hfw) * buffer(hfw)
        num = 1
        finished = False
        do while (.not.finished)
          avgwidth = avgwidth + 1
          if ((buffer(hfw - avgwidth).gt.-999.0) .and. \
              (buffer(hfw + avgwidth).gt.-999.0)) then
            totalwgt = \
              totalwgt + filter(hfw - avgwidth) + filter(hfw + avgwidth)
            total = \
              total + (filter(hfw - avgwidth) * buffer(hfw - avgwidth)) + \
              (filter(hfw + avgwidth) * buffer(hfw + avgwidth))
            num = num + 2
          else
            finished = True
          end if
          if (avgwidth.eq.hfw) then
            finished = True
          end if
        end do

        if (method.eq."EdgeP") then
          if (num.eq.1) then
            if (buffer(hfw + 1).lt.-999.0) then
              smthdata(t) = (buffer(hfw) + buffer(hfw - 1)) / 2.0
            else
              smthdata(t) = (buffer(hfw) + buffer(hfw + 1)) / 2.0
            end if
          else
            smthdata(t) = total / totalwgt
          end if
        else
          if (i.eq.iterations) then
            minval = hfw
          else
            minval = 0
          end if
          if (num.gt.minval) then
            smthdata(t) = total / totalwgt
          else
            smthdata(t) = -999.99
          end if
        end if  ; if(method)
      end if
      do j = 1, (nwin - 1)
        buffer(j-1) = buffer(j)
      end do
      idx_ed = idx_ed+1
      if (idx_ed.gt.ntim) then
        buffer(nwin-1) = -999.99
      else
        buffer(nwin-1) = hold(idx_ed-1)
      end if
    end do  ; loop(t)
  end do  ; loop(i)
  xAnnAnom(:) = (/smthdata(:)/)

  ; Create an informational attribute:
  xAnnAnom@info = "smoothed annual anomalies: CCMOzoneModelData.ncl"
  xAnnAnom@_FillValue = -999.99

end


; #############################################################################
undef("clmMon2clmDayn")
function clmMon2clmDayn(x:numeric,
                        retOrder:integer,
                        opt:integer)
;
; Arguments
;
;
; Return value
;
;
; Description
;    Take a monthly climatology and make a daily climatology
;
;
; Caveats
;
; Modification history
;    20190605-righi_mattia: ported from v1.
;
local dNames, dimx, rank, X, midMon, day
begin
  if (.not.(retOrder.eq.0 .or. retOrder.eq.1)) then
    print("clmMon2clmDay: retOrder must be 0 or 1,  retOrder=" + retOrder)
    exit
  end if

  dNames = getvardims(x)
  if (any(ismissing(dNames))) then
    print("clmMon2clmDay: named dimensions required:" + dNames)
    exit
  end if

  dimx = dimsizes(x)
  if (dimx(0).ne.12) then
    print("clmMon2clmDay: leftmost dimension must be size=12: SIZE=" + dimx(0))
    exit
  end if

  rank = dimsizes(dimx)
  if (rank.gt.4) then
    print("clmMon2clmDay: currently max of 4 dimension supported: rank=" + \
          rank)
    exit
  end if

  ; transfer to work array, if necessary, reorder array
  if (rank.eq.1) then
    X = x
  end if
  if (rank.eq.2) then
    X = x($dNames(1)$|:, $dNames(0)$|:)
  end if
  if (rank.eq.3)
    X = x($dNames(1)$|:, $dNames(2)$|:, $dNames(0)$|:)
  end if
  if (rank.eq.4)
    X = x($dNames(1)$|:, $dNames(2)$|:, $dNames(3)$|:, $dNames(0)$|:)
  end if

  ; mid day of each month
  if (isatt(opt, "midmon")) then
    if (dimsizes(opt@midMon).eq.12) then
      midMon = opt@midMon
    else
      print("clmMon2clmDay: midMon required to be size 12: size=" + \
            dimsizes(opt@midMon))
      exit
    end if
  else
    midMon = (/15.5, 45, 74.5, 105, 135.5, 166, 196.5, \
              227.5, 258, 288.5, 319, 349.5/)
  end if
  midMon@long_name = "middle of month"

  day    = ispan(0, 364, 1)  ; use 0 => 364 for interpolation
  day!0  = "day"

  Z      = linint1_Wrap(midMon, X, True, day, 0)
  Z@info = "NCL: clmMon2clmDay"

  day    = ispan(1, 365, 1)  ; use 1 => 365 for coord variable
  day@long_name = "day of year: no leap"
  day@units     = "1=Jan 1,  32=Feb 1,  ...,  365-Dec 31"
  Z!(rank-1) = "day"
  Z&day      =  day

  if (retOrder.eq.1) then
    return(Z)
  end if

  if (retOrder.eq.0) then
    if (rank.eq.1) then
      return(Z)
    end if
    if (rank.eq.2) then
      return(Z(day|:, $dNames(1)$|:))
    end if
    if (rank.eq.3) then
      return(Z(day|:, $dNames(1)$|:, $dNames(2)$|:))
    end if
    if (rank.eq.4) then
      return(Z(day|:, $dNames(1)$|:, $dNames(2)$|:, $dNames(3)$|:))
    end if
  end if

end
