; #############################################################################
; ROUTINES FOR CREATING LEGENDS AS SEPARATE PLOTS
; #############################################################################
; Please consider using or extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
; Contents:
;
;    function legend_lines
;    function legend_markers
;
; #############################################################################

load "$diag_scripts/../interface_scripts/constants.ncl"
load "$diag_scripts/../interface_scripts/logging.ncl"

; #############################################################################
undef("legend_lines")
function legend_lines(wks_in[1],
                      source,
                      varname[1]:string)
;
; Arguments
;     wks_in: workstation ("graphic" object or default will be used)
;     source: data to be plotted, either passed directly or via netCDF file
;             legend styles and strings are passed as attributes data@
;     varname: variable name, needed for netCDF files with multiple variables
;
;     The following attributes of the input data are evaluated:
;     source@diag_script: name(s) of the calling script(s)
;     source@colors: vector of colors
;     source@thicks: vector of line thicknesses
;     source@dashes: vector of dash styles
;     source@annots: vector of label strings
;     source@nrow: number of rows (optional)
;     source@ncol: number of columns (optional)
;
; Return value
;    wks: workstation with legend
;
; Description
;    Creates an extra plot with a legend, specified by labels and line styles.
;    It will be saved in the outfile directory and returned as a workstation
;
; Modification history
;    20140326-gottschaldt_klaus-dirk: written based on code by frank_franziska
;
local data, defaults, var, diag_script, colors, annots, dashes, thicks
begin

  funcname = "legend_lines"
  scriptname = "diag_scripts/shared/plot/legends.ncl"
  enter_msg(scriptname, funcname)

  ; Get data, either directly or via netCDF file
  if (typeof(source) .eq. "string") then
    ; Function in ~/interface_scripts/auxiliary.ncl
    data = ncdf_read(source, varname)  ; can deal with defaults
  else
    data = source
    copy_VarMeta(source, data)
  end if

  ; Retrieve basic metadata from data
  defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
  if(any(varname .eq. defaults)) then  ; Use default
    var = att2var(data, "var")      ; interface_scripts/auxiliary.ncl
  else                                ; Use explicitely given name
    var = varname
  end if
  diag_script = att2var(data, "diag_script")
  colors = att2var(data, "colors")
  thicks = att2var(data, "thicks")
  dashes = att2var(data, "dashes")
  annots = att2var(data, "annots")

  ; Attributes for polyline and text
  gsres = True
  txres = True
  copy_VarAtt_sel(data, txres, (/"tx"/))
  ; txres@txAngleF = -90.0

  ; Determine dimensions of plot
  dim_LAB = dimsizes(annots)

  ; Columns & rows
  if (isatt(data, "ncol") .and. isatt(data, "nrow")) then
    col = data@ncol
    row = data@nrow
  else if (isatt(data, "ncol") .and. .not.isatt(data, "nrow")) then
    col = data@ncol
    row = dim_LAB / col
    if (row.lt.tofloat(dim_LAB) / tofloat(col)) then
      row = row + 1
    end if
  else if(.not. isatt(data, "ncol") .and. isatt(data, "nrow")) then
    row = data@nrow
    col = dim_LAB / row
    if (col.lt.tofloat(dim_LAB) / tofloat(row)) then
      col = col + 1
    end if
  else  ; default
    col = 3
    row = dim_LAB / col
    if (row.lt.tofloat(dim_LAB) / tofloat(col)) then
      row = row + 1
    end if
  end if
  end if
  end if

  ; Determines anchor points of lines and labels
  xleg = new((/col/), float)
  yleg = new((/row/), float)

  xleg(0) = 1. / (1.1 * col)
  do j = 1, col - 1
    xleg(j) = xleg(j - 1) + 1. / col
  end do

  yleg(0) = 1 - 1. / (2 * row)
  do j = 1, row - 1
    yleg(j) = yleg(j - 1) - 1. / row
  end do

    ; Attributes of workstation
;    file_type = True
;    file_type@wkOrientation = "portrait"
;    file_type@wkPaperHeightF = 4. * col
;    file_type@wkPaperWidthF = 1.5 * row
;    file_type = "default"
;    absFontH = 0.3 ; absolute Font Height
;    txres@txFontHeightF = absFontH / file_type@wkPaperWidthF
;    txres@txFontHeightF = absFontH / (1.5 * row)
;    txres@txFontAspectF = 12. * (1. / row)

  ; Check if a valid wks has been provided, otherwise invoke default
  ; attributes of file_type will be used by get_wks
  wks_legend = get_wks(wks_in, diag_script, var+"_legend")

  ; Loop over rows (nrow) and columns (col) to draw lines and text
  i = 0
  do r = 0, row - 1
    do c = 0, col - 1
      if i.lt.(dim_LAB) then

        ; Set line color and other optional styles
        gsres@gsLineColor = colors(i)
        gsres@gsLineDashPattern = dashes(i)
        gsres@gsLineThicknessF = thicks(i)

        ; Draw line and add text
        gsn_polyline_ndc(wks_legend, (/yleg(r), yleg(r)/), \
                         (/xleg(c), xleg(c) - 0.05/), gsres)
        gsn_text_ndc(wks_legend, annots(i), yleg(r), \
                     xleg(c)-0.17, txres)
        i = i + 1
      end if
    end do
  end do

  ; Advance frame
  frame(wks_legend)

  ; Outfile name
;    if (isatt(wks_legend, "fullname")) then
;        plot@outfile = wks_legend@fullname
;    else
;        plot@outfile = wks_legend@name
;        log_debug("warning: wks@fullname missing, consider to use " +\
;                    "get_wks to open wks")
;    end if
;    log_info(" Wrote " + plot@outfile)
  log_info(" Wrote " + wks_legend@fullname)

  leave_msg(scriptname, funcname)
  return(wks_legend)

end

; #############################################################################
undef("legend_markers")
function legend_markers(wks_in[1],
                        source,
                        varname[1]:string)
;
; Arguments
;     wks: workstation ("graphic" object or default will be used)
;     source: data to be plotted, either passed directly or via netCDF file
;             legend styles and strings are passed as attributes data@
;     varname: variable name, needed for netCDF files with multiple variables
;
;     The following attributes of the input data are evaluated:
;     source@diag_script: name(s) of the calling script(s)
;     source@colors: vector of colors
;     source@annots: vector of label strings
;     source@nrow: number of rows (optional)
;     source@ncol: number of columns (optional)
;
; Return value
;    wks: workstation with legend
;
; Description
;    Creates an extra plot with a legend, specified by labels and line styles.
;    It will be saved in the outfile directory and returned as a workstation
;
; Modification history
;    20140326-wenzel_sabrina: written based on legend_lines.
;
local data, defaults, var, diag_script, colors, annots, dashes, thicks
begin

  funcname = "legend_markers"
  scriptname = "diag_scripts/shared/plot/legends.ncl"
  enter_msg(scriptname, funcname)

  ; Get data, either directly or via netCDF file
  if (typeof(source) .eq. "string") then
    ; Function in ~/interface_scripts/auxiliary.ncl
    data = ncdf_read(source, varname)  ; can deal with defaults
  else
    data = source
    copy_VarMeta(source, data)
  end if

  ; Retrieve basic metadata from data
  defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
  if(any(varname .eq. defaults)) then
    var = att2var(data, "var")
  else
    var = varname
  end if
  diag_script = att2var(data, "diag_script")
  colors = att2var(data, "colors")
  thicks = att2var(data, "thicks")
  annots = att2var(data, "annots")

  ; Attributes for polyline and text
  gsres = True
  txres = True
  copy_VarAtt_sel(data, txres, (/"tx"/))
  ; txres@txAngleF = -90.0

  ; Determine dimensions of plot
  dim_LAB = dimsizes(annots)

  ; Columns & rows
  if (isatt(data, "ncol") .and. isatt(data, "nrow")) then
    col = data@ncol
    row = data@nrow
  else if (isatt(data, "ncol") .and. .not.isatt(data, "nrow")) then
    col = data@ncol
    row = dim_LAB / col
    if (row.lt.tofloat(dim_LAB) / tofloat(col)) then
      row = row + 1
    end if
  else if(.not. isatt(data, "ncol") .and. isatt(data, "nrow")) then
    row = data@nrow
    col = dim_LAB / row
    if (col.lt.tofloat(dim_LAB) / tofloat(row)) then
      col = col + 1
    end if
  else  ; default
    col = 3
    row = dim_LAB / col
    if (row.lt.tofloat(dim_LAB) / tofloat(col)) then
      row = row + 1
    end if
  end if
  end if
  end if

  ; Determines anchor points of lines and labels
  xleg = new((/col/), float)
  yleg = new((/row/), float)

  xleg(0) = 1. / (1.1 * col)
  do j = 1, col - 1
    xleg(j) = xleg(j - 1) + 1. / col
  end do

  yleg(0) = 1 - 1. / (2 * row)
  do j = 1, row - 1
    yleg(j) = yleg(j - 1) - 1. / row
  end do

  ; Check if a valid wks has been provided, otherwise invoke default
  ; attributes of file_type will be used by get_wks
  wks_legend = get_wks(wks_in, diag_script, var+"_legend")

  ; Loop over rows (nrow) and columns (col) to draw lines and text
  i = 0
  do r = 0, row - 1
    do c = 0, col - 1
      if i.lt.(dim_LAB) then

        ; Set marker color and other optional styles
        gsres@gsMarkerColor      = colors(i)
        gsres@gsMarkerIndex      = thicks(i)
        gsres@gsMarkerSizeF      = 0.015
        gsres@gsMarkerThicknessF = 2

        ; Draw marker and add text
        gsn_polymarker_ndc(wks_legend, (/yleg(r), yleg(r)/), \
                           (/xleg(c), xleg(c) - 0.05/), gsres)
        gsn_text_ndc(wks_legend, annots(i), yleg(r), \
                     xleg(c)-0.17, txres)
        i = i + 1
      end if
    end do
  end do

  ; Advance frame
  frame(wks_legend)

  ; log_info(" Wrote " + wks_legend@fullname)

  leave_msg(scriptname, funcname)
  return(wks_legend)

end
