; #############################################################################
; PLOT SCRIPTS FOR CONTOUR PLOTS
; #############################################################################
; Please consider using or extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
;    contour: wrapper for gsn_csm_contour
;
; #############################################################################

load "$diag_scripts/../interface_scripts/auxiliary.ncl"
load "$diag_scripts/../interface_scripts/logging.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"

; #############################################################################
undef("contourplot")
function contourplot(wks_in[1],
                     source,
                     varname[1] : string)
;
; Arguments
;    wks    : workstation (graphic object or default will be used).
;    source : data to be plotted or a NetCDF filename with data.
;    varname: variable name in the file.
;
; Source prototype
;    source[*,*]
;
; Return value
;    A graphic variable.
;
; Description
;    Wrapper for gsn_csm_contour (http://www.ncl.ucar.edu/Document/Graphics/
;    Interfaces/gsn_csm_contour.shtml)
;    Opens default wks, if not provided as argument of type "graphic".
;    Defines default resources, which are overridden by argument res.
;    Creates plot, according to wks & res.
;
; Caveats
;    Selection of defaults for res almost arbitrary.
;    Please check results of all scripts that use this routine if modifying
;    the defaults!
;    Input via netCDF not yet implemented.
;
; Modification history
;    20131217-vanulft_bert: written, largely copied from countour_map.ncl.
;
local data, defaults, diag_script, funcname, res, res_new, scriptname, \
  source, varname, wks, wks_in
begin

  funcname = "contourplot"
  scriptname = "diag_scripts/shared/plot/contourplot.ncl"
  enter_msg(scriptname, funcname)

  ; Get data, either directly or via netCDF file
  if(typeof(source) .eq. "string") then
    data = ncdf_read(source, varname)
  else
    data = source
    copy_VarMeta(source, data)
  end if

  ; Retrieve basic metadata from data
  defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
  if (any(varname .eq. defaults)) then
    ; Use default
    var = att2var(data, "var")
  else
    ; Use explicitly given name
    var = varname
  end if
  diag_script = att2var(data, "diag_script")

  ; Check if a valid wks has been provided, otherwise invoke default
  wks = get_wks(wks_in, diag_script, var)

  ; Set default resources
  res = True
  res@gsnDraw            = True                        ; Draw
  res@gsnFrame           = True                        ; Advance frame
  res@cnFillPalette      = "ncl_default"               ; Set colormap
  res@tiMainString       = "Output of " + diag_script  ; Main title
  res@gsnLeftString      = ""                          ; Annotation
  res@trYReverse         = False                       ; reverse the Y-axis?
  if (isatt(var, "long_name")) then
    res@gsnRightString = var@long_name
  else
    res@gsnRightString = var
  end if
  if (isatt(var, "units")) then
    res@gsnRightString = res@gsnRightString + " in " + var@units
  end if

  ; Override defaults with "res_" attributes of "data"
  res_new = att2var(data, "res_")
  copy_VarMeta(res_new, res)

  ; Create plot
  plot = gsn_csm_contour(wks, data, res)

  if (isatt(wks, "fullname")) then
    plot@outfile = wks@fullname
  else
    plot@outfile = wks@name
    error_msg("w", "contourplot.ncl", funcname, "wks@fullname " + \
              "missing, consider to use get_wks to open wks")
  end if

  if (isatt(res, "gsnDraw")) then
    if (res@gsnDraw .eq. True) then
      log_info(" Wrote " + plot@outfile)
    end if
  else
    log_info(" Wrote " + plot@outfile)
  end if

  leave_msg(scriptname, funcname)
  return(plot)

end
