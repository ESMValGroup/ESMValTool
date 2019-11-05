; #############################################################################
; PLOT SCRIPTS FOR ZONAL MEAN PROFILES
; #############################################################################
; Please consider using or extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
;    function zonalmean_profile
;
; #############################################################################

load "$diag_scripts/shared/set_operators.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"

; #############################################################################
undef("zonalmean_profile")
function zonalmean_profile(wks_in[1],
                           source,
                           varname[1]: string)
;
; Arguments
;    wks_in: workstations (graphic object or default will be used).
;    source: data to be plotted or a NetCDF filename with data.
;    varname: variable name in the file.
;
; Source prototype
;    source[*,*]
;    source!0 = plev
;    source!1 = lat
;
; Return value
;    A graphic variable
;
; Description
;    Draws a pressure-latitude plot
;
; Caveats
;
; Modification history
;    20131210-frank_franziska: written.
;
local funcname, scriptname, wks, wks_in, data, source, res, atts, base, varname
begin

  funcname = "zonalmean_profile"
  scriptname = "diag_scripts/shared/plot/zonalmean_profile.ncl"
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
  if (any(varname .eq.  defaults)) then
    var = att2var(data, "var")
    DIAG_SCRIPT = att2var(data, "diag_script")
  else
    var = varname
  end if

  ; Check if a valid wks has been provided, otherwise invoke default
  wks = get_wks(wks_in, DIAG_SCRIPT, var)

  ; Define default plot resources
  res = True
  res@gsnDraw          = False
  res@gsnFrame         = False
  res@cnLineLabelsOn   = False
  res@cnLinesOn        = True
  res@cnFillOn         = True

  ; Set left label
  if (isatt(data, "long_name")) then
    res@gsnLeftString = data@long_name
  else
    res@gsnLeftString = var
  end if
  if (isatt(data, "units")) then
    res@gsnLeftString = \
      res@gsnLeftString + " [" + format_units(data@units) + "]"
  end if

  ; Set right label
  if (isatt(data, "time_avg")) then
    res@gsnRightString = data@time_avg
  else
    res@gsnRightString = ""
  end if

  res@lbLabelAutoStride    = True
  res@tmLabelAutoStride    = False
  res@tmYRMode             = "Automatic"
  res@tmXTOn               = False
  res@tiYAxisString        = "Pressure [hPa]"

  ; Override defaults with "res_" attributes of "data"
  res_new = att2var(data, "res_")
  copy_VarMeta(res_new, res)

  ; Create plot
  plot = gsn_csm_pres_hgt(wks, data, res)
  plot@outfile = wks@name

  leave_msg(scriptname, funcname)
  return(plot)

end
