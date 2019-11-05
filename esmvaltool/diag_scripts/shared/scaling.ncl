; #############################################################################
; GENERAL ROUTINES FOR UNIT CONVERSION
; #############################################################################
; Please consider using of extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
; Contents:
;    function convert_units
;
; #############################################################################

load "$diag_scripts/../interface_scripts/logging.ncl"

; #############################################################################
undef("convert_units")
function convert_units(var:numeric,
                       units_to[1]:string)
;
; Arguments
;    var: a numeric field on any dimensionality, it must have a units
;         attribute.
;    units_to: a string, specifying the units to be converted to.
;
; Return value
;    A numeric field of the same dimensionality of var.
;
; Description
;    Converts units of var from the one specified by the units attribute to
;    the one given by the units_to argument. An error message is issued if
;    the requested unit conversion is not available.
;
; Caveats
;    This function shall be used to perform unit conversion before plotting.
;    For unit conversion to be performed before data processing within the
;    diagnostic (e.g., non-linear conversions), it is recommended to use
;     variable_defs.
;
; References
;
; Modification history
;    20150216-righi_mattia: written.
;
local funcname, scriptname
begin

  funcname = "convert_units"
  scriptname = "diag_scripts/shared/scaling.ncl"
  enter_msg(scriptname, funcname)

  ; Define output (copy metadata)
  out = var

    ; Check for units attribute
  if (.not.isatt(var, "units")) then
    error_msg("f", scriptname, funcname, "The input variable must " + \
              "have a units attribute")
  else
    units_from = var@units
  end if

  if (units_from.eq."1") then
    if (units_to.eq."g/kg") then
      out = out * 1000.
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
    if (units_to.eq."ppmv") then
      out = out * 1.e-6
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
    if (units_to.eq."ppm9") then
      out = out * 1.e-9
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
  end if

  if (units_from.eq."1e-6") then
    if (any(units_to.eq.(/"ppm", "ppmv", "umol/mol"/))) then
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
  end if

  if (units_from.eq."1e-9") then
    if (any(units_to.eq.(/"ppm", "ppmv", "umol/mol"/))) then
      out = out * 1.e-3
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
    if (any(units_to.eq.(/"ppb", "ppbv", "ppmb", "nmol/mol"/))) then
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
    if (any(units_to.eq.(/"ppt", "pptv", "pmol/mol"/))) then
      out = out * 1.e3
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
  end if

  if (units_from.eq."mole mole-1") then
    if (any(units_to.eq.(/"ppm", "ppmv", "umol/mol"/))) then
      out = out * 1.e6
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
    if (any(units_to.eq.(/"ppb", "ppbv", "nmol/mol"/))) then
      out = out * 1.e9
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
    if (any(units_to.eq.(/"ppt", "pptv", "pmol/mol"/))) then
      out = out * 1.e12
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
  end if

  if (units_from.eq."kg kg-1") then
    if (any(units_to.eq.(/"ng kg-1", "ng/kg"/))) then
      out = out * 1.e12
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
  end if

  if (units_from.eq."kg m-2 s-1") then
    if (any(units_to.eq.(/"PgC y-1", "GtC y-1"/))) then
      out = out * 3600. * 24. * 365 / 1.e12
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
    if (any(units_to.eq.(/"mm d-1", "mm/day"/))) then
      out = out * 24 * 3600
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
  end if

  if (units_from.eq."kg m-2") then
    if (any(units_to.eq.(/"PgC", "GtC"/))) then
      out = out * 1000 / 1.e15
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
    if (units_to.eq."mm") then  ; change only the label
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
  end if

  if (units_from.eq."mol m-2 s-1") then
    if (any(units_to.eq.(/"PgC/y"/))) then
      out = out * 60 * 60 * 24 * 365 * 1000 * 12. * 1.e-3 / 1.e15
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
    if (any(units_to.eq.(/"kg m-2 s-1"/))) then
      out = out * 12. * 1.e-3
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
  end if

  if (units_from.eq."m-3") then
    if (units_to.eq."cm-3") then
      out = out * 1.e-6
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
  end if

  if (units_from.eq."m-3-STP") then
    if (units_to.eq."cm-3-STP") then
      out = out * 1.e-6
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
  end if

  if (units_from.eq."kg m-3") then
    if (any(units_to.eq.(/"ug m-3", "ug/m3"/))) then
      out = out * 1.e9
      out@units = str_sub_str(units_to, "u", "~F33~m~F21~")
      leave_msg(scriptname, funcname)
      return(out)
    end if
  end if

  if (units_from.eq."K") then
    if (units_to.eq."degC") then
      out = out - 273.15
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
  end if

  if (units_from.eq."W m-2") then
    if (units_to.eq."-W m-2") then  ; change only the sign, not the label
      out = -1. * out
      leave_msg(scriptname, funcname)
      return(out)
    end if
    if (units_to.eq."mm d-1") then
      out = out / 2.5e6 * 24 * 3600.
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
  end if

  if (units_from.eq."kg m-2 s-1") then
    if (units_to.eq."mm d-1") then
      out = out * 24 * 3600.
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
  end if

  if (units_from.eq."Pa") then
    if (units_to.eq."hPa") then
      out = out / 100.
      out@units = units_to
      leave_msg(scriptname, funcname)
      return(out)
    end if
  end if

  error_msg("f", scriptname, funcname, "conversion from " + units_from + \
            " to " + units_to + " not defined")

end
