;;#############################################################################
;; GENERAL ROUTINES FOR UNIT CONVERSION
;;#############################################################################
;; Please consider using of extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    function convert_units
;;    function scalfac
;;    function scale_units
;;    function rescale
;;
;;#############################################################################

load "./interface_scripts/logging.ncl"

;;#############################################################################
undef("convert_units")
function convert_units(var:numeric,
                       units_to[1]:string)
;;
;; Arguments
;;    var: a numeric field on any dimensionality, it must have a units
;;         attribute.
;;    units_to: a string, specifying the units to be converted to.
;;
;; Return value
;;    A numeric field of the same dimensionality of var.
;;
;; Description
;;    Converts units of var from the one specified by the units attribute to
;;    the one given by the units_to argument. An error message is issued if
;;    the requested unit conversion is not available.
;;
;; Caveats
;;    This function shall be used to perform unit conversion before plotting.
;;    For unit conversion to be performed before data processing within the
;;    diagnostic (e.g., non-linear conversions), it is recommended to use
;     variable_defs.
;;
;; References
;;
;; Modification history:
;;    20150216-A_righ_ma: written.
;;
local funcname, scriptname
begin

    funcname = "convert_units"
    scriptname = "diag_scripts/shared/scaling.ncl"
    enter_msg(scriptname, funcname)

    ;; Define output (copy metadata)
    out = var

    ;; Check for units attribute
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

    error_msg("f", scriptname, funcname, "conversion from " + units_from + \
              " to " + units_to + " not defined")

end

;;#############################################################################
undef("scalfac")
function scalfac(u:numeric,
                 digits[1]:numeric)
;;
;; Arguments
;;    u: numeric field of any dimension containing the data to be scaled.
;;    digits: number of desired relevant digits left of the decimal point.
;;
;; Return value
;;    Scaling factor of type float, which would need to be multiplied to the
;;    original data to bring them in the desired range (needed e.g. for
;;    annotation).
;;
;; Description
;;    Calculates a factor for scaling the input data to a range with the
;;    desired number of relevant digits left of the decimal point.
;;    The sign is counted as one of those digits if required for negative
;;    values.
;;
;; Caveats
;;
;; References
;;
;; Modification history:
;;    20140220-A_gott_kl: written.
;;
local funcname, scriptname, result, u, digits, d, log10_max_abs_u
begin

    funcname = "scalfac"
    scriptname = "diag_scripts/shared/scaling.ncl"
    enter_msg(scriptname, funcname)

    ;; consider minus sign for number of digits
    if(abs(min(u)) .gt. 0.1 * max(u)) then
        d = toint(digits) - 1
    else
        d = toint(digits)
    end if

    ;; digit gymnastics
    max_abs_u = max(abs(u))
    if(max_abs_u.gt.0.) then
        log10_max_abs_u = floor(log10(max_abs_u)) - 1
        result = 10 ^ (-log10_max_abs_u + d - 2)
    else
        result = 1.
    end if
    log_debug(" digits=" + digits + " ; " + min(u) + " < data < " + \
              max(u) + " ; " + result * min(u) + " < data*" + result + \
              " < " + result * max(u))

    leave_msg(scriptname, funcname)
    return(result)

end

;;#############################################################################
undef("scale_units")
function scale_units(u:numeric,
                     digits[1]:numeric)
;;
;; Arguments
;;    u: numeric field of any dimension containing the data to be scaled.
;;    digits: number of desired relevant digits left of the decimal point
;;            (needed e.g. for annotation).
;;
;; Return value
;;    Scaling factor of type float, which would need to be multiplied to the
;;    original data to bring them in the desired range.
;;    The string of the new unit is added as attribute.
;;    If the "units" string contains function codes, the code used will be
;;    added in s@FuncCode.
;;
;; Description
;;    Calculates a factor for scaling the input data to a range with the
;;    desired number of relevant digits left of the decimal point.
;;    The result is further scaled to match a "nice" unit (if available)
;;    depending on the original unit attribute of u.
;;    Format codes: http://www.ncl.ucar.edu/Applications/text.shtml
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20140220-A_gott_kl: written.
;;
local funcname, scriptname, u, digits, s
begin

    funcname = "scale_units"
    scriptname = "diag_scripts/shared/scaling.ncl"
    enter_msg(scriptname, funcname)

    if(.not.all(ismissing(u))) then
        s = scalfac(u, digits)  ; function in scaling.ncl
    else
        error_msg("w", scriptname, funcname, "all input values missing, " + \
                  "no rescaling in scale_units")
        s = 1.
    end if

    if(isatt(u, "units")) then
        ;; mole mole-1 ########################################################
        if (u@units.eq."mole mole-1") then
            if(s.ge.0.5e5 .and. s.le.5.e7) then
                s = 1.e6
    ;; s@units = ":V1::F33:m:F0:mol:H-90::V-6::F21:______:H-80::V-30::F0:mol"
                s@FuncCode = ":"
                s@units = ":F33:m:F:mol/mol"  ; ppmv
            else if(s.gt.5.e7 .and. s.le.5.e10) then
                s = 1.e9
    ;; s@units = ":F0::V1:nmol:H-90::V-6::F21:______:H-80::V-30::F0:mol"
                ;; s@FuncCode = ":"
                s@units = "nmol/mol"  ; ppbv
            else if(s.gt.5.e10 .and. s.le.5.e13) then
                s = 1.e12
    ;; s@units = ":F0::V1:pmol:H-90::V-6::F21:______:H-80::V-30::F0:mol"
                ;; s@FuncCode = ":"
                s@units = "pmol/mol"  ; pptv
            end if
            end if
            end if
        end if  ; mole mole-1

        ;; DU #################################################################
        if (u@units.eq."DU") then  ; leave as is
            s = 1.
            sunits = "DU"
        end if  ; DU

        ;; Specific Humidity ##################################################
        if(u@long_name.eq."Specific Humidity") then  ; units is "1"
            if(s.le.50.) then
                s = 1.
                s@units = "kg/kg"
            else if(s.gt.50. .and. s.le.50000.) then
                s = 1.e3
                s@units = "g/kg"
            else if(s.gt.50000.) then
                s = 1.e6
                s@units = "mg/kg"
            end if
            end if
            end if
        end if  ; Specific Humidity

        ;; insert other units here ############################################

        ;; Default ############################################################
        if(.not.isatt(s, "units")) then
            if(s.eq.1.) then
                s@units = u@units
            else
                s@units = tostring(s) + " * " + u@units
            end if
        end if  ; Default

    else

        if(s.ne.1.) then
            s@units = "/" + tostring(s)
        else
            s@units = "1"
        end if

    end if

    log_info(" scalfac=" + s + "  ;  units: " + s@units)

    leave_msg(scriptname, funcname)
    return(s)

end

;;#############################################################################
undef("rescale")
function rescale(u:numeric,
                 s[1]:numeric)
;;
;; Arguments
;;    u: numeric field of any dimension containing the data to be scaled.
;;    s: scaling factor, e.g. as output by function scale_units.
;;
;; Return value
;;    Field of the same dimensions and type of u, rescaled to the units
;;    of s.
;;
;; Description
;;    Rescales the input field according to scaling factor s.
;;    Metadata of u are kept, but the units attribute will be changed according
;;    to s@units.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20140221-A_gott_kl: written.
;;
local funcname, scriptname, u, s, result
begin

    funcname = "rescale"
    scriptname = "diag_scripts/shared/scaling.ncl"
    enter_msg(scriptname, funcname)

    result = u * s
    copy_VarMeta(u, result)
    if(isatt(s, "units")) then
        result@units = s@units
    end if
    if(isatt(s, "FuncCode")) then
        result@FuncCode = s@FuncCode
    end if

    leave_msg(scriptname, funcname)
    return(result)

end
