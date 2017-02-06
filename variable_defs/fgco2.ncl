;
;  Requires: fgco2:*2*s
;
variable_info = True
variable_info@derived = True 
variable_info@units = "kg m-2 s-1"
variable_info@plot_units = "PgC y-1"

;;#############################################################################
load "interface_scripts/data_handling.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef("calculate")
function calculate(index[1]:integer,
                   variable[1]:string,
                   field_type[1]:string)
;;
;; Arguments
;;    index: index to current infile defined in the 
;;           'interface_data/ncl.interface'-file
;;    variable: a string with current variable
;;    field_type: a string with field type classification
;;
;; Return value
;;    data_new: logical
;;
;; Description
;;    Corrects fgco2 for coastal regions. If available, the ocean mask for the
;;    given model is used, otherwise the NCL mask is applied (for regular grids
;;    only.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20151029-A_righ_ma: written based on a previous version by wenz_sa.
;;
local funcname, scriptname, verbosity, result, tmp, maskok, f, a, sftof, \
    FGCO2, w, areacello, lat, lon
begin

    funcname = "calculate"
    scriptname = "variable_defs/FGCO2.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    ;; Read original variable
    result = read_data(index, "fgco2", "*2*s")
    tmp = extract_data_extended(index, "fgco2" , result, -1, 0, 0)

    maskok = False
    f = addfile(interface_get_fullpath("fgco2", field_types(0), index), "r")
    if (isfilevar(f, "sftof")) then
        sftof = f->sftof
        sftof = sftof / 100.  ; [%] --> [1]
        maskok = True
    else
        if (iscoord(tmp, "lat") .and. iscoord(tmp, "lon")) then
            if (dimsizes(dimsizes(tmp&lat)).eq.1 .and. \
                dimsizes(dimsizes(tmp&lon)).eq.1) then
                f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
                a = f->LSMASK
                sftof = byte2flt(landsea_mask(a, tmp&lat, tmp&lon))
                sftof = where(sftof.eq.0., 1., 0.)
                info_output("Ocean-mask file not available, using NCL " + \ 
                            "default mask", verbosity, 1)
                maskok = True
            end if
        end if
    end if
         
    if (.not.maskok) then
        error_msg("f", scriptname, funcname, "ocean-mask (sftof) not " + \
                  "available for this model, cannot calculate FGCO2")
    end if
    sftof = where(sftof.eq.0, sftof@_FillValue, sftof)

    ;; Apply mask
    FGCO2 = tmp
    FGCO2 = FGCO2 * conform(tmp, sftof, (/1, 2/))

    data_new = True
    add_data_var(index, data_new, FGCO2, "FGCO2")

    ;; Copy over irregular grid info if available
    if (isfilevar(f, "areacello")) then
        areacello = f->areacello
        lat = f->lat
        lon = f->lon
        w=addfile(interface_get_fullpath("FGCO2", field_types(0), index), "w")
        w->lat = lat
        w->lon = lon
        w->areacello = areacello
        delete(w)
    end if

    leave_msg(scriptname, funcname, 4)
    return(data_new)

end
