;;#############################################################################
; CONTOURS AND VECTORS ON A POLARSTEREOGRAPHIC MAP #
;;#############################################################################
;; Please consider using or extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;;    function vector_scalar_map_polar: wrapper for
;;                                      gsn_csm_vector_scalar_map_polar
;;
;;#############################################################################

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/logging.ncl"

load "./diag_scripts/shared/plot/aux_plotting.ncl"

;;#############################################################################
undef("vector_scalar_map_polar")
function vector_scalar_map_polar(wks_in[1],
                                 sourceu,
                                 sourcev,
                                 sourcec,
                                 res : logical,
                                 varnameu[1] : string,
                                 varnamev[1] : string,
                                 varnamec[1] : string)
;;
;; Arguments
;;    wks     : workstation, must be passed - no default used.
;;    sourceu : u vector data to be plotted.
;;    sourcev : v vector data to be plotted.
;;    sourcec : contour data to be plotted.
;;    res     : diag_script-specific resources passed from diag_script.
;;    varnameu: u vector variable name in the file.
;;    varnamev: v vector variable name in the file.
;;    varnamec: contour vector variable name in the file.
;;
;; Source prototype
;;    source[*,*]
;;    source!0 = lat
;;    source!1 = lon
;;
;; Return value
;;    A graphic variable.
;;
;; Description
;;    Wrapper for gsn_csm_vector_scalar_map_polar (http://www.ncl.ucar.edu/
;;    Document/Graphics/Interfaces/gsn_csm_vector_scalar_map_polar.shtml)
;;    Combines local resources and accepts diag_script-specific resources
;;    Creates plot, according to wks & res
;;
;; Caveats
;;    Selection of defaults for res almost arbitrary
;;    Please check results of all scripts that use this routine if modifying
;;    the defaults!
;;
;; Modification history
;;    20141016-A_vanu_be: adapted from contour_map_polar.ncl.
;;
local areaplot, datac, datau, datav, funcname, scriptname, sourcec, sourceu, \
    sourcev, res, loc_res, varnamec, varnameu, varnamev, wks_in
begin

    funcname = "vector_scalar_map_polar"
    scriptname = "diag_scripts/shared/plot/vector_scalar_map_polar.ncl"
    enter_msg(scriptname, funcname)

    ;; Get data, either directly or via netCDF file
    if(typeof(sourceu) .eq. "string") then
        datau = ncdf_read(sourceu, varnameu)
    else
        datau = sourceu
        copy_VarMeta(sourceu, datau)
    end if
    if(typeof(sourcev) .eq. "string") then
        datav = ncdf_read(sourcev, varnamev)
    else
        datav = sourcev
        copy_VarMeta(sourcev, datav)
    end if
    if(typeof(sourcec) .eq. "string") then
        datac = ncdf_read(sourcec, varnamec)
    else
        datac = sourcec
        copy_VarMeta(sourcec, datac)
    end if

    ;; additional plot resources to be added to res from diag_script
    loc_res                       = True
    loc_res                       = res
    loc_res@gsnDraw               = False
    loc_res@gsnFrame              = False       ; Don't advance frame.
    loc_res@gsnLeftStringOrthogonalPosF = 0.1   ; shift the string up a bit
    loc_res@cnFillOn              = True        ; Turn on color
    loc_res@lbLabelBarOn          = False       ; Turn on later in panel
    loc_res@vcVectorDrawOrder     = "PostDraw"  ; draw vectors last

    ;; specify the hemisphere
    if (diag_script_info@region .eq. "Antarctic") then
        loc_res@gsnPolar          = "SH"
        if (isatt(diag_script_info, "max_lat")) then
            loc_res@mpMaxLatF     = diag_script_info@max_lat
        else
            loc_res@mpMaxLatF     = 0.0
        end if
    else if (diag_script_info@region .eq. "Arctic")
        loc_res@gsnPolarNH        = True
        if (isatt(diag_script_info, "min_lat")) then
            loc_res@mpMinLatF     = diag_script_info@min_lat
        else
            loc_res@mpMinLatF     = 0.0
        end if
    else
        error_msg("f", "vector_scalar_map_polar.ncl", funcname, \
                  "no (valid) region selection: " + \
                  diag_script_info@region)
    end if
    end if

    ;; Draw plot
    areaplot = gsn_csm_vector_scalar_map_polar(wks_in, \
                                               datau, datav, datac, \
                                               loc_res)

    ;; outfile name
    if (isatt(wks_in, "fullname")) then
        outfile = wks_in@fullname
    else
        outfile = wks_in@name
    end if
    log_info(" Wrote " + outfile)

    leave_msg(scriptname, funcname)
    return(areaplot)

end
