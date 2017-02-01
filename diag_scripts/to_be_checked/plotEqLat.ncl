; zonlnplot.ncl
;***************************************
; NCL script
; plot_zonal_lines.ncl
; Mark Stevens, Sept 2001
;***************************************
load "./plot_type/functions_contrib.ncl"
load "./plot_type/calc_seasonal.ncl"
;***************************************

begin
; compare    = getenv("COMPARE")
; plot_type  = getenv("PLOTTYPE")
; color_type = getenv("COLORTYPE")
; season     = getenv("SEASON")
; version    = getenv("DIAG_VERSION")
    wkdir      = getenv("ESMValTool_wrk_dir")
; time_stamp = getenv("TIMESTAMP")
; case_names = getenv("CASENAMES")
; ncdf_mode  = getenv("NCDF_MODE")

    plot_type  = "PS"
    color_type = "COLOR"
    wkdir      = wkdir + "/Equivalent_latitude/"
    time_stamp = False
    version    = "beta"

    do imod = 0, dimsizes(models@name) - 1
        A0 = extract_data(imod, data, -1, 0, 0)
        print(A0)
        delete(A0)
    end do
end
