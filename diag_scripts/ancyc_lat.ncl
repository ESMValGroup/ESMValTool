;##############################################################################
;# ANNUAL CYCLE CONTOUR PLOTS OF ZONAL MEANS                                  #
;##############################################################################
;# Author: Klaus-Dirk Gottschaldt (DLR, Germany)                              #
;# ESMVal project PI: Veronika Eyring (DLR, Germany)                          #
;##############################################################################
;#                                                                            #
;# Description:                                                               #
;#    *                                                                       #
;# Required diag_script_info attributes (diagnostics specific):               #
;#    *                                                                       #
;# Required variable_info attributes (variable specific):                     #
;#    *                                                                       #
;# Caveats:                                                                   #
;# Modification history:                                                      #
;#    * 20151027-A_laue_ax moved call to 'write_references' to the beginning  #
;#                         of the code                                        #
;#    * 20150325-A_laue_ax modified reference tags used for acknowledgements  #
;#                        (projects, observations, etc)                       #
;#    * 20140328 A_gott_kl written based on anncycplot.ncl                    #
;#                                                                            #
;##############################################################################
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/regridding.ncl"
load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/ensemble.ncl"

;; fetch general parameters
begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    ;; References Information
    write_references(diag_script,      \  ; script name
                     "A_gott_kl",      \  ; authors
                     "",               \  ; contributors
                     "D_0000",         \  ; diag_references
                     "",               \  ; obs_references
                     (/"P_esmval"/))          ; proj_references

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    do idx_mod = 0, dim_MOD - 1
        data1 = read_data(idx_mod, var0, field_type0)
    end do
end

;; Load variable_info
begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end


;; Check for special obs
begin
    dim_SPECIAL = 0
    if(isatt(diag_script_info, "obsfiles")) then
        flist = diag_script_info@obsfiles
        dim_SPECIAL = dimsizes(flist)
    end if
end

;; main code
begin
    ;; Check valid plot fields
    plot_fields = (/"T2Ms", "T2Ds", "C2Ms", "C2Ds"/)
    ind_f = inlist(field_type0, plot_fields)
    if (all(ismissing(ind_f))) then
        print("fatal:can't use field, " + field_type0 + \
              " for this diag_script")
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)

    ;; determine destination grid (regridding.ncl)
    idx = ispan(0, dim_MOD - 1, 1)
    dest_grd = find_destination_grid(idx, var0, field_type0, "coarsest")

    ;; default time range & annots
    y1 = toint(models@start_year)
    y2 = toint(models@end_year)
    annots = project_style(diag_script_info, "annots")

    ;; Loop over models to regrid to a reference grid and create
    ;;     mask for missing values
    ifile = -1  ; counter and flag for special obs
    do idx_mod = 0, dim_MOD + dim_SPECIAL - 1

        if(idx_mod.le.dim_MOD - 1) then  ; models via nml
            ;; Extract data
            A0 = extract_data_extended(idx_mod, var0, data1, -1, 0, 0)
            ;; Use overlapping time period or as given in nml?
            if(isatt(diag_script_info, "range_option")) then
                y1(idx_mod) = get_start_year(idx_mod)
                y2(idx_mod) = get_end_year(idx_mod)
                temp = time_operations(A0, y1(idx_mod), y2(idx_mod), \
                                       "extract", "", 0)
                delete(A0)
                A0 = temp
                delete(temp)
            end if  ; Default: use as in nml
            ;; annual cycle (Note: Returns time dimension "month")
            A1 = time_operations(A0, -1, -1, "average", "monthlyclim", True)
        else  ; special obs
            ;; read data
            ifile = idx_mod - dim_MOD
            info_output(" reading " + flist(ifile), verbosity, 2)
            f = addfile(flist(ifile), "r")
            A0 = f->$var0$
            ;; conform dimensions
            dims = (/12, dimsizes(A0&lat), dimsizes(dest_grd&lon)/)
            if(isatt(diag_script_info, "case")) then
                A1 = conform_dims(dims, \
                                  A0(diag_script_info@case, :, :), (/0, 1/))
            else
                A1 = conform_dims(dims, A0(:, :), (/0, 1/))
            end if
            A1!0 = "month"
            A1!1 = "lat"
            A1!2 = "lon"
            A1&month = A0&month
            A1&lat = A0&lat
            A1&lon = dest_grd&lon
            ;; extend fields for annotation
            temp1 = array_append_record(y1, diag_script_info@y1(ifile), 0)
            temp2 = array_append_record(y2, diag_script_info@y2(ifile), 0)
            temp3 = array_append_record(annots, \
                                        diag_script_info@obsname(ifile), 0)
            delete([/y1, y2, annots/])
            y1 = temp1
            y2 = temp2
            annots = temp3
            delete([/temp1, temp2, temp3/])
        end if

        ;; initialize target fields ...
        if(idx_mod.eq.0) then
            ;; ... for regridded data
            if(.not.isatt(A1, "_FillValue")) then
                Fill = default_fillvalue(typeof(A1))
            else
                Fill = A1@_FillValue
            end if
            var_grd = new((/dim_MOD + dim_SPECIAL, 12, dimsizes(dest_grd&lat),\
                          dimsizes(dest_grd&lon)/), typeof(A1), Fill)
            var_grd!0 = "models"
            var_grd!1 = "month"
            var_grd!2 = "lat"
            var_grd!3 = "lon"
            var_grd&month = A1&month
            var_grd&lat = dest_grd&lat
            var_grd&lon = dest_grd&lon

            ;; ... for the mask
            mask0 = var_grd(0, :, :, :)
            mask0(:, :, :) = 1
        end if

        ;; regrid to ref
        A1 = where(ismissing(A1), Fill, A1)
        A1@_FillValue = Fill  ; ensure consistent _FillValue
        var_grd(idx_mod, :, :, :) = \
            area_hi2lores_Wrap(A1&lon, A1&lat, A1, \
            True, 1, dest_grd&lon, dest_grd&lat, False)

        ;; Build the mask: missing in at least one model = missing in the mask
        mask0 = where(ismissing(var_grd(idx_mod, :, :, :)), \
                      mask0@_FillValue, mask0)

        ;; clean up
        delete([/A0, A1/])
    end do  ; idx_MOD, for regridding and mask

    ;; create mask to collect missings over time
    temp = mask0(0, :, :)
    do i = 1, 11
        temp = where(ismissing(mask0(i, :, :)), mask0(i, :, :), 1)
    end do
    mask1 = conform_dims(dimsizes(mask0), temp, (/1, 2/))
    delete(temp)
    info_output(diag_script + ": Mean value calculated after regridding " + \
        "(month, lon, lat) to coarsest grid, masking all points that " + \
        "are missing in at least one model or month, considering area " + \
        "weighting and days per month", verbosity, 1)

    ;; Output dir
    plot_dir = getenv("ESMValTool_plot_dir")
    diag_script_base = basename(diag_script)
    output_dir = get_output_dir(plot_dir, diag_script_base)

    ;; General settings for this plot
    res = True
    res@cnFillOn = True
    res@cnLinesOn = True
    ;res@cnFillPalette =
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = variable_info@acp_cnlvl
    res@txFont = diag_script_info@font
    res@tmXBLabelFont = diag_script_info@font
    res@tmYLLabelFont = diag_script_info@font
    res@tmXBMode = "Explicit"
    res@tmXBValues = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)
    res@tmXBLabels = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", "O", \
                      "N", "D"/)
    if(isatt(diag_script_info, "lbLabelBarOn")) then
        res@lbLabelBarOn = diag_script_info@lbLabelBarOn
    else
        res@lbLabelBarOn = True
    end if
    res@lbLabelFont = diag_script_info@font
    res@lbTitleFont = diag_script_info@font
    res@lbTitleFontHeightF = 0.01  ; default: 0.025
    res@lbTitleOffsetF = -0.25  ; default: 0.03 (negative values move it down)
    res@lbTopMarginF = -0.1  ; default: 0.05 (negative values move it up)
    res@lbTitleJust = "BottomCenter"  ; default: CenterCenter
    ;res@lbTitleString = variable_info@long_name

    ;; modify aspect ratio
    res@vpKeepAspect = False
    res@vpHeightF = 0.3
    res@vpWidthF = 0.6

    ;; Loop over regridded models for plotting
    do idx_mod = 0, dim_MOD + dim_SPECIAL - 1
        ifile = idx_mod - dim_MOD

        ;; calculate mean
        A2 = where(ismissing(mask1), mask1, var_grd(idx_mod, :, :, :))
        copy_VarMeta(mask0, A2)
        area_mean = area_operations(A2, -90, 90, 0, 360, "average", True)
        dpm = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
        time_mean = dim_avg_wgt_Wrap(area_mean, dpm, 1)
        if(time_mean.ge.1000) then
            fmt = tostring(toint(time_mean))
        else if(time_mean.lt.1000 .and. time_mean.ge.100) then
            fmt = sprintf("%5.1f", time_mean)
        else if(time_mean.lt.100 .and. time_mean.ge.10) then
            fmt = sprintf("%5.2f", time_mean)
        else if(time_mean.lt.10 .and. time_mean.ge.1) then
            fmt = sprintf("%5.3f", time_mean)
        else
            fmt = tostring(time_mean)
        end if
        end if
        end if
        end if

        ;; Define output workstation
        file_type = getenv("ESMValTool_output_file_type")
        if(ismissing(file_type)) then
            file_type = "PS"
        end if

        aux_title_info = ""
        if(ifile.lt.0) then
            output_filename = interface_get_figure_filename(diag_script_base,\
                var0, field_type0, aux_title_info, idx_mod)
        else
            output_filename = diag_script_base + "_" + var0 + "_" + \
                field_type0 + \
                "_" + annots(idx_mod) + "_" + tostring(y1(idx_mod)) + "_" + \
                tostring(y1(idx_mod)) + "." + file_type
        end if
        output_file_path = output_dir + output_filename
        wks = gsn_open_wks(file_type, output_file_path)

        ;; Data to plot: zonal mean
        A3 = dim_avg_n_Wrap(var_grd(idx_mod, :, :, :), 2)
        delete(A3&month)
        A3!0 = "time"
        A3&time = ispan(1, 12, 1)

        ;; model specific ressources
        res@tiMainString = annots(idx_mod)
        res@gsnLeftString = variable_info@short_name + \
            " [" + variable_info@units + "]"
        res@gsnRightString = fmt + " " + variable_info@units
;        res@gsnRightString = tostring(y1(idx_mod)) + "-" + \
;            tostring(y2(idx_mod)) + ":  " + fmt + " " + \
;            variable_info@units

        ;; Plotting
        map = gsn_csm_lat_time(wks, A3(lat|:, time|:), res)
        print("gv " + output_file_path + "." + file_type)

        ;; clean up
        delete([/wks, A2, A3, time_mean/])

    end do  ; idx_MOD, for plotting

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
