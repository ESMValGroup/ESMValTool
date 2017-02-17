;;#############################################################################
;; ANNUAL CYCLE PLOTS FOR SEA ICE AREA AND EXTENT
;;#############################################################################
;; Author: Daniel Senftleben (DLR, Germany)
;; ESMVal project
;;#############################################################################
;;
;; Description
;;    Creates an annual cycle (line) plot of sea ice area and extent.
;;
;; Required variable_info attributes (variable specific)
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    multi_model_mean: additional plot MMM and standard deviation.
;;    region: region to be plotted ( "Arctic" or "Antarctic").
;;    range_option: select specific or intersecting time range.
;;
;; Optional diag_script_info attributes (diagnostics specific)
;;    EMs_in_lg: create legend label for each individual ensemble member
;;               (Default: "True")
;;
;; Caveats
;;    CMOR units expected: areacello in m2, sic in %.
;;
;; Modification history
;;    * 20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                          of the code
;;    * 20150325-A_laue_ax: modified reference tags used for acknowledgements
;;                          (projects, observations, etc)
;;    * 20140312-A_righ_ma: adapted to new time_operations structure.
;;    * 20140129-A_senf_da: modularized to fit to new structure.
;;    * 20130603-A_gott_kl: implemented into ./diag_scripts.
;;    * 20130405-A_gott_kl: written based on previous code by Melanie Braeu.
;;
;;#############################################################################

load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/ensemble.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/set_operators.ncl"

load "./plot_scripts/ncl/aux_plotting.ncl"
load "./plot_scripts/ncl/xy_line.ncl"

begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    ;; References Information
    write_references(diag_script, \
                     "A_gott_kl", \
                     (/"A_hueb_mi", "A_eval_ma", "A_brae_me", "A_senf_da"/), \
                     "D_0001", \
                     "", \
                     (/"P_embrace", "P_climval"/))

    do idx_mod = 0, dim_MOD - 1
        data = read_data(idx_mod, var0, field_type0)
    end do
end

; Load var_info
begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin
    ;; Get environment variables
    file_type = getenv("ESMValTool_output_file_type")
    if(ismissing(file_type)) then
        file_type = "PS"
    end if

    ;; Check field type
    plot_fields = (/"T2Ms"/)
    ind_f = inlist(field_type0, plot_fields)
    if (all(ismissing(ind_f))) then
        error_msg("f", diag_script, "", "can't use field, " + field_type0)
    end if
    delete(plot_fields)
    delete(ind_f)

    ;; Select colors and other plotting attributes
    ;; See ./diag_scripts/lib/ncl/style.ncl
    colors = project_style(diag_script_info, "colors")
    dashes = project_style(diag_script_info, "dashes")
    thicks = project_style(diag_script_info, "thicks")
    annots = project_style(diag_script_info, "annots")
    avgstd = project_style(diag_script_info, "avgstd")
    ;; Select colors and other plotting attributes for multi-model mean
    if (diag_script_info@multi_model_mean .eq. "y") then
        ;; project_style evaluates metadata of variable "models"
        temp = models

        ;; -> keep original "models" in "temp" and restore later
        copy_VarMeta(models, temp)
        delete(models)

        ;; Use "models" to pass on attribute names
        models = getvaratts(temp)
        do i = 0, dimsizes(models) - 1
            ;; Define all original attributes again, but empty
            models@$models(i)$ = ""
        end do
        models@name = "model_mean"

        ;; See ./diag_scripts/lib/ncl/style.ncl
        colors_mmm = project_style(diag_script_info, "colors")
        dashes_mmm = project_style(diag_script_info, "dashes")
        thicks_mmm = project_style(diag_script_info, "thicks")
        annots_mmm = project_style(diag_script_info, "annots")
        avgstd_mmm = project_style(diag_script_info, "avgstd")
        delete(models)
        models = temp   ; Restore original "models"
        copy_VarMeta(temp, models)
        delete(temp)
    end if

    ;; Output file name specifications
    ;; (passsed to get_wks for final output file name)
    outfile_ext  = "extent_" + var0 + "_" + diag_script_info@region
    outfile_area = "area_" + var0 + "_" + diag_script_info@region

    ;; Start processing
    ;; Array for collecting sea-ice extent for all models
    val_ext = new((/dimsizes(annots), 12/), float)
    val_ext!1 = "month"
    val_ext&month = int2dble(ispan(1, 12, 1))

    ;; Array for collecting sea-ice area for all models
    val_area = new((/dimsizes(annots), 12/), float)
    val_area!1 = "month"
    val_area&month = int2dble(ispan(1, 12, 1))

    ;; Loop over models
    do imod = 0, dimsizes(models@name) - 1

        info_output("processing " + annots(imod) + \
                    "  ############################################", \
                    verbosity, 2)

        ;; See ./interface_scripts/data_handling.ncl
        var = extract_data(imod, data, -1, 0, 0)
        if (var@units .eq. "%") then
            var = var / 100.  ; Subsequent code expects %/100
            var@units = "%/100"
        end if

        ;; Further processing depends on whether var is on a lat-lon grid
        if(all(iscoord(var, (/"lat", "lon"/)))) then
            info_output("SeaIce_ancyc.ncl: " + var0 + \
                        " is on a lat-lon grid", verbosity, 3)
            area_temp = map_area(var&lat, var&lon)
            ;; Make 2D to be compatible with irregular grid handling below
            area_lat = conform_dims(dimsizes(area_temp), var&lat, 0)
            area_lon = conform_dims(dimsizes(area_temp), var&lon, 1)
            if (isatt(var, "_FillValue")) then
                if (isatt(area_temp, "_FillValue")) then
                    area_temp = where(area_temp .eq. area_temp@_FillValue, \
                                      var@_FillValue, area_temp)
                end if
                area_temp@_FillValue = var@_FillValue
            end if
        else  ; Irregular grid
            info_output("SeaIce_ancyc.ncl: var is not on a lat-lon grid", \
                        verbosity, 3)
            ;; Get area data & lat-lon from another file
            ;; (model specific: ~/reformat_scripts folder)
            area_lon = extract_area(imod, data, "lon", var0)
            area_lat = extract_area(imod, data, "lat", var0)
            area_temp = extract_area(imod, data, "areacello", var0)
        end if

        ;; Determine time range (years)
        start_year = get_start_year(imod)
        end_year   = get_end_year(imod)
        if (diag_script_info@range_option.eq.0) then
            range = "different years (" + \
                sprinti("%0.4i", min(toint(models@start_year))) + "-" + \
                sprinti("%0.4i", max(toint(models@end_year))) + ")"
        else if (diag_script_info@range_option.eq.1) then
            range = sprinti("%0.4i", start_year) + " - " + \
                sprinti("%0.4i", end_year)
        end if
        end if

        ;; Time axis
        var_temp = time_operations(var, start_year, end_year, "extract", "", 0)
        delete(var)
        var = var_temp
        copy_VarMeta(var_temp, var)
        delete(var_temp)

        ;; Select region -> name corresponds to the pole region, but entire
        ;; hemisphere will be evaluated and convert areas to million sqkm
        if (diag_script_info@region .eq. "Antarctic") then
            ;; Set area_temp to 0 in NH
            area_temp = where(area_lat .lt. 0.0, area_temp / 1.e12, 0.0)
        else if (diag_script_info@region .eq. "Arctic")
            ;; Set area_temp to 0 in SH
            area_temp = where(area_lat .gt. 0.0, area_temp / 1.e12, 0.0)
        else
            error_msg("f", diag_script, "", "no (valid) region selection: " + \
                      diag_script_info@region)
        end if
        end if

        ;; Calculate sea ice area
        var_area = var * conform(var, area_temp, (/1, 2/))
        copy_VarMeta(var, var_area)

        ;; Calculate sea ice extent
        extent_temp = where(var .ge. 0.15, 1., 0.)
        var_ext = conform(var, area_temp, (/1, 2/)) * extent_temp
        delete(extent_temp)
        copy_VarMeta(var, var_ext)
        delete(var)

        ;; Calculate monthly mean for every lon/lat coordinate
        var_months_ext = time_operations(var_ext, start_year, end_year, \
                                         "average", "monthlyclim", True)
        var_months_area = time_operations(var_area, start_year, end_year, \
                                          "average", "monthlyclim", True)

        ;; Loop over months
        do imon = 0, 11

            ;; Sum up over lon and lat and put into global array
            temp_ext = dim_sum_Wrap(var_months_ext(imon, :, :))
            temp_area = dim_sum_Wrap(var_months_area(imon, :, :))
            val_ext(imod, imon)  = tofloat(dim_sum_Wrap(temp_ext))
            val_area(imod, imon) = tofloat(dim_sum_Wrap(temp_area))
            delete([/temp_ext, temp_area/])

        end do  ; Loop over months

        ;; Clean up
        delete([/var_ext, var_area, area_lon, area_lat, \
               area_temp, var_months_ext, var_months_area/])

    end do  ; imod

    ;; Calculate standard deviation of models
    if (diag_script_info@multi_model_mean .ne. "y") then
        ;; define anyway, because fields are needed as parameters for xy_line
        val_ext_stddev  = 0
        val_area_stddev = 0
    else
        ;; Mean, stddev, -1 * stddev, +1 * stddev
        val_ext_stddev  = new((/4, 12/), float)
        val_area_stddev = new((/4, 12/), float)
        ;; See ./diag_scripts/lib/ncl/style.ncl for which data shall be
        ;; included in the statistics
        temp = ind(avgstd .eq. 0)
        do imon = 0, 11
            ;; For extent
            val_ext_stddev(0, imon) = \
                dim_avg_n_Wrap(val_ext(temp, imon), (/0/))
            val_ext_stddev(1, imon) = \
                dim_stddev_n_Wrap(val_ext(temp, imon), (/0/))
            val_ext_stddev(2, imon) = \
                val_ext_stddev(0, imon) - val_ext_stddev(1, imon)
            val_ext_stddev(3, imon) = \
                val_ext_stddev(0, imon) + val_ext_stddev(1, imon)

            ;; For area
            val_area_stddev(0, imon) = \
                dim_avg_n_Wrap(val_area(temp, imon), (/0/))
            val_area_stddev(1, imon) = \
                dim_stddev_n_Wrap(val_area(temp, imon), (/0/))
            val_area_stddev(2, imon) = \
                val_area_stddev(0, imon) - val_area_stddev(1, imon)
            val_area_stddev(3, imon) = \
                val_area_stddev(0, imon) + val_area_stddev(1, imon)
        end do
        delete(temp)
    end if

    ;; Some attributes for extent plot
    val_ext&month@units = "Multi-year monthly mean"
    val_ext@units = "10~S~6~N~ km~S~2~N~"
    val_ext@long_name = "Sea Ice Extent"

    ;; Some attributes for area plot
    val_area&month@units = "Multi-year monthly mean"
    val_area@units = "10~S~6~N~ km~S~2~N~"
    val_area@long_name = "Sea Ice Area"

;;*********** PLOTTING ************

    ;; Create plot variables
    wks_ext  = get_wks(0, diag_script, outfile_ext)
    wks_area = get_wks(0, diag_script, outfile_area)

    ;; Define plot resources
    res             = True
    res@gsnPaperOrientation = "landscape"

    ;; Min and max values on x-axis
    res@trXMinF = min(val_ext&month) - 0.05 * \
        (max(val_ext&month) - min(val_ext&month))
    res@trXMaxF = max(val_ext&month) + 0.25 * \
        (max(val_ext&month) - min(val_ext&month))

    ;res@vpXF     = 0.18               ; Control X position of first plot.
    ;res@vpYF     = 0.85
    res@tmXBMode        = "Explicit"
    res@tmXBValues      = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)
    res@tmXBLabels      = (/"J", "F", "M", "A", "M", "J",\
                           "J", "A", "S", "O", "N", "D"/)
    res@tiXAxisString = range

    ;; Resources specific for plot_ext
    res_ext = True
    copy_VarMeta(res, res_ext)
    res_ext@tiMainString  = val_ext&month@units + " " + \
        diag_script_info@region + " " + val_ext@long_name
    res_ext@tiYAxisString = val_ext@units

    ;; Resources specific for plot_area
    res_area = True
    copy_VarMeta(res, res_area)

    ;; Min val_area on y-axis
    res_area@tiMainString  = val_area&month@units + " " + \
        diag_script_info@region + " " + val_area@long_name
    res_area@tiYAxisString = val_area@units

    ;; Draw the plots (procedure in ./plot_scripts/ncl/xy_line.ncl)
    debuginfo = False
    xy_line(wks_ext, val_ext, val_ext&month, val_ext_stddev, res_ext, \
            debuginfo)
    xy_line(wks_area, val_area, val_area&month, val_area_stddev, res_area, \
            debuginfo)

    info_output(" gv " + wks_ext@fullname + " &", verbosity, 1)
    info_output(" gv " + wks_area@fullname + " &", verbosity, 1)

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
