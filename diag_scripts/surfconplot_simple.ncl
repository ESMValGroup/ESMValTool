;##############################################################################
;# CONTOUR PLOTS OF 2D FIELDS                                                 #
;# Author: Mark Stevens (NCAR, USA), Richard Neale (NCAR, USA)                #
;# ESMVal project
;##############################################################################
;#                                                                            #
;# Description:                                                               #
;#    Creates longitude vs latitude 2D-field contour plots (set5) for the     #
;#    AMWG diagnostics package (http://www.cgd.ucar.edu/cms/diagnostics/)     #
;#    The original script has been modified for use within the CCMVal         #
;#    diagnostic software package                                             #
;#    The CCMVal script script has been furthered developed as a part of      #
;#    the ESMValTool                                                          #
;#                                                                            #
;# Required variable_info attributes (variable specific):                     #
;#    scp_dcnlvl: vector of contour levels covering range of values for       #
;#                difference field                                            #
;#    scp_cnlvl: vector of contour levels for original field                  #
;#    refModel: reference model name (used with scp_dcnlvl)                   #
;#    rgb_file: color map filename                                            #
;#                                                                            #
;# Optional variable_info attributes (variable specific):                     #
;#    seasons: (optional) array of for averaging. (DJF, MAM, JJA, SON or ANN  #
;#             for annual mean). If not set, use default                      #
;#             (/"DJF", "JJA", "ANN"/)                                        #
;#                                                                            #
;# Caveats:                                                                   #
;#                                                                            #
;# Modification history:                                                      #
;#    * 200801: ERA40 reanalysis                                              #
;#    * 200303: ISCCP fields                                                  #
;#    * 200101: updated with significance tests                               #
;#                                                                            #
;##############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;; Temporary file updated by Python on the fly
load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/calc_seasonal.ncl"
load "./diag_scripts/lib/ncl/set_operators.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/functions_contrib.ncl"
load "./diag_scripts/uniform_coordinates.ncl"

begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    print("++++++++++++++++++++++++++++++++++++++++++")
    print("plot " + diag_script + " (var: " + variables(0) + ")")
    print("++++++++++++++++++++++++++++++++++++++++++")

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

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
    ;; 'output_file_type' if fetched from ncl.interface
    if (isdefined("output_file_type")) then
        file_type = output_file_type
    else
        file_type = "PS"
    end if


    plot_fields = (/"T2Ms"/)
    ind_f = inlist(field_type0, plot_fields)
    if (all(ismissing(ind_f))) then
        error_msg("f", diag_script, "", "can't use field, " + field_type0)
    end if
    delete(plot_fields)
    delete(ind_f)

    ;;--------------------------
    ;; Check diag_script info
    ;;--------------------------
    if(isatt(diag_script_info, "seasons")) then
        seasons = diag_script_info@seasons
    else
        error_msg("f", diag_script, "", "plot requires diag_script_info " + \
                  "attribute 'seasons'")
    end if

    if(isatt(diag_script_info, "ref_model")) then
        ref_model = diag_script_info@ref_model
    else
        ref_model = "None"
    end if

    if(isatt(diag_script_info, "cn_levels_mean")) then
        cn_levels_mean = diag_script_info@cn_levels_mean
    else
        error_msg("f", diag_script, "", "plot requires diag_script_info " + \
                  "attribute 'cn_levels_mean'")
    end if

    if(isatt(diag_script_info, "cn_levels_mean_diff")) then
        cn_levels_mean_diff = diag_script_info@cn_levels_mean_diff
    end if

    ;; Output dir
    ;; 'plot_dir' is fetched from ncl.interface
    diag_script_base = basename(diag_script)
    output_dir = get_output_dir(plot_dir, diag_script_base)

    ;; ------------------------- GLOBAL CE PLOTS ------------------------
    ;; common resources for global contour plots
    res = True
    res@mpCenterLonF = 180
    res@mpGeophysicalLineThicknessF = 1.0
    res@gsnDraw = False           ; Do not draw plot
    res@gsnFrame = False          ; Do not advance frame
    res@gsnTickMarksOn = False
    res@tiMainFont = "Helvetica"
    res@vpXF = 0.05
    res@vpWidthF = 1.00
    res@vpHeightF = 0.70
    res@mpPerimOn = True
    res@mpPerimDrawOrder = "PostDraw"
    res@mpFillOn = False
    res@lbLabelBarOn = True
    res@lbTitleOn = True
    res@lbBoxMinorExtentF = 0.15
    res@lbLabelStride = 1
    res@txFontHeightF = 0.020
    res@tiMainFontHeightF = 0.025
    res@lbOrientation = "Vertical"
    res@lbLabelFontHeightF = 0.015
    res@lbTitleFontHeightF = 0.018
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnFillOn = True
    res@cnLinesOn = False
    res@cnMissingValFillColor = "Background"
    res@cnLineLabelsOn = False

    if (isatt(diag_script_info, "colormap")) then
        colormap = diag_script_info@colormap
        res@cnFillPalette     = colormap        ; Set colormap
        res@cnSpanFillPalette = True            ; Use full range of color map
    end if

    title_ref = "_ref"
    idx_ref = ind(models@name .eq. ref_model)
    if (ismissing(idx_ref)) then
        idx_ref = -1
        title_ref = ""
    end if

    ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ;;loop through seasons
    do s = 0, dimsizes(seasons) - 1
        season = seasons(s)

        ;; ***************************************************************
        ;; ------------------ Open files for plots -----------------------
        aux_title_info = season + title_ref
        output_filename = interface_get_figure_filename(diag_script_base,\
                                                        var0,\
                                                        field_type0,\
                                                        aux_title_info,\
                                                        -1)
        wks = gsn_open_wks(file_type, output_dir + output_filename)

        if (idx_ref .ne. -1) then
            ;; (b) difference of each model from 'reference'
            plot = new(dim_MOD - 1, "graphic")

            ;; CASE 1
            case1 = models@name(idx_ref) + " " + models@case_name(idx_ref)
            A0 = extract_data(idx_ref, data, -1, 0, 0)
            uniform_coordinates(A0, field_type0)

            lat1 = A0&lat
            nlat1 = dimsizes(lat1)
            lon1 = A0&lon
            if (lon1(0) .lt. 0.0) then
                lon1 = lon1 + 180.0
                A0&lon = (/lon1/)
            end if

            ;; Handle odd number of latitudes explicitly
            if (mod(nlat1, 2) .eq. 1) then
                gw1_tmp = latGauWgt(nlat1 - 1, "lat", \
                                    "gaussian weights", "none")
                gw1 = new(nlat1, float)
                gw1(0:(nlat1/2)) = gw1_tmp(0:(nlat1/2))
                gw1(((nlat1/2) + 1):(nlat1 - 1)) = \
                    gw1_tmp((nlat1/2):(nlat1 - 2))
                delete(gw1_tmp)
            else
                ;; An even number of latitudes is ok
                gw1 = latGauWgt(nlat1, "lat", "gaussian weights", "none")
            end if
            gw1!0 = "lat"
            gw1&lat = lat1

            A = calc_SEASON(A0, season)
            A&lat@units = "degrees_north"
            A&lon@units = "degrees_east"
            delete(A0)

            ;; compute area weighted mean of the fields
            gblmean_A = wgt_areaave(A, gw1, 1.0, 0)

            jmod = 0
            do imod = 0, dim_MOD - 1
                if (imod .ne. idx_ref) then
                    ;; CASE 2
                    case2 = models@name(imod) + " " + models@case_name(imod)
                    B0 = extract_data(imod, data, -1, 0, 0)
                    uniform_coordinates(B0, field_type0)

                    lat2 = B0&lat
                    nlat2 = dimsizes(lat2)
                    lon2 = B0&lon
                    if (lon2(0) .lt. 0.0) then
                        lon2 = lon2 + 180.0
                        B0&lon = (/lon2/)
                    end if

                    ;; gw2 = inptr2->gw
                    if (mod(nlat2, 2) .eq. 1) then
                        ;; no gw in file
                        gw2_tmp = latGauWgt(nlat2 - 1, "lat", \
                                            "gaussian weights", "none")
                        gw2 = new(nlat2, float)
                        gw2(0:(nlat2/2)) = gw2_tmp(0:(nlat2/2))
                        gw2(((nlat2/2) + 1):(nlat2 - 1)) = \
                            gw2_tmp((nlat2/2):(nlat2 - 2))
                        delete(gw2_tmp)
                    else
                        ;; no gw in file
                        gw2 = latGauWgt(nlat2, "lat", \
                                        "gaussian weights", "none")
                    end if
                    gw2!0 = "lat"
                    gw2&lat = lat2

                    B = calc_SEASON(B0, season)
                    B&lat@units = "degrees_north"
                    B&lon@units = "degrees_east"
                    delete(B0)

                    ;; compute area weighted mean of the fields
                    gblmean_B = wgt_areaave(B, gw2, 1.0, 0)

                    ;; compute RMSE of the two fields after interpolating to
                    ;; the same grid if needed, and the mean of the
                    ;; difference of the interpolated fields
                    if (nlat1 .eq. nlat2) then  ; same grids
                        gblrmse = wgt_arearmse(A, B, gw1, 1.0, 0)
                        C = A
                        C = (/A - B/)
                        gblmean = wgt_areaave(C, gw1, 1.0, 0)
                    else  ; different grids
                        if (nlat1 .gt. nlat2) then  ; interp A to nlat2
                            tmp = B
                            tmp = linint2(lon1, lat1, A, True, lon2, lat2, 0)
                            gblrmse = wgt_arearmse(tmp, B, gw2, 1.0, 0)
                            C = B
                            C = (/tmp - B/)
                            gblmean = wgt_areaave(C, gw2, 1.0, 0)
                        else  ; interp B to nlat1
                            tmp = A
                            tmp = linint2(lon2, lat2, B, True, lon1, lat1, 0)
                            gblrmse = wgt_arearmse(A, tmp, gw1, 1.0, 0)
                            C = A
                            C = (/A - tmp/)
                            gblmean = wgt_areaave(C, gw1, 1.0, 0)
                        end if
                        delete(tmp)
                    end if

                    ;; ---------------- GLOBAL CE PLOTS ---------------------
                    mind = min(C)
                    maxd = max(C)
                    res@tiMainString = case1 + " - " + case2
                    res@gsnLeftString = "mean = " + sprintf("%6.2f", gblmean)
                    res@gsnCenterString = "rmse = " + sprintf("%6.2f", gblrmse)
                    res@gsnRightString = variable_info@units

                    res@lbTitleString = "Min = " + sprintf("%6.2f", mind) + \
                        " Max = " + sprintf("%6.2f", maxd)
                    res@cnLevels = cn_levels_mean_diff

                    ;; plot differences
                    plot(jmod) = gsn_csm_contour_map_ce(wks, C, res)
                    delete(res@cnLevels)
                    delete(res@tiMainString)
                    delete(res@gsnLeftString)
                    delete(res@gsnRightString)
                    delete(res@gsnCenterString)

                    delete(mind)
                    delete(maxd)

                    delete(lat2)
                    delete(nlat2)
                    delete(lon2)
                    delete(gw2)
                    delete(B)
                    delete(gblmean_B)

                    delete(C)
                    delete(gblrmse)

                    jmod = jmod + 1
                end if
            end do  ; imod

            delete(lat1)
            delete(nlat1)
            delete(lon1)
            delete(gw1)
            delete(A)
            delete(gblmean_A)

        else
            ;; (a) full field
            plot = new(dim_MOD, "graphic")

            do imod = 0, dim_MOD - 1
                ;; CASE 1
                case1 = models@name(imod) + " " + models@case_name(imod)
                if ((stringtoint(models@end_year(imod)) - \
                     stringtoint(models@start_year(imod))) .gt. 0) then
                    yrs_ave1 = models@start_year(imod) + "-" + \
                        models@end_year(imod)
                end if
                A0 = extract_data(imod, data, -1, 0, 0)
                uniform_coordinates(A0, field_type0)

                lat1 = doubletofloat(A0&lat)
                nlat1 = dimsizes(lat1)

                ;; gw1 = inptr1->gw
                if (mod(nlat1, 2) .eq. 1) then
                    ;; no gw in file
                    gw1_tmp = latGauWgt(nlat1 - 1, \
                                        "lat", "gaussian weights", "none")
                    gw1 = new(nlat1, float)
                    gw1(0:(nlat1/2)) = gw1_tmp(0:(nlat1/2))
                    gw1(((nlat1/2) + 1):(nlat1 - 1)) = \
                        gw1_tmp((nlat1/2):(nlat1 - 2))
                    delete(gw1_tmp)
                else
                    ;; no gw in file
                    gw1 = latGauWgt(nlat1, "lat", "gaussian weights", "none")
                end if
                gw1!0 = "lat"
                gw1&lat = lat1

                A = calc_SEASON(A0, season)
                A&lat@units = "degrees_north"
                A&lon@units = "degrees_east"
                delete(A0)

                ;; compute area weighted mean of the fields
                gblmean_A = wgt_areaave(A, gw1, 1.0, 0)

                ;; ---------------- GLOBAL CE PLOTS ---------------------
                min1 = min(A)
                max1 = max(A)

                res@cnLevels = cn_levels_mean
                res@lbTitleString = "Min = " + sprintf("%6.2f", min1) \
                    + " Max = " + sprintf("%6.2f", max1)
                res@cnRasterModeOn = False

                if (isvar("yrs_ave1")) then
                    res@tiMainString = case1 + " (yrs " + yrs_ave1 + ")"
                    delete(yrs_ave1)
                else
                    res@tiMainString = case1
                end if
                res@gsnCenterString = "mean= " + sprintf("%6.2f", gblmean_A)

                if (isatt(variable_info, "long_name")) then
                    res@gsnLeftString = variable_info@long_name
                end if

                if (isatt(variable_info, "units")) then
                    res@gsnRightString = variable_info@units
                end if

                ;; plot model data
                plot(imod) = gsn_csm_contour_map_ce(wks, A, res)

                delete(res@lbTitleString)
                delete(res@tiMainString)
                delete(res@gsnCenterString)

                delete(min1)
                delete(max1)

                delete(lat1)
                delete(nlat1)
                delete(gw1)
                delete(A)
                delete(gblmean_A)
            end do  ; imod

        end if

        ;; common resources for panel plots
        pan = True
        pan@gsnMaximize = True
        pan@gsnFrame = False
        pan@gsnPaperOrientation = "portrait"
        pan@gsnPanelYWhiteSpacePercent = 2
        pan@txFontHeightF = 0.018
        pan@txString = season

        n_p = dimsizes(plot) / 2
        if (mod(dimsizes(plot), 2) .ne. 0) then
            n_p = n_p + 1
        end if
        gsn_panel(wks, plot, (/n_p, 2/), pan)    ; panel plot
        frame(wks)

        delete(plot)
        delete(pan)

    end do

    delete(res)
    delete(seasons)
    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
