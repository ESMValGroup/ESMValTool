;##############################################################################
;# ANNUAL CYCLE CONTOUR PLOTS OF ZONAL MEANS                                  #
;##############################################################################
;# Author: Mark Stevens (NCAR, USA)                                           #
;# ESMVal project
;##############################################################################
;#                                                                            #
;# Description:                                                               #
;#    * Creates annual cycle contour plots of zonal means (set8) for the AMWG #
;#      diagnostics package (http://www.cgd.ucar.edu/cms/diagnostics/)        #
;#    * The original script has been modified for use within the CCMVal       #
;#      diagnostic software package                                           #
;#    * The CCMVal script script has been furthered developed as a part of    #
;#      the ESMValTool                                                        #
;#                                                                            #
;# Required variable_info attributes (variable specific):                     #
;#    * acp_dcnlvl: vector of contour levels covering range of values for     #
;#                  difference field                                          #
;#    * acp_cnlvl: vector of contour levels for original field                #
;#    * refModel: reference model name (used with acp_dcnlvl)                 #
;#    * colormap: named color map                                             #
;#                                                                            #
;# Required diag_script_info attributes (diagnostics specific):               #
;#    * refModel: this is the reference model. If present, differences with   #
;#      the other models are displayed. If omitted, all models are displayed  #
;#      individually.                                                         #
;#    * rgb_file: color map filename                                          #
;#                                                                            #
;# Caveats:                                                                   #
;#                                                                            #
;# Modification history:                                                      #
;#                                                                            #
;##############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

;; Temporary file updated by Python on the fly
load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/calc_cycle.ncl"
load "./diag_scripts/calc_seasonal.ncl"
load "./diag_scripts/lib/ncl/set_operators.ncl"
load "./diag_scripts/uniform_coordinates.ncl"
load "./diag_scripts/write_infobox.ncl"

;*********************************************************************

begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output("plot SET8 - " + diag_script + " (var: " + \
                variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

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
    plot_fields = (/"T2Ms", "T2Ds", "T1Ms", "C2Ms", "C2Ds", "C1Ms"/)
    ind_f = inlist(field_type0, plot_fields)
    if (all(ismissing(ind_f))) then
        error_msg("f", diag_script, "", "cant' use field " + field)
    end if
    delete(plot_fields)
    delete(ind_f)

    supported_info_transforms = (/"latrange", "lonrange"/)

    ;; Output dir
    ;; 'plot_dir' is fetched from ncl.interface
    diag_script_base = basename(diag_script)
    output_dir = get_output_dir(plot_dir, diag_script_base)

    compare = "MODEL"

    ;; 'output_file_type' if fetched from ncl.interface
    if (isdefined("output_file_type")) then
        file_type = output_file_type
    else
        file_type = "PS"
    end if

    file_type = str_lower(file_type)

;******************************************************************
;check the input parameters

    if (isatt(diag_script_info, "refModel")) then
        idx_ref = ind(models@name .eq. diag_script_info@refModel)
        if (all(ismissing(idx_ref))) then
            error_msg("f", diag_script, "", "refModel is not included in " + \
                      "model list (namelist)")
        end if
        if (isatt(diag_script_info, "acp_dcnlvl")) then
            cntrs = diag_script_info@acp_dcnlvl
        else
            error_msg("f", diag_script, "", "plot requires variable " + \
                      "attribute 'acp_dcnlvl'")
        end if
        p_ref = "ref_"
    else
        idx_ref = -1
        if (isatt(variable_info, "acp_cnlvl")) then
            cntrs = variable_info@acp_cnlvl
        else
            error_msg("f", diag_script, "", "plot requires variable " + \
                      "attribute 'acp_cnlvl'")
        end if
        p_ref = ""
    end if

    if(isvar("latrange")) then
        delete(latrange)
    end if
    if(isatt(diag_script_info, "ac_latrange")) then
        latrange = diag_script_info@ac_latrange
    end if

    if(isvar("lonrange")) then
        delete(lonrange)
    end if
    if(isatt(diag_script_info, "ac_lonrange")) then
        lonrange = diag_script_info@ac_lonrange
    end if

    if (isatt(diag_script_info, "rgb_file")) then
        color_type = "COLOR"
        rgb_file = "./diag_scripts/lib/ncl/rgb/" + diag_script_info@rgb_file
    else
        color_type =  "MONO"
        rgb_file = ""
    end if

;******************************************************************
; common resources for contour plots
    res                        = True
    res@gsnDraw                = False
    res@gsnFrame               = False
    res@txFontHeightF          = 0.018
    res@sfXArray               = ispan(0, 12, 1)
    res@tiMainFontHeightF      = 0.022
    res@tmXBMode               = "Explicit"
    res@tmXBValues             = ispan(0, 12, 1)
    res@tmXBLabels = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", \
                       "O", "N", "D", "J"/)
    if (color_type .eq. "COLOR") then
        res@cnFillOn  = True
        res@cnLinesOn = False
        res@lbTitleOn = True
        res@lbLabelFontHeightF = 0.018
        res@lbTitleFontHeightF = 0.02
        res@lbBoxMinorExtentF = 0.18
        res@cnLevelSelectionMode = "ExplicitLevels"
        res@cnLevels = cntrs
        res@lbLabelStride = 2
        res@lbOrientation = "Horizontal"
    else
        res@cnLevelSelectionMode = "ManualLevels"
        res@cnFillOn  = False
        res@cnLinesOn = True
        res@cnLineLabelsOn = True
        res@cnLineLabelInterval = 1
        res@cnLineLabelFontHeightF =  0.02
        res@cnInfoLabelFontHeightF =  0.02
        cnlvls = 15
        dcnlvls = 12
    end if
    res@gsnLeftString = variable_info@long_name
    if(isvar("lonrange") .and. \
       inlist("lonrange", supported_info_transforms)) then
        res@gsnLeftString = res@gsnLeftString + " (" + lonrange(0) \
                                              + "<lon<" + lonrange(1) + ")"
        res@gsnLeftStringFontHeightF = 14
    end if
    res@gsnRightString = variable_info@units

;;----------------------- Open files for plots ----------------------
    if (color_type .eq. "COLOR") then
        filename_aux = "c"
    else
        filename_aux = ""
    end if

    output_filename = interface_get_figure_filename(diag_script_base, \
                                                    var0, \
                                                    field_type0,\
                                                    filename_aux, -1)
    output_filename = output_filename + "." + file_type

    wks = gsn_open_wks(file_type, output_dir + output_filename)

    if(isatt(variable_info, "colormap")) then
        gsn_define_colormap(wks, variable_info@colormap)
    end if

    if (idx_ref .ne. -1) then
        ;; (b) difference of each model from 'reference'
        plot = new(dim_MOD - 1, "graphic")

        ;; CASE 1
        case1 = models@name(idx_ref) + " " + models@case_name(idx_ref)
        A0 = extract_data(idx_ref, data, -1, 0, 0)
        uniform_coordinates(A0, field_type0)

        ;; convert to T2Ms, C2Ms
        if ((field_type0 .eq. "T2Ds") .or. (field_type0 .eq. "C2Ds"))then
            tmp = day_to_month(A0, field_type0)
            delete(A0)
            A0 = tmp
            delete(tmp)
        end if

        ;; Select a longitude sub-set to plot
        if (iscoord(A0, "lon") .and. isvar("lonrange") \
            .and. inlist("lonrange", supported_info_transforms)) then
            Atmp = A0(time|:, lat|:, {lon|lonrange(0):lonrange(1)})
            delete(A0)
            A0 = Atmp
            delete(Atmp)
        end if
        A = calc_cycle(A0)
        A&lat@units = "degrees_north"
        lat1 = A&lat
        nlat1 = dimsizes(lat1)
        delete(A0)

        jmod = 0
        do imod = 0, dim_MOD - 1
            if (imod .ne. idx_ref) then
                ;; CASE 2
                case2 = models@name(imod) + " " + models@case_name(imod)
                B0 = extract_data(imod, data, -1, 0, 0)
                uniform_coordinates(B0, field_type0)

                ;; convert to T2Ms, C2Ms
                if ((field_type0 .eq. "T2Ds") .or. \
                    (field_type0 .eq. "C2Ds")) then
                    tmp = day_to_month(B0, field_type0)
                    delete(B0)
                    B0 = tmp
                    delete(tmp)
                end if

                ;; Select a longitude sub-set to plot
                if (iscoord(B0, "lon") .and. isvar("lonrange") \
                    .and. inlist("lonrange", supported_info_transforms)) then
                    Btmp = B0(time|:, lat|:, {lon|lonrange(0):lonrange(1)})
                    delete(B0)
                    B0 = Btmp
                    delete(Btmp)
                end if
                B = calc_cycle(B0)
                B&lat@units = "degrees_north"
                lat2 = B&lat
                nlat2 = dimsizes(lat2)
                delete(B0)

                ;;---------------------------------------------------------
                ;; difference plot of time(x) vs. latitude(y)
                ;; check for different number of latitudes and then
                ;; lineary interpolate to the smaller number of latitudes

                if (nlat1 .ne. nlat2) then  ; lat grids different
                    if (nlat1 .gt. nlat2) then
                        tmp = linint1(lat1, A, False, lat2, 0)  ; a(time, lat)
                        C = B  ; copy dims, coords
                        C = (/tmp - B/)  ; get diff values
                        res@sfYArray = lat2
                    else
                        tmp = linint1(lat2, B, False, lat1, 0)  ; b(time, lat)
                        C = A  ; copy dims, coords
                        C = (/A - tmp/)  ; get diff values
                        res@sfYArray = lat1
                    end if
                    delete(tmp)
                else  ; same grid latitudes
                    C = A  ; copy dims, coords
                    C = (/A - B/)  ; get diff values
                    res@sfYArray = lat1
                end if

                mind = min(C)
                maxd = max(C)

                res@tiMainString = case1 + " - " + case2
                if (color_type .eq. "COLOR") then
                    res@lbLabelStride = 1
                    res@lbTitleString = "MIN = " + sprintf("%6.2e", mind) \
                                        + " MAX = " + sprintf("%6.2e", maxd)
                else
                    mnmxint = nice_mnmxintvl(mind, maxd, dcnlvls, False)
                    res@cnMinLevelValF = mnmxint(0)
                    res@cnMaxLevelValF = mnmxint(1)
                    res@cnLevelSpacingF = mnmxint(2)
                end if

                if(isvar("latrange") .and. \
                   inlist("latrange", supported_info_transforms))then
                    lat1!0 = "lat"
                    lat1&lat = lat1
                    delete(res@sfYArray)
                    res@sfYArray = lat1({lat|latrange(0):latrange(1)})
                    plot(imod) = gsn_csm_lat_time( \
                        wks, C({lat|latrange(0):latrange(1)}, time|:), res)
                else
                    plot(jmod) = gsn_csm_lat_time(wks, C(lat|:, time|:), res)
                end if
                if (color_type .eq. "MONO") then
                    plot(jmod) = ZeroNegDashLineContour(plot(jmod))
                else
                    delete(res@lbTitleString)
                end if

                delete(mind)
                delete(maxd)
                delete(lat2)
                delete(nlat2)
                delete(B)
                delete(C)
                delete(res@sfYArray)

                jmod = jmod + 1
            end if
        end do  ; imod

        delete(lat1)
        delete(nlat1)
        delete(A)

    else
        ;; (a) full field
        plot = new(dim_MOD, "graphic")

        do imod = 0, dim_MOD - 1
            ;; CASE 1
            case1 = models@name(imod) + " " + models@case_name(imod)
            if ((stringtoint(models@end_year(imod)) \
                - stringtoint(models@start_year(imod))) .gt. 0) then
                yrs_ave1 = models@start_year(imod) + "-" + \
                    models@end_year(imod)
            end if
            A0 = extract_data(imod, data, -1, 0, 0)
            uniform_coordinates(A0, field_type0)

            ;; convert to T2Ms, C2Ms
            if ((field_type0 .eq. "T2Ds") .or. (field_type0 .eq. "C2Ds"))then
                tmp = day_to_month(A0, field_type0)
                delete(A0)
                A0 = tmp
                delete(tmp)
            end if

            ;; Select a longitude sub-set to plot
            if (iscoord(A0, "lon") .and. isvar("lonrange") \
                .and. inlist("lonrange", supported_info_transforms)) then
                Atmp = A0(time|:, lat|:, {lon|lonrange(0):lonrange(1)})
                delete(A0)
                A0 = Atmp
                delete(Atmp)
            end if
            A = calc_cycle(A0)
            A&lat@units = "degrees_north"
            lat1 = A&lat
            delete(A0)

            ;;------------------------------------------------------------
            ;; case contour plots of time(x) vs. latitude(y)
            min1 = min(A)
            max1 = max(A)

            if (color_type .eq. "COLOR") then
                res@lbLabelStride = 2
                res@lbTitleString = "MIN = " + sprintf("%6.2e", min1) \
                                    + " MAX = " + sprintf("%6.2e", max1)
            else
                mnmxint = nice_mnmxintvl(min1, max1, cnlvls, False)
                res@cnMinLevelValF = mnmxint(0)
                res@cnMaxLevelValF = mnmxint(1)
                res@cnLevelSpacingF = mnmxint(2)
            end if

            if (isvar("yrs_ave1")) then
                res@tiMainString = case1 + " (yrs " + yrs_ave1 + ")"
                delete(yrs_ave1)
            else
                res@tiMainString = case1
            end if
            if(isvar("latrange") .and. inlist("latrange", \
                                              supported_info_transforms)) then
                lat1!0 = "lat"
                lat1&lat = lat1
                res@sfYArray = lat1({lat|latrange(0):latrange(1)})
                plot(imod) = gsn_csm_lat_time(\
                    wks, A({lat|latrange(0):latrange(1)}, time|:), res)
            else
                res@sfYArray = lat1
                plot(imod) = gsn_csm_lat_time(wks, A(lat|:, time|:), res)
            end if
            if (color_type .eq. "MONO") then
                plot(imod) = ZeroNegDashLineContour(plot(imod))
            else
                delete(res@lbTitleString)
            end if
            delete(res@sfYArray)
            delete(res@tiMainString)

            delete(min1)
            delete(max1)
            delete(A)
            delete(lat1)
        end do  ; imod

    end if

    ;; Panel plots
    pan = True
    pan@gsnMaximize = True
    pan@gsnFrame = False
    pan@gsnPaperOrientation = "portrait"
    n_p = dimsizes(plot) / 3
    if (mod(dimsizes(plot), 3) .ne. 0) then
        n_p = n_p + 1
    end if
    gsn_panel(wks, plot, (/n_p, 3/), pan)
    frame(wks)

    delete(plot)
    delete(res)
    delete(pan)

    info_output("wrote " + output_dir + output_filename, verbosity, 1)
    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
