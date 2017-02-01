;*********************************************************************
; seacycplot.ncl
;*********************************************************************
; NCL script
; plot_seas_cycle.ncl (set 10)
; Mark Stevens April 2002
;*********************************************************************
; The original NCL script, plot_seas_cycle.ncl, was written to make...
; annual cycle line plots of global means (set10)
; for the AMWG diagnostics package.
; http://www.cgd.ucar.edu/cms/diagnostics/
;
; The NCL script described above has been modified for use with the CCMVal
; diagnostic software package.
; (PIs: Andrew Gettelman & Veronika Eyring, co-I: Hisako Shiona)
;*********************************************************************
load "./plot_type/functions_contrib.ncl"
load "./plot_type/calc_seasonal.ncl"
load "./plot_type/uniform_coordinates.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
;*********************************************************************

begin
    print("++++++++++++++++++++++++++++++++++++++++++")
    print("plot SET10 - seacycplot (var: " + variable + ")")
    print("++++++++++++++++++++++++++++++++++++++++++")

    field_number = getenv("ESMValTool_ORIG_FIELD_NUMBER")
    plot_fields = (/"T2Ms", "T2Ds", "C2Ms"/)
    ind_f = ind(plot_fields .eq. field_number)
    if (all(ismissing(ind_f))) then
        print("fatal:can't use field, " + field_number + " for this plot_type")
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)

    wkdir = plot_dir + "/seacycplot/"
    compare = "MODEL"
    file_type = "PS"
    color_type = "COLOR"

    ;; The 'models' array has been loaded in main.ncl
    dim_MOD = dimsizes(models@name)

    ;;******************************************************************
    ;;--------------------- Open files for plots -----------------------
    if (color_type .eq. "COLOR") then
        if (variable .eq. "ICEFRAC" .or. variable .eq. "SNOWHICE" .or. \
            variable .eq. "SNOWHLND" .or. variable .eq. "PRECST") then
            wksNP = gsn_open_wks(file_type, wkdir + "seacycplot_" + variable + "_NP_c")
            wksSP = gsn_open_wks(file_type, wkdir + "seacycplot_" + variable + "_SP_c")
        else
            wks = gsn_open_wks(file_type, wkdir + "seacycplot_" + variable + "_c")
        end if
    else
        wks = gsn_open_wks(file_type, wkdir + "seacycplot_" + variable)
    end if

    ;;******************************************************************
    ;;------------------------------------------------------------------
    months = ispan(0, 12, 1)
    months!0 = "month"
    gblmean = new((/dim_MOD, 13/), float)
    gblmean!0 = "line"
    gblmean!1 = "month"
    nhmean = new((/dim_MOD, 13/), float)
    shmean = new((/dim_MOD, 13/), float)
    nhmean!0 = "line"
    nhmean!1 = "month"
    shmean!0 = "line"
    shmean!1 = "month"
    lgLabels = new((/dim_MOD/), string)

    do imod = 0, dim_MOD - 1
        ;; CASE 1  A(time, lat, lon)
        case1 = models@name(imod) + " " + models@case_name(imod)
        if ((stringtoint(models@end_date(imod)) \
                  - stringtoint(models@start_date(imod))) .gt. 0) then
            yrs_ave1 = models@start_date(imod) + "-" + models@end_date(imod)
            lgLabels(imod) = case1 + " (" + yrs_ave1 + ")"
        else
            lgLabels(imod) = case1
        end if
        A0 = extract_data(imod, data, -1, 0, 0)
        uniform_coordinates(A0, field_number)

        ;; convert to T2Ms
        if (field_number .eq. "T2Ds") then
            tmp = day_to_month(A0, field_number)
            delete(A0)
            A0 = tmp
            delete(tmp)
        end if
        if isStrSubset(field_number, "T") then
            tmp_clm = clmMonTLL(A0)
            delete(A0)
            A0 = tmp_clm
            delete(tmp_clm)
        end if

        A0&lat@units = "degrees_north"
        lat1 = A0&lat
        nlat1 = dimsizes(lat1)

        ;;gw1 = inptr1->gw
        if (mod(nlat1, 2) .eq. 1) then
            ;; no gw in file
            gw1_tmp = latGauWgt(nlat1 - 1, "lat", "gaussian weights", "none")
            gw1 = new(nlat1, float)
            gw1(0:(nlat1/2)) = gw1_tmp(0:(nlat1/2))
            gw1(((nlat1/2) + 1):(nlat1 - 1)) = gw1_tmp((nlat1/2):(nlat1 - 2))
            delete(gw1_tmp)
        else
            ;; no gw in file
            gw1 = latGauWgt(nlat1, "lat", "gaussian weights", "none")
        end if
        gw1!0 = "lat"
        gw1&lat = lat1

        ;;***************************************************************
        do m = 0, 11         ; loop over months
            A = A0(m, :, :)

            if (variable .eq. "ICEFRAC" .or. variable .eq. "SNOWHICE" .or. \
                variable .eq. "SNOWHLND" .or. variable .eq. "PRECST") then
                if (variable .eq. "ICEFRAC") then
                    SLAT = True
                    if (SLAT) then
                        ;; finite volume lat/lon grid
                        ;;if (isfilevar(inptr1(m), "slat")) then
                        shmean(imod, m) = iceareaFIX(A, 0)
                        nhmean(imod, m) = iceareaFIX(A, 1)
                    else      ; gaussian lat/lon grid
                        shmean(imod, m) = iceareaGAU(A, 0)
                        nhmean(imod, m) = iceareaGAU(A, 1)
                    end if
                end if
                if (variable .eq. "PRECST") then
                    shmean(imod, m) = wgt_areaave(A({-90:-30}, :), gw1({-90:-30}), 1., 0)
                    nhmean(imod, m) = wgt_areaave(A({30:90}, :), gw1({30:90}), 1., 0)
                end if
                if (variable .eq. "SNOWHICE") then
                    shmean(imod, m) = wgt_areaave(A({-90:-30}, :), gw1({-90:-30}), 1., 0)
                    nhmean(imod, m) = wgt_areaave(A({30:90}, :), gw1({30:90}), 1., 0)
                end if
                if (variable .eq. "SNOWHLND") then
                    shmean(imod, m) = wgt_areaave(A({-90:-1}, :), gw1({-90:-1}), 1., 0)
                    nhmean(imod, m) = wgt_areaave(A({1:90}, :), gw1({1:90}), 1., 0)
                end if
            else
                if (compare .eq. "OBS" .and. variable .eq. "PRECT_TRMM") then
                    gblmean(imod, m) = wgt_areaave(A({-38:38}, :), gw1({-38:38}), 1., 0)
                else
                    gblmean(imod, m) = wgt_areaave(A, gw1, 1., 0)
                end if
            end if

            delete(A)
        end do       ; loop over month

        ;; wrap around for plot
        if (variable .eq. "ICEFRAC" .or. variable .eq. "SNOWHICE" .or. \
            variable .eq. "SNOWHLND" .or. variable .eq. "PRECST") then
            nhmean(imod, 12) = nhmean(imod, 0)
            shmean(imod, 12) = shmean(imod, 0)
        else
            gblmean(imod, 12) = gblmean(imod, 0)
        end if

        delete(A0)
        delete(lat1)
        delete(nlat1)
        delete(gw1)
    end do  ; imod

    ;;******************************************************************
    colors = CCMVal2_ColorTable(models@name)
    dashs = CCMVal2_DashTable(models@name)
    LineThick = new((/dim_MOD/), integer)
    LineThick(:) = 2
    LineMode = new((/dim_MOD/), string)
    LineMode(:) = "MarkLines"
    ;;------------------------------------------------------------------
    ;; common plot resources
    res = True
    res@gsnMaximize = True
    res@gsnPaperOrientation = "portrait"
    res@tiMainString = "Annual Cycle Global Mean Climatology"
    res@tiMainFontHeightF = 0.022
    res@tiMainFont = "helvetica-bold"
    res@tiYAxisFontHeightF = 0.015
    if (color_type .eq. "COLOR") then
        res@xyLineColors = colors
    else
        res@xyMonoLineColor = True
        res@xyLineColor = "black"
    end if
    res@xyLineThicknesses = LineThick
    res@xyDashPatterns = dashs
    res@xyMarkLineModes = LineMode
    res@xyMarker = 16
    res@xyMarkerSizeF = 0.01
    res@pmLegendDisplayMode = "Always"
    res@pmLegendSide = "Bottom"
    res@pmLegendWidthF = 0.15
    res@pmLegendOrthogonalPosF = -.018
    res@pmLegendParallelPosF = 0.5
    if (dim_MOD .le. 4) then
        res@pmLegendHeightF = 0.08
        res@lgLabelFontHeightF = 0.016
    else
        if (dim_MOD .le. 8) then
            res@pmLegendHeightF = 0.15
            res@lgLabelFontHeightF = 0.012
        else
            res@pmLegendHeightF = 0.30
            res@lgLabelFontHeightF = 0.008
        end if
    end if
    res@lgPerimOn = True
    res@txFontHeightF = 0.020
    res@tmXBMode = "Explicit"
    res@tmXBLabelFontHeightF = 0.015
    res@tmXBLabels = (/"Jan", "Feb", "Mar", "Apr", "May", "Jun", \
                       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan"/)
    res@tmXBValues = months
    res@xyExplicitLegendLabels = lgLabels

    ;;-------------------------------------------------------------------------
    if (variable .eq. "ICEFRAC" .or. variable .eq. "SNOWHICE" .or. \
        variable .eq. "SNOWHLND" .or. variable .eq. "PRECST") then
        if (var .eq. "ICEFRAC") then
            res@tiYAxisString = "Area x 10~S~6~N~ km~S~2~N~"
            res@gsnLeftString = "Northern Hemisphere Sea-ice area"
            plotNP = gsn_csm_xy(wksNP, months, nhmean, res)
            res@gsnLeftString = "Southern Hemisphere Sea-ice area"
            plotSP = gsn_csm_xy(wksSP, months, shmean, res)
        end if
        if (variable .eq. "PRECST") then
            res@tiYAxisString = "mm/day"
            res@gsnLeftString = "NH snowfall rate"
            plotNP = gsn_csm_xy(wksNP, months, nhmean, res)
            res@gsnLeftString = "SH snowfall rate"
            plotSP = gsn_csm_xy(wksSP, months, shmean, res)
        end if
        if (variable .eq. "SNOWHICE") then
            res@tiYAxisString = "centimeters"
            res@gsnLeftString = "NH depth of snow on seaice"
            plotNP = gsn_csm_xy(wksNP, months, nhmean, res)
            res@gsnLeftString = "SH depth of snow on seaice"
            plotSP = gsn_csm_xy(wksSP, months, shmean, res)
        end if
        if (variable .eq. "SNOWHLND") then
            res@tiYAxisString = "centimeters"
            res@gsnLeftString = "NH depth of snow on land"
            plotNP = gsn_csm_xy(wksNP, months, nhmean, res)
            res@gsnLeftString = "SH depth of snow on land"
            plotSP = gsn_csm_xy(wksSP, months, shmean, res)
        end if
    else
        res@tiYAxisString = info@units
        res@gsnLeftString = info@long_name
        plot = gsn_csm_xy(wks, months, gblmean, res)
    end if

    delete(months)
    delete(nhmean)
    delete(shmean)
    delete(gblmean)
    delete(res)
    if (variable .eq. "ICEFRAC" .or. variable .eq. "SNOWHICE" .or. \
        variable .eq. "SNOWHLND" .or. variable .eq. "PRECST") then
        delete(plotNP)
        delete(plotSP)
    else
        delete(plot)
    end if
    delete(lgLabels)
    delete(colors)
    delete(dashs)
    delete(LineThick)
    delete(LineMode)

end
