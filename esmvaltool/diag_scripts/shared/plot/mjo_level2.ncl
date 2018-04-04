
;; Temporary file updated by Python on the fly


load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/logging.ncl"

load "./diag_scripts/shared/set_operators.ncl"
load "./diag_scripts/shared/statistics.ncl"

load "./diag_scripts/aux/SAMonsoon/SAMonsoon.ncl"

undef("mjo_wave_freq_plot")
function mjo_wave_freq_plot(wks[1]:graphic, \
                             wf[*][*]:numeric,\
                          dummy[*][*]:graphic,\
                        idx_mod[1]:numeric,\
                            opt[1]:logical,\
                             di[1]:logical)
; Description:
;   Place the 2D-model fields stored in the 'storage_vault' on a set of paneled
;   figures. First entry is always the reference data set. The various logcials
;   carry switches for the plot routine behaviour.
local res, NW, fMin, fMax, day1, fline1, day2, fline2, gsres, tres
begin
    ;; Set plot setting
    res = True                   ; plot mods desired
    res@gsnFrame = False
    res@gsnDraw = False
    res@cnFillOn = True          ; turn on color
    res@gsnSpreadColors = True   ; use full range of colormap
    res@lbLabelAutoStride = True
    res@cnLinesOn = False        ; no contour lines
    res@cnLineLabelsOn = False
    res@cnLevelFlags = False
    res@cnInfoLabelOn = False
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = (/0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.05/)
    res@lbLabelBarOn = False
    res@gsnLeftString = ""
    res@gsnRightString = ""
    res@tiXAxisString = ""

    gsres = True
    gsres@gsLineDashPattern = 11

    tres = True
    tres@txFontHeightF = 0.0175

    NW = 6
    fMin = -0.05
    fMax = 0.05
    day1 = di@min_days
    fline1 = 1. / day1
    day2 = di@max_days
    fline2 = 1. / day2

    if(isatt(opt, "gsnLeftString")) then
        res@gsnLeftString = opt@gsnLeftString
    end if
    if(isatt(opt, "gsnCenterString")) then
        res@gsnCenterString = opt@gsnCenterString
    end if
    if(isatt(opt, "gsnRightString")) then
        res@gsnRightString = opt@gsnRightString
    end if

    plot = gsn_csm_contour(wks, wf({0:NW}, {fMin:fMax}), res)

    dummy(idx_mod, 0) = gsn_add_polyline(wks,\
                                         plot,\
                                        (/fline1, fline1/),\
                                        (/ 0., NW/),\
                                         gsres)

    dummy(idx_mod, 1) = gsn_add_text(wks,\
                                     plot,\
                                    (day1 + "d"),\
                                     fline1 + 0.005,\
                                     0.93 * NW,\
                                     tres)

    dummy(idx_mod, 2) = gsn_add_polyline(wks,\
                                         plot,\
                                        (/fline2, fline2/),\
                                        (/ 0., NW/),\
                                         gsres)

    dummy(idx_mod, 3) = gsn_add_text(wks,\
                                     plot,\
                                    (day2 + "d"),\
                                     fline2 + 0.005,\
                                     0.93 * NW,\
                                     tres)
    return(plot)
end

; -----------

undef("addHorVertLinesCross_extended")
function addHorVertLinesCross_extended(wks[1]    : graphic,\
                                      plot[1]    : graphic,\
                                        nw[1],\
                                      dumy[*][*] : graphic,\
                                   idx_mod[1]    : numeric)
; freq [y] axis:  Add horizontal lines that explicitly
;                 print time in days. This assumes the units
;                 of the freq axis are "cpd" [cycles per day]
local gsres, txres, xx, dely, m, nwl, nwr
begin
    gsres = True
    gsres@gsLineDashPattern = 1

    nwl = -nw + 3.5   ; left
    nwr = nw          ; right
    dumy(idx_mod, 0) = gsn_add_polyline(wks, plot, (/0, 0/), (/0.0, 0.5 /), gsres)
    dumy(idx_mod, 1) = gsn_add_polyline(wks, plot, (/nwl, nwr/), (/1./80, 1./80/), gsres)
    dumy(idx_mod, 2) = gsn_add_polyline(wks, plot, (/nwl, nwr/), (/1./20, 1./20/), gsres)
    dumy(idx_mod, 3) = gsn_add_polyline(wks, plot, (/nwl, nwr/), (/1./10, 1./10/), gsres)
    dumy(idx_mod, 4) = gsn_add_polyline(wks, plot, (/nwl, nwr/), (/1./5, 1./5/), gsres)
    dumy(idx_mod, 5) = gsn_add_polyline(wks, plot, (/nwl, nwr/), (/1./3, 1./3/), gsres)

    txres = True
    txres@txJust = "CenterLeft"
    txres@txFontHeightF = 0.013

    xx = -nw + 0.3
    dely = 0.000   ; yy
    dumy(idx_mod, 6) = gsn_add_text(wks, plot, "3 days", xx, (1./3 + dely), txres)
    dumy(idx_mod, 7) = gsn_add_text(wks, plot, "5 days", xx, (1./5 + dely), txres)
    dumy(idx_mod, 8) = gsn_add_text(wks, plot, "10 days", xx, (1./10 + dely), txres)
    dumy(idx_mod, 9) = gsn_add_text(wks, plot, "20 days", xx, (1./20 + dely), txres)
    dumy(idx_mod, 10) = gsn_add_text(wks, plot, "80 days", xx, (1./80 + dely), txres)

    return(plot)
end

;-----------------------------------------------

undef("mjo_cross_spectra_plot")
function mjo_cross_spectra_plot(wks[1]       : graphic, \
                                STC[*][*][*] : numeric,\
                        spectraType[1]       : string,\
                             dumdcs[*][*]    : graphic,\
                               dums[*][*]    : graphic,\
                            idx_mod[1]       : numeric,\
                                opt[1]       : logical)
local res
begin
    res = True             ; plot mods desired
    res@gsnDraw = False
    res@gsnFrame = False

    res@cnFillOn = True            ; turn on color
    res@cnFillMode = "RasterFill"  ; match WMO Clivar

    res@gsnSpreadColors = True     ; use full range of colormap
    if(opt .and. isatt(opt, "gsnSpreadColorStart")) then
        res@gsnSpreadColorStart = opt@gsnSpreadColorStart
    end if

    if(opt .and. isatt(opt, "gsnSpreadColorEnd")) then
        res@gsnSpreadColorEnd = opt@gsnSpreadColorEnd
    end if

    res@cnLinesOn = False
    res@cnLineLabelsOn = False
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = 0.05
    res@cnMaxLevelValF = 0.65         ; correlation^2 = 0.8
    res@cnLevelSpacingF = 0.05
    res@cnInfoLabelOn = False
    res@lbLabelBarOn = False          ; no individual label bars

    if(.not. opt .or. .not. isatt(opt, "pltPhase") .or. opt@pltPhase) then
        plotPhase = True
        res@vcRefMagnitudeF = 1.0              ; define vector ref mag
        res@vcRefLengthF = 0.01                ; define length of vec ref
        res@vcRefAnnoOrthogonalPosF = -1.0     ; move ref vector
        res@vcRefAnnoArrowLineColor = "black"  ; change ref vector color
        res@vcMinDistanceF = 0.0075            ; thin out vectors
        res@vcMapDirection = False
        res@vcRefAnnoOn = False                ; do not draw
        res@gsnScalarContour = True            ; contours desired
    else
        plotPhase = False
    end if

    ;;  res@gsnLeftString = "Coh^2: Symmetric"
    if(isatt(opt, "gsnLeftString")) then
        res@gsnLeftString = opt@gsnLeftString
    end if

    res@gsnRightString = "10% = "\
                         + sprintf("%3.2f", STC@prob_coh2(2))\
                         + "   5% = "\
                         + sprintf("%3.2f", STC@prob_coh2(4))

    if(opt .and. isatt(opt, "pltZonalWaveNumber")) then
        nWavePlt = opt@pltZonalWaveNumber
    else
        nWavePlt = 15     ; default
    end if

    ;;---------------------------------------------------------------
    ;; dispersion: curves
    ;;---------------------------------------------------------------
    rlat = 0.0
    Ahe = (/50., 25., 12./)
    nWaveType = 6
    nPlanetaryWave = 50
    nEquivDepth = dimsizes(Ahe)
    Apzwn = new((/nWaveType, nEquivDepth, nPlanetaryWave/), "double", 1e20)
    Afreq = Apzwn
    genDispersionCurves(nWaveType,\
                        nEquivDepth,\
                        nPlanetaryWave,\
                        rlat,\
                        Ahe,\
                        Afreq,\
                        Apzwn)

    ;;---------------------------------------------------------------
    ;; dispersion curve and text plot resources
    ;;---------------------------------------------------------------
    dcres = True
    dcres@gsLineThicknessF = 2.0
    dcres@gsLineDashPattern = 0

    txres = True
    txres@txPerimOn = True
    txres@txFontHeightF = 0.013
    txres@txBackgroundFillColor = "Background"

    ;;---------------------------------------------------------------
    ;; plot symmetric data
    ;;---------------------------------------------------------------
    if(spectraType .eq. "Symmetric") then
        n = 8

        c2s = STC(n, :, {-nWavePlt:nWavePlt})
        c2s@_FillValue = 1e20
        c2s(0, :) = c2s@_FillValue
        c2s = where(c2s .lt. 0.05, c2s@_FillValue, c2s)  ; mask

        n = 12
        phs1 = STC(n, :, {-nWavePlt:nWavePlt})
        phs1@long_name = "symmetric phase-1"
        phs1@_FillValue = c2s@_FillValue
        phs1(0, :) = phs1@_FillValue
        phs1 = where(c2s .lt. 0.05, phs1@_FillValue, phs1)  ; mask

        n = 14
        phs2 = STC(n, :, {-nWavePlt:nWavePlt})
        phs2@long_name = "symmetric phase-2"
        phs2@_FillValue = c2s@_FillValue
        phs2(0, :) = phs2@_FillValue
        phs2 = where(c2s .lt. 0.05, phs2@_FillValue, phs2)  ; mask

        if(opt .and. isatt(opt, "pltProb")) then
            np = ind(c2s@prob_coh2 .eq. opt@pltProb)
            if(.not. ismissing(np)) then
                c2s = where(c2s .lt. STC@prob_coh2(np), c2s@_FillValue, c2s)
                phs1 = where(ismissing(c2s), phs1@_FillValue, phs1)
                phs2 = where(ismissing(c2s), phs2@_FillValue, phs2)
            end if
        end if

        if(opt .and. isatt(opt, "coh2Cutoff")) then
            c2s = where(c2s .lt. opt@coh2Cutoff, c2s@_FillValue, c2s)
        end if

        if(opt .and. isatt(opt, "phaseCutoff")) then
            phs1 = where(c2s .lt. opt@phaseCutoff, phs1@_FillValue, phs1)
            phs2 = where(c2s .lt. opt@phaseCutoff, phs2@_FillValue, phs2)
        end if

        if(opt .and. isatt(opt, "elimIsoVals") .and. .not.opt@elimIsoVals) then
            print("mjo_cross_plot: no values eliminated")
        else
            mjo_elimIsolatedValues(c2s, phs1, phs2, 1)
            mjo_elimIsolatedValues(c2s, phs1, phs2, 2)
        end if

        if(plotPhase) then
            scl_one = sqrt(1./(phs1^2 + phs2^2))
            phs1 = scl_one*phs1
            phs2 = scl_one*phs2
            plot = gsn_csm_vector_scalar(wks, phs1, phs2, c2s, res)
        else
            plot = gsn_csm_contour(wks, c2s, res)
        end if

        plot = addHorVertLinesCross_extended(wks, plot, nWavePlt, dums, idx_mod)

        dumdcs(idx_mod, 0) = gsn_add_polyline(wks, plot, Apzwn(3, 0, :), Afreq(3, 0, :), dcres)
        dumdcs(idx_mod, 1) = gsn_add_polyline(wks, plot, Apzwn(3, 1, :), Afreq(3, 1, :), dcres)
        dumdcs(idx_mod, 2) = gsn_add_polyline(wks, plot, Apzwn(3, 2, :), Afreq(3, 2, :), dcres)
        dumdcs(idx_mod, 3) = gsn_add_polyline(wks, plot, Apzwn(4, 0, :), Afreq(4, 0, :), dcres)
        dumdcs(idx_mod, 4) = gsn_add_polyline(wks, plot, Apzwn(4, 1, :), Afreq(4, 1, :), dcres)
        dumdcs(idx_mod, 5) = gsn_add_polyline(wks, plot, Apzwn(4, 2, :), Afreq(4, 2, :), dcres)
        dumdcs(idx_mod, 6) = gsn_add_polyline(wks, plot, Apzwn(5, 0, :), Afreq(5, 0, :), dcres)
        dumdcs(idx_mod, 7) = gsn_add_polyline(wks, plot, Apzwn(5, 1, :), Afreq(5, 1, :), dcres)
        dumdcs(idx_mod, 8) = gsn_add_polyline(wks, plot, Apzwn(5, 2, :), Afreq(5, 2, :), dcres)

        dumdcs(idx_mod, 9) = gsn_add_text(wks, plot, "Kelvin", 11.5, .40, txres)
        dumdcs(idx_mod, 10) = gsn_add_text(wks, plot, "n = 1 ER", -10.7, .07, txres)
        dumdcs(idx_mod, 11) = gsn_add_text(wks, plot, "n = 1 IG", -3.0, .45, txres)
        dumdcs(idx_mod, 12) = gsn_add_text(wks, plot, "h = 50", -14.0, .78, txres)
        dumdcs(idx_mod, 13) = gsn_add_text(wks, plot, "h = 25", -14.0, .60, txres)
        dumdcs(idx_mod, 14) = gsn_add_text(wks, plot, "h = 12", -14.0, .46, txres)

    ;;---------------------------------------------------------------
    ;; plot asymmetric data
    ;;---------------------------------------------------------------
    else if(spectraType .eq. "Asymmetric") then
        n = 9
        if(isatt(opt, "gsnLeftString")) then
            res@gsnLeftString = opt@gsnLeftString
        end if

        res@gsnRightString = "10% = "\
                             + sprintf("%3.2f", STC@prob_coh2(2))\
                             + " "\
                             + "  5% = "\
                             + sprintf("%3.2f", STC@prob_coh2(4))
        c2a = STC(n, :, {-nWavePlt:nWavePlt})
        c2a@_FillValue = 1e20
        c2a(0, :) = c2a@_FillValue
        c2a = where(c2a .lt. 0.05, c2a@_FillValue, c2a)  ; mask

        n = 13
        pha1 = STC(n, :, {-nWavePlt:nWavePlt})
        pha1@long_name = "asymmetric phase-1"
        pha1@_FillValue = c2a@_FillValue
        pha1(0, :) = pha1@_FillValue

        n = 15
        pha2 = STC(n, :, {-nWavePlt:nWavePlt})
        pha2@long_name = "asymmetric phase-2"
        pha2@_FillValue = c2a@_FillValue
        pha2(0, :) = pha2@_FillValue
        pha2 = where(c2a .lt. 0.05, pha2@_FillValue, pha2)  ; mask

        if(opt .and. isatt(opt, "pltProb")) then
            np = ind(c2a@prob_coh2 .eq. opt@pltProb)
            if(.not. ismissing(np)) then
                c2a = where(c2a .lt. STC@prob_coh2(np), c2a@_FillValue, c2s)
                pha1 = where(ismissing(c2a), pha1@_FillValue, pha1)
                pha2 = where(ismissing(c2a), pha2@_FillValue, pha2)
            end if
        end if

        if(opt .and. isatt(opt, "coh2Cutoff")) then
            c2a = where(c2a .lt. opt@coh2Cutoff, c2s@_FillValue, c2s)
        end if
        if(opt .and. isatt(opt, "phaseCutoff")) then
            pha1 = where(c2a .lt. opt@phaseCutoff, pha1@_FillValue, pha1)
            pha2 = where(c2a .lt. opt@phaseCutoff, pha2@_FillValue, pha2)
        end if

        if(opt .and. isatt(opt, "elimIsoVals") .and. .not.opt@elimIsoVals) then
            mjo_elimIsolatedValues(c2a, pha1, pha2, 2)
            mjo_elimIsolatedValues(c2a, pha1, pha2, 1)
        end if

        if(plotPhase) then
            scl_one = sqrt(1. / (pha1 ^ 2 + pha2 ^ 2))
            pha1 = scl_one * pha1
            pha2 = scl_one * pha2
            plot = gsn_csm_vector_scalar(wks, pha1, pha2, c2a, res)
        else
            plot = gsn_csm_contour(wks, c2a, res)
        end if

        plot = addHorVertLinesCross_extended(wks, plot, nWavePlt, dums, idx_mod)

        dumdcs(idx_mod, 0) = gsn_add_polyline(wks, plot, Apzwn(0, 0, :), Afreq(0, 0, :), dcres)
        dumdcs(idx_mod, 1) = gsn_add_polyline(wks, plot, Apzwn(0, 1, :), Afreq(0, 1, :), dcres)
        dumdcs(idx_mod, 2) = gsn_add_polyline(wks, plot, Apzwn(0, 2, :), Afreq(0, 2, :), dcres)
        dumdcs(idx_mod, 3) = gsn_add_polyline(wks, plot, Apzwn(1, 0, :), Afreq(1, 0, :), dcres)
        dumdcs(idx_mod, 4) = gsn_add_polyline(wks, plot, Apzwn(1, 1, :), Afreq(1, 1, :), dcres)
        dumdcs(idx_mod, 5) = gsn_add_polyline(wks, plot, Apzwn(1, 2, :), Afreq(1, 2, :), dcres)
        dumdcs(idx_mod, 6) = gsn_add_polyline(wks, plot, Apzwn(2, 0, :), Afreq(2, 0, :), dcres)
        dumdcs(idx_mod, 7) = gsn_add_polyline(wks, plot, Apzwn(2, 1, :), Afreq(2, 1, :), dcres)
        dumdcs(idx_mod, 8) = gsn_add_polyline(wks, plot, Apzwn(2, 2, :), Afreq(2, 2, :), dcres)

        dumdcs(idx_mod, 10) = gsn_add_text(wks, plot, "MRG", -10.0, .15, txres)
        dumdcs(idx_mod, 11) = gsn_add_text(wks, plot, "n = 2 IG", -3.0, .58, txres)
        dumdcs(idx_mod, 12) = gsn_add_text(wks, plot, "n = 0 EIG", 6.5, .40, txres)
        dumdcs(idx_mod, 13) = gsn_add_text(wks, plot, "h = 50", -10.0, .78, txres)
        dumdcs(idx_mod, 14) = gsn_add_text(wks, plot, "h = 25", -10.0, .63, txres)
        dumdcs(idx_mod, 15) = gsn_add_text(wks, plot, "h = 12", -10.0, .51, txres)

        end if
    end if

    resP = True
    resP@gsnMaximize = True
    resP@gsnPanelLabelBar = True
    resP@lbLabelStride = 2                   ; every other one
    resP@cnLabelBarEndLabelsOn = True
    resP@cnLabelBarEndStyle = "ExcludeOuterBoxes"

    resP@cnLevelSelectionMode = res@cnLevelSelectionMode
    resP@cnMinLevelValF = res@cnMinLevelValF
    resP@cnMaxLevelValF = res@cnMaxLevelValF    ; correlation^2 = 0.8
    resP@cnLevelSpacingF = res@cnLevelSpacingF

    return(plot)
end
;---------------------------------
undef("mjo_ceof_plot")
function mjo_ceof_plot(wks[1]    : graphic, \
                      ceof[*][*] : numeric,\
               pcv_eof_olr[1]    : numeric,\
              pcv_eof_u850[1]    : numeric,\
              pcv_eof_u200[1]    : numeric,\
                       opt[1]    : logical)
local res, rts, var0, lon, plot
begin
    ;;************************************************
    ;; Multivariate EOF plots
    ;;************************************************
    rts = True
    rts@gsnDraw = False       ; don't draw yet
    rts@gsnFrame = False      ; don't advance frame yet
    rts@gsnScale = True       ; force text scaling

    rts@vpHeightF = 0.40      ; Changes the aspect ratio
    rts@vpWidthF = 0.85
    rts@vpXF = 0.10        ; change start locations
    rts@vpYF = 0.75        ; the plot
    rts@xyLineThicknesses = (/2, 2, 2/)
    rts@xyLineColors = (/"black", "red", "green"/)
    rts@xyDashPatterns = (/0.0, 0.0, 0.0/)    ; line patterns

    rts@trYMaxF = diag_script_info@trYMaxF_ceof
    rts@trYMinF = diag_script_info@trYMinF_ceof

    rts@gsnYRefLine = 0.                 ; reference line

    rts@pmLegendDisplayMode = "Always"   ; turn on legend
    rts@pmLegendSide = "Top"             ; Change location of
    rts@pmLegendParallelPosF = 0.86      ; move units right
    rts@pmLegendOrthogonalPosF = -0.50   ; move units down
    rts@pmLegendWidthF = 0.15            ; Change width and
    rts@pmLegendHeightF = 0.15           ; height of legend.
    rts@lgLabelFontHeightF = 0.0175

    var0 = variables(0)

    if(var0 .eq. "pr-mmday") then
        rts@xyExplicitLegendLabels = (/"Pre:  " + sprintf("%4.1f", pcv_eof_olr) + "%",\
                                       "U850: " + sprintf("%4.1f", pcv_eof_u850) + "%",\
                                       "U200: " + sprintf("%4.1f", pcv_eof_u200) + "%" /)
    end if

    if(var0 .eq. "rlut") then
        rts@xyExplicitLegendLabels = (/"OLR:  " + sprintf("%4.1f", pcv_eof_olr) + "%",\
                                       "U850: " + sprintf("%4.1f", pcv_eof_u850) + "%",\
                                       "U200: " + sprintf("%4.1f", pcv_eof_u200) + "%"/)
    end if

    if(isatt(opt, "gsnLeftString")) then
        rts@gsnLeftString = opt@gsnLeftString
    end if

    if(isatt(opt, "gsnCenterString")) then
        rts@gsnCenterString = opt@gsnCenterString
    end if

    if(isatt(opt, "gsnRightString")) then
        rts@gsnRightString = opt@gsnRightString
    end if

    lon = ceof&lon
    plot = gsn_csm_xy(wks, lon, ceof(:, :), rts)
    return(plot)
end

;----------------------

undef("mjo_life_cycle_plot")
procedure mjo_life_cycle_plot(x[*][*][*] : numeric,\
                              u[*][*][*] : numeric,\
                              v[*][*][*] : numeric,\
                        idx_mod[1]       : numeric,\
                             di[1]       : logical)
local amres,           ang,               angBnd,\
    aux_title_info,  diag_script_base,  f,\
    field_type0,     iday,              imon,\
    mjo_indx,        nDays,             netcdf_dir,\
    nn,              nPhase,            nt,\
    output_dir,      output_filename,   input_netcdf,\
    pc1,             pc2,               phaBnd,\
    phase,           pLabel,            plot,\
    pltSubTitle,     r2d,               res,\
    resP,            time,              txid,\
    txres,           uAvg,              var0,\
    vAvg,            wks,               xAvg,\
    ymdhms


begin
    ;***********************************************************
    ; Open PC components file created in 'mjo_14.ncl'
    ;***********************************************************
    var0 = variables(0)
    field_type0 = field_types(0)

    netcdf_dir = plot_dir + "netcdf/"
    aux_title_info = "MJO_PC_INDEX"
    input_netcdf = interface_get_figure_filename("mjo_life_cycle_plot",\
                                                  var0,\
                                                  field_type0,\
                                                  aux_title_info,\
                                                  idx_mod)
    f = addfile(netcdf_dir + input_netcdf + ".nc", "r")

    ;;***********************************************************
    ;; Find the indices corresponding to the start/end times
    ;;      on the MJO_PC_INDEX.nc file
    ;;***********************************************************
    ;;  TIME = f->time                          ; days since ...
    ;;  YMD = cd_calendar(TIME, -2)            ; entire(time, 6)
    ;
    ;;  iStrt = ind(YMD .eq. twStrt)              ; index start
    ;;  iLast = ind(YMD .eq. twLast)              ; index last
    ;;  delete(TIME)
    ;;  delete(YMD)
    ;;***********************************************************
    ;; Read the data for the desired period
    ;;***********************************************************
    pc1 = f->PC1(:)
    pc2 = f->PC2(:)
    mjo_indx = f->MJO_INDEX(:)

    time = pc1&time
    ymdhms = cd_calendar(time, 0)

    imon = floattoint(ymdhms(:, 1))   ; convenience
    iday = floattoint(ymdhms(:, 2))   ; subscripts must be integer
    ;;***********************************************************
    ;; Place each array into an appropriate array
    ;;***********************************************************
    phase = (/202.5, 247.5, 292.5, 337.5,\   ; phi: Matt Wheeler's f77 code
               22.5,  67.5, 112.5, 157.5/)
    nPhase = dimsizes(phase)
    phaBnd = 180.0 / nPhase
    angBnd = new((/2, nPhase/), "float")
    angBnd(0, :) = phase - phaBnd
    angBnd(1, :) = phase + phaBnd

    r2d = 180. / (4. * atan(1.0))
    ang = atan2(pc2, pc1) * r2d     ; phase space
    nn = ind(ang .lt. 0)
    ang(nn) = ang(nn) + 360       ; make 0 to 360

    nDays = new(nPhase, "integer")
    pLabel = "P" + ispan(1, nPhase, 1) + ": "

    res = True
    res@gsnDraw = False          ; don't draw yet
    res@gsnFrame = False         ; don't advance frame yet
    res@gsnSpreadColors = True   ; spread out color table

    res@mpFillOn = False        ; turn off map fill
    res@mpMinLatF = di@latrange_life(0)   ; zoom in on map
    res@mpMaxLatF = di@latrange_life(1)
    res@mpCenterLonF = 210.
    res@cnFillOn = True          ; turn on color fill
    res@cnLinesOn = False        ; True is default
    res@cnLineLabelsOn = False   ; True is default
    res@lbLabelBarOn = False     ; turn off individual lb's
    res@gsnScalarContour = True  ; contour 3rd array
    res@gsnMajorLatSpacing = 15
    res@gsnMajorLonSpacing = 60
    res@tmXBLabelFontHeightF = 0.01
    res@tmYLLabelFontHeightF = 0.01

    res@cnLevelSelectionMode = "ManualLevels"

    if(var0 .eq. "rlut") then
        res@cnMinLevelValF = -40    ; -98 mnmxint(0)
        res@cnMaxLevelValF = 40     ;   80 mnmxint(1)
        res@cnLevelSpacingF = 5     ;   20 mnmxint(2)
    end if

    if(var0 .eq. "pr-mmday") then
        res@cnMinLevelValF = -8     ; -100 mnmxint(0)
        res@cnMaxLevelValF = 8      ;   80 mnmxint(1)
        res@cnLevelSpacingF = 1     ;   20 mnmxint(2)
    end if

    res@vcMinDistanceF = 0.01              ; thin the vector density
    res@vcRefMagnitudeF = 2.0              ; define vector ref mag
    res@vcRefLengthF = 0.025               ; define length of vec ref
    res@vcRefAnnoOrthogonalPosF = -1.0     ; move ref vector
    res@vcRefAnnoArrowLineColor = "black"  ; change ref vector color
    res@vcRefAnnoArrowUseVecColor = False  ; don't use vec color for ref

                                            ; panel plot only resources
    resP = True                    ; modify the panel plot
    resP@gsnMaximize = True        ; large format
    resP@gsnPanelLabelBar = True   ; add common colorbar
    resP@lbLabelAutoStride = True  ; auto stride on labels
    resP@lbLabelFontHeightF = 0.01
    resP@gsnPanelBottom = 0.05     ; add some space at bottom
    resP@pmLabelBarWidthF = 0.8    ; label bar width
    resP@pmLabelBarHeightF = 0.05
    resP@gsnPanelFigureStringsFontHeightF = 0.0125  ; bit larger than default

    amres = True
    amres@amOrthogonalPosF = 0.75
    amres@amJust = "CenterCenter"

    ;;*******************************************
    ;; Loop over each phase
    ;;*******************************************
    res@gsnLeftString = ""
    res@gsnRightString = ""

    ;; Output dir
    ;; 'plot_dir' is fetched from ncl.interface
    diag_script_base = basename(diag_script)
    output_dir = config_user_info@plot_dir

    pltSubTitle = "Anomalous: OLR, U850, V850"

    do nSeason = 1, 2
        aux_title_info = model_info@model(idx_mod) + "_"\
                         + di@season_life_cycle(nSeason - 1)
        output_filename = interface_get_figure_filename(diag_script_base,\
                                                        var0,\
                                                        field_type0,\
                                                        aux_title_info,\
                                                        -1)
        wks = gsn_open_wks(output_file_type, output_dir + output_filename)
        gsn_define_colormap(wks, "ViBlGrWhYeOrRe")
        plot = new(nPhase, graphic)              ; create graphic array

        txres = True
        txres@txFontHeightF = 0.01
        txid = gsn_create_text(wks, pltSubTitle, txres)

        if(nSeason .eq. 1) then
            resP@txString = model_info@model(idx_mod) + ": "\
                            + model_info@start_year(idx_mod) + "-"\
                            + model_info@end_year(idx_mod)\
                            + ": May to Oct"
        else
            resP@txString = model_info@model(idx_mod) + ": "\
                            + model_info@start_year(idx_mod) + "-"\
                            + model_info@end_year(idx_mod)\
                            + ": Nov to Apr"
        end if
        do n = 0, nPhase - 1
            ;; All times where the following criteria are met
            if(nSeason .eq. 1) then
                nt = ind(mjo_indx .gt. 1.0\
                         .and.(imon .ge. 5 .and. imon .le. 10)\
                         .and. ang .ge. angBnd(0, n)\
                         .and. ang .lt. angBnd(1, n))
            else
                nt = ind(mjo_indx .gt. 1.0\
                         .and.(imon .ge. 11 .or. imon .le. 4)\
                         .and. ang .ge. angBnd(0, n) \
                         .and. ang .lt. angBnd(1, n))
            end if
            if(.not. all(ismissing(nt))) then
                xAvg = dim_avg_n_Wrap(x(nt, :, :), 0)
                uAvg = dim_avg_n_Wrap(u(nt, :, :), 0)
                vAvg = dim_avg_n_Wrap(v(nt, :, :), 0)

                nDays(n) = dimsizes(nt)

                res@tmXBLabelsOn = False      ; do not draw lon labels
                res@tmXBOn = False            ; lon tickmarks
                if(n .eq. (nPhase - 1)) then
                    res@tmXBLabelsOn = True   ; draw lon labels
                    res@tmXBOn = True         ; tickmarks
                end if

                plot(n) = gsn_csm_vector_scalar_map_ce(wks, uAvg, vAvg, xAvg, res)
            end if
            delete(nt)     ; will change next iteration
        end do

        resP@gsnPanelFigureStrings = pLabel + nDays
        gsn_panel(wks, plot, (/nPhase, 1/), resP)     ; now draw as one plot

        delete(wks)
        delete(plot)
    end do
end
