;;#############################################################################
;; DIAGNOSTIC SCRIPT for reproducing IPCC ch. 9 fig. 9.31
;; Author: Axel Lauer (DLR, Germany)
;; CRESCENDO project
;;#############################################################################
;;
;; Description
;;    Calculates the two leading Empirical Orthogonal Functions (EOFs) of
;;    outgoing longwave radiation (OLR). The 20- to 100-day filtered OLR
;;    from observations and each of the CMIP5 historical simulations is
;;    projected on these two leading EOFs to obtain MJO Principal Component
;;    time series. A scatterplot shows the maximum positive correlation
;;    between the resulting MJO principal components (PC) and the time lag at
;;    which it occurred for all winters (November to March). The maximum
;;    positive correlation is an indication of the coherence with which the
;;    MJO convection propagates from the Indian Ocean to the Maritime
;;    Continent/western Pacific, and the time lag is approximately one
;;    fourth of the period of the MJO (Constructed following Sperber and
;;    Kim, 2012.)
;;
;; Required diag_script_info attributes (diagnostics specific)
;;     none
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;     none
;;
;; Required variable_info attributes (variable specific)
;;     none
;;
;; Optional variable_info attributes (variable specific)
;;     none
;;
;; Caveats
;;     none
;;
;; Modification history
;;    20170309-A_laue_ax: written.
;;
;;#############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"

load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"

load "./plot_scripts/ncl/xy_line.ncl"

begin

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(diag_script, "", 2)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    ref_model = var_attr_ref(0)  ; produces an error if ref_model is not set in the namelist
    info_output("reference model = " + ref_model, verbosity, 1)
    idx_ref = get_ref_model_idx(models, var_attr_ref(0))

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    ;; Read model data
    data = True
    do imod = 0, dim_MOD - 1
        data = read_data(imod, var0, field_type0)
    end do

end

begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin
    write_references(diag_script,                  \
                     (/"A_laue_ax"/),              \  ; authors
                     (/""/),                       \  ; contributors
                     (/"D_sperber12asl"/),         \  ; diag_references
                     (/""/),                       \  ; obs_references
                     (/"P_crescendo"/))               ; proj_references

    ;; Basename of diag_script
    diag_script_base = basename(diag_script)

    ;; Define file type
    file_type = getenv("ESMValTool_output_file_type")
    if (ismissing(file_type)) then
        file_type = "ps"
    end if

    ;; Check field type
    plot_fields = (/"T2Ds"/)
    if (field_type0.ne.plot_fields) then
        error_msg("f", diag_script, "", "can't use field " + field_type + \
                  " for this diag_script")
    end if
    delete(plot_fields)

    ; define region of interest
    latmin = -20.0
    latmax = 20.0
    lonmin = 45.0  ;0.0
    lonmax = 240.0 ;360.0

    ; output plot directory
    plot_dir = getenv("ESMValTool_plot_dir")
    output_dir = get_output_dir(plot_dir, diag_script_base)

    rmax = new(dim_MOD, float)
    rmaxlag = new(dim_MOD, float)

    ;; Loop over models
    do imod = 0, dim_MOD - 1

        info_output("Processing " + models@name(imod), verbosity, 0)

        ; extract model corresponding data
        var  = extract_data_extended(imod, var0, data, -1, 0, 0)

        ; extract region
        regvar = area_operations(var, latmin, latmax, \
                 lonmin, lonmax, "extract", False)
        delete(var)

        ; compute Lanczos weights
        ihp = 2      ; 2 = band-pass filter
        nWgt = 201
        sigma = 1.0 ; Lanczos sigma
        fca = 1.0 / 100.0    ; start freq
        fcb = 1.0 / 20.0     ; end freq
        wgt = filwgts_lanczos(nWgt, ihp, fca, fcb, sigma)

        ymd = cd_calendar(regvar&time, -2)

        if (isatt(ymd, "calendar")) then
            if (ymd@calendar .eq. "proleptic_gregorian") then
                ymd@calendar = "gregorian"
            end if
        end if

        yyyyddd = yyyymmdd_to_yyyyddd(ymd)

        ClmDay = clmDayTLL(regvar, yyyyddd)  ; mean annual cycle (daily resolution)
        no_of_harmonics = 3
        ClmDay_sm = smthClmDayTLL(ClmDay, no_of_harmonics)
        anom = calcDayAnomTLL(regvar, yyyyddd, ClmDay_sm)  ; daily anomalies

        delete(ClmDay)
        delete(ClmDay_sm)
        delete(regvar)
        delete(yyyyddd)
        delete(ymd)

        ; apply bandpass filter (time dimension = leftmost dimension)
        xbp = wgt_runave_Wrap(anom(lat|:, lon|:, time|:), wgt, 0)

        ; remove means of band pass series: *not* necessary
        xbp = dim_rmvmean(xbp)

        delete(anom)

;        ; extract winter months, i.e. November through March
;
;        date = cd_calendar(xbp&time, 0)
;        month = date(:, 1)
;        idx = ind(month.eq.1 .or. month.eq.2 .or. month.eq.3 .or. month.eq.11 .or. month.eq.12)
;        timeseries = xbp(:, :, idx)
;
;        delete(date)
;        delete(idx)
;        delete(month)

        ; calculate fir two EOFs
;        eof = eofunc_Wrap(timeseries, 2, False)
;        eof_ts = eofunc_ts_Wrap(timeseries, eof, False)
        eof = eofunc_Wrap(xbp, 2, False)
        eof_ts = eofunc_ts_Wrap(xbp, eof, False)

        delete(xbp)
;        delete(timeseries)

        ; calculate cross-correlation between the resulting MJO principal components 1 + 2
        maxlag = 20
        result01 = esccr(eof_ts(0, :), eof_ts(1, :), maxlag)
        result10 = esccr(eof_ts(1, :), eof_ts(0, :), maxlag)

        ccr = new((/2 * maxlag + 1/), float)
        ccr(maxlag:)  = result01
        ccr(0:maxlag) = result10(::-1)       ; reverse order

        lag = ispan(-maxlag, maxlag, 1)      ; days
        imaxcorr = maxind(ccr)

        rmaxlag(imod) = lag(imaxcorr)
        rmax(imod) = ccr(imaxcorr)

        delete(lag)
        delete(ccr)
        delete(result01)
        delete(result10)

;        system ("rm debug.nc")
;        debugfile = addfile("debug.nc", "c")
;        debugfile->eof = eof
;        debugfile->eof_ts = eof_ts

        ; save EOFs from observation (ref_model) for plotting

        if (imod.eq.idx_ref) then
            eof1 = eof(0, :, :)
            eof2 = eof(1, :, :)
            eof1@vpercent = sprintf("%5.2f", eof@pcvar(0)) + "%"
            eof2@vpercent = sprintf("%5.2f", eof@pcvar(1)) + "%"
        end if

        delete(eof)
        delete(eof_ts)

    end do  ; imod-loop

    ; ===========================================================================================
    ; plot results - based on script "Fig9_31.ncl" from Hongmei Li (MPI, Germany)
    ; ===========================================================================================

    ;********************
    ; plot parameters
    ;********************

    outfile = output_dir + diag_script_base + "_" + var0
    wks = gsn_open_wks(file_type, outfile)
    gsn_define_colormap(wks, "BlueWhiteOrangeRed")  ; "psl"
    plot_L = new(2, graphic)

    ;********************
    ; first plot
    ;********************

    res                 = True           ; plot mods desired
    res@gsnDraw         = False
    res@gsnFrame        = False          ; don't advance frame
    res@cnFillOn        = True           ; color fill
    res@cnLineLabelsOn  = False          ; turn off contour lines
    res@cnLinesOn       = False
    res@mpMinLatF       = -20.
    res@mpMaxLatF       = 20.
    res@mpMinLonF       = 45.
    res@mpMaxLonF       = 240.
    res@mpCenterLonF    = 180.
    res@gsnAddCyclic    = False
    res@gsnLeftString   = " "            ; no left string
    res@gsnRightString  = " "            ; no right string
    res@lbLabelBarOn = False
    res@tmYLLabelFontHeightF  = 0.016
    res@tmXBLabelFontHeightF  = 0.016
    res@gsnStringFontHeightF  = 0.016
    res@tmXBMajorOutwardLengthF = 0.009
    res@tmYLMajorOutwardLengthF = 0.009
    res@tmXBMinorOutwardLengthF = 0.0045
    res@tmYLMinorOutwardLengthF = 0.0045
    res@tmXBMajorLengthF = 0.009
    res@tmYLMajorLengthF = 0.009
    res@tmXBMinorLengthF = 0.0045
    res@tmYLMinorLengthF = 0.0045

    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels = (/-0.07, -0.05, -0.03, -0.01, 0., 0.01, 0.03, 0.05, 0.07/)

    ; this controls the size and location of the first plot
    res@mpShapeMode     = "FreeAspect"
    res@vpWidthF        = 0.6
    res@vpHeightF       = 0.3
 
    res@gsnLeftString   = "(a) EOF1"
    res@gsnRightString  = eof1@vpercent ; title
    plot_L(0) = gsn_csm_contour_map_ce(wks, eof1, res)
    res@gsnLeftString    = "(b) EOF2"
    res@gsnRightString  = eof2@vpercent ; title
    plot_L(1) = gsn_csm_contour_map_ce(wks, eof2, res)

    pres_L = True
    pres_L@gsnFrame = False
    pres_L@gsnPanelLabelBar= True             ; label bar on panel
    pres_L@lbOrientation   = "Horizontal"     ; vertical label bar
    pres_L@lbLabelFontHeightF  = 0.014
    pres_L@gsnPanelRight = 0.55
    gsn_panel(wks, plot_L, (/2, 1/), pres_L)

    ;*****************************************
    ; Third plot
    ;*****************************************

    if (isatt(diag_script_info, "styleset")) then
    else
        diag_script_info@styleset = "CMIP5"
    end if

    colortab = project_style(diag_script_info, "colors")
    markertab = project_style(diag_script_info, "markers")

    delete(res)
    res                    = True
    res@gsnDraw            = False       ; Don't draw plot
    res@gsnFrame           = False       ; Don't advance frame
    ;
    ; Don't use gsnMaximize. Instead, we will control
    ; exactly where plot goes.
    ;
    res@vpXF              = 0.1
    res@vpYF              = 0.9
    res@vpWidthF          = 0.4
    res@vpHeightF         = 0.63

    ;---Some calculations to get min/max for Y axis.

    ;---Resources for X/Y axes.
    res@trYMinF            =  0.               ; min value on y-axis
    res@trYMaxF            =  0.8              ; max value on y-axis
    res@trXMinF            =  -20.             ; min value on x-axis
    res@trXMaxF            =  20.              ; max value on x-axis

    ;---XY marker resources
    res@xyMarkLineMode     = "Markers"
    res@xyMarkerSizeF      = 0.014
    res@xyMarkerThicknessF = 3.
    res@xyMonoMarkerColor  = False          ; Allow different colors for markers
    res@xyMarkerColors     = colortab
    res@xyMarkers          = markertab

    ;---Tickmark resources
    res@tmXBMode         = "Manual"
    res@tmYLMode         = "Manual"
    res@tmXBTickSpacingF = 5.
    res@tmYLTickSpacingF = 0.2

    ;---Title resources
    res@tiMainOffsetXF     = -0.058
    res@tiMainFontHeightF  = 0.015
    res@tiMainString       = ""
    res@tiXAxisString      = "Lag of Maximum Positive Correlation (Day)"
    res@tiYAxisString      = "Maximum Positive Correlation (Rmax)"
    res@tmYLLabelFontHeightF  = 0.016
    res@tmXBLabelFontHeightF  = 0.016
    res@tiXAxisFontHeightF    = 0.016
    res@tiYAxisFontHeightF    = 0.016
    res@gsnStringFontHeightF  = 0.016
    res@tmXBMajorOutwardLengthF = 0.00
    res@tmYLMajorOutwardLengthF = 0.00
    res@tmXBMinorOutwardLengthF = 0.00
    res@tmYLMinorOutwardLengthF = 0.00
    res@tmXBMajorLengthF = 0.009
    res@tmYLMajorLengthF = 0.009
    res@tmXBMinorLengthF = 0.0045
    res@tmYLMinorLengthF = 0.0045

    ;----------------------------------------------------------------------
    ; Create the scatter plot, but don't draw it yet.
    ;----------------------------------------------------------------------

    res@gsnLeftString = "(c) Maximum Positive Correlation"
    plot_R = gsn_csm_xy(wks, transpose((/rmaxlag, rmaxlag/)), \
                      transpose((/rmax, rmax/)), res)

    ;----------------------------------------------------------------------
    ; Draw some markers and text.
    ;----------------------------------------------------------------------

    labels = models@name
    ypos2  = fspan(0.78, 0.39, 22)

    mkres               = True         ; Marker resources
    txres               = True         ; Text resources
    txres@txFontHeightF = 0.01
    txres@txJust        = "CenterLeft"

    do i = 0, dim_MOD - 1
        mkres@gsMarkerThicknessF = 3.
        mkres@gsMarkerSizeF      = 0.012
        mkres@gsMarkerIndex      = markertab(i)
        mkres@gsMarkerColor      = colortab(i)
        gsn_polymarker_ndc(wks, 0.67, ypos2(i), mkres)
        gsn_text_ndc(wks, labels(i), 0.69, ypos2(i), txres)
    end do

    pres_R  = True
    pres_R@gsnFrame = False
    pres_R@gsnPanelLabelBar= False                ; label bar on panel
    pres_R@gsnPanelLeft = 0.56
    pres_R@gsnPanelBottom = 0.055
    gsn_panel(wks,plot_R,(/1,1/),pres_R)

    frame(wks)             ; advance frame after extra line drawn       

    ; ===========================================================================================

    leave_msg(diag_script, "", 2)

end
