;;#############################################################################
;; DIAGNOSTIC SCRIPT for reproducing IPCC ch. 9 fig. 9.31
;; Author: Axel Lauer (DLR, Germany)
;; CRESCENDO project
;;#############################################################################
;;
;; Description
;;    Calculates the two leading Empirical Orthogonal Functions (EOFs) of
;;    the observed outgoing longwave radiation (OLR). The 20- to 100-day
;;    filtered OLR from observations and each of the CMIP5 historical simulations is
;;    projected on these two leading EOFs to obtain MJO Principal Component
;;    time series. A scatterplot shows the maximum positive correlation
;;    between the resulting MJO principal components 1 + 2 and the time lag at
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
;; Required variable attributes (defined in namelist)
;;    ref_model: name of reference data set (observations)
;;
;; Caveats
;;     no selection of winters with strong MJO, all winters (Nov-Mar)
;;     will be used for the calculations
;;
;; Modification history
;;    20170622-A_laue_ax: added tags for reporting
;;    20170310-A_laue_ax: written.
;;
;;#############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"

load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/aux/perfmetrics/aux_perfmetrics.ncl"

load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"

load "./plot_scripts/ncl/xy_line.ncl"

load "./diag_scripts/lib/ncl/meta_data.ncl"

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

    climofiles = new(dim_MOD, string)

    do ii = 0, dim_MOD - 1
        climofiles(ii) = \
            interface_get_inpaths(ii) + "/" \
            + interface_get_infile(variables(0), field_types(0), ii)
    end do

    alltags = array_append_record(tags, (/"PT_geo", "ST_eof", "DM_reg"/), 0)
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
        error_msg("f", diag_script, "", "can't use field " + field_type0 + \
                  " for this diag_script")
    end if
    delete(plot_fields)

    ; output plot directory
    plot_dir = getenv("ESMValTool_plot_dir")
    output_dir = get_output_dir(plot_dir, diag_script_base)

    rmax = new(dim_MOD, float)
    rmaxlag = new(dim_MOD, float)

    ; the reference model (observations) has to be processed first
    ; --> create list of all models with the reference model as the first entry

    all_models = ispan(0, dim_MOD -1, 1)
    idx = ind(all_models .ne. idx_ref)
    sorted_models = new(dim_MOD, integer)
    sorted_models(0) = idx_ref
    sorted_models(1:dim_MOD - 1) = all_models(idx)
    delete(idx)
    delete(all_models)

    ; loop over models
    do iloop = 0, dim_MOD - 1

        imod = sorted_models(iloop)

        info_output("Processing " + models@name(imod), verbosity, 0)

        ; get grid definition of reference data

        if (imod .eq. idx_ref) then  ; this is always the first model in the list
            grid = find_destination_grid(idx_ref, var0, field_type0, "coarsest")
        end if

        ; ========================================================
        ; ====================== regridding ======================
        ; ========================================================

        opt = True
        opt@plot_type = "regridded"
        opt@time_avg = "ref"

        regr_path = clim_diag_file_path(var0, field_type0, imod, diag_script, opt)

        ; do not process data if already processed
        if (isfilepresent(regr_path) .and. \
            getenv("ESMValTool_force_calc").eq."False") then
            pre_var = ncdf_read(regr_path, var0)
            dnames = getVarDimNames(pre_var)
        else
            info_output("Regridding " + models@name(imod), verbosity, 2)
            ; extract model corresponding data
            data = read_data(imod, var0, field_type0)
            var = extract_data_extended(imod, var0, data, -1, 0, 0)
            dnames = getVarDimNames(var)
    
            ; determine dimensions to regrid
            xregrid = False  ; lon
            yregrid = False  ; lat
            dnames = getVarDimNames(var)
            if (any(dnames.eq."lon")) then
                xregrid = True
                if (dimsizes(var&lon).eq.dimsizes(grid&lon)) then
                    if (max(abs(var&lon - grid&lon)).lt.1.0d-5) then
                        xregrid = False
                    end if
                end if
            end if
            if (any(dnames.eq."lat")) then
                yregrid = True
                if (dimsizes(var&lat).eq.dimsizes(grid&lat)) then
                    if (max(abs(var&lat - grid&lat)).lt.1.0d-5) then
                        yregrid = False
                    end if
                end if
            end if

           ; regrid in the horizontal (if necessary)
           if (xregrid.and.yregrid) then
                var_hreg = linint2_Wrap(var&lon, var&lat, var, True, \
                                        grid&lon, grid&lat, 0)
            else if (xregrid.and..not.yregrid) then
                dx = ind(dnames.eq."lon")
                var_hreg = linint1_n_Wrap(var&lon, var, True, grid&lon, 0, dx)
                info_output("Latitude regridding not required", verbosity, 4)
            else if (.not.xregrid.and.yregrid) then
                dx = ind(dnames.eq."lat")
                var_hreg = linint1_n_Wrap(var&lat, var, False, grid&lat, 0, dx)
                info_output("Longitude regridding not required", verbosity, 4)
            else if (.not.xregrid .and. .not.yregrid) then
                var_hreg = var
                info_output("Horizontal regridding not required", verbosity, 4)
            end if
            end if
            end if
            end if

            delete(data)
            delete(var)

            pre_var = var_hreg
            delete(var_hreg)

            ;; output to NetCDF
            pre_var@var = var0
            pre_var@diag_script = (/diag_script/)
            pre_var@ncdf = regr_path
            ncdf_outfile = ncdf_write(pre_var, regr_path)
        end if

        ; ========================================================

        ; extract region (0-360 longitude seems to be standard)
        regvar = area_operations(pre_var, -20.0, 20.0, 0.0, 360.0, "extract", False)
        delete(pre_var)

        ; the bandpass filering and calculation of the EOFs is taken from
        ; function "calc_daily_univari_eof" in diag_scripts/aux/mjo/mjo.ncl

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

        ; extract winter months, i.e. November through March

        date = cd_calendar(anom&time, 0)
        month = date(:, 1)
        idx = ind(month.eq.1 .or. month.eq.2 .or. month.eq.3 .or. month.eq.11 .or. month.eq.12)
        anomts = anom(idx, :, :)

        delete(date)
        delete(idx)
        delete(month)
        delete(anom)

        ; apply bandpass filter (time dimension = leftmost dimension)
        xbp = wgt_runave_Wrap(anomts(lat|:, lon|:, time|:), wgt, 0)

        ; remove means of band pass series: *not* necessary
        xbp = dim_rmvmean(xbp)

        delete(anomts)

        ; EOFs are only calculated for the reference data (observations)
        ; (the reference dataset is always the first dataset)

        if (.not.isvar("eof_obs")) then
            ; calculate first two EOFs
            eof_obs = eofunc_Wrap(xbp, 2, False)
        end if

        ; the 20-100 day band-pass filtered daily OLR from each model is
        ; projected onto the observed modes (Sperber and Kim, Atmos. Sci. Let., 2012)

        eof_ts = eofunc_ts_Wrap(xbp, eof_obs, False)

        delete(xbp)

;        system ("rm debug.nc")
;        debugfile = addfile("debug.nc", "c")
;        debugfile->eof_obs = eof_obs
;        debugfile->eof_ts = eof_ts
;        ts = timeseries(time|:, lat|:, lon|:)
;        debugfile->ts = ts
;        exit

;        delete(timeseries)

        ; calculate cross-correlation between the resulting MJO principal components 1 + 2
        maxlag = 20
        result_p = esccr(eof_ts(1, :), eof_ts(0, :), maxlag)
        result_n = esccr(eof_ts(0, :), eof_ts(1, :), maxlag)

        delete(eof_ts)

        ccr = new((/2 * maxlag + 1/), float)
        ccr(0:maxlag-1) = result_n(1:maxlag:-1)  ; "negative lag", -1 reverses order
        ccr(maxlag:)    = result_p(0:maxlag)     ; "positive lag"

        lag = ispan(-maxlag, maxlag, 1)      ; days
        imaxcorr = maxind(ccr)
        rmaxlag(imod) = lag(imaxcorr)
        rmax(imod) = ccr(imaxcorr)

        delete(lag)
        delete(ccr)
        delete(result_n)
        delete(result_p)

    end do  ; iloop

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
    res@mpMinLatF       = -20.0
    res@mpMaxLatF       = 20.0
    res@mpMinLonF       = 45.0
    res@mpMaxLonF       = 240.0
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
    res@gsnRightString  = sprintf("%4.1f", eof_obs@pcvar(0)) + "%" ; title
    plot_L(0) = gsn_csm_contour_map_ce(wks, eof_obs(0, :, :), res)
    res@gsnLeftString    = "(b) EOF2"
    res@gsnRightString  = sprintf("%4.1f", eof_obs@pcvar(1)) + "%" ; title
    plot_L(1) = gsn_csm_contour_map_ce(wks, eof_obs(1, :, :), res)

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

    ; add meta data to plot (for reporting)

    caption = "(a, b) The two leading Empirical Orthogonal Functions (EOFs) of variable " + variables(0) \
              + " (" + models@name(idx_ref) + "). The scatterplot (c) shows the maximum positive correlation " \
              + "between the resulting MJO Principal Components and the time lag at which it occurred for all" \
              + " winters (November to March). Similar to IPCC AR5, fig. 9.31."
    id = diag_script + "_" + variables(0)

    ESMValMD(outfile, alltags, caption, id, variables(0), models@name, climofiles)

    leave_msg(diag_script, "", 2)

end
