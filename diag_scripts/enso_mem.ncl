;; ENSO_MEM
;;#############################################################################
;; Author: Axel Lauer (DLR, Germany)
;; PROJECT-NAME CRESCENDO
;;#############################################################################
;; Description
;;    Calculates maximum entropy power spectum for the Nino3 region (5S-5N, 150W-90W)
;;    from monthly tas anomalies with the climatological annaul cycle removed.
;;    The maximum entropy method (mem) is taken from Numerical Recipes in C
;;    (sections 13.6 and 13.7).
;;
;; Required diag_script_info attributes (diagnostics specific)
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;
;; Required variable_info attributes (variable specific)
;;
;; Optional variable_info attributes (variable specific)
;;    none
;;
;; Caveats
;;    none
;;
;; Modification history
;;    20170622-A_laue_ax: added tags for reporting
;;    20170315-A_laue_ax: written.
;;
;;#############################################################################

; #####################################
; # load external NCL code, if needed #
; #####################################

; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL

load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"

load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/style.ncl"

load "./diag_scripts/aux/perfmetrics/aux_perfmetrics.ncl"

load "./diag_scripts/aux/enso_mem/evlmem.ncl"
load "./diag_scripts/aux/enso_mem/memcof.ncl"

load "./diag_scripts/lib/ncl/meta_data.ncl"

begin
    ; ###########################################################
    ; # Fetch general parameters, set in namelist_cloud.xml,    #
    ; # passed via environment variables by python code         #
    ; ###########################################################

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)

    ; 'diag_script' and 'variables' are fetched from the
    ; above 'ncl.interface' file

    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    ; ###########################################
    ; # references                              #
    ; ###########################################

    write_references(diag_script,      \  ; script name
                     "A_laue_ax",      \  ; authors
                     "",               \  ; contributors
                     "",               \  ; diag_references
                     "",               \  ; obs_references
                     (/"P_crescendo"/))   ; proj_references

    ; 'models', 'variables' and 'field_types' are fetched from the
    ; above 'ncl.interface' file

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

    alltags = array_append_record(tags, (/"PT_pro", "ST_spectrum", "DM_reg"/), 0)
end

begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin
    ;; Basename of diag_script
    diag_script_base = basename(diag_script)

    ;; Define file type
    file_type = getenv("ESMValTool_output_file_type")
    if (ismissing(file_type)) then
        file_type = "ps"
    end if

    ;; Check field type
    if (field_type0.ne."T2Ms") then
        error_msg("f", diag_script, "", "can't use field " + field_type + \
                  " for this diag_script")
    end if

    ;; define region
    region = "Nino3"
    latmin = -5.0
    latmax = 5.0
    lonmin = 210.0
    lonmax = 270.0

    ;; Output plot directory
    plot_dir = getenv("ESMValTool_plot_dir")
    output_dir = get_output_dir(plot_dir, diag_script_base)

    ; results from maximum entropy method (mem)
    steps = 100
    mem_period = new((/steps/), float)
    mem_power = new((/steps, dim_MOD/), float)

    ;; Loop over models
    do imod = 0, dim_MOD - 1

        info_output("Processing " + models@name(imod), verbosity, 0)

        yr0 = stringtoint(models@start_year(imod))
        yr1 = stringtoint(models@end_year(imod))

        ;; Extract model corresponding data
        var  = extract_data_extended(imod, var0, data, -1, 0, 0)

        ; calculate domain average
        regmean = area_operations(var, latmin, latmax, lonmin, lonmax, "average", True)
        ; remove trend

        dimname = getvardims(regmean)
        regmean_dtrend = dtrend_msg(regmean&$dimname$, regmean, False, False)

        regmean_dtrend!0 = regmean!0
        regmean_dtrend&$dimname$ = regmean&$dimname$

        ; calculate climatological seasonal cycle
        seasclim = time_operations(regmean_dtrend, -1, -1, "average", "monthlyclim", True)

        ; subtract climatological seasonal cycle
        dseas = regmean_dtrend
        do i = 0, 11
            dseas(i::12) = dseas(i::12) - seasclim(i)
        end do

        delete(regmean)
        delete(regmean_dtrend)
        delete(seasclim)
        delete(var)

        ; ------------------------------------------------------------------------
        ; estimation of the linear prediction (LP) coefficients d
        ; ------------------------------------------------------------------------

        m = 20  ; default
        xms = 0.0
        d = new((/m/), float)

        memcof(dseas, m, xms, d)

        ; ------------------------------------------------------------------------------
        ; estimate power spectrum as a function of frequency times the sampling interval
        ; (sampling interval = 1 month)
        ; ------------------------------------------------------------------------------

        tmin = 2.0    ; min. period (months)
        tmax = 400.0  ; max. period (months)
        delta = log10(tmax / tmin) / tofloat(steps - 1)
        t0 = log10(tmin)

        do i = 0, steps - 1
            mem_period(i) = 10.0^(t0 + tofloat(i) * delta)  ; months
            mem_power(i, imod) = tofloat(evlmem(1.0 / mem_period(i), d, m, xms))
        end do

        delete(d)

        ; ------------------------------------------------------------------------

        delete(dseas)

    end do  ; imod-loop

    ; ==================================================================================
    ; plot results - based on script "Fig9_35.ncl" from Hongmei Li (MPI, Germany)
    ; ==================================================================================

    ;************************************************
    ; plotting parameters
    ;************************************************

    outfile = output_dir + diag_script_base + "_" + var0 + "_" + region
    wks = gsn_open_wks(file_type, outfile)

    res                   = True        ; plot mods desired
    res@gsnDraw           = False
    res@gsnFrame          = False
    res@tiMainString      = ""          ; add title
    res@vpHeightF= 0.55                 ; change aspect ratio of plot
    res@vpWidthF = 0.75
    res@trYMinF  = 0.001                ; min value on x-axis
    res@trYMaxF  = 1000.                ; max value on x-axis
    res@trXMinF  =  2.                  ; min value on y-axis
    res@trXMaxF  =  400. 
    res@trXLog = True
    res@trYLog = True
    res@tmXBMinorPerMajor = 8
    res@tmYLMinorPerMajor = 8
    res@gsnStringFont = "helvetica"
    res@gsnStringFontHeightF = 0.018
    res@tiMainFont = "helvetica"
    res@tiMainFontHeightF  = 0.02
    res@tiMainOffsetXF = -0.18
    res@tiMainOffsetYF = -0.01
    res@tiYAxisFontHeightF = 0.02
    res@tiXAxisFontHeightF = 0.02
    res@gsnStringFontHeightF = 0.02
    res@tmXBLabelFontHeightF = 0.018
    res@tmYLLabelFontHeightF = 0.018
    res@tiYAxisString = "Power"
    res@tiXAxisString = "Period (Months)"
    res@tmXBMajorOutwardLengthF = 0.006
    res@tmYLMajorOutwardLengthF = 0.006
    res@tmXBMinorOutwardLengthF = 0.003
    res@tmYLMinorOutwardLengthF = 0.003
    res@tmXBMajorLengthF = 0.006
    res@tmYLMajorLengthF = 0.006
    res@tmXBMinorLengthF = 0.003
    res@tmYLMinorLengthF = 0.003

    ; note a similiar resource is xyLineThicknessF, which will effect all
    ; the lines in the array.

    if (isatt(diag_script_info, "styleset")) then
    else
        diag_script_info@styleset = "CMIP5"
    end if

    lcolor = project_style(diag_script_info, "colors")
    ldash = project_style(diag_script_info, "dashes")
    lthick = project_style(diag_script_info, "thicks")

    res@tiMainString    = ""        ; add title
    res@xyLineThicknessF = lthick(0)
    res@xyLineColor  = lcolor(0)
    res@xyDashPattern = ldash(0)

    plot = gsn_csm_xy(wks, mem_period, mem_power(:, 0), res) ; create plot

    do i = 1, dim_MOD - 1
       res@xyLineThicknessF = lthick(i)
       res@xyLineColor  = lcolor(i)
       res@xyDashPattern = ldash(i)
       plot_tmp = gsn_csm_xy(wks, mem_period, mem_power(:, i), res) ; create plot
       overlay(plot, plot_tmp)
       delete(plot_tmp)
    end do

    res_lines                   = True        ; polyline mods desired
    res_lines@gsLineColor       = "gray"   ; line color
    res_lines@gsLineThicknessF  = 2.           ; line thicker
    res_lines@gsLineDashPattern = 1.   ; dash pattern

    txres               = True
    txres@txFontColor   = 5
    txres@txFontHeightF = 0.016            ; Set the font height
    txres@txAngleF      = 90.

    xx = (/24.,24./)
    yy = (/0.001,1000./)
    dum0 = gsn_add_polyline(wks, plot, xx,yy, res_lines)
    text0 = gsn_add_text(wks, plot, "2 Years", 27., 200., txres)

    xx = (/84.,84./)
    dum1 = gsn_add_polyline(wks, plot, xx, yy, res_lines)
    text1 = gsn_add_text(wks, plot, "7 Years", 97., 200., txres)
 
    ;***********************************************
    ; legend resources
    ;***********************************************

    lgres                    = True
    lgres@lgDashIndexes      = ldash
    lgres@lgLineColors       = lcolor
    lgres@lgMonoLineThickness = False
    lgres@lgLineThicknesses  = lthick       ; legend line thickness = contour line thickness
    lgres@lgItemType         = "Lines"      ; show lines only (default)
    lgres@lgLabelFontHeightF = .05          ; set the legend label font thickness
    lgres@vpWidthF           = 0.25         ; width of legend (NDC)
    lgres@vpHeightF          = 0.8 * tofloat(dim_MOD) / 30.0  ; height of legend (NDC)
    lgres@lgPerimColor       = "gray"       ; draw the box perimeter in orange
    lgres@lgPerimThicknessF  = 1.0          ; thicken the box perimeter
    labels                   = models@name
    lbid = gsn_create_legend(wks, dim_MOD, labels, lgres)         ; create legend

    amres = True
    amres@amParallelPosF   = 0.75           ; move legend to the right
    amres@amOrthogonalPosF = 0.22           ; move the legend down
    annoid1 = gsn_add_annotation(plot, lbid, amres)            ; add legend to plot

    resP = True
    resP@gsnMaximize = True
    resP@gsnPaperOrientation = "portrait"
    resP@gsnPaperMargin =  0.8
    resP@gsnPanelYF = (/.9,.49/)      ; Adjust bottommost plots.

    gsn_panel(wks,plot,(/2,1/),resP)

    ; ==================================================================================

    delete(mem_power)
    delete(mem_period)

    ; add meta data to plot (for reporting)

    caption = "Maximum entropy power spectra of variable " + variables(0) + " averaged over the NINO3 region" \
              + " (5N to 5S, 150W to 90W), similar to IPCC AR5, fig. 9.35."
    id = diag_script + "_" + variables(0)
    contrib_authors = "A_laue_ax"

    ESMValMD(outfile+"."+file_type, alltags, caption, id, variables(0), models@name, climofiles, diag_script, contrib_authors)

end

