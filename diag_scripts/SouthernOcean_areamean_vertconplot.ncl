;;#############################################################################
;; SouthernOcean_areamean_vertconplot.ncl
;; Author: Bert van Ulft (KNMI, the Netherlands)
;; EMBRACE project
;;#############################################################################
;; Description
;;    Calculate average over a subdomain and create a hovmoller like plot
;;
;; Required diag_script_info attributes (diagnostic specific)
;;    ref_model   : reference model as specified by annots
;;    styleset    : plot style
;;    West        : Western limit for area
;;    East        : Eastern limit for area
;;    South       : Southern limit for area
;;    North       : Northern limit for area
;;    Region      : name for selected region
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    grid_min    : min for plot scale
;;    grid_max    : max for plot scale
;;    grid_step   : step size for plot scale
;;    colormap    : e.g., WhiteBlueGreenYellowRed, rainbow
;;    dbg_plots   : create plots for all intermediate steps
;;
;; Caveats
;;    It is assumed masks vary per level, but are constant in time
;;    Year of reference must be 0 or same as other models, months must match
;;    days, hours, minutes, seconds are not checked
;;
;; Modification history
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;;    20140416-A_vanu_be: written.
;;
;;#############################################################################

;; Load addition ncl code
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;; Load Temporary file with current namelist configuration
load "./interface_data/ncl.interface"

;; Load auxiliary NCL routines
load "./diag_scripts/lib/ncl/ensemble.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/aux/SAMonsoon/SAMonsoon.ncl"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./plot_scripts/ncl/contourplot.ncl"
load "./plot_scripts/ncl/contour_maps.ncl"

;;#############################################################################
;; Collect information
;;#############################################################################
begin
    nl = str_get_nl()
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(diag_script, "", 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    write_references(diag_script,           \  ; script name
                     "A_vanu_be",           \  ; authors
                     (/""/),                \  ; contributors
                     "D_0000",              \  ; diag_references
                     "",                    \  ; obs_references
                     (/"P_embrace"/))          ; proj_acknowledgements

    do imod = 0, dim_MOD - 1
        data = read_data(imod, var0, field_type0)
    end do

    ;; Check field type
    if (field_type0 .ne. "T3M") then
        error_msg("f", diag_script, "", "can't use field of type " + \
                  field_type0 + ".")
    end if
end

;;#############################################################################
;; Load var_info
;;#############################################################################
begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

;;#############################################################################
;; Prepare for regridding and plotting
;;#############################################################################
begin
    ;; Prepare annotation
    annots = project_style(diag_script_info, "annots")

    ;; Reference model index
    idx_ref = get_ref_model_idx(models, var_attr_ref(0))
    diag_script_info@ref_model = models@name(idx_ref)

    ;; Check if reference model is available in annots
    if (.not. isatt(diag_script_info, "ref_model")) then
        error_msg("f", diag_script, "", "missing " + \
                  "diag_script_info@ref_model " + \
                  "-> define in config file " + diag_script_cfg + ".")
    end if
    iref = ind(annots .eq. diag_script_info@ref_model)

    ;; if ref_model not properly defined -> issue a list of valid strings
    if (ismissing(iref)) then
        errstr = "diag_script_info@ref_model = " + \
                 diag_script_info@ref_model + " is not valid." + nl + \
                 "Please use one of the following strings in " + \
                 "config file " + diag_script_cfg + ":" + nl + \
                 str_join(str_insert(annots, \
                          "diag_script_info@ref_model = ", 0), nl)
        error_msg("f", diag_script, "", errstr)
    else if (dimsizes(iref) .ne. 1) then  ; Ref_model is ambiguous
        errstr = "diag_script_info@ref_model = " + \
                 diag_script_info@ref_model + " is ambiguous." + nl + \
                 "Please use a unique string from the following " + \
                 "list in config file " + diag_script_cfg + ":" + nl + \
                 str_join(str_insert(annots, \
                          "diag_script_info@ref_model = ", 0), nl)
        error_msg("f", diag_script, "", errstr)
    else  ; Put reference model in first position to
          ; have it available for subsequent comparisons
        nannots = dimsizes(annots)
        ref_first = ispan(0, nannots - 1, 1)  ; Only valid if iref = 0
        if (iref .ne. 0) then                 ; reorder otherwise
            ref_first(0) = iref
            ref_first(1:iref) = ispan(0, iref - 1, 1)
            if (iref .ne. nannots - 1) then
                ref_first(iref + 1:nannots - 1) = \
                    ispan(iref + 1, nannots - 1, 1)
            end if
        end if
    end if
    end if

    ;; For regridding: define domain, grid size
    rgWest    =   0.5d
    rgEast    = 359.5d
    rgSouth   = -89.5d
    rgNorth   =  89.5d
    dLon      =   1.0
    dLat      =   1.0
    nLon      = doubletoint(round((rgEast - rgWest) / dLon, 0) + 1)
    nLat      = doubletoint(round((rgNorth - rgSouth) / dLat, 0) + 1)
    if (rgWest + (nLon - 1) * dLon .ne. rgEast .or. \
        rgSouth + (nLat - 1) * dLat .ne. rgNorth) then
        error_msg("f", diag_script, "", "definition of domain for " + \
                  "regridding is inconsistent.")
    end if
    strDomain = rgWest + "to" + rgEast + "_" + \
                rgSouth + "to" + rgNorth + "_" + \
                dLat + "x" + dLon

    ;; Must open the workstation before the loop
    name_string = var0
    wks = get_wks(0, diag_script, name_string)

    grid_min    = -0.14  ; default min in diff plot
    grid_max    =  0.14  ; default max in diff plot
    grid_step   =  0.02  ; default step size in diff plot
    colormap    = "temp_diff_18lev"
    dbg_plots   = False
    if (isatt(diag_script_info, "grid_min")) then
        grid_min = tofloat(diag_script_info@grid_min)
    end if
    if (isatt(diag_script_info, "grid_max")) then
        grid_max = tofloat(diag_script_info@grid_max)
    end if
    if (isatt(diag_script_info, "grid_step")) then
        grid_step = tofloat(diag_script_info@grid_step)
    end if
    if (isatt(diag_script_info, "colormap")) then
        colormap = diag_script_info@colormap
        gsn_define_colormap(wks, colormap)
    end if
    if (isatt(diag_script_info, "dbg_plots")) then
        dbg_plots = diag_script_info@dbg_plots
    end if
end

;;#############################################################################
;; loop over models (starting with reference) and:
;; - read data
;; - regrid
;; - calculate area mean
;; - interpolate to reference vertical profile
;; - calculate difference with reference
;; - plot
;;#############################################################################
begin
    do iannots = 0, dimsizes(annots) - 1
        ;; Make sure to start with reference model
        imod = ref_first(iannots)
        info_output("processing " + models@name(imod) \
                                  + "_" + models@case_name(imod) \
                                  + "_" + models@ensemble(imod), \
                    verbosity, 5)

        ;; Read all data values for this model
        ;; (interface_scripts/data_handling.ncl)
        var_raw = extract_data(imod, data, -1, 0, 0)

        ;; Get 2d lat & lon coordinates
        var = make_latlon2D(var0, imod, data, var_raw)
        delete(var_raw)

        ;; part of file names for regridding grid and weights
        stem      = interface_get_inpaths(imod) + "/ESMF_regrid_" + \
                    models@project(imod) + "_" + \
                    models@name(imod) + "_"
        ending    = var0 + "_" + strDomain + "_lev"

        ;; array to store regridded data
        var_regrid = new((/dimsizes(var&time), dimsizes(var&lev), nLat, nLon/), float)
        var_regrid!1 = "lev"
        var_regrid&lev = var&lev

        ;; Regridding per level because ESMF mask can only be 2D
        ;; ESMF_regrid may issue a warning
        ;; (e.g. if Dst is not covered by Src grid)
        ;; -> set "exit_on_warning" to False (in main nml) to continue to
        ;; the next diag
        ilev_min = 0
        do ilev = ilev_min, dimsizes(var&lev) - 1
            info_output(diag_script + ": regridding level " + ilev, \
                        verbosity, 4)
            if (all(ismissing(var(0, ilev, :, :)))) then
                info_output("  all points missing, skipping regridding", \
                            verbosity, 4)
            else
                Opt                = True
                Opt@DstLLCorner    = (/rgSouth, rgWest /)
                Opt@DstURCorner    = (/rgNorth, rgEast /)
                Opt@DstGridType    = dLat + "x" + dLon  ; Destination grid
                Opt@ForceOverwrite = True
                Opt@PrintTimings   = True
                Opt@Debug          = False  ; setting Debug to True causes
                                            ; ESMF to hang after some levels

                ;; Create a mask to filter out missing values
                Opt@SrcMask2D      = where(ismissing(var(0, ilev, :, :)), 0, 1)

                ;; "bilinear" ; "patch" ; patch may be better but slower
                Opt@InterpMethod   = "patch"

                ;; Setting Opt@SrcRegional=True avoids artefacts
                ;; from incomplete Src grids
                Opt@SrcRegional    = True  ; Activate if only part of the
                                           ; sphere is given

                ;; Save regrid files for later use: careful with naming!
                Opt@SrcFileName    = stem + "Src_" + ending + ilev + ".nc"
                Opt@DstFileName    = stem + "Dst_" + ending + ilev + ".nc"
                Opt@WgtFileName    = stem + "Wgt_" + ending + ilev + ".nc"

                ;; reuse regrid files if present
                Opt@SkipSrcGen     = False
                Opt@SkipDstGen     = False
                Opt@SkipWgtGen     = False
                if (isfilepresent_esmval(Opt@SrcFileName)) then
                    Opt@SkipSrcGen = True
                end if
                if (isfilepresent_esmval(Opt@DstFileName)) then
                    Opt@SkipDstGen = True
                end if
                if (isfilepresent_esmval(Opt@WgtFileName)) then
                    Opt@SkipWgtGen = True
                    temp = ESMF_regrid_with_weights(var(:, ilev, :, :),\
                                                    Opt@WgtFileName,\
                                                    Opt)
                else
                    temp = ESMF_regrid(var(:, ilev, :, :), Opt)
                end if
                var_regrid(:, ilev, :, :) = temp
                delete(temp)
                delete(Opt)
            end if
        end do

        ;; Create contour plot for surface to check regridding went ok
        if (dbg_plots) then
            N = dimsizes(var_regrid)
            data1 = time_operations(var_regrid(:, N(1) - 1, :, :), \
                                    -1, -1, "average", "annualclim", \
                                    True)
            data1@var                = var0
            data1@diag_script        = (/diag_script/)
            data1@res_gsnMaximize    = True      ; Use full page for the plot
            data1@res_cnFillOn       = True      ; Color plot desired
            data1@res_cnLineLabelsOn = False     ; Contour labels
            data1@res_cnLinesOn      = False     ; Contour lines
            data1@res_cnFillPalette  = colormap  ; Set colormap
            map = contour_map(wks, data1, \
                              "toplevel regridded" + \
                              models@name(imod) + " " + variables(0))
            delete(data1)
        end if

        ;; Calculate area weighted average
        ;; (diag_scripts/lib/ncl/latlon.ncl)
        latmin    =  diag_script_info@South
        latmax    =  diag_script_info@North
        lonmin    =  diag_script_info@West
        lonmax    =  diag_script_info@East
        var_prof = area_operations(var_regrid, latmin, latmax, lonmin, \
                                   lonmax, "average", True)
        ;; plot area weighted average
        if (dbg_plots) then
            data1 = var_prof
            data1@var                = var0
            data1@diag_script        = (/diag_script/)
            data1@res_gsnMaximize    = True      ; Use full page for the plot
            data1@res_cnFillOn       = True      ; Color plot desired
            data1@res_cnLineLabelsOn = False     ; Contour labels
            data1@res_cnLinesOn      = False     ; Contour lines
            data1@res_cnFillPalette  = colormap  ; Set colormap
            map = contourplot(wks, data1, \
                              models@name(imod) + " " + variables(0) + \
                              " for " + diag_script_info@region)
            delete(data1)
        end if

        ;; Just keep ref_model OR plot difference to ref_model
        if (iannots .eq. 0) then  ; ref_model -> keep, but don't plot
            var_prof_ref = var_prof
            copy_VarMeta(var_prof, var_prof_ref)
            ntime_ref = dimsizes(var_prof_ref&time)
            utime_ref = cd_calendar(var_prof_ref&time, 0)
            yy_ref    = tointeger(utime_ref(:, 0))
            mm_ref    = tointeger(utime_ref(:, 1))
        else  ; plot difference to ref_model
            ;; interpolate to reference vertical profile
            var_prof_vint = linint1_n_Wrap(var_prof&lev,\
                                           var_prof,\
                                           False,\
                                           var_prof_ref&lev,\
                                           0,\
                                           1)

            ;; subtract reference, time dimensions must have the same length OR
            ;; time of var must be a multiple of that of the reference
            var_prof_diff = var_prof_vint
            copy_VarMeta(var_prof_vint, var_prof_diff)
            ntime_var  = dimsizes(var_prof_vint&time)
            utime_var  = cd_calendar(var_prof_vint&time, 0)
            yy_var     = tointeger(utime_var(:, 0))
            mm_var     = tointeger(utime_var(:, 1))
            monfreq    = 3                    ; show label every monfreq months
            tmYLMjVals = new(floattoint(ntime_var / monfreq) + 1, double)
            tmYLMnVals = new(ntime_var, double)
            tmYLLbls   = new(floattoint(ntime_var / monfreq) + 1, string, "")
            if (mod(ntime_var, ntime_ref) .eq. 0) then
                ; loop over time steps
                itmYL = 0
                do itime = 0, ntime_var - 1, ntime_ref
                    iend = itime + ntime_ref - 1

                    ; check dates match and create y label array
                    do itc = 0, ntime_ref - 1
                        if (yy_ref(itc) .ne. 0 .and. \
                            yy_ref(itc) .ne. yy_var(itime + itc)) then
                            error_msg("f", diag_script, "", \
                                      "year of reference (" + \
                                      yy_ref(itc) + ") is not 0 " + \
                                      "and is not equal to year of " + \
                                      "the current model (" + \
                                      yy_var(itime + itc) + ").")
                        end if
                        if (mm_ref(itc) .ne. mm_var(itime + itc)) then
                            error_msg("f", diag_script, "", \
                                      "month of reference (" + \
                                      mm_ref(itc) + ") is not " + \
                                      "equal to month of current " + \
                                      "model (" + \
                                      mm_var(itime + itc) + ").")
                        end if
                        if (mod(mm_var(itime + itc) - 1, monfreq) .eq. 0) then
                            tmYLMjVals(itmYL) = var_prof_vint&time(itime + itc)
                            tmYLLbls(itmYL)   = \
                                sprinti("%0.4i", yy_var(itime + itc)) + \
                                "-" + \
                                sprinti("%0.2i", mm_var(itime + itc))
                            itmYL = itmYL + 1
                        end if
                        tmYLMnVals(itime + itc) = \
                            var_prof_vint&time(itime + itc)
                    end do

                    ; subtract
                    var_prof_diff(itime:iend, :) = \
                        var_prof_vint(itime:iend, :) - var_prof_ref

                end do
            else
                error_msg("f", diag_script, "", \
                          "unsupported combination of number of " + \
                          "time steps. Number of time steps of " + \
                          "current model (" + models@name(imod) + \
                          ") must be a multiple of that of the " + \
                          "reference model.")
            end if

            ;; plot vint
            if (dbg_plots) then
                data1 = var_prof_vint
                data1@var                = var0 + "vint"
                data1@diag_script        = (/diag_script/)
                data1@res_gsnMaximize    = True         ; Full page for plot
                data1@res_cnFillOn       = True         ; Color plot desired
                data1@res_cnLineLabelsOn = False        ; Contour labels
                data1@res_cnLinesOn      = False        ; Contour lines
                data1@res_cnFillPalette  = colormap     ; Set colormap
                data1@res_tmYLMode       = "Explicit"   ; time on Y explicit
                data1@res_tmYLValues     = tmYLMjVals   ; time on y axis
                data1@res_tmYLLabels     = tmYLLbls     ; time on Y labels
                data1@res_tmYLLabelFontHeightF = 0.015  ; font size of y axis
                map = contourplot(wks, data1, \
                                  models@name(imod) + " " + \
                                  variables(0) + " for " + \
                                  diag_script_info@region)
                delete(data1)
            end if

            ;; plot diff; put time on x, depth on y, with surface on top
            data1 = transpose(var_prof_diff)
            data1@var = "~F33~D~F~ " + var0 + \
                " (" + var_prof_diff@units + ")"
            data1@diag_script = (/diag_script/)
            data1@res_tiMainString   = "~F33~D~F~ of " + \
                var@long_name + " for " + diag_script_info@region
            data1@res_tiXAxisString  = "date (yyyy-mm)"
            data1@res_tiYAxisString  = "ocean depth (m)"  ; assume metres
            data1@res_trYReverse     = True         ; reverse y-axis
            data1@res_gsnLeftString  = annots(imod) + " - " + \
                diag_script_info@ref_model
            data1@res_gsnMaximize    = True         ; Use full page for plot
            data1@res_lbOrientation  = "Vertical"
            data1@res_cnFillOn       = True         ; Color plot desired
            data1@res_cnLineLabelsOn = False        ; Contour labels
            data1@res_cnLinesOn      = False        ; Contour lines
            data1@res_cnFillPalette  = colormap     ; Set colormap
            data1@res_cnLevelSelectionMode = "ManualLevels"  ; "ExplicitLevels"
            data1@res_cnMinLevelValF = grid_min     ; min of scale
            data1@res_cnMaxLevelValF = grid_max     ; max of scale
            data1@res_cnLevelSpacingF = grid_step   ; step of scale
            data1@res_tmXBMode       = "Explicit"   ; time on Y explicit
            data1@res_tmXBValues     = tmYLMjVals   ; time on y axis
            data1@res_tmXBLabels     = tmYLLbls     ; time on Y explicit labels
            data1@res_tmXBMinorValues = tmYLMnVals  ; minor ticks
            data1@res_tmXBLabelFontHeightF = 0.014  ; font size of y axis
            data1@res_tmXBLabelAngleF = 90          ; rotate labels
            ;data1@res_tmXBLabelDeltaF = 2.0        ; move labels away
            data1@res_tmXBLabelJust  = "CenterRight"  ; position
            data1@res_tmXMajorGrid   = True         ; show major grid lines
            data1@res_tmXMajorGridThicknessF = 1.0  ; thickness of major
            map = contourplot(wks, data1, "dummy_for_var")
            delete(data1)
            delete([/var_prof_vint, var_prof_diff/])
            delete([/utime_var, yy_var, mm_var/])
            delete([/tmYLMjVals, tmYLMnVals, tmYLLbls/])
        end if

        ;; clean up
        delete([/var, var_regrid, var_prof/])
        if (isfilepresent_esmval("./PET0.RegridWeightGen.Log")) then
            system("rm ./PET0.RegridWeightGen.Log")
        end if
    end do
end

;;#########################################################################
;; References Information
;;#########################################################################
begin
    leave_msg(diag_script, "", 4)
end
