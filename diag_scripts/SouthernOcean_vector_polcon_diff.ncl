;;#############################################################################
;; SouthernOcean_vector_polcon_diff.ncl
;;#############################################################################
;; Author: Bert van Ulft (KNMI, the Netherlands)
;; EMBRACE project
;;#############################################################################
;; Description
;;    Transfers all data to a common grid
;;    Draws absolute difference between two data sets for var(0)
;;    Overlays vectors for two data sets from var(1) and var(2)
;;
;; Required diag_script_info attributes (diagnostic specific)
;;    styleset    : plot style
;;    month       : selected month or annual mean
;;    range_option: use models own or intersection of time range (0|1)
;;    region      : name for the region (Arctic|Antarctic)
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    max_lat     : for Antarctic, maximum latitude to plot
;;    min_lat     : for Arctic, minimum latitude to plot
;;    max_hori    : max allowed nr of columns on a panel page (horizontal)
;;    max_vert    : max allowed nr of rows on a panel page (vertical)
;;    grid_min    : min for plot scale
;;    grid_max    : max for plot scale
;;    grid_step   : step size for plot scale
;;    colormap    : e.g., WhiteBlueGreenYellowRed, rainbow
;;    showunits   : display units in figure title
;;    vcmagnitude : magnitude, larger value, shorter arrows
;;    vcmindist   : controls density of vectors, larger value, less arrows
;;    vccolormod  : colour of vectors for current model
;;    vccolorref  : colour of vectors for reference model
;;
;; Caveats
;;    Only for Antarctic/Southern hemisphere
;;
;; Modification history
;;    * 20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                          of the code
;;    * 20141014-A_vanu_be: written, largely based on SeaIce_polcon_diff
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

load "./plot_scripts/ncl/aux_plotting.ncl"
load "./plot_scripts/ncl/vector_scalar_map_polar.ncl"

;;#############################################################################
;; Collect information
;;#############################################################################
begin
    nl = str_get_nl()
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(diag_script, "", 4)
    diagsc_base = systemfunc("basename " + diag_script)

    ;; print info, part 1
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script, verbosity, 1)

    ;; check that we have three variables
    dim_VAR = dimsizes(variables)
    if (dim_VAR .ne. 3) then
        error_msg("f", diag_script, "", "must have 3 input " + \
                  "variables , (e.g. tos, u and v) but got " + \
                  dim_VAR + ".")
    end if
    varc = variables(0)     ; variable for contours
    varu = variables(1)     ; u(-like) variable, used with v for vectors
    varv = variables(2)     ; v(-like) variable, used with u for vectors

    write_references(diag_script,           \  ; script name
                     "A_vanu_be",           \  ; authors
                     (/""/),                \  ; contributors
                     "D_0000",              \  ; diag_references
                     "",                    \  ; obs_references
                     (/"P_embrace"/))          ; proj_references

    ;; print info, part 2
    info_output("    contours for: " + varc + "(model) - " + \
                varc + "(ref_model)", verbosity, 1)
    info_output("    vectors from: " + varu + " & " + varv + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    ;; Check field types
    field_type0 = field_types(0)
    if (field_type0 .ne. "T2Ms" .and. field_type0 .ne. "TO2Ms") then
        error_msg("f", diag_script, "", "can't use field of type " + \
                  field_type0 + " for diag_script " + diag_script + ".")
    end if
    if (.not. all(field_types(1:) .eq. field_type0)) then
        error_msg("f", diag_script, "", "field types of all " + \
                  "variables must be the same for diag_script " + \
                  diag_script + " but they are not: " + nl + \
                  str_join(str_insert(variables + ": " + field_types, \
                           "field_type for var ", 0), nl))
    end if

    ;; Get environment variables
    file_type = getenv("ESMValTool_output_file_type")
    if(ismissing(file_type)) then
        file_type = "PS"
    end if

    ;; Load var_info
    vardeffile = "interface_data/" + varc + "_info.tmp"
    loadscript(vardeffile)
    vardeffile = "interface_data/" + varu + "_info.tmp"
    loadscript(vardeffile)
    vardeffile = "interface_data/" + varv + "_info.tmp"
    loadscript(vardeffile)

    ;; check all required diag_script_info attributes have been set
    required_attributes = (/"styleset",\
                            "month",\
                            "range_option",\
                            "region"/)
    exit_if_missing_atts(diag_script_info, required_attributes)

    ;; Only allow Antarctic as region
    if (diag_script_info@region .ne. "Antarctic") then
        error_msg("f", diag_script, "", "currently the only " + \
                  "supported region is Antarctic, because the u " + \
                  "and v vectors are relative to the model grid " + \
                  "and are not rotated to the lat-lon grid. " + \
                  "On the Southern hemisphere the grid is " + \
                  "approximately a regular lat-lon grid.")
    end if
end

;;#############################################################################
;; Prepare for analysis and plotting
;;#############################################################################
begin
    ;; Prepare annotation
    annots = project_style(diag_script_info, "annots")
    smonths = (/"January", "February", "March", "April", "May", "June", \
                "July", "August", "September", "October", "November", \
                "December"/)

    ;; Time indication in output file name
    if (diag_script_info@month .eq. "A") then
        i2month = "annual"
    else
        i2month = smonths(toint(diag_script_info@month) - 1)
    end if

    ;; Define output file name
    name_string = varc + "_" + varu + "_" + varv + "_" + \
                  diag_script_info@region + "_" + i2month + "_"
    if (diag_script_info@range_option .eq. 1) then
        name_string = name_string + \
                      sprinti("%0.4i", max(toint(models@start_year))) + \
                      "-" + \
                      sprinti("%0.4i", min(toint(models@end_year)))
    else if (diag_script_info@range_option .eq. 0) then
        name_string = name_string + \
                      sprinti("%0.4i", min(toint(models@start_year))) + \
                      "-" + \
                      sprinti("%0.4i", max(toint(models@end_year)))
    else
        error_msg("f", diag_script, "", "range_option not set " + \
                  "correctly - must be 1 or 0")
    end if
    end if

    ;; Must open the workstation before the loop
    wks = get_wks(0, diag_script, name_string)

    ;; define default color scale
    grid_min    = -0.2        ; default min in diff plot
    grid_max    =  0.2        ; default max in diff plot
    grid_step   =  0.02       ; default step size in diff plot
    colormap    = "ncl_default"
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
    end if
    gsn_define_colormap(wks, colormap)
    nsteps = round((grid_max - grid_min) / grid_step, 3) + 1
    cnLevels = fspan(grid_min, grid_max, nsteps)

    ;;  Properties of vectors
    vcmagnitude =  0.2        ; magnitude of vectors
    vcmindist   =  0.04       ; controls density of vectors
    vccolormod  =  "green4"   ; color of vectors for current model
    vccolorref  =  "black"    ; color of vectors for reference model
    if (isatt(diag_script_info, "vcmagnitude")) then
        vcmagnitude = diag_script_info@vcmagnitude
    end if
    if (isatt(diag_script_info, "vcmindist")) then
        vcmindist = diag_script_info@vcmindist
    end if
    if (isatt(diag_script_info, "vccolormod")) then
        vccolormod = diag_script_info@vccolormod
    end if
    if (isatt(diag_script_info, "vccolorref")) then
        vccolorref = diag_script_info@vccolorref
    end if

    ;; Panel parameters
    nvert = 1  ; Default
    nhori = 1  ; Default
    if (isatt(diag_script_info, "max_vert")) then
        ;; Maximum allowed number of plots per page (vertical)
        nvert = toint(diag_script_info@max_vert)
    end if
    if (isatt(diag_script_info, "max_hori")) then
        ;; Maximum allowed number of plots per page (horizontal)
        nhori = toint(diag_script_info@max_hori)
    end if
    npages = toint(ceil(tofloat(dimsizes(annots) - 1) / \
                        tofloat(nvert * nhori)))

    ;; Define array for collecting plots
    if (npages .gt. 1) then
        ;; Size of all stamps is identical only if every page contains
        ;; the same number of plots
        plots = new(npages * nvert * nhori, graphic)
    else
        plots = new(dimsizes(annots) - 1, graphic)
    end if

    ;; Reference model index
    idx_ref = get_ref_model_idx(models, var_attr_ref(0))
    diag_script_info@ref_model = models@name(idx_ref)

    ;; ref_model not properly defined -> issue a list of valid strings
    iref = ind(annots .eq. diag_script_info@ref_model)
    if (ismissing(iref)) then
        errstr = "diag_script_info@ref_model = " + \
                 diag_script_info@ref_model + " is not valid." + nl + \
                 "Please use one of the following strings in " + \
                 "the xml file " + diag_script_cfg + ":" + nl + \
                 str_join(str_insert(annots, \
                          "diag_script_info@ref_model = ", 0), nl)
        error_msg("f", diag_script, "", errstr)
    else if (dimsizes(iref).ne.1) then  ; Ref_model is ambiguous
        errstr = "diag_script_info@ref_model = " + \
                 diag_script_info@ref_model + " is ambiguous." + nl + \
                 "Please use a unique string from the following " + \
                 "list in the xml file " + diag_script_cfg + ":" + nl + \
                 str_join(str_insert(annots, \
                          "diag_script_info@ref_model = ", 0), nl)
        error_msg("f", diag_script, "", errstr)
    else
        ;; Put reference model in first position to have it available
        ;; for subsequent comparisons
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
end

;;#############################################################################
;; Read data
;;#############################################################################
begin
    dim_MOD = dimsizes(models@name)
    do idx_mod = 0, dim_MOD - 1
        datac = read_data(idx_mod, varc, field_type0)
        datau = read_data(idx_mod, varu, field_type0)
        datav = read_data(idx_mod, varv, field_type0)
    end do
end

;;#############################################################################
;; Start processing
;;#############################################################################
begin
    ;; Loop over models
    do iannots = 0, dimsizes(annots) - 1

        ;; Make sure to start with reference model
        imod = ref_first(iannots)

        info_output("processing " + annots(imod), verbosity, 2)
        ;; See ./interface_scripts/data_handling.ncl
        extrc_raw = extract_data(imod, datac, -1, 0, 0)
        extru_raw = extract_data(imod, datau, -1, 0, 0)
        extrv_raw = extract_data(imod, datav, -1, 0, 0)

        ;; store and check units
        if (iannots.eq.0) then
            units = extrc_raw@units
        end if
        if (extru_raw@units .ne. units .or. \
            extrv_raw@units .ne. units) then
            error_msg("f", diag_script, "", "units mismatch " + \
                      "between models and/or variables" + nl + \
                      "units      : " + units + nl + \
                      "extru@units: " + extru_raw@units + nl + \
                      "extrv@units: " + extrv_raw@units)
        end if

        ;; Show units in title (default False)
        strUnits = ""
        if (isatt(diag_script_info, "showunits")) then
            strUnits = " (" + units + ")"
        end if

        ;; Get 2d lat & lon coordinates for variables
        extrc = make_latlon2D(varc, imod, datac, extrc_raw)
        extru = make_latlon2D(varu, imod, datau, extru_raw)
        extrv = make_latlon2D(varv, imod, datav, extrv_raw)
        delete([/extrc_raw, extru_raw, extrv_raw/])

        ;; Time axis
        start_year = get_start_year(imod)
        end_year   = get_end_year(imod)
        if (diag_script_info@month.eq."A") then
            extrc_ann = time_operations(extrc, start_year, end_year, \
                                       "average", "annualclim", True)
            extru_ann = time_operations(extru, start_year, end_year, \
                                       "average", "annualclim", True)
            extrv_ann = time_operations(extrv, start_year, end_year, \
                                       "average", "annualclim", True)
        else
            extrc_ann = time_operations(extrc, start_year, end_year, \
                                       "average", diag_script_info@month, True)
            extru_ann = time_operations(extru, start_year, end_year, \
                                       "average", diag_script_info@month, True)
            extrv_ann = time_operations(extrv, start_year, end_year, \
                                       "average", diag_script_info@month, True)
        end if
        range = start_year + " - " + end_year
        delete([/extrc, extru, extrv/])
        extrc = extrc_ann
        extru = extru_ann
        extrv = extrv_ann
        copy_VarMeta(extrc_ann, extrc)
        copy_VarMeta(extru_ann, extru)
        copy_VarMeta(extrv_ann, extrv)
        delete([/extrc_ann, extru_ann, extrv_ann/])

        ;; Regridding of variables
        ;; ESMF_regrid may issue a warning
        ;; (e.g. if Dst is not covered Src grid)
        ;; -> set "exit_on_warning" to False (in main nml) to continue to the
        ;; next diag
        ;; general settings
        Opt                   = True
        Opt@DstLLCorner       = (/-89.75d,   0.75d /)
        Opt@DstURCorner       = (/ 89.75d, 359.75d /)
        Opt@SrcMask2D         = where(ismissing(extrc), 0, 1)
        Opt@DstGridType       = "1deg"      ; Destination grid
        Opt@ForceOverwrite    = True
        Opt@PrintTimings      = True
        Opt@Debug             = True
        Opt@InterpMethod      = "patch"  ; "bilinear"; "patch"; "conserve"
        Opt@SrcRegional       = True     ; Activate if only part of the sphere
                                         ; is given, avoids artefacts from
                                         ; incomplete Src grids
        stem   = interface_get_inpaths(imod) + "/ESMF_regrid_" + \
                 models@project(imod) + "_" + models@name(imod) + "_"
        ending = diag_script_info@region + "_" + \
                 diag_script_info@month + "_" + \
                 start_year + "-" + end_year + ".nc"

        ;; regrid contour variable
        ;; Save regrid files for later use: careful with naming!
        ;; reuse regrid files if present
        Opt@WgtFileName       = stem + "Wgt_" + varc + "_" + ending
        Opt@SrcFileName       = stem + "Src_" + varc + "_" + ending
        Opt@DstFileName       = stem + "Dst_" + varc + "_" + ending
        Opt@SkipSrcGen        = False
        Opt@SkipDstGen        = False
        Opt@SkipWgtGen        = False
        if (isfilepresent_esmval(Opt@SrcFileName)) then
            Opt@SkipSrcGen    = True
        end if
        if (isfilepresent_esmval(Opt@DstFileName)) then
            Opt@SkipDstGen    = True
        end if
        if (isfilepresent_esmval(Opt@WgtFileName)) then
            Opt@SkipWgtGen    = True
            extrc_regrid      = ESMF_regrid_with_weights(extrc, Opt@WgtFileName, Opt)
        else
            extrc_regrid      = ESMF_regrid(extrc, Opt)
        end if

        ;; regrid u variable
        ;; Save regrid files for later use: careful with naming!
        ;; reuse regrid files if present
        Opt@WgtFileName       = stem + "Wgt_" + varu + "_" + ending
        Opt@SrcFileName       = stem + "Src_" + varu + "_" + ending
        Opt@DstFileName       = stem + "Dst_" + varu + "_" + ending
        Opt@SkipSrcGen        = False
        Opt@SkipDstGen        = False
        Opt@SkipWgtGen        = False
        if (isfilepresent_esmval(Opt@SrcFileName)) then
            Opt@SkipSrcGen    = True
        end if
        if (isfilepresent_esmval(Opt@DstFileName)) then
            Opt@SkipDstGen    = True
        end if
        if (isfilepresent_esmval(Opt@WgtFileName)) then
            Opt@SkipWgtGen    = True
            extru_regrid      = ESMF_regrid_with_weights(extru, Opt@WgtFileName, Opt)
        else
            extru_regrid      = ESMF_regrid(extru, Opt)
        end if

        ;; regrid v variable
        ;; Save regrid files for later use: careful with naming!
        ;; reuse regrid files if present
        Opt@WgtFileName       = stem + "Wgt_" + varv + "_" + ending
        Opt@SrcFileName       = stem + "Src_" + varv + "_" + ending
        Opt@DstFileName       = stem + "Dst_" + varv + "_" + ending
        Opt@SkipSrcGen        = False
        Opt@SkipDstGen        = False
        Opt@SkipWgtGen        = False
        if (isfilepresent_esmval(Opt@SrcFileName)) then
            Opt@SkipSrcGen    = True
        end if
        if (isfilepresent_esmval(Opt@DstFileName)) then
            Opt@SkipDstGen    = True
        end if
        if (isfilepresent_esmval(Opt@WgtFileName)) then
            Opt@SkipWgtGen    = True
            extrv_regrid      = ESMF_regrid_with_weights(extrv, Opt@WgtFileName, Opt)
        else
            extrv_regrid      = ESMF_regrid(extrv, Opt)
        end if

        ;; Some attributes for area plot
        if (diag_script_info@month .eq. "A") then
            ;; Weighted annual mean
            extrc_regrid@units = "annual mean " + range
        else
            ;; Select a month
            units_str = \
                smonths(floattoint(stringtofloat(diag_script_info@month)) - 1)\
                + " mean " + range
            extrc_regrid@units = units_str
            extru_regrid@units = units_str
            extrv_regrid@units = units_str
        end if
        extrc_regrid@long_name = diag_script_info@region + " " + \
                                 extrc@long_name + strUnits
        extru_regrid@long_name = diag_script_info@region + " " + \
                                 extru@long_name
        extrv_regrid@long_name = diag_script_info@region + " " + \
                                 extrv@long_name

        ;; Clean up after regridding
        delete([/extrc, extru, extrv, Opt/])

        ;; Start plotting
        if (iannots .eq. 0) then
            ;; ref_model -> keep, plot vectors only
            extrc_ref = extrc_regrid
            extru_ref = extru_regrid
            extrv_ref = extrv_regrid
            copy_VarMeta(extrc_regrid, extrc_ref)
            copy_VarMeta(extru_regrid, extru_ref)
            copy_VarMeta(extrv_regrid, extrv_ref)
            delete([/extrc_regrid, extru_regrid, extrv_regrid/])
        else
            ;; other model, plot difference to ref_model and vectors
            ;; Enforce identical _FillValue for extrc_ref & extrc_regrid
            if (extrc_regrid@_FillValue .ne. extrc_ref@_FillValue) then
                extrc_regrid = where(ismissing(extrc_regrid), \
                                     extrc_ref@_FillValue, extrc_regrid)
                extrc_regrid@_FillValue = extrc_ref@_FillValue
            end if

            ;; Difference to ref_model
            extrc = extrc_regrid - extrc_ref
            copy_VarMeta(extrc_regrid, extrc)
            extrc = where(ismissing(extrc_ref), extrc_ref@_FillValue, extrc)
            extrc = where(ismissing(extrc_regrid), extrc_ref@_FillValue, extrc)

            extrc@long_name = "~F33~D~F~ of " + extrc@long_name
            delete(extrc_regrid)

            ;; plot vectors for this model and contours of diff with reference
            ;; Define diag-specific resources to be passed to plot routine
            ;; common resources
            res                     = True
            res@vcRefMagnitudeF     = vcmagnitude         ; vector magnitude
            res@vcRefLengthF        = 0.05                ; size of ref vector
            res@vcMinDistanceF      = vcmindist           ; controls density of vectors
            res@vcGlyphStyle        = "CurlyVector"       ; turns on curly vectors
            res@vcRefAnnoArrowUseVecColor = True          ; use colour of vector for reference vector
            res@vcRefAnnoZone       = 2                   ; specifies location of reference vector
            res@vcRefAnnoString2    = res@vcRefMagnitudeF + " " + units  ; text below reference vector
            res@vcRefAnnoString2On  = True                ; show second string
            res@vcRefAnnoPerimOn    = False               ; no box around reference vector

            ;; resources for contour plot of difference and vectors of current model
            resmod                  = res
            resmod@tiMainString     = annots(imod) + " - " + \
                                      diag_script_info@ref_model
            resmod@gsnScalarContour = True                ; vectors over scalars
            resmod@gsnSpreadColors  = True                ; use full colormap
            resmod@cnLevelSelectionMode = "ExplicitLevels"  ; define all levels
            resmod@cnLevels         = cnLevels            ; explicit levels
            resmod@cnLinesOn        = False               ; no contour lines
            resmod@cnLineLabelsOn   = False               ; no contour line labels
            resmod@vcLineArrowColor = vccolormod          ; colour of vectors
            resmod@vcRefAnnoString1 = annots(imod)        ; text above reference vector

            ;; resources for plot of vectors for reference
            resref                  = res
            resref@gsnLeftString    = ""                  ; don't show titles
            resref@gsnRightString   = ""                  ; don't show titles
            resref@gsnDraw          = False               ; don't draw it
            resref@gsnFrame         = False               ; Don't advance frame.
            resref@vcLineArrowColor = vccolorref          ; colour of vectors
            resref@vcRefAnnoString1 = diag_script_info@ref_model  ; text above reference vector

            ;; make plots
            plot = vector_scalar_map_polar(wks,\
                                           extru_regrid,\
                                           extrv_regrid,\
                                           extrc, \
                                           resmod, \
                                           varu,\
                                           varv,\
                                           varc)

            plot_ref = gsn_csm_vector(wks,\
                                      extru_ref,\
                                      extrv_ref,\
                                      resref)

            ;; combine plots
            overlay(plot, plot_ref)
            plots(iannots - 1)      = plot

            ;; Clean up
            delete([/extrc, extru_regrid, extrv_regrid, res, resmod, resref/])

        end if  ; ref_model?

        if (isfilepresent_esmval("./PET0.RegridWeightGen.Log")) then
            system("rm ./PET0.RegridWeightGen.Log")
        end if
    end do
end

;;#############################################################################
;; Draw the panel
;;#############################################################################
begin
    pres = True
    pres@lbLabelStride = 2
    pres@gsnPanelLabelBar = True  ; turn on panel label bar

    outfile = panelling(wks, plots, nvert, nhori, pres)

    info_output(" gv " + outfile + " &", verbosity, 1)
end

;;#############################################################################
;; References Information
;;#############################################################################
begin
    leave_msg(diag_script, "", 4)
end
