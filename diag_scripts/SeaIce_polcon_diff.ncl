;;#############################################################################
;; POLAR-STEREOGRAPHIC CONTOUR PLOTS OF SEA ICE AREA DIFFERENCE
;; Author: Daniel Senftleben (DLR, Germany)
;; ESMVal project
;;#############################################################################
;;
;; Description
;;    Transfers all data to a common 1°x1° Gaussian grid
;;    Draws contour maps of absolute difference in sea ice concentration between
;;      all models and the specified reference model (e.g., observations)
;;
;; Required variable_info attributes (variable specific)
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    month: selected month or annual mean
;;    region: region to be plotted ("Arctic" or "Antarctic")
;;    contour_extent: optionally draw a red contour line at 0.15
;;    max_lat, min_lat: plot ranges for SH & NH
;;    ref_model: reference model as specified by annots
;;    range_option: select specific or intersecting time range
;;    dst_grid: path to destination grid description file for 'cdo remapdis'
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    max_vert, max_hori: allowed dimensions of a panel page
;;                        (number of plots per column/row)
;;    PanelTop: tune to get full title of uppermost row (1 = no top margin)
;;    grid_min, grid_max, grid_step, grid_center: contour parameters
;;    colormap: NCL colour map name
;;    PanelLabelBar: use single label bar per page
;;    showunits: display units in figure title
;;    multi_model_mean: only plot multi-model-mean minus ref_model
;;
;; Caveats
;;    No masking of the apparent line of missing data between the two grid
;;      poles yet (for a workaround, see
;;      http://www.ncl.ucar.edu/Applications/Scripts/ice_4.ncl).
;;    When using multi-model-mean option, make sure to exclude all unwanted 
;;      files (e.g., observations that are not the reference) from the namelist.
;;
;; Modification history
;;    20170512-A_senf_da: Added reporting
;;    20170120-A_senf_da: Moved regridding to cdo_remapdis (latlon.ncl)
;;    20160713-A_senf_da: Included multi-model-mean option
;;    20151027-A_laue_ax: Moved call to 'write_references' to the beginning
;;                        of the code
;;    20150605-A_vanu_be: Additional attributes for more formatting control
;;    20150325-A_laue_ax: Modified reference tags used for acknowledgements
;;                        (projects, observations, etc)
;;    20141023-A_vanu_be: Externalized make_latlon2D + var specific weights
;;    20140623-A_senf_da: res now passed to plot script as attributes of var
;;    20140312-A_righ_ma: Adapted to new time_operations structure
;;    20140114-A_senf_da: Modularized to fit to new structure
;;    20130613-A_gott_kl: Implemented
;;
;;#############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/calc_seasonal.ncl"
load "./diag_scripts/lib/ncl/ensemble.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/style.ncl"

load "./plot_scripts/ncl/aux_plotting.ncl"
load "./plot_scripts/ncl/contour_maps.ncl"

load "./diag_scripts/lib/ncl/meta_data.ncl"


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
    write_references(diag_script,    \  ; script name
                     "A_gott_kl",    \  ; authors
                     "A_senf_da",    \  ; contributors
                     "D_0001",       \  ; diag_references
                     "",             \  ; obs_references
                     (/"P_embrace", "P_climval"/))   ; proj_acknowledgements

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
;    climo_dir    = getenv("ESMValTool_climo_dir")

    file_type = getenv("ESMValTool_output_file_type")
    if(ismissing(file_type)) then
        file_type = "PS"
    end if

    ;; Check field type
    plot_fields = (/"T2Ms"/)
    ind_f = inlist(field_type0, plot_fields)
    if (all(ismissing(ind_f))) then
        error_msg("f", diag_script, "", "can't use  field " + field_type0)
    end if
    delete(plot_fields)
    delete(ind_f)

    ;; Prepare annotation
    annots = project_style(diag_script_info, "annots")
    smonths = (/"January", "February", "March", "April", "May", "June", \
                "July", "August", "September", "October", "November", \
                "December"/)
    i2months = (/"01", "02", "03", "04", "05", "06", "07", \
                 "08", "09", "10", "11", "12"/)

    ;; Output file name
    if (diag_script_info@month.eq."A") then
        i2month = "annual"
    else
        i2month = smonths(toint(diag_script_info@month) - 1)
    end if
    ;; Specify additional elements for file name
    if (diag_script_info@range_option .eq. 1) then
        name_string = var0 + "_" + diag_script_info@region + "_" + i2month + \
            "_" + sprinti("%0.4i", max(toint(models@start_year))) + "-"\
            + sprinti("%0.4i", min(toint(models@end_year)))
    else if (diag_script_info@range_option .eq. 0) then
        name_string = var0 + "_" + diag_script_info@region + "_" + i2month + \
            "_" + sprinti("%0.4i", min(toint(models@start_year))) + "-"\
            + sprinti("%0.4i", max(toint(models@end_year)))
    else
        error_msg("f", diag_script, "", "range_option not set correctly, " + \
                  "must be 1 or 0")
    end if
    end if

    ;; Must open the workstation before the loop
    wks = get_wks(0, diag_script, name_string)

    ;; Panel parameters
    nvert  = 1  ; Default
    nhori  = 1  ; Default
    if (isatt(diag_script_info, "max_vert")) then
        ;; Maximum allowed number of plots per page (vertical)
        nvert = toint(diag_script_info@max_vert)
    end if
    if (isatt(diag_script_info, "max_hori")) then
        ;; Maximum allowed number of plots per page (horizontal)
        nhori = toint(diag_script_info@max_hori)
    end if
    npages = toint(ceil(tofloat(dimsizes(annots) - 1)/tofloat(nvert * nhori)))

    ;; Check further diag_script_info attributes
    MMM = False
    if (diag_script_info@multi_model_mean.eq."y") then
        MMM    = True
        nvert  = 1
        nhori  = 1
        npages = 1
    end if

    ;; individual label bar (default) or 1 per panel
    lblBarOn = True
    PanelLabelBar = False
    if (isatt(diag_script_info, "PanelLabelBar")) then
        ;; turn on panel label bar?
        if (diag_script_info@PanelLabelBar) then
            lblBarOn = False
            PanelLabelBar = True
        end if
    end if

    ;; Define array for collecting plots
    if (npages .gt. 1) then
        ;; Size of all stamps is identical only if every page contains
        ;; the same number of plots
        plots = new(npages * nvert * nhori, graphic)
    else
        plots = new(dimsizes(annots) - 1, graphic)
    end if

;   start processing ##########################################
    ;; Check if reference model is available in annots
    if (.not. isatt(diag_script_info, "ref_model")) then
        error_msg("f", diag_script, "", "missing " + \
                  "diag_script_info@ref_model -> define in config file")
    end if
    iref = ind(annots .eq. diag_script_info@ref_model)
    ;; ref_model not properly defined -> issue a list of valid strings
    if (ismissing(iref)) then
        error_msg("f", diag_script, "", "diag_script_info@ref_model = " + \
                  diag_script_info@ref_model + " is not valid.")
    else if (dimsizes(iref).ne.1) then  ; Ref_model is ambiguous
        error_msg("f", diag_script, "", "diag_script_info@ref_model = " + \
                  diag_script_info@ref_model + " is ambiguous")
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

    ;; Loop over models
    do iannots = 0, dimsizes(annots) - 1

        ;; Make sure to start with reference model
        imod = ref_first(iannots)

        info_output("processing " + annots(imod) \
                    + "  ############################################", \
                    verbosity, 2)
        ;; See ./interface_scripts/data_handling.ncl
        var_raw = extract_data(imod, data, -1, 0, 0)
        if (var0.eq."sic") then
            if (var_raw@units.eq."%") then
                var_raw = var_raw / 100.    ; Subsequent code expects %/100
                var_raw@units = "%/100"
            end if
        end if

        ;; Show units in title (default False)
        strUnits = ""
        if (isatt(diag_script_info, "showunits")) then
            strUnits = " (" + var_raw@units + ")"
        end if

        ;; Get 2d lat & lon coordinates
        var = make_latlon2D(var0, imod, data, var_raw)
        delete(var_raw)

        ;; Time axis
        start_year = get_start_year(imod)
        end_year   = get_end_year(imod)
        if (diag_script_info@month.eq."A") then
            var_ann = time_operations(var, start_year, end_year, \
                                      "average", "annualclim", True)
        else
            var_ann = time_operations(var, start_year, end_year, \
                                      "average", diag_script_info@month, True)
        end if
        range = start_year + " - " + end_year

        delete(var)

        var = var_ann
        copy_VarMeta(var_ann, var)
        delete(var_ann)

        ;; Regridding
        var_regrid = cdo_remapdis(var0, var, imod, data, \
                                  diag_script_info@dst_grid, "")
        if (var_regrid!0.eq."lon") then  ; reorder dimensions
            var_regrid := var_regrid(lat|:, lon|:)
        end if        

        ;; Some attributes for area plot
        if (diag_script_info@month.eq."A") then  ; Weighted annual mean
            var_regrid@units = "annual mean " + range
        else  ; Select a month
            var_regrid@units = smonths(\
                floattoint(stringtofloat(diag_script_info@month)) - 1) \
                + " mean " + range
        end if
        var_regrid@long_name = diag_script_info@region + \
                               " " + var@long_name
        delete(var)

        ;; Just keep ref_model OR plot difference to ref_model
        if (iannots .eq. 0) then  ; ref_model -> keep, but don't plot

            var_ref = var_regrid
            copy_VarMeta(var_regrid, var_ref)
            delete(var_regrid)

        else  ; plot difference to ref_model

            if (MMM) then
                if (.not.isdefined("all_models")) then
                    all_models = new((/dimsizes(annots), \
                                       dimsizes(var_regrid&lat), \
                                       dimsizes(var_regrid&lon)/), \
                                      typeof(var_regrid))
                    all_models!0 = "models"
                end if
                all_models(iannots, :, :) = var_regrid
                delete(var_regrid)
                if (iannots.eq.dimsizes(annots) - 1) then
                    var_regrid := dim_avg_n_Wrap(all_models, 0)
                else
                    continue
                end if
            end if

            ;; Enforce identical _FillValue for var_ref & var_regrid
            if (var_regrid@_FillValue .ne. var_ref@_FillValue) then
                var_regrid = where(ismissing(var_regrid), \
                                   var_ref@_FillValue, var_regrid)
                var_regrid@_FillValue = var_ref@_FillValue
            end if

            ;; Difference to ref_model
            var = var_regrid - var_ref
            copy_VarMeta(var_regrid, var)
            var = where(ismissing(var_ref), var_ref@_FillValue, var)
            var = where(ismissing(var_regrid), var_ref@_FillValue, var)
            var@long_name = "~F33~D~F~ of " + var@long_name + strUnits
            var_perc = 100. * (var_regrid - var_ref) / \
                where(var_ref .le. 0., var_ref@_FillValue, var_ref)
            copy_VarMeta(var_regrid, var_perc)
            var_perc@long_name = "% ~F33~D~F~ of " + var@long_name
            delete(var_regrid)

            ;; *********** PLOTTING ************
            ;; Set contour parameters
            grid_min    = -1.   ; default
            grid_max    =  1.   ; default
            grid_step   =  0.2  ; default
            grid_center =  0.0  ; default
            if (isatt(diag_script_info, "grid_min")) then
                grid_min = tofloat(diag_script_info@grid_min)
            else
                if (var0.eq."sit") then
                    error_msg("w", diag_script, "", "For sea ice thickness," + \
                              " grid min and max values " + \
                              " should be specified in the cfg file.")
                end if
            end if
            if (isatt(diag_script_info, "grid_max")) then
                grid_max = tofloat(diag_script_info@grid_max)
            end if
            if (isatt(diag_script_info, "grid_step")) then
                grid_step = tofloat(diag_script_info@grid_step)
            end if
            if (isatt(diag_script_info, "grid_center")) then
                grid_center = tofloat(diag_script_info@grid_center)
            end if

            ;; Define color table
            ;; "BlueYellowRed"  ;"BlWhRe"
            colormap = "BlueWhiteOrangeRed"
            if (isatt(diag_script_info, "colormap")) then
                colormap = diag_script_info@colormap
            end if
            gsn_define_colormap(wks, colormap)
            ;; Adjust gsnSpreadColor* (see function gsn_range in
            ;; ./diag_scripts/lib/ncl/style.ncl)

            ;; 255 colors, transition at 98
            gsn_range = gsnColorRange(grid_min, grid_max, \
                                      grid_step, grid_center, 255, 98)

            ;; Define diag-specific resources to be passed to plot routine
            if (MMM) then
                title_model = "MMM"
            else
                title_model = annots(imod)
            end if
            var@res_tiMainString          = title_model + " - " + \
                                            diag_script_info@ref_model
            var@res_gsnRightString        = var@units
            var@res_lbLabelBarOn          = lblBarOn  ; individual bars
            var@res_cnLevelSelectionMode  = "ManualLevels"  ; "ExplicitLevels"
            var@res_cnMinLevelValF        = grid_min
            var@res_cnMaxLevelValF        = grid_max
            var@res_cnLevelSpacingF       = grid_step
            var@res_gsnSpreadColors       = True    ; use full colormap
            var@res_gsnSpreadColorStart   = gsn_range@ColorStart
            var@res_gsnSpreadColorEnd     = gsn_range@ColorEnd
            ;var@res_cnLevels              = ispan(0, 10, 1)/10.
            var@res_cnLinesOn             = False
            var@res_cnLineLabelsOn        = False
            ;var@res_gsnStringFontHeightF  = 15.  ;tune to increase font size
            var@res_gsnDraw               = False  ; Draw and Frame is
            var@res_gsnFrame              = False  ; called in panelling

            ; Pass additional meta data
            var@diag_script = (/diag_script/)

            ; Plot
            plots(iannots - 1) = contour_map_polar(wks, var, variables(0))

            ;; Clean up
            delete(var)

        end if  ; ref_model?

        if (isfilepresent_esmval("./PET0.RegridWeightGen.Log")) then
            system("rm ./PET0.RegridWeightGen.Log")
        end if

    end do  ; iannots

    ;; Draw the panel
    pres = True
    pres@lbLabelStride = 2
    pres@gsnPanelLabelBar = PanelLabelBar  ; 1 bar per panel?

    if (MMM) then
        outfile = panelling(wks, plots(dimsizes(plots)-1), nvert, nhori, pres)
        model_string = "multi-model-mean"
    else
        outfile = panelling(wks, plots, nvert, nhori, pres)
        model_string = "CMIP model"
    end if

    ;; Reporting
    alltags = array_append_record(tags, (/"DM_polar", "PT_polar", \
                                          "ST_mean"/), 0)
    caption = "Polar-stereographic contour map of difference in sea ice concentration between " + model_string + " and " + diag_script_info@ref_model  + ". TODO: Problem with OBS* classes - climofile not passed correctly. for OBS* classes, interface_get_infile already returns path + file! Reason yet unknown."
    if (npages.eq.1 .or. MMM) then  ; only 1 file: no sequence number added
        plot_fname = wks@fullname
        id = "#idSeaIcePolcon"
        climofiles = interface_get_inpaths(0) + "/" + \
                     interface_get_infile(var0, field_type0, 0)
        if (dimsizes(models@name).gt.1) then
            do imod=1, dimsizes(models@name)-1
;;TODO: Problem with OBS* classes: for OBS* classes, interface_get_infile 
;;      already returns path + file! Reason yet unknown.
                climofiles := array_append_record( climofiles, \
                    interface_get_inpaths(imod) + "/" + \
                    interface_get_infile(var0, field_type0, imod), 0 )
            end do
        end if
        contrib_authors = (/"A_gott_kl", "A_senf_da"/)
        ESMValMD(plot_fname, alltags, caption, id, var0, \
                 models@name, climofiles, diag_script, contrib_authors)
    else  ; Multiple files, each with its own sequence number
        ; Split wks@fullname to later add number sequence before '.png'
        plotname_split = str_split(wks@fullname, ".")
        plotname = ""
        do isplit=0, dimsizes(plotname_split)-2  ; reassemble except '.png'
            plotname = plotname + plotname_split(isplit) + "."
        end do
        ; Add number sequence before '.png' and call reporting procedure
        do ipage=1, npages
            seq_num = sprinti("%0.6i", ipage)
            plot_fname = plotname + seq_num + "." + \
                         plotname_split(dimsizes(plotname_split)-1)
            id = "#idSeaIcePolcon" + seq_num
            if (ipage.ne.npages) then
                ind_models_shown = ispan( (ipage-1)*nvert*nhori, \
                                          ipage*nvert*nhori-1, 1 )
            else
                ind_models_shown = ispan( (ipage-1)*nvert*nhori, \
                                          dimsizes(plots)-1, 1 )
            end if
            ind_ref_model = ind(models@name .eq. diag_script_info@ref_model)
            ind_models_shown_exclref = ind_models_shown( get1Dindex_Collapse( \
                                       ind_models_shown, ind_ref_model) )
            climofiles = interface_get_inpaths(ind_ref_model) + "/" + \
                         interface_get_infile(var0, field_type0, ind_ref_model)
            do imod=0, dimsizes(ind_models_shown)-1
;;TODO: Problem with OBS* classes: for OBS* classes, interface_get_infile 
;;      already returns path + file! Reason yet unknown.
                climofiles := array_append_record( climofiles, \
                    interface_get_inpaths(ind_models_shown_exclref(imod)) + \
                    "/" + interface_get_infile(var0, field_type0, \
                    ind_models_shown_exclref(imod)), 0 )
            end do
            contrib_authors = (/"A_gott_kl", "A_senf_da"/)
            ESMValMD(plot_fname, alltags, caption, id, var0, \
                     models@name(ind_models_shown), climofiles, diag_script, \
                     contrib_authors)
            delete([/ seq_num, plot_fname, id, ind_models_shown, \
                      ind_models_shown_exclref, climofiles /])
        end do
    end if

    info_output(" gv " + outfile + " &", verbosity, 1)

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)

end
