;; CLOUDS
;;#############################################################################
;; Author: Axel Lauer (DLR, Germany)
;; PROJECT-NAME EMBRACE
;;#############################################################################
;; Description
;;    Calculates annual means of 2-d cloud parameters as well as the
;;    multi-model mean for comparison with a reference data set.
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    grid:            grid for regridding (coarsest, finest, ref, 1x1) in
;;                     case calculation of the multi-model mean or difference
;;                     plots are requested
;;                     (see optional diag_script_info attributes below)
;;    ref_model:       name of reference data set
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    calcmm:          include calculating of the multi-model mean
;;    embracesetup:    True = 2 plots per line, False = 4 plots per line (default)
;;    showdiff         calculate and plot differences
;;    timemean:        time averaging - "seasonal" = DJF, MAM, JJA, SON), 
;;                                      "annual" = annual mean
;;
;; Required variable_info attributes (variable specific)
;;    long_name:       description of variable
;;    map_Levels:      contour levels for plotting
;;    units:           variable units
;;
;; Optional variable_info attributes (variable specific)
;;    none
;;
;; Caveats
;;    none
;;
;; Modification history
;;    20160901-A_laue_ax: added regridding option 1 deg x 1 deg
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;;    20150415-A-laue_ax: written.
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
load "./plot_scripts/ncl/contour_maps.ncl"

load "./diag_scripts/aux/perfmetrics/aux_perfmetrics.ncl"


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
                     "D_lauer13jclim", \  ; diag_references
                     "",               \  ; obs_references
                     (/"P_embrace"/))     ; proj_references

    ; 'models', 'variables' and 'field_types' are fetched from the
    ; above 'ncl.interface' file

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    flag_mean = False
    if (isatt(diag_script_info, "calcmm")) then
        if (diag_script_info@calcmm) then
            flag_mean = True
        end if
    end if

    flag_diff = False
    if (isatt(diag_script_info, "showdiff")) then
        if (diag_script_info@showdiff) then
            flag_diff = True
        end if
    end if

    embracesetup = False
    if (isatt(diag_script_info, "embrace_setup")) then
        if (diag_script_info@embrace_setup) then
            embracesetup = True
        end if
    end if

    if (isatt(diag_script_info, "projection")) then
        projection = diag_script_info@projection
        perim      = False
    else
        projection = "CylindricalEquidistant"
        perim      = True
    end if

    ; time averaging: at the moment, only "annualclim" and "seasonalclim"
    ; are supported

    timemean = "annualclim"  ; default
    numseas = 1              ; default
    season = (/"DJF", "MAM", "JJA", "SON"/)

    if (isatt(diag_script_info, "timemean")) then
        timemean = diag_script_info@timemean
        if (timemean.eq."seasonalclim") then
            numseas = 4
        end if
    end if

    panel_labels = True  ; default
    if (isatt(diag_script_info, "panel_labels")) then
        panel_labels = diag_script_info@panel_labels
    end if

    treat_var_as_error = False
    if (isatt(diag_script_info, "treat_var_as_error")) then
        treat_var_as_error = diag_script_info@treat_var_as_error
    end if

    ref_ind = -1  ; set to invalid value

    if ((flag_mean).or.(flag_diff)) then
        ; check for reference model definition
        if (.not.(isatt(diag_script_info, "ref_model"))) then
            error_msg("f", diag_script, "", "no reference model is selected")
        end if

        ; set reference model

        ref_ind = ind(models@name.eq.diag_script_info@ref_model(0))
        ref_model = diag_script_info@ref_model
        if (ismissing(ref_ind)) then
            if (diag_script_info@ref_model(0).eq."n.a.") then
                ref_ind = 0
            else
                error_msg("f", diag_script, "", \
                      "no adequate reference model provided")
            end if
        end if

        ; set regridding method

        if (.not.(isatt(diag_script_info, "grid"))) then
            error_msg("w", diag_script, "", \
                      "Regridding method not defined. " + "Using 1x1")
            selected_grid = "1x1"
        else
            selected_grid = diag_script_info@grid
        end if

        ; define the target grid for regridding

        models_ind = ispan(0, dim_MOD - 1, 1)
        info_output("Setting target grid:", verbosity, 2)
        if (selected_grid.eq."ref") then
            grid = find_destination_grid(ref_ind, \
                                         var0, field_type0, "coarsest")
        else if (selected_grid.eq."1x1") then
            grid     = new((/180, 360/), float)
            grid!0   = "lat"
            grid!1   = "lon"
            grid&lat = flt2dble(fspan(-89.5,  89.5, 180))  ; CMOR standard
            grid&lon = flt2dble(fspan(0.5, 359.5, 360))    ; CMOR standard
            grid&lat@standard_name = "latitude"
            grid&lon@standard_name = "longitude"
            grid&lat@long_name     = "latitude"
            grid&lon@long_name     = "longitude"
            grid&lat@units         = "degrees_north"
            grid&lon@units         = "degrees_east"
            grid&lat@axis          = "Y"
            grid&lon@axis          = "X"
        else
            grid = find_destination_grid(models_ind, \
                                         var0, field_type0, selected_grid)
        end if
        end if
        delete(models_ind)

        if (flag_mean) then
            ; add entry "multi-model mean"

            dim_MOD = dim_MOD + 1
            modelsatname = array_append_record((/"multi-model mean"/), \
                                               models@name, 0)
            delete(models@name)
            models@name = modelsatname
            delete(modelsatname)

            projectsatname = array_append_record((/"mmm"/), models@project, 0)
            delete(models@project)
            models@project = projectsatname
            delete(projectsatname)

            ; move "multi-model mean" entry to the end of the array

            dumstr = models@name(0)
            do i = 1, dim_MOD - 1
                models@name(i - 1) = models@name(i)
            end do
            mm_ind = dim_MOD - 1
            models@name(mm_ind) = dumstr
        end if

        opt = True
        opt@plot_type = "regridded"
        opt@time_avg = selected_grid

        nummod = 0
    end if

    if (.not.flag_mean) then
        mm_ind = -1
    end if

end

; load variable_info

begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin
    ; ###########################################
    ; # get data and average time               #
    ; ###########################################

    maps = new((/dim_MOD, 4/), graphic)  ; collect individual maps in a graphic array

    do imod = 0, dim_MOD - 1

        if (isvar("data1")) then
            delete(data1)
        end if

        if (flag_mean.and.(imod.eq.mm_ind)) then
            info_output("processing multi-model mean", verbosity, 0)
            data1 = mmdata / nummod
            copy_VarMeta(mmdata, data1)
        end if

        if (imod.ne.mm_ind) then
            info_output("processing " + models@name(imod), verbosity, 0)

            if (isvar("A0")) then
                delete(A0)
            end if

            data = read_data(imod, var0, field_type0)
            A0 = extract_data(imod, data, -1, 0, 0)
            info_output(" ++++++++++++++ Variable summary of original data " \
                        + "(A0): ", verbosity, 5)
            if (verbosity .ge. 4) then
                printVarSummary(A0)
            end if

            ; check dimensions

            dims = getvardims(A0)
            if (dimsizes(dims) .lt. 2) then
                error_msg("f", diag_script, "", dimsizes(dims) + \
                          " dimensions, need 2 or 3")
            end if
            idx = ind(dims .eq. "lat")
            if (ismissing(idx)) then
                error_msg("f", diag_script, "", "no lat dimension")
            end if
            idx = ind(dims .eq. "lon")
            if (ismissing(idx)) then
                error_msg("f", diag_script, "", "no lon dimension")
            end if

            ; average over time
            ; (function in ./diag_scripts/lib/ncl/statistics.ncl)

            ; if variable is an error variable, we have to square it before
            ; averaging and then calculate the square-root afterwards

            if (treat_var_as_error) then
                info_output(" ++++++++++++++ Treating variable as error " + \
                            "variable when averaging ", verbosity, 0)
                A0 = A0 * A0
            end if

            data1 = time_operations(A0, -1, -1, "average", timemean, True)

            if (treat_var_as_error) then
                data1 = sqrt(data1)
            end if

            info_output(" ++++++++++++++ Variable summary of time averaged" + \
                        " data (data1): ", verbosity, 5)

        end if

        if (verbosity .ge. 5) then
            printVarSummary(data1)
        end if

        ; =====================================================================
        ; =============================== regridding ==========================
        ; =====================================================================

        ; always regrid (for calculating correlation coefficients)

;        ; only needed if multi-model mean and/or differences are requested

        if ((flag_mean.or.flag_diff).and.(imod.ne.mm_ind)) then
;        if (imod.ne.mm_ind) then

            info_output("Regridding " + models@name(imod), verbosity, 2)
            regr_path = clim_diag_file_path(var0, field_type0, imod, \
                                            diag_script, opt)

            ; Do not process data if already processed

            if (isfilepresent(regr_path) .and. \
                getenv("ESMValTool_force_calc").eq."False") then

                info_output("Reading in file = " + regr_path, verbosity, 1)
                var_hreg = ncdf_read(regr_path, var0)
            else

                ; Determine dimensions to regrid

                xregrid = False  ; lon
                yregrid = False  ; lat
                dnames = getVarDimNames(data1)
                if (any(dnames.eq."lon")) then
                    xregrid = True
                    if (dimsizes(data1&lon).eq.dimsizes(grid&lon)) then
                        if (max(abs(data1&lon - grid&lon)).lt.1.0d-5) then
                            xregrid = False
                        end if
                    end if
                end if
                if (any(dnames.eq."lat")) then
                    yregrid = True
                    if (dimsizes(data1&lat).eq.dimsizes(grid&lat)) then
                        if (max(abs(data1&lat - grid&lat)).lt.1.0d-5) then
                            yregrid = False
                        end if
                    end if
                end if

                ; regrid in the horizontal (if necessary)

                if (xregrid.and.yregrid) then
                    var_hreg = linint2_Wrap(data1&lon, data1&lat, data1, \
                                            True, grid&lon, grid&lat, 0)
                else if (xregrid.and..not.yregrid) then
                    dx = ind(dnames.eq."lon")
                    var_hreg = linint1_n_Wrap(data1&lon, data1, True, \
                                              grid&lon, 0, dx)
                    info_output("Latitude regridding not required", \
                                verbosity, 4)
                else if (.not.xregrid.and.yregrid) then
                    dx = ind(dnames.eq."lat")
                    var_hreg = linint1_n_Wrap(data1&lat, data1, False, \
                                              grid&lat, 0, dx)
                    info_output("Longitude regridding not required", \
                                verbosity, 4)
                else if (.not.xregrid .and. .not.yregrid) then
                    var_hreg = data1
                    info_output("Horizontal regridding not required", \
                                verbosity, 4)
                end if
                end if
                end if
                end if

                ; output to NetCDF

                var_hreg@var = var0
                var_hreg@diag_script = (/diag_script/)
                var_hreg@ncdf = regr_path
                ncdf_outfile = ncdf_write(var_hreg, regr_path)
            end if

            if (flag_mean) then
                ; skip reference model (observations) and multi-model mean entry

                if ((imod.ne.ref_ind).and.(imod.ne.mm_ind)) then
                    nummod = nummod + 1
                    if (nummod.eq.1) then  ; first element
                        mmdata = var_hreg
                        if (numseas.eq.1) then
                            copy_VarMeta(grid, mmdata)
                        else
                            do is = 0, numseas - 1
                                copy_VarMeta(grid, mmdata(is, :, :))
                            end do
                        end if
                    else
                        mmdata = mmdata + var_hreg
                    end if
                end if
            end if  ; flag_mean

            if (flag_diff) then
                if (imod.eq.0) then
                    tmp = dimsizes(grid)
                    if (dimsizes(dimsizes(data1)).eq.2) then
                        diff = new((/dim_MOD, tmp(0), tmp(1)/), float)
                        do j = 0, dim_MOD - 1
                            copy_VarMeta(grid, diff(j, :, :))
                        end do
                    else
                        diff = new((/dim_MOD, numseas, tmp(0), tmp(1)/), float)
                        do j = 0, dim_MOD - 1
                            do is = 0, numseas - 1
                                copy_VarMeta(grid, diff(j, is, :, :))
                            end do
                        end do
                    end if
                    delete(tmp)
                    diff@res = True
                end if
                if (dimsizes(dimsizes(var_hreg)).eq.2) then
                    diff(imod, :, :) = var_hreg
                else
                    diff(imod, :, :, :) = var_hreg
                end if
            end if  ; flag_diff

        end if  ; flag_mean or flag_diff

     end do  ; imod-loop

     do imod = 0, dim_MOD - 1

        if (isvar("data1")) then
            delete(data1)
        end if

        if (flag_mean.and.(imod.eq.mm_ind)) then
            info_output("processing multi-model mean", verbosity, 0)
            data1 = mmdata / nummod
            copy_VarMeta(mmdata, data1)
        end if

        if (imod.ne.mm_ind) then
            if (isvar("A0")) then
                delete(A0)
            end if

            data = read_data(imod, var0, field_type0)
            A0 = extract_data(imod, data, -1, 0, 0)

            ; average over time

            data1 = time_operations(A0, -1, -1, "average", timemean, True)
        end if

        ; ###########################################
        ; # Style dependent annotation              #
        ; ###########################################
        ; retrieve unique strings describing the data
        ; function in ./diag_scripts/lib/ncl/style.ncl

        annots = project_style(diag_script_info, "annots")

        ; ###########################################
        ; # plot ressources                         #
        ; ###########################################

        data1@res_cnFillOn          = True         ; Color plot desired
        data1@res_cnLineLabelsOn    = False        ; Contour lines
        data1@res_cnLinesOn         = False

        ; colors
        ; http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml

;        data1@res_cnFillPalette           = diag_script_info@colormap

        ; annotation

;        data1@res_tiMainString            = MyParam  ; Main title
        data1@res_tiMainOn                 = False
;        data1@res_gsnLeftString           = annots(imod)
;        data1@res_gsnLeftStringFontHeightF = 0.015
;        data1@res_gsnLeftString           = models@name(imod)
        data1@res_cnLevelSelectionMode     = "ExplicitLevels"

        data1@res_mpOutlineOn              = True
        data1@res_mpFillOn                 = False

        ; variable specific plotting settings

        if (isatt(variable_info, "map_ref_Levels")) then
            if (isatt(data1, "res_cnLevels")) then
                delete(data1@res_cnLevels)
            end if
            data1@res_cnLevels = variable_info@map_ref_Levels
        end if

        if (var0.eq."pr-mmday") then
            if (isatt(data1, "res_cnLevels")) then
                delete(data1@res_cnLevels)
            end if
            data1@res_cnLevels = fspan(0.5, 10, 20)
        end if

        if (var0.eq."prc-mmh") then
            if (isatt(data1, "res_cnLevels")) then
                delete(data1@res_cnLevels)
            end if
            data1@res_cnLevels = fspan(0.025, 0.5, 20)
            pal = read_colormap_file("diag_scripts/lib/ncl/rgb/ipcc-tas.rgb")
;            pal = read_colormap_file("diag_scripts/lib/ncl/rgb/qcm3.rgb")
            data1@res_cnFillColors             = pal
        end if

        if (var0.eq."lwp") then
            delete(data1@res_cnLevels)
            data1@res_cnLevels                 = variable_info@map_ref_Levels * 0.001
            data1@res_mpOutlineOn              = False
            data1@res_mpFillOn                 = True
            data1@res_mpLandFillColor          = "Black"
            pal = read_colormap_file("diag_scripts/lib/ncl/rgb/qcm3.rgb")
            data1@res_cnFillColors             = pal
        end if

        if (var0.eq."tas-degC") then
            pal = read_colormap_file("diag_scripts/lib/ncl/rgb/ipcc-tas.rgb")
            data1@res_cnFillColors             = pal
        end if

        if (var0.eq."clt") then
            if (isatt(data1, "res_cnLevels")) then
                delete(data1@res_cnLevels)
            end if
            data1@res_cnLevels = fspan(5, 100, 20)
        end if

        if (var0.eq."clivi") then
            delete(data1@res_cnLevels)
            data1@res_cnLevels                 = variable_info@map_ref_Levels * 0.001
        end if

        data1@res_lbLabelBarOn             = False
        data1@res_gsnRightString           = ""

        data1@res_mpFillDrawOrder          = "PostDraw"     ; draw map last
        data1@res_cnMissingValFillColor    = "Gray"

        ; no tickmarks and no labels

        data1@res_tmYLLabelsOn             = False
        data1@res_tmYLOn                   = False
        data1@res_tmYRLabelsOn             = False
        data1@res_tmYROn                   = False
        data1@res_tmXBLabelsOn             = False
        data1@res_tmXBOn                   = False
        data1@res_tmXTLabelsOn             = False
        data1@res_tmXTOn                   = False
;        data1@res_gsnRightString           = variable_info@long_name + \
;                                             " in " + variable_info@units
        data1@res_cnInfoLabelOn            = False     ; turn off cn info label
        data1@res_mpPerimOn                = perim     ; draw line around map

        data1@res_mpGridLineColor   = -1
        data1@res_mpGridAndLimbOn   = True

        ; specified in ./nml/cfg_clouds/cfg_clouds.ncl

        data1@res_mpProjection             = diag_script_info@projection

        ; set explicit contour levels

        if (isatt(diag_script_info, "explicit_cn_levels")) then
            data1@res_cnLevelSelectionMode = "ExplicitLevels"
            data1@res_cnLevels = diag_script_info@explicit_cn_levels
        end if

        ; ###########################################
        ; # other Metadata: diag_script, var        #
        ; ###########################################
        ; add to data1, as attributes without prefix

        if (isatt(data1, "diag_script")) then  ; add to existing entries
            temp = data1@diag_script
            delete(data1@diag_script)
            data1@diag_script = array_append_record(temp, (/diag_script/), 0)
            delete(temp)
        else  ; add as new attribute
            data1@diag_script = (/diag_script/)
        end if
        data1@var = var0  ; Overwrite existing entry
        if (isatt(variable_info, "long_name")) then
            data1@var_long_name = variable_info@long_name
        end if
        if (isatt(variable_info, "units")) then
            data1@var_units     = variable_info@units
        end if

        corr = new((/numseas/), float)
        gavg = new((/numseas/), float)

        if (all(ismissing(data1))) then
            data1@res_gsnRightString = ""
            data1@res_gsnLeftString = ""
        else
            if (numseas.gt.1) then
                do is = 0, numseas - 1
                    if (flag_diff) then
                        corr(is) = calculate_metric(diff(ref_ind, is, :, :), diff(imod, is, :, :), "correlation")
                    else
                        corr(is) = -999.
                    end if
                    gavg(is) = area_operations(data1(is, :, :), -90., 90., 0., 360., "average", True)
                end do
            else
                if (flag_diff) then
                    corr(0) = calculate_metric(diff(ref_ind, :, :), diff(imod, :, :), "correlation")
                else
                    corr(0) = -999.
                end if
                gavg(0) = area_operations(data1, -90., 90., 0., 360., "average", True)
            end if
        end if

        data1@res_gsnLeftStringFontHeightF  = min((/0.02, 0.015 * 6.0 / tofloat((dim_MOD + 1) / 2)/))
        data1@res_gsnRightStringFontHeightF = min((/0.02, 0.015 * 6.0 / tofloat((dim_MOD + 1) / 2)/))

        ; ###########################################
        ; # create the plot                         #
        ; ###########################################

        data1@res_gsnDraw           = False  ; Do not draw yet
        data1@res_gsnFrame          = False  ; Don't advance frame.

        ; function in aux_plotting.ncl

        if (imod.eq.0) then
            wks = get_wks("dummy_for_wks", diag_script, var0)
;            drawNDCGrid(wks)  ; debugging option
        end if

        if (numseas.gt.1) then
            do is = 0, numseas - 1
                if (corr(is).ge.-1.) then
                    data1@res_gsnRightString = "corr = " + sprintf("%6.3f", corr(is))
                    data1@res_gsnLeftString = "mean = " + sprintf("%6.3f", gavg(is))
                end if
                maps(imod, is) = contour_map(wks, data1(is, :, :), var0)
            end do
        else
            if (corr(0).ge.-1.) then
                data1@res_gsnRightString = "corr = " + sprintf("%6.3f", corr(0))
                data1@res_gsnLeftString = "mean = " + sprintf("%6.3f", gavg(0))
            end if
            maps(imod, 0) = contour_map(wks, data1, var0)
        end if

        if (imod.eq.ref_ind) then  ; remove correlation string for reference model
            data1@res_gsnRightString = ""
        end if
    end do  ; imod-loop (models)

    ; sort plots if needed (observations go first)

    plottmp = ispan(0, dim_MOD - 1, 1)
    plotind = plottmp

    ; if there are observational datasets, move the corresponding plots to the first lines
    ; of the panel plot

    j = 0
    do i = 0, dimsizes(plottmp) - 1
        if (i.eq.ref_ind) then
            plotind(j) = plottmp(i)
            j = j + 1
        else if (plottmp(i) .lt. dimsizes(models@project)) then
            if (isStrSubset(str_lower(models@project(plottmp(i))), "obs")) then
               plotind(j) = plottmp(i)
               j = j + 1
            end if
        end if
        end if
    end do
;    do i = 0, dimsizes(plottmp) - 1
;        if ((isStrSubset(str_lower(models@project(plottmp(i))), "obs")).or.(i.eq.ref_ind)) then
;        else
;            plotind(j) = plottmp(i)
;            j = j + 1
;        end if
;    end do

    pres                                  = True      ; needed to override
                                                      ; panelling defaults
    pres@gsnPanelLabelBar                 = True      ; add common colorbar
    if (panel_labels) then
        pres@gsnPanelFigureStrings        = annots(plotind)    ; print model name on
                                                               ; each panel
    end if
    pres@gsnPanelFigureStringsFontHeightF = min((/0.016, 0.008 * 6.0 / tofloat((dim_MOD + 1) / 2)/))
    pres@lbLabelFontHeightF               = min((/0.02, 0.01 * 6.0 / tofloat((dim_MOD + 1) / 2)/))
    pres@lbAutoManage                     = False
    pres@lbTopMarginF                     = 0.1
    pres@lbTitleOn                        = True
    pres@lbTitleFontHeightF               = min((/0.02, 0.01 * 6.0 / tofloat((dim_MOD + 1) / 2)/))
    pres@lbTitlePosition                  = "Bottom"
    pres@lbTitleString                    = variable_info@long_name + " (" + \
                                            variable_info@units + ")"
    pres@lbPerimOn                        = False     ; draw line around label
                                                      ; bar area
;    pres@gsnPanelTop                      = 0.96
;    pres@gsnPanelBottom                   = 0.012   ; add space for common
;                                                    ; colorbar (default = 0.0)
    pres@gsnPanelCenter                   = False
    pres@pmLabelBarOrthogonalPosF         = -0.01  ; shift label bar a bit to
                                                   ; the bottom
;    pres@gsnPanelDebug                    = True

    if (embracesetup) then
    if (numseas.gt.1) then
        do is = 0, numseas - 1
            pres@txString   = season(is)
            outfile = panelling(wks, maps(plotind, is), (dim_MOD + 3) / 4, 4, pres)
            info_output(" Wrote " + outfile, verbosity, 1)
        end do
    else
;        outfile = panelling(wks, maps(:, 0), (dim_MOD + 3) / 4, 4, pres)
;        outfile = panelling(wks, maps(:, 0), (dim_MOD + 1) / 2, 2, pres)

        pres@gsnPanelRowSpec = True                   ; tell panel what order to plt
;        pres@gsnMaximize       = True
        pres@gsnPanelYWhiteSpacePercent = 5
        pres@gsnPanelXWhiteSpacePercent = 5
        if (isatt(diag_script_info, "PanelTop")) then
            top = tofloat(diag_script_info@PanelTop)
        else
            top = 0.99  ; default
        end if
        pres@gsnPanelTop = top

        if (isvar("plotsperline")) then
            delete(plotsperline)
        end if

        plotsperline = new((dim_MOD + 1) / 2, integer)
        plotsperline = 2

        if ((isStrSubset(str_lower(models@project(plotind(0))), "obs")).and. \
            .not.(isStrSubset(str_lower(models@project(plotind(1))), "obs"))) then
            plotsperline(0) = 1
        end if

        if (sum(plotsperline).gt.dimsizes(plotind)) then
            plotsperline(dimsizes(plotsperline) - 1) = 1
        end if

        if (sum(plotsperline).lt.dimsizes(plotind)) then
            xadd = 1
            xtmp = array_append_record(plotsperline, xadd, 0)
            delete(plotsperline)
            plotsperline = xtmp
            delete(xtmp)
        end if

        gsn_panel(wks, maps(plotind, 0), plotsperline, pres)

;        info_output(" Wrote " + outfile, verbosity, 1)
    end if
    else  ; if embracesetup
    if (numseas.gt.1) then
        do is = 0, numseas - 1
            pres@txString   = season(is)
            outfile = panelling(wks, maps(:, is), (dim_MOD + 3) / 4, 4, pres)
            info_output(" Wrote " + outfile, verbosity, 1)
        end do
    else
        outfile = panelling(wks, maps(:, 0), (dim_MOD + 3) / 4, 4, pres)
        info_output(" Wrote " + outfile, verbosity, 1)
    end if
    end if

    ; ===================================================================================

    if (.not.flag_diff) then
        exit
    end if

    pres@lbTitleString                    = "~F33~D~F21~" + variable_info@long_name + " (" + \
                                            variable_info@units + ")"

    do imod = 0, dim_MOD - 1
        if (imod.eq.ref_ind) then
            continue
        end if

        if (numseas.gt.1) then
;            do is = 0, numseas - 1
;                diff@res_gsnLeftString = "rmsd = " + calculate_metric(diff(imod, is, :, :), diff(ref_ind, is, :, :), "RMSD")
;            end do
            diff(imod, :, :, :) = diff(imod, :, :, :) - diff(ref_ind, :, :, :)
            diff@res_gsnLeftString  = ""
            diff@res_gsnRightString = ""
        else
            rmsd = calculate_metric(diff(imod, :, :), diff(ref_ind, :, :), "RMSD")
            diff@res_gsnRightString = "rmsd = " + sprintf("%6.3f", rmsd)
            diff(imod, :, :) = diff(imod, :, :) - diff(ref_ind, :, :)
            if (all(ismissing(diff(imod, :, :)))) then
                diff@res_gsnLeftString = ""
                diff@res_gsnRightString = ""
            else
                gavg = area_operations(diff(imod, :, :), -90., 90., 0., 360., "average", True)
                diff@res_gsnLeftString = "bias = " + sprintf("%6.3f", gavg)
            end if
        end if

        ; -------------------------------------------------------------------------------

        ; ###########################################
        ; # plot ressources                         #
        ; ###########################################

        diff@res_gsnLeftStringFontHeightF  = min((/0.02, 0.015 * 6.0 / tofloat((dim_MOD + 1) / 2)/))
        diff@res_gsnRightStringFontHeightF = min((/0.02, 0.015 * 6.0 / tofloat((dim_MOD + 1) / 2)/))

        diff@res_tiMainOn      = False

        diff@res_cnFillOn          = True         ; Color plot desired
        diff@res_cnLineLabelsOn    = False        ; Contour lines
        diff@res_cnLinesOn         = False

        ; colors
        ; http://www.ncl.ucar.edu/Document/Graphics/color_table_gallery.shtml

;        diff@res_cnFillPalette           = diag_script_info@colormap

        ; annotation

;        diff@res_tiMainString            = MyParam  ; Main title
;        diff@res_tiMainOn                 = False
;        diff@res_gsnLeftString           = annots(imod)
;        diff@res_gsnLeftStringFontHeightF = 0.015
;        diff@res_gsnLeftString           = models@name(imod)
        diff@res_cnLevelSelectionMode     = "ExplicitLevels"

        diff@res_mpOutlineOn              = True
        diff@res_mpFillOn                 = False

        ; variable specific plotting settings

        if (isatt(variable_info, "map_diff_Levels")) then
            if (isatt(diff, "res_cnLevels")) then
                delete(diff@res_cnLevels)
            end if
            diff@res_cnLevels = variable_info@map_diff_Levels
        end if

        if ((var0.eq."pr-mmday").or.(var0.eq."prc-mmh")) then
            if (imod.eq.ref_ind) then
                if (isatt(diff, "res_cnLevels")) then
                    delete(diff@res_cnLevels)
                end if
                if (isatt(diff, "res_cnFillColors")) then
                    delete(diff@res_cnFillColors)
                end if
                if (isvar("pal")) then
                    delete(pal)
                end if
                diff@res_cnLevels = fspan(0.5, 7, 14)
;                pal = read_colormap_file("diag_scripts/lib/ncl/rgb/rainbow.rgb")
;                diff@res_cnFillPalette = pal
                diff@res_cnFillPalette = "wgne15"
                diff@res_lbOrientation = "vertical"
            else
                if (isatt(diff, "res_cnLevels")) then
                    delete(diff@res_cnLevels)
                end if
                if (isatt(diff, "res_cnFillColors")) then
                    delete(diff@res_cnFillColors)
                end if
                if (isvar("pal")) then
                    delete(pal)
                end if
;                diff@res_cnLevels = fspan(-3, 3, 13)
                diff@res_cnLevels = ispan(-30, 30, 5) * 0.1
                pal = read_colormap_file("diag_scripts/lib/ncl/rgb/ipcc-precip-delta.rgb")
                diff@res_cnFillColors = pal
                diff@res_lbOrientation = "horizontal"
            end if
            if (var0.eq."prc-mmh") then
                diff@res_cnLevels = diff@res_cnLevels * 0.05
            end if
        end if

        if (var0.eq."tas-degC") then
            if (isatt(diff, "res_cnFillColors")) then
                delete(diff@res_cnFillColors)
            end if
            if (isvar("pal")) then
                delete(pal)
            end if
            pal = read_colormap_file("diag_scripts/lib/ncl/rgb/ipcc-tas-delta.rgb")
            diff@res_cnFillPalette = pal
        end if

        if (var0.eq."lwp") then
            delete(diff@res_cnLevels)
            diff@res_cnLevels                 = variable_info@map_diff_Levels * 0.001
            diff@res_mpOutlineOn              = False
            diff@res_mpFillOn                 = True
            diff@res_mpLandFillColor          = "Black"
            pal = read_colormap_file("diag_scripts/lib/ncl/rgb/qcm3.rgb")
            diff@res_cnFillColors             = pal
        end if

        if (var0.eq."clt") then
            if (isatt(diff, "res_cnLevels")) then
                delete(diff@res_cnLevels)
            end if
            diff@res_cnLevels = fspan(-25, 25, 11)
        end if

        if (var0.eq."clivi") then
            delete(diff@res_cnLevels)
            diff@res_cnLevels                 = variable_info@map_diff_Levels * 0.001
        end if

        if (imod.eq.ref_ind) then
            diff@res_lbLabelBarOn = True
        else
            diff@res_lbLabelBarOn = False
        end if

        ; map attributes

        diff@res_mpFillDrawOrder          = "PostDraw"     ; draw map last
        diff@res_cnMissingValFillColor    = "Gray"

        ; no tickmarks and no labels

        diff@res_tmYLLabelsOn             = False
        diff@res_tmYLOn                   = False
        diff@res_tmYRLabelsOn             = False
        diff@res_tmYROn                   = False
        diff@res_tmXBLabelsOn             = False
        diff@res_tmXBOn                   = False
        diff@res_tmXTLabelsOn             = False
        diff@res_tmXTOn                   = False
;        diff@res_gsnRightString           = variable_info@long_name + \
;                                             " in " + variable_info@units
        diff@res_cnInfoLabelOn            = False     ; turn off cn info label

        ; specified in ./nml/cfg_clouds/cfg_clouds.ncl

        diff@res_mpProjection             = diag_script_info@projection

        ; set explicit contour levels

        if (isatt(diag_script_info, "explicit_cn_levels")) then
            diff@res_cnLevelSelectionMode = "ExplicitLevels"
            diff@res_cnLevels = diag_script_info@explicit_cn_levels
        end if

        ; ###########################################
        ; # other Metadata: diag_script, var        #
        ; ###########################################
        ; add to data1, as attributes without prefix

        diff@var = var0  ; Overwrite existing entry
        if (isatt(variable_info, "long_name")) then
            diff@var_long_name = variable_info@long_name
        end if
        if (isatt(variable_info, "units")) then
            diff@var_units     = variable_info@units
        end if

        ; ###########################################
        ; # create the plot                         #
        ; ###########################################

        diff@res_gsnDraw           = False  ; Do not draw yet
        diff@res_gsnFrame          = False  ; Don't advance frame.

        ; -------------------------------------------------------------------------------

        if (numseas.gt.1) then
            do is = 0, numseas - 1
                maps(imod, is) = contour_map(wks, diff(imod, is, :, :), var0)
            end do
        else
            maps(imod, 0) = contour_map(wks, diff(imod, :, :), var0)
        end if

    end do  ; imod-loop (models)

    if (isvar("plottmp")) then
        delete(plottmp)
    end if

    if (isvar("plotind")) then
        delete(plotind)
    end if
 
    plottmp = ind(ispan(0, dim_MOD - 1, 1).ne.ref_ind)
    plotind = plottmp

    ; if there is a second observational dataset, move the corresponding plot to the first line
    ; of the panel plot

    j = 0
    do i = 0, dimsizes(plottmp) - 1
        if (isStrSubset(str_lower(models@project(plottmp(i))), "obs")) then
            plotind(j) = plottmp(i)
            j = j + 1
        end if
    end do
    do i = 0, dimsizes(plottmp) - 1
        if (isStrSubset(str_lower(models@project(plottmp(i))), "obs")) then
        else
            plotind(j) = plottmp(i)
            j = j + 1
        end if
    end do

    if (isatt(pres, "gsnPanelFigureStrings")) then
        delete(pres@gsnPanelFigureStrings)
    end if
    if (panel_labels) then
        pres@gsnPanelFigureStrings            = annots(plotind)
    end if

    if (embracesetup) then
    if (numseas.gt.1) then
        do is = 0, numseas - 1
            pres@txString   = season(is)
            outfile = panelling(wks, maps(plotind, is), (dim_MOD + 3) / 4, 4, pres)
            info_output(" Wrote " + outfile, verbosity, 1)
        end do
    else
;        outfile = panelling(wks, maps(plotind, 0), (dim_MOD + 3) / 4, 4, pres)
        pres@gsnPanelRowSpec = True                   ; tell panel what order to plt
;        pres@gsnMaximize       = True
        pres@gsnPanelYWhiteSpacePercent = 5
        pres@gsnPanelXWhiteSpacePercent = 5
        if (isatt(diag_script_info, "PanelTop")) then
            top = tofloat(diag_script_info@PanelTop)
        else
            top = 0.99  ; default
        end if
        pres@gsnPanelTop = top

        if (isvar("plotsperline")) then
            delete(plotsperline)
        end if

        plotsperline = new(dim_MOD / 2, integer)
        plotsperline = 2

        if ((isStrSubset(str_lower(models@project(plotind(0))), "obs")).and. \
            .not.(isStrSubset(str_lower(models@project(plotind(1))), "obs"))) then
            plotsperline(0) = 1
        end if

        if (sum(plotsperline).gt.dimsizes(plotind)) then
            plotsperline(dimsizes(plotsperline) - 1) = 1
        end if

        if (sum(plotsperline).lt.dimsizes(plotind)) then
            xadd = 1
            xtmp = array_append_record(plotsperline, xadd, 0)
            delete(plotsperline)
            plotsperline = xtmp
            delete(xtmp)
        end if

        gsn_panel(wks, maps(plotind, 0), plotsperline, pres)

;        outfile = panelling(wks, maps(plotind, 0), (dim_MOD + 1) / 2, 2, pres)
;        info_output(" Wrote " + outfile, verbosity, 1)
    end if
    else  ; embracesetup = False
    if (numseas.gt.1) then
        do is = 0, numseas - 1
            pres@txString   = season(is)
            outfile = panelling(wks, maps(:, is), (dim_MOD + 3) / 4, 4, pres)
            info_output(" Wrote " + outfile, verbosity, 1)
        end do
    else
        outfile = panelling(wks, maps(:, 0), (dim_MOD + 3) / 4, 4, pres)
        info_output(" Wrote " + outfile, verbosity, 1)
    end if
    end if  ; end if embracesetup

    ; ===================================================================================

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
