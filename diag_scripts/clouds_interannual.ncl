;; CLOUDS
;;#############################################################################
;; Author: Axel Lauer (DLR, Germany)
;; PROJECT-NAME EMBRACE
;;#############################################################################
;; Description
;;    Calculates the interannual variability estimated as the temporal standard
;;    deviation calculated from monthly mean anomalies after subtracting the
;;    climatological mean seasonal cycle.
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    grid:            grid for regridding (coarsest, finest, ref) in case
;;                     calculation og the multi-model mean is requested
;;                     (see optional diag_script_info attributes below)
;;    ref_model:       name of reference data set
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    calcmm:          include calculating of the multi-model mean
;;    colormap:        e.g., WhiteBlueGreenYellowRed, rainbow
;;    ncdf:            enable to output to netCDF; either use "default" or give
;;                     a full file name
;;    projection:      map projection, e.g., Mollweide, Mercator
;;    styleset:        "CMIP5", "DEFAULT"
;;
;; Required variable_info attributes (variable specific)
;;    none
;;
;; Optional variable_info attributes (variable specific)
;;    none
;;
;; Caveats
;;    none
;;
;; Modification history
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

load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
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

    if (flag_mean) then
        ; check for reference model definition
        if (.not.(isatt(diag_script_info, "ref_model"))) then
            error_msg("f", diag_script, "", "no reference model is selected")
        end if

        ; set reference model

        ref_ind = ind(models@name.eq.diag_script_info@ref_model(0))
        ref_model = diag_script_info@ref_model
        if (ismissing(ref_ind)) then
            error_msg("f", diag_script, "", \
                      "no adequate reference model provided")
        end if

        ; set regridding method

        if (.not.(isatt(diag_script_info, "grid"))) then
            error_msg("w", diag_script, "", \
                      "Regridding method not defined. " + "Using coarsest")
            selected_grid = "coarsest"
        else
            selected_grid = diag_script_info@grid
        end if

        ; defined the target grid for regridding

        models_ind = ispan(0, dim_MOD - 1, 1)
        info_output("Setting target grid:", verbosity, 2)
        if (selected_grid.eq."ref") then
            grid = find_destination_grid(ref_ind, \
                                         var0, field_type0, "coarsest")
        else
            grid = find_destination_grid(models_ind, \
                                         var0, field_type0, selected_grid)
        end if
        delete(models_ind)

        ; add entry "multi-model mean"

        dim_MOD = dim_MOD + 1
        modelsatname = array_append_record((/"multi-model mean"/), \
                                           models@name, 0)
        delete(models@name)
        models@name = modelsatname
        delete(modelsatname)

        ; move "multi-model mean" entry to the end of the array

        dumstr = models@name(0)
        do i = 1, dim_MOD - 1
            models@name(i - 1) = models@name(i)
        end do
        mm_ind = dim_MOD - 1
        models@name(mm_ind) = dumstr

        opt = True
        opt@plot_type = "regridded"
        opt@time_avg = selected_grid

        nummod = 0
    else
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

    maps = new(dim_MOD, graphic)  ; collect individual maps in a graphic array

    do imod = 0, dim_MOD - 1

        if (isvar("data1")) then
            delete(data1)
        end if

        if (imod.eq.mm_ind) then
            info_output("processing multi-model mean", verbosity, 0)
            data1 = mmdata / nummod
            copy_VarMeta(mmdata, data1)
        else
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

            ; subtract climatological seasonal cycle from time series

            if (isvar("timeseries")) then
                delete(timeseries)
            end if

            timeseries = calcMonAnomTLL(A0, time_operations(A0, -1, -1, \
                                        "average", "monthlyclim", True))

            ; calculate temporal standard deviation for each grid cell

            data1 = dim_stddev_n_Wrap(timeseries, 0)

            ; normalize standard deviation and convert to percent

            if (isvar("mean")) then
                delete(mean)
            end if
            mean = time_operations(A0, -1, -1, "average", "annualclim", True)
            mean = where(abs(mean).lt.1.0e-4, mean@_FillValue, mean)  ; replace
                                        ; "epsilon" values with missing value
            data1 = 100.0 * data1 / abs(mean)

        end if

        if (verbosity .ge. 5) then
            printVarSummary(data1)
        end if

        ; =====================================================================
        ; =============================== regridding ==========================
        ; =====================================================================

        ; only needed if multi-model mean is requested

        if (flag_mean.and.(imod.ne.mm_ind)) then

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

                zregrid = False  ; plev
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

            ; skip reference model (observations) and multi-model mean entry

            if ((imod.ne.ref_ind).and.(imod.ne.mm_ind)) then
                nummod = nummod + 1
                if (nummod.eq.1) then  ; first element
                    mmdata = var_hreg
                    copy_VarMeta(grid, mmdata)
                else
                    mmdata = mmdata + var_hreg
                end if
            end if

        end if  ; flag_mean

        ; ###########################################
        ; # Style dependent annotation              #
        ; ###########################################
        ; retrieve unique strings describing the data
        ; function in ./diag_scripts/lib/ncl/style.ncl

        annots = project_style(diag_script_info, "annots")

        ; ###########################################
        ; # plot ressources                         #
        ; ###########################################

        data1@res_gsnMaximize       = True         ; Use full page for the plot
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
        data1@res_gsnLeftStringFontHeightF = 0.015
;        data1@res_gsnLeftString           = models@name(imod)
        data1@res_cnLevelSelectionMode     = "ExplicitLevels"

        data1@res_mpOutlineOn              = True
        data1@res_mpFillOn                 = False

        ; variable specific plotting settings

        data1@res_cnLevels                 = ispan(5, 100, 5)

        if (var0.eq."lwp") then
            data1@res_mpOutlineOn          = False
            data1@res_mpFillOn             = True
            data1@res_mpLandFillColor      = "Black"
            pal = read_colormap_file("diag_scripts/lib/ncl/rgb/qcm3.rgb")
            data1@res_cnFillColors         = pal
        end if

        data1@res_lbLabelBarOn             = False
        data1@res_gsnRightString           = ""

        ; map attributes

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
;        data1@res_gsnRightString           = var0 + \
;                                             " in " + variable_info@units
        data1@res_cnInfoLabelOn            = False     ; turn off cn info label

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
        data1@var_units     = "%"

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

        maps(imod) = contour_map(wks, data1, var0)
    end do  ; imod-loop (models)

    pres                                  = True      ; needed to override
                                                      ; panelling defaults
    pres@gsnPanelLabelBar                 = True      ; add common colorbar
    pres@gsnPanelFigureStrings            = annots    ; print model name on
                                                      ; each panel
    pres@gsnPanelFigureStringsFontHeightF = 0.007
    pres@lbLabelFontHeightF               = 0.01
    pres@lbAutoManage                     = False
    pres@lbTopMarginF                     = 0.1
    pres@lbTitleOn                        = True
    pres@lbTitleFontHeightF               = 0.009
    pres@lbTitlePosition                  = "Bottom"
    pres@lbTitleString                    = "~F8~s~F21~" + var0 + " (%)"
    pres@lbPerimOn                        = False     ; draw line around label
                                                      ; bar area
;    pres@gsnPanelTop                      = 0.96
;    pres@gsnPanelBottom                   = 0.012   ; add space for common
;                                                    ; colorbar (default = 0.0)
    pres@gsnPanelCenter                   = False
    pres@pmLabelBarOrthogonalPosF         = -0.01  ; shift label bar a bit to
                                                   ; the bottom
;    pres@gsnPanelDebug                    = True

    outfile = panelling(wks, maps, (dim_MOD + 3) / 4, 4, pres)

    info_output(" Wrote " + outfile, verbosity, 1)

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
