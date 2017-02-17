;;#############################################################################
;; DIAGNOSTIC SCRIPT FOR AEROSOL STATION DATA
;; Author: Mattia Righi (DLR, Germany)
;; ESMVal project
;##############################################################################
;;
;; Description
;;    Plots surface-level concentrations and aerosol optical depth and comapres
;;    with station data.
;;    Path to observational data is provided in the cfg namelists.
;;    Data are plotted as time series (model/obs vs. time), scatter plot (model
;;    vs. obs) and contour plot (model data with observations overlaid as
;;    markers).
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    network: name of the station network for observational data
;;             (AERONET, CASTNET, EANET, EMEP or IMPROVE).
;;    time_avg: time average (monthly or yearly).
;;    datapath: path to the observational data (observations must be of type
;;              type T0M and must contain the global attributes stname, stlat,
;;              stlon and stalt.
;;    scatter_log: if True, use log axes in the scatter plot.
;;    legend_outside: save legend in a separate file.
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    styleset: project-specific styleset.
;;
;; Required variable_info attributes (variable specific)
;;    contours_abs: array of contour levels for the absolute value map
;;                         plot
;;
;; Caveats
;;    Input units for aerosol concentrations must be kg/m3 (as per CMOR
;;    standard). Conversion to microg/m3 is performed within the script.
;;
;; Modification history
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;;    20150325-A_laue_ax: modified reference tags used for acknowledgements
;;                        (projects, observations, etc)
;;    20140528-A_righ_ma: added extra check for data availability.
;;    20140324-A_righ_ma: improved modularization.
;;    20140114-A_righ_ma: adapted for using time_operations.
;;    20131212-A_righ_ma: adapted to the new tool structure.
;;    20130619-A_righ_ma: adapted to standard namelist.
;;    20130529-A_righ_ma: written.
;;
;;#############################################################################

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/scaling.ncl"
load "./diag_scripts/aux/aerosol/aux_aerosol.ncl"

load "./plot_scripts/ncl/xy_line.ncl"
load "./plot_scripts/ncl/scatterplot.ncl"
load "./plot_scripts/ncl/contour_maps.ncl"

begin

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, \
                verbosity, 2)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    ;; Reference information
    write_references(diag_script,\
                     "A_righ_ma", \  ; authors
                     "A_kais_ch", \  ; contributors
                     "D_0000",    \  ; diag_ref
                     (/"E_improve", "E_castnet", "E_emep", "E_eanet", "E_aeronet"/), \  ; obs_ref
                     (/"P_esmval", "P_dlrveu"/) \  ; proj_acknow
                     )

    ;; Check for valid fields
    field_type0 = field_types(0)
    if (field_type0 .ne. "T2Ms") then
        error_msg("f", diag_script, "", "type " + field_type0 + \
                  " not supported")
    end if

    ;; Check required diag_script_info attributes
    req_atts = (/"network", "time_avg", "datapath", "scatter_log", \
                 "legend_outside"/)
    exit_if_missing_atts(diag_script_info, req_atts)

    ;; Read data
    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    do idx_mod = 0, dim_MOD - 1
        data = read_data(idx_mod, var0, field_type0)
    end do

end

;; Load var_info
begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

;; Main
begin

    ;; Check for data availability
    if (all(ismissing(systemfunc("ls " + diag_script_info@datapath + \
                                 "/OBS_" + diag_script_info@network + \
                                 "_ground_*_T0M_" + var0 + ".nc")))) then
        info_output(diag_script_info@network + \
                    " observational data for " + var0 + \
                    " not available, skipping", verbosity, 1)
        exit
    end if

    ;; File type
    valid_types = (/"PS", "ps", "PDF", "pdf"/)  ; Support for multiple pages
    file_type = getenv("ESMValTool_output_file_type")
    if (ismissing(file_type)) then
        file_type = "PS"
    end if
    if (dimsizes(models@name) .gt. 1 .and. all(file_type .ne. valid_types)) then
        info_output(file_type + " format does not support multiple pages, " \
                    + "switching to ps", verbosity, 1)
        file_type = "PS"
    end if

    ;; Print summary info
    info_output("Network: " + diag_script_info@network, verbosity, 2)
    info_output("Variable: " + var0, verbosity, 2)
    info_output("Models:", verbosity, 2)
    do mID = 0, dimsizes(models@name) - 1
        info_output("  " + models@name(mID), verbosity, 2)
    end do

    ;; Create time-series of co-located model and observations for all stations
    info_output("Reading data for stations and models", verbosity, 2)
    obser_model_pts = obs_model_timeseries(min(toint(models@start_year)), \
                                           max(toint(models@end_year)), \
                                           data, var0)
    obser_model_pts@var = var0
    obser_model_pts@diag_script = (/diag_script/)

    ;; Check for data availability
    if (all(ismissing(obser_model_pts))) then
        info_output(diag_script_info@network + \
                    " observational data not available in the selected " + \
                    "time period, skipping", verbosity, 1)
        exit
    end if

    ;; Calculate yearly mean (if required)
    if (diag_script_info@time_avg.eq."yearly") then
        tmp = time_operations(obser_model_pts(time|:, model|:, station|:), \
                              -1, -1, "average", "yearly", True)
        delete(obser_model_pts)
        obser_model_pts = tmp(model|:, year|:, station|:)
        delete(tmp)
        years = toint(obser_model_pts&year)
    else
        years = cd_calendar(obser_model_pts&time, -1) / 100
    end if

    ;; Optional output to NetCDF
    if (getenv("ESMValTool_write_netcdf").eq."True") then
        new_path = get_output_dir(get_work_dir(), basename(diag_script))
        new_path = new_path + basename(diag_script) + "_" + \
            diag_script_info@network + "_" + var0 + ".nc"
        obser_model_pts@ncdf = new_path
        ncdf_outfile = ncdf_write(obser_model_pts, new_path)
    end if

    ;; For aerosol concentrations, convert units [kg/m3] --> [ug/m3]
    if (isatt(variable_info, "plot_units")) then
        obser_model_pts = convert_units(obser_model_pts, \
                                        variable_info@plot_units)
    end if

    ;; Define path and names of plot files
    plot_dir = getenv("ESMValTool_plot_dir")
    diag_script_base = basename(diag_script)
    output_dir = get_output_dir(plot_dir, diag_script_base)

    ;; Time-series plot ######################################################
    info_output("Time-series plot range: " + min(years) + "-" + max(years), \
                verbosity, 2)
    plot_file = diag_script_info@network + "_timeseries_" + \
        diag_script_info@time_avg + "_" + min(years) + "-" + max(years) + \
        "_" + var0
    plot_file = output_dir + plot_file

    ;; Define line properties (first element is used for observations)
    annots = array_append_record(diag_script_info@network, \
                                 project_style(diag_script_info, "annots"), 0)
    colors = array_append_record("black", \
                                 project_style(diag_script_info, "colors"), 0)
    thicks = array_append_record(2., \
                                 project_style(diag_script_info, "thicks"), 0)
    dashes = array_append_record(0, \
                                 project_style(diag_script_info, "dashes"), 0)
    wks = gsn_open_wks(file_type, plot_file)
    gsn_define_colormap(wks, array_append_record((/"white", "black"/), \
                        colors, 0))

    ;; Loop over stations
    do stID = 0, dimsizes(obser_model_pts&station) - 1

        ;; Timeseries for this station
        data_arr = obser_model_pts(:, :, stID)

        ;; Skip if not observations available
        if (all(ismissing(data_arr(0, :)))) then
            delete(data_arr)
            continue
        end if

        ;; Attach station info as attribute
        data_arr@stname = str_squeeze(str_get_field( \
            obser_model_pts&station(stID), 1, ";"))
        data_arr@stlat  = tofloat(str_get_field( \
            obser_model_pts&station(stID), 2, ";"))
        data_arr@stlon  = tofloat(str_get_field( \
            obser_model_pts&station(stID), 3, ";"))
        data_arr@stalt  = tofloat(str_get_field( \
            obser_model_pts&station(stID), 4, ";"))

        ;; Resources
        data_arr@res_xyLineColors        = colors
        data_arr@res_xyLineThicknesses   = thicks
        data_arr@res_xyDashPatterns      = dashes

        ;; Call plot funciton
        plot = timeseries_station(wks, data_arr, var0)

        delete(data_arr)
        draw(plot)
        frame(wks)

    end do

    ;; Create legend file
    if (diag_script_info@legend_outside) then
        lg_outfile = plot_file + "_legend"
        styles = True
        styles@colors = colors
        styles@dashes = dashes
        styles@thicks = thicks
        create_legend_lines(annots, styles, lg_outfile, "lines")
    end if

    delete(annots)
    delete(colors)
    delete(dashes)
    delete(thicks)
    delete(wks)

    info_output(" gv " + plot_file + "." + file_type, verbosity, 1)

    ;; Scatter plot ##########################################################
    info_output("Scatter-plot range: " + min(years) + "-" + max(years), \
                verbosity, 2)
    plot_file = diag_script_info@network + "_scatter_" + \
        diag_script_info@time_avg + "_" + min(years) + "-" + max(years) + \
        "_" + var0
    plot_file = output_dir + plot_file

    ;; Style options
    colors = project_style(diag_script_info, "colors")
    annots = project_style(diag_script_info, "annots")

    wks = gsn_open_wks(file_type, plot_file)
    gsn_define_colormap(wks, array_append_record((/"white", "black"/), \
                        colors, 0))

    ;; Define arrays
    obser = ndtooned(obser_model_pts(0, :, :))
    obser = where(obser .eq. 0., obser@_FillValue, obser)
    data_arr = new((/2, dimsizes(obser)/), float)
    data_arr!0 = "models"
    data_arr(0, :) = obser
    copy_VarAtts(obser_model_pts, data_arr)

    ;; Loop over models
    do mID = 0, dimsizes(models@name) - 1

        data_arr&models = (/diag_script_info@network, annots(mID)/)

        ;; Transform to 1-D arrays
        model = ndtooned(obser_model_pts(mID + 1, :, :))
        model@units = obser_model_pts@units

        ;; Filter out zeros and missing (skip if all missing)
        model = where(ismissing(obser), model@_FillValue, model)
        if (all(ismissing(model))) then
            delete(model)
            continue
        end if

        ;; Data
        data_arr(1, :) = model

        ;; Resources
        data_arr@res_xyMarkerColor = colors(mID)

        ;; Call plot function
        plot = scatterplot(wks, data_arr, var0, True, True)
        draw(plot)
        frame(wks)
        delete(model)

    end do
    delete(data_arr)
    delete(wks)
    if (isdefined("plot")) then
        delete(plot)
    end if
    delete(obser)
    delete(colors)

    info_output(" gv " + plot_file + "." + file_type, verbosity, 1)

    ;; Map plot ###############################################################
    plot_file = diag_script_info@network + "_map_annualclim_" \
        + min(years) + "-" + max(years) + "_" + var0
    plot_file = output_dir + plot_file
    wks = gsn_open_wks(file_type, plot_file)
    gsn_define_colormap(wks, variable_info@colortable_abs)

    do mID = 0, dimsizes(models@name) - 1

        ;; Calculate time average over consistent time ranges
        tmp = extract_data(mID, data, -1, 0, 0)
        date = cd_calendar(tmp&time, 0)
        t1 = max((/min(toint(date(:, 0))), min(years)/))
        t2 = min((/max(toint(date(:, 0))), max(years)/))
        model = time_operations(tmp, t1, t2, "average", "annualclim", False)
        if (diag_script_info@time_avg.eq."monthly") then
            obser = time_operations(obser_model_pts(0, :, :), t1, t2, \
                                    "average", "annualclim", False)
        else
            obser = dim_avg_n_Wrap(obser_model_pts(0, {t1:t2}, :), 0)
        end if
        delete(tmp)
        delete(date)
        info_output("Map plot range " + models@name(mID) + ": " + \
                    t1 + "-" + t2, verbosity, 2)

        ;; Convert units
        if (isatt(variable_info, "plot_units")) then
            model = convert_units(model, variable_info@plot_units)
        end if

        ;; Plot resources
        model@res_cnLinesOn  =  False
        model@res_cnFillMode =  "AreaFill"
        model@res_cnLevelSelectionMode = "ExplicitLevels"
        model@res_cnLevels = variable_info@contours_abs
        model@res_tiMainString = var0
        model@res_lbBoxLinesOn = True
        if (isatt(diag_script_info, "minlon").and. \
            isatt(diag_script_info, "maxlon").and. \
            isatt(diag_script_info, "minlat").and. \
            isatt(diag_script_info, "maxlat")) then
            ;; Set new map boundaries, make sure aspect ratio is preserved
            dlat = abs(diag_script_info@maxlat - diag_script_info@minlat)
            dlon = abs(diag_script_info@maxlon - diag_script_info@minlon)
            if (dlon.gt.2 * dlat) then
                model@res_mpMinLonF = diag_script_info@minlon
                model@res_mpMaxLonF = diag_script_info@maxlon                
                model@res_mpMinLatF = diag_script_info@minlat                
                model@res_mpMaxLatF = dlon/2. + diag_script_info@minlat
            else
                model@res_mpMinLatF = diag_script_info@minlat
                model@res_mpMaxLatF = diag_script_info@maxlat
                model@res_mpMinLonF = diag_script_info@minlon
                model@res_mpMaxLonF = 2*dlat + diag_script_info@minlon
            end if
            if (model@res_mpMaxLonF.gt.180) then
                model@res_mpCenterLonF = \
                    0.5 * (model@res_mpMinLonF + model@res_mpMaxLonF)
            end if
            model@res_pmTickMarkDisplayMode = "Always"
            model@res_tiMainOffsetYF = -0.035
        end if
        model@res_gsnRightString = t1 + "-" + t2
        model@res_tiMainString = var0 + " (" + annots(mID) + " - " + \
            diag_script_info@network + ")"
        plot = contour_map_ce(wks, model, var0)

        ;; Overlay observations as filled markers
        stat_lat  = tofloat(str_get_field(obser&station, 2, ";"))
        stat_lon  = tofloat(str_get_field(obser&station, 3, ";"))

        add_markers_to_map(wks, plot, stat_lat, stat_lon, obser)

        draw(plot)
        frame(wks)
        delete(model)
        delete(obser)

    end do

    info_output(" gv " + plot_file + "." + file_type, verbosity, 1)

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 2)

end
