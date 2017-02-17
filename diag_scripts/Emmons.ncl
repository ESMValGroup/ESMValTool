;##############################################################################
;# COMPARING TRACE GAS PROFILES OF Emmons et al. (2000) TO MODELS             #
;# Author: Klaus Gottschaldt (DLR, Germany)                                   #
;# ESMVal project                                                             #
;##############################################################################
;#                                                                            #
;# Description:                                                               #
;#    * Chooses observational profile data that match some selection criteria #
;#      (time, region, species) from a directory (specified in cfg_Emmons)    #
;#    * The same selection criteria are applied to all models specified in    #
;#      namelist_Emmons (as far as possible)                                  #
;#    * One plot "mixing ratio vs pressure altitude" is created for each      #
;#      observational file, with all trace gas profiles regridded to          #
;#      the observational levels                                              #
;#    * Whiskers for min/stddev/max or percentiles (5,25,75,95)%              #
;#      are added to the observational (reference) profile                    #
;#    * Output is via the panelling routine: Choose max_vert = 1 and          #
;#      max_hori = 1 in cfg_Emmons to get individual plots                    #
;#    * Reference: Emmons et al.: Data composites of airborne observation of  #
;#                 tropospheric ozone and its precursors,                     #
;#                 J. Geophys. Res., 105, 20497-20538, 2000.                  #
;#                                                                            #
;# Required variable_info attributes (variable specific): none                #
;#                                                                            #
;# Required diag_script_info attributes (diagnostics specific):               #
;#    * See cfg_Emmons for details. Mandatory attributes are:                 #
;#    * @obsfiles         ;obs files to be scanned (wildcards allowed)        #
;#    * @obsvar_mean or @obsvar_median (depending on selection of @quantity)  #
;#         ;name of the base parameter in the obs files                       #
;#    * @obsatt_campaign  ;for plot annotation                                #
;#    * @obsatt_location  ;for plot annotation                                #
;#    * @obsatt_period    ;for selecting matching time of year in model data  #
;#    * @obsatt_latmin    ;for selecting matching region in the model data    #
;#    * @obsatt_latmax    ;for selecting matching region in the model data    #
;#    * @obsatt_lonmin    ;for selecting matching region in the model data    #
;#    * @obsatt_lonmax    ;for selecting matching region in the model data    #
;#                                                                            #
;# Caveats:                                                                   #
;#    * statistics for single model modus not implemented yet:                #
;#      look for "single" to see an outline on how to do it                   #
;#    * only mean fully implemented yet,                                      #
;#      see comments near int2p_n_Wrap to add median                          #
;#    * performance metrics not implemented yet                               #
;#                                                                            #
;# Modification history:                                                      #
;#    * 20151027-A_laue_ax: moved call to 'write_references' to the beginning #
;#                        of the code                                         #
;#    * 20150325-A_laue_ax: modified reference tags used for acknowledgements #
;#                        (projects, observations, etc)                       #
;#    * 20140226 written by Klaus Gottschaldt                                 #
;##############################################################################

;##############################################################################
;# load external NCL code, if needed                                          #
;##############################################################################

;; Temporary file with current namelist configuration
;; (updated by Python on the fly)
load "./interface_data/ncl.interface"

; Auxiliary NCL routines
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
      ;; for int2p_n_Wrap, array_append_record
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"
load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/set_operators.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./plot_scripts/ncl/xy_line.ncl"
;load "./plot_scripts/ncl/legends.ncl"
load "./diag_scripts/aux/Emmons/aux_Emmons.ncl"
load "./diag_scripts/lib/ncl/scaling.ncl"  ; for scalfac
load "./plot_scripts/ncl/aux_plotting.ncl"  ; for get_wks & get_outfile_name

begin
    ;;#########################################################################
    ;;# Fetch general parameters, set in namelist_Emmons.xml,                 #
    ;;# passed via environment variables by python code                       #
    ;;#########################################################################
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)

    ;; 'diag_script' and 'variables' are fetched from the
    ;; above 'ncl.interface' file
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    ;; 'models', 'variables' and 'field_types' are fetched from the
    ;; above 'ncl.interface' file
    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    ;;#########################################################################
    ;;# References Information                                                #
    ;;#########################################################################
    write_references(diag_script,        \  ; script name
                     "A_gott_kl",        \  ; authors
                     "",                 \  ; contributors
                     "D_0000",           \  ; diag_references
                     "E_emmons00jgr",             \  ; obs_references
                     (/"P_esmval"/))          ; proj_references

    do idx_mod = 0, dim_MOD - 1
        data = read_data(idx_mod, var0, field_type0)
    end do
end

begin
    ;; Load variable_info
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin
    ;;#########################################################################
    ;;# Select files with observational data,                                 #
    ;;# according to diag_script_info                                         #
    ;;#########################################################################
    do iobs = 0, dimsizes(diag_script_info@obsfiles) - 1
        info_output("Scanning " + diag_script_info@obsfiles(iobs) + \
                    " for valid observations:", verbosity, 5)
        if(iobs.eq.0) then
            flist = systemfunc("ls " + diag_script_info@obsfiles(iobs))
        else
            temp1 = flist
            delete(flist)
            temp2 = systemfunc("ls " + diag_script_info@obsfiles(iobs))
            flist = array_append_record(temp1, temp2, 0)
            delete(temp1)
            delete(temp2)
        end if
    end do
    ;; Loop over all files and compare
    ;; diag_script_info@roi to @obsatt_region
    ;; & diag_script_info@obsatt_period to @poi.
    ;; Skip files that do not match.
    do ifile = 0, dimsizes(flist) - 1
        info_output(" reading " + flist(ifile), verbosity, 2)
        f = addfile(flist(ifile), "r")

        ;; check target variable
        vars = getfilevarnames(f)
        test_var = var0 + diag_script_info@obsvar_mean

        if(.not.inlist(test_var, vars)) then
            info_output(test_var + " not in " + flist(ifile) +\
                        ": skipping", verbosity, 5)
            flist(ifile) = "skipped"
        end if
        delete(vars)

        ;; check period
        if(flist(ifile).ne."skipped") then
            period_in = ncdf_att(f, diag_script_info@obsatt_period)
            period_dat = get_period(period_in)
                         ;; function in aux/Emmons/aux_Emmons.ncl
                         ;; returns array (/3,2/)
                         ;; with strings yyyy, mm, dd
            if(isatt(diag_script_info, "poi")) then
                period_of_interest = diag_script_info@poi
            else  ; default poi for obs: use everything
                period_of_interest = (/1900, 2100/)
                info_output("diag_script_info@poi not defined: " +\
                            "using everything", verbosity, 5)
            end if
            result_poi = poi((/period_dat(0, 0), period_dat(0, 1)/),\
                             period_of_interest)  ; aux/Emmons/aux_Emmons.ncl
            if(isatt(diag_script_info, "poi_match")) then
                poi_match = diag_script_info@poi_match
            else  ; default: obs period and poi must have some overlap
                poi_match = "overlap"
            end if
            if(.not.inlist(diag_script_info@poi_match, result_poi)) then
                    ;; interface_scripts/auxiliary.ncl
                info_output("  period of interest not matching:" +\
                            " skipping", verbosity, 5)
                flist(ifile) = "skipped"
            end if
            delete(result_poi)
        end if

        ;; check region
        if(flist(ifile).ne."skipped") then
            latmin = ncdf_att(f, diag_script_info@obsatt_latmin)
            latmax = ncdf_att(f, diag_script_info@obsatt_latmax)
            lonmin = ncdf_att(f, diag_script_info@obsatt_lonmin)
            lonmax = ncdf_att(f, diag_script_info@obsatt_lonmax)
            latlon_dat = tofloat((/latmin, latmax, lonmin, lonmax/))
            if(isatt(diag_script_info, "roi")) then
                latlon_roi = tofloat(diag_script_info@roi)
            else  ; default: global
                latlon_roi = (/-90., 90., 0., 360./)
                info_output("diag_script_info@roi not defined: " +\
                            "using global domain", verbosity, 5)
            end if
            result_roi = roi(latlon_dat, latlon_roi)
                         ;; diag_scripts/lib/ncl/latlon.ncl
            if(isatt(diag_script_info, "roi_match")) then
                roi_match = diag_script_info@roi_match
            else  ; default: center of the obs region must be in roi
                roi_match = "center"
            end if
            if(.not.inlist(roi_match, result_roi)) then
                    ;; interface_scripts/auxiliary.ncl
                info_output("  region of interest not matching:" +\
                            " skipping", verbosity, 5)
                flist(ifile) = "skipped"
            end if  ; roi
            delete([/result_roi/])
        end if

        delete([/f/])
    end do  ; flist

    ;;#########################################################################
    ;;# Are there matching obs files?                                         #
    ;;#########################################################################
    imatching = ind(flist.ne."skipped")
    if(.not.all(ismissing(imatching))) then
        files2use = flist(imatching)
        info_output(" Using " + dimsizes(files2use) + " of " +\
                    dimsizes(flist) + " files.", verbosity, 5)

        ;;#####################################################################
        ;;# Set up wks and panelling, plot file name                          #
        ;;#####################################################################
        ;; Do not overwrite existing plots by appending a new number
        if(isatt(diag_script_info, "aux_info")) then
            aux_info = diag_script_info@aux_info
        else
            aux_info = sprinti("%0.4i", diag_script_info@poi(0)) + "-"\
                + sprinti("%0.4i", diag_script_info@poi(1))
        end if
        i = 0
        name_string = var0 + "_" + aux_info
        name_string_n = name_string
        fullname = get_outfile_name(diag_script, name_string_n)
        do while(isfilepresent(fullname))
            name_string_n = name_string + "_" + tostring(i)
            fullname = get_outfile_name(diag_script, name_string_n)
            i = i + 1
        end do
        delete([/fullname, i/])

        ;; wks for profile plots
        wks = get_wks("dummy_for_wks", diag_script, name_string_n)

        ;; legend
        if(isatt(diag_script_info, "legend")) then
            legend = diag_script_info@legend
        else
            legend = "separate"
        end if
        if(legend.eq."separate") then  ; wks for separate legend
            ;wks_leg = get_wks("dummy_for_wks", diag_script, \
            ;                  name_string_n + "_legend")
            base = get_outfile_name(diag_script, name_string_n)
            base = str_sub_str(base, get_file_suffix(base, 0), "")
            name_string_leg = base + "_legend"
            delete(base)
        end if

        ;; panel parameters: number of plots determined by files2use
        nvert = 1  ; default
        nhori = 1  ; default
        if (isatt(diag_script_info, "max_vert")) then
            nvert = toint(diag_script_info@max_vert)
        end if
        if (isatt(diag_script_info, "max_hori")) then
            nhori = toint(diag_script_info@max_hori)
        end if
        npages = toint(ceil(tofloat(dimsizes(files2use)) \
                       / tofloat(nvert * nhori)))

        ;; define array for collecting plots
        if (npages .gt. 1) then
        ;; size of all stamps is identical only if every page
        ;; contains the same number of plots
            plots = new(npages * nvert * nhori, graphic)
        else
            plots = new(dimsizes(files2use), graphic)
        end if

        ;; diag_script-specific resources to be passed to plot routine
        res                       = True

        ;;#####################################################################
        ;;# Process matching observational files                              #
        ;;#####################################################################
        ;; Target quantity: mean or median?
        targetq = "mean"  ; default
        if(isatt(diag_script_info, "quantity")) then
            targetq = diag_script_info@quantity
        end if
        info_output(" Target quantity = " + targetq, verbosity, 10)
        if(targetq.eq."mean") then
            quantity = (/"N", "mean", "stddev", "min", "max"/)
        else if(targetq.eq."median") then
            quantity = (/"N", "median", "5", "25", "75", "95"/)
            error_msg("f", diag_script, "", "median not fully implemented yet")
        else
            error_msg("f", diag_script, "", "diag_script_info@quantity = " + \
                      targetq + " not implemented. Disable to " + \
                      "use the default (mean)")
        end if
        end if

        ;; loop over matching files
        do ifile = 0, dimsizes(files2use) - 1
            info_output(" processing " + files2use(ifile), \
                        verbosity, 2)
            f = addfile(files2use(ifile), "r")

            ;; fetch (again) obs region and period
            latmin = ncdf_att(f, diag_script_info@obsatt_latmin)
            latmax = ncdf_att(f, diag_script_info@obsatt_latmax)
            lonmin = ncdf_att(f, diag_script_info@obsatt_lonmin)
            lonmax = ncdf_att(f, diag_script_info@obsatt_lonmax)
            period_in = ncdf_att(f, diag_script_info@obsatt_period)
            period_dat = get_period(period_in)  ; aux_Emmons

            ;; fetch observational pressure axis
            ;;  (to be used as reference)
            plev = f->plev
            plev_status = isMonotonic(plev)
            if(plev_status.ne.-1) then
                error_msg("f", diag_script, "", "plev not monotonically " + \
                          "decreasing: status=" + plev_status)
            end if

            ;; Create field var to hold results
            ;; (plev is target coordinate, also for models).
            ;; Dimensions: plev, dim_MOD+obs, number of possible
            ;;   diag_script_info@obsvar_*:
            ;;   N,mean,stddev,min,max OR N,5,25,median,75,95
            var = new((/dimsizes(plev), dim_MOD + 1, dimsizes(quantity)/),\
                      float)
            var!0 = "plev"
            var&plev = plev
            var!1 = "models"
            ;; refmodel included via "if", but hardcoded here for "else"
            if(isatt(diag_script_info, "annots")) then
                var&models = diag_script_info@annots
            else
                annots = project_style(diag_script_info, "annots")
                ;var&models = array_append_record("ref_" + \
                ;             diag_script_info@quantity, annots, 0)
                var&models = array_append_record("Emmons", annots, 0)
                delete(annots)
            end if
            var!2 = "quantity"
            var&quantity = quantity
            if(.not.isatt(var, "_FillValue")) then
                var@_FillValue = default_fillvalue("float")
            end if
            var = var@_FillValue

            ;; fetch observational data
            ;; and store as first (index 0 hardcoded) entry in var
            do iq = 0, dimsizes(quantity) - 1
                q = "obsvar_" + quantity(iq)
                if(isatt(diag_script_info, q)) then
                    temp = var0 + diag_script_info@$q$
                    var(:, 0, iq) = f->$temp$
                end if
            end do
            delete([/temp, iq, q/])

            ;; exclude data < 0
            temp = ndtooned(var(:, 0, :))
            itemp = ind(temp.lt.0.)
            if(.not.all(ismissing(temp))) then
                info_output("info: Assigning zero to " + dimsizes(itemp) + \
                            " negative observational values.", verbosity, 2)
                var(:, 0, :) = where(var(:, 0, :).lt.0., 0., var(:, 0, :))
            end if
            delete([/itemp, temp/])

            ;;#################################################################
            ;;# Process model data                                            #
            ;;#################################################################
            do imod = 0, dim_MOD - 1  ; loop over models

                info_output("processing " + models@name(imod), verbosity, 2)

                A0 = extract_data(imod, data, -1, 0, 0)
                     ;; interface_scripts/data_handling.ncl

                if (verbosity .ge. 10) then
                    info_output(" ++++++++++++++ Variable summary " +\
                                "of original data (A0): ", \
                                verbosity, 10)
                    printVarSummary(data)
                    printVarSummary(A0)
                end if

                ;; Check dimensions
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

                ;; Extract region corresponding to obs ...
                if(targetq.eq."mean") then
                    ;; ... and do area weighted average
                    profile2d = area_operations(A0, latmin, latmax, lonmin,\
                                                lonmax, "average", True)
                                ;; latlon.ncl
                else if(targetq.eq."median") then
                    ;; ... and calculate median
;; see https://www.ncl.ucar.edu/Applications/Scripts/box_5.ncl
;; for percentile calculation
                    error_msg("f", diag_script, "", "median not fully " + \
                              "implemented yet")
                else
                    error_msg("f", diag_script, "", \
                              "diag_script_info@quantity = " + \
                              targetq + " not implemented. Disable to " + \
                              "use the default (mean)")
                end if
                end if
                ;; Extract & average time steps that match obs time of year
                ;; Weighting depends on time step size of the data
                ;;   (see function toy_weights in aux_Emmons.ncl).
                ;; period_dat = (/3,2/) with strings yyyy, mm, dd
                start = period_dat(1:2, 0)
                final = period_dat(1:2, 1)
                wgt_time = toy_weights(start, final, profile2d&time)
                ;; averaging all non-missing values
                ;; weights are being normalized automaticlly
                profile1d = dim_avg_wgt_n_Wrap(profile2d, wgt_time, 1, 0)

;; For single model modus it might be desired to get
;; statistical information
;; (N,mean,stddev,min,max,5,25,median,75,95) for model data.
;; To implement this, area & time averages must be disabled
;; above. Instead the following might be a strategy to get
;; this functionality:
;; * field = area_operations(A0,latmin,latmax,lonmin,lonmax,"extract",True)
;; * wgt_area = map_area(field&lat, field&lon)  ... area weights
;; * wgt_time ...time weights
;; * consider wgt_area * wgt_time as number of "measurements" per cell
;;   and do statistics accordingly

                ;; adjust missing values to obs
                profile1d = where(ismissing(profile1d),\
                                  var@_FillValue, profile1d)
                profile1d@_FillValue = var@_FillValue

                ;; check consistency of obs pressure levels
                if(isMonotonic(profile1d&plev).eq.1) then  ; need reversing
                    profile1d = profile1d(::-1)  ; must be decreasing
                else if(isMonotonic(profile1d&plev).eq.0) then  ; crash
                    error_msg("f", diag_script, "", "plev_dat is not " + \
                              "monotonic")
                end if
                end if

                ;; interpolate to obs pressure levels
                if(isatt(diag_script_info, "interpolation")) then
                    interpol_flag = diag_script_info@interpolation
                    info_output("Interpolation set by diag_script_info",\
                                verbosity, 5)
                else  ; default
                    valid_plev = plev(ind(.not.ismissing(plev)))
                    nvalid = dimsizes(valid_plev)
                    if(max(profile1d&plev).gt.valid_plev(1) .and.\
                       min(profile1d&plev).lt.valid_plev(nvalid - 2)) then
                        interpol_flag = -1
                        info_output("Linear interpolation with extrapolation",\
                                    verbosity, 5)
                    else
                        interpol_flag = 1
                        info_output("Linear interpolation, no extrapolation",\
                                    verbosity, 5)
                    end if
                end if
                var(:, imod + 1, 1) = int2p_n_Wrap(profile1d&plev, profile1d,\
                                                   plev, interpol_flag, 0)
                ;; Fill only mean/median column of quantity coordinate (yet).
                ;; Other columns might be needed for single model modus or
                ;; performance metrics.

                ;; clean up
                delete([/A0, profile1d, profile2d, wgt_time/])

            end do  ; imod

            ;;#################################################################
            ;;# Add general Metadata to var                                   #
            ;;#################################################################
            ;; Add as attributes without prefix
            if(isatt(var, "diag_script")) then  ; Add to existing entries
                temp = var@diag_script
                delete(var@diag_script)
                var@diag_script = array_append_record(temp, (/diag_script/), 0)
                delete(temp)
            else  ; Add as new attribute
                var@diag_script = (/diag_script/)
            end if
            var@var = var0  ; Overwrite existing entry
            if(isatt(variable_info, "long_name")) then
                var@var_long_name = variable_info@long_name
            end if
            if(isatt(variable_info, "units")) then
                var@var_units = variable_info@units
            end if
            if(isatt(diag_script_info, "ptop")) then
                var@ptop = diag_script_info@ptop
            end if
            if(isatt(diag_script_info, "pbot")) then
                var@pbot = diag_script_info@pbot
            end if
            if(isatt(diag_script_info, "zoom")) then
                var@zoom = diag_script_info@zoom
            end if
            if(isatt(diag_script_info, "font")) then
                var@font = diag_script_info@font  ; used for several
            else                                  ; plot ressources
                var@font = 25
            end if
            if(isatt(variable_info, "FuncCode")) then
                var@FuncCode = variable_info@FuncCode
            end if
            if(isatt(variable_info, "short_name")) then
                var@short_name = variable_info@short_name
            end if
            var@Refmodel = var&models(0)  ; as hardcoded throughout this script
            var@legend = legend

            ;;#################################################################
            ;;# Update Performance metrics                                    #
            ;;#################################################################
            ;; to be implemented

            ;; Convert units if appropriate
            if (isatt(variable_info, "plot_units")) then
                var(:, :, 1:) = \
                    convert_units(var(:, :, 1:), variable_info@plot_units)
            end if

            ;; fetch obs parameters for annotation; styles etc.
            campaign = ncdf_att(f, diag_script_info@obsatt_campaign)
            location = ncdf_att(f, diag_script_info@obsatt_location)
            period = ncdf_att(f, diag_script_info@obsatt_period)

            ;; pass specific plot ressources via attributes of var
            var@res_tiMainString      = location
            var@res_gsnLeftString     = campaign
            var@res_gsnRightString    = period
            var@res_gsnDraw           = False  ; Do not draw yet
            var@res_gsnFrame          = False  ; Don't advance frame.

            ;; set style parameters:
            ;; reference model (Emmons) must be included if taken
            ;; from diag_script_info, but is added (hardcoded) here otherwise
            if(isatt(diag_script_info, "colors")) then
                var@res_xyLineColors = diag_script_info@colors
            else
                colors = project_style(diag_script_info, "colors")
                var@res_xyLineColors = array_append_record("black", colors, 0)
                delete(colors)
            end if
            if(isatt(diag_script_info, "dashes")) then
                var@res_xyDashPatterns = diag_script_info@dashes
            else
                dashes = project_style(diag_script_info, "dashes")
                var@res_xyDashPatterns = array_append_record(0, dashes, 0)
                delete(dashes)
            end if
            if(isatt(diag_script_info, "thicks")) then
                var@res_xyLineThicknesses = diag_script_info@thicks
            else
                thicks = tofloat(project_style(diag_script_info, "thicks"))
                var@res_xyLineThicknesses = \
                    array_append_record(max(thicks) * 2, thicks, 0)
                delete(thicks)
            end if

            ;;#################################################################
            ;;# Optional output to netCDF                                     #
            ;;#################################################################
            if(isatt(diag_script_info, "ncdf")) then
                var@ncdf = diag_script_info@ncdf
                ncdf_outfile = ncdf_write(var, diag_script_info@ncdf)
                  ;; Function in ~/interface_scripts/auxiliary.ncl
            end if

            ;;#################################################################
            ;;# Draw ONE separate legend, using var attributes                #
            ;;#################################################################
            if(legend.eq."separate" .and. ifile.eq.0) then
                leg = True
                leg@txFontQuality     = "High"
                leg@txFont            = var@font
                leg@txFontHeightF     = 0.02  ; used in
                leg@diag_script = diag_script
                leg@annots = var&models
                leg@colors = var@res_xyLineColors
                leg@dashes = var@res_xyDashPatterns
                leg@thicks = var@res_xyLineThicknesses
;                plot = legend_lines(wks_leg, leg, "leg")
                ;pres = True
                ;outfile = panelling(wks_leg, plot, 1, 1, pres)
                ;info_output(" Wrote " + outfile, verbosity, 1)
            end if

            ;;#################################################################
            ;;# Draw the plot                                                 #
            ;;#################################################################
            plots(ifile) = profile_plev(wks, var, var0)

            ;; clean up
            delete([/f, plev, var/])

        end do  ; ifile (files2use)
        delete([/files2use/])

        ;;#####################################################################
        ;;# Draw the panel(s)                                                 #
        ;;#####################################################################
        pres = True
        outfile = panelling(wks, plots, nvert, nhori, pres)
        info_output(" Wrote " + outfile, verbosity, 1)

        create_legend_lines(leg@annots, leg, name_string_leg, "lines")

    else  ; no matching obs files

        info_output("No observational files for " + var0 + " found in "\
                    + str_concat(diag_script_info@obsfiles + ", "), \
                    verbosity, 2)
        info_output("Consider to modify poi & roi in the config " +\
                    "file (cfg_Emmons).", verbosity, 2)
        info_output("No plots have been produced for this " +\
                    "diagnostic.", verbosity, 2)

    end if  ; matching obs files or not

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)

end
