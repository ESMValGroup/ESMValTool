;;#############################################################################
;; SouthernOcean_polcon.ncl
;; Author: Daniel Senftleben (DLR, Germany)
;; Contributors: Bert van Ulft (KNMI, the Netherlands)
;; ESMVal project
;;#############################################################################
;; Description
;;    Uses original grid (no regridding).
;;    Draws sea ice concentration color coded and a red contour line at 0.15
;;    (extent).
;;    Panels multiple plots: starts a new page after max_vert*max_hori
;;    Adapted from the SeaIce Polar contor plot
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    month: selected month or annual mean.
;;    region: region to be plotted ( "Arctic" or "Antarctic").
;;    contour_extent: optionally draw a red contour line at 0.15.
;;    max_lat, min_lat: plot ranges for SH & NH.
;;    range_option: select specific or intersecting time range.
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    max_vert, max_hori: allowed dimensions of a panel page
;;                        (number of plots per column/row)
;;    PanelTop: tune to get full title of uppermost row (1 = no top margin)
;;    grid_min, grid_max, grid_step: min, max and step size for color scale
;;    colormap: NCL colour map name
;;    PanelLabelBar: use single label bar per page
;;    showunits: display units in figure title
;;
;; Caveats
;;    no masking of the apparent line of missing data between the two grid
;;    poles yet (for a workaround, see
;;    http://www.ncl.ucar.edu/Applications/Scripts/ice_4.ncl).
;;
;; Modification history
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;;    20150605-A_vanu_be: additional attributes for more formatting control
;;    20150325-A_laue_ax: modified reference tags used for acknowledgements
;;                        (projects, observations, etc)
;;    20140623-A_senf_da: res now passed to plot script as attributes of var.
;;    20140312-A_righ_ma: adapted to new time_operations structure.
;;    20131216-A_vanu_be: generalized a bit so it can plot other vars.
;;    20131212-A_senf_da: modularized to fit into new structure.
;;    20130604-A_gott_kl: implemented.
;;
;;#############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

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
    write_references(diag_script,       \  ; script name
                     "A_gott_kl",          \  ; authors
                     "A_senf_da",          \  ; contributors
                     "D_0001",            \  ; diag_references
                     "",                \  ; obs_references
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
    file_type = getenv("ESMValTool_output_file_type")
    if(ismissing(file_type)) then
        file_type = "PS"
    end if

    ;; check field type
    plot_fields = (/"T2Ms"/)
    ind_f = inlist(field_type0, plot_fields)
    if (all(ismissing(ind_f))) then
        error_msg("f", diag_script, "", "can't use field, " + fields)
    end if
    delete(plot_fields)
    delete(ind_f)

    ;; prepare annotation
    annots = project_style(diag_script_info, "annots")
    smonths = (/"January", "February", "March", "April", \
                "May", "June", "July", "August", \
                "September", "October", "November", "December"/)

    i2months = (/"01", "02", "03", "04", "05", "06",\
                "07", "08", "09", "10", "11", "12"/)

    ;; output file name
    if (diag_script_info@month .eq. "A") then
        i2month = "annual"
    else
        i2month = smonths(toint(diag_script_info@month) - 1)
    end if

    ;;*** prepare plotting ***
    ;; Specify additional elemens for file name
    ;; function get_wks transfers full output file name
    ;; including path via wks@fullname
    if (diag_script_info@range_option .eq. 1) then
        name_string = var0 + "_" + diag_script_info@region + "_" + i2month + \
            "_" + sprinti("%0.4i", max(toint(models@start_year))) + "-" \
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

    ;; open the workstation
    wks = get_wks(0, diag_script, name_string)

    ;; Set contour parameters
    grid_min    =  0.   ; default
    grid_max    =  1.   ; default
    grid_step   =  0.1  ; default
    if (isatt(diag_script_info, "grid_min")) then
        grid_min = tofloat(diag_script_info@grid_min)
    end if
    if (isatt(diag_script_info, "grid_max")) then
        grid_max = tofloat(diag_script_info@grid_max)
    end if
    if (isatt(diag_script_info, "grid_step")) then
        grid_step = tofloat(diag_script_info@grid_step)
    end if
    nsteps = round((grid_max - grid_min) / grid_step, 3) + 1
    cnLevels = fspan(grid_min, grid_max, nsteps)

    ;; Define color table
    ;; "WhiteBlue", "BlueWhiteOrangeRed", "BlueYellowRed" ;"BlWhRe"
    if (isatt(diag_script_info, "colormap")) then
        colormap = diag_script_info@colormap
        gsn_define_colormap(wks, colormap)
        ColorStart = 2
        ColorEnd = -1
    else
        colormap = "WhiteBlue"
        gsn_define_colormap(wks, "WhiteBlue")
        gsn_reverse_colormap(wks)            ; Reverse the color map.
        colors = gsn_retrieve_colormap(wks)  ; retrieve color map for editing
        colors(3, :) = (/1., 0., 0./)  ; red
        colors(4, :) = (/0., 0., 1./)  ; blue
        colors(5, :) = (/1., 1., 1./)  ; black
        gsn_define_colormap(wks, colors)
        ColorStart = 10
        ColorEnd = 240
    end if

    ;; panel parameters
    nvert = 1  ; default
    nhori = 1  ; default
    if (isatt(diag_script_info, "max_vert")) then
        ;; Maximum allowed number of plots per page (vertical)
        nvert = toint(diag_script_info@max_vert)
    end if
    if (isatt(diag_script_info, "max_hori")) then
        ;; Maximum allowed number of plots per page (horizontal)
        nhori = toint(diag_script_info@max_hori)
    end if
    npages = toint(ceil(tofloat(dimsizes(annots)) / tofloat(nvert * nhori)))

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

    ;; define array for collecting plots
    if (npages .gt. 1) then
        ;; size of all stamps is identical only if every page
        ;; contains the same number of plots
        plots = new(npages * nvert * nhori, graphic)
    else
        plots = new(dimsizes(annots), graphic)
    end if

;   start processing ##########################################
    ;; loop over models
    do imod = 0, dimsizes(annots) - 1

        info_output("processing " + annots(imod) + \
                    "  #####################################", \
                    verbosity, 3)

        ;; see ./interface_scripts/data_handling.ncl
        var = extract_data(imod, data, -1, 0, 0)
        if (var@units.eq."%") then
            var = var / 100.  ; Subsequent code works fine without rescaling
            var@units = "%/100"  ; But annotation is hardcoded for %/100
        end if

        ;; Show units in title (default False)
        strUnits = ""
        if (isatt(diag_script_info, "showunits")) then
            strUnits = " (" + var@units + ")"
        end if

        ;; get 2d lat & lon coordinates (see e.g.
        ;; http://www.ncl.ucar.edu/Applications/Scripts/ice_3.ncl),
        ;; but attach them later
        if(all(iscoord(var, (/"lat", "lon"/)))) then
            info_output("SeaIce_polcon.ncl: var is on a lat-lon grid", \
                        verbosity, 3)
            ;; order coordinates
            temp = var(time|:, lon|:, lat|:)
            delete(var)
            var = temp
            delete(temp)
            ;; create 2d lat & lon
            area_temp = map_area(var&lat, var&lon)

            ;; make 2D to be compatible with irregular grid handling below
            area_lat = conform_dims((/dimsizes(var&lon), dimsizes(var&lat)/), \
                                    var&lat, 1)
            area_lon = conform_dims((/dimsizes(var&lon), dimsizes(var&lat)/), \
                                    var&lon, 0)
            area_lat!0 = "i"
            area_lat&i = var&lon
            area_lat!1 = "j"
            area_lat&j = var&lat
            var!1 = "i"
            var!2 = "j"
            copy_VarMeta(area_lat, area_lon)
            delete(area_temp)
        else  ; Irregular grid
            info_output("SeaIce_polcon.ncl: var is not on a lat-lon grid", \
                        verbosity, 3)
            ;; Get lat - lon from another file
            ;; (model specific: ~/reformat_scripts folder)
            area_lon = extract_area(imod, data, "lon", var0)
            area_lat = extract_area(imod, data, "lat", var0)
        end if
        var@lat2d = area_lat
        var@lon2d = area_lon

        ;; time axis
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

        ;; Some attributes for area plot
        if diag_script_info@month .eq. "A" then  ; weighted annual mean
            var_ann@units = "annual mean " + range
        else  ; select a month
            imon = stringtofloat(diag_script_info@month)
            var_ann@units = smonths(floattoint(imon) - 1) + " mean " + range
        end if
        var_ann@long_name = diag_script_info@region\
                          + " "\
                          + var@long_name\
                          + strUnits

        ;; resources to be passed to plot routine
        ; as attributes to var_ann with prefix "res_"
        var_ann@res_lbLabelBarOn          = lblBarOn    ; individual bars
        var_ann@res_cnLevelSelectionMode  = "ExplicitLevels"
        var_ann@res_cnLevels              = cnLevels
        var_ann@res_gsnSpreadColors       = True        ; use full colormap
        var_ann@res_gsnSpreadColorStart   = ColorStart  ; color to start
        var_ann@res_gsnSpreadColorEnd     = ColorEnd    ; color to end
        var_ann@res_cnLinesOn             = False
        var_ann@res_cnLineLabelsOn        = False
        var_ann@res_tiMainString          = annots(imod)
        var_ann@res_gsnRightString        = var_ann@units
        var_ann@res_gsnDraw               = False       ; Draw and Frame is
        var_ann@res_gsnFrame              = False       ; called in panelling

        ; pass additional meta data
        var_ann@diag_script = (/diag_script/)

        ;; draw the plot
        plots(imod) = contour_map_polar(wks, var_ann, variables(0))

        ;; clean up
        delete(var)
        delete(area_lon)
        delete(area_lat)
        delete(var_ann)

    end do  ; imod

    ;; Draw the panel
    pres = True
    pres@gsnPanelLabelBar = PanelLabelBar  ; 1 bar per panel?
    outfile = panelling(wks, plots, nvert, nhori, pres)

    info_output(" gv " + outfile + " &", verbosity, 1)

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)

end
