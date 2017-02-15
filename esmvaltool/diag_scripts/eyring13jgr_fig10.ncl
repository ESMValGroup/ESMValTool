;##############################################################################
;# eyring13jgr_fig10.ncl                                                      #
;##############################################################################
;# Author: Irene Cionni (ENEA, Italy)                                         #
;# ESMVal project PI: Veronika Eyring (DLR, Germany)                          #
;##############################################################################
;#                                                                            #
;# Description:                                                               #
;#     Calculates trends in annual mean near-global (82.5°S-82.5°N) ozone at  #
;#     50 hPa and temperature at 100 hPa, September-October-November-December #
;#     (SOND) ozone at 50 hPa over Antarctica (60°S-90°S), October-November-  #
;#     December-January (ONDJ) temperature at 100 hPa over Antarctica         #
;#     (60°S-90°S), DJF SH jet position at 850 hPa, and DJF upper             #
;#     tropospheric tropical (30°S-30°N) temperatures at 250 hPa. The trends  #
;#     are calculated over 1979-2005 for the past and over 2006-2050 for the  #
;#     future.                                                                #
;#                                                                            #
;# Required diag_script_info attributes (diagnostics specific):               #
;#     e13fig10_latrange: min. and max. latitude of the regions (n-element    #
;#         array of 2-element pairs, e.g., (/(/-30,30/)/));                   #
;#         one pair of latitudes is required for each season (see below)      #
;#     styleset: style set ("DEFAULT", "CMIP5")                               #
;#     e13fig10_season: seasons (n-element array of strings, "ANN", "JAN",    #
;#         "FEB", "MAR", "DJF", "SON", etc.)                                  #
;#     e13fig07_lev: vertical level (in hPa)                                  #
;#     plot_number_string: string used for plot labeling / sub-figure         #
;#         (e.g., “(a)”)                                                      #
;#                                                                            #
;# Optional diag_script_info attributes (diagnostics specific):               #
;#     list_chem_mod: list of models in the group “chem” (array of strings,   #
;#         default = (/”All”))                                                #
;#     list_nochem_mod: list of models in the group “nochem” (array of        #
;#         strings, default = (/"None"/))                                     #
;#                                                                            #
;# Caveats:                                                                   #
;#                                                                            #
;# Modification history:                                                      #
;#    * 20151209-A-laue_ax: code clean-up, added call to write_references,    #
;#                          several bugfixes for processing zonal wind (ua)   #
;#                          and handling missing variable attribute "units"   #
;#    * 201510??-A-cion_ir: written                                           #
;#                                                                            #
;##############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"
load "./diag_scripts/lib/ncl/set_operators.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/regridding.ncl"
load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/ensemble.ncl"

;************************************************************************
;************************************************************************

begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output("plot- " + diag_script + ".ncl (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    dim_MOD = dimsizes(models@name)
    var0 = variables(0)
    field_type0 = field_types(0)  
     print(diag_script_info) 
    colors=project_style(diag_script_info, "colors")
    dashes=project_style(diag_script_info, "dashes")
    thicks=tofloat(project_style(diag_script_info, "thicks")) 
    do idx_mod = 0, dim_MOD - 1
        data0 = read_data(idx_mod,var0, field_type0)
    end do
end

begin
    vardeffile = "interface_data/" + variables + "_info.tmp"
    loadscript(vardeffile)
end
begin
    ;references Information
    write_references(diag_script,      \  ; script name
                     "A_cion_ir",      \  ; authors
                     "",               \  ; contributors
                     "D_eyring13jgr",  \  ; diag_references
                     "",               \  ; obs_references
                     "P_esmval")          ; proj_references

    dim_EXP = 1  ; dimsizes of different experiments
    label_exp = models@experiment(0)
    do n = 1, dimsizes(models@experiment) - 1
        if (models@experiment(n).ne.models@experiment(n - 1)) then
            dim_EXP = dim_EXP + 1
            ltmp = array_append_record(label_exp, models@experiment(n), 0)
            delete(label_exp)
            label_exp = ltmp
            delete(ltmp)
        end if
    end do

    enssize = dimsizes(ind(models@experiment(0).eq.models@experiment))
    yearsize = max(toint(models@end_year) - toint(models@start_year) + 1)
    dum = new((/dimsizes(diag_script_info@e13fig10_season), \
                (dim_EXP) * 3/), graphic)
    error_bar = new((/dimsizes(diag_script_info@e13fig10_season), \
                (dim_EXP) * 3/), graphic)

    if ((isatt(diag_script_info, "list_chem_mod"))) then
        list_chem_mod = diag_script_info@list_chem_mod
    else
        list_chem_mod = models@name
    end if

    if ((isatt(diag_script_info,"list_nochem_mod"))) then
        list_nochem_mod = diag_script_info@list_nochem_mod
    else
        list_nochem_mod = "None"
    end if

    if (isvar("latrange")) then
        delete(latrange)
    end if

    if (isatt(diag_script_info, "e13fig10_latrange")) then
        latrange = diag_script_info@e13fig10_latrange
        if (dimsizes(dimsizes(latrange)).le.1) then
            latrange1 = new((/1, 2/), typeof(latrange))
            latrange1(0, :) = latrange
            delete(latrange)
            latrange = latrange1
            delete(latrange1)
        end if
    else
        print("fatal:plot requires diag_script_info attribute 'e13fig10_latrange'")
        status_exit(1)
    end if

    if (isatt(diag_script_info, "e13fig10_season")) then
        season = diag_script_info@e13fig10_season
    else
        if (isStrSubset(field_type0, "z") .eq. False) then
            print("fatal:plot requires diag_script_info attribute 'e13fig10_season'")
            status_exit(1)
        end if
    end if

    if (isatt(diag_script_info, "e13fig10_lev")) then
        lev = toint(diag_script_info@e13fig10_lev)
    end if

    if ((isatt(diag_script_info, "plot_number"))) then
        plot_num = new(dimsizes(diag_script_info@e13fig10_season), "string")
        plot_num(:) = diag_script_info@plot_number
    else
        plot_num = integertochar(ispan(97, 122, 1))
    end if

    diag_script_base = basename(diag_script)
    output_dir = get_output_dir(plot_dir, diag_script_base)
    aux_title_info = plot_num(0)
    output_filename = interface_get_figure_filename(diag_script_base,\
        var0, lev, aux_title_info, -1)
    file_type = getenv("ESMValTool_output_file_type")
    if (ismissing(file_type)) then
        file_type = "ps"
    end if
    file_type = str_lower(file_type)
    color_type = "COLOR"
    output_file_path = output_dir + output_filename
    wks = gsn_open_wks(file_type, output_file_path)
    plot = new(dimsizes(diag_script_info@e13fig10_season), graphic)

    ;--------------------------------------------------------------------------
    ;--------------------------------------------------------------------------

    do ipan = 0, dimsizes(season) - 1  ; panel loop
        chem_t       = new((/dim_EXP, dimsizes(list_chem_mod) + 1, yearsize/), "float")
        chem_trend   = new((/dim_EXP, dimsizes(list_chem_mod) + 1/), "float")
        nochem_t     = new((/dim_EXP, dimsizes(list_nochem_mod) + 1, yearsize/), "float")
        nochem_trend = new((/dim_EXP, dimsizes(list_nochem_mod) + 1/), "float")
;        trend        = new((dim_EXP) * 3, "float")
;        trend_std    = new((dim_EXP) * 3, "float")

        ;----------------------------chem------------------------------------

        iexp = 0
        nochem = 0
        chem = 0

        time  = new((/dim_EXP, yearsize/), float)

        do imod = 0, dim_MOD - 1
            A0 = extract_data(imod, data0, -1, 0, 0)
            if (max(A0&plev).ge.100000) then
                A0&plev = A0&plev / 100
            end if
            ttmp = ispan(toint(models@start_year(imod)), \
                toint(models@end_year(imod)), 1)

            if (var0.ne."ua") then
                dummy1 = area_operations(A0(time|:, \
                    {plev|lev(ipan)}, lat|:, lon|:), latrange(ipan, 0), \
                    latrange(ipan, 1), 0., 360., "average", True)
            else
                ; A-laue_ax: calculate average over longitudes
                ; bugfix: the selection of the vertical level reduces the dimension
                ;         of the array by one ---> the last dimension (longitudes)
                ;         is therefore at index=2 and not at index=3
                tmp_pos = dim_avg_n_Wrap(A0(time|:, {plev|lev(ipan)}, \
                    {lat|latrange(ipan, 0):latrange(ipan, 1)}, lon|:), 2)
;                    {lat|latrange(ipan, 0):latrange(ipan, 1)}, lon|:), 3)
                ; A-laue_ax (bugfixes):
                ;   1. dummy_ua is not defined --> use ua instead
                ;   2. copy time dimension from tmp_pos for later call to time_operations
;                dummy1 = new(dimsizes(dummy_ua&time), typeof(A0&lat))
                dummy1 = new(dimsizes(A0&time), typeof(A0&lat))
                copy_VarCoords_1(tmp_pos, dummy1)
                do tt = 0, dimsizes(A0&time) - 1
                    ; A-laue_ax (bugfix): inserted missing time index tt
;                    dummy1(tt) = tmp_pos&lat(maxind(tmp_pos))
                    dummy1(tt) = tmp_pos&lat(maxind(tmp_pos(tt,:)))
                end do
                delete(tmp_pos)
            end if

            ; convert to ppmv
            if (var0.eq."tro3") then
                dummy1 = dummy1 * 10^-3
                variable_info@units = "ppm"
            end if

            if (var0.eq."vmro3") then
                dummy1 = dummy1 * 10^6
                variable_info@units = "ppm"
            end if

            if ((isatt(variable_info, "units"))) then
            else
                if (var0.eq."ta") then
                    variable_info@units = "K"
                end if
            end if

            ; ua is used to calculate the jet position (= latitude of
            ; the zonally averaged wind maximum) --> use "degree" as unit
            ; string
            if (var0.eq."ua") then
                variable_info@units = "deg"
            end if

            delete(A0)

            if (season(ipan).ne."yearly") then
                tmp = time_operations(dummy1, \
                    toint(models@start_year(imod)), toint(models@end_year(imod)), \
                    "extract", season(ipan), True)
                dummy2 = time_operations(tmp, toint(models@start_year(imod)), \
                    toint(models@end_year(imod)), "average", "yearly", True)
                delete(tmp)
            else
                dummy2 = time_operations(dummy1, toint(models@start_year(imod)), \
                    toint(models@end_year(imod)), "average", season(ipan), True)
            end if

            dummy3 = dummy2 - dim_avg_n(dummy2, 0)

            if (imod.ge.1) then
                if (models@experiment(imod).ne.models@experiment(imod - 1)) then
                    iexp = iexp + 1
                    nochem = 0
                    chem = 0
                    tetmp = ispan(toint(models@start_year(imod)), \
                        toint(models@end_year(imod)), 1)
                    time(iexp, :dimsizes(ttmp) - 1) = tetmp
                end if
            end if

            if (any(models@name(imod).eq.list_chem_mod)) then
                chem_t(iexp, chem,:dimsizes(dummy3) - 1) = tofloat(dummy3)
                chem_trend(iexp, chem) = tofloat(regline(ttmp, dummy3) * 10) 
                chem = chem + 1
            end if

            if (any(models@name(imod).eq.list_nochem_mod)) then
                nochem_t(iexp, nochem, :dimsizes(dummy3) - 1) = tofloat(dummy3)
                nochem_trend(iexp, nochem) = tofloat(regline(ttmp, dummy3) * 10)
                nochem = nochem + 1
            end if

            delete(dummy1)
            delete(dummy2)
            delete(dummy3)
            delete(ttmp)
        end do  ; loop over models

        if (nochem.gt.0) then
            trend     = new((dim_EXP) * 3, "float")
            trend_std = new((dim_EXP) * 3, "float")
            ti = ispan(1, (dim_EXP) * 3, 1)
        else
            trend     = new((dim_EXP), "float")
            trend_std = new((dim_EXP), "float")
            ti = ispan(1, (dim_EXP), 1)
        end if

        ;-------------------------------------------------------

        do iexp = 0, dim_EXP - 1
            chem_t(iexp, dimsizes(list_chem_mod), :) = \
                tofloat(dim_avg_n(chem_t(iexp, \
                :dimsizes(list_chem_mod) - 1, :), 0))
            chem_trend(iexp, dimsizes(list_chem_mod)) = \
                tofloat(regline(time(iexp, :), chem_t(iexp, \
                dimsizes(list_chem_mod), :)) * 10)

            if (nochem.gt.0) then
                nochem_t(iexp, dimsizes(list_nochem_mod), :) = \
                    tofloat(dim_avg_n(nochem_t(iexp, \
                    :dimsizes(list_nochem_mod) - 1, :), 0))
                nochem_trend(iexp, dimsizes(list_nochem_mod)) = \
                    tofloat(regline(time(iexp, :), chem_t(iexp, \
                    dimsizes(list_nochem_mod), :)) * 10) 
            end if 
        end do  ; loop over experiments

        do iexp = 0, dim_EXP - 1
            if (nochem.gt.0) then
                a = array_append_record(chem_trend(iexp, \
                    0:dimsizes(chem_trend(iexp, :)) - 2), \
                    nochem_trend(iexp, 0:dimsizes(nochem_trend(iexp, :)) - 2), 0)
                trend(iexp * 3) = avg(a)
                trend(iexp * 3 + 1) = dim_avg_n(chem_trend(iexp, \
                    0:dimsizes(chem_trend(iexp, :)) - 2), 0)
                trend(iexp * 3 + 2) = dim_avg_n(chem_trend(iexp, \
                    0:dimsizes(nochem_trend(iexp, :)) - 2), 0) 
                trend_std(iexp * 3) = stddev(a)
                trend_std(iexp * 3 + 1) = dim_stddev_n(chem_trend(iexp, \
                    0:dimsizes(chem_trend(iexp, :)) - 2), 0)
                trend_std(iexp * 3 + 2) = dim_stddev_n(chem_trend(iexp, \
                    0:dimsizes(nochem_trend(iexp, :)) - 2), 0)
            else
                a = chem_trend(iexp, 0:dimsizes(chem_trend(iexp, :)) - 2)
                trend(iexp) = avg(a)
                trend_std(iexp) = stddev(a)
            end if

            delete(a)
        end do

        delete(time)

        ;--------PLOT---------------------------------------

        if ((latrange(ipan, 0).eq.-82.5).and.(latrange(ipan, 1).eq.82.5)) then
            latstring = "NG"
        else
            if (latrange(ipan, 0).ge.0.) then
                latstring = "(" + latrange(ipan, 0) + "N-"
            else if (latrange(ipan, 0).lt.0.) then
                latstring = "(" + -1 * latrange(ipan, 0) + "S-"
            end if
            end if

            if (latrange(ipan, 1).ge.0.) then
                latstring = latstring + latrange(ipan, 1) + "N)"
            else if (latrange(ipan, 1).lt.0.) then
                latstring = latstring + -1 * latrange(ipan, 1) + "S)"
            end if
            end if
        end if

        res                   = True                      ;  plot mods desired

        if ((isatt(variable_info, "units"))) then
            res@tiYAxisString = "(" + variable_info@units + "/dec)"  ; add units title
        else
            res@tiYAxisString = "(per dec)"
        end if

        if (season.eq."yearly") then
            season1 = "ANN"
        else
            season1 = season
        end if

        res@tiMainString      = "(" + plot_num + ") " + lev + " hPa "\
            + season + " " + latstring + " " + var0 + " trend "  ; title
        res@gsnFrame          = False       ; don't draw yet
        res@gsnDraw           = False       ; don't advance frame
        res@xyMarkLineMode    = "Markers"   ; choose to use markers
        res@xyMarkers         =  16         ; choose type of marker  
        res@xyMarkerColor     = "Black"     ; Marker color
        res@xyMarkerSizeF     = 0.015       ; Marker size (default 0.01)
        res@tiMainFontHeightF = 0.024
        res@trYMinF = 2 * min(trend)  ; set Y-axis min
        res@trYMaxF = 2 * max(trend)  ; set Y-axis max
        res@trXMinF = 0               ; set X-axis min
        res@trXMaxF = dimsizes(trend) + 1
        res@vpWidthF = 0.8            ; Change the aspect ratio, but 
        res@vpHeightF = 0.4           ; make plot as large as possible.
        res@gsnYRefLine = 0.0
        res@gsnYRefLineColor = "Grey"
        res@tmXBMinorOn = False
        res@tmXTMinorOn = False
        res@tmXBLabelsOn = False

        plot(ipan) = gsn_csm_xy(wks, (/ti(0), ti(0)/), (/trend(0), trend(0)/), res)

        mkres = True
        mkres@gsMarkerSizeF = 0.015
        mkres@gsMarkerThicknessF = 3.
        color_m = (/"Black", "Blue", "lightskyblue", "orange", "red", \
            "Grey", "green"/)

        do iexp = 0, dim_EXP - 1
            mkres@gsMarkerColor = color_m(iexp)
            mkres@gsMarkerIndex = 16
            if (nochem.gt.0) then
                dum(ipan, iexp * 3) = gsn_add_polymarker(wks, plot(ipan), \
                    ti(iexp * 3), trend(iexp * 3), mkres)
                if (nochem.gt.0) then
                    mkres@gsMarkerIndex = 6
                    dum(ipan, iexp * 3 + 1) = gsn_add_polymarker(wks, \
                        plot(ipan), ti(iexp * 3 + 1), trend(iexp * 3 + 1), mkres)
                    mkres@gsMarkerIndex = 7
                    dum(ipan, iexp * 3 + 2) = gsn_add_polymarker(wks, \
                        plot(ipan), ti(iexp * 3 + 2), trend(iexp * 3 + 2), mkres)
                end if
            else
                dum(ipan, iexp) = gsn_add_polymarker(wks, plot(ipan), \
                    ti(iexp), trend(iexp), mkres)
            end if
        end do

        ;---------add error bars----------------------------

        polyres   = True  ; marker resources
        color_lin = (/"black", "black", "black", "blue", "blue", "blue", \
            "lightskyblue", "lightskyblue", "lightskyblue", "orange", \
            "orange", "orange", "red", "red", "red", "green"/)

        do t = 0, dimsizes(trend) - 1
            if (nochem.gt.0) then
                polyres@gsLineColor = color_lin(t)
            else
                polyres@gsLineColor = color_m(t)
            end if   

            error_bar(ipan, t) = gsn_add_polyline(wks, plot(ipan), \
                (/ti(t), ti(t)/), (/trend(t) - trend_std(t), \
                trend(t) + trend_std(t)/), polyres)
        end do

        delete(nochem_trend)
        delete(nochem_t)
        delete(chem_trend)
        delete(chem_t)
        delete(trend)
        delete(trend_std)
    end do  ; ipan

    ;-------------labels------------------------------

    if (nochem.gt.0) then
        lbres                    = True          ; labelbar only resources
        lbres@vpWidthF           = 0.1           ; labelbar width
        lbres@vpHeightF          = 0.1           ; labelbar height
        lbres@lbBoxMajorExtentF  = 0.15          ; puts space between color boxes
        lbres@lbMonoFillPattern  = True          ; Solid fill pattern
        lbres@lbLabelFontHeightF = 0.015         ; font height. default is small
        lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
        lbres@lbPerimOn          = False
        xpos = (/0.05, 0.25, 0.45/)
        index_p = (/16, 7, 6/)
        labels_1 = (/"All", "nochem", "chem"/)

        txres               = True         ; Text resources
        txres@txFontHeightF = 0.015
        txres@txJust        = "CenterLeft"
        mkres               = True         ; Marker resources

        do i = 0, 2
            mkres@gsMarkerThicknessF = 3.5
            mkres@gsMarkerSizeF      = 0.02
            mkres@gsMarkerIndex      = index_p(i)
            mkres@gsMarkerColor      = "Black"
            gsn_polymarker_ndc(wks, xpos(i), 0.2, mkres)
            gsn_text_ndc(wks, labels_1(i), xpos(i) + 0.02, 0.2, txres)
        end do
    end if

    xpos1 = (/0.01, 0.16, 0.26, 0.36, 0.46, 0.56, 0.66, 0.76/)
    color_p = (/"Black", "Blue", "lightskyblue", "Orange", "Red", "Green", "Grey"/)
    labels_2 = label_exp

    lbres                    = True          ; labelbar only resources
    lbres@vpWidthF           = 0.1           ; labelbar width
    lbres@vpHeightF          = 0.1           ; labelbar height
    lbres@lbBoxMajorExtentF  = 0.15          ; puts space between color boxes
    lbres@lbMonoFillPattern  = True          ; Solid fill pattern
    lbres@lbLabelFontHeightF = 0.015         ; font height. default is small
    lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
    lbres@lbPerimOn          = False

    do i = 0, dim_EXP - 1
        lbres@lbFillColors  = color_p(i)
        lbres@gsMarkerIndex = 16
        gsn_labelbar_ndc(wks, 1, labels_2(i), xpos1(i), 0.1, lbres)
    end do

    ;; common resources for panel plots 
    pan = True 
    pan@gsnMaximize = True
    pan@gsnFrame    = False
    pan@gsnPaperOrientation = "portrait" 
    pan@gsnPanelXWhiteSpacePercent = 5
    pan@gsnPanelYWhiteSpacePercent = 5
    gsn_panel(wks, plot, (/dimsizes(season), 1/), pan)  ; panel plot               
    frame(wks)

    ;--------------------------------------------------

end
