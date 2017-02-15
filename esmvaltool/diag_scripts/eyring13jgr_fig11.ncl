;##############################################################################
;# eyring13jgr_fig11.ncl                                                      #
;##############################################################################
;# Author: Irene Cionni (ENEA, Italy)                                         #
;# ESMVal project PI: Veronika Eyring (DLR, Germany)                          #
;##############################################################################
;#                                                                            #
;# Description:                                                               #
;#    * Trend relationship (scatter plot) between 2 variables ("variable 0"   #
;#      and "variable 1"); models are classified by experiment name.          #
;#                                                                            #
;# Required diag_script_info attributes (diagnostics specific):               #
;#      styleset: style set ("DEFAULT", "CMIP5")                              #
;#      e13fig11_V0_units: unit label for "variable 0" (x-axis) (string)      #
;#      e13fig11_V1_units: unit label for "variable 1" (y-axis) (string)      #
;#      e13fig11_V0_latrange: min. and max. latitude of the region for        #
;#          "variable 0"                                                      #
;#      e13fig11_V1_latrange: min. and max. latitude of the region for        #
;#          "variable 1"                                                      #
;#      e13fig11_V0_season: season for "variable 0" (e.g., "yearly")          #
;#      e13fig11_V1_season: season for "variable 1" (e.g., "yearly")          #
;#      e13fig10_V0_lev: vertical level (in hPa) for "variable 0"             #
;#      e13fig10_V1_lev: vertical level (in hPa) for "variable 1"             #
;#      plot_number: string used for plot labeling / sub-figure               #
;#          (e.g., "(a)")                                                     #
;#      e13fig11_XMin: min. x-value for plotting                              #
;#      e13fig11_XMax: max. x-value for plotting                              #
;#      e13fig11_YMin: min. y-value for plotting                              #
;#      e13fig11_YMax: max. y-value for plotting                              #
;#                                                                            #
;# Optional diag_script_info attributes (diagnostics specific):               #
;#      list_chem_mod: list of models in the group “chem” (array of strings,  #
;#          default = empty)                                                  #
;#      list_nochem_mod: list of models in the group “nochem” (array of       #
;#          strings, default = empty)                                         #
;#                                                                            #
;# Caveats:                                                                   #
;#                                                                            #
;# Modification history:                                                      #
;#    * 20151210-A-laue_ax: code clean-up, added call to write_references     #
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

;************************************************

begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output("plot 11 - " + diag_script + ".ncl (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    var1 = variables(1)
    field_type0 = field_types(0)
    field_type1 = field_types(1)
    colors=project_style(diag_script_info, "colors")
    dashes=project_style(diag_script_info, "dashes")
    thicks=tofloat(project_style(diag_script_info, "thicks")) 
    
    do idx_mod = 0, dim_MOD - 1
        data0 = read_data(idx_mod, var0, field_type0)
    end do
    do idx_mod = 0, dim_MOD - 1
        data1 = read_data(idx_mod, var1, field_type1)
    end do
end

begin
    var0deffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(var0deffile)
    var1deffile = "interface_data/" + var1 + "_info.tmp"
    loadscript(var1deffile)
end

begin
    ;references Information
    write_references(diag_script,      \  ; script name
                     "A_cion_ir",      \  ; authors
                     "",               \  ; contributors
                     "D_eyring13jgr",  \  ; diag_references
                     "",               \  ; obs_references
                     "P_esmval")          ; proj_references

    if ((isatt(diag_script_info, "plot_number"))) then
        plot_num = new(dimsizes(diag_script_info@e13fig11_V0_season), "string")
        plot_num(:) = diag_script_info@plot_number
    else
        plot_num = integertochar(ispan(97, 122, 1))
    end if

    plot_dir = getenv("ESMValTool_plot_dir")
    diag_script_base = basename(diag_script)
    output_dir = get_output_dir(plot_dir, diag_script_base)
    aux_title_info = plot_num(0)
    output_filename = interface_get_figure_filename(diag_script_base, \
        var0, var1, aux_title_info, -1)
    file_type = getenv("ESMValTool_output_file_type")

    if (ismissing(file_type)) then
        file_type = "ps"
    end if

    file_type = str_lower(file_type)

    color_type = "COLOR"
    output_file_path = output_dir + output_filename
    wks = gsn_open_wks (file_type, output_file_path)
    plot = new(dimsizes(diag_script_info@e13fig11_V0_season), graphic)

    ;-----------------------------------------------------------------

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

    ;----------------------------------------------------------------------

    if ((isatt(diag_script_info, "list_chem_mod"))) then
        list_chem_mod = diag_script_info@list_chem_mod
    else
        list_chem_mod = models@name 
    end if

    if ((isatt(diag_script_info, "list_nochem_mod"))) then
        list_nochem_mod = diag_script_info@list_nochem_mod
    else
        list_nochem_mod = "none"
    end if

    chem_t = new((/dim_EXP, 2, dimsizes(list_chem_mod) + 1, yearsize/), "float")
    chem_trend = new((/dim_EXP, 2, dimsizes(list_chem_mod) + 1/), "float")
    nochem_t = new((/dim_EXP, 2, dimsizes(list_nochem_mod) + 1, yearsize/), "float")
    nochem_trend = new((/dim_EXP, 2, dimsizes(list_nochem_mod) + 1/), "float")

    if (isatt(diag_script_info, "e13fig11_V0_season").and. \
        isatt(diag_script_info, "e13fig11_V1_season")) then
        season_V0 = diag_script_info@e13fig11_V0_season
        season_V1 = diag_script_info@e13fig11_V1_season
    else
        if (isStrSubset(field_type0, "z") .eq. False) then
            print("fatal:plot requires diag_script_info attribute 'e13fig11_season' and 'e13fig11_V1_season'")
            status_exit(1)
        end if
    end if

    if (isvar("latrange_V0")) then
        delete(latrange_V0)
    end if

    if (isatt(diag_script_info, "e13fig11_V0_latrange")) then
        latrange_V0 = new((/dimsizes(season_V0), 2/), "float")    
        latrange_V0 = diag_script_info@e13fig11_V0_latrange
    else
        print("fatal:plot requires diag_script_info attribute 'e13fig11_V0_latrange'")
        status_exit(1)
    end if

    if (isvar("latrange_V1")) then
        delete(latrange_V1)
    end if

    if (isatt(diag_script_info, "e13fig11_V1_latrange")) then
        latrange_V1 = new((/dimsizes(season_V1), 2/), "float") 
        latrange_V1 = diag_script_info@e13fig11_V1_latrange
    else
        print("fatal:plot requires diag_script_info attribute 'e13fig11_V1_latrange'")
        status_exit(1)
    end if

    if (isatt(diag_script_info, "e13fig11_V0_lev")) then
        V0_lev = toint(diag_script_info@e13fig11_V0_lev)
    end if

    if (isatt(diag_script_info, "e13fig11_V1_lev")) then
        V1_lev = toint(diag_script_info@e13fig11_V1_lev)
    end if

    ;-----------------------------------------------------------------------------
    do ipan = 0, dimsizes(season_V0) - 1  ; panel loop
        ;----------------------------chem------------------------------------
        iexp = 0
        nochem = 0
        chem = 0

        timeh = ispan(1979, 2005, 1)
        timef = ispan(2006, 2050, 1)
        time  = new((/dim_EXP, yearsize/), float)

        do imod = 0, dim_MOD - 1
            A0 = extract_data(imod, data0, -1, 0, 0)
            ; convert ozone (tro3, vmro3) to ppmv
            if (var0.eq."tro3") then
                A0=A0 * 10^-3
            end if
            if (var0.eq."vmro3") then
                A0 = A0 * 10^6
            end if
            A1 = extract_data(imod, data1, -1, 0, 0)
            if (var1.eq."tro3") then
                A1 = A1 * 10^-3
            end if 
            if (var1.eq."vmro3") then
                A0 = A0 * 10^6
            end if

            if (max(A0&plev).ge.100000) then
                A0&plev = A0&plev / 100
            end if

            if (max(A1&plev).ge.100000) then
                A1&plev = A1&plev / 100
            end if

            ttmp = ispan(toint(models@start_year(imod)), toint(models@end_year(imod)), 1)

            if (var0.ne."ua") then
                dummy1_ta = area_operations(A0(time|:, {plev|V0_lev}, lat|:, lon|:), \
                    latrange_V0(ipan, 0), latrange_V0(ipan, 1), 0., 360., "average", True)
            else
                tmp_pos = dim_avg_n_Wrap(A0(time|:, {plev|V0_lev}, \
                    {lat|latrange_V0(ipan, 0):latrange_V0(ipan, 1)}, lon|:), 2)
                dummy1_ta = new(dimsizes(A0&time), typeof(A0&lat))
                do tt = 0, dimsizes(A0&time) - 1
                    dummy1_ta(tt) = tmp_pos&lat(maxind(tmp_pos(tt, :)))
                end do 
                dummy1_ta!0 = "time"
                dummy1_ta&time = A0&time
                delete(tmp_pos)
            end if

            if (var1.ne."ua") then      
                dummy1_O3 = area_operations(A1(time|:, {plev|V1_lev}, lat|:,lon|:), \
                    latrange_V1(ipan, 0), latrange_V1(ipan, 1), 0., 360., "average", True)
                dummy1_O3 = dummy1_O3
            else
                tmp_pos = dim_avg_n_Wrap(A0(time|:, {plev|V1_lev}, \
                    {lat|latrange_V1(ipan, 0):latrange_V1(ipan, 1)}, lon|:), 2)
                dummy1_O3 = new(dimsizes(A0&time), typeof(A0&lat))
                do tt = 0, dimsizes(A0&time) - 1
                    dummy1_O3(tt) = tmp_pos&lat(maxind(tmp_pos(tt, :)))
                end do 
                dummy1_O3!0 = "time"
                dummy1_O3&time = A0&time
 
                delete(tmp_pos)
            end if

            delete(A0)
            delete(A1)

            if (season_V0(ipan).ne."yearly") then
                tmp_ta = time_operations(dummy1_ta, toint(models@start_year(imod)), \
                    toint(models@end_year(imod)), "extract", season_V0(ipan), True)
                tmp_O3 = time_operations(dummy1_O3, toint(models@start_year(imod)), \
                    toint(models@end_year(imod)), "extract", season_V1(ipan), True)
                dummy2_ta = time_operations(tmp_ta, toint(models@start_year(imod)), \
                    toint(models@end_year(imod)), "average", "yearly", True)
                dummy2_O3 = time_operations(tmp_O3, toint(models@start_year(imod)), \
                    toint(models@end_year(imod)), "average", "yearly", True)
                delete(tmp_ta)
                delete(tmp_O3)
            else
                dummy2_ta = time_operations(dummy1_ta, toint(models@start_year(imod)), \
                    toint(models@end_year(imod)), "average", season_V0(ipan), True)
                dummy2_O3 = time_operations(dummy1_O3, toint(models@start_year(imod)), \
                    toint(models@end_year(imod)), "average", season_V1(ipan), True)
            end if

            dummy3_ta = dummy2_ta-dim_avg_n(dummy2_ta, 0)
            dummy3_O3 = dummy2_O3-dim_avg_n(dummy2_O3, 0)
        
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
                chem_t(iexp, 0, chem,:dimsizes(dummy3_ta) - 1) = tofloat(dummy3_ta)
                chem_t(iexp, 1, chem,:dimsizes(dummy3_O3) - 1) = tofloat(dummy3_O3)
                chem_trend(iexp, 0, chem) = tofloat(regline(ttmp, dummy3_ta) * 10)
                chem_trend(iexp, 1, chem) = tofloat(regline(ttmp, dummy3_O3) * 10)
                chem = chem + 1
            end if

            if (any(models@name(imod).eq.list_nochem_mod)) then
                nochem_t(iexp, 0, nochem, :dimsizes(dummy3_ta) - 1) = tofloat(dummy3_ta)
                nochem_t(iexp, 1, nochem, :dimsizes(dummy3_O3) - 1) = tofloat(dummy3_O3)
                nochem_trend(iexp, 0, nochem) = tofloat(regline(ttmp, dummy3_ta) * 10) 
                nochem_trend(iexp, 1, nochem) = tofloat(regline(ttmp, dummy3_O3) * 10)
                nochem = nochem + 1
            end if

            delete(dummy1_ta)
            delete(dummy2_ta)
            delete(dummy3_ta)
        
            delete(dummy1_O3)
            delete(dummy2_O3)
            delete(dummy3_O3)
            delete(ttmp)

        end do  ; loop over seasons (var0)

        do iexp = 0, dim_EXP - 1
            chem_t(iexp, 0, dimsizes(list_chem_mod), :) = \
                tofloat(dim_avg_n(chem_t(iexp, 0, :dimsizes(list_chem_mod) - 1, :), 0))
            chem_t(iexp, 1, dimsizes(list_chem_mod), :) = \
                tofloat(dim_avg_n(chem_t(iexp, 1, :dimsizes(list_chem_mod) - 1, :), 0))
            chem_trend(iexp, 0, dimsizes(list_chem_mod)) = \
                tofloat(regline(time(iexp, :), chem_t(iexp, 0, dimsizes(list_chem_mod), :)) * 10) 
            chem_trend(iexp, 1, dimsizes(list_chem_mod)) = \
                tofloat(regline(time(iexp, :), chem_t(iexp, 1, dimsizes(list_chem_mod), :)) * 10)

            nochem_t(iexp, 0, dimsizes(list_nochem_mod), :) = \
                tofloat(dim_avg_n(nochem_t(iexp, 0, :dimsizes(list_nochem_mod) - 1, :), 0))
            nochem_t(iexp, 1, dimsizes(list_nochem_mod), :) = \
                tofloat(dim_avg_n(nochem_t(iexp, 1, :dimsizes(list_nochem_mod) - 1, :), 0))
            nochem_trend(iexp,0,dimsizes(list_nochem_mod)) = \
                tofloat(regline(time(iexp, :), chem_t(iexp, 0, dimsizes(list_nochem_mod), :)) * 10)
            nochem_trend(iexp, 1, dimsizes(list_nochem_mod)) = \
                tofloat(regline(time(iexp, :), chem_t(iexp, 1, dimsizes(list_nochem_mod), :)) * 10)
        end do

        delete(timef)
        delete(timeh)
        delete(time)

        ;-------correlation----------------------------------------------------------------------

        if (any(isStrSubset(models@experiment(0), "historical").eq.True)) then
;            c_nochem_trend_h = new((/2, dim_MOD + 1/), "float")
;            c_chem_trend_h   = new((/2, dim_MOD + 1/), "float")
            c_tot_trend_h    = new((/2, (dimsizes(list_chem_mod) + \
                dimsizes(list_nochem_mod) + 2)/), "float")
            c_nochem_trend_f = new((/2, (dimsizes(list_nochem_mod) + 1) * \
                (dim_EXP - 1)/), "float")
            c_chem_trend_f   = new((/2, (dimsizes(list_chem_mod) + 1) * \
                (dim_EXP - 1)/), "float")
            c_tot_trend_f    = new((/2, (dimsizes(c_nochem_trend_f(0, :)) + \
                dimsizes(c_chem_trend_f(0, :)))/), "float")

            c_nochem_trend_h = (/nochem_trend(0, :, :)/)
            c_chem_trend_h   = (/chem_trend(0, :, :)/)
            c_tot_trend_h(0, :) = array_append_record(c_nochem_trend_h(0, :), \
                c_chem_trend_h(0, :), 0)
            c_tot_trend_h(1, :) = array_append_record(c_nochem_trend_h(1, :), \
                c_chem_trend_h(1, :), 0)

            c_nochem_trend_f(0, :dimsizes(list_nochem_mod)) = (/nochem_trend(1, 0, :)/)
            c_nochem_trend_f(1, :dimsizes(list_nochem_mod)) = (/nochem_trend(1, 1, :)/)
            c_chem_trend_f(0, :dimsizes(list_chem_mod))   = (/chem_trend(1, 0, :)/)
            c_chem_trend_f(1, :dimsizes(list_chem_mod))   = (/chem_trend(1, 1, :)/)

            do iexp = 2, dim_EXP - 1
                c_nochem_trend_f(0, dimsizes(nochem_trend(iexp - 1, 0, :)) * \
                    (iexp - 1):dimsizes(nochem_trend(iexp - 1, 0, :)) * \
                    (iexp - 1) + dimsizes(nochem_trend(iexp, 0, :)) - 1) = \
                    nochem_trend(iexp, 0, :)
              c_nochem_trend_f(1, dimsizes(nochem_trend(iexp - 1, 1, :)) * \
                  (iexp -1 ):dimsizes(nochem_trend(iexp - 1, 1, :)) * \
                  (iexp - 1) + dimsizes(nochem_trend(iexp, 1, :)) - 1) = \
                  nochem_trend(iexp, 1, :)
              c_chem_trend_f(0, dimsizes(chem_trend(iexp - 1, 0, :)) * \
                  (iexp - 1):dimsizes(chem_trend(iexp - 1, 0, :)) * \
                  (iexp - 1) + dimsizes(chem_trend(iexp, 0, :)) - 1) = \
                  chem_trend(iexp, 0, :)
              c_chem_trend_f(1, dimsizes(chem_trend(iexp - 1, 1, :)) * \
                  (iexp - 1):dimsizes(chem_trend(iexp - 1, 1, :)) * \
                  (iexp - 1) + dimsizes(chem_trend(iexp, 1, :)) - 1) = \
                  chem_trend(iexp, 1, :)
            end do
        else
            c_nochem_trend_f(0, :dimsizes(list_nochem_mod)) = (/nochem_trend(0, 0, :)/)
            c_nochem_trend_f(1, :dimsizes(list_nochem_mod)) = (/nochem_trend(0, 1, :)/)
            c_chem_trend_f(0, :dimsizes(list_chem_mod))     = (/chem_trend(0, 0, :)/)
            c_chem_trend_f(1, :dimsizes(list_chem_mod))     = (/chem_trend(0, 1, :)/)

            do iexp = 1, dim_EXP - 1
                c_nochem_trend_f(0, dimsizes(nochem_trend(iexp - 1, 0, :)) * \
                    (iexp - 1):dimsizes(nochem_trend(iexp - 1, 0, :)) * \
                    (iexp - 1) + dimsizes(nochem_trend(iexp, 0, :)) - 1) = \
                    nochem_trend(iexp, 0, :)
                c_nochem_trend_f(1, dimsizes(nochem_trend(iexp - 1, 1, :)) * \
                    (iexp - 1):dimsizes(nochem_trend(iexp - 1, 1, :)) * \
                    (iexp - 1) + dimsizes(nochem_trend(iexp, 1, :)) - 1) = \
                    nochem_trend(iexp, 1, :)
                c_chem_trend_f(0, dimsizes(chem_trend(iexp - 1, 0, :)) * \
                    (iexp - 1):dimsizes(chem_trend(iexp - 1, 0, :)) * \
                    (iexp - 1) + dimsizes(chem_trend(iexp, 0, :)) - 1) = \
                    chem_trend(iexp, 0, :)
                c_chem_trend_f(1, dimsizes(chem_trend(iexp - 1, 1, :)) * \
                    (iexp - 1):dimsizes(chem_trend(iexp - 1, 1, :)) * \
                    (iexp - 1) + dimsizes(chem_trend(iexp, 1, :)) - 1) = \
                    chem_trend(iexp, 1, :)
            end do
        end if

        c_tot_trend_f(0, :) = array_append_record(c_nochem_trend_f(0, :), c_chem_trend_f(0, :), 0)
        c_tot_trend_f(1, :) = array_append_record(c_nochem_trend_f(1, :), c_chem_trend_f(1, :), 0)

        ;____________ccc______________________________________________________

        rc_nochem_trend_h = escorc(c_nochem_trend_h(0, :), c_nochem_trend_h(1, :))
        rc_chem_trend_h = escorc(c_chem_trend_h(0, :), c_chem_trend_h(1, :))
        rc_tot_trend_h = escorc(c_tot_trend_h(0, :), c_tot_trend_h(1, :))

        rc_nochem_trend_f = escorc(c_nochem_trend_f(0, :), c_nochem_trend_f(1, :))
        rc_chem_trend_f = escorc(c_chem_trend_f(0, :), c_chem_trend_f(1, :))
        rc_tot_trend_f = escorc(c_tot_trend_f(0, :), c_tot_trend_f(1, :))

        ;____________ccc______________________________________________________

;        t_trend_nochem = rc_trend_nochem * ((dimsizes(c_trend_nochem(0, :)) - 2) \
;            / (1 - rc_trend_nochem^2))^0.5
;        t_trend_chem = rc_trend_chem * ((dimsizes(c_trend_chem(0, :)) - 2) \
;            / (1 - rc_trend_chem^2))^0.5
;        t_trend_tot = rc_trend_tot * ((dimsizes(c_trend_tot(0, :)) - 2) \
;            / (1 - rc_trend_tot^2))^0.5

;        t_trend_nochem_h = rc_trend_nochem_h * ((dimsizes(c_trend_nochem_h(0, :)) - 2) \
;            / (1 - rc_trend_nochem_h^2))^0.5
;        t_trend_chem_h = rc_trend_chem_h * ((dimsizes(c_trend_chem_h(0, :)) - 2) \
;            / (1 - rc_trend_chem_h^2))^0.5
;        t_trend_tot_h = rc_trend_tot_h * ((dimsizes(c_trend_tot_h(0, :)) - 2) \
;            / (1 - rc_trend_tot_h^2))^0.5

;        t_trend_nochem_f = rc_trend_nochem_f * ((dimsizes(c_trend_nochem_f(0, :)) - 2) \
;            / (1 - rc_trend_nochem_f^2))^0.5
;        t_trend_chem_f = rc_trend_chem_f * ((dimsizes(c_trend_chem_f(0, :)) - 2) \
;            / (1 - rc_trend_chem_f^2))^0.5
;        t_trend_tot_f=rc_trend_tot_f * ((dimsizes(c_trend_tot_f(0, :)) - 2) \
;            / (1 - rc_trend_tot_f^2))^0.5

;        rt_trend_nochem = rtest(rc_trend_nochem, dimsizes(c_trend_nochem(1, :)), 0)
;        rt_trend_chem = rtest(rc_trend_chem, dimsizes(c_trend_chem(1, :)), 0)
;        rt_trend_tot = rtest(rc_trend_tot, dimsizes(c_trend_tot(1, :)), 0)
      
;        rt_trend_nochem_h = rtest(rc_trend_nochem_h, dimsizes(c_trend_nochem_h(1, :)), 0)
;        rt_trend_chem_h = rtest(rc_trend_chem_h, dimsizes(c_trend_chem_h(1, :)), 0)
;        rt_trend_tot_h = rtest(rc_trend_tot_h, dimsizes(c_trend_tot_h(1, :)), 0)

;        rt_trend_nochem_f = rtest(rc_trend_nochem_f, dimsizes(c_trend_nochem_f(1, :)), 0)
;        rt_trend_chem_f = rtest(rc_trend_chem_f, dimsizes(c_trend_chem_f(1, :)), 0)
;        rt_trend_tot_f = rtest(rc_trend_tot_f, dimsizes(c_trend_tot_f(1, :)), 0)

        ; plot 
        ;-----------------------------------------------------------------------------
        ; scatter
        p_chem   = new((/dim_EXP, dimsizes(list_chem_mod)/), "graphic")
        p_nochem = new((/dim_EXP, dimsizes(list_nochem_mod)/), "graphic")
        ;-----------------------------------------------------------------------------

        res                    = True
        res@gsnDraw            = False       ; Don't draw plot
        res@gsnFrame           = False       ; Don't advance frame

        if (isatt(diag_script_info, "e13fig11_Xmin")) then
            res@trXMinF = diag_script_info@e13fig11_Xmin
        end if

        if (isatt(diag_script_info, "e13fig11_Xmax")) then
            res@trXMaxF = diag_script_info@e13fig11_Xmax
        end if

        if (isatt(diag_script_info, "e13fig11_Ymin")) then
            res@trYMinF = diag_script_info@e13fig11_Ymin
        end if

        if (isatt(diag_script_info, "e13fig11_Ymax")) then
            res@trYMaxF = diag_script_info@e13fig11_Ymax
        end if

;        res@trYMinF            = -0.1              ; min value on y-axis
;        res@trYMaxF            = 0.1               ; max value on y-axis
;        res@trXMinF            = -0.5              ; min value on x-axis
;        res@trXMaxF            = 0.5               ; max value on x-axis
        res@xyMarker           = 6 
        res@xyMarkLineMode     = "Markers"
        res@xyMarkerSizeF      = 0.01
        res@xyMarkerThicknessF = 3.
        res@xyMarkerColor      = "White"
        res@tiMainString       = "("+plot_num(ipan)+")"
        res@tiMainPosition     = "Left"
        res@gsnYRefLine        = 0.0                        
        res@gsnXRefLine        = 0.0

        if (latrange_V0(ipan, 0).eq.-82.5.and.latrange_V0(ipan, 1).eq.82.5) then
            latstring1 = "NG"
        else
            if (latrange_V0(ipan, 0).ge.0.) then
                latstring1 = "(" + latrange_V0(ipan, 0) + "N"
            else if (latrange_V0(ipan, 0).lt.0.) then
                latstring1 = "(" + (-latrange_V0(ipan, 0)) + "S"
            end if
            end if
            if (latrange_V0(ipan, 1).ge.0.) then
                latstring1 = latstring1 + latrange_V0(ipan, 1) + "N)"
            else if (latrange_V0(ipan, 1).lt.0.) then
                latstring1 = latstring1 + (-latrange_V0(ipan, 1)) + "S)"
            end if
            end if
        end if

        if (latrange_V1(ipan, 0).eq.-82.5.and.latrange_V1(ipan, 1).eq.82.5) then
            latstring2 = "NG"
        else
            if (latrange_V1(ipan, 0).ge.0.) then
                latstring2 = "(" + latrange_V1(ipan, 0) + "N"
            else if (latrange_V1(ipan, 0).lt.0.) then
                latstring2 = "(" + (-latrange_V1(ipan, 0)) + "S"
            end if
            end if
            if (latrange_V1(ipan, 1).ge.0.) then
                latstring2 = latstring2+latrange_V1(ipan, 1)+"N)"
            else if (latrange_V1(ipan, 1).lt.0.) then
                latstring2 = latstring2 + (-latrange_V1(ipan, 1)) + "S)"
            end if
            end if
        end if

        if (season_V0.eq."yearly") then
            season_l_V0 = "ANN"
        else
            season_l_V0 = season_V0
        end if

        if (season_V1.eq."yearly") then
            season_l_V1 = "ANN"
        else
            season_l_V1 = season_V1
        end if

        res@tiXAxisString     = V0_lev + " hPa " + season_l_V0 + " " + latstring1 + \
            " " + var0 + " trend (" + diag_script_info@e13fig11_V0_units + "/dec)"
        res@tiYAxisString     = V1_lev + " hPa " + season_l_V1 + " " + latstring2 + \
            " " + var1 + " trend (" + diag_script_info@e13fig11_V1_units+"/dec) "
        res@tiMainFontHeightF = 0.03

        plot(ipan) = gsn_csm_xy(wks, chem_trend(:, 0, :), chem_trend(:, 1, :), res) 

        color_p = (/"Black", "Blue", "lightskyblue", "Orange", "Red", "Grey", "Green"/)

        do iexp = 0, dim_EXP - 1
            print(iexp) 
            print(chem_trend(iexp, 0, :))
            print(chem_trend(iexp, 1, :))
            do ichem = 0, dimsizes(list_chem_mod) - 1
                ;----------add chem trend---------------------------------
                mkres = True
                mkres@gsMarkerIndex = 6
                mkres@gsMarkerColor = color_p(iexp)

                p_chem(iexp, ichem) = gsn_add_polymarker(wks, plot(ipan), \
                    chem_trend(iexp, 0, ichem), chem_trend(iexp, 1, ichem), mkres)
            end do

            ;----------add nochem trend-------------------------------

            do inochem = 0, dimsizes(list_nochem_mod) - 1
                mkres@gsMarkerIndex = 7
                mkres@gsMarkerColor = color_p(iexp)
                p_nochem(iexp, inochem) = gsn_add_polymarker(wks, plot(ipan), \
                    nochem_trend(iexp, 0, inochem), nochem_trend(iexp, 1, inochem), mkres)
            end do
        end do

        txres               = True         ; Text resources
        txres@txFontHeightF = 0.015
        txres@txJust        = "CenterLeft"
        dummy1 = gsn_add_text(wks, plot(ipan), "CMIP5 all =" + \
            sprintf("%3.2f", rc_tot_trend_h) + "*, " + sprintf("%3.2f", rc_tot_trend_f), \
            -1.8, 0.9, txres)

        print("-----------------------------------------------------")
        print("chem CMIP5 historical=" + sprintf("%3.2f", rc_chem_trend_h))
        print("nochem CMIP5 historical=" + sprintf("%3.2f", rc_nochem_trend_h))
        print("CMIP5 all historical=" + sprintf("%3.2f", rc_tot_trend_h))
        print("chem CMIP5 rcp=" + sprintf("%3.2f", rc_chem_trend_f))
        print("nochem CMIP5 rcp=" + sprintf("%3.2f", rc_nochem_trend_f))
        print("CMIP5 all rcp=" + sprintf("%3.2f", rc_tot_trend_f))
        print("-----------------------------------------------------")

    end do  ; loop over seasons
 
    ;-------------labels------------------------------

;    xpos     = (/0.05, 0.25, 0.45/)
;    index_p  = (/6, 4, 7/)
;    labels_1 = (/"nochem CMIP5", "chem CMIP5","CCMVal-2"/)
;    mkres    = True         ; Marker resources
;    do i = 0, 1
;        mkres@gsMarkerThicknessF = 3.5
;        mkres@gsMarkerSizeF      = 0.02
;        mkres@gsMarkerIndex      = index_p(i)
;        mkres@gsMarkerColor      = "Black"
;        gsn_polymarker_ndc(wks, xpos(i), 0.07, mkres)
;        gsn_text_ndc(wks,labels_1(i), xpos(i) + 0.02, 0.07, txres)
;    end do

    lbres                    = True          ; labelbar only resources
    lbres@vpWidthF           = 0.1           ; labelbar width
    lbres@vpHeightF          = 0.1           ; labelbar height
    lbres@lbBoxMajorExtentF  = 0.15          ; puts space between color boxes
    lbres@lbMonoFillPattern  = True          ; Solid fill pattern
    lbres@lbLabelFontHeightF = 0.015         ; font height. default is small
    lbres@lbLabelJust        = "CenterLeft"  ; left justify labels
    lbres@lbPerimOn          = False
    xpos = (/0.05, 0.25, 0.45, 0.65/)
    index_p = (/6, 4, 7/)
    labels_1 = (/"CMIP5 all", "nochem CMIP5", "chem CMIP5", "CCMVal-2"/)

    mkres               = True         ; Marker resources

    do i = 0, 1
        mkres@gsMarkerThicknessF = 3.5
        mkres@gsMarkerSizeF      = 0.02
        mkres@gsMarkerIndex      = index_p(i)
        mkres@gsMarkerColor      = "Black"
        gsn_polymarker_ndc(wks, xpos(i), 0.07, mkres)
        gsn_text_ndc      (wks, labels_1(i), xpos(i) + 0.02, 0.07, txres)
    end do

    xpos1 = (/0.01, 0.16, 0.26, 0.36, 0.46, 0.56, 0.66, 0.76/)
    color_p = (/"Black", "Blue", "lightskyblue", "Orange", "Red", "Green", "Grey"/)
    labels_2 = label_exp

    do i = 0, dim_EXP - 1
        lbres@lbFillColors   = color_p(i)
         lbres@gsMarkerIndex = 16
         gsn_labelbar_ndc(wks, 1, labels_2(i), xpos1(i), 0.07, lbres)
    end do

    ;; common resources for panel plots 
    pan = True 
    pan@gsnMaximize = True
    pan@gsnFrame    = False
    pan@gsnPaperOrientation = "portrait" 
    pan@gsnPanelXWhiteSpacePercent = 5
    pan@gsnPanelYWhiteSpacePercent = 10
    gsn_panel(wks, plot, (/3, 2/), pan)   ; panel plot               
    frame(wks)

    ;--------------------------------------------------

end
