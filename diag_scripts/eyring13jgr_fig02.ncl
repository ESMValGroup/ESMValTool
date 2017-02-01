;##############################################################################
;# TIME SERIES PLOT OF SELECTED SEASONS and  LATITUDES BANDs                  #
;##############################################################################
;# Author: Irene Cionni (ENEA, ITALY)                                         #
;# ESMVal project                                                             #
;##############################################################################
;#                                                                            #
;# Description:                                                               #
;#    Calculates time series of area-weighted total ozone from 1960-2005      #
;#    averaged over given domains and seasons.                                #
;#                                                                            #
;# Required diag_script_info attributes (diagnostics specific):               #
;#    * e13fig02_latrange: extremes of latitude band                          #
;#    * e13fig02_season: ANN=Annual, DJF=December, January,February, ect      #
;#    * e13fig02_lev: selction of level if original fiels is T3*, is          #
;#      not required if variables already in T2*s fiels                       #
;#                                                                            #
;# Optional diag_script_info attributes (diagnostics specific):               #
;#    * list_chem_mod: names of models to apply a MMM, "All" if not need      #
;#    * list_nochem_mod: names of models to apply a MMM, "None" if not need   #
;#    * e13fig02_legend: switch on or off the full legend.Opt "True","False"  # 
;#    * e13fig02_legend_MMM: switch on or off the legend of MMM groups        #
;#                                                                            #
;# Required variable_info attributes (variable specific):                     #
;#                                                                            #
;# Caveats:                                                                   #
;#                                                                            #
;# Modification history:                                                      #
;#    20151203-A-laue_ax: code clean-up, moved call to write_references to    #
;#                        beginning of code                                   #
;#    201510??-A-cion_ir: written                                             #
;#                                                                            #
;##############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

load "./diag_scripts/lib/ncl/misc_function.ncl"

load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"
load "./diag_scripts/lib/ncl/set_operators.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/style.ncl"

;************************************************
;************************************************

begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(" " + diag_script + ".ncl (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    dim_MOD1 = dimsizes(models@name)
    var0 = variables(0)
    field_type0 = field_types(0)

    colors = project_style(diag_script_info, "colors")
    dashes = project_style(diag_script_info, "dashes")
    thicks = tofloat(project_style(diag_script_info, "thicks"))
    markers = project_style(diag_script_info, "markers")
    flag_mod = where(models@project.ne."OBS", 1, 0)
    index_mod = ind(flag_mod.gt.0)
    index_obs = ind(flag_mod.eq.0)
    dim_MOD = dimsizes(index_mod)
    dim_OBS = 0
    if (.not.all(ismissing(index_obs))) then
        dim_OBS = dimsizes(index_obs)
        list_obs_mod = models@name(index_obs)
    end if

    do idx_mod = 0, dim_MOD1 - 1
        data0 = read_data(idx_mod,var0, field_type0)
    end do
end

begin
    vardeffile = "interface_data/" + variables + "_info.tmp"
    loadscript(vardeffile)
end

begin
    ;; References Information
    write_references(diag_script,      \  ; script name
                     "A_cion_ir",      \  ; authors
                     "",               \  ; contributors
                     "D_eyring13jgr",  \  ; diag_references
                     "E_niwa",         \  ; obs_references
                     "P_esmval")          ; proj_references

    plot_fields = (/"T3M", "T2Ms"/)
    ind_f = inlist(field_type0, plot_fields)
    if (all(ismissing(ind_f))) then
        print("fatal:can't use field, " + field_type0 + " for diag_script " \
              + diag_script)
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)
    
    plot_dir = getenv("ESMValTool_plot_dir")
    diag_script_base = basename(diag_script)
    output_dir = get_output_dir(plot_dir, diag_script_base)
    file_type = getenv("ESMValTool_output_file_type")
    
    if(ismissing(file_type)) then
        file_type = "ps"
    end if
    file_type = str_lower(file_type)
    color_type = "COLOR"
    
     aux_title_info = ""
     output_filename = interface_get_figure_filename(diag_script_base,\
                var0, field_type0, aux_title_info, -1)
     output_file_path = output_dir + output_filename
   

    wks  = gsn_open_wks (file_type, output_file_path)
    plot = new(dimsizes(diag_script_info@e13fig02_season), "graphic")
    plot_num = integertochar(ispan(97, 122, 1))

    ;-----------------------------------------------------------------

    dim_EXP = 1  ; dimsizes of different experiments
    label_exp = models@experiment(0)
    do n = 1, dimsizes(models@experiment) - 1
        if (models@project(n).ne."OBS") then
            if (models@experiment(n).ne.models@experiment(n - 1)) then
                dim_EXP = dim_EXP + 1
                ltmp = array_append_record(label_exp, models@experiment(n), 0)
                delete(label_exp)
                label_exp = ltmp
                delete(ltmp)
            end if
        end if
    end do

    enssize = dimsizes(ind(models@experiment(0).eq.models@experiment))
    yearsize = max(toint(models@end_year)) - min(toint(models@start_year)) + 1
    dum = new((/dimsizes(diag_script_info@e13fig02_season), (dim_EXP) * 3/), graphic)
    
    if ((isatt(diag_script_info, "list_chem_mod"))) then
        list_chem_mod = diag_script_info@list_chem_mod
        attrib_chem = "_CHEM"
        if (list_chem_mod.eq."All") then
            delete(list_chem_mod)
            list_chem_mod = models@name(index_mod)
            attrib_chem = " "
        end if
    else 
        list_chem_mod = models@name 
        attrib_chem = " "
    end if

    if ((isatt(diag_script_info, "list_nochem_mod"))) then
        list_nochem_mod = diag_script_info@list_nochem_mod
    end if
    if (isvar("latrange")) then
        delete(latrange)
    end if
    if (isatt(diag_script_info, "e13fig02_latrange")) then
        latrange = diag_script_info@e13fig02_latrange
    else
        print("fatal:plot requires diag_script_info attribute 'e13fig02_latrange'")
        status_exit(1)
    end if
    if (isatt(diag_script_info, "e13fig02_season")) then
        season = diag_script_info@e13fig02_season
    else
        print("fatal:plot requires diag_script_info attribute 'e13fig02_season'")
    end if
    if (isatt(diag_script_info, "e13fig02_lev")) then
        lev = toint(diag_script_info@e13fig02_lev)
    end if
    if ((isatt(diag_script_info, "e13fig02_XMin"))) then
        XMin = diag_script_info@e13fig02_XMin
    else
        XMin= new(dimsizes(season), "float")
    end if
    if ((isatt(diag_script_info, "e13fig02_XMax"))) then
        XMax = diag_script_info@e13fig02_XMax
    else
        XMax=new(dimsizes(season), "float")
    end if
    if ((isatt(diag_script_info, "e13fig02_XMin"))) then
        YMin = diag_script_info@e13fig02_YMin
    else
        YMin=new(dimsizes(season), "float")
    end if
    if ((isatt(diag_script_info, "e13fig02_XMax"))) then
        YMax = diag_script_info@e13fig02_YMax
    else
        YMax=new(dimsizes(season), "float")
    end if

    color_fill_exp=(/(/0.90, 0.90, 0.90/), \
                     (/0.105, 0.67, 0.8/), \
                     (/0.88, 0.88, 0.88/), \
                     (/0.66, 0.66, 0.66/), \
                     (/0.55, 0.55, 0.55/), \
                     (/0.33, 0.33, 0.33/), \
                     (/0.22, 0.22, 0.22/)/)
  
    index_mark=(/4, 6, 7, 5, 12/)
    plot_poly = new(dimsizes(diag_script_info@e13fig02_season) * dim_EXP, "graphic")
    if (dim_OBS.gt.0) then
        plot_obs = new(dimsizes(diag_script_info@e13fig02_season) * dim_OBS, "graphic")
    end if
    plot_line = new(dimsizes(diag_script_info@e13fig02_season) * dim_EXP \
                    * (dimsizes(list_nochem_mod) + 1), "graphic")

    ;----------------------------------------------------------------------------

    ip = 0
    iobs1 = 0
    in = 0
    let = (/"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", \
            "o", "p", "q", "o", "s", "t", "u", "v", "w", "x", "y", "z"/)
    iobs2 = 0

    do ipan = 0, dimsizes(season) - 1  ; panel loop
        month=(/"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", \
                "SEP", "OCT", "NOV", "DEC"/)
        ind_month = season(ipan)  
        do im = 0, 11
            if (season(ipan).eq.month(im)) then
                delete(ind_month) 
                ind_month = im + 1
            end if
        end do

        chem_t   = new((/dim_EXP, dimsizes(list_chem_mod) + 1, yearsize/), "float")
        nochem_t = new((/dim_EXP, dimsizes(list_nochem_mod) + 1, yearsize/), "float")

        chem_t_error   = new((/dim_EXP, yearsize/), "float")
        nochem_t_error = new((/dim_EXP, yearsize/), "float")

        if (dim_OBS.gt.0) then
            obs_t       = new((/dim_OBS,yearsize/), "float")
            colors_obs  = new((/dim_OBS/), typeof(colors))
            markers_obs = new((/dim_OBS/), typeof(dashes))
            thicks_obs  = new((/dim_OBS/), typeof(thicks))
            lines_obs   = new((/dim_OBS/), typeof(colors))
        end if
    
        colors_chem   = new((/dim_EXP, dimsizes(list_chem_mod) + 1/), typeof(colors))
        colors_nochem = new((/dim_EXP, dimsizes(list_nochem_mod) + 1/), typeof(colors))

        dashes_chem   = new((/dim_EXP, dimsizes(list_chem_mod) + 1/), typeof(dashes))
        dashes_nochem = new((/dim_EXP, dimsizes(list_nochem_mod) + 1/), typeof(dashes))

        thicks_chem   = new((/dim_EXP, dimsizes(list_chem_mod) + 1/), typeof(thicks))
        thicks_nochem = new((/dim_EXP, dimsizes(list_nochem_mod) + 1/), typeof(thicks))

        lines_chem    = new((/dim_EXP, dimsizes(list_chem_mod) + 1/), typeof(colors))
        lines_nochem  = new((/dim_EXP, dimsizes(list_nochem_mod) + 1/), typeof(colors))

        ;----------------------------chem------------------------------------

        iexp = 0
        nochem = 0
        chem = 0
        iobs = 0

        time_chem = new((/dim_EXP, dimsizes(list_chem_mod) + 1, yearsize/), float)

        time_c = ispan(min(toint(models@start_year)), max(toint(models@end_year)), 1)
        time_chem!0 = "exp"
        time_chem!2 = "time"
        time_chem&time = time_c
 
        time_nochem = new((/dim_EXP, dimsizes(list_nochem_mod) + 1, yearsize/), float)

        time_c = ispan(min(toint(models@start_year)), max(toint(models@end_year)), 1)
        time_nochem!0 = "exp"
        time_nochem!2 = "time"
        time_nochem&time = time_c
        if (dim_OBS.gt.0) then
            time_obs = new((/dim_OBS, yearsize/), float)
           time_o=ispan(min(toint(models@start_year)), max(toint(models@end_year)), 1)
           time_obs!1 = "time"
           time_obs&time = time_o
           obs_t!1 = "time"
           obs_t&time = time_o
        end if

        chem_t!2 = "time"
        nochem_t!2 = "time"
        chem_t&time = time_c
        nochem_t&time = time_c

        do imod = 0, dim_MOD1 - 1
            A0 = extract_data(imod, data0, -1, 0, 0)
            if (isatt(diag_script_info, "e13fig02_lev")) then
                if (max(A0&plev).ge.100000) then
                    A0&plev = A0&plev / 100
                end if
                A1 = A0(time|:, {plev|lev(ipan)}, lat|:, lon|:)
                delete(A0)
                A0 = A1
                delete(A1)
            end if
            dummy1 = area_operations(A0(time|:, lat|:, lon|:), latrange(ipan, 0), \
                                     latrange(ipan, 1), 0., 360., "average", True)
            if (var0.eq."tro3") then
               dummy1 = dummy1 * 10^-3
            end if

            delete(A0)

            if (season(ipan).ne."ANN") then
                tmp = time_operations(dummy1, toint(models@start_year(imod)),  \
                                      toint(models@end_year(imod)), "extract", \
                                      tostring(ind_month), True)
                dummy2 = time_operations(tmp, toint(models@start_year(imod)),  \
                                      toint(models@end_year(imod)), "average", \
                                      "yearly", True)
                delete(tmp)
            else
                dummy2 = time_operations(dummy1, toint(models@start_year(imod)), \
                                      toint(models@end_year(imod)), "average",  \
                                      "yearly", True)
            end if

            if (imod.ge.1) then
                if (models@project(imod).ne."OBS") then   
                    if (models@experiment(imod).ne.models@experiment(imod - 1)) then
                        iexp = iexp + 1
                        nochem = 0
                        chem = 0
                    end if
                end if
            end if

            if (any(models@name(imod).eq.list_chem_mod)) then
                chem_t(iexp, chem, {toint(models@start_year(imod)):toint(models@end_year(imod))}) = \
                    tofloat(dummy2)
                time_chem(iexp, chem, {toint(models@start_year(imod)):toint(models@end_year(imod))}) = \
                    ispan(toint(models@start_year(imod)), toint(models@end_year(imod)), 1)
                colors_chem(iexp, chem) = colors(imod)
                dashes_chem(iexp, chem) = dashes(imod)
                thicks_chem(iexp, chem) = 2.
                lines_chem(iexp, chem) = "Lines"
                chem = chem + 1
            end if

            if (any(models@name(imod).eq.list_nochem_mod)) then
                nochem_t(iexp, nochem, {toint(models@start_year(imod)):toint(models@end_year(imod))}) = \
                    tofloat(dummy2)
                time_nochem(iexp, nochem, {toint(models@start_year(imod)):toint(models@end_year(imod))}) = \
                    ispan(toint(models@start_year(imod)), toint(models@end_year(imod)), 1)
                colors_nochem(iexp, nochem) = colors(imod)
                dashes_nochem(iexp, nochem) = dashes(imod)
                thicks_nochem(iexp, nochem) = 1.
                lines_nochem(iexp, nochem) = "Lines"
                nochem = nochem + 1
            end if

            if (dim_OBS.gt.0) then
                if (any(models@name(imod).eq.list_obs_mod)) then
                    obs_t(iobs, {toint(models@start_year(imod)):toint(models@end_year(imod))}) = \
                        tofloat(dummy2)
                    time_obs(iobs, {toint(models@start_year(imod)):toint(models@end_year(imod))}) = \
                        ispan(toint(models@start_year(imod)), toint(models@end_year(imod)), 1)
                    colors_obs(iobs) = colors(imod)
                    markers_obs(iobs) = markers(imod)
                    thicks_obs(iobs) = 3.
                    lines_obs(iobs) = "Markers"
                    iobs = iobs + 1
                end if
            end if

            delete(dummy1)
            delete(dummy2)
        end do  ; (do imod)

        do iexp = 0, dim_EXP - 1
            chem_t(iexp, dimsizes(list_chem_mod), :) = \
                tofloat(dim_avg_n(chem_t(iexp, :dimsizes(list_chem_mod) - 1, :), 0))
            chem_t_error(iexp, :) = \
                tofloat(dim_stddev_n(chem_t(iexp, :dimsizes(list_chem_mod) - 1, :), 0))
            if (chem.gt.0) then
                time_chem(iexp, dimsizes(list_chem_mod), \
                    {toint(min(time_chem(iexp, :dimsizes(list_chem_mod) - 1, :))):toint(max(time_chem(iexp, :dimsizes(list_chem_mod) - 1, :)))}) = \
                    ispan(toint(min(time_chem(iexp, :dimsizes(list_chem_mod) - 1, :))), \
                    toint(max(time_chem(iexp, :dimsizes(list_chem_mod) - 1, :))) ,1)
            end if

            dashes_chem(iexp,dimsizes(list_chem_mod)) = 0
            thicks_chem(iexp,dimsizes(list_chem_mod)) = 3.
            lines_chem(iexp,dimsizes(list_chem_mod)) = "Lines"

            if (((label_exp).eq."historical").or.((label_exp).eq."refC1")) then
                colors_chem(iexp,dimsizes(list_chem_mod)) = "black"
            end if
            if ((label_exp).eq."rcp26") then
                colors_chem(iexp, dimsizes(list_chem_mod)) = "(/0., 0., 1./)"
            end if
            if ((label_exp).eq."rcp45") then
                colors_chem(iexp, dimsizes(list_chem_mod)) = "(/0, 4757, 0.7373, 1./)"
            end if
            if ((label_exp).eq."rcp60") then
                colors_chem(iexp, dimsizes(list_chem_mod)) = "(/1., 0.5098, 0.1765/)"
            end if
            if ((label_exp).eq."rcp85") then
                colors_chem(iexp, dimsizes(list_chem_mod)) = "(/1., 0., 0./)"
            end if

           ;******no chem**********************************

            nochem_t(iexp, dimsizes(list_nochem_mod), :) = \
                tofloat(dim_avg_n(nochem_t(iexp, :dimsizes(list_nochem_mod) - 1, :), 0))
            nochem_t_error(iexp, :) = \
                tofloat(dim_stddev_n(nochem_t(iexp, :dimsizes(list_nochem_mod) - 1, :), 0))
            if (nochem.gt.0) then
                time_nochem(iexp, dimsizes(list_nochem_mod), \
                    {toint(min(time_nochem(iexp, :dimsizes(list_nochem_mod) - 1, :))):toint(max(time_nochem(iexp, :dimsizes(list_nochem_mod) - 1, :)))}) = \
                    ispan(toint(min(time_nochem(iexp, :dimsizes(list_nochem_mod) - 1, :))), \
                    toint(max(time_nochem(iexp, :dimsizes(list_nochem_mod) - 1, :))), 1)
            end if

            dashes_nochem(iexp, dimsizes(list_nochem_mod)) = 1
            thicks_nochem(iexp, dimsizes(list_nochem_mod)) = 3
            lines_nochem(iexp, dimsizes(list_nochem_mod)) = "Lines"

            if (((label_exp).eq."historical").or.((label_exp).eq."refC1")) then
                colors_nochem(iexp, dimsizes(list_nochem_mod)) = "black"
            end if
            if ((label_exp).eq."rcp26") then
                colors_nochem(iexp, dimsizes(list_nochem_mod)) = "(/0., 0., 1./)"
            end if
            if ((label_exp).eq."rcp45") then
                colors_nochem(iexp,dimsizes(list_nochem_mod)) = "(/0, 4757, 0.7373, 1./)"
            end if
            if ((label_exp).eq."rcp60") then
                colors_nochem(iexp, dimsizes(list_nochem_mod)) = "(/1., 0.5098, 0.1765/)"
            end if
            if ((label_exp).eq."rcp85") then
                colors_nochem(iexp, dimsizes(list_nochem_mod)) = "(/1., 0., 0./)"
            end if
        end do

        ;--------PLOT---------------------------------------

        ti = ispan(1, (dim_EXP) * 3, 1)

        if (latrange(ipan, 0).ge.0.) then
            latstring1 = latrange(ipan,0) + "N"
        else
            latstring1 = -latrange(ipan,0) + "S"
        end if

        if (latrange(ipan, 1).ge.0.) then
            latstring2 = latrange(ipan, 1) + "N"
        else
            latstring2 = -latrange(ipan,1) + "S"
        end if

        res = True  ; plot mods desired
        res@tiYAxisString = variable_info@units  ; add units title
        res@tiXAxisString = "Year"               ; y-axis label
        res@tiMainString  = let(ipan) + ") " + variable_info@long_name \
                            + " " + season(ipan) + " " + latstring1 + "-" + latstring2
        res@gsnFrame      = False                ; don't draw yet
        res@gsnDraw       = False

        if (.not.ismissing(YMax(ipan))) then
            res@trYMaxF = YMax(ipan)
        else
            vv = (/max(chem_t), max(nochem_t)/)
            res@trYMaxF = max((/max(chem_t), max(nochem_t)/))
        end if

        if (.not.ismissing(YMin(ipan))) then
            res@trYMinF = YMin(ipan)
        else
            res@trYMinF = min((/min(chem_t), min(nochem_t)/))
        end if

        if (.not.ismissing(XMax(ipan))) then
            res@trXMaxF = XMax(ipan)
        end if

        if (.not.ismissing(XMin(ipan))) then
            res@trXMinF = XMin(ipan)
        end if

        res@vpHeightF         = 0.4                    ; change aspect ratio of plot
        res@vpWidthF          = 0.5
        res@xyLineThicknesses = ndtooned(thicks_chem)
        res@xyLineColors      = ndtooned(colors_chem)  ; line color
        res@xyDashPatterns    = ndtooned(dashes_chem)
        res@xyMarkLineModes   = ndtooned(lines_chem)

        plot(ipan) = gsn_csm_xy(wks, onedtond(ndtooned(time_chem(:, :, :)),  \
                     (/dim_EXP * (dimsizes(list_chem_mod) + 1), yearsize/)), \
                     onedtond(ndtooned(chem_t),                              \
                     (/dim_EXP * (dimsizes(list_chem_mod) + 1), yearsize/)), res)

        ;*************** add error*****************

        dim_dummy = dimsizes(time_chem(0, dimsizes(list_chem_mod), :))
        xp = new((/dim_EXP, 2 * dim_dummy/), typeof(time_chem))
        yp = new((/dim_EXP, 2 * dim_dummy/), typeof(chem_t))

        do k = 0, dim_dummy - 1
            xp(:, k) = time_chem(:, dimsizes(list_chem_mod), k)
            xp(:, dim_dummy + k) = time_chem(:, dimsizes(list_chem_mod), dim_dummy - 1 - k)
            yp(:,k) = chem_t(:, dimsizes(list_chem_mod), k) - chem_t_error(:, k)
            yp(:, dim_dummy + k) = chem_t(:, dimsizes(list_chem_mod), dim_dummy-1-k) \
                                   + chem_t_error(:, dim_dummy - 1 - k)
        end do

        do iexp = 0, dim_EXP - 1
            if (label_exp(iexp).ne."OBS") then
                gsres = True
                gsres@gsFillColor = color_fill_exp(iexp, :)
                gsres@tfPolyDrawOrder = "Draw"
                plot_poly(ip) = gsn_add_polygon(wks, plot(ipan), xp(iexp,:), yp(iexp,:), gsres)
                ip = ip + 1
            end if
        end do

        ;************add nochem**************************

        if (nochem.gt.0) then
            do iexp = 0, dim_EXP - 1
                do imod = 0, dimsizes(list_nochem_mod)
                   plres = True
                   plres@gsLineColor = colors_nochem(iexp, imod)
                   plres@gsLineDashPattern = dashes_nochem(iexp, imod)
                   plres@gsLineThicknessF = thicks_nochem(iexp, imod)
                   plot_line(in) = gsn_add_polyline(wks, plot(ipan), \
                                   time_nochem(iexp, imod, :), nochem_t(iexp, imod, :), plres)
                   in = in + 1
                end do
            end do
        end if

        ;*************** add obs *******************

         if (dim_OBS.gt.0) then
             do iobs1 = 0, dim_OBS - 1
                 pmres = True
                 pmres@gsMarkerSizeF = 0.01
                 pmres@gsMarkerThicknessF = thicks_obs(iobs1)
                 pmres@gsMarkerColor = colors_obs(iobs1)
                 pmres@gsMarkerIndex = markers_obs(iobs1)
                 plot_obs(iobs2) = gsn_add_polymarker(wks,plot(ipan), \
                                   time_obs(iobs1, :), obs_t(iobs1, :), pmres)
                 iobs2 = iobs2 + 1
            end do
        end if

        ;****************************************** 

        delete(nochem_t)
        delete(chem_t)
        delete(xp)
        delete(yp)      
        delete(ind_month)
    end do  ; loop over all "seasons"

    Pres = True
    Pres@gsnFrame = False
    Pres@gsnPanelBottom = 0.2

    gsn_panel(wks, plot, (/tointeger(dimsizes(season) / 2), 2/), Pres)
 
    if (dim_OBS.gt.0) then
        add_legenda(wks, list_obs_mod, colors_obs, markers_obs, lines_obs, thicks_obs, "0.35")
    end if

    if (diag_script_info@e13fig02_legend_MMM.eq."True") then
        chem_att = new(dim_EXP, "string")
        nochem_att = new(dim_EXP, "string")
        do iexp = 0, dim_EXP - 1
            chem_att(iexp) = attrib_chem
            nochem_att(iexp) = "_NOCHEM"
        end do

        if (nochem.gt.0) then    
            add_legenda(wks, label_exp + nochem_att, colors_nochem(:, dimsizes(list_nochem_mod)), \
                dashes_nochem(:, dimsizes(list_nochem_mod)), \
                lines_nochem(:, dimsizes(list_nochem_mod)),  \
                thicks_nochem(:, dimsizes(list_nochem_mod)), "0.08")
        end if
    end if

    ;******************a second page********************

    if (diag_script_info@e13fig02_legend.eq."True") then
        add_legenda(wks, list_chem_mod, colors_chem(0, :dimsizes(list_chem_mod) - 1), \
            dashes_chem(0, :dimsizes(list_chem_mod) - 1), \
            lines_chem(0, :dimsizes(list_chem_mod) - 1),  \
            thicks_chem(0, :dimsizes(list_chem_mod) - 1), "0.3")

        if (nochem.gt.0) then
            add_legenda(wks, list_nochem_mod, colors_nochem(0, :dimsizes(list_nochem_mod) - 1), \
                dashes_nochem(0, :dimsizes(list_nochem_mod) - 1), \
                lines_nochem(0, :dimsizes(list_nochem_mod) - 1),  \
                 thicks_nochem(0, :dimsizes(list_nochem_mod) - 1), "0.05")
        end if
    end if

    frame(wks)

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
