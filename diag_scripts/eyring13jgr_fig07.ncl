;##############################################################################
;# eyring13jgr_fig07.ncl                                                      #
;##############################################################################
;# Author: Irene Cionni (ENEA, Italy)                                         #
;# ESMVal project PI: Veronika Eyring (DLR, Germany)                          #
;##############################################################################
;#                                                                            #
;# Description:                                                               #
;#     Cacluates differences in vertically resolved ozone between the 2090s   #
;#     ("period2") and 2000s ("period1") averaged over given domains and      #
;#     seasons.                                                               #
;#                                                                            #
;# Required diag_script_info attributes (diagnostics specific):               #               
;#     e13fig07_latrange: min. and max. latitude of the regions (n-element    #
;#         array of 2-element pairs, e.g., (/(/-90,90/), (/-90,-60/)/));      #
;#         one pair of latitudes is required for each season (see below)      #
;#     styleset: style set ("DEFAULT", "CMIP5")                               #
;#     e13fig07_season: seasons (n-element array of strings, "ANN", "JAN",    #
;#         "FEB", "MAR", "DJF", "SON", etc.)                                  #
;#     e13fig07_XMin: min. x-values for plotting (n-element array, e.g.,      #
;#         (/-2., -2./)); array is required to have the same number of        #
;#         elements as "seasons" (see above)                                  #
;#     e13fig07_XMax: max. x-values for plotting (n-element array, e.g.,      #
;#         (/2., 12./)); array is required to have the same number            #
;#         of elements as "seasons" (see above)                               #
;#     e13fig07_period1: start and end year of “period1” (= 2000s), e.g.,     #
;#         (/2000., 2009/)                                                    #
;#     e13fig07_period2: start and end year of “period2” (= 2090s), e.g.,     #
;#         (/2090., 2099/)                                                    #
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
;#    * 20151208-A-laue_ax: code clean-up, added call to write_references     #
;#    * 201510??-A-cion_ir: written                                           #
;#                                                                            #
;##############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

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
load "./diag_scripts/lib/ncl/misc_function.ncl"

;************************************************

begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output("plot 07 - " + diag_script + ".ncl (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    dim_MOD = dimsizes(models@name)
    var0 = variables(0)
    field_type0 = field_types(0)
    print(diag_script_info)
    colors=project_style(diag_script_info,"colors")
    dashes=project_style(diag_script_info, "dashes")
    thicks=tofloat(project_style(diag_script_info, "thicks"))

    do idx_mod = 0, dim_MOD - 1
        data0 = read_data(idx_mod, var0, field_type0)
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

    plot_fields = (/"T3M", "T2Mz"/)
    ind_f = inlist(field_type0, plot_fields)
    if (all(ismissing(ind_f))) then
        print("fatal:can't use field, " + field_type0 + " for diag_script " + diag_script)
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

    output_filename = interface_get_figure_filename(diag_script_base, \
                var0, field_type0, "", -1)
    output_file_path = output_dir + output_filename

    color_type = "COLOR"
    wks  = gsn_open_wks (file_type, output_file_path)
    plot = new(dimsizes(diag_script_info@e13fig07_season), "graphic")
    plot_num = integertochar(ispan(97, 122, 1))

    ;-----------------------------------------------------------------

    dim_EXP = 1  ; dimsizes of different experiments
    label_exp = models@experiment(0)

    do n = 1, dimsizes(models@experiment) - 1
        if (models@experiment(n).ne.models@experiment(n - 1)) then
            dim_EXP = dim_EXP + 1
            ltmp = array_append_record(label_exp,models@experiment(n), 0)
            delete(label_exp)
            label_exp = ltmp
            delete(ltmp)
        end if
    end do

    enssize = num(ind(models@experiment.eq.models@experiment(0)))
    yearsize = max(toint(models@end_year)) - min(toint(models@start_year)) + 1
    dum = new((/dimsizes(diag_script_info@e13fig07_season), (dim_EXP) * 3/), graphic)
;    error_bar = new((/dimsizes(diag_script_info@e13fig07_season), (dim_EXP) * 3/), graphic)

    if ((isatt(diag_script_info,"list_chem_mod"))) then
        list_chem_mod = diag_script_info@list_chem_mod
        attrib_chem = "_CHEM"
        if (list_chem_mod.eq."All") then
            delete(list_chem_mod)
            list_chem_mod = models@name(:enssize) 
            attrib_chem = " "
        end if
    else 
        list_chem_mod = models@name 
        attrib_chem = " "
    end if

    print(list_chem_mod)

    if ((isatt(diag_script_info, "list_nochem_mod"))) then
        list_nochem_mod = diag_script_info@list_nochem_mod
    end if

    if ((isatt(diag_script_info, "e13fig07_mod_plot"))) then
        mod_plot = diag_script_info@e13fig07_mod_plot
    else
        mod_plot = "MMM"
    end if

    if (isvar("latrange")) then
        delete(latrange)
    end if

    if (isatt(diag_script_info, "e13fig07_latrange")) then
        latrange = diag_script_info@e13fig07_latrange
    else
        print("fatal:plot requires diag_script_info attribute 'e13fig07_latrange'")
        status_exit(1)
    end if

    if (isatt(diag_script_info, "e13fig07_season")) then
        season = diag_script_info@e13fig07_season
    else
        if (isStrSubset(field_type0, "z") .eq. False) then
            print("fatal:plot requires diag_script_info attribute 'e13fig07_season'")
            status_exit(1)
        end if
    end if

    if (isatt(diag_script_info, "e13fig07_period1")) then
        period1 = toint(diag_script_info@e13fig07_period1)
    else
        print("fatal:plot requires diag_script_info attribute 'e13fig07_period1'")
    end if

    if (isatt(diag_script_info, "e13fig07_period2")) then
        period2 = toint(diag_script_info@e13fig07_period2)
    else     
        print("fatal:plot requires diag_script_info attribute 'e13fig07_period2'")
    end if 

    if ((isatt(diag_script_info, "e13fig07_XMax"))) then
        XMax = diag_script_info@e13fig07_XMax
    else
        XMax = new(dimsizes(season), "float")
    end if

    if ((isatt(diag_script_info, "e13fig07_XMin"))) then
        XMin = diag_script_info@e13fig07_XMin
    else
        XMin=new(dimsizes(season), "float")
    end if
 
    if ((isatt(diag_script_info, "e13fig07_YMin"))) then
        YMin = diag_script_info@e13fig07_YMin
    else
        YMin=new(dimsizes(season), "float")
    end if

    if ((isatt(diag_script_info, "e13fig07_YMax"))) then
        YMax = diag_script_info@e13fig07_YMax
    else
        YMax=new(dimsizes(season), "float")
    end if
   
    color_fill_exp = (/(/0.90, 0.90, 0.90/), \
                       (/0.105, 0.67, 0.8/), \
                       (/0.88, 0.88, 0.88/), \
                       (/0.66, 0.66, 0.66/), \
                       (/0.55, 0.55, 0.55/), \
                       (/0.33, 0.33, 0.33/), \
                       (/0.22, 0.22, 0.22/)/)

    index_mark = (/4, 6, 7, 5, 12/)
    plot_poly = new(dimsizes(diag_script_info@e13fig07_season) * dim_EXP, \
                "graphic")
    plot_line = new(dimsizes(diag_script_info@e13fig07_season) * dim_EXP \
                * (dimsizes(list_nochem_mod) + 1), "graphic")

    ;----------------------------------------------------------------------------

    ip = 0
    in = 0
    let = (/"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", \
            "n", "o", "p", "q", "o", "s", "t", "u", "v", "w", "x", "y", "z"/)

    do ipan = 0, dimsizes(season) - 1  ; panel loop
        month = (/"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", \
                  "SEP", "OCT", "NOV", "DEC"/)
        ind_month = season(ipan)
        do im = 0, 11
            if (season(ipan).eq.month(im))
                delete(ind_month) 
                ind_month = im + 1
            end if
        end do  

        lev = (/1000, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 80, \
                70, 50, 30, 20, 15, 10, 7, 5, 3, 2, 1.5, 1./)
        lev!0 = "lev"
        lev&lev = lev
        lev@units = "hPa"
        lev@long_name = "Pressure"

        dim_lev = dimsizes(lev)
        chem_t   = new((/dim_EXP, dimsizes(list_chem_mod) + 1, yearsize, \
                         dim_lev/), "float")
        nochem_t = new((/dim_EXP, dimsizes(list_nochem_mod) + 1, yearsize, \
                         dim_lev/), "float")

        chem_t_error = new((/dim_EXP, yearsize, dim_lev/), "float")
        nochem_t_error = new((/dim_EXP, yearsize, dim_lev/), "float")

        colors_chem   = new((/dim_EXP, dimsizes(list_chem_mod) + 1/), \
                        typeof(colors)) 
        colors_nochem = new((/dim_EXP, dimsizes(list_nochem_mod) + 1/), \
                        typeof(colors))

        dashes_chem   = new((/dim_EXP,dimsizes(list_chem_mod)+1/), typeof(dashes))
        dashes_nochem = new((/dim_EXP,dimsizes(list_nochem_mod)+1/), typeof(dashes))
 
        thicks_chem   = new((/dim_EXP, dimsizes(list_chem_mod) + 1/), typeof(thicks))
        thicks_nochem = new((/dim_EXP, dimsizes(list_nochem_mod) + 1/), typeof(thicks))

        lines_chem    = new((/dim_EXP, dimsizes(list_chem_mod) + 1/), typeof(colors))
        lines_nochem  = new((/dim_EXP, dimsizes(list_nochem_mod) + 1/), typeof(colors))

        ;----------------------------chem------------------------------------

        iexp = 0
        nochem = 0
        chem = 0

        time_chem = new((/dim_EXP, dimsizes(list_chem_mod) + 1, yearsize/), float)
        time_c = ispan(min(toint(models@start_year)), max(toint(models@end_year)), 1)
        time_chem!0 = "exp"
        time_chem!2 = "time"
        time_chem&time = (/time_c/)

        time_nochem = new((/dim_EXP, dimsizes(list_nochem_mod) + 1, yearsize/), float)

        time_c=ispan(min(toint(models@start_year)),max(toint(models@end_year)),1)
        time_nochem!0 = "exp"
        time_nochem!2 = "time"
        time_nochem&time = (/time_c/)

        chem_t!2 = "time"
        nochem_t!2 = "time"
        chem_t&time = time_c
        nochem_t&time = time_c

        do imod = 0, dim_MOD - 1
            A0 = extract_data(imod, data0, -1, 0, 0)
            if (var0.eq."tro3") then
                var0@long_name = "O3"
            end if
            if (var0.eq."vmro3") then
                A0 = A0 * 10^9
                var0@long_name = "O3"
            end if

            if (max(A0&plev).ge.100000) then
                llev=A0&plev / 100.
                llev@units = "hPa"
                A0&plev = llev
                delete(llev)
            end if

            dummy1 = area_operations(A0(time|:, plev|:, lat|:, lon|:), \
                     latrange(ipan, 0), latrange(ipan, 1), 0., 360., "average", True)
            if (var0.eq."tro3") then
                dummy1 = dummy1 * 10^-3
            end if
            delete(A0)

            if (season(ipan).ne."ANN") then
                tmp = time_operations(dummy1, toint(models@start_year(imod)), \
                      toint(models@end_year(imod)), "extract", ind_month, True)
                dummy2 = time_operations(tmp, toint(models@start_year(imod)), \
                      toint(models@end_year(imod)), "average", "yearly", True)
                delete(tmp)
            else
                dummy2 = time_operations(dummy1, toint(models@start_year(imod)), \
                      toint(models@end_year(imod)), "average", "yearly", True)
            end if

            if (imod.ge.1) then
                if (models@experiment(imod).ne.models@experiment(imod - 1)) then
                    iexp = iexp + 1
                    nochem = 0
                    chem = 0
                end if
            end if

            if (any(models@name(imod).eq.list_chem_mod)) then
                chem_t(iexp, chem, \
                    {toint(models@start_year(imod)):toint(models@end_year(imod))}, :)\
                    = tofloat(int2p_n_Wrap(dummy2&plev, dummy2, lev, 2, 1))
                time_chem(iexp, chem, \
                    {toint(models@start_year(imod)):toint(models@end_year(imod))}) \
                    = ispan(toint(models@start_year(imod)), toint(models@end_year(imod)), 1)
                colors_chem(iexp, chem) = colors(imod)
                dashes_chem(iexp, chem) = dashes(imod)
                thicks_chem(iexp, chem) = 2.
                lines_chem(iexp, chem) = "Lines"
                chem = chem + 1
            end if

            if (any(models@name(imod).eq.list_nochem_mod)) then
                nochem_t(iexp, nochem, \
                    {toint(models@start_year(imod)):toint(models@end_year(imod))}, :) \
                    = tofloat(int2p_n_Wrap(dummy2&plev, dummy2, lev, 2, 1))
                time_nochem(iexp, nochem, \
                    {toint(models@start_year(imod)):toint(models@end_year(imod))}) \
                    = ispan(toint(models@start_year(imod)), toint(models@end_year(imod)), 1)
                colors_nochem(iexp, nochem) = colors(imod)
                dashes_nochem(iexp, nochem) = dashes(imod)
                thicks_nochem(iexp, nochem) = 1.
                lines_nochem(iexp, nochem) = "Lines"
                nochem = nochem + 1
            end if

            delete(dummy1)
            delete(dummy2)
        end do  ; loop over models

        do iexp = 0, dim_EXP - 1
            chem_t(iexp,dimsizes(list_chem_mod), :, :) = \
                tofloat(dim_avg_n(chem_t(iexp, :dimsizes(list_chem_mod) - 1, :, :), 0))
            chem_t_error(iexp, :, :) = \
                tofloat(dim_stddev_n(chem_t(iexp, :dimsizes(list_chem_mod) - 1, :, :), 0))
            if (chem.gt.0) then
                time_chem(iexp, dimsizes(list_chem_mod), \
                    {toint(min(time_chem(iexp, :dimsizes(list_chem_mod) - 1, \
                    :))):toint(max(time_chem(iexp,:dimsizes(list_chem_mod) - 1, :)))}) \
                    = ispan(toint(min(time_chem(iexp, :dimsizes(list_chem_mod) - 1, :))), \
                    toint(max(time_chem(iexp, :dimsizes(list_chem_mod) - 1, :))), 1)
            end if
            dashes_chem(iexp, dimsizes(list_chem_mod)) = 0 
            thicks_chem(iexp, dimsizes(list_chem_mod)) = 3.
            lines_chem(iexp, dimsizes(list_chem_mod)) = "Lines"
            if ((label_exp(iexp)).eq."historical") then
                colors_chem(iexp, dimsizes(list_chem_mod)) = "black"
            end if
            if ((label_exp(iexp)).eq."rcp26") then
                colors_chem(iexp, dimsizes(list_chem_mod)) = "(/0., 0., 1./)"
            end if
            if ((label_exp(iexp)).eq."rcp45") then
                colors_chem(iexp, dimsizes(list_chem_mod)) = "(/0, 4757, 0.7373, 1./)"
            end if
            if ((label_exp(iexp)).eq."rcp60") then
                colors_chem(iexp, dimsizes(list_chem_mod)) = "(/1., 0.5098, 0.1765/)"
            end if
            if ((label_exp(iexp)).eq."rcp85") then
                colors_chem(iexp, dimsizes(list_chem_mod)) = "(/1., 0., 0./)"
            end if 

            ;******no chem**********************************

            nochem_t(iexp, dimsizes(list_nochem_mod), :, :) = \
                tofloat(dim_avg_n(nochem_t(iexp, :dimsizes(list_nochem_mod) - 1, :, :), 0))
            nochem_t_error(iexp, :, :) = \
                tofloat(dim_stddev_n(nochem_t(iexp, :dimsizes(list_nochem_mod) - 1, :, :), 0))
            if (nochem.gt.0) then 
                time_nochem(iexp, dimsizes(list_nochem_mod), \
                    {toint(min(time_nochem(iexp, :dimsizes(list_nochem_mod) - 1, \
                    :))):toint(max(time_nochem(iexp, :dimsizes(list_nochem_mod) - 1, :)))}) \
                    = ispan(toint(min(time_nochem(iexp, :dimsizes(list_nochem_mod) - 1, :))), \
                    toint(max(time_nochem(iexp, :dimsizes(list_nochem_mod) - 1, :))), 1)
            end if
            dashes_nochem(iexp, dimsizes(list_nochem_mod)) = 1
            thicks_nochem(iexp, dimsizes(list_nochem_mod)) = 3
            lines_nochem(iexp, dimsizes(list_nochem_mod)) = "Lines"
            if ((label_exp(iexp)).eq."historical") then
                colors_nochem(iexp, dimsizes(list_nochem_mod)) = "black"
            end if
            if ((label_exp(iexp)).eq."rcp26") then
                colors_nochem(iexp, dimsizes(list_nochem_mod)) = "(/0., 0., 1./)"
            end if
            if ((label_exp(iexp)).eq."rcp45") then
                colors_nochem(iexp, dimsizes(list_nochem_mod)) = "(/0, 4757, 0.7373, 1./)"
            end if
            if ((label_exp(iexp)).eq."rcp60") then
                colors_nochem(iexp, dimsizes(list_nochem_mod)) = "(/1., 0.5098, 0.1765/)"
            end if
            if ((label_exp(iexp)).eq."rcp85") then
                colors_nochem(iexp, dimsizes(list_nochem_mod)) = "(/1., 0., 0./)"
            end if
        end do  ; loop over experiments

        if (label_exp(0).eq."historical") then
            do iexp = 1, dim_EXP - 1
                do imod = 0, dimsizes(list_chem_mod) 
                    chem_t(iexp, imod, :, :) = \
                        (/dim_avg_n_Wrap((/chem_t(0, imod, :, :), \
                        chem_t(iexp, imod, :, :)/), 0)/)
                end do
            end do
            istart = 1
        else
            istart = 0  
        end if

        ;***********period average****************************

        print(period1)
        period1_chem_t = dim_avg_n_Wrap(chem_t(istart:, :, \
            {period1(0):period1(1)}, :), 2)
        period1_nochem_t = dim_avg_n_Wrap(nochem_t(istart:, :, \
            {period1(0):period1(1)}, :), 2)
        period2_chem_t = dim_avg_n_Wrap(chem_t(istart:, :, \
            {period2(0):period2(1)}, :), 2)
        period2_nochem_t = dim_avg_n_Wrap(nochem_t(istart:, :, \
            {period2(0):period2(1)}, :), 2)
        diff_chem = period2_chem_t - period1_chem_t
        printVarSummary(diff_chem) 
        diff_nochem = period2_nochem_t - period1_nochem_t
        diff_chem = (/diff_chem * (conform(diff_chem, lev, 2) / 10)/)
        diff_nochem = (/diff_nochem * (conform(diff_nochem, lev, 2) / 10)/)

        ;***************************************************

        ;--------PLOT---------------------------------------

        ti = ispan(1, (dim_EXP) * 3, 1)

        if (latrange(ipan, 0).ge.0.) then
            latstring1 = latrange(ipan,0) + "N"
        else
            latstring1 = -latrange(ipan,0) + "S"
        end if
       
        if (latrange(ipan, 1).ge.0.) then
            latstring2 = latrange(ipan,1) + "N"
        else 
            latstring2 = -latrange(ipan, 1) + "S"
        end if

        res                   = True  ; plot mods desired
        res@tmYLMode          = "Explicit"
        res@tmYLValues        = (/1,3,10,30,100,300,500,1000/)
        res@tmYLLabels        = (/1,3,10,30,100,300,500,1000/)
        res@tmLabelAutoStride = True
        res@trYAxisType       = "LogAxis"
        res@trYReverse        = True
        res@tmXMajorGrid      = True                   ; implement x grid
        res@tmXMajorGridThicknessF = 1.0               ; 2.0 is default
        res@tmXMajorGridLineDashPattern = 2            ; select short dash lines
        res@tiXAxisString     = var0@long_name + " (mPa/Pa)" ; add units title
        res@tiYAxisString     = "Pressure (hPa)"       ; y-axis label      
        res@tiMainString      = let(ipan) + ")  " + season(ipan) \
                                + " " + latstring1 + "-" + latstring2
        res@gsnFrame          = False                  ; don't draw yet
        res@gsnDraw           = False  

        if (.not.ismissing(XMax(ipan))) then
            res@trXMaxF = XMax(ipan)
        else
            vv=(/max(diff_chem), max(diff_nochem)/)
            res@trXMaxF = max((/max(diff_chem), max(diff_nochem)/))
        end if

        if (.not.ismissing(XMin(ipan))) then
            res@trXMinF = XMin(ipan)
        else
            res@trXMinF = min((/min(diff_chem), min(diff_nochem)/))
        end if

        res@vpHeightF = 0.4                    ; change aspect ratio of plot
        res@vpWidthF  = 0.5

        print(onedtond(ndtooned(diff_chem), (/(dim_EXP - 1) * \
            (dimsizes(list_chem_mod) + 1), dimsizes(lev)/)))

        if (mod_plot.eq."All") then
            res@xyLineThicknesses = ndtooned(thicks_chem(istart:, :))
            res@xyLineColors = ndtooned(colors_chem(istart:, :))  ; line color
            res@xyDashPatterns = ndtooned(dashes_chem(istart:, :))
            plot(ipan) = gsn_csm_xy(wks, onedtond(ndtooned(diff_chem), \
                (/(dim_EXP - istart) * (dimsizes(list_chem_mod) + 1),  \
                dimsizes(lev)/)), lev, res)
        else
            if (mod_plot.eq."MMM") then
                res@xyLineThicknesses = thicks_chem(istart:, dimsizes(list_chem_mod))
                res@xyLineColors = colors_chem(istart:, dimsizes(list_chem_mod))  ; line color
                res@xyDashPatterns = dashes_chem(istart:, dimsizes(list_chem_mod))
                plot(ipan) = gsn_csm_xy(wks, \
                    onedtond(ndtooned(diff_chem(:, dimsizes(list_chem_mod), :)), \
                    (/(dim_EXP - istart), dimsizes(lev)/)), lev, res)
            end if
            if (mod_plot.eq."Models") then
                res@xyLineThicknesses = ndtooned(thicks_chem(istart:, \
                    dimsizes(list_chem_mod) - 1))
                res@xyLineColors = ndtooned(colors_chem(istart:, \
                    dimsizes(list_chem_mod) - 1))   ; line color
                res@xyDashPatterns = ndtooned(dashes_chem(istart:, \
                    dimsizes(list_chem_mod) - 1))
                plot(ipan) = gsn_csm_xy(wks, onedtond(ndtooned(diff_chem(:, \
                    :list_chem_mod - 1, :)), (/(dim_EXP - istart) * \
                    (dimsizes(list_chem_mod)), dimsizes(lev)/)), lev, res)
            end if
        end if

        ;************add nochem**************************

        if (nochem.gt.0) then
            do iexp = istart, dim_EXP - 1
                if (mod_plot.eq."All") then
                    do imod = 0, dimsizes(list_nochem_mod)
                        plres = True
                        plres@gsLineColor = colors_nochem(iexp, imod)
                        plres@gsLineDashPattern = dashes_nochem(iexp, imod)
                        plres@gsLineThicknessF = thicks_nochem(iexp, imod)
                        plot_line(in) = gsn_add_polyline(wks, plot(ipan), \
                            diff_nochem(iexp, imod, :), lev, plres)
                        in = in + 1
                    end do
                else
                    if (mod_plot.eq."MMM") then
                        plres = True
                        plres@gsLineColor = colors_nochem(iexp, dimsizes(list_nochem_mod))
                        plres@gsLineDashPattern = dashes_nochem(iexp, dimsizes(list_nochem_mod))
                        plres@gsLineThicknessF = thicks_nochem(iexp, dimsizes(list_nochem_mod))
                        plot_line(in) = gsn_add_polyline(wks, plot(ipan), \
                            diff_nochem(iexp, dimsizes(list_nochem_mod), :), lev, plres)
                        in = in + 1
                    end if
                    if (mod_plot.eq."Models") then
                        do imod = 0, dimsizes(list_nochem_mod) - 1
                            plres = True
                            plres@gsLineColor = colors_nochem(iexp, imod)
                            plres@gsLineDashPattern = dashes_nochem(iexp, imod)
                            plres@gsLineThicknessF = thicks_nochem(iexp, imod)
                            plot_line(in) = gsn_add_polyline(wks, plot(ipan),\
                                diff_nochem(iexp, imod,:), lev, plres)
                            in = in + 1
                        end do
                    end if
                end if 
            end do  ; loop over experiments
        end if
    end do  ; loop over seasons

    Pres                = True
    Pres@gsnFrame       = False
    Pres@gsnPanelBottom = 0.4

    gsn_panel(wks, plot, (/tointeger(dimsizes(season)/ 2), 2/), Pres)

    if (mod_plot.eq."All") then
        add_legenda(wks, list_chem_mod, colors_chem(0, :dimsizes(list_chem_mod) - 1), \
            dashes_chem(0, :dimsizes(list_chem_mod) - 1), \
            lines_chem(0, :dimsizes(list_chem_mod) - 1),  \
            thicks_chem(0, :dimsizes(list_chem_mod)) - 1, "0.38")

        if (nochem.gt.0) then
            add_legenda(wks, list_nochem_mod, \
                colors_nochem(0, :dimsizes(list_nochem_mod) - 1), \
                dashes_nochem(0, :dimsizes(list_nochem_mod) - 1), \
                lines_nochem(0, :dimsizes(list_nochem_mod) - 1), \
                thicks_nochem(0, :dimsizes(list_nochem_mod) - 1), "0.3")
        end if

        add_legenda(wks, label_exp(istart:), colors_chem(istart:, \
            dimsizes(list_chem_mod)), dashes_chem(istart:, \
            dimsizes(list_chem_mod)), lines_chem(istart:, \
            dimsizes(list_chem_mod)), thicks_chem(istart, \
            dimsizes(list_chem_mod)), "0.38")

        if (nochem.gt.0) then
            add_legenda(wks, label_exp(istart:), \
                colors_nochem(istart:, dimsizes(list_nochem_mod)), \
                dashes_nochem(istart:, dimsizes(list_nochem_mod)), \
                lines_nochem(istart, dimsizes(list_nochem_mod)), \
                thicks_nochem(istart, dimsizes(list_nochem_mod)), "0.3")
        end if
    end if

    if (mod_plot.eq."MMM") then
        add_legenda(wks, label_exp(istart:), \
            colors_chem(istart:, dimsizes(list_chem_mod)), \
            dashes_chem(istart:, dimsizes(list_chem_mod)), \
            lines_chem(istart:, dimsizes(list_chem_mod)), \
            thicks_chem(istart:, dimsizes(list_chem_mod)), "0.38")

        if (nochem.gt.0) then
            add_legenda(wks, label_exp(istart:), \
                colors_nochem(istart:, dimsizes(list_nochem_mod)), \
                dashes_nochem(istart:, dimsizes(list_nochem_mod)), \
                lines_nochem(istart, dimsizes(list_nochem_mod)), \
                thicks_nochem(istart, dimsizes(list_nochem_mod)), "0.3")
        end if
    end if

    if (mod_plot.eq."Models") then
        add_legenda(wks, list_chem_mod, \
            colors_chem(0, :dimsizes(list_chem_mod) - 1), \
            dashes_chem(0, :dimsizes(list_chem_mod) - 1), \
            lines_chem(0, :dimsizes(list_chem_mod) - 1), \
            thicks_chem(0, :dimsizes(list_chem_mod) - 1), "0.38")

        if (nochem.gt.0) then
            add_legenda(wks, list_nochem_mod, \
                colors_nochem(0, :dimsizes(list_nochem_mod) - 1), \
                dashes_nochem(0, :dimsizes(list_nochem_mod) - 1), \
                lines_nochem(0, :dimsizes(list_nochem_mod) - 1), \
                thicks_nochem(0, :dimsizes(list_nochem_mod) - 1), "0.3")
        end if
    end if

    frame(wks)
end
