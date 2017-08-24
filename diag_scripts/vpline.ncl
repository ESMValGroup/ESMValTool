;#############################################################################
;# Vertical PROFILES                                 
;# ESMVal project                                                             
;#############################################################################
;#                                                                            
;# Description: Produces vertical profiles according to Eyring et al. 2006
;#              Figure 5 upper panels following eyring06jgr_fig05.ncl
;#                                                                            
;# Required variable_info attributes (variable specific):
;#          units
;#                                                                            
;# Required diag_script_info attributes (diagnostics specific):   
;#        e06fig05_multimean
;#        e06fig05_levmax 
;#        e06fig05_levmin
;#        e06fig05_XMax_latsel 
;#        e06fig05_XMin_latsel 
;#        e06fig05_latsel
;#        e06fig05_monthsel 
;#        e06fig05_levsel                   
;#                                                                            
;# Caveats:      
;#                                                              
;# Modification history:                                                      
;#    * 20170824-A_bock_ls added tags for reporting
;#    * 20160926-A_laue_ax: 1) removed erroneous, hard-wired "fixing" of vertical
;#                          coordinate (plev) and introduced fix-files
;#                          (in reformat_scripts/fixes) for BDBP and ESACCI-OZONE
;#                          2) bugfix: multi-model mean only calculated from models
;#                          even if more than one reference data set is used
;#                          (new attribute: diag_script_info@models_to_skip to 
;#                          specify list of data sets to skip when calculating
;#                          multi-model mean)
;#    * 20160309-A_wenz_sa: written based on eyring06jgr_fig05.ncl
;#############################################################################
load "./diag_scripts/lib/ncl/misc_function.ncl"

load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/set_operators.ncl"
load "./diag_scripts/lib/ncl/scaling.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/regridding.ncl"
load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/ensemble.ncl"
load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/misc_function.ncl"
load "./diag_scripts/calc_seasonal.ncl"
load "./diag_scripts/uniform_coordinates.ncl"
load "./plot_scripts/ncl/aux_plotting.ncl"

load "./diag_scripts/lib/ncl/meta_data.ncl"

begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output("plot "+diag_script+" (var: "+variables(0)+")",verbosity,1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)


    ;; References Information
    ;write_references(diag_script,      \  ; script name
    ;                 "A_wenz_sa",      \  ; authors
    ;                 "",               \  ; contributors
    ;                 "D_esacci",  \  ; diag_references
    ;                 "E_esacci",       \  ; obs_references
    ;                 "P_cmug")            ; proj_references

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)
   
    colors = project_style(diag_script_info,"colors")
    dashes = project_style(diag_script_info, "dashes")
    thicks = tofloat(project_style(diag_script_info, "thicks")) 
 
    do idx_mod = 0, dim_MOD - 1
        data = read_data(idx_mod, var0, field_type0)
    end do
end

begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin

    plot_fields = (/"T3M", "T2Mz"/)
    ind_f = inlist(field_type0, plot_fields)
    if (all(ismissing(ind_f))) then
        error_msg("f", diag_script, "", \
                  "can't use field, " + field_type0 + \
                  " for diag_script " + diag_script)
    end if

    delete(plot_fields)
    delete(ind_f)
    if(isatt(variable_info,"short_name")) then
        varname=variable_info@short_name
    else
       varname=var0
    end if
    plot_dir = getenv("ESMValTool_plot_dir")
    diag_script_base = basename(diag_script)
    output_dir = get_output_dir(plot_dir, diag_script_base)

    file_type = getenv("ESMValTool_output_file_type")
    if(ismissing(file_type)) then
        file_type = "ps"
    end if
    file_type = str_lower(file_type)

    color_type = "COLOR"

    ;; Check requried diag_script_info attributes
    req_atts1 = (/"e06fig05_multimean", "e06fig05_levmax", "e06fig05_levmin", \
                  "e06fig05_latsel", "e06fig05_monthsel", "e06fig05_XMax", "e06fig05_XMin"/)
    if (field_type0 .ne. "T0M") then
        exit_if_missing_atts(diag_script_info, req_atts1)
    end if
    
    latsel    = diag_script_info@e06fig05_latsel
    monthsel     = diag_script_info@e06fig05_monthsel
    XMax      = diag_script_info@e06fig05_XMax
    XMin      = diag_script_info@e06fig05_XMin
    Max_lev   = diag_script_info@e06fig05_levmax
    Min_lev   = diag_script_info@e06fig05_levmin
    multimean = diag_script_info@e06fig05_multimean

    if (dimsizes(monthsel).ne.dimsizes(latsel)) then
        error_msg("f", diag_script, "", \
                  "plot requires dimension 'e06fig05_monthsel' " + \
                  "equal last dimension 'e06fig05_latsel' ")       
    end if

    if (isatt(diag_script_info,"models_to_skip")) then
        refModels = diag_script_info@models_to_skip
        dim_OBS   = dimsizes(refModels)
        idx_ref   = new(dim_MOD - dim_OBS, integer)
        i = 0
        do imod = 0, dim_MOD - 1
            if (any(models@name(imod) .eq. refModels)) then
                if (multimean) then
                    print("skipped model " + models@name(imod) + " in multi-model mean calculation")
                end if
            else
                idx_ref(i) = imod
                i = i + 1
            end if
        end do
    else
        refModels = "noref"
        dim_OBS   = 0
        idx_ref   = ind(models@name.eq.models@name)
    end if

    if (multimean) then
        dimTot = dim_MOD+1
    else
        dimTot = dim_MOD
    end if

    if(isatt(variable_info, "units")) then
        UNITS1 = variable_info@units
    end if

    diag_script_base = basename(diag_script)
    aux_title_info   = str_join(monthsel,"_")
   
    ;; Define some variables
    nPanel = dimsizes(monthsel)

    plot_num = (/"a", "b", "c", "d", "e", "f", "g", "h", "i", "l", \
                 "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "z"/)

    rad  = 4.0 * atan(1.0) / 180.0
    lev0 = 100.*(/0.1, 0.2, 0.3, 0.5, 1, 1.5, 2., 3., 5., 7., 10., \
                  15., 20., 30., 50., 70., 80., 90., 100., 115., 120., \
                  150., 170., 200., 250., 300., 400., 500., 700., 850., 1000./)
    nlat   = 73
    lat0   = latGlobeF(nlat, "lat", "latitude", "degrees_north")
    dimLev = dimsizes(lev0)
   
    tmpData = new((/dim_MOD+dim_OBS, 12,dimsizes(lev0),dimsizes(lat0),2/), "float")
    tmpData!0 = "line"
    tmpData!1 = "month"
    tmpData!2 = "lev"
    tmpData!3 = "lat"
    tmpData!4 = "stat"

    Xdata = new((/dimTot,12,dimsizes (lev0),dimsizes(lat0),2/), "float")
    Xdata!0 = "line"
    Xdata!1 = "month"
    Xdata!2 = "lev"
    Xdata!3 = "lat"
    Xdata!4 = "stat"
    Xdata_month = ispan(1,12,1)
    Xdata_plev  = lev0
    Xdata_lat   = lat0
   
    Xref        = new((/nPanel,12,dimsizes(lev0),2/), "float")
    plotMod     = new(dimTot,typeof(models@name))
    plotColor   = new(dimTot,typeof(models@name))
    plotDash    = new(dimTot,typeof(models@name))
    plotScatters = new(dimTot,typeof(models@name))
    LineScatters = new(dimTot,typeof(models@name)) 
    LineTick     = new(dimTot,"float") 

    plotColor(:dim_MOD-1)    = colors
    plotDash(:dim_MOD-1)     = dashes
    plotScatters(:dim_MOD-1) = 0.
    LineScatters(:dim_MOD-1) = "Lines"
    LineTick(:dim_MOD-1)     = 2.

    ;; MAIN------------------------------------
    do imod = 0, dim_MOD - 1
        A0 = extract_data(imod, data, -1, 0, 0)

        ;; Convert units if appropriate
        ;if (isatt(variable_info, "plot_units")) then
        ;    A0 = convert_units(A0, variable_info@plot_units)
        ;    UNITS = variable_info@plot_units
        ;else
        ;    UNITS = UNITS1
        ;end if

        if (A0@units .eq. "mole mole-1" .or. A0@units .eq. "1e-9")then
           A0    = A0*10^6
           UNITS = "ppmv"
        else
           UNITS = UNITS1
        end if

        uniform_coordinates(A0, field_type0)
        plev_n = A0&plev 

        if (field_type0 .eq. "T3M") then
            tmp    = clmMonTLL(dim_avg_n_Wrap(A0,3))
            tmp_st = stdMonTLL(dim_avg_n_Wrap(A0,3))
        end if
        if (field_type0 .eq. "T2Mz") then
            tmp    = clmMonTLL(A0)
            tmp_st = stdMonTLL(A0)
        end if
        tmp2          = new((/12, dimsizes(tmp&plev),dimsizes(tmp&lat),2/),typeof(tmp))
        tmp2(:,:,:,0) = tmp
        tmp2(:,:,:,1) = tmp_st
        tmp2          = where(tmp2.ge.10^20,tmp2@_FillValue,tmp2) 

        if(plev_n(1).lt.plev_n(0)) then
            plev_nn = plev_n(::-1)
            delete(plev_n)
            plev_n = plev_nn
            tmp_n  = tmp(:,::-1,:)   
            delete(tmp)
            tmp = tmp_n
            delete([/tmp_n,plev_nn/])               
        end if

        Aint        = int2p_n_Wrap(plev_n,tmp,lev0,2,1)
        Aint_st     = int2p_n_Wrap(plev_n,tmp_st,lev0,2,1)
        Aint_lat    = linint1_n_Wrap(Aint&lat,Aint,False,lat0,0,2)
        Aint_lat_st = linint1_n_Wrap(Aint_st&lat,Aint_st,False,lat0,0,2)
        tmpData(imod,:,:,:,0) = (/Aint_lat/)
        tmpData(imod,:,:,:,1) = (/Aint_lat_st/)
        tmpData = where(tmpData.ge.10^20,tmpData@_FillValue,tmpData)
        delete(A0)

        Xdata(imod,:,:,:,:)         = tmpData(imod,:,:,:,:) 
        plotMod(imod)               = models@name(imod)
; A_laue_ax+
;       add ensemble member to model name in legend if ensemble is not r1i1p1
        if (models@ensemble(imod).ne."r1i1p1") then
            if (isStrSubset(models@ensemble(imod), "r") .and. \
                isStrSubset(models@ensemble(imod), "i") .and. \
                isStrSubset(models@ensemble(imod), "p")) then
                plotMod(imod) = models@name(imod) + " (" + models@ensemble(imod) + ")"
            end if
        end if
; A_laue_ax-
        Xdata@$plotMod(imod)$       = tmp2
        Xdata_month@$plotMod(imod)$ = tmp2&month
        Xdata_plev@$plotMod(imod)$  = tmp2&plev
        Xdata_lat@$plotMod(imod)$   = tmp2&lat

        ;Xdata@plotMod(imod)=tmp2
        delete([/tmp,tmp2,tmp_st,Aint,Aint_st,Aint_lat,Aint_lat_st,plev_n/])         
             
   end do
    
    month_n = month_sel(monthsel)

   if (multimean) then
       Xdata(dimTot-1,:,:,:,:) = (/dim_avg_n_Wrap(Xdata(idx_ref,:,:,:,:),0)/)

       plotMod(dimTot-1)      = "MMM"
       plotColor(dimTot-1)    = "Red"
       plotDash(dimTot-1)     = 0
       LineScatters(dimTot-1) = "Lines"
       LineTick(imod)         = 3.

       MMM       = (/dim_avg_n_Wrap(Xdata(idx_ref,:,:,:,0),0)/)
       X_MMM_std = (/dim_stddev_n(Xdata(idx_ref,:,:,:,0),0)/)      
       tmp2_MMM  = new((/12, dimsizes(lev0),dimsizes(lat0),2/),typeof(MMM))

       tmp2_MMM(:,:,:,0) = MMM
       tmp2_MMM(:,:,:,1) = X_MMM_std
       Xdata@MMM         = tmp2_MMM
       delete(tmp2_MMM)
       Xdata_month@MMM   = ispan(1,12,1)
       Xdata_plev@MMM    = lev0
       Xdata_lat@MMM     = lat0
   end if
   
   plot_d      = new(dimsizes(latsel),"graphic")
   plot_xp_R_d = new(dimsizes(latsel),"graphic")  
   plot_xp_M_d = new(dimsizes(latsel),"graphic") 

   ;; Define workstation
   if isatt(models,"experiment") then
       outfile = var0 + "_" + models@experiment(0) + "_" + models@start_year(0) + "-" + models@end_year(0)
   else
       outfile = var0 + "_" + models@start_year(0) + "-" + models@end_year(0)
   end if
   wks = get_wks(0, diag_script, outfile)

   do ipanel=0, dimsizes(latsel)-1
       res=True   
       res@gsnDraw           = False
       res@gsnFrame          = False 
       res@trYReverse        = True    
       res@trXMaxF           = XMax(ipanel)
       res@trXMinF           = XMin(ipanel)
       res@trYMaxF           = Max_lev
       res@trYMinF           = Min_lev
       res@tmYLMode          = "Explicit"
       res@tmYLValues        = (/0,1,2,5,7,10,30,50,100,200,500,850/)
       res@tmYLLabels        = (/0,1,2,5,7,10,30,50,100,200,500,850/)
       res@tmLabelAutoStride = True
       res@vpHeightF         = 0.6
       res@vpWidthF          = 0.5
       res@tmLabelAutoStride = True

       res@tiYAxisString     = "Pressure (hPa)"
       res@tiXAxisString     = varname+" "+"("+UNITS+")" 
       if (latsel(ipanel).le.0)
           latsel1 = -latsel(ipanel)+"S" 
       end if
       
       if (latsel(ipanel).ge.0)
           latsel1 = latsel(ipanel)+"N"
       end if
       
       if (latsel(ipanel).eq.0)
           latsel1 = "EQ"
       end if
       
       res@tiMainString       = varname+" "+latsel1+" "+monthsel(ipanel)      
       res@xyMarkLineModes    = LineScatters(0)
       res@trYAxisType        = "LogAxis"

       res@xyLineThicknessF   = 2.0
       res@xyLineDashSegLenF  = 0.13  
       res@xyLineColors       = plotColor(0)
       res@xyDashPatterns     = plotDash(0)
       res@xyMarkers          = plotDash(0) 
       res@xyMarkerSizeF      = 0.004
       res@xyMarkerThicknessF = 1.8

       tmp = Xdata@$plotMod(0)$
       tmp = where(tmp.ge.10^20,Xdata@_FillValue,tmp) 
       tmp!0 = "month"
       tmp!1 = "plev"
       tmp!2 = "lat"
       tmp!3 = "stat"
       tmp&month = Xdata_month@$plotMod(0)$
       tmp&plev  = Xdata_plev@$plotMod(0)$/100.
       tmp&lat   = Xdata_lat@$plotMod(0)$
       plot_d(ipanel) = gsn_csm_xy (wks,tmp(month|month_n(ipanel)-1,plev|:,\
                                        {lat|latsel(ipanel)},stat|0),tmp&plev,res)
       delete(tmp)

       do  imodi = 1,dimTot-1 

           imod = dimTot-imodi  
           tmp  = Xdata@$plotMod(imod)$
; A-laue_ax+
           tmp@_FillValue = Xdata@_FillValue
; A-laue_ax-
           tmp  = where(tmp.ge.10^20,tmp@_FillValue,tmp)
           tmp1   = tmp(:,:,:,0)
           tmp_st = tmp(:,:,:,1)  
           delete(tmp)
           tmp    = tmp1
           delete(tmp1)
 
           tmp!0 = "month"
           tmp!1 = "plev"
           tmp!2 = "lat"
     
           tmp&month = Xdata_month@$plotMod(imod)$
           tmp&plev  = Xdata_plev@$plotMod(imod)$/100.
           tmp&lat   = Xdata_lat@$plotMod(imod)$
     
           tmp_st!0 = "month"
           tmp_st!1 = "plev"
           tmp_st!2 = "lat"
                    
           tmp_st&plev = Xdata_plev@$plotMod(imod)$/100.
           tmp_st&lat  = Xdata_lat@$plotMod(imod)$

           if (LineScatters(imod).eq."Markers")  then

               xp_R = new((/2*dimsizes(tmp&plev)/),typeof(tmp))
               yp_R = new((/2*dimsizes(tmp&plev)/),typeof(tmp)) 

               if (latsel(ipanel).gt.0) then
                   ltmp    = tmp(month|month_n(ipanel)-1,plev|:,{lat|(latsel(ipanel)+0.1)})
                   ltmp_st = tmp_st(month|month_n(ipanel)-1,plev|:,{lat|(latsel(ipanel)+0.1)})
               else
                   ltmp    = tmp(month|month_n(ipanel)-1,plev|:,{lat|(latsel(ipanel)-0.1)})
                   ltmp_st = tmp_st(month|month_n(ipanel)-1,plev|:,{lat|(latsel(ipanel)-0.1)})
               end if

               do j = 0,dimsizes(tmp&plev)-1
                   xp_R(j) = ltmp(j) + ltmp_st(j) 
                   yp_R(j) = tofloat(tmp&plev(j))
                   yp_R(dimsizes(tmp&plev)+j) = tofloat(tmp&plev(dimsizes(tmp&plev)-1-j))      
                   xp_R(dimsizes(tmp&plev)+j) = ltmp(dimsizes(tmp&plev)-1-j)-ltmp_st(dimsizes(tmp&plev)-1-j)
               end do
               gsres = True    
               gsres@gsFillColor = "Snow3"
     
               plot_xp_R_d(ipanel) = gsn_add_polygon(wks,plot_d(ipanel),xp_R,yp_R,gsres)

               lplev = tmp&plev
               lcol  = plotColor(imod)
               lsca  = plotDash(imod)

               add_scatt(wks,plot_d(ipanel),ltmp,lplev, lcol,lsca ) 
               
               delete([/xp_R,yp_R,gsres,ltmp,ltmp_st,lplev,lcol,lsca/]) 

           else
               xp_M = new (2*dimsizes(tmp&plev),typeof(tmp))
               yp_M = new (2*dimsizes(tmp&plev),typeof(tmp))

               do j = 0,dimsizes(tmp&plev)-1         
                   xp_M(j) = tmp(month|month_n(ipanel)-1,plev|j,{lat|latsel(ipanel)}) + \
                             tmp_st(month|month_n(ipanel)-1,plev|j,{lat|latsel(ipanel)})
                   yp_M(j) = tofloat(tmp&plev(j))
                   yp_M(dimsizes(tmp&plev)+j) = tofloat(tmp&plev(dimsizes(tmp&plev)-1-j))      
                   xp_M(dimsizes(tmp&plev)+j) = tmp(month|month_n(ipanel)-1,\
                                                    plev|(dimsizes(tmp&plev)-1-j),{lat|latsel(ipanel)})-\
                                                tmp_st(month|month_n(ipanel)-1, \
                                                       plev|(dimsizes(tmp&plev)-1-j),{lat|latsel(ipanel)})
               end do
               if(plotMod(imod).eq."MMM") then
                   gsres = True
                   gsres@gsFillColor = "cyan"
                   gsres@tfPolyDrawOrder   = "Predraw"
                   plot_xp_M_d(ipanel) = gsn_add_polygon(wks,plot_d(ipanel),xp_M,yp_M,gsres)
               end if
               add_line(wks,plot_d(ipanel),tmp(month|month_n(ipanel)-1,plev|:,{lat|latsel(ipanel)}),\
                                           tmp&plev, plotColor(imod),plotDash(imod) )
               delete([/xp_M,yp_M/])          
           end if
           delete([/tmp,tmp_st/]) 
       end do
       
       ;draw(plot_d(ipanel))
       ;frame(wks)
   end do
;   delete(res)

   ;; Call legend-plotting routine
   Pres = True
;   Pres@gsnFrame       = False
   Pres@gsnPanelBottom = 0.2
   Pres@gsnPanelFigureStrings=plot_num(0:dimsizes(latsel)-1) 
   Pres@gsnPanelFigureStringsFontHeightF = 0.02
   Pres@gsnPanelFigureStringsPerimOn = "False"
   Pres@gsnPanelRowSpec = True
   Pres@gsnPanelXWhiteSpacePercent=5.
   Pres@gsnPanelYWhiteSpacePercent=3.
   Pres@amJust="TopLeft"

   ;; Call panel routine
   add_legenda(wks,plotMod,plotColor,tofloat(plotDash),LineScatters, LineTick,"0.4")
   gsn_panel(wks,plot_d,dimsizes(latsel),Pres)      

;**************************************************************
   ; add meta data to plot (for reporting)

   climofiles = new(dim_MOD, string)
   do ii = 0, dim_MOD - 1
       climofiles(ii) = \
           interface_get_inpaths(ii) + "/" \
           + interface_get_infile(variables(0), field_types(0), ii)
   end do

   alltags = array_append_record(tags, (/"DM_reg", "ST_clim", "PT_clim"/), 0)

   if (multimean) then
       allmodelnames = array_append_record(models@name, (/"multi-model-mean"/), 0)
   end if

   caption = "Vertical profiles of variable " + variables(0) \
         + " with given latitude and month."

   id = diag_script + "_" + variables(0)

   contrib_authors = "A_wenz_sa"

   ESMValMD(output_dir+outfile, alltags, caption, id, variables(0), allmodelnames, climofiles, diag_script, contrib_authors)
;**************************************************************

   info_output(" gv " + wks@fullname + " &", verbosity, 1)

 end
