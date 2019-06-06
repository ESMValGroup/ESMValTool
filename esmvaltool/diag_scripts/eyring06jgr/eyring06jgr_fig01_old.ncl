;;#############################################################################
;; EYRING06JGR_FIG01
;; Author: Irene Cionni (ENEA, Italy), Veronika Eyring (DLR, Germany)
;; PROJECT-NAME CCMI
;;#############################################################################
;; Description
;;    Vertical profile climatological mean bias,
;;    climatological mean for selected seasons and latitudinal region
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    @e06fig01_latrange:  min lat and max lat where variable is averaged
;;                         (i.e. (/60,90/))
;;    @e06fig01_season:    season when variable is averaged
;;                         (i.e. "DJF", "MAM", "JJA", "SON")
;;    @e06fig01_multimean: calculate multi-model mean (i.e. False/True)
;;    @e06fig01_refModel:  name of the reference model or observation for
;;                         bias calculation (e.g. (/"ERA-Interim"/))
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    @e06fig01_XMax:     max limit X axis
;;    @e06fig01_XMin:     min limit X axis
;;    @e06fig01_levmax:   max limit Y axis
;;    @e06fig01_levmin:   min limit Y axis
;;
;; Required variable_info attributes (variable specific)
;;
;; Caveats
;;
;; Modification history
;;    20180320-A_laue_ax: added tags for reporting and call to write_references
;;    20170315-A_cion_ir: written
;;#############################################################################
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/style.ncl"

load "$diag_scripts/shared/plot/misc_function.ncl"


begin
   enter_msg(DIAG_SCRIPT, "")
   var0 = variable_info[0]@short_name
   info0 = select_metadata_by_name(input_file_info, var0)
   dim_MOD1 = ListCount(info0)
   if (isatt(variable_info[0], "reference_dataset")) then
     refname = variable_info[0]@reference_dataset
   end if
   names = metadata_att_as_array(info0, "dataset")
   projects = metadata_att_as_array(info0, "project")
   exps = metadata_att_as_array(info0, "exp")
   ensembles = metadata_att_as_array(info0, "ensemble")

   log_info("++++++++++++++++++++++++++++++++++++++++++")
   log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
   log_info("++++++++++++++++++++++++++++++++++++++++++")
  
   ; Set default values for optional diag_script_info attributes

   set_default_att(diag_script_info, "e13fig01_season", "ANN")
   set_default_att(diag_script_info, "styleset", "CMIP5")
   
   if (.not.isatt(diag_script_info, "e06fig01_latrange")) then
      error_msg("f", DIAG_SCRIPT, "", "no latitude bands " + \
              "(e06fig01_latrange) defined.")
   end if  
    
   colors = project_style(info0, diag_script_info, "colors")
   dashes = project_style(info0, diag_script_info, "dashes")
   thicks = tofloat(project_style(info0, diag_script_info, "thicks"))
   markers = project_style(info0, diag_script_info, "markers") 


  flag_mod = where(projects.ne."OBS" .and. projects.ne."obs4mips", 1, 0)
  index_mod = ind(flag_mod.gt.0)
  index_obs = ind(flag_mod.eq.0)
  dim_MOD = dimsizes(index_mod)
  dim_OBS = 0

  if (.not.all(ismissing(index_obs))) then
    dim_OBS = dimsizes(index_obs)
    list_obs_mod = names(index_obs)
  end if

  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

end
;begin
;    vardeffile = "interface_data/" + var0 + "_info.tmp"
;    loadscript(vardeffile)
;end
begin
  aux_title_info = ""
  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring06f01_" + var0)
  output_filename = wks@fullname

  plot = new(dimsizes(diag_script_info@e06fig01_season), "graphic")
  plot_num = integertochar(ispan(97, 122, 1))
;******************************************************************
 do imod = 0, dim_MOD1 - 1
      ; Extract data
      A0 = read_data(info0[imod])
      printVarSummary(A0)
end do
end 
      dummy1 = area_operations(A0(time|:, lat|:, lon|:), latrange(ipan, 0), \
                               latrange(ipan, 1), 0., 360., "average", True)

    A = extract_data(0, data, -1, 0, 0)
    rank = dimsizes(dimsizes(A))
    delete(A) 
    if(isvar("clima")) then
       delete(clima)
    end if 
    if(isvar("clima_files")) then
       delete(clima_files)
    end if    
    dim_OBS=0

   if ((isatt(diag_script_info,"e06fig01_climatology")).and. (isatt(diag_script_info,"e06fig01_climatology_files"))) then
       clima= diag_script_info@e06fig01_climatology
       clima_files=diag_script_info@e06fig01_climatology_files
       dim_OBS=dimsizes(clima)
       ;dim_OBS=0
   end if
 ;   iob=1
    if (isatt(diag_script_info, "e06fig01_refModel").or.(isatt(diag_script_info, "e06fig01_refObs"))) then
      ;print(diag_script_info@e06fig01_refModel)
        do imod = 0, dim_MOD - 1
            if (names(imod) .eq. diag_script_info@e06fig01_refModel) then
                iref = imod
                ref_mod=names(imod)
                print(iref)
            end if
        end do
     else
        print("fatal:plot requires diag_script_info attribute 'e06fig01_refModel' or 'e06fig01_refObs' ")
        status_exit(1)    
    end if

    if(rank .gt. 1) then
        if(isvar("latrange")) then
            delete(latrange)
        end if
        if(isatt(diag_script_info, "e06fig01_latrange")) then
            latrange = diag_script_info@e06fig01_latrange
        else
            print("fatal:plot requires diag_script_info attribute 'e06fig01_latrange'")
            status_exit(1)
        end if

        if(isvar("lonrange")) then
            delete(lonrange)
        end if
        if(isatt(diag_script_info, "e06fig01_season")) then
            seasons = diag_script_info@e06fig01_season
         else
            if (isStrSubset(field_type0, "z") .eq. False) then
                print("fatal:plot requires diag_script_info attribute 'e06fig01_season'")
                status_exit(1)
            end if
        end if
        if (dimsizes(seasons).ne.dimsizes(latrange(:,0))) then
               print("fatal:plot requires dimension 'e06fig01_season' equal last dimension 'e06fig01_latrange' ")
               status_exit(1)       
        end if

    end if  
    
    if(isatt(diag_script_info, "levmax")) then
        Max_lev = diag_script_info@e06fig01_levmax
    else
        Max_lev=1000.
    end if
    if(isatt(diag_script_info, "levmin")) then
        Min_lev = diag_script_info@e06fig01_levmin
    else
        Min_lev=1.
    end if
    if(isatt(diag_script_info, "XMax")) then
        XMax = diag_script_info@e06fig01_XMax
    else
        XMax=20.
    end if
    if(isatt(diag_script_info, "XMin")) then
        XMin = diag_script_info@e06fig01_XMin
    else
        XMin=-30. 
    end if
    multimean = "False"
    if(isatt(diag_script_info, "e06fig01_multimean")) then
        multimean = diag_script_info@e06fig01_multimean
    end if
    if (multimean.eq."True")
        dimTot=dim_MOD+dim_OBS
    else
        dimTot=dim_MOD+dim_OBS-1
    end if
    if(isatt(variable_info, "units")) then
        UNITS = variable_info@units
    else
        UNITS= "K"
    end if
   
   if(isatt(variable_info, "long_name")) then
        LONG_NAME = variable_info@long_name
    else
        LONG_NAME = var0
    end if
;    diag_script_base = basename(diag_script)
;    aux_title_info=str_join(seasons,"_")
;
;    output_filename = interface_get_figure_filename(diag_script_base,\
;                                              var0,\
;                                              field_type0,\
;                                              aux_title_info,\
;                                              iref)
     ;print(output_filename)
     
    nPanel = dimsizes(seasons)

    plot_num = (/"a", "b", "c", "d", "e", "f", "g", "h", "i", "l", \
                 "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "z"/)

    rad = 4.0 * atan(1.0) / 180.0

    A0 = extract_data(iref, data, -1, 0, 0)
    ;printVarSummary(A0)
    lev0 = A0&plev
    dimLev = dimsizes(lev0)
    delete(A0)
   
    tmpData = new((/dim_MOD+dim_OBS, nPanel,dimsizes(lev0),2/), "float")
    tmpData!0 = "line"
    tmpData!1 = "season"
    tmpData!2 = "plev" 
    tmpData!3 = "stat"


    Xdata = new((/dimTot, nPanel,dimsizes(lev0),2/), "float")
    Xdata!0 = "line"
    Xdata!1 = "season"
    Xdata!2 = "plev"
    Xdata!3 = "stat"
    Xdata&plev=lev0
    Xref= new((/ nPanel,dimsizes(lev0),2/), "float")
    plotMod=new(dimTot,typeof(names))
    plotColor=new(dimTot,typeof(names))
    plotDash=new(dimTot,typeof(names))
    plotThicks=new(dimTot,"float")
    plotScatters=new(dimTot,typeof(names))
    LineScatters=new(dimTot,typeof(names)) 
    ik=0

    do imod = 0, dim_MOD - 1

       A0 = extract_data(imod, data, -1, 0, 0)
       ;uniform_coordinates(A0, field_type0)
       plev_n=A0&plev  
       if (max(plev_n).gt.1000.) then
          plev_n=plev_n
       else
          plev_n=plev_n*100.
       end if  
       Aint=int2p_n_Wrap(plev_n,A0,lev0,2,1)
       delete(plev_n)    
       delete(A0)
       A0=Aint         
       delete(Aint)
       do iseas=0,dimsizes(seasons)-1
          if (rank .gt. 1) then

             if (field_type0 .eq. "T3M") then
               ;  dummy1=area_operations(A0(time|:,plev|:,lat|:,lon|:),latrange(iseas,0),latrange(iseas,1),0.,360.,"average",True)
                 dummy2=dim_avg_Wrap(A0)
                 coslat   = cos(A0&lat * rad)
                 coslat!0 = "lat"
                 coslat&lat = A0&lat
                 ;printVarSummary(A0)

                dummy1 = dim_avg_wgt_Wrap(dummy2(:,:,{latrange(iseas,0):latrange(iseas,1)}), coslat({latrange(iseas,0):latrange(iseas,1)}), 1)
                delete(dummy2)
                delete(coslat)
             else
             coslat   = cos(A0&lat * rad)
             coslat!0 = "lat"
             coslat&lat = A0&lat 
             ;printVarSummary(A0)
             
             dummy1 = dim_avg_wgt_Wrap(A0(:,:,{latrange(iseas,0):latrange(iseas,1)}), coslat({latrange(iseas,0):latrange(iseas,1)}), 1)
                ;dummy1=area_operations(A0(time|:,plev|:,lat|:),latrange(iseas,0),latrange(iseas,1),0.,360.,"average",True) ;not sure it's working
             delete(coslat)
             end if
             if (seasons(iseas).ne."ANN") then
                 tmp = time_operations(dummy1,toint(models@start_year(imod)),toint(models@end_year(imod)),"extract",seasons(iseas),True)
                 dummy2 = time_operations(tmp,toint(models@start_year(imod)),toint(models@end_year(imod)),"average","yearly",True)
                 delete(tmp)
                 delete(dummy1)
             else
                 dummy2 = time_operations(dummy1,toint(models@start_year(imod)),toint(models@end_year(imod)),"average","yearly",True)
             end if
             ;tmplt = dim_avg_Wrap(month_to_season(A0(time |:, plev | :, \
             ;                        {lat | latrange(iseas,0):latrange(iseas,1)}, \
             ;                        {lon | :}),seasons(iseas)))
             ; tmpz = dim_avg_wgt_Wrap(tmplt(:, :, :), coslat({latrange(iseas,0):latrange(iseas,1)}), 1)
            
             tmp = dim_avg_n_Wrap(dummy2,0)
             tmp_st = dim_stddev_n_Wrap(dummy2,0)
             delete(dummy2)
             tmpData(imod,iseas,:,0) = (/tmp/)
             tmpData(imod,iseas,:,1) = (/tmp_st/)
             delete(tmp)
             delete(tmp_st)
          end if  ; rank = 1 loop
       end do
       delete(A0)

       if (imod.ne.iref) then
          Xdata(ik,:,:,:)= (/tmpData(imod,:,:,:)/) 
          plotMod(ik)=models@name(imod)
          plotColor(ik)=colors(imod)
          plotDash(ik)=dashes(imod)
          plotThicks(ik)=thicks(imod)
          plotScatters(ik)=0.
          LineScatters(ik)="Lines"
          ik=ik+1
       else
          Xref=tmpData(imod,:,:,:)
      end if

    end do
;********************************************************
;print(dim_OBS)

  if (dim_OBS.ge.1) then   
     do iobs = 0, dim_OBS - 1
        ;print( clima_files(iobs))
        f_obs = addfile(clima_files(iobs),"r")
        A0 = f_obs->$var0$
        delete(f_obs)
        ;printVarSummary(A0)
        ;uniform_coordinates(A0, "T2Mz")
        plev_n=A0&plev
        if (max(A0&plev).gt.1000.) then
           plev_n=plev_n
        else
           plev_n=plev_n*100.
        end if     
        Aint=int2p_n_Wrap(plev_n,A0,lev0,2,1)
        delete(plev_n)    
        ;printVarSummary(Aint)
        delete(A0)
        A0=Aint 
        rank_l= dimsizes(dimsizes(A0))       
        ;print(rank_l)        
        delete(Aint)
        do iseas=0,dimsizes(seasons)-1
           ;print(rank_l)           
           ;print(iseas)
           if(rank_l .gt. 1) then
             ;printVarSummary(A0)
              ;print(A0&lat)
   
              coslat   = cos(A0&lat * rad)
              ;print(coslat)
               
              ;print(A0&lat)              
              coslat!0 = "lat"
              coslat&lat = A0&lat     
              
                 ;printVarSummary((A0(time |:, plev |:, \
                 ;                          {lat | latrange(iseas,0):latrange(iseas,1)})))
               
                 tmp = dim_avg_wgt_Wrap(month_to_season(A0(time |:, plev |:, \
                                           {lat | latrange(iseas,0):latrange(iseas,1)}),seasons(iseas)),coslat({latrange(iseas,0):latrange(iseas,1)}),1)
   ;  tmp = dim_avg_Wrap(month_to_season(A0(time |:, plev |:, \
   ;                                        {lat | latrange(iseas,0):latrange(iseas,1)}),seasons(iseas)))
            ;print(dimsizes(tmp))
 
             if(dimsizes(dimsizes(tmp)).gt.1) then
                 tmp1=dim_avg_n_Wrap(tmp,0)
              else
                 tmp1=tmp
              end if 
             ; printVarSummary(tmp)
             ; printVarSummary(tmpData)
            ;  print(dim_MOD+iobs)
               
              tmpData(dim_MOD+iobs,iseas,:,0) = (/tofloat(tmp1)/)
           end if  ; rank = 1 loop
        end do
        delete (tmp)
        delete(tmp1)
        delete(A0)
        delete(coslat)
        
        Xdata(dim_MOD-1+iobs,:,:,:)= tmpData(dim_MOD+iobs,:,:,:) 
        plotMod(dim_MOD-1+iobs)=clima(iobs)
        plotColor(dim_MOD-1+iobs)="Black"
        plotDash(dim_MOD-1+iobs)=iobs+2
        plotScatters(dim_MOD-1+iobs)=iobs+2        
        if(clima(iobs).eq."UKMO") then
          plotDash(dim_MOD-1+iobs)=2
          plotScatters(dim_MOD-1+iobs)=2
        end if 
        if(clima(iobs).eq."NCEP") then
          plotDash(dim_MOD-1+iobs)=3
          plotScatters(dim_MOD-1+iobs)=3
        end if        
        LineScatters(dim_MOD-1+iobs)="Markers"   
     end do
  end if

   if (multimean.eq."True") then
       Xdata(dimTot-1,:,:,:)=(/dim_avg_n(Xdata(0:dim_MOD-2,:,:,:),0)/)
       plotMod(dimTot-1)="MMM"
       plotColor(dimTot-1)="Black"
       plotDash(dimTot-1)=0
       LineScatters(dim_MOD-1+iobs)="Lines"
   end if
    X_MMM_std=(/dim_stddev_n(Xdata(0:dim_MOD-2,:,:,:),0)/)
      
   Xdata=(/Xdata-conform_dims(dimsizes(Xdata),Xref,(/1,2,3/))/)
  
   Ydata= lev0/100.

   plot_d=new(dimsizes(seasons),"graphic")
    plot_xp_R_d=new(dimsizes(seasons),"graphic")  
    plot_xp_M_d=new(dimsizes(seasons),"graphic") 
    xp_R= new((/2*dimsizes(lev0),dimsizes(seasons)/),typeof(Xref))
      yp_R=new((/2*dimsizes(lev0),dimsizes(seasons)/),typeof(Xref))

      xp_M= new((/2*dimsizes(lev0),dimsizes(seasons)/),typeof(Xref))
      yp_M= new((/2*dimsizes(lev0),dimsizes(seasons)/),typeof(Xref)) 
      do j=0,dimsizes(lev0)-1
      
         xp_R(j,:) = -((Xref(:,j,1)))
         yp_R(j,:) = tofloat(lev0(j)/100.)
         yp_R(dimsizes(lev0)+j,:) = tofloat(lev0(dimsizes(lev0)-1-j)/100)      
         xp_R(dimsizes(lev0)+j,:) = (Xref(:,(dimsizes(lev0)-1-j),1))
      end do
     if (multimean.eq."True") then 
        do j=0,dimsizes(lev0)-1
    ;        printVarSummary(X_MMM_std)
            xp_M(j,:) = Xdata(dimTot-1,:,j,0)-((X_MMM_std(:,j,0)))
            yp_M(j,:) = tofloat(lev0(j)/100.)
            yp_M(dimsizes(lev0)+j,:) = tofloat(lev0(dimsizes(lev0)-1-j)/100.)
            xp_M(dimsizes(lev0)+j,:) = Xdata(dimTot-1,:,j,0)+(X_MMM_std(:,(dimsizes(lev0)-1-j),0))
         end do
      end if
   ;print(file_type)
   ;print(output_dir)

   output_filename = output_dir + diag_script_base + "_" + var0 + "." + file_type
   wks = gsn_open_wks(file_type, output_filename)

  do ipanel=0, dimsizes(seasons)-1
      res                   = True
      res@gsnDraw           = False
      res@gsnFrame          = False
      res@tfPolyDrawOrder   = "Predraw"
      
      res@trXMaxF           = XMax
      res@trXMinF           = XMin
      res@trYMaxF           = Max_lev
      res@trYMinF           = Min_lev
      res@tmYLMode          = "Explicit"
      res@tmYLValues        = (/1,3,10,30,100,300,1000/)
      res@tmYLLabels        = (/1,3,10,30,100,300,1000/)
      res@tmLabelAutoStride = True

      res@vpHeightF         = 0.2                    ; change aspect ratio of plot
      res@vpWidthF          = 1./dimsizes(seasons)

      res@tiYAxisString     = "Pressure (hPa)"       ; y-axis label
      res@tiXAxisString     = "Model "+LONG_NAME+" bias "+"("+UNITS+")"     ; y-axis label
         
      res@tiMainString      = seasons(ipanel)+" "+lat_names(latrange(ipanel,0))+"- "+lat_names(latrange(ipanel,1))
      res@xyMarkLineModes    = LineScatters
      ;print(LineScatters)
     
      res@trYAxisType       = "LogAxis"
      res@trYReverse        = True
      res@xyLineThicknessF= 2.0
      res@xyLineDashSegLenF=0.13
      res@xyLineColors      = plotColor                 ; line color
      res@xyDashPatterns    = plotDash
      res@xyMarkers         =plotDash 
      res@xyMarkerSizeF     =0.004
      res@xyMarkerThicknessF =1.8
      res@tmXMajorGrid      = True                   ; implement x grid
      res@tmXMajorGridThicknessF = 1.0               ; 2.0 is default
      res@tmXMajorGridLineDashPattern = 2            ; select short dash lines


      ;printVarSummary(Xdata)
      Xdata&plev=lev0
      Ydata!0="plev"
      Ydata&plev=lev0
      plot_d(ipanel) = gsn_csm_xy (wks,Xdata(line|:,season|ipanel,{plev|80000:},stat|0),Ydata({plev|80000:}),res)

      gsres = True
    
      gsres@gsFillColor = "Snow3"
 ;gsres@gsFillOpacityF = 0.5
      plot_xp_R_d(ipanel) = gsn_add_polygon(wks,plot_d(ipanel),xp_R(:,ipanel),yp_R(:,ipanel),gsres)
      ;gsres@gsFillOpacityF = 0.5 
      if (multimean.eq."True") then
        gsres@gsFillColor       = "cyan"                 ; color chosen
        plot_xp_M_d(ipanel) = gsn_add_polygon(wks,plot_d(ipanel),xp_M(:,ipanel),yp_M(:,ipanel),gsres)
        
        
      end if
   end do
;; Call legend-plotting routine
   Pres = True
   Pres@gsnFrame       = False
   Pres@gsnPanelBottom = 0.2
   Pres@gsnPanelFigureStrings=plot_num(0:dimsizes(seasons)-1) 
   Pres@gsnPanelFigureStringsFontHeightF=0.025
   Pres@gsnPanelFigureStringsPerimOn="False"
   Pres@amJust="TopLeft"
;; Call panel routine
   gsn_panel(wks,plot_d,(/2,dimsizes(seasons)/2/),Pres)         
   add_legenda(wks,plotMod,plotColor,tofloat(plotDash),LineScatters,plotThicks,"0.2")

   frame(wks)

    ; *************************************************************************************************
    ; add meta data to plots (for reporting)

    caption = "Climatological mean biases of "+  variables(0) + " for (top) 60-90N and " \
              + "(bottom) 60-90S for the (left) winter and (right) spring seasons. " \
              + "Similar to Eyring et al. (2006), fig. 1."
    id = diag_script + "_" + variables(0)
    contrib_authors = "A_cion_ir"

    if (isatt(diag_script_info, "e06fig01_multimean")) then
        if (diag_script_info@e06fig01_multimean.eq."True") then
            modelsatname = array_append_record((/"multi-model-mean"/), \
                                               models@name, 0)
            delete(models@name)
            models@name = modelsatname
            delete(modelsatname)
        end if
    end if

    ESMValMD(output_filename, alltags, caption, id, variables(0), models@name, climofiles, \
             diag_script, contrib_authors)

    info_output("wrote " + output_filename, verbosity, 1)
end
