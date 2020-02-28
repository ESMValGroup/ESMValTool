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
;;    @e06fig12a_start_year: year when start the climatology calculation
;;    @e06fig12a_start_year: year when end the climatology calculation
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

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/style.ncl"

begin
   enter_msg(DIAG_SCRIPT, "")
   var0 = variable_info[0]@short_name
   info0 = select_metadata_by_name(input_file_info, var0)
   dim_MOD1 = ListCount(info0)
   names = metadata_att_as_array(info0, "dataset")
   projects = metadata_att_as_array(info0, "project")
   exps = metadata_att_as_array(info0, "exp")
   ensembles = metadata_att_as_array(info0, "ensemble")
   y1 = metadata_att_as_array(info0, "start_year")
   y2 = metadata_att_as_array(info0, "end_year") 
   if (isatt(variable_info[0], "reference_dataset")) then
     refname = variable_info[0]@reference_dataset
     do imod = 0, dimsizes(names) - 1
         if (names(imod) .eq. refname) then
            iref = imod
          end if
      end do       
   else
     error_msg("f", DIAG_SCRIPT, "", "no reference dataset  " + \
              "(variable_info[0]@reference_dataset) needs to be defined.")
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
   print(LONG_NAME)

   log_info("++++++++++++++++++++++++++++++++++++++++++")
   log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
   log_info("++++++++++++++++++++++++++++++++++++++++++")
  
   ; Set default values for optional diag_script_info attributes

   set_default_att(diag_script_info, "e06fig01_season", "ANN")
   set_default_att(diag_script_info, "styleset", "CMIP5")
   set_default_att(diag_script_info, "e06fig01_latmin", -90.)
   set_default_att(diag_script_info, "e06fig01_latmax", 90.)
   set_default_att(diag_script_info, "e06fig01_multimean", False)
   set_default_att(diag_script_info, "e06fig01_start_year", max(y1))
   set_default_att(diag_script_info, "e06fig01_end_year",min(y2))
   
   seasons=diag_script_info@e06fig01_season
   latrange=new((/dimsizes(seasons),2/),"float")
   latrange(:,0)=diag_script_info@e06fig01_latmin
   latrange(:,1)=diag_script_info@e06fig01_latmax
   seasons=diag_script_info@e06fig01_season


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

   Max_lev=new((/dimsizes(seasons)/),"float")
   Min_lev=new((/dimsizes(seasons)/),"float")
   XMax=new((/dimsizes(seasons)/),"float")
   XMin=new((/dimsizes(seasons)/),"float")   

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

   multimean = diag_script_info@e06fig01_multimean
   if (multimean)
       dimTot=dim_MOD+dim_OBS
   else
       dimTot=dim_MOD+dim_OBS-1
   end if

   year_start=toint(diag_script_info@e06fig01_start_year)
   year_end=toint(diag_script_info@e06fig01_end_year)

  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

end
begin
  aux_title_info = ""
  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring06f01_" + var0)
  output_filename = wks@fullname
  nPanel = dimsizes(seasons)
  plot = new(nPanel, "graphic")
  plot_num = integertochar(ispan(97, 122, 1))
  rad = 4.0 * atan(1.0) / 180.0
  A0 = read_data(info0[iref]) 
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


;******************************************************************
  do imod = 0, dim_MOD1 - 1
    do ipan=0,dimsizes(seasons)-1
      ; Extract data
      A0 = read_data(info0[imod])
      if(all(isdim(A0,(/"time","plev","lat","lon"/)))) then        

        dummy1 = area_operations(A0(time|:,plev|:, lat|:, lon|:), latrange(ipan, 0), \
                               latrange(ipan, 1), 0., 360., "average", True)
           coslat   = cos(A0&lat * rad)
           coslat!0 = "lat"
           coslat&lat = A0&lat
           delete(coslat)
           delete(A0)
       else
          if(all(isdim(A0,(/"time","plev","lat"/)))) then
             coslat   = cos(A0&lat * rad)
             coslat!0 = "lat"
             coslat&lat = A0&lat

             dummy1 = dim_avg_wgt_Wrap(A0(:,:,{latrange(ipan,0):latrange(ipan,1)}), coslat({latrange(ipan,0):latrange(ipan,1)}), 1)
             delete(coslat)
             delete(A0)
       
          end if
       end if
       if (seasons(ipan).ne."ANN") then
                 tmp = time_operations(dummy1,year_start,year_end,"extract",seasons(ipan),True)
                 dummy2 = time_operations(tmp,year_start,year_end,"average","yearly",True)
                 delete(tmp)
                 delete(dummy1)
       else
                 dummy2 = time_operations(dummy1,year_start,year_end,"average","yearly",True)
       end if
       tmp = dim_avg_n_Wrap(dummy2,0)
       tmp_st = dim_stddev_n_Wrap(dummy2,0)
       delete(dummy2)
       tmpData(imod,ipan,:,0) = (/tmp/)
       tmpData(imod,ipan,:,1) = (/tmp_st/)
       delete(tmp)
       delete(tmp_st)
     end do ; end seasons loop
     if (imod.ne.iref) then
          Xdata(ik,:,:,:)= (/tmpData(imod,:,:,:)/) 
          plotMod(ik)=names(imod)
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
    if (multimean) then
       Xdata(dimTot-1,:,:,:)=(/dim_avg_n(Xdata(index_mod,:,:,:),0)/)
       plotMod(dimTot-1)="MMM"
       plotColor(dimTot-1)="Black"
       plotDash(dimTot-1)=0
       LineScatters(dimTot-1)="Lines"
   end if
   X_MMM_std=(/dim_stddev_n(Xdata(index_mod,:,:,:),0)/)
      
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
   if (multimean) then 
       do j=0,dimsizes(lev0)-1
    ;        printVarSummary(X_MMM_std)
            xp_M(j,:) = Xdata(dimTot-1,:,j,0)-((X_MMM_std(:,j,0)))
            yp_M(j,:) = tofloat(lev0(j)/100.)
            yp_M(dimsizes(lev0)+j,:) = tofloat(lev0(dimsizes(lev0)-1-j)/100.)
            xp_M(dimsizes(lev0)+j,:) = Xdata(dimTot-1,:,(dimsizes(lev0)-1-j),0)+(X_MMM_std(:,(dimsizes(lev0)-1-j),0))
         end do
    end if
;*************************plotting*******************************************************

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

     ; res@vpHeightF         = 0.4                    ; change aspect ratio of plot
     ; res@vpWidthF          = 1./dimsizes(seasons)

      res@tiYAxisString     = "Pressure (hPa)"       ; y-axis label
      res@tiXAxisString     = "Model "+LONG_NAME+" bias "+"("+UNITS+")"     ; y-axis label

      res@tiMainString      = seasons(ipanel)+" "+lat_names(latrange(ipanel,0))+"- "+lat_names(latrange(ipanel,1))
      res@xyMarkLineModes    = LineScatters

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


      Xdata&plev=lev0
      Ydata!0="plev"
      Ydata&plev=lev0
      plot_d(ipanel) = gsn_csm_xy (wks,Xdata(line|:,season|ipanel,{plev|80000:},stat|0),Ydata({plev|80000:}),res)

      gsres = True

      gsres@gsFillColor = "Snow3"
      plot_xp_R_d(ipanel) = gsn_add_polygon(wks,plot_d(ipanel),xp_R(:,ipanel),yp_R(:,ipanel),gsres)
      if (multimean) then
        gsres1=True
        gsres1@gsFillColor       = "cyan"                 ; color chosen
        gsres1@gsFillOpacityF = 0.2
        plot_xp_M_d(ipanel) = gsn_add_polygon(wks,plot_d(ipanel),xp_M(:,ipanel),yp_M(:,ipanel),gsres1)


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
   if (dimsizes(seasons).gt.1) then
      gsn_panel(wks,plot_d,(/2,dimsizes(seasons)/2/),Pres)
   else

      gsn_panel(wks,plot_d,(/1,1/),Pres)
   end if
     print(plotDash)
   
    if (dimsizes (plotMod).gt.21) then
        frame(wks)
    end if
     
       add_legenda(wks,plotMod,plotColor,tofloat(plotDash),LineScatters,plotThicks,"0.2")

   frame(wks)
;**********************common output file**************************
     nc_filename = work_dir + "eyring06f01_" + var0 + ".nc"
     Xdata@var = var0
     Xdata@diag_script="eyring06f01"

     printVarSummary(Xdata)
     Xdata&stat=(/"mean","stddev"/)
     Xdata&season=seasons
     Xdata&line=plotMod
     printVarSummary(Xdata) 
     nc_outfile = ncdf_write(Xdata, nc_filename)


end

