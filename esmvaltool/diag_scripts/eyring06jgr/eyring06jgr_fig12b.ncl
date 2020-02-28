;;#############################################################################
;; EYRING06JGR_FIG12B
;; Author:  Irene Cionni (ENEA, Italy) , Veronika Eyring (DLR, Germany)
;; PROJECT-NAME CCMI
 ;;#############################################################################
;; Description
;; Time series of monthly mean chemical tracer at selected latitude 
;; Add more bullets if required ;;
;; Required diag_script_info attributes (diagnostics specific)
;;
;;  @e06fig12b_latmin:  min lat  where variable is averaged
;;  @e06fig12b_latmax:  max lat where variable is averaged
;;  @e06fig12b_month:    season  when variable is averaged
;;(i.e. "1","2","3","10")
;;  @e06fig12b_multimean= (i.e. False/True)
;; Optional diag_script_info attributes (diagnostic specific)
;;  @e06fig12b_XMin : X axis min
;;  @e06fig12b_XMax : X axis max
;;  @e06fig12b_YMin : Y Axis min
;;  @e06fig12b_YMax : Y Axis max
;;  e06fig12b_levstring : string for the title about the level selected in the preprocessing


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
           iref=0 
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

   log_info("++++++++++++++++++++++++++++++++++++++++++")
   log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
   log_info("++++++++++++++++++++++++++++++++++++++++++")
  
   ; Set default values for optional diag_script_info attributes

   set_default_att(diag_script_info, "e06fig12b_month", "ANN")
   set_default_att(diag_script_info, "styleset", "CMIP5")
   set_default_att(diag_script_info, "e06fig12b_latmin", -90.)
   set_default_att(diag_script_info, "e06fig12b_latmax", 90.)
   set_default_att(diag_script_info, "e06fig12b_multimean", False)
   set_default_att(diag_script_info, "e06fig12b_xMin",1960)
   set_default_att(diag_script_info, "e06fig12b_xMax",2001)
   set_default_att(diag_script_info, "e06fig12b_yMin",0)
   set_default_att(diag_script_info, "e06fig12b_yMax",4)     
   set_default_att(diag_script_info, "e06fig12b_levstring"," ")
   seasons=diag_script_info@e06fig12b_month
   month_names=(/"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"/)
   seasons_name=month_names(toint(seasons)-1)

   latrange=new((/dimsizes(seasons),2/),"float")
   latrange(:,0)=diag_script_info@e06fig12b_latmin
   latrange(:,1)=diag_script_info@e06fig12b_latmax

   lev0=diag_script_info@e06fig12b_levstring
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

   XMax = diag_script_info@e06fig12b_xMax
   XMin = diag_script_info@e06fig12b_xMin
   YMax = diag_script_info@e06fig12b_yMax
   YMin = diag_script_info@e06fig12b_yMin

   multimean = diag_script_info@e06fig12b_multimean
   if (multimean)
      dimTot=dim_MOD+dim_OBS+1
   else
      dimTot=dim_MOD+dim_OBS
   end if


  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

end
begin
  aux_title_info = ""
  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring06f12b_" + var0)
  output_filename = wks@fullname
  nPanel = dimsizes(seasons)
  plot = new(nPanel, "graphic")
  plot_num = integertochar(ispan(97, 122, 1))
  rad = 4.0 * atan(1.0) / 180.0
  min_y1=min(y1)
  max_y2=max(y2)
  years=ispan(min_y1,max_y2,1)
  years!0="years"
  years&years=years  
  nyears=dimsizes(years)
  tmpData = new((/dim_MOD1, nPanel,nyears,2/), "float")
  tmpData!0 = "line"
  tmpData!1 = "season"
  tmpData!2 = "years" 
  tmpData!3 = "stat"
  tmpData&years=years

  Xdata = new((/dimTot, nPanel,nyears,2/), "float")
  Xdata!0 = "line"
  Xdata!1 = "season"
  Xdata!2 = "years"
  Xdata!3 = "stat"
  Xdata&years=years
  Xref= new((/ nPanel,nyears,2/), "float")
 
  Xref!0 = "season"
  Xref!1 = "years"
  Xref!2 = "stat"
  Xref&years=years
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
        
      if(all(isdim(A0,(/"time","lat","lon"/)))) then        

           dummy1 = area_operations(A0(time|:, lat|:, lon|:), latrange(ipan, 0), \
                               latrange(ipan, 1), 0., 360., "average", True)
           coslat   = cos(A0&lat * rad)
           coslat!0 = "lat"
           coslat&lat = A0&lat
           delete(coslat)
           delete(A0)
      else
          if(all(isdim(A0,(/"time","lat"/)))) then
             coslat   = cos(A0&lat * rad)
             coslat!0 = "lat"
             coslat&lat = A0&lat
             dummy1 = dim_avg_wgt_Wrap(A0(:,{latrange(ipan,0):latrange(ipan,1)}), coslat({latrange(ipan,0):latrange(ipan,1)}), 1)
             delete(coslat)
             delete(A0)
       
          end if
      end if
      if (seasons(ipan).ne."ANN") then
              dummy2=time_operations(dummy1,y1(imod),y2(imod),"extract",seasons(ipan),True)
              ;3 year smoothing
              delete(dummy1)
      else
                   
              dummy2=month_to_annual( dummy1, 1)          
      end if
      tmpData(imod,ipan,{y1(imod):y2(imod)},0) = (/dummy2/)
      delete(dummy2)
   end do ; end seasons loop
   Xdata(imod,:,:,:)= (/tmpData(imod,:,:,:)/) 
;**********UNITS CHANGE TO PPMV***********************************
   if ((UNITS.eq."10e6").or.(UNITS.eq."ppmv")) then
        Xdata=(/Xdata*10^3/)
        tmpData=(/tmpData*10^3/)
        Xref=(/Xref*10^3/)
        UNITS="ppbv"
        ;print("Units changed")
   end if

;*****************************************************************
   plotMod(imod)=names(imod)
   plotColor(imod)=colors(imod)
   plotDash(imod)=dashes(imod)
   plotThicks(imod)=thicks(imod)
   plotScatters(imod)=markers(imod) 
   if (isatt(variable_info[0], "reference_dataset")) then
       if (imod.eq.iref) then
          Xref=tmpData(imod,:,:,:)
          plotDash(imod)=     plotScatters(imod)
          LineScatters(imod)="Markers"
       end if
   end if 
 end do ;end loop model

 plotScatters(index_mod)=0
 LineScatters(index_obs)="Markers"
 LineScatters(index_mod)="Lines"
 
 if (multimean) then
     Xdata(dimTot-1,:,:,0)=(/dim_avg_n(Xdata(index_mod,:,:,0),0)/)
     plotMod(dimTot-1)="MMM"
     plotColor(dimTot-1)="Black"
     plotDash(dimTot-1)=0
     LineScatters(dimTot-1)="Lines"
     plotScatters(dimTot-1)=0
     Xdata(dimTot-1,:,:,1)=(/dim_stddev_n(Xdata(index_mod,:,:,0),0)/)
     X_MMM_std=(/dim_stddev_n(Xdata(index_mod,:,:,0),0)/) 
 end if

 plot_d=new(dimsizes(seasons),"graphic")
 plot_xp_R_d=new(dimsizes(seasons),"graphic")  
 plot_xp_M_d=new(dimsizes(seasons),"graphic") 
  
 
 xp_R= new((/2*dimsizes(Xdata&years),dimsizes(seasons)/),typeof(Xdata))
 yp_R= new((/2*dimsizes(Xdata&years),dimsizes(seasons)/),typeof(Xdata))

 xp_M= new((/2*dimsizes(Xdata&years),dimsizes(seasons)/),typeof(Xdata))
 yp_M= new((/2*dimsizes(Xdata&years),dimsizes(seasons)/),typeof(Xdata)) 
 if (multimean) then 
     do j=0,dimsizes(years)-1        
        yp_M(j,:) = Xdata(dimTot-1,:,j,0)-Xdata(dimTot-1,:,j,1)
        xp_M(j,:) = tofloat(years(j))
        xp_M(dimsizes(years)+j,:) = tofloat(years(dimsizes(years)-1-j))
        yp_M(dimsizes(years)+j,:) = Xdata(dimTot-1,:,(dimsizes(years)-1-j),0)+Xdata(dimTot-1,:,(dimsizes(years)-1-j),1)
     end do
  end if
 if (isatt(variable_info[0], "reference_dataset")) then
    do j=0,dimsizes(years)-1
       yp_R(j,:) = Xref(:,j,0)-Xref(:,j,1)
       xp_R(j,:) = tofloat(years(j))
       xp_R(dimsizes(years)+j,:) = tofloat(years(dimsizes(years)-1-j))
       yp_R(dimsizes(years)+j,:) = Xref(:,(dimsizes(years)-1-j),0)+Xref(:,(dimsizes(years)-1-j),1)
    end do
 end if
;*************************plotting*******************************************************

  do ipanel=0, dimsizes(seasons)-1
      res=True
      res@gsnDraw           = False
      res@gsnFrame          = False
      res@tfPolyDrawOrder   = "Predraw"
      
      res@trXMaxF           = XMax
      res@trXMinF           = XMin
      res@trYMaxF           = YMax(ipanel)
      res@trYMinF           = YMin(ipanel)
      res@tmLabelAutoStride = True
      res@vpHeightF         = 0.6                    ; change aspect ratio of plot
      res@vpWidthF          = 0.6
      res@tiYAxisString     = var0+" ("+UNITS+")"     ; y-axis label
      res@tiMainFontHeightF=  0.023 
      res@tiMainString      = seasons_name(ipanel)+" "+lat_names(latrange(ipanel,0))+"-"+lat_names(latrange(ipanel,1))+" at "+lev0
      res@xyMarkLineModes    = LineScatters(:)
      res@xyLineThicknessF= 4.0
      res@xyLineDashSegLenF=0.13
      res@xyLineColors      = plotColor(:)                 ; line color
      res@xyDashPatterns    = plotDash(:)
      res@xyMarkers         =plotScatters(:)
      res@xyMarkerSizeF     =0.01
      res@xyMarkerThicknessF =1.8
      plot_d(ipanel) = gsn_csm_xy (wks,years,Xdata(line|:,season|ipanel,years|:,stat|0),res)

      gsres = True
      if (isatt(variable_info[0], "reference_dataset")) then 
         gsres@gsFillColor = "Snow3"
        plot_xp_R_d(ipanel) = gsn_add_polygon(wks,plot_d(ipanel),xp_R(:,ipanel),yp_R(:,ipanel),gsres)

      end if
      if (multimean) then
         gsres@gsFillColor       = "cyan"                 ; color chosen
         gsres@gsFillOpacityF = 0.2  
         plot_xp_M_d(ipanel) = gsn_add_polygon(wks,plot_d(ipanel),xp_M(:,ipanel),yp_M(:,ipanel),gsres)
      end if
     

end do


   ;; Call legend-plotting routine
  Pres = True
  Pres@gsnFrame       = False
  Pres@gsnPanelBottom = 0.2
  if (dimsizes(seasons).gt.1) then
      gsn_panel(wks,plot_d,(/dimsizes(seasons),1/),Pres)
  else

      gsn_panel(wks,plot_d,(/1,1/),Pres)
  end if
        frame(wks)
  plotDash(index_obs)=plotScatters(index_obs)

  plotThicks=plotThicks+2
  add_legenda(wks,plotMod,plotColor,tofloat(plotDash),LineScatters,plotThicks,"0.2")

  frame(wks)
    ;**********************common output file**************************
     nc_filename = work_dir + "eyring06f12b_" + var0 + ".nc"
     ;Xdata(line|:,season|ipanel,years|:,stat|0)
     Xdata@var = var0
     Xdata@diag_script="eyring06f12b"
     Xdata&stat=(/"mean","std"/)
     Xdata&season=seasons_name
     Xdata&line=plotMod
     Xdata&years=years
     nc_outfile = ncdf_write(Xdata, nc_filename)

end

