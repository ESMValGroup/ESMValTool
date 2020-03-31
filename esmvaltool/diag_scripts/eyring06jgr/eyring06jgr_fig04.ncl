;;#############################################################################
;; EYRING06JGR_FIG04
;; Author:  Irene Cionni (ENEA, Italy) , Veronika Eyring (DLR, Germany)
;; PROJECT-NAME CCMI
 ;;#############################################################################
;; Description
;; Time series of monthly mean anomalies respect 
;; Additional description of the diagnostic
;; cliamtological mean for selected month and latitudinal region
;; Add more bullets if required ;;
;; Required diag_script_info attributes (diagnostics specific)
;;
;;  @e06fig04_latmin:  min lat  where variable is averaged
;;  @e06fig04_latmax:  max lat where variable is averaged
;;  (i.e. (/60,90/))
;;  @e06fig04_season:    season  when variable is averaged
;;(i.e. "DJF","MAM","JJA","SON")
;;  @e06fig04_multimean= (i.e. False/True)
;;  e06fig04_reference_yearmax: min year for clculating the mena for the anomalies 
;;  e06fig04_reference_yearmax: max year for calculating  the mena for the anomalies
;;  @e06fig04_show_trend= (i.e. True/False) ; calculate the trend and write the value in the legend
;;  @e06fig04_trend_yearmin:min year for the trend
;;  @e06fig04_trend_yearmax: max year for the trend
;; Optional diag_script_info attributes (diagnostic specific)
;;  @e06fig04_XMin : X axis min
;;  @e06fig04_XMax : X axis max
;;  @e06fig04_YMin : Y Axis min
;;  @e06fig04_YMax : Y Axis max

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

   set_default_att(diag_script_info, "e06fig04_season", "ANN")
   set_default_att(diag_script_info, "styleset", "CMIP5")
   set_default_att(diag_script_info, "e06fig04_latmin", -90.)
   set_default_att(diag_script_info, "e06fig04_latmax", 90.)
   set_default_att(diag_script_info, "e06fig04_multimean", False)
   set_default_att(diag_script_info, "e06fig04_reference_yearmin",max(y1))
   set_default_att(diag_script_info, "e06fig04_reference_yearmax",min(y2))
   set_default_att(diag_script_info, "e06fig04_show_trend","False")   
   set_default_att(diag_script_info, "e06fig04_xMin",1960)
   set_default_att(diag_script_info, "e06fig04_xMax",2001)
   set_default_att(diag_script_info, "e06fig04_yMin",-20)
   set_default_att(diag_script_info, "e06fig04_yMax",20)     
   year_ref_max=toint(diag_script_info@e06fig04_reference_yearmax)
   year_ref_min=toint(diag_script_info@e06fig04_reference_yearmin)

   seasons=diag_script_info@e06fig04_season
   latrange=new((/dimsizes(seasons),2/),"float")
   latrange(:,0)=diag_script_info@e06fig04_latmin
   latrange(:,1)=diag_script_info@e06fig04_latmax
   seasons=diag_script_info@e06fig04_season
   trend = diag_script_info@e06fig04_show_trend

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

   XMax = diag_script_info@e06fig04_xMax
   XMin = diag_script_info@e06fig04_xMin
   YMax = diag_script_info@e06fig04_yMax
   YMin = diag_script_info@e06fig04_yMin

   multimean = diag_script_info@e06fig04_multimean
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
  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring06f04_" + var0)
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
    print(imod)
    do ipan=0,dimsizes(seasons)-1
      print(ipan)
      print(latrange(ipan, 0))
      print(latrange(ipan, 1))
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

              dummy2a=month_to_season (dummy1, seasons(ipan))
              ;3 year smoothing
              dummy2=runave (dummy2a, 3, 0)
              dummy2_avg = time_operations(dummy1,year_ref_min,year_ref_max,"average",seasons(ipan),True)
              delete(dummy2a) 
              delete(dummy1)
      else
                   
             ;dummy2 = time_operations(dummy1,y1(imod),y2(imod),"extract","yearly",True)
              dummy2=month_to_annual( dummy1, 1)          
              dummy2_avg = time_operations(dummy1,year_ref_min,year_ref_max,"average","annualclim",True) 
      end if
      tmp = dummy2
      ano_tmp=(/dummy2-dummy2_avg/)
      delete(dummy2)
      delete(dummy2_avg)
      tmpData(imod,ipan,{y1(imod):y2(imod)},0) = (/ano_tmp/)
      tmpData(imod,ipan,{y1(imod):y2(imod)},1) = (/tmp/)
      delete(tmp)
      delete(ano_tmp)
     end do ; end seasons loop
          Xdata(imod,:,:,:)= (/tmpData(imod,:,:,:)/) 
          plotMod(imod)=names(imod)
          plotColor(imod)=colors(imod)
          plotDash(imod)=dashes(imod)
          plotThicks(imod)=thicks(imod)
   
 end do ;end loop model

 plotScatters(index_mod)=0
 plotScatters(index_obs)=1
 LineScatters(index_obs)="Markers"
 LineScatters(index_mod)="Lines"

 if (multimean) then
     Xdata(dimTot-1,:,:,:)=(/dim_avg_n(Xdata(index_mod,:,:,:),0)/)
     plotMod(dimTot-1)="MMM"
     plotColor(dimTot-1)="Black"
     plotDash(dimTot-1)=0
     LineScatters(dimTot-1)="Lines"
 end if
 X_MMM_std=(/dim_stddev_n(Xdata(index_mod,:,:,:),0)/)
 if (trend.eq."True") then
     year_Min_trend = diag_script_info@e06fig04_trend_yearmin
     year_Max_trend = diag_script_info@e06fig04_trend_yearmax 
     trend_Xdata=regCoef(years({year_Max_trend:year_Min_trend}),Xdata(line|:,season|:,{years|year_Max_trend:year_Min_trend},stat|0))
     plotMod1=plotMod+" trend"
 end if        

  plot_d=new(dimsizes(seasons),"graphic")
  plot_xp_R_d=new(dimsizes(seasons),"graphic")  
  plot_xp_M_d=new(dimsizes(seasons),"graphic") 
  
  xp_M= new((/2*dimsizes(Xdata&years),dimsizes(seasons)/),typeof(Xdata))
  yp_M= new((/2*dimsizes(Xdata&years),dimsizes(seasons)/),typeof(Xdata)) 
  if (multimean.eq."True") then 
     do j=0,dimsizes(years)-1        
        yp_M(j,:) = Xdata(dimTot-1,:,j,0)-((X_MMM_std(:,j,0)))
        xp_M(j,:) = tofloat(years(j))
        xp_M(dimsizes(years)+j,:) = tofloat(years(dimsizes(years)-1-j))
        yp_M(dimsizes(years)+j,:) = Xdata(dimTot-1,:,(dimsizes(years)-1-j),0)+(X_MMM_std(:,(dimsizes(years)-1-j),0))
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
      res@vpHeightF         = 0.3                    ; change aspect ratio of plot
      res@vpWidthF          = 0.9
      res@tiYAxisString     = "Anomaly "+"("+UNITS+")"     ; y-axis label
      res@tiMainFontHeightF=  0.023 
      res@tiMainString      = seasons(ipanel)+" "+lat_names(latrange(ipanel,0))+"-"+lat_names(latrange(ipanel,1))
      res@xyMarkLineModes    = LineScatters
      res@xyLineThicknessF= 3.0
      res@xyLineDashSegLenF=0.13
      res@xyLineColors      = plotColor                 ; line color
      res@xyDashPatterns    = plotDash
      res@xyMarkers         =plotDash 
      res@xyMarkerSizeF     =0.004
      res@xyMarkerThicknessF =1.8
      if (trend.eq."True") then
         plotMod1=plotMod+" =  "+sprintf("%5.2f",10*trend_Xdata(:,ipanel))
         res@pmLegendDisplayMode    = "Always"            ; turn on legend
         res@pmLegendSide           = "Top"               ; Change location of
         res@pmLegendParallelPosF   = 1.3                ; move units right
         res@pmLegendOrthogonalPosF = -1.2                ; move units down
         res@pmLegendWidthF         = 0.16                ; Change width and
         res@pmLegendHeightF        = 0.3                ; height of legend.
         res@lgLabelFontHeightF     = .018                 ; change font height
         res@lgPerimOn              = False 
         res@xyExplicitLegendLabels = plotMod1         ; explicit labels
         delete(plotMod1)
      end if
      plot_d(ipanel) = gsn_csm_xy (wks,years,Xdata(line|:,season|ipanel,years|:,stat|0),res)
      add_line(wks,plot_d(ipanel),(/years(0),years(dimsizes(years)-1)/),(/0.,0./),"gray",1) 
      gsres = True
    
      gsres@gsFillColor = "Snow3"
     
      if (multimean.eq."True") then
         gsres@gsFillColor       = "cyan"                 ; color chosen
   
         plot_xp_M_d(ipanel) = gsn_add_polygon(wks,plot_d(ipanel),xp_M(:,ipanel),yp_M(:,ipanel),gsres)
         print(xp_M(:,ipanel))
         print(yp_M(:,ipanel))
        
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
      gsn_panel(wks,plot_d,(/dimsizes(seasons),1/),Pres)
   else

      gsn_panel(wks,plot_d,(/1,1/),Pres)
   end if
     print(plotDash)
   
    if (dimsizes (plotMod).gt.21) then
        frame(wks)
    end if
     plotThicks=plotThicks+1
     add_legend(wks,plotMod,plotColor,tofloat(plotDash),LineScatters,plotThicks,"0.2")

     frame(wks)
end

