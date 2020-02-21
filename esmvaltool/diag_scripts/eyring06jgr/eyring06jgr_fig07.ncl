;EYRING06JGR_FIG07
;; Author: Irene Cionni (ENEA, Italy), Veronika Eyring (DLR, Germany)
;; PROJECT-NAME CCMI
;;#############################################################################
;; Description
;;   Climatological mean annual cycle calculated from
;;   the monthly mean temperature or H2O mixing ratio,
;;    climatological mean for selected levels and latitudinal region
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    @e06fig07_multimean: calculate multi-model mean (i.e. False/True)
;;    @e06fig07_refModel:  name of the reference model or observation for
;;                         bias calculation (e.g. (/"ERA-Interim"/))
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    @e06fig07_Ymax:   max limit Y axis
;;    @e06fig07_Ymin:   min limit Y axis
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
   print(names)
   projects = metadata_att_as_array(info0, "project")
   print(projects)
   exps = metadata_att_as_array(info0, "exp")
   print(exps)
   ensembles = metadata_att_as_array(info0, "ensemble")
   print("ensembles")
   print(metadata_att_as_array(info0, "start_year"))
   print(metadata_att_as_array(info0, "end_year"))
   y1 =(/ metadata_att_as_array(info0, "start_year")/)
   y2 =(/ metadata_att_as_array(info0, "end_year") /)
   log_info(y1)
   log_info(y2)
   

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

   print(refname)
   print(iref)
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
   print(UNITS)
   print(LONG_NAME)

   log_info("++++++++++++++++++++++++++++++++++++++++++")
   log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
   log_info("++++++++++++++++++++++++++++++++++++++++++")
  
   ; Set default values for optional diag_script_info attributes

   set_default_att(diag_script_info, "styleset", "CMIP5")
   set_default_att(diag_script_info, "e06fig07_Ymin", 185.)
   set_default_att(diag_script_info, "e06fig07_Ymax", 200.)
   set_default_att(diag_script_info, "e06fig07_multimean", False)
   
   set_default_att(diag_script_info, "e06fig07_lev", "100hPa")
   set_default_att(diag_script_info, "e06fig07_latrange", "10S-10N")

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

         
  YMax = diag_script_info@e06fig07_Ymax
  YMin = diag_script_info@e06fig07_Ymin
  sellev=diag_script_info@e06fig07_lev
  sellat=diag_script_info@e06fig07_latrange
    multimean = diag_script_info@e06fig07_multimean
    print(multimean)
    if (multimean)
        dimTot=dim_MOD+dim_OBS+1
    else
        dimTot=dim_MOD+dim_OBS
    end if

   print(dimTot)
   print(dim_MOD)
   print(dim_OBS)
      print(dim_MOD1)
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
  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring06f07_" + var0)
  output_filename = wks@fullname
  nPanel = 1
  plot = new(nPanel, "graphic")
  plot_num = integertochar(ispan(97, 122, 1))
  rad = 4.0 * atan(1.0) / 180.0
   
  tmpData = new((/dim_MOD+dim_OBS, nPanel,12,2/), "double")
  tmpData!0 = "line"
  tmpData!1 = "season"
  tmpData!2 = "month"
  tmpData!3 = "stat"

  XXdata = new((/dimTot, nPanel,12,2/), "double")
  XXdata!0 = "line"
  XXdata!1= "season"
  XXdata!2 = "month"
  XXdata!3 = "stat"

  XXref= new((/ nPanel,12,2/), "double")
  XXref!0= "season"
  XXref!1 = "month"
  XXdata!2 = "stat"



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
    do ipan=0,nPanel-1
      ; Extract data
      A0 = read_data(info0[imod])
  ;    Aref = read_data(info0[iref])
      dummy1 = time_operations(A0,y1(imod),y2(imod),"average","monthlyclim",True)
      dummy2_b=tofloat(A0)
      date= yyyymm_time(y1(imod),y2(imod),"integer") 
      yyyy  = date/100 
      mm  = date-yyyy*100                ; mmdd = yyyymmdd%10000
      weights = days_in_month(yyyy, mm)

      dummy2 = dummy2_b(0:11 )  ; Copy metadata
      do im = 0, 11
          dummy2(im ) = \
               dim_stddev_wgt_Wrap(dummy2_b(im::12), weights(im::12), 1)
      end do
     delete(A0)
     delete(weights)
     delete(date)
     delete(yyyy)
     delete(mm)
     delete(dummy2_b)
     dummy2!0=dummy1!0 
      ; a= clmMon2clmDayn(dummy1,0,0)
      ;printVarSummary(a)
     ; tmpData(imod,ipan,:,:,0) = (/ clmMon2clmDayn(dummy1,0,0)/)
       
      ;tmpData(imod,ipan,:,:,1) = (/clmMon2clmDayn(dummy2,0,0)/)
      ;print(tmpData(imod,ipan,:,:,0))
      printVarSummary(dummy1)
      printVarSummary(dummy2)
      tmpData(imod,ipan,:,0) = (/dummy1/)
      tmpData(imod,ipan,:,1) = (/dummy2/)
 
      delete(dummy1)
      delete(dummy2)
   end do ; end seasons loop
 end do
 ;X_MMM_std=(/dim_stddev_n(Xdata(index_mod,:,:,:),0)/)
 XXdata(:dim_MOD1-1,:,:,:)=(/tmpData/)

 plotMod(:dim_MOD1-1)=names
 plotColor(:dim_MOD1-1)=colors
 plotDash(:dim_MOD1-1)=dashes
 plotThicks(:dim_MOD1-1)=thicks
 plotScatters(:dim_MOD1-1)=0.
 LineScatters(index_mod)="Lines"
 LineScatters(index_obs)="Markers"
  XXref= (/tmpData(iref,:,:,:)/)
         
  XXdata!0 = "line"
  XXdata!1= "season"
  XXdata!2 = "month"
  XXdata!3 = "stat"
  printVarSummary(XXdata)

  ;°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
  if (multimean) then
             XXdata(dimTot-1,:,:,0)=(/dim_avg_n(XXdata(index_mod,:,:,0),0)/)
             plotMod(dimTot-1)="MMM"
             plotColor(dimTot-1)="Black"
             plotDash(dimTot-1)=0.
             LineScatters(dimTot-1)="Lines"
  end if
; print(Xdata(line|:,season|:,day|:,stat|0))

X_MMM_std=(/dim_stddev_n(XXdata(index_mod,:,:,:),0)/)
  YData=ispan(1,12,1)
  YData!0="month"
  YData&month=YData
  plot_xp_M_d=new((nPanel),"graphic") 
  
  xp_M= new((/2*dimsizes(YData),nPanel/),typeof(XXdata))
  yp_M= new((/2*dimsizes(YData),nPanel/),"double") 

  if (multimean.eq."True") then 
     do j=0,dimsizes(YData)-1
                xp_M(j,:) = XXdata(dimTot-1,:,j,0)-((X_MMM_std(:,j,0)))
                yp_M(j,:) =YData(j)
                yp_M(dimsizes(YData)+j,:) = YData(dimsizes(YData)-1-j)
            
                xp_M(dimsizes(YData)+j,:) = XXdata(dimTot-1,:,(dimsizes(YData)-1-j),0)+(X_MMM_std(:,(dimsizes(YData)-1-j),0))
     end do
   end if
   xp_R= new((/2*dimsizes(YData),nPanel/),typeof(XXdata))
   yp_R= new((/2*dimsizes(YData),nPanel/),"double")
   ;°°°°°°
   if (isatt(variable_info[0], "reference_dataset")) then
      do j=0,dimsizes(YData)-1
          xp_R(j,:) = XXref(:,j,0)-XXref(:,j,1)
          yp_R(j,:) =YData(j)
          yp_R(dimsizes(YData)+j,:) = YData(dimsizes(YData)-1-j)
          xp_R(dimsizes(YData)+j,:) = XXref(:,(dimsizes(YData)-1-j),0)+XXref(:,(dimsizes(YData)-1-j),1)
      end do
;print(xp_R)
;print(yp_R)
   end if

     
   
;   Ydata= lev0/100.

   plot_d=new(nPanel,"graphic")
   plot_xp_R_d=new(nPanel,"graphic")  
   plot_xp_M_d=new(nPanel,"graphic") 
   plot_R=new(nPanel,"graphic")
;   xp_R= new((/2*dimsizes(lev0),dimsizes(seasons)/),typeof(Xref))
;   yp_R=new((/2*dimsizes(lev0),dimsizes(seasons)/),typeof(Xref))
;
;   xp_M= new((/2*dimsizes(lev0),dimsizes(seasons)/),typeof(Xref))
;   yp_M= new((/2*dimsizes(lev0),dimsizes(seasons)/),typeof(Xref)) 
;   do j=0,dimsizes(lev0)-1   
;         xp_R(j,:) = -((Xref(:,j,1)))
;         yp_R(j,:) = tofloat(lev0(j)/100.)
;         yp_R(dimsizes(lev0)+j,:) = tofloat(lev0(dimsizes(lev0)-1-j)/100)      
;         xp_R(dimsizes(lev0)+j,:) = (Xref(:,(dimsizes(lev0)-1-j),1))
;    end do
;    if (multimean) then 
;        do j=0,dimsizes(lev0)-1
;    ;        printVarSummary(X_MMM_std)
;            xp_M(j,:) = Xdata(dimTot-1,:,j,0)-((X_MMM_std(:,j,0)))
;            yp_M(j,:) = tofloat(lev0(j)/100.)
;            yp_M(dimsizes(lev0)+j,:) = tofloat(lev0(dimsizes(lev0)-1-j)/100.)
;            xp_M(dimsizes(lev0)+j,:) = Xdata(dimTot-1,:,(dimsizes(lev0)-1-j),0)+(X_MMM_std(:,(dimsizes(lev0)-1-j),0))
;         end do
;    end if
;   printVarSummary(Xdata)
;   printVarSummary(tmpData)
;   printVarSummary(Xref)
;*************************plotting*******************************************************
   ;mid_month=(/15.5, 45  , 74.5,105  ,135.5,166 ,196.5,227.5,258  ,288.5,319  ,349.5/)
    month_names=(/"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"/)
 ;print(Xdata(line|:,season|:,day|:,stat|0))
do ipanel=0, nPanel-1
      res                   = True
      res@gsnDraw           = False
      res@gsnFrame          = False
      res@tfPolyDrawOrder   = "Predraw"
      res@trXMinF           = 0 
      res@trXMaxF           = 12
      res@tmXBMode          = "Explicit"
      res@tmXBValues        = ispan(1,12,1) 
      res@tmXBLabels        = month_names
      res@trYMaxF           = YMax 
      res@trYMinF           = YMin
      res@tmLabelAutoStride = True
      res@vpHeightF         = 0.35                    ; change aspect ratio of plot
      res@vpWidthF          = 0.8

     ; res@vpHeightF         = 0.4                    ; change aspect ratio of plot
     ; res@vpWidthF          = 1./dimsizes(seasons)
      res@tiMainString      = sellev+" "+LONG_NAME+" at "+sellat 
      res@tiYAxisString     = LONG_NAME+"[ "+UNITS+" ]"       ; y-axis label
      res@tiXAxisString     = "Month"  
      res@xyMarkLineModes    = LineScatters
      ;print(LineScatters)

      res@xyLineThicknessF= 2.0
      res@xyLineDashSegLenF=0.13
      res@xyLineColors      = plotColor                 ; line color
      res@xyDashPatterns    = plotDash
      res@xyMarkers         =plotScatters
      res@xyMarkerSizeF     =0.01
      res@xyMarkerThicknessF =1.8
      ;res@tmXMajorGrid      = True                   ; implement x grid
      ;res@tmXMajorGridThicknessF = 1.0               ; 2.0 is default
      ;res@tmXMajorGridLineDashPattern = 2            ; select short dash lines
      printVarSummary(XXdata)
      XXdata!0 = "line"
      XXdata!1= "season"
      XXdata!2 = "month"
      XXdata!3 = "stat" 
      print(XXdata(line|:,season|ipanel,month|:,stat|0))
      plot_d(ipanel) = gsn_csm_xy (wks,YData,XXdata(line|:,season|ipanel,month|:,stat|0),res)
      if (isatt(variable_info[0], "reference_dataset")) then
           gsres = True      
           gsres@gsFillOpacityF =0.5
           gsres@gsFillColor = "Gray"
 ;gsres@gsFillOpacityF = 0.5
           resR=True
            
           resR@gsLineColor      = colors(iref)                 ; line color
           resR@gsLineDashPattern=dashes(iref)
           resR@gsFontThicknessF=thicks(iref)
      ;      plot_R(ipanel)=gsn_add_polyline(wks,plot_d(ipanel),YData,Xref(season|ipanel,day|:,stat|0),resR)
            plot_xp_R_d(ipanel) = gsn_add_polygon(wks,plot_d(ipanel),yp_R(:,ipanel),xp_R(:,ipanel),gsres)
      end if
      ;gsres@gsFillOpacityF = 0.5
      if (multimean) then
        gsres1=True
        gsres1@gsFillColor       = "cyan"                 ; color chosen
        gsres1@gsFillOpacityF = 0.2
        plot_xp_M_d(ipanel) = gsn_add_polygon(wks,plot_d(ipanel),yp_M(:,ipanel),xp_M(:,ipanel),gsres1)
      end if
   end do
   ;; Call legend-plotting routine
   Pres = True
   Pres@gsnFrame       = False
   Pres@gsnPanelBottom = 0.2
   ;Pres@gsnPanelFigureStrings=plot_num(0:nPanel-1)
   ;Pres@gsnPanelFigureStringsFontHeightF=0.025
   ;Pres@gsnPanelFigureStringsPerimOn="False"
   Pres@amJust="TopLeft"
;; Call panel routine
   if (nPanel.gt.1) then
      gsn_panel(wks,plot_d,(/2,nPanel/2/),Pres)
   else

      gsn_panel(wks,plot_d,(/1,1/),Pres)
   end if
     print(plotDash)
   
    if (dimsizes (plotMod).gt.21) then
        frame(wks)
    end if
            ;       plotMod(ikk)=names(imod)
            ;           plotColor(ikk)=colors(imod)
            ;           plotDash(ikk)=dashes(imod)
            ;           plotThicks(ikk)=thicks(imod)
;Lin=new(dimsizes(names),"string")
;Lin(:)="Lines"
   
    plotDash(iref)=plotScatters(iref)
    add_legenda(wks,plotMod,plotColor,tofloat(plotDash),LineScatters,plotThicks,"0.2")

   frame(wks)
end

