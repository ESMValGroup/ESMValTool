;;#############################################################################
;; EYRING06JGR_FIG03
;; Author: Irene Cionni (ENEA, Italy), Veronika Eyring (DLR, Germany)
;; PROJECT-NAME CCMI
;;#############################################################################
;; Description
;;    Scatter plot var1 (X axis)  vs. var2  (Y axis)
;;    climatological mean for selected month and latitudinal region
;;
;; Required diag_script_info attributes (diagnostics specific)
;;   @e06fig03X_lev: lev where variable var1 is selected  ["10000","10000"]
;;   @e06fig03X_month:   month when variable is averaged ["JF","JA"]
;;   @e06fig03X_latmin:  lat min where the var1 is averaged [40,-80]
;;   @e06fig03X_latmax:  lat max where the var1 is averaged [80,-40]
;;   @e06fig03Y_lev: lev where variable var2 is selected  ["5000","5000"]
;;   @e06fig03Y_month:   month when var2 is averaged ["FM","AS"]
;;   @e06fig03Y_latmin:  lat min where the var2 is averaged [60,-90]
;;   @e06fig03Y_latmax:  lat max where the var2 is averaged [90,-60]
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    @e06fig03X_Max:     max limit X axis
;;    @e06fig03X_Min:     min limit X axis
;;    @e06fig03Y_Max:     max limit Y axis
;;    @e06fig03Y_Min:     min limit Y axis
;;
;; Required variable_info attributes (variable specific)
;;
;; Caveats
;;
;; Modification history
;;    20190228-A_cion_ir: written
;;#############################################################################
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/style.ncl"

begin
   enter_msg(DIAG_SCRIPT, "")
   var0 = variable_info[0]@short_name
   info0 = select_metadata_by_name(input_file_info, var0)

   var1 = variable_info[1]@short_name
   info1 = select_metadata_by_name(input_file_info, var1)

   dim_MOD1 = ListCount(info0)
   names = metadata_att_as_array(info0, "dataset")
   projects = metadata_att_as_array(info0, "project")
   exps = metadata_att_as_array(info0, "exp")
   ensembles = metadata_att_as_array(info0, "ensemble")
   y1 = metadata_att_as_array(info0, "start_year")
   y2 = metadata_att_as_array(info0, "end_year")  
   print(y1)
; y1_var0 = variable_info[0]@reference_dataset@start_year
  ; y2_var0 = variable_info[0]@reference_dataset@end_year 
  ; y1_var1 = variable_info[1]@reference_dataset@start_year
  ; y2_var1 = variable_info[1]@reference_dataset@end_year
   if (isatt(variable_info[0], "reference_dataset")) then
     refname = variable_info[0]@reference_dataset
     do imod = 0, dimsizes(names) - 1
         if (names(imod) .eq. refname) then
            iref = imod
          end if
      end do       
   else
;     error_msg("f", DIAG_SCRIPT, "", "no reference dataset  " + \
;              "(variable_info[0]@reference_dataset) needs to be defined.")
        iref=0
  end if

print(variable_info[0])
print(variable_info[1])
   if(isatt(variable_info[0], "units")) then
        UNITSX = variable_info[0]@units
   end if
   
   if(isatt(variable_info[1], "long_name")) then
        LONG_NAMEX = variable_info[0]@long_name
   end if
      
   if(isatt(variable_info[1], "units")) then
        UNITSY = variable_info[1]@units
   end if

   if(isatt(variable_info[0], "long_name")) then
        LONG_NAMEY = variable_info[1]@long_name
   end if
   log_info("++++++++++++++++++++++++++++++++++++++++++")
   log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
   log_info("++++++++++++++++++++++++++++++++++++++++++")
  
   ; Set default values for optional diag_script_info attributes

   set_default_att(diag_script_info, "styleset", "CMIP5")
   set_default_att(diag_script_info, "e06fig03X_lev", 5000.)
   set_default_att(diag_script_info, "e06fig03X_month", "JF")
   set_default_att(diag_script_info, "e06fig03X_latmin", 40.)
   set_default_att(diag_script_info, "e06fig03X_latmax", 80.)
   set_default_att(diag_script_info, "e06fig03X_Min", 0.)
   set_default_att(diag_script_info, "e06fig03X_Max", 25.)
   set_default_att(diag_script_info, "e06fig03_start_year", max(y1))
   set_default_att(diag_script_info, "e06fig03_end_year", min(y2))

   set_default_att(diag_script_info, "e06fig03Y_lev", 5000.)
    set_default_att(diag_script_info, "e06fig03Y_month", "FM")
   set_default_att(diag_script_info, "e06fig03Y_latmin", 60.)
   set_default_att(diag_script_info, "e06fig03Y_latmax", 90.)
      set_default_att(diag_script_info, "e06fig03Y_Min", 200.)
   set_default_att(diag_script_info, "e06fig03Y_Max", 225.)
   set_default_att(diag_script_info, "e06fig03_multimean", False)
   
   seasonX=diag_script_info@e06fig03X_month
   seasonY=diag_script_info@e06fig03Y_month
  if (dimsizes(seasonX).ne.dimsizes(seasonY)) then
       error_msg("f", DIAG_SCRIPT, "", "dimension @e06fig03X_month needs to be equal to @e06fig03Y_month  ")
  end if


   levX=tofloat(diag_script_info@e06fig03X_lev)
   levY=tofloat(diag_script_info@e06fig03Y_lev)
   levY=levY/100.
   levX=levX/100.

 
   latminX=tofloat(diag_script_info@e06fig03X_latmin)
   latmaxX=tofloat(diag_script_info@e06fig03X_latmax)

   latminY=tofloat(diag_script_info@e06fig03Y_latmin)
   latmaxY=tofloat(diag_script_info@e06fig03Y_latmax)
print(levX)
print(levY)
print(latminX)
print(latmaxX)
print(latminY)
print(latmaxY)
  XMin=tofloat(diag_script_info@e06fig03X_Max) 
  XMAx=tofloat(diag_script_info@e06fig03X_Min)  
    
  YMin=tofloat(diag_script_info@e06fig03Y_Max)    
  YMAx=tofloat(diag_script_info@e06fig03Y_Min)  
  start_year=toint(diag_script_info@e06fig03_start_year)
  end_year=toint(diag_script_info@e06fig03_end_year)

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



   multimean = diag_script_info@e06fig03_multimean
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
  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring06f03_" + var0)
  output_filename = wks@fullname
  nPanel = dimsizes(seasonX)
  
  plot = new(nPanel, "graphic")
  plot_num = integertochar(ispan(97, 122, 1))
  rad = 4.0 * atan(1.0) / 180.0
  A0 = read_data(info0[iref]) 
  years=ispan(toint(start_year),toint(end_year),1)
 nyears=dimsizes(years)  
  print(nyears) 
  tmpDataX = new((/dim_MOD+dim_OBS,nPanel,nyears,2/), "double")
  tmpDataX!0 = "line"
  tmpDataX!1 = "season"
  tmpDataX!2 = "years"
  tmpDataX!3 = "stat"

  tmpDataY = new((/dim_MOD+dim_OBS, nPanel,nyears,2/), "double")
  tmpDataY!0 = "line"
  tmpDataY!1 = "season"
  tmpDataY!2 = "years"
  tmpDataY!3 = "stat"
  
  tmpDataY2 = new((/dim_MOD+dim_OBS, nPanel,nyears,2/), "double")
  tmpDataY2!0 = "line"
  tmpDataY2!1 = "season"
  tmpDataY2!2 = "years"
  tmpDataY2!3 = "stat"

  XdataX = new((/dimTot, nPanel,nyears,2/), "double")
  XdataX!0 = "line"
  XdataX!1 = "season"
  XdataX!2 = "years"
  XdataX!3 = "stat"

  XdataY = new((/dimTot, nPanel,nyears,2/), "double")
  XdataY!0 = "line"
  XdataY!1 = "season"
  XdataY!2 = "years"
  XdataY!3 = "stat"

  XdataY2 = new((/dimTot, nPanel,nyears,2/), "double")
  XdataY2!0 = "line"
  XdataY2!1 = "season"
  XdataY2!2 = "years"
  XdataY2!3 = "stat"
  
  Xref= new((/ nPanel,nyears,2/), "double")
  Yref=new((/ nPanel,nyears,2/), "double")
  plotMod=new(dimTot,typeof(names))
  plotColor=new(dimTot,typeof(names))
  plotDash=new(dimTot,typeof(names))
  plotThicks=new(dimTot,"float")
  plotScatters=new(dimTot,typeof(names))
  LineScatters=new(dimTot,typeof(names)) 
  ik=0


;******************************************************************
  do imod = 0, dim_MOD1 - 1
     do ipan=0,dimsizes(seasonX)-1
      ; Extract data
       AX = read_data(info0[imod])
      AY = read_data(info1[imod])  
      printVarSummary(AX)
      printVarSummary(AY)  
       if(all(isdim(A0,(/"time","plev","lat","lon"/)))) then        
                   dummy1X = area_operations(AX(time|:,{plev|levX(ipan)}, lat|:, lon|:), latminX(ipan), \
                               latmaxX(ipan), 0., 360., "average", True)
                   dummy1Y = area_operations(AY(time|:,{plev|levY(ipan)}, lat|:, lon|:), latminY(ipan), \
                               latmaxY(ipan), 0., 360., "average", True)
                delete(AX)
                delete(AY)

       else
          if(all(isdim(AX,(/"time","plev","lat"/)))) then
             coslat   = cos(AX&lat * rad)
             coslat!0 = "lat"
             coslat&lat = AX&lat
             print(levX(ipan))
             pp=AX&plev
             if (pp.gt.10000) then
              pp=pp/100.           
             AX&plev=pp
              end if
             delete(pp)
             pp=AY&plev
             if (pp.gt.10000) then
              pp=pp/100.
             AY&plev=pp
              end if

           
             dummy1X =dim_avg_wgt_Wrap(AX(time|:,{plev|levX(ipan)}, {lat|latminX(ipan):latmaxX(ipan)}), coslat({latminX(ipan):latmaxX(ipan)}), 1)
             delete (AX) 
             delete(coslat)  
                   coslat   = cos(AY&lat * rad)
             coslat!0 = "lat"
             coslat&lat = AY&lat
             print(levY(ipan))
             print(AY&plev)
             dummy1Y =dim_avg_wgt_Wrap(AY(time|:,{plev|levY(ipan)}, {lat|latminY(ipan):latmaxY(ipan)}), coslat({latminY(ipan):latmaxY(ipan)}), 1)
             delete(coslat) 
             delete(AY)
             
       
          end if
       end if
       print(dummy1X)
       if (seasonX(ipan).ne."ANN") then
                 dummy2X1 = time_operations(dummy1X,start_year,end_year,"extract",seasonX(ipan),True)
                 dummy2Y1 = time_operations(dummy1Y,start_year,end_year,"extract",seasonY(ipan),True)
                 dimseasX= dimsizes(tochar(seasonX(ipan)))
                 dummy2Xf=onedtond(dummy2X1,(/dimseasX,nyears/))   
                 dummy2X=dim_avg_n_Wrap(dummy2Xf,0)
                 dimseasY= dimsizes(tochar(seasonY(ipan)))
                 dummy2Yf=onedtond(dummy2Y1,(/dimseasY,nyears/))         
                 dummy2Y=dim_avg_n_Wrap(dummy2Yf,0)
                 time1=dummy2X1&time(::2)
                 dummy2X!0="time"
                 dummy2X&time=time1
                 time2=dummy2Y1&time(::2)
                 dummy2Y!0="time"
                 dummy2Y&time=time2
                 delete(time1)
                 delete(time2)
                 delete(dummy2X1)
                 delete(dummy2Y1)  
                 delete(dummy2Xf)
                 delete(dummy2Yf)
      else
                   dummy2X = time_operations(dummy1X,start_year,end_year,"extract","yearly",True)
                   dummy2Y = time_operations(dummy1Y,start_year,end_year,"extract","yearly",True)
       
       end if
       delete(dummy1X)
              delete(dummy1Y)
       printVarSummary(dummy2X)
       printVarSummary(dummy2Y)
       print(dummy2X&time)
       print(dummy2Y&time)
       ;Xtmp = dim_avg_n_Wrap(dummy2X,0)
       ;Xtmp_st = dim_stddev_n_Wrap(dummy2X,0)
       ;delete(dummy2X)
       ;Ytmp = dim_avg_n_Wrap(dummy2Y,0)
      ; Ytmp_st = dim_stddev_n_Wrap(dummy2Y,0)
      ; delete(dummy2Y)

       tmpDataX(imod,ipan,:,0) = dummy2X
       ;tmpDataX(imod,ipan,:,1) = 

       tmpDataY(imod,ipan,:,0) = dummy2Y
       ;tmpDataY(imod,ipan,1) = (/Ytmp_st/)
  
       ; **************************regression line************************************
         ind1=ind(.not.ismissing(dummy2X))
         ind2=ind(.not.ismissing(dummy2Y))
         if((.not.all(ismissing(ind1))).or.(.not.all(ismissing(ind1)))) then
              rc=regline(dummy2X(ind1),dummy2Y(ind2))
              tmpDataY2(imod,ipan,:,0)=rc*dummy2X+ rc@yintercept
              delete(dummy2X)
              delete(dummy2Y)
              delete(rc)
              delete(ind1)
              delete(ind2) 
         end if
     end do ; end seasons loop
     XdataX(imod,:,:,:)= (/tmpDataX(imod,:,:,:)/) 
     XdataY(imod,:,:,:)= (/tmpDataY(imod,:,:,:)/)  
     XdataY(imod,:,:,:)= (/tmpDataY2(imod,:,:,:)/)

     plotMod(imod)=names(imod)
     plotColor(imod)=colors(imod)
     plotDash(imod)=dashes(imod)
     plotThicks(imod)=2.
     plotScatters(imod)=markers(imod)
     LineScatters(imod)="Markers"
     
     if (imod.eq.iref) then 
          Xref=tmpDataX(imod,:,:,:)
          Yref=tmpDataY(imod,:,:,:)
          Yref2=tmpDataY2(imod,:,:,:) 
          plotDash(imod)=16. 
          plotThicks(imod)=2.
          plotScatters(imod)=16.
          LineScatters(imod)="Markers"
     end if
   
  end do ;end model loop
  if (multimean) then
       XdataX(dimTot-1,:,:,0)=(/dim_avg_n(XdataX(index_mod,:,:,0),0)/)
       XdataY(dimTot-1,:,:,0)=(/dim_avg_n(XdataY(index_mod,:,:,0),0)/)
       do ipan=0,nPanel-1        
             ind1=ind(.not.ismissing(XdataX(dimTot-1,ipan,:,0)))
             ind2=ind(.not.ismissing(XdataY(dimTot-1,ipan,:,0)))
             rc=regline(XdataX(dimTot-1,ipan,ind1,0),dummy2Y(dimTot-1,ipan,ind2,0))
             XdataY2(dimTot-1,ipan,:,0)=rc*XdataX(dimTot-1,ipan,ind1,0)+ rc@yintercept
            delete(rc)
            delete(ind1)
            delete(ind2)
        end do

  

          plotThicks(dimTot-1)=5
       plotMod(dimTot-1)="MMM"
       plotColor(dimTot-1)="Black"
       plotDash(dimTot-1)=0.
       plotScatters(dimTot-1)=0.
       LineScatters(dimTot-1)="Lines"
  end if


   plot_d=new(nPanel,"graphic")
   line=new(nPanel*dimTot,"graphic")
;   plot_xp_R_d=new(dimsizes(seasons),"graphic")  
;   plot_xp_M_d=new(dimsizes(seasons),"graphic") 
;   xp_R= new((/2*dimsizes(lat0),dimsizes(seasons)/),typeof(Xref))
;   yp_R=new((/2*dimsizes(lat0),dimsizes(seasons)/),typeof(Xref))

;   xp_M= new((/2*dimsizes(lat0),dimsizes(seasons)/),typeof(Xref))
;   yp_M= new((/2*dimsizes(lat0),dimsizes(seasons)/),typeof(Xref)) 
;   do j=0,dimsizes(lat0)-1   
;         yp_R(j,:) = (Xref(:,j,0)-Xref(:,j,1))
;         xp_R(j,:) = tofloat(lat0(j))
;         xp_R(dimsizes(lat0)+j,:) = tofloat(lat0(dimsizes(lat0)-1-j))      
;         yp_R(dimsizes(lat0)+j,:) = (Xref(:,(dimsizes(lat0)-1-j),0)+Xref(:,(dimsizes(lat0)-1-j),1))
;   end do
;   if (multimean) then 
;        do j=0,dimsizes(lat0)-1
;            yp_M(j,:) = Xdata(dimTot-1,:,j,0)-((X_MMM_std(:,j,0)))
;            xp_M(j,:) = tofloat(lat0(j))
;            xp_M(dimsizes(lat0)+j,:) = tofloat(lat0(dimsizes(lat0)-1-j))
;            yp_M(dimsizes(lat0)+j,:) = Xdata(dimTot-1,:,(dimsizes(lat0)-1-j),0)+(X_MMM_std(:,(dimsizes(lat0)-1-j),0))
;;        end do
;   end if
;*************************plotting*******************************************************
   ij=0
    do ipanel=0, dimsizes(seasonX)-1
     res=True
     res@gsnDraw=False
     res@gsnFrame=False
     res@xyMarkerColors=plotColor
     res@xyMarkLineMode     = (/"Markers"/)
     res@xyMarkers=plotScatters
     res@trXMaxF=XMax(ipanel)
     res@trXMinF=XMin(ipanel)
     res@trYMaxF=YMax(ipanel)
     res@trYMinF=YMin(ipanel)
            res@tiYMainString      ="Mean "+LONG_NAMEY+" "+"("+UNITSY+") at "+levY(ipanel)+"hPa "+seasonY(ipanel) 
      res@tiXMainString      ="Mean "+LONG_NAMEX+" "+"("+UNITSX+") at "+levX(ipanel)+"hPa "+seasonX(ipanel) 
      res@xyMarkLineMode    = "Markers"

      res@xyLineColors      = plotColor                 ; line color
      res@xyMarkers         =plotScatters
      res@xyMarkerSizeF     =0.01
      res@xyMarkerThicknessF =1.8

      plot_d(ipanel) = gsn_csm_xy (wks,XdataX(:,ipanel,:,0),XdataY(:,ipanel,:,0),res)
      do imod=0,dimTot-1
         resL = True
         resL@gsLineThicknessF=    plotThicks(imod)
         resL@gsLineColor      = plotColor(imod)
        resL@gsLineDashPattern     = plotDash(imod)
      line(ij)=gsn_add_polyline(wks, plot_d(ipanel), \
                                 XdataX(imod,ipanel,:,0), \
                                 XdataY2(imod,ipanel,:,0), resL)

         ij=ij+1 
      end do 
     


  ;    gsres = True
  ;    gsres@gsFillColor = "Snow3"
  ;    if (isatt(variable_info[0], "reference_dataset")) then
  ;        plot_xp_R_d(ipanel) = gsn_add_polygon(wks,plot_d(ipanel),xp_R(:,ipanel),yp_R(:,ipanel),gsres)
  ;    end if
  ;    if (multimean) then
  ;      gsres1=True
  ;      gsres1@gsFillColor       = "cyan"                 ; color chosen
  ;      gsres1@gsFillOpacityF = 0.2
  ;      plot_xp_M_d(ipanel) = gsn_add_polygon(wks,plot_d(ipanel),xp_M(:,ipanel),yp_M(:,ipanel),gsres1)

   ;   end if
   end do ;end loop panel
   ;; Call legend-plotting routine
   Pres = True
   Pres@gsnFrame       = False
   Pres@gsnPanelBottom = 0.2
   Pres@amJust="TopLeft"
;; Call panel routine

      gsn_panel(wks,plot_d,(/1,nPanel/),Pres)
     frame(wks)
   ; the legend is included in a second figure 
     
    add_legenda(wks,plotMod,plotColor,plotScatters,LineScatters,plotThicks+2.,"0.2")
    frame(wks)
;**********************common output file var0**************************
     nc_filename1 = work_dir + "eyring06f03_" + var0 + ".nc"
     XdataX@var = var0
     XdataX@diag_script="eyring06f03"
     ;Xdata(line|:,season|(ipanel),lat|:,stat|:)
     XdataX&stat=(/"mean","avg"/)
     XdataX&season=seasonX
     XdataX&line=plotMod
     XdataX&years=years
     nc_outfile = ncdf_write(XdataX, nc_filename1)
;***************************************************************
;**********************common output file var1**************************
     nc_filename2 = work_dir + "eyring06f03_" + var1 + ".nc"
     XdataY@var = var1
     XdataY@diag_script="eyring06f03"
     ;Xdata(line|:,season|(ipanel),lat|:,stat|:)
     XdataY&stat=(/"mean","avg"/)
     XdataY&season=seasonX
     XdataY&line=plotMod
     XdataY&years=years
     nc_outfile2 = ncdf_write(XdataY, nc_filename2)
;**********************common output file var1 interpolated**************************
     nc_filename3 = work_dir + "eyring06f03_" + var1 + "interpolated.nc"
     XdataY2@var = var1
     XdataY2@diag_script="eyring06f03"
     ;Xdata(line|:,season|(ipanel),lat|:,stat|:)
     XdataY2&stat=(/"mean","avg"/)
     XdataY2&season=seasonX
     XdataY2&line=plotMod
     XdataY2&years=years
     nc_outfile3 = ncdf_write(XdataY2, nc_filename3)





end

