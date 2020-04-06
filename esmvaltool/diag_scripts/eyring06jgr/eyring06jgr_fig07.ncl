; #############################################################################
;SCRIPT TO REPRODUCE FIG. 7 OF EYRING ET AL. (JGR, 2006)
; Author: Irene Cionni (ENEA, Italy), Veronika Eyring (DLR, Germany)
; PROJECT-NAME CCMI
;#############################################################################
; Description
;   Climatological mean annual cycle calculated from
;   the monthly mean temperature or H2O mixing ratio,
;    climatological mean for selected levels and latitudinal region
;
; Required diag_script_info attributes (diagnostics specific)
;    multimean: calculate multi-model mean (i.e. False/True)
;    refModel:  name of the reference model or observation for
;                         (e.g. (/"ERA-Interim"/))
;
; Optional diag_script_info attributes (diagnostic specific)
;    Ymax:   max limit Y axis
;    Ymin:   min limit Y axis
;    start_year: year when start the climatology calculation
;    start_year: year when end the climatology calculation
;
; Required variable_info attributes (variable specific)
;
;; Caveats
;;
;; Modification history
;;    20190324-A_cion_ir: cleaned code 
;;    20180320-A_laue_ax: added tags for reporting and call to write_references
;;    20170315-A_cion_ir: written
;;#############################################################################
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/style.ncl"

begin
   enter_msg(DIAG_SCRIPT, "")
   ; Get variables and datasets  
   var0 = variable_info[0]@short_name
   info0 = select_metadata_by_name(input_file_info, var0)
   dim_MOD1 = ListCount(info0)
   names = metadata_att_as_array(info0, "dataset")
   projects = metadata_att_as_array(info0, "project")
   exps = metadata_att_as_array(info0, "exp")
   ensembles = metadata_att_as_array(info0, "ensemble")
   print(metadata_att_as_array(info0, "start_year"))
   print(metadata_att_as_array(info0, "end_year"))
   y1 =(/ metadata_att_as_array(info0, "start_year")/)
   y2 =(/ metadata_att_as_array(info0, "end_year") /)

   ; Check for reference dataset definition
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

   ; Check variable attributes
   if(isatt(variable_info, "units")) then
     units = variable_info@units
   else
     units = "K"
   end if

   if(isatt(variable_info, "long_name")) then
     long_name = variable_info@long_name
   else
     long_name = var0
   end if 

   ; Set default values for optional diag_script_info attributes

   set_default_att(diag_script_info, "styleset", "CMIP5")
   set_default_att(diag_script_info, "Ymin", 185.)
   set_default_att(diag_script_info, "Ymax", 200.)
   set_default_att(diag_script_info, "multimean", False)
   
   set_default_att(diag_script_info, "lev", "100hPa")
   set_default_att(diag_script_info, "latrange", "10S-10N")
   set_default_att(diag_script_info, "start_year", max(y1))
   set_default_att(diag_script_info, "end_year",min(y2))

   flag_mod = where(projects.ne."OBS".and. projects.ne."OBS6" .and. projects.ne."obs4mips", 1, 0)
   index_mod = ind(flag_mod.gt.0)
   index_obs = ind(flag_mod.eq.0)
   dim_MOD = dimsizes(index_mod)
   dim_OBS = 0

   if (.not.all(ismissing(index_obs))) then
     dim_OBS = dimsizes(index_obs)
     list_obs_mod = names(index_obs)
   end if

         
   YMax = diag_script_info@Ymax
   YMin = diag_script_info@Ymin
   sellev = diag_script_info@lev
   sellat = diag_script_info@latrange
   multimean = diag_script_info@multimean
   
   if (multimean)
     dimTot=dim_MOD+dim_OBS+1
   else
     dimTot=dim_MOD+dim_OBS
   end if

   year_start=toint(diag_script_info@start_year)
   year_end=toint(diag_script_info@end_year)
  
   ; Styles
   colors = project_style(info0, diag_script_info, "colors")
   dashes = project_style(info0, diag_script_info, "dashes")
   thicks = tofloat(project_style(info0, diag_script_info, "thicks"))
   markers = project_style(info0, diag_script_info, "markers")
   print(names+" "+colors)
   ; Create output directories
   system("mkdir -p " + config_user_info@work_dir)
   system("mkdir -p " + config_user_info@plot_dir)

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

end

begin
  nPanel = 1
  plot = new(nPanel, "graphic")
  plot_num = integertochar(ispan(97, 122, 1))
  rad = acos(-1.) / 180.0
 
  tmpData = new((/dim_MOD + dim_OBS, nPanel, 12, 2/), "double")
  tmpData!0 = "line"
  tmpData!1 = "season"
  tmpData!2 = "month"
  tmpData!3 = "stat"

  XXdata = new((/dimTot, nPanel, 12, 2/), "double")
  XXdata!0 = "line"
  XXdata!1 = "season"
  XXdata!2 = "month"
  XXdata!3 = "stat"

  XXref= new((/ nPanel, 12, 2/), "double")
  XXref!0 = "season"
  XXref!1 = "month"
  XXdata!2 = "stat"



  plotMod = new(dimTot, typeof(names))
  plotColor = new(dimTot, typeof(names))
  plotDash = new(dimTot, typeof(names))
  plotThicks = new(dimTot, "float")
  plotScatters = new(dimTot, typeof(names))
  LineScatters = new(dimTot, typeof(names)) 
  
  ik=0
  ; Start datasets loop
  do imod = 0, dim_MOD1 - 1
    
    ; Start datasets loop at the present nPanel=!
    do ipan=0,nPanel-1
    ; Extract data
      A0 = read_data(info0[imod])
    
      if (var0.eq."hus") then
           A0 = (/A0 / 0.62/)
           units = "ppmv"
           A0 = (/A0 * 10 ^ 6/)
      end if
      dummy1 = time_operations (A0, year_start, year_end, "average",\
                                "monthlyclim", True)
      dummy2_b=time_operations (A0,year_start,year_end,"extract",\
                                "", True)
      date = yyyymm_time ( year_start, year_end, "integer") 
      yyyy  = date / 100 
      mm  = date - yyyy * 100               
      weights = days_in_month(yyyy, mm)
      dummy2 = dummy2_b ( 0 : 11 )  
      do im = 0, 11
        dummy2(im ) = \
             dim_stddev_wgt_Wrap(dummy2_b(im :: 12), \
                                   weights(im :: 12), 1)
      end do
      delete(A0)
      delete(weights)
      delete(date)
      delete(yyyy)
      delete(mm)
      delete(dummy2_b)
      dummy2!0 = dummy1!0 
      tmpData(imod, ipan, :, 0) = (/dummy1/)
      tmpData(imod, ipan, :, 1) = (/dummy2/)
      delete(dummy1)
      delete(dummy2)
    end do ; end nPanel loop
  end do ;end Dataset loop

  XXdata(:dim_MOD1-1, :, :, :) = (/tmpData/)

  plotMod(:dim_MOD1-1) = names
  plotColor(:dim_MOD1-1) = colors
  plotDash(:dim_MOD1-1) = dashes
  plotThicks(:dim_MOD1-1) = thicks
  plotScatters(:dim_MOD1-1) =markers
  print(index_mod)
  print(index_obs)
  LineScatters(index_mod) = "Lines"
  LineScatters(index_obs) = "Markers"
  XXref = (/tmpData(iref, :, :, :)/)
         
  XXdata!0 = "line"
  XXdata!1 = "season"
  XXdata!2 = "month"
  XXdata!3 = "stat"

  ; Calculate multimodel mean
  if (multimean) then
    XXdata(dimTot-1, :, :, 0) = (/dim_avg_n(XXdata(index_mod, :, :, 0), 0)/)
    plotMod(dimTot-1) = "MMM"
    plotColor(dimTot-1) = "Black"
    plotDash(dimTot-1) = 0.
    LineScatters(dimTot-1) = "Lines"
    plotScatters(dimTot-1) =0.
    plotThicks(dimTot-1) = 4.
  end if
  X_MMM_std=(/dim_stddev_n(XXdata(index_mod,:,:,0),0)/)
  
  YData=ispan(1,12,1)
  YData!0="month"
  YData&month=YData
 
  plot_xp_M_d=new((nPanel),"graphic") 
  
  xp_M= new((/2*dimsizes(YData),nPanel/),typeof(XXdata))
  yp_M= new((/2*dimsizes(YData),nPanel/),"double") 

  if (multimean.eq."True") then 
    do j=0,dimsizes(YData)-1
      xp_M(j, :) = XXdata(dimTot-1, :, j, 0)-((X_MMM_std(:, j)))
      yp_M(j, :) = YData(j)
      yp_M(dimsizes(YData)+j, :) = YData(dimsizes(YData)-1-j)
      xp_M(dimsizes(YData)+j, :) = XXdata(dimTot-1, :, (dimsizes(YData)-1-j), \
                                          0) + (X_MMM_std(:, (dimsizes(YData) \
                                          -1-j)))
    end do
  end if

  xp_R = new((/2 * dimsizes(YData), nPanel/), typeof(XXdata))
  yp_R = new((/2 * dimsizes(YData), nPanel/), "double")
  
  if (isatt(variable_info[0], "reference_dataset")) then
    do j = 0, dimsizes(YData)-1
      xp_R(j, :) = XXref(:, j, 0) - XXref(:, j, 1)
      yp_R(j, :) = YData(j)
      yp_R(dimsizes(YData)+j, :) = YData(dimsizes(YData)-1-j)
      xp_R(dimsizes(YData)+j, :) = XXref(:, (dimsizes(YData)-1-j), 0) + XXref(:,                                         (dimsizes(YData)-1-j), 1)
    end do
  end if
  
  ; Plotting
  plotpath = config_user_info@plot_dir + "eyring06f07_" + var0
  wks = gsn_open_wks(file_type, plotpath)
  plotpath1 = config_user_info@plot_dir + "eyring06f07_legend_" + var0
  wks1 = gsn_open_wks(file_type, plotpath1)


  plot_d = new(nPanel, "graphic")
  plot_xp_R_d = new(nPanel, "graphic")  
  plot_xp_M_d = new(nPanel, "graphic") 
  plot_R = new(nPanel, "graphic")
  month_names = (/"JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG",\
                  "SEP", "OCT", "NOV", "DEC"/)

  do ipanel=0, nPanel-1
    res = True
    res@gsnDraw = False
    res@gsnFrame = False
    res@tfPolyDrawOrder = "Predraw"
    res@trXMinF = 0 
    res@trXMaxF = 12
    res@tmXBMode = "Explicit"
    res@tmXBValues = ispan(1, 12, 1) 
    res@tmXBLabels = month_names
    res@trYMaxF = YMax 
    res@trYMinF = YMin
    res@tmLabelAutoStride = True
    res@vpHeightF = 0.35                    
    res@vpWidthF = 0.8
    res@tiMainString = sellev + " " + long_name + " at " + sellat 
    res@tiYAxisString = long_name + "[ " + units + " ]"
    res@tiXAxisString = "Month"  
    res@xyMarkLineModes = LineScatters
    res@xyLineThicknessF = 2.0
    res@xyLineDashSegLenF =0.13
    res@xyLineColors = plotColor                 ; line color
    res@xyMarkerColors = plotColor  
    res@xyDashPatterns = plotDash
    res@xyMarkers = plotScatters
    res@xyMarkerSizeF = 0.01
    res@xyMarkerThicknessF = 1.8
     print(LineScatters)
     print(plotColor)
     print(plotDash)
     print(plotScatters)
 
    plot_d(ipanel) = gsn_csm_xy (wks, YData, XXdata(line|:, season|ipanel,\
                                 month|:, stat|0), res)
    if (isatt(variable_info[0], "reference_dataset")) then
      gsres = True      
      gsres@gsFillOpacityF = 0.5
      gsres@gsFillColor = "Gray"
;      resR = True
;      resR@gsLineColor = colors(iref)                 
 ;     resR@gsLineDashPattern = dashes(iref)
;      resR@gsFontThicknessF = thicks(iref)
      plot_xp_R_d(ipanel) = gsn_add_polygon(wks, plot_d(ipanel), yp_R(:, \
                                            ipanel), xp_R(:, ipanel), gsres)
      end if
      if (multimean) then
        gsres1 = True
        gsres1@gsFillColor = "cyan"                 ; color chosen
        gsres1@gsFillOpacityF = 0.2
        plot_xp_M_d(ipanel) = gsn_add_polygon(wks, plot_d(ipanel), yp_M(:, \
                                              ipanel), xp_M(:, ipanel), gsres1)
      end if
   end do
   
   ; Call legend-plotting routine
   Pres = True
   Pres@gsnFrame = False
   Pres@gsnPanelBottom = 0.2
   Pres@amJust = "TopLeft"

   if (nPanel.gt.1) then
      gsn_panel(wks,plot_d,(/2,nPanel/2/),Pres)
   else
      gsn_panel(wks,plot_d,(/1,1/),Pres)
   end if

   plotDash(index_obs) = plotScatters(index_obs)
     
   frame(wks)
  print(plotMod+" "+plotColor+" "+plotDash+" "+LineScatters+" "+plotThicks)

   add_legend(wks1, plotMod, plotColor, tofloat(plotDash), LineScatters,\
                plotThicks, "0.9") 
      
   frame(wks1)
   
;common output file
   nc_filename = config_user_info@work_dir + "eyring06f07_" + var0 + ".nc"
   Xdata = XXdata(line|:, season|0, month|:, stat|:)
   Xdata@var = var0
   Xdata@diag_script ="eyring06f07"
   Xdata&stat = (/"mean","stddev"/)
   Xdata&line = plotMod
   Xdata&month = month_names
   nc_outfile = ncdf_write(Xdata, nc_filename)

; Call provenance logger
  log_provenance(nc_filename, \
                 plotpath + "." + file_type, \
                 "Annual cycle of " + var0 + " bias", \
                 (/"mean", "stddev"/), \
                 (/"trop"/), \
                 "seas", \
                 (/"cionni_irene", "eyring_veronika"/), \
                 (/"eyring06jgr"/), \
                 metadata_att_as_array(info0, "filename"))

end
