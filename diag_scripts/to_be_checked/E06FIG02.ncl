;##################################################################################
;# This NCL script reads in CF compliant netCDF CCMVal model output as specified  #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and produces FIGURE 2 of                                                       #
;# Eyring et al., Assessment of temperature, trace species and ozone in           #
;# chemistry-climate model simulations of the recent past, J. Geophys. Res., 111, #
;# D22308, doi:10.1029/2006JD007327, 2006.                                        #
;# (hereinafter "E06")                                                            #
;##################################################################################
;# Authors and Contact:                                                           #
;#                                                                                #
;# Irene Cionni (irene.cionni@enea.it) & Veronika Eyring (veronika.eyring@dlr.de) #
;##################################################################################
;# The NCL script has been developed for use within the CCMVal diagnostic tool    #
;# (PIs: Andrew Gettelman & Veronika Eyring, co-I: Hisako Shiona)                 #
;##################################################################################
;  Description: This script is designed to plot the time series of the            #
;  climatological mean of the zero zonal wind lines at the poles from CCMs and to #
;  compare then to reanalysis data to diagnose the timing of the vortex breakup.  #
;  Variables plotted can also be output to netcdf files.                          #
;  The number of panels in the plot is determined by the number of months         #
;  selected in the info array (info@fig02_first_month ...)                        #
;##################################################################################
;# INPUT/OUTPUT:                                                                  #
;# 1. Reads in                                                                    #
;#    (a) '${WKDIR}/var_att/ta_att.ncl'   (${WKDIR} specified in main namelist)   #
;#    (b) '${WKDIR}/namelist_EO6'         (specification of model names etc.)     #
;#    (c) '${WKDIR}/diag_att/E06Diag.att' (specifies plot type)                   #
;#    (d) Model input files               (as defined in namelists)               #
;# 2. produces                                                                    #
;#    (a) '${diri_out}/E06FIG02.ps'       (plot as ps-file)                       #
;#    (b) '${diri_out}/E06Fig02/'         (plot variables as CF netCDF file)      #
;##################################################################################
; This script takes input from the '${WKDIR}/var_att/ua_att.ncl file
; info attributes required:
; fig02_first_month     array of start months of climatology
; fig02_last_month      array of end months of climatology
; fig02_lat_sel         array of latitudes for plots (one panel per latitude).
;                       Note: his array must have the same dimensions as
;                       fig02_first_month and fig02_last_month
; fig02_refModel        name of the model or observation file as in the namelist_E06
;                       used to calculate the bias and the standard deviation.
;                       The standard deviation is plotted as shadow line.
;**********************************************************************************
load "./plot_type/CCMVal_FUNCTION/new_function.ncl"
load "./plot_type/CCMVal_FUNCTION/Climatology_CCMVal.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
load "./plot_type/CCMVal_FUNCTION/proc_E06FIG02.ncl"
load "./plot_type/get_field_number.ncl"

begin
; fields supported T3M T3D T2Mz T2Dz
   print("Begin ./plot_type/E06FIG02.ncl")
   diri_out = getenv("ESMValTool_plot_dir")
   project_name = getenv("ESMValTool_project")
   OPTIONS_FILE = getenv("ESMValTool_write_plot_vars")
 
   ;******************************************************************
   ;check the input parameters & translate them
   year_min=stringtointeger(models@start_date)
   year_max=stringtointeger(models@end_date)

   if (isatt(info,"fig02_refModel")) then
      OBS=info@fig02_refModel
      idx_ref = ind(models@name .eq. info@fig02_refModel)
      if (all(ismissing(idx_ref))) then
         print("fatal:refModel isn't included in model list (namelist)")
         status_exit(1)
      end if
   else
      print("fatal:no 'refModel'")
      status_exit(1)
   end if

   if(isatt(info,"fig02_first_month")) then
      first_month = info@fig02_first_month
      first_mon_sel = month_sel(first_month)
   else
      print("fatal:plot requires variable attribute 'fig02_first_month'")
      status_exit(1)
   end if

   if(isatt(info,"fig02_last_month")) then
      last_month = info@fig02_last_month
      last_mon_sel = month_sel(last_month)
   else
      print("fatal:plot requires variable attribute 'fig02_last_month'")
      status_exit(1)
   end if

   if(isvar("lat_sel")) then
      delete(lat_sel)
   end if
   if(isvar("lat_sel_name")) then
      delete(lat_sel_name)
   end if
   if(isatt(info,"fig02_lat_sel")) then
      lat_sel = info@fig02_lat_sel
      lat_sel_name = lat_names(lat_sel)
   else
      print("fatal:plot requires variable attribute 'fig02_lat_sel'")
      status_exit(1)
   end if

   ;******************************************************************

   Obs_model = models@name
   dim_MOD=dimsizes(Obs_model)
   colors= new(dim_MOD, "string")
   dashs = new(dim_MOD, "integer")

   names_month = (/"January","February","March","April","May","June","July","August","September","October","November","December"/)

   ;******************************************************************
   ; read obs file and extract data


   if (.not.all(ismissing(idx_ref))) then

      E40 = extract_data(idx_ref,data,-1,0,0)
      field = get_field_number(E40)

      v_sel_E40r = new((/dimsizes(lat_sel),3,12/),"double")

      v_sel_ens = new((/dimsizes(lat_sel),dim_MOD,12/),"double")
      v_sel_std = new((/dimsizes(lat_sel),3,dim_MOD,12/),"double")
   end if
   ;******************************************************************
   month_name_r = new((/dimsizes(lat_sel),12/),"string")
   v_sel_E40 = new((/dimsizes(lat_sel),3,12/),"double")

   mont = (/0,1,2,3,4,5,6,7,8,9,10,11/)
   month = (/1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5/)

   do k=0,dimsizes(lat_sel)-1
      jmod=0
      if (.not.all(ismissing(idx_ref))) then
         Obs_model(jmod)=models@name(idx_ref)
         v_sel_E40(k,:,:) = figure2(E40,field,"ERA_40",lat_sel(k),year_min(idx_ref),year_max(idx_ref))
         v_sel_E40!0 = "sel_lat"

         v_sel_E40r(k,:,0:12-first_mon_sel(k)) = (/v_sel_E40(sel_lat|k,case|:,month|first_mon_sel(k)-1:11)/)
         v_sel_E40r(k,:,12-first_mon_sel(k)+1:11) = (/v_sel_E40(sel_lat|k,case|:,month|0:first_mon_sel(k)-2)/)

         month_name_r(k,0:12-first_mon_sel(k)) = (/names_month(first_mon_sel(k)-1:11)/)
         month_name_r(k,12-first_mon_sel(k)+1:11) = (/names_month(0:first_mon_sel(k)-2)/)

         v_sel_E40r!0 = "sel_lat"
         v_sel_E40r!2 = "month"
         v_sel_E40r&month = mont
         v_sel_ens(k,jmod,:) = (/v_sel_E40r(k,0,:)/)
         jmod=1
         delete(E40)
      end if
      do imod = 0,dim_MOD-1
         if (imod.ne.idx_ref) then
            Obs_model(jmod)=models@name(imod)
            print(models@name(imod))
            A0 = extract_data(imod,data,-1,0,0)
            v_sel = figure2(A0,field_number,models@name(imod),lat_sel(k),year_min(imod),year_max(imod))
            v_sel_ens(k,jmod,0:12-first_mon_sel(k)) = v_sel(case|0,month|first_mon_sel(k)-1:11)
            v_sel_ens(k,jmod,12-first_mon_sel(k)+1:11) = v_sel(case|0,month|0:first_mon_sel(k)-2)
            v_sel_std(k,:,jmod,0:12-first_mon_sel(k)) = v_sel(case|:,month|first_mon_sel(k)-1:11)
            v_sel_std(k,:,jmod,12-first_mon_sel(k)+1:11) = v_sel(case|:,month|0:first_mon_sel(k)-2)
            delete(A0)
            delete(v_sel)
            jmod=jmod+1
         end if
      end do
   end do

   v_sel_std!0 = "sel_lat"
   v_sel_std!1 = "case"
   v_sel_std!2 = "model"
   v_sel_std!3 = "month"
   v_sel_std@units = "hPa"
   v_sel_std&model = Obs_model

   v_sel_ens!0 = "sel_lat"
   v_sel_ens!1 = "model"
   v_sel_ens!2 = "month"
   v_sel_ens@units = "hPa"
   v_sel_ens&model = Obs_model

   if (project_name.eq."CCMVal2") then
      colors = CCMVal2_ColorTable(Obs_model)
      dashs = CCMVal2_DashTable(Obs_model)
   else
      colors = CCMVal1_ColorTable(Obs_model)
      dashs = CCMVal1_DashTable(Obs_model)
   end if

   ;****************  plot  ******************************************
   plot = new(dimsizes(lat_sel),"graphic")
   plot_ERA = new(dimsizes(lat_sel),"graphic")
   lbid  = new(dimsizes(lat_sel),"graphic")

   wks = gsn_open_wks ("ps",diri_out+"/E06FIG02")

   ;***************************************************************
   num_mon = new((/dimsizes(lat_sel)/),"integer")

   do k=0,dimsizes(lat_sel)-1

      if (first_mon_sel(k).gt.last_mon_sel(k)) then
         num_mon(k) = 12-first_mon_sel(k)+last_mon_sel(k)
      else
         num_mon(k) = (last_mon_sel(k)-first_mon_sel(k))
      end if
      if (.not.all(ismissing(idx_ref))) then
         xp = new((/2*num_mon(k)/),typeof(v_sel_E40r))
         yp = new((/2*num_mon(k)/),typeof(v_sel_E40r))
         do j=0,num_mon(k)-1
            xp(j) = ((v_sel_E40r(k,1,j)))
            yp(j) = month(j)
            yp(num_mon(k)+j) = month(num_mon(k)-j)
            xp(num_mon(k)+j) = (v_sel_E40r(k,2,(num_mon(k)-j)))
         end do
      end if
      ;**************  OUTPUT FILES  *********************************
      if (OPTIONS_FILE.eq."True") then
         fili = (project_name+"_"+field_number+"_"+variable+".nc")
         f_out = diri_out+"/E06FIG02/E06FIG02_"+lat_sel_name(k)+"_"+first_month(k)+"_"+fili
         delete(fili)

         system("/bin/rm -f "+f_out)
         fout = addfile(f_out,"c")
         setfileoption(fout,"DefineMode",True)

         ;============================================================
         ; create global attributes of the file
         ;============================================================
         fAtt               = True            ; assign file attributes
         fAtt@title         = "Climatological-mean data (month)"
         ; fAtt@source_file = f
         fAtt@creation_date = systemfunc("date")
         fAtt@latitude      = lat_sel_name(k)
         do im=0,dimsizes(Obs_model)-1
            str=unique_string("model_")
            fAtt@$str$= Obs_model(im)
         end do

         fAtt@comment       = "This file contains the output of FIG02.The  values of Reference model/Observation are also included"
         ; copy attributes from the source file
         fileattdef(fout,fAtt)            ; copy file attributes

         dimNames = (/"model","month","case"/)
         dimSizes = (/dimsizes(Obs_model),dimsizes(month),3/)
         dimUnlim = (/False,False,False/)
         filedimdef(fout,dimNames,dimSizes,dimUnlim)

         ;============================================================
         var_name = "Pressure"
         var_name3 = "Pressure_plus_minus_std"

         filevardef(fout,"month",typeof(month),"month")
         filevardef(fout,var_name,typeof(v_sel_ens),(/"model","month"/))
         filevardef(fout,var_name3,typeof(v_sel_std),(/"case","model","month"/))
         if (.not.all(ismissing(idx_ref))) then
            var_name2= "Std_Pressure_"+info@fig02_refModel
            filevardef(fout,var_name2,typeof(v_sel_E40r),(/"case","month"/))
         end if

         setfileoption(fout,"DefineMode",False)

         fout->month      = (/month/)
         fout->$var_name$ = (/v_sel_ens(sel_lat|k,model|:,month|:)/)
         fout->$var_name3$ = (/v_sel_std(sel_lat|k,case|:,model|:,month|:)/)
         if (.not.all(ismissing(idx_ref))) then
            fout->$var_name2$ = (/v_sel_E40r(k,:,:)/)
         end if
      else
         print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
      end if

      ;***************************************************************
      res                      = True
      res@gsnDraw              = False
      res@gsnFrame             = False
      res@trYAxisType          = "LogAxis"
      res@trYReverse           = True
      res@tfPolyDrawOrder      = "Predraw"
      res@trXMaxF              = month(num_mon(k))
      res@trXMinF              = month(0)
      res@trYMaxF              = 100
      res@trYMinF              = 0.1
      res@tmYLMode             = "Explicit"
      res@tmYLValues           = (/.1,.3,1,3,10,30,100/)
      res@tmYLLabels           = (/.1,.3,1,3,10,30,100/)
      res@tmLabelAutoStride    = True

      res@tmXBMode             = "Explicit"
      num_mesi = ispan(1,num_mon(k)+1,1)
      res@tmXBValues           = num_mesi
      res@tmXBLabels           = month_name_r(k,0:num_mon(k))

      res@vpHeightF            = 0.7/dimsizes(lat_sel)   ; change aspect ratio of plot
      res@vpWidthF             = 1.5*0.7/dimsizes(lat_sel)

      res@tiYAxisString        = "Pressure (hPa)"        ; y-axis label
      res@tiMainString         = "Transition to easterlies at "+ lat_sel_name(k)     ; y-axis label
      res@tiXAxisString        = "Date"
      res@trYAxisType          = "LogAxis"
      res@trYReverse           = True
      res@xyLineColors         = colors                  ; line color
      res@xyDashPatterns       = dashs

      plot(k) = gsn_csm_xy (wks,month(0:num_mon(k)),v_sel_ens(k,:,0:num_mon(k)),res)

      gsres                    = True
      gsres@gsFillColor        = "Snow3"
      plot_ERA(k) = gsn_add_polygon(wks,plot(k),yp(:),xp(:),gsres)

      lgres                    = True
      lgres@lgMonoDashIndex    = False
      lgres@lgLineColors       = colors(::-1)
      lgres@lgLineDashSegLenF  = 0.07
      lgres@lgDashIndexes      = dashs(::-1)
      lgres@vpWidthF           = 0.21
      lgres@vpHeightF          = 0.015+0.7/dimsizes(lat_sel)
      lgres@lgLabelFontHeightF = 0.13
      labels = new(dimsizes(Obs_model),"string")
      labels = Obs_model(::-1)

      nitems = dimsizes(labels)
      lbid(k) = gsn_create_legend(wks,nitems,labels,lgres)

      amres                    = True
      amres@amParallelPosF     = 0.73
      amres@amOrthogonalPosF   = 0.0
      annoid = gsn_add_annotation(plot(k),lbid(k),amres)

      delete(num_mesi)
      delete(res)
      delete(xp)
      delete(yp)
   end do

   Pres = True
   gsn_panel(wks,plot,(/dimsizes(lat_sel)+1,1/),Pres)

   delete(first_month)
   delete(lat_sel_name)
   delete(v_sel_E40r)
   delete(v_sel_E40)
   delete(nitems)
   delete(dim_MOD)
   delete(first_mon_sel)
   delete(month)
   delete(last_mon_sel)
   delete(colors)
   delete(plot_types)
   delete(lat_sel)
   delete(mont)
   delete(k)
   delete(j)
   delete(Pres)
   delete(i)
   delete(labels)
   delete(lgres)
   delete(v_sel_ens)
   delete(gsres)
   delete(imod)
   delete(amres)
   delete(num_mon)
   delete(names_month)
   delete(month_name_r)
   delete(Obs_model)
   delete(year_min)
   delete(dashs)

   print ("=================== FINISH (E06FIG02) ====================")

end
;##################################################################################
