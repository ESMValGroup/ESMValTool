;##################################################################################
;# This NCL script reads in CF compliant netCDF CCMVal model output as specified  #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and produces FIGURE 5d,e and FIGURE 10 of                                      #
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
;# Description:                                                                   #
;# Figure 5b. This script is designed to plot latitudinal profiles (climatological#
;# (monthly mean) zonal mean) for specific variables (e.g. mean age of air, CH4,  #
;# H2O, HCl and O3) on selected pressure levels and where possible compare CCM    #
;# output with observations. Variables plotted can also be output to netcdf files.#
;# The number of panels in the plot is determined by the number of pressure       #
;# levels selected in the info array (info@fig05_lev_sel)                         #
;##################################################################################
;# INPUT/OUTPUT:                                                                  #
;# 1. Reads in                                                                    #
;#    (a)'${WKDIR}/var_att/${var_name}_att.ncl' (${WKDIR} specified in main       #
;# 				                 namelist)                        #
;#                                               ${var_name} is the name          #
;#                                               of the var to process            #
;#    (b)'${WKDIR}/namelist_EO6'         (specification of model names etc.)      #
;#    (c)'${WKDIR}/diag_att/E06Diag.att' (specifies plot type)                    #
;#    (d) Model input files               (as defined in namelists)               #
;# 2. produces                                                                    #
;#    (a)'${diri_out}/E06FIG05_${var_name}.ps'    (plot as ps-file)               #
;#    (b)'${diri_out}/E06FIG15/'         (plot variables as CF netCDF file)       #
;##################################################################################
; This script takes input from the CH4_att.ncl
; (H2O_att.ncl, HCl_att.ncl, O3_att.ncl) files
; info attributes required:
; fig05_lat_sel         array of latitudes for profile plots (Fig 5a)
; fig05_month_lat_sel   array of months for profile plots (Fig 5a). Note: this
;                       array must have the same dimensions as fig05_lat_sel
; fig05_lev_sel         array of pressure level for latitudinal transects (Fig 5b)
; fig05_month_lev_sel   array of months for transect plots (Fig 5b). Note: this
;                       array must have the same dimensions as fig05_lev_sel
; fig05_xmax_lat        (optional) x axis maximum for profile plots (Fig 5a)
; fig05_xmin_lat        (optional) x axis minimum for profile plots (Fig 5a)
; fig05_ymax_lat        (optional) y axis maximum for profile plots (Fig 5a)
; fig05_ymin_lat        (optional) y axis minimum for profile plots (Fig 5a)
; fig05_xmax_lev        (optional) x axis maximum for latitude transect plots (Fig 5b)
; fig05_xmin_lev        (optional) x axis minimum for latitude transect plots (Fig 5b)
; fig05_ymax_lev        (optional) y axis maximum for latitude transect plots (Fig 5b)
; fig05_ymin_lev        (optional) y axis minimum for latitude transect plots (Fig 5b)
; fig05_obs_name        (optional) observation name as climatological file
; fig05_obs_file        (optional) file name of file containing observation
;                       processed data (mean) for plot
; fig05_obs_stdv_file   (optional) file name of file containing observation
;                       processed data (stdv) for plot
; fig05_refModel        name of the model or observation file as in the namelist_E06
;                       used to calculate the standard deviation.
;                       The standard deviation is plotted as shadow line.
;                       Note: fig05_refModel must to be switch off if fig05_obs_file
;                       is already specified
;
;**********************************************************************************
load "./plot_type/CCMVal_FUNCTION/new_function.ncl"
load "./plot_type/CCMVal_FUNCTION/Climatology_CCMVal.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
load "./plot_type/CCMVal_FUNCTION/proc_E06FIG05.ncl"

begin
; fields supported T3M T3D T2Mz T2Dz
   print("Begin ./plot_type/E06FIG05B.ncl")
   ;*********** get enviromental varibles ****************************

   diri_out = getenv("WKDIR")
   project_name = getenv("PROJECT_NAME")
   OPTIONS_FILE = getenv("WRITE_PLOT_VARS")

   ;******************************************************************
   ;check the input parameters

   year_min=stringtointeger(models@start_date)
   year_max=stringtointeger(models@end_date)

   if(isvar("lev_sel")) then
      delete(lev_sel)
   end if
   if(isvar("month_names_lev")) then
      delete(month_names_lev)
   end if
   if(isatt(info,"fig05_lev_sel") .and. isatt(info,"fig05_month_lev_sel")) then
      lev_sel = info@fig05_lev_sel
      month_names_lev = info@fig05_month_lev_sel
   else
      print("fatal:plot requires variable attributes 'fig05_lev_sel' and 'fig05_month_lev_sel'")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if

   if(isatt(info,"units")) then
      UNITS = info@units
   end if

   ; optional input parameters
   if(isatt(info,"fig05_xmax_lev")) then
      xMax_lev = info@fig05_xmax_lev
   end if
   if(isatt(info,"fig05_xmin_lev")) then
      xMin_lev = info@fig05_xmin_lev
   end if
   if(isatt(info,"fig05_ymax_lev")) then
      yMax_lev = info@fig05_ymax_lev
   end if
   if(isatt(info,"fig05_ymin_lev")) then
      yMin_lev = info@fig05_ymin_lev
   end if

   if(isatt(info,"fig05_obs_name") .and. isatt(info,"fig05_obs_file") \
                             .and. isatt(info,"fig05_obs_stdv_file")) then
      OBS = info@fig05_obs_name
      dim_OBS=dimsizes(OBS)
      obs_filename = info@fig05_obs_file
      obs_stdv_file = info@fig05_obs_stdv_file
      idx_ref=-999
   else

      if(isatt(info,"fig05_obs_name_mean_age") .and. isatt(info,"fig05_obs_file_mean_age") \
                                .and. isatt(info,"fig05_obs_stdv_file_mean_age")) then
         OBS = info@fig05_obs_name_mean_age
         print(info@fig05_obs_file_mean_age)
         print(info@fig05_obs_stdv_file_mean_age)
         print(OBS)
         dim_OBS=dimsizes(OBS)
         obs_filename = info@fig05_obs_file_mean_age
         obs_stdv_file = info@fig05_obs_stdv_file_mean_age
         idx_ref=-999
      else
         if (isatt(info,"fig05_refModel")) then
            print(info@fig05_refModel)

            idx_ref = ind(models@name .eq. info@fig05_refModel)
            if (all(ismissing(idx_ref))) then
               idx_ref=-999
               print("refModel isn't included in model list (namelist)")
               OBS = (/"NaN"/)
               dim_OBS=0
            else
               OBS =info@fig05_refModel
               dim_OBS=dimsizes(OBS)
            end if
         else
            OBS = (/"NaN"/)
            dim_OBS=0
         end if
      end if
   end if
   ;******************************************************************
   dim_MOD = dimsizes(models@name)
   if (OBS.ne."NaN")then
      if  (isatt(info,"fig05_obs_file")) then
         Obs_mod = new(dim_MOD+1,"string")
         Obs_mod(dim_MOD) = OBS
         Obs_mod(:dim_MOD-1) = models@name
      else
         if  (isatt(info,"fig05_obs_file_mean_age")) then
            Obs_mod = new(dim_MOD+1,"string")
            Obs_mod(dim_MOD) = OBS
            Obs_mod(:dim_MOD-1) = models@name
         else
            if (isatt(info,"fig05_refModel")).and.(.not.all(ismissing(idx_ref))) then
               idx_ref = ind(models@name .eq. info@fig05_refModel)
               OBS=info@fig05_refModel
               Obs_mod = new(dim_MOD,"string")
               Obs_mod(dim_MOD-1) = OBS
               jmod=0
               do imod=0,dim_MOD-1
                  if (imod.ne.idx_ref) then
                     Obs_mod(jmod) = models@name(imod)
                     jmod=jmod+1
                  end if
               end do
            else
               dim_OBS = 0
               Obs_mod = new(dim_MOD,"string")
               Obs_mod = models@name
            end if
         end if
      end if
   else
      idx_ref=-999
      Obs_mod = new(dim_MOD,"string")
      Obs_mod= models@name
   end if

   ;******************************************************************
   ;translate the input parameter :
   ;1) month_sel from misc_function.ncl associate the correspondent number for each month
   ;2) lat_names from misc_function.ncl associate a string at any latitude value

   if (isdefined("lev_sel").eq.True).and.(isdefined("month_names_lev").eq.True) then
      if ((dimsizes(lev_sel)).ne.(dimsizes(month_names_lev))) then
         print("fatal:Error!!!!The number of months must to be the same of the number of lev selected!!!!!")
         print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
         print("fatal:!No files or plots are been produced!")
         print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
         status_exit(1)
      end if
      month_n_lev = month_sel(month_names_lev)
   end if

   ;******************************************************************
   ; The function CCMVal_Color and CCMVal_Dash  in misc_function.ncl
   ; associate at any model a color and dash pattern
   ;******************************************************************
   if (project_name.eq."CCMVal2") then
      color = CCMVal2_ColorTable(Obs_mod)
      dash = CCMVal2_DashTable(Obs_mod)
   else
      color = CCMVal1_ColorTable(Obs_mod)
      dash = CCMVal1_DashTable(Obs_mod)
   end if

   plot_num = (/"a","b","c","d","e","f","g","h","i","l","m","n","o","p","q","r","s","t","u","v","z"/)

   ;******************************************************************
   ;***********************  LEV   ***********************************
   if((isdefined("lev_sel"))) then

      wks = gsn_open_wks ("ps",diri_out+"E06FIG05B_"+variable)

      plot_OBS_lev = new((/dimsizes(lev_sel)/),"graphic")
      plot_lev = new((/dimsizes(lev_sel)/),"graphic")

   ;***************************************************************
   ;################ LOOP ON LEVELS ###############################
      do k=0,dimsizes(lev_sel)-1
         if (dim_OBS.ne.0) then
            if (isatt(info,"fig05_obs_file")) then
               f_obs = addfile(obs_filename,"r")
               v_wgt = f_obs->$variable$      ;Observation by HALOE are already processed:mean
               delete(f_obs)
               f_obs = addfile(obs_stdv_file,"r")
               vName = variable + "@STDEV"
               v_wgt_STD = f_obs->$vName$     ;Observation by HALOE  are already processed:standard deviation
               delete(f_obs)
               v_wgt!0 = "month"
               v_wgt!1 = "plev"
               v_wgt!2 = "lat"
               v_wgt_STD!0 = "month"
               v_wgt_STD!1 = "plev"
               v_wgt_STD!2 = "lat"
               if ((v_wgt@units.eq."m^3/m^3").and.(variable.eq."HCl")) then
                  v_wgt=v_wgt*10^9
                  v_wgt_STD=v_wgt_STD*10^9
                  v_wgt@units="ppbv"
                  v_wgt_STD@units="ppbv"
               end if
               if (v_wgt@units.eq."m^3/m^3") then
                  v_wgt=v_wgt*10^6
                  v_wgt_STD=v_wgt_STD*10^6
                  v_wgt@units="ppm"
                  v_wgt_STD@units="ppm"
               end if
               n_start = 0
            else
               if (isatt(info,"fig05_obs_file_mean_age")) then
                  f_obs = addfile(obs_filename,"r")
                  v_wgt = f_obs->$variable$      ;Observation by OBS are already processed:mean
                  delete(f_obs)
                  f_obs = addfile(obs_stdv_file,"r")
                  print(obs_stdv_file)
                  vName = variable + "_std"
                  v_wgt_STD = f_obs->$vName$     ;Observation by OBS  are already processed:standard deviation
                  delete(f_obs)
                  n_start=0
               else
                  A0 = extract_data(models,data,idx_ref,-1,0,0)
                  v_wgt = figure5(A0,field_number,models@name(idx_ref),year_min(idx_ref),year_max(idx_ref))
                  v_wgt_STD = figure5_std(A0,field_number,models@name(idx_ref),year_min(idx_ref),year_max(idx_ref))
                  delete(A0)
                  n_start = 1
               end if
            end if
         else
            imod = 0
            A0 = extract_data(models,data,imod,-1,0,0)
            v_wgt = figure5(A0,field_number,models@name(imod),year_min(imod),year_max(imod))
            delete(A0)
            n_start = 1
         end if

         plev = v_wgt&plev
         lat = v_wgt&lat

         ;************************************************************
         hres                  = True
         hres@gsnDraw          = False
         hres@gsnFrame         = False

         ;************************************************************
         ;YAxis
         if (isvar("yMax_lev").eq.True) then
            hres@trYMaxF       = yMax_lev
         else
            if (variable.eq."CH4") then
               hres@trYMaxF    = 2.0
            end if
            if (variable.eq."H2O") then
               hres@trYMaxF    = 10.0
            end if
            if (variable.eq."HCl") then
               hres@trYMaxF    = 3.5
            end if
            if (variable.eq."mean_age") then
               hres@trYMaxF    = 8
            end if
            if (variable.eq."O3") then
               hres@trYMaxF    = 5.0
            end if
         end if

         if (isvar("yMin_lev").eq.True) then
            hres@trMinF        = yMin_lev
         else
            if (variable.eq."CH4") then
               hres@trYMinF    = 0.4
            end if
            if (variable.eq."H2O") then
               hres@trYMinF    = 0.0
            end if
            if (variable.eq."HCl") then
               hres@trYMinF    = 0.0
            end if
            if (variable.eq."mean_age") then
               hres@trYMinF    = 0.0
            end if
            if (variable.eq."O3") then
               hres@trYMinF    = 0.0
            end if
         end if

         ;************************************************************
         ;XAxis
         if (isvar("xMin_lev").eq.True) then
            hres@trXMinF       = xMin_lev
         else
            hres@trXMinF       = -90.
         end if
         if (isvar("xMax_lev").eq.True) then
            hres@trXMaxF       = xMin_lev
         else
            hres@trXMaxF       = 90.
         end if

         ;************************************************************
         hres@tmXBMode         = "Explicit"
         hres@tmXBValues       = (/-90.,-60.,-30.,0.,30.,60.,90./)
         hres@tmXBLabels       = (/-90.,-60.,-30.,0.,30.,60.,90./)
         hres@tmLabelAutoStride = True

         hres@vpHeightF        = 0.4                    ; change aspect ratio of plot
         hres@vpWidthF         = 0.8/dimsizes(lev_sel)
         hres@tiYAxisString    = variable+" ("+UNITS+")" ; y-axis label
         hres@tiXAxisString    = " Latitude "           ; y-axis label
         hres@tiMainString     = "("+plot_num(k)+") "+variable+" "+lev_sel(k)+" hPa "+month_names_lev(k)
         hres@tfPolyDrawOrder  = "Predraw"

         if (dim_OBS.ne.0) then

            if (isatt(info,"fig05_obs_file")) then
               hres@xyMarkLineModes   = "Markers"                ; choose which have markers
               hres@xyMarkers         =  16                      ; choose type of marker
               hres@xyMarkerColor     = "Gray0"                  ; Marker color
               hres@xyMarkerSizeF     = 0.01                     ; Marker size (default 0.01)
            end if
            if (isatt(info,"fig05_obs_file")) then
               xp = new((/2*dimsizes(v_wgt&lat)/),typeof(v_wgt))
               yp = new((/2*dimsizes(v_wgt&lat)/),typeof(v_wgt))
               do j=0,dimsizes(v_wgt&lat)-1
                  xp(j) = v_wgt(month|month_n_lev(k),lat|j,{plev|lev_sel(k)})+v_wgt_STD(month|month_n_lev(k),lat|j,{plev|lev_sel(k)})
                  yp(j) = v_wgt&lat(j)
                  yp(dimsizes(v_wgt&lat)+j) = v_wgt&lat(dimsizes(v_wgt&lat)-1-j)
                  xp(dimsizes(v_wgt&lat)+j) = v_wgt(month|month_n_lev(k),lat|(dimsizes(v_wgt&lat)-1-j),{plev|lev_sel(k)})- v_wgt_STD(month|month_n_lev(k),lat|(dimsizes(v_wgt&lat)-1-j),{plev|lev_sel(k)})
               end do

               plot_lev(k) = gsn_csm_xy(wks,lat,v_wgt(month|month_n_lev(k),lat|:,{plev|lev_sel(k)}),hres)

               gsres = True
               gsres@gsFillColor = "Snow3"
               plot_OBS_lev(k) = gsn_add_polygon(wks,plot_lev(k),yp(:),xp(:),gsres)

               delete(xp)
               delete(yp)
               delete(v_wgt_STD)
            else
               if (isatt(info,"fig05_obs_file_mean_age")) then
                  hres@xyMarkLineModes   = "Markers"                ; choose which have markers
                  hres@xyMarkers         =  16                      ; choose type of marker
                  hres@xyMarkerColor     = "Gray0"                  ; Marker color
                  hres@xyMarkerSizeF     = 0.01
                  plot_lev(k) = gsn_csm_xy(wks,lat,v_wgt(lat|:,{plev|lev_sel(k)}),hres)
                  v_wgt1_STD=new((/2,dimsizes(v_wgt&lat)/),typeof(v_wgt))
                  v_wgt1_STD(0,:)=v_wgt(lat|:,{plev|lev_sel(k)})+v_wgt_STD(lat|:,{plev|lev_sel(k)})
                  v_wgt1_STD(1,:)=v_wgt(lat|:,{plev|lev_sel(k)})-v_wgt_STD(lat|:,{plev|lev_sel(k)})

                  v_wgt1_STD!0="lat"
                  latit=new((/2,dimsizes(lat)/),typeof(lat))
                  latit(0,:)=lat
                  latit(1,:)=lat
                  do j=0,dimsizes(v_wgt&lat)-1
                     das=0
                     if (.not.any(ismissing(v_wgt1_STD(:,j)))) then
                        add_line(wks,plot_lev(k),latit(:,j),v_wgt1_STD(:,j),"Gray0",das)
                     end if
                  end do
                  delete(v_wgt1_STD)
               else
                  hres@xyLineColors     = color(imod)            ; line color
                  hres@xyDashPatterns   = dash(imod)
                  plot_lev(k) = gsn_csm_xy(wks,lat,v_wgt(month|month_n_lev(k),lat|:,{plev|lev_sel(k)}),hres)
               end if
            end if

            ;************************************************************
            ; OUTPUT FILE
            if (OPTIONS_FILE.eq."True") then
               if (isatt(info,"fig05_refModel")).and.(.not.all(ismissing(idx_ref))) then
                  imod=0
                  fili = (project_name+"_"+models@case_name(imod)+"_"+models@name(imod)+"_"+models@ensemble(imod)+"_"+field_number+"_"+variable+".nc")
                  f_out = diri_out+"/E06FIG05B/E06FIG05B_"+month_names_lev(k)+"_"+lev_sel(k)+"hPa_"+fili
               else
                  fili = (project_name+"_"+OBS+"_"+variable+".nc")
                  if (isatt(info,"fig05_obs_file_mean_age")) then
                     f_out = diri_out+"/E06FIG05B/E06FIG05B_"+lev_sel(k)+"hPa_"+fili
                  else
                     f_out = diri_out+"/E06FIG05B/E06FIG05B_"+month_names_lev(k)+"_"+lev_sel(k)+"hPa_"+fili
                  end if
               end if
               system("/bin/rm -f "+f_out)
               fout = addfile(f_out,"c")
               setfileoption(fout,"DefineMode",True)

               ;=========================================================
               ; create global attributes of the file
               ;=========================================================
               fAtt               = True            ; assign file attributes
               fAtt@title         = "Climatological-mean data (lat)"
               ; fAtt@source_file = f
               fAtt@creation_date = systemfunc("date")
               fAtt@month         = month_names_lev(k)
               if (isatt(info,"fig05_refModel")).and.(.not.all(ismissing(idx_ref))) then
                  imod=0
                  fAtt@model         = models@name(imod)
               else
                  fAtt@model         = OBS
               end if
               fAtt@level         = lev_sel(k)
               ; copy attributes from the source file
               fileattdef(fout,fAtt)            ; copy file attributes

               dimNames = (/"lat"/)
               dimSizes = (/dimsizes(lat)/)
               dimUnlim = (/False/)
               filedimdef(fout,dimNames,dimSizes,dimUnlim)

               ;=========================================================
               var_name = variable

               filevardef(fout,"lat",typeof(lat),"lat")
               filevardef(fout,var_name,typeof(v_wgt),(/"lat"/))
               setfileoption(fout,"DefineMode",False)

               fout->lat        = (/lat/)
               if (isatt(info,"fig05_obs_file_mean_age")) then
                  fout->$var_name$ = (/v_wgt(lat|:,{plev|lev_sel})/)
               else
                  fout->$var_name$ = (/v_wgt(month|month_n_lev(k),lat|:,{plev|lev_sel})/)
               end if
               print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
               print("!")
               print("!  Figure E06FIG05B_"+variable+".ps is produced in "+diri_out)
               print("!")
               print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            else
               print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
            end if
         end if  ;if (dim_OBS.ne.0) then

         delete(hres)
         delete(v_wgt)
         delete(lat)

         ;####################### LOOP ON MODEL  #####################
         do imod=n_start,dim_MOD-1
            print(models@name(imod))
            A0 = extract_data(models,data,imod,-1,0,0)
            v_wgt = figure5(A0,field_number,models@name(imod),year_min(imod),year_max(imod))

            lat = v_wgt&lat
            delete(A0)

            if (project_name.eq."CCMVal2") then
               col = CCMVal2_ColorTable(models@name(imod))
               das = CCMVal2_DashTable(models@name(imod))
            else
               col = CCMVal1_ColorTable(models@name(imod))
               das = CCMVal1_DashTable(models@name(imod))
            end if

            add_line(wks,plot_lev(k),lat,v_wgt(month|month_n_lev(k),lat|:,{plev|lev_sel(k)}),col,das)
            delete(col)
            delete(das)
            ;************************************************************
            ; OUTPUT FILE
            if (OPTIONS_FILE.eq."True") then
               fili = (project_name+"_"+models@case_name(imod)+"_"+models@name(imod)+"_"+models@ensemble(imod)+"_"+field_number+"_"+variable+".nc")
               f_out = diri_out+"/E06FIG05B/E06FIG05B_"+month_names_lev(k)+"_"+lev_sel(k)+"hPa_"+fili

               system("/bin/rm -f "+f_out)
               fout = addfile(f_out,"c")
               setfileoption(fout,"DefineMode",True)

               ;=========================================================
               ; create global attributes of the file
               ;=========================================================
               fAtt               = True            ; assign file attributes
               fAtt@title         = "Climatological-mean data (lat)"
               ; fAtt@source_file = f
               fAtt@creation_date = systemfunc("date")
               fAtt@month         = month_names_lev(k)
               fAtt@model         = models@name(imod)
               fAtt@level         = lev_sel(k)
               ; copy attributes from the source file
               fileattdef(fout,fAtt)            ; copy file attributes
               dimNames = (/"lat"/)
               dimSizes = (/dimsizes(lat)/)
               dimUnlim = (/False/)
               filedimdef(fout,dimNames,dimSizes,dimUnlim)
               ;=========================================================
               var_name = variable
               filevardef(fout,"lat",typeof(lat),"lat")
               filevardef(fout,var_name,typeof(v_wgt) ,(/"lat"/))
               setfileoption(fout,"DefineMode",False)
               fout->lat        = (/lat/)
               fout->$var_name$ = (/v_wgt(month|month_n_lev(k),lat|:,{plev|lev_sel})/)
               print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
               print("!")
               print("!  Figure E06FIG05B_"+variable+".ps is produced in "+diri_out)
               print("!")
               print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            else
               print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
            end if

            delete(lat)
            delete(v_wgt)
         ;###############################################################
         end do                  ;END LOOP MODEL
      end do                     ;END LOOP LEVEL
      ;###############################################################

      dim_lev_sel = dimsizes(lev_sel)

      ;***************************************************************
      Pres                 = True
      Pres@gsnFrame        = False
      Pres@gsnPanelBottom  = 0.33
      Pres@gsnPanelXWhiteSpacePercent = 10
      n_p = dim_lev_sel/3
      if (mod(dim_lev_sel,3) .ne. 0) then
         n_p = n_p+1
      end if
      gsn_panel(wks,plot_lev,(/n_p,3/),Pres)

      add_legenda(wks,Obs_mod,color,dash,"bottom")
      print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("!")
      print("!  Figure E06FIG05B_"+variable+".ps is produced in "+diri_out)
      print("!")
      print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   else
      print ("No plot created due to info@fig05_lev_sel")
   end if          ;END IF LEVEL

   delete(lev_sel)
   delete(color)
   delete(dim_MOD)
   delete(plot_num)
   delete(year_max)
   delete(n_start)
   delete(dim_lev_sel)
   delete(k)
   delete(Pres)
   delete(i)
   delete(imod)
   delete(plev)
   delete(month_n_lev)
   delete(models)
   delete(year_min)
   delete(month_names_lev)
   delete(dash)
   delete(Obs_mod)
;   list_vars()
   print ("=================== FINISH (E06FIG05B) ===================")

end
;##################################################################################
