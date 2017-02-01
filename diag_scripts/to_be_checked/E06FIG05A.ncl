;##################################################################################
;# This NCL script reads in CF compliant netCDF CCMVal model output as specified  #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and produces FIGURE 5a of                                                      #
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
;# Figure 5a. This script is designed to plot vertical profiles (climatological   #
;# (monthly mean) zonal meana) for specific species  (e.g. CH4, H2O, HCl and O3)  #
;# at selected latitudes and where possible compare CCM output with measurements. #
;# Variables plotted can also be output to netcdf files.                          #
;# The number of panels in the plot is determined by the number of latitudes      #
;# selected in the info array (info@fig05_lat_sel)                                #
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
; fig05_xmax_lat        (optional) x axis maximum for profile plots (Fig 5a)
; fig05_xmin_lat        (optional) x axis minimum for profile plots (Fig 5a)
; fig05_ymax_lat        (optional) y axis maximum for profile plots (Fig 5a)
; fig05_ymin_lat        (optional) y axis minimum for profile plots (Fig 5a)
; fig05_obs_name        (optional) observation name
; fig05_obs_file        (optional) file name of file containing observation
;                       processed data (mean) for plot
; fig05_obs_stdv_file   (optional) file name of file containing observation
;                       processed data (stdv) for plot
; fig05_altitude        (optional) y-axis option, pressure or altitude
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
    print("Begin ./plot_type/E06FIG05A.ncl")
   ;*********** get enviromental varibles ****************************
   diri_out = getenv("ESMValTool_plot_dir")
   project_name = getenv("ESMValTool_project")
   OPTIONS_FILE = getenv("ESMValTool_write_plot_vars")

   ;******************************************************************
   ;check the input parameters
   year_min=stringtointeger(models@start_date)
   year_max=stringtointeger(models@end_date)

   if(isvar("lat_sel")) then
      delete(lat_sel)
   end if
   if(isvar("month_names_lat")) then
      delete(month_names_lat)
   end if

   if(isatt(info,"fig05_lat_sel") .and. isatt(info,"fig05_month_lat_sel")) then
      lat_sel = info@fig05_lat_sel
      month_names_lat = info@fig05_month_lat_sel
      do lsel=0,dimsizes(lat_sel)-1
         if (lat_sel(lsel).gt.0.).and.(lat_sel(lsel).lt.90.) then
            lat_sel(lsel)=lat_sel(lsel)+0.1
         end if
         if (lat_sel(lsel).lt.0.).and.(lat_sel(lsel).gt.-90.) then
            lat_sel(lsel)=lat_sel(lsel)-0.1
         end if
      end do
   else
      print("plot requires variable attributes 'fig05_lat_sel' and 'fig05_month_lat_sel'")
      print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("!No files or plots are been produced!")
      print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   end if

   if(isatt(info,"units")) then
      UNITS = info@units
   end if

   ; optional input parameters
   if(isatt(info,"fig05_xmax_lat")) then
      xMax_lat = info@fig05_xmax_lat
   end if
   if(isatt(info,"fig05_xmin_lat")) then
      xMin_lat = info@fig05_xmin_lat
   end if
   if(isatt(info,"fig05_ymax_lat")) then
      yMax_lat = info@fig05_ymax_lat
   end if
   if(isatt(info,"fig05_ymin_lat")) then
      yMin_lat = info@fig05_ymin_lat
   end if
   if(isatt(info,"fig05_altitude")) then
      altitude = info@fig05_altitude
   else
      altitude = False
   end if

   if(isatt(info,"fig05_obs_name") .and. isatt(info,"fig05_obs_file") \
                                .and. isatt(info,"fig05_obs_stdv_file")) then
      OBS = info@fig05_obs_name
      dim_OBS=dimsizes(OBS)
      obs_filename = info@fig05_obs_file
      obs_stdv_file = info@fig05_obs_stdv_file
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

   ;******************************************************************
   dim_MOD = dimsizes(models@name)
   if (OBS.ne."NaN")then
      if (isatt(info,"fig05_obs_file")) then
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
   else
      idx_ref=-999
      Obs_mod = new(dim_MOD,"string")
      Obs_mod= models@name
   end if

   ;******************************************************************
   ;translate the input parameter :
   ;1) month_sel from misc_function.ncl associate the correspondent number for each month
   ;2) lat_names from misc_function.ncl associate a string at any latitude value

   if ((isdefined("lat_sel").eq.True).and.(isdefined("month_names_lat").eq.True)) then
      if ((dimsizes("lat_sel")).ne.(dimsizes("month_names_lat"))) then
         print("fatal:Error!!!!The number of months must to be the same of the number of latitude selected!!!!!")
         print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
         print("fatal:!No files or plots are been produced!")
         print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
         status_exit(1)
      end if
      month_n_lat = month_sel(month_names_lat)
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
   ;***********************  LAT   ***********************************
   if ((isdefined("lat_sel").eq.True)) then

      lat_name = lat_names(info@fig05_lat_sel)
      if altitude then
         wks = gsn_open_wks ("ps",diri_out+"/E06FIG05Aalt_"+variable)
      else
         wks = gsn_open_wks ("ps",diri_out+"/E06FIG05A_"+variable)
      end if

      plot_OBS_lat = new((/dimsizes(lat_sel)/),"graphic")
      plot_lat = new((/dimsizes(lat_sel)/),"graphic")

      ;***************************************************************
      ;################ LOOP ON LATITUDES ############################
      do k=0,dimsizes(lat_sel)-1

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
               ;v_wgt&month = ispan(1,12,1)
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
               ;v_wgt_STD&month = ispan(1,12,1)
               n_start = 0
            else
               A0 = extract_data(idx_ref,data,-1,0,0)
               v_wgt = figure5(A0,field_number,models@name(idx_ref),year_min(idx_ref),year_max(idx_ref))
               v_wgt_STD = figure5_std(A0,field_number,models@name(idx_ref),year_min(idx_ref),year_max(idx_ref))
               delete(A0)
               n_start = 1
            end if
         else
            imod = 0
            A0 = extract_data(imod,data,-1,0,0)
            v_wgt = figure5(A0,field_number,models@name(imod),year_min(imod),year_max(imod))
            delete(A0)
            n_start = 1
         end if
         plev = v_wgt&plev
         lat = v_wgt&lat

         ;********************** plot features ***********************
         hres                  = True
         hres@gsnDraw          = False
         hres@gsnFrame         = False

         ;************************************************************
         ;XAxis
         if (isvar("xMax_lat").eq.True) then
            hres@trMaxF        = xMax_lat
         else
            ;print("automatic selection, trMaxF")
            if (variable.eq."CH4") then
               hres@trXMaxF    = 2.0
            end if
            if (variable.eq."Cly") then
               hres@trXMaxF    = 4.0
            end if
            if (variable.eq."H2O") then
               hres@trXMaxF    = 8.0
            end if
            if (variable.eq."HCl") then
               hres@trXMaxF    = 4.0
            end if
            if (variable.eq."O3") then
               if (lat_sel(k).ge.-30).and.(lat_sel(k).le.30) then
                  hres@trXMaxF = 12.0
               else
                  hres@trXMaxF = 7.0
               end if
            end if
         end if

         if (isvar("xMin_lat").eq.True) then
            hres@trMinF        = xMin_lat
         else
            ;print("automatic selection, trMinF")
            if (variable.eq."CH4") then
               hres@trXMinF    = 0.0
            end if
            if (variable.eq."Cly") then
               hres@trXMinF    = 0.0
            end if
            if (variable.eq."H2O") then
               hres@trXMinF    = 1.0
            end if
            if (variable.eq."HCl") then
               hres@trXMinF    = 0.0
            end if
            if (variable.eq."O3") then
               hres@trXMinF    = 0.0
            end if
         end if

         ;************************************************************
         ;YAxis
         if altitude then
            if (isvar("yMax_lat").eq.True) then
               hres@trMaxF        = yMax_lat
            else
               hres@trYMaxF       = 50.
            end if
            if (isvar("yMin_lat").eq.True) then
               hres@trMinF        = yMin_lat
            else
               hres@trYMinF       = 10.
            end if
            hres@tmYLMode         = "Explicit"
            hres@tmYLValues       = (/10,20,30,40,50/)
            hres@tmYLLabels       = (/10,20,30,40,50/)
            hres@tiYAxisString    = "Height (km)"         ; y-axis label
         else
            if (isvar("yMax_lat").eq.True) then
               hres@trMaxF        = yMax_lat
            else
               hres@trYMaxF       = 200.
            end if
            if (isvar("yMin_lat").eq.True) then
               hres@trMinF        = yMin_lat
            else
               hres@trYMinF       = 1.
            end if
            hres@tmYLMode         = "Explicit"
            hres@tmYLValues       = (/1,2,3,5,7,10,20,30,50,70,100,200/)
            hres@tmYLLabels       = (/1,2,3,5,7,10,20,30,50,70,100,200/)
            hres@tiYAxisString    = "Pressure (hPa)"         ; y-axis label
            hres@trYAxisType      = "LogAxis"
            hres@trYReverse       = True
         end if
         ;************************************************************
         hres@vpHeightF        = 0.5                      ; change aspect ratio of plot
         hres@vpWidthF         = 1./dimsizes(lat_sel)
         hres@tiXAxisString    = variable+" ("+UNITS+")"  ; y-axis label
         hres@tiMainString     = "("+plot_num(k)+") "+variable+" "+lat_name(k)+" "+month_names_lat(k)
         hres@tfPolyDrawOrder  = "Predraw"
         if (dim_OBS.ne.0) then
            if (isatt(info,"fig05_obs_file")) then
               hres@xyMarkLineModes   = "Markers"                ; choose which have markers
               hres@xyMarkers         =  16                      ; choose type of marker
               hres@xyMarkerColor     = "Gray0"                  ; Marker color
               hres@xyMarkerSizeF     = 0.01                     ; Marker size (default 0.01)
            else
               hres@xyLineColors     = color(idx_ref)              ; line color
               hres@xyDashPatterns   = dash(idx_ref)
            end if

            xp = new((/2*dimsizes(v_wgt&plev)/),typeof(v_wgt))
            yp = new((/2*dimsizes(v_wgt&plev)/),typeof(v_wgt))

            do j=0,dimsizes(v_wgt&plev)-1
               xp(j) = v_wgt(month|month_n_lat(k)-1,{lat|lat_sel(k)},plev|j)-v_wgt_STD(month|month_n_lat(k)-1,{lat|lat_sel(k)},plev|j)
               yp(j) = v_wgt&plev(j)
               yp(dimsizes(v_wgt&plev)+j) = v_wgt&plev(dimsizes(v_wgt&plev)-1-j)
               xp(dimsizes(v_wgt&plev)+j) = v_wgt(month|month_n_lat(k)-1,{lat|lat_sel(k)},plev|(dimsizes(v_wgt&plev)-1-j))+v_wgt_STD(month|month_n_lat(k)-1,{lat|lat_sel(k)},plev|(dimsizes(v_wgt&plev)-1-j))
            end do

            plot_lat(k) = gsn_csm_xy(wks,v_wgt(month|month_n_lat(k)-1,{lat|lat_sel(k)},plev|:),plev,hres)

            gsres = True
            gsres@gsFillColor = "Snow3"
            plot_OBS_lat(k) = gsn_add_polygon(wks,plot_lat(k),xp(:),yp(:),gsres)

            delete(xp)
            delete(yp)
         else
            ;if Observation are not included then the base plot should be the plot of the first model
            hres@xyLineColors     = color(imod)              ; line color
            hres@xyDashPatterns   = dash(imod)
            plot_lat(k) = gsn_csm_xy(wks,v_wgt(month|month_n_lat(k)-1,{lat|lat_sel(k)},plev|:),plev,hres)
         end if
         ;************************************************************
         ; OUTPUT FILE
         if (OPTIONS_FILE.eq."True") then
            if (dim_OBS.ne.0) then
               fili = (project_name+"_"+info@fig05_obs_name+"_"+variable+".nc")
            else
               fili = (project_name+"_"+models@case_name(imod)+"_"+models@name(imod)+"_"+models@ensemble(imod)+"_"+field_number+"_"+variable+".nc")
            end if
            f_out = diri_out+"/E06FIG05A/E06FIG05A_"+month_names_lat(k)+"_"+lat_name(k)+"_"+fili

            system("/bin/rm -f "+f_out)
            fout = addfile(f_out,"c")
            setfileoption(fout,"DefineMode",True)

            ;=========================================================
            ; create global attributes of the file
            ;=========================================================
            fAtt               = True            ; assign file attributes
            fAtt@title         = "Climatological-mean data (lev)"
            ; fAtt@source_file = f
            fAtt@creation_date = systemfunc("date")
            fAtt@month         = month_names_lat(k)

            if (dim_OBS.ne.0) then
               fAtt@model         = OBS
            else
               fAtt@model         = models@name(imod)
            end if
            fAtt@latitude      = lat_name(k)
            ; copy attributes from the source file
            fileattdef(fout,fAtt)            ; copy file attributes

            dimNames = (/"plev"/)
            dimSizes = (/dimsizes(plev)/)
            dimUnlim = (/False/)
            filedimdef(fout,dimNames,dimSizes,dimUnlim)

            ;=========================================================
            var_name = variable

            filevardef(fout,"plev",typeof(plev),"plev")
            filevardef(fout,var_name,typeof(v_wgt),(/"plev"/))
            if (dim_OBS.ne.0) then
               var_name2=var_name+"@STD"
               filevardef(fout,var_name2,typeof(v_wgt),(/"plev"/))
            end if
            setfileoption(fout,"DefineMode",False)

            fout->plev       = (/plev/)
            fout->$var_name$ = (/v_wgt(month|month_n_lat(k)-1,{lat|lat_sel(k)},plev|:)/)
            if (dim_OBS.ne.0) then
               fout->$var_name2$=(/v_wgt_STD(month|month_n_lat(k)-1,{lat|lat_sel(k)},plev|:)/)
               delete(v_wgt_STD)
            end if
            delete(fout)
            delete(fAtt)
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print("!")
            print("!  Files are been produced in "+diri_out+"/E06FIG05A/  !!")
            print("!")
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
         else
            if (dim_OBS.ne.0) then
               delete(v_wgt_STD)
            end if
            print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
         end if

         delete(hres)
         delete(v_wgt)
         delete(plev)

         ;####################### LOOP ON MODEL  #####################
         do imod=n_start,dim_MOD-1

            if (imod.ne.idx_ref) then
               A0 = extract_data(imod,data,-1,0,0)
               v_wgt = figure5(A0,field_number,models@name(imod),year_min(imod),year_max(imod))
               plev = v_wgt&plev
               delete(A0)
               if (project_name.eq."CCMVal2") then
                  col = CCMVal2_ColorTable(models@name(imod))
                  das = CCMVal2_DashTable(models@name(imod))
               else
                  col = CCMVal1_ColorTable(models@name(imod))
                  das = CCMVal1_DashTable(models@name(imod))
               end if

               add_line(wks,plot_lat(k),v_wgt(month|month_n_lat(k)-1,{lat|lat_sel(k)},plev|:),plev,col,das)
               delete(col)
               delete(das)
               ;************************************************************
               ; OUTPUT FILE
               if (OPTIONS_FILE.eq."True") then
                  fili = (project_name+"_"+models@case_name(imod)+"_"+models@name(imod)+"_"+models@ensemble(imod)+"_"+field_number+"_"+variable+".nc")
                  f_out = diri_out+"/E06FIG05A/E06FIG05A_"+month_names_lat(k)+"_"+lat_name(k)+"_"+fili

                  system("/bin/rm -f "+f_out)
                  fout = addfile(f_out,"c")
                  setfileoption(fout,"DefineMode",True)

                  ;=========================================================
                  ; create global attributes of the file
                  ;=========================================================
                  fAtt               = True            ; assign file attributes
                  fAtt@title         = "Climatological-mean data (lev)"
                ; fAtt@source_file = f
                  fAtt@creation_date = systemfunc("date")
                  fAtt@month         = month_names_lat(k)
                  fAtt@model         = models@name(imod)
                  fAtt@latitude      = lat_name(k)
                  ; copy attributes from the source file
                  fileattdef(fout,fAtt)            ; copy file attributes

                  dimNames = (/"plev"/)
                  dimSizes = (/dimsizes(plev)/)
                  dimUnlim = (/False/)
                  filedimdef(fout,dimNames,dimSizes,dimUnlim)

                  ;=========================================================
                  var_name = variable

                  filevardef(fout,"plev",typeof(plev),"plev")
                  filevardef(fout,var_name,typeof(v_wgt),(/"plev"/))
                  setfileoption(fout,"DefineMode",False)

                  fout->plev       = (/plev/)
                  fout->$var_name$ = (/v_wgt(month|month_n_lat(k)-1,{lat|lat_sel(k)},plev|:)/)
                  print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                  print("!")
                  print("!  Files are been produced in "+diri_out+"/E06FIG05A/  !!")
                  print("!")
                  print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
               else
                  print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
               end if

               delete(v_wgt)
               delete(plev)
            end if

         ;###############################################################
         end do                ;END LOOP ON MODEL
      end do                   ;END LOOP ON LATITUDE
      ;###############################################################

      dim_lat_sel = dimsizes(lat_sel)

      ;***************************************************************
      Pres                 = True
      Pres@gsnFrame        = False
      Pres@gsnPanelBottom  = 0.33
      Pres@gsnPanelXWhiteSpacePercent = 10
      n_p = dim_lat_sel/3
      if (mod(dim_lat_sel,3) .ne. 0) then
         n_p = n_p+1
      end if
      gsn_panel(wks,plot_lat,(/n_p,3/),Pres)

      add_legenda(wks,Obs_mod,color,dash,"bottom")
      print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("!")
      print("!  Figure E06FIG05A_"+variable+".ps is produced in "+diri_out)
      print("!")
      print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   else
      dim_lat_sel = 0
      print ("No plot created due to info@fig05_lat_sel")
   end if            ;END IF LATITUDE
   delete(month_n_lat)
   delete(dim_lat_sel)
   delete(lat)
   delete(color)
   delete(dim_MOD)
   delete(month_names_lat)
   delete(plot_num)
   delete(year_max)
   delete(n_start)
   delete(lat_sel)
   delete(k)
   delete(Pres)
   delete(i)
   delete(imod)
   delete(models)
   delete(year_min)
   delete(dash)
   delete(Obs_mod)
   delete(lat_name)

   print ("=================== FINISH (E06FIG05A) ===================")

end
;##################################################################################
