;##################################################################################
;# This NCL script reads in CF compliant netCDF CCMVal model output as specified  #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and produces FIGURE 12 of                                                      #
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
;# Figure 12a. This script is designed to plot profiles of Cly (or other species  #
;# such as CH4, H2O, HCl, and O3) mixing ratios (climatological (monthly mean)    #
;# zonal mean) at selected latitudes and where possible compares CCM output with  #
;# measurements. Variables plotted can also be output to netcdf files.            #
;# The number of panels in the plot is determined by the number of latitudes      #
;# selected in the info array (info@fig12_lat_sel).                               #
;#                                                                                #
;# Figure 12b: This script is designed to plot time series of Cly (or other       #
;# species such as CH4, H2O,  HCl and O3) mixing ratios (zonal mean monthly mean) #
;# at selected latitudes and pressure levels and where possible, compare CCM      #
;# output with measurements. Variables plotted can also be output to netcdf files.#
;# The number of panels in the plot is determined by the number of latitudes      #
;# selected in the info array (info@fig12_lat_sel).                               #
;;#################################################################################
;# INPUT/OUTPUT:                                                                  #
;# 1. Reads in                                                                    #
;#    (a) '${WKDIR}/var_att/Cly_att.ncl'   (${WKDIR} specified in main namelist)  #
;#    (b) '${WKDIR}/namelist_EO6'         (specification of model names etc.)     #
;#    (c) '${WKDIR}/diag_att/E06Diag.att' (specifies plot type)                   #
;#    (d) Model input files               (as defined in namelists)               #
;# 2. produces                                                                    #
;#    (a) '${diri_out}/E06FIG12.ps'       (plot as ps-file)                       #
;#    (b) '${diri_out}/E06FIG12/'         (plot variables as CF netCDF file)      #
;##################################################################################
; This script takes input from the Cly_att.ncl (CH4_att.ncl, H2O_att.ncl,
; HCl_att.ncl, O3_att.ncl) files
; info attributes required:
; fig12_lat_sel         array of latitudes for profile plots (Fig 12a) and
;                       time series plots (Fir 12b)
; fig12_month_lat_sel   array of months for profile plots (Fig 12a). Note:
;                       this array must have the same dimensions as
;                       fig012_lat_sel
; fig12_lev_sel         array of pressure levels for time series plots (Fig 12b).
;                       Note: this array must have the same dimensions as
;                       fig012_lat_sel
; fig12_month_lev_sel   array of months for time series plots (Fig 12b). Note:
;                       this array must have the same dimensions as
;                       fig012_lat_sel
; fig12_xmax_var        (optional) x axis maximum for profile plots (Fig 12a)
; fig12_xmin_var        (optional) x axis minimum for profile plots (Fig 12a)
; fig12_ymax_var        (optional) y axis maximum for profile plots (Fig 12a)
; fig12_ymin_var        (optional) y axis minimum for profile plots (Fig 12a)
; fig12_ymax_lev        (optional) y axis maximum for time series plots (Fig 12b)
; fig12_ymin_lev        (optional) y axis minimum for time series plots (Fig 12b)
; fig12_xmax_year       (optional) start year for time series plots (Fig 12b)
; fig12_xmin_year       (optional) end year for time series plots (Fig 12b)
;
;**********************************************************************************
load "./plot_type/CCMVal_FUNCTION/new_function.ncl"
load "./plot_type/CCMVal_FUNCTION/Climatology_CCMVal.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
load "./plot_type/CCMVal_FUNCTION/proc_E06FIG12.ncl"
load "./plot_type/CCMVal_FUNCTION/proc_E06FIG05.ncl"

begin
; fields supported T2Mz T2Dz
   print("Begin ./plot_type/E06FIG12.ncl")
   ;*********** get enviromental varibles ****************************
   print("Begin ./plot_type/E06FIG08.ncl")
   diri_out = getenv("ESMValTool_plot_dir")
   project_name = getenv("ESMValTool_project")
   OPTIONS_FILE = getenv("ESMValTool_write_plot_vars")

   dim_MOD = dimsizes(models@name)

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
   if(isatt(info,"fig12_lat_sel") .and. isatt(info,"fig12_month_lat_sel")) then
      lat_sel = info@fig12_lat_sel
      month_names_lat = info@fig12_month_lat_sel
      do lsel=0,dimsizes(lat_sel)-1
         if (lat_sel(lsel).gt.0.).and.(lat_sel(lsel).lt.90.) then
            lat_sel(lsel)=lat_sel(lsel)+0.1
         end if
         if (lat_sel(lsel).lt.0.).and.(lat_sel(lsel).gt.-90.) then
            lat_sel(lsel)=lat_sel(lsel)-0.1
         end if
      end do
   else
      print("fatal:plot requires variable attributes 'fig12_lat_sel' and 'fig12_month_lat_sel'")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if
   if(isvar("lev_sel")) then
      delete(lev_sel)
   end if
   if(isvar("month_names_lev")) then
      delete(month_names_lev)
   end if
   if(isatt(info,"fig12_lev_sel") .and. isatt(info,"fig12_month_lev_sel")) then
      lev_sel = info@fig12_lev_sel
      month_names_lev = info@fig12_month_lev_sel
   else
      print("fatal:plot requires variable attributes 'fig12_lev_sel' and 'fig12_month_lev_sel'")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if

   if(isatt(info,"units")) then
      UNITS = info@units
   end if

   ; optional input parameters
   if(isatt(info,"fig12_xmax_var")) then
      xMax_var = info@fig12_xmax_var
   end if
   if(isatt(info,"fig12_xmin_var")) then
      xMin_var = info@fig12_xmin_var
   end if
   if(isatt(info,"fig12_ymax_var")) then
      yMax_var = info@fig12_ymax_var
   end if
   if(isatt(info,"fig12_ymin_var")) then
      yMin_var = info@fig12_ymin_var
   end if
   if(isatt(info,"fig12_ymax_lev")) then
      yMax_lev = info@fig12_ymax_lev
   end if
   if(isatt(info,"fig12_ymin_lev")) then
      yMin_lev = info@fig12_ymin_lev
   end if
   if(isatt(info,"fig12_xmax_year")) then
      xMax_year = info@fig12_xmax_year
   end if
   if(isatt(info,"fig12_xmin_year")) then
      xMin_year = info@fig12_xmin_year
   end if
   if(isatt(info,"fig12_obs_name") .and. isatt(info,"fig12_obs_file") \
                                .and. isatt(info,"fig12_obs_stdv_file")) then
      OBS = info@fig12_obs_name
      dim_OBS=dimsizes(OBS)
      obs_filename = info@fig12_obs_file
      obs_stdv_file = info@fig12_obs_stdv_file
      idx_ref=-999
   end if

   ;***********************************************************************
   ;translate the input parameter :
   ;1) month_sel from misc_function.ncl associate the correspondent number for each month
   ;2) lat_names from misc_function.ncl associate a string at any latitude value

   if ((isdefined("lat_sel").eq.True).and.(isdefined("month_names_lat").eq.True)) then
      if ((dimsizes("lat_sel")).ne.(dimsizes("month_names_lat"))) then
         print("fatal:Error!!!!The number of months must to be the same of the number of latitude selected!!!!!")
         status_exit(1)
      end if
      month_n_lat = month_sel(month_names_lat)
   end if

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

   if ((dimsizes("lat_sel")).ne.(dimsizes("lev_sel"))) then
      print("fatal:Error!!!!The number of lev must to be the same of the number of latitude selected!!!!!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if

   zlev_sel = (-7*log(lev_sel/10^3))

   if ((isdefined("lat_sel").eq.True)) then
      lat_name = lat_names(info@fig12_lat_sel)
   end if

   ;******************************************************************
   ;The function CCMVal_Color and CCMVal_Dash  in misc_function.ncl
   ; associate at any model a color and dash pattern
   ;******************************************************************
   if (project_name.eq."CCMVal2") then
      color = CCMVal2_ColorTable(models@name)
      dash = CCMVal2_DashTable(models@name)
   else
      color = CCMVal1_ColorTable(models@name)
      dash = CCMVal1_DashTable(models@name)
   end if

   wks = gsn_open_wks ("ps",diri_out+"/E06FIG12_"+variable)
   plot_num = (/"a","b","c","d","e","f","g","h","i","l","m","n","o","p","q","r","s","t","u","v","z"/)

;*********************************************************************
;***********************  LAT   **************************************
   if ((isdefined("lat_sel").eq.True)) then

      plot_lat = new((/dimsizes(lat_sel)/),"graphic")

      ;***************************************************************
      ;################# LOOP ON LATITUDES ###########################
      do k=0,dimsizes(lat_sel)-1

         imod = 0
         A0 = extract_data(imod,data,-1,0,0)
         v_wgt = figure5(A0,field_number,models@name(imod),year_min(imod),year_max(imod))
         plev = v_wgt&plev
         zlev = (-7*log(plev/10^3))
         zlev@units = "Km"
         zlev@long_name = "Height"
         v_wgt&plev = zlev
         ;lat = v_wgt&lat
         delete(A0)

         ;********************** plot features ***********************
         hres                  = True
         hres@gsnDraw          = False
         hres@gsnFrame         = False

         ;************************************************************
         ;XAxis
         if (isvar("xMax_var").eq.True) then
            hres@trXMaxF       = xMax_var
         else
            ;print("automatic selection, trXMaxF")
            if (variable.eq."CH4") then
               hres@trXMaxF    = 2.0
            end if
            if (variable.eq."H2O") then
               hres@trXMaxF    = 8.0
            end if
            if (variable.eq."HCl") then
               hres@trXMaxF    = 4.0
            end if
            if (variable.eq."Cly") then
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

         if (isvar("xMin_var").eq.True) then
            hres@trXMinF       = xMin_var
         else
            ;print("automatic selection, trXMinF")
            if (variable.eq."CH4") then
               hres@trXMinF    = 0.0
            end if
            if (variable.eq."H2O") then
               hres@trXMinF    = 1.0
            end if
            if (variable.eq."HCl") then
               hres@trXMinF    = 0.0
            end if
            if (variable.eq."Cly") then
               hres@trXMinF    = 0.0
            end if
            if (variable.eq."O3") then
               hres@trXMinF    = 0.0
            end if
         end if

         ;************************************************************
         ;YAxis
         if (isvar("yMin_lev").eq.True) then
            hres@trYMinF       = yMin_lev
         else
            hres@trYMinF       = 1.0
         end if
         if (isvar("yMax_lev").eq.True) then
            hres@trYMaxF       = yMax_lev
         else
            hres@trYMaxF       = 200.
         end if

         ;************************************************************
         hres@tmYLMode         = "Explicit"
         hres@tmYLValues       = (/10,20,30,40,50/)
         hres@tmYLLabels       = (/10,20,30,40,50/)
         hres@tmLabelAutoStride = True

         hres@tiYAxisString    = "Height (km)"    ; y-axis label
         hres@tiXAxisString    = variable+" ("+UNITS+")"   ; y-axis label
         hres@tiMainString     = "("+plot_num(k)+") "+variable+" "+lat_name(k)+" "+month_names_lat(k)
         ;hres@trYAxisType = "LogAxis"
         ;hres@trYReverse = True
         hres@xyLineColors     = color(imod)                    ; line color
         hres@xyDashPatterns   = dash(imod)

         plot_lat(k) = gsn_csm_xy(wks,v_wgt(month|month_n_lat(k)-1,{lat|lat_sel(k)},plev|:),zlev,hres)
         delete(hres)

         ;************************************************************
         ; OUTPUT FILE
         if (OPTIONS_FILE.eq."True") then
            fili = (project_name+"_"+models@case_name(imod)+"_"+models@name(imod)+"_"+models@ensemble(imod)+"_"+field_number+"_"+variable+".nc")
            f_out = diri_out+"/E06FIG12/E06FIG12_"+month_names_lat(k)+"_"+lat_name(k)+"_"+fili

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

            dimNames = (/"zlev"/)
            dimSizes = (/dimsizes(zlev)/)
            dimUnlim = (/False/)
            filedimdef(fout,dimNames,dimSizes,dimUnlim)

            ;===================================================================
            var_name = variable

            filevardef(fout,"zlev",typeof(zlev),"zlev")
            filevardef(fout,var_name,typeof(v_wgt),(/"zlev"/))

            setfileoption(fout,"DefineMode",False)

            fout->zlev       = (/zlev/)
            fout->$var_name$ = (/v_wgt(month|month_n_lat(k)-1,{lat|lat_sel(k)},plev|:)/)
            delete(fout)
            delete(fAtt)
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print("!")
            print("!  Files "+f_out+" are been produced in "+diri_out+"/E06FIG12/  !!")
            print("!")
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
         else
            print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
         end if

         delete(v_wgt)
         delete(plev)
         delete(zlev)
         ; delete(lat)

         ;####################### LOOP ON MODEL  #####################
         do imod=1,dim_MOD-1

            A0 = extract_data(imod,data,-1,0,0)
            v_wgt = figure5(A0,field_number,models@name(imod),year_min(imod),year_max(imod))
            plev = v_wgt&plev
            zlev=  (-7*log(plev/10^3))
            zlev@units = "Km"
            zlev@long_name = "Height"
            v_wgt&plev = zlev
            ;lat=v_wgt&lat
            delete(A0)

            add_line(wks,plot_lat(k),v_wgt(month|month_n_lat(k)-1,{lat|lat_sel(k)},plev|:),zlev,color(imod),dash(imod))

            ;************************************************************
            ; OUTPUT FILE
            if (OPTIONS_FILE.eq."True") then
               fili = (project_name+"_"+models@case_name(imod)+"_"+models@name(imod)+"_"+models@ensemble(imod)+"_"+field_number+"_"+variable+".nc")
               f_out = diri_out+"/E06FIG12/E06FIG12_"+month_names_lat(k)+"_"+lat_name(k)+"_"+fili

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

               dimNames = (/"zlev"/)
               dimSizes = (/dimsizes(zlev)/)
               dimUnlim = (/False/)
               filedimdef(fout,dimNames,dimSizes,dimUnlim)

               ;=========================================================
               var_name = variable

               filevardef(fout,"zlev",typeof(zlev),"zlev")
               filevardef(fout,var_name,typeof(v_wgt),(/"zlev"/))

               setfileoption(fout,"DefineMode",False)

               fout->zlev       = (/zlev/)
               fout->$var_name$ = (/v_wgt(month|month_n_lat(k)-1,{lat|lat_sel(k)},plev|:)/)
               print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
               print("!")
               print("!  Files "+f_out+" are been produced in "+diri_out+"/E06FIG012/  !!")
               print("!")
               print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            else
               print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
            end if
            delete(v_wgt)
            delete(plev)
            delete(zlev)
            ; delete(lat)

      ;###############################################################
         end do                ;END LOOP ON MODEL
      end do                   ;END LOOP ON LATITUDE
      ;###############################################################

      dim_lat_sel = dimsizes(lat_sel)
   else
      dim_lat_sel = 0
   end if            ;END LOOP LATITUDE

;*********************************************************************
;***********************  LEV   **************************************
   if((isdefined("lev_sel"))) then

      plot_lev = new((/dimsizes(lev_sel)/),"graphic")

      ;***************************************************************
      ;################ LOOP ON LEVELS ###############################
      do k=0,dimsizes(lev_sel)-1

         imod = 0
         A0 = extract_data(imod,data,-1,0,0)
         v_wgt = figure12(A0,field_number,month_n_lev(k),models@name(imod))
         ;lat = v_wgt&lat
         delete(A0)

         ;************************************************************
         hres                  = True
         hres@gsnDraw          = False
         hres@gsnFrame         = False

         ;************************************************************
         ;YAxis
         if (isvar("yMax_var").eq.True) then
            hres@trYMaxF       = yMax_var
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
               hres@trYMaxF    = 8.0
            end if
            if (variable.eq."O3") then
               hres@trYMaxF    = 5.0
            end if
            if (variable.eq."Cly") then
               hres@trYMaxF    = 4.0
            end if
         end if

         if (isvar("yMin_var").eq.True) then
            hres@trYMinF       = yMin_var
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
              hres@trYMinF     = 0.0
            end if
            if (variable.eq."O3") then
               hres@trYMinF    = 0.0
            end if
            if (variable.eq."Cly") then
               hres@trYMinF    = 0.0
            end if
         end if

         ;************************************************************
         ;XAxis
         xaxis_year = v_wgt&year
         if (typeof(xaxis_year).eq."double") then
            xaxis_year_flt = dble2flt(xaxis_year)
         else
            xaxis_year_flt = xaxis_year
         end if
         if (isvar("xMin_year").eq.True) then
            hres@trXMinF       = xMin_year
         else
            hres@trXMinF       = min(xaxis_year_flt)
         end if
         if (isvar("xMax_year").eq.True) then
            hres@trXMaxF       = xMax_year
         else
            hres@trXMaxF       = max(xaxis_year_flt)
         end if

         ;************************************************************
         hres@tiYAxisString    = variable+" ("+UNITS+")"    ; y-axis label
         hres@tiXAxisString    = " Year "    ; y-axis label
         hres@tiMainString     = "("+plot_num(k+dim_lat_sel)+") "+variable+" "+lev_sel(k)+"hPa "+lat_name(k)+" "+month_names_lev(k)
         hres@xyLineColors     = color(imod)                    ; line color
         hres@xyDashPatterns   = dash(imod)
         ;x=dim_avg_Wrap(v_wgt({plev|lev_sel(k)},year|:,{lat|lat_sel(k)-12.5:lat_sel(k)+12.5}))

         plot_lev(k) = gsn_csm_xy(wks,v_wgt&year,v_wgt({lat|lat_sel(k)},{plev|lev_sel(k)},year|:),hres)
         ;plot_lev(k) = gsn_csm_xy(wks,v_wgt&year,x,hres)
         ;delete(x)
         delete(hres)

         ;************************************************************
         ; OUTPUT FILE
         if (OPTIONS_FILE.eq."True") then
            fili = (project_name+"_"+models@case_name(imod)+"_"+models@name(imod)+"_"+models@ensemble(imod)+"_"+field_number+"_"+variable+".nc")
            f_out = diri_out+"/E06FIG12/E06FIG12_"+month_names_lev(k)+"_"+lat_name(k)+lev_sel(k)+"hPa_"+fili

            system("/bin/rm -f "+f_out)
            fout = addfile(f_out,"c")
            setfileoption(fout,"DefineMode",True)

            ;=========================================================
            ; create global attributes of the file
            ;=========================================================
            fAtt               = True            ; assign file attributes
            fAtt@title         = "Time series data (year)"
            ; fAtt@source_file = f
            fAtt@creation_date = systemfunc("date")
            fAtt@month         = month_names_lev(k)
            fAtt@model         = models@name(imod)
            fAtt@level         = lev_sel(k)
            fAtt@latitude      = lat_sel(k)
            ; copy attributes from the source file
            fileattdef(fout,fAtt)            ; copy file attributes

            dimNames = (/"year"/)
            dimSizes = (/dimsizes(v_wgt&year)/)
            dimUnlim = (/False/)
            filedimdef(fout,dimNames,dimSizes,dimUnlim)

            ;=========================================================
            var_name = variable

            filevardef(fout,"year",typeof(v_wgt&year),"year")
            filevardef(fout,var_name,typeof(v_wgt),(/"year"/))

            setfileoption(fout,"DefineMode",False)

            fout->year       = (/v_wgt&year/)
            fout->$var_name$ = (/v_wgt(year|:,{lat|lat_sel},{plev|lev_sel})/)
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print("!")
            print("!  Files "+f_out+" are been produced in "+diri_out+"/E06FIG12/  !!")
            print("!")
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
         else
            print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
         end if

          delete(v_wgt)

         ;####################### LOOP ON MODEL  #####################
         do imod=1,dim_MOD-1

            A0 = extract_data(imod,data,-1,0,0)
            v_wgt = figure12(A0,field_number,month_n_lev(k),models@name(imod))
            ;lat=v_wgt&lat
            delete(A0)
            ;x=dim_avg_Wrap(v_wgt({plev|lev_sel(k)},year|:,{lat|lat_sel(k)-12.5:lat_sel(k)+12.5}))
            ;add_line(wks,plot_lev(k),v_wgt&year,x,color(imod),dash(imod))
            ;delete(x)
            add_line(wks,plot_lev(k),v_wgt&year,v_wgt({lat|lat_sel(k)},{plev|lev_sel(k)},year|:),color(imod),dash(imod))

            ;************************************************************
            ; OUTPUT FILE
            if (OPTIONS_FILE.eq."True") then
               fili = (project_name+"_"+models@case_name(imod)+"_"+models@name(imod)+"_"+models@ensemble(imod)+"_"+field_number+"_"+variable+".nc")
               f_out = diri_out+"/E06FIG12/E06FIG12_"+month_names_lev(k)+"_"+lat_name(k)+lev_sel(k)+"hPa_"+fili

               system("/bin/rm -f "+f_out)
               fout = addfile(f_out,"c")
               setfileoption(fout,"DefineMode",True)

               ;=========================================================
               ; create global attributes of the file
               ;=========================================================
               fAtt               = True            ; assign file attributes
               fAtt@title         = "Time series data (year)"
               ; fAtt@source_file = f
               fAtt@creation_date = systemfunc("date")
               fAtt@month         = month_names_lev(k)
               fAtt@model         = models@name(imod)
               fAtt@latitude      = lat_sel(k)
               fAtt@level         = lev_sel(k)
               ; copy attributes from the source file
               fileattdef(fout,fAtt)            ; copy file attributes

               dimNames = (/"year"/)
               dimSizes = (/dimsizes(v_wgt&year)/)
               dimUnlim = (/False/)
               filedimdef(fout,dimNames,dimSizes,dimUnlim)

               ;=========================================================
               var_name = variable

               filevardef(fout,"year",typeof(v_wgt&year),"year")
               filevardef(fout,var_name,typeof(v_wgt),(/"year"/))

               setfileoption(fout,"DefineMode",False)

               fout->year       = (/v_wgt&year/)
               fout->$var_name$ = (/v_wgt(year|:,{lat|lat_sel},{plev|lev_sel})/)
               print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
               print("!")
               print("!  Files "+f_out+" are been produced in "+diri_out+"/E06FIG12/  !!")
               print("!")
               print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            else
               print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
            end if

            delete(v_wgt)

         ;###############################################################
         end do                  ;END LOOP MODEL
         if (isatt(info,"fig12_obs_name") .and. isatt(info,"fig12_obs_file") \
                                .and. isatt(info,"fig12_obs_stdv_file")) then
            f_obs = addfile(info@fig12_obs_file,"r")
            v_wgt = f_obs->$variable$      ;Observation by HALOE are already processed:mean
            delete(f_obs)
            f_obs = addfile(info@fig12_obs_stdv_file,"r")
            vName = variable + "@STD"
            v_wgt_STD = f_obs->$vName$     ;Observation by HALOE  are already processed:standard deviation
            delete(f_obs)
            v_wgt!0 = "time"
            v_wgt!1 = "plev"
            v_wgt!2 = "lat"
            ;v_wgt&month = ispan(1,12,1)
            v_wgt_STD!0 = "time"
            v_wgt_STD!1 = "plev"
            v_wgt_STD!2 = "lat"
            v_wgt0 = figure12(v_wgt,field_number,month_n_lev(k),OBS)
            v_wgt0_STD = figure12(v_wgt_STD,field_number,month_n_lev(k),OBS)
            v_wgt1_STD=new((/2,dimsizes(v_wgt0&year)/),typeof(v_wgt_STD))
            tempo=new((/2,dimsizes(v_wgt&time)/),typeof(v_wgt&time))
            indn_lat=ind_nearest_coord(info@fig12_lat_sel(k),v_wgt0&lat,0)
            if (.not.ismissing(indn_lat)) then
               print("Note.The nearest coordinates in the Observation file"+info@fig12_obs_file+" is:"+ v_wgt0&lat(indn_lat))
               mark_color="black"
               mark_ind=7
               add_scatt(wks,plot_lev(k),v_wgt0&year,v_wgt0(lat|indn_lat,{plev|lev_sel(k)},year|:),mark_color,mark_ind)
               ; add_scatt(wks,plot_lev(k),v_wgt0&year,x,mark_color,mark_ind)
            end if
            ; OUTPUT FILE
            if (OPTIONS_FILE.eq."True") then
               fili = (project_name+"_"+info@fig12_obs_name+"_"+variable+".nc")
               f_out = diri_out+"/E06FIG12/E06FIG12_"+month_names_lev(k)+"_"+lat_name(k)+lev_sel(k)+"hPa_"+fili
               system("/bin/rm -f "+f_out)
               fout = addfile(f_out,"c")
               setfileoption(fout,"DefineMode",True)
               ;=========================================================
               ; create global attributes of the file
               ;========================================================
               fAtt               = True            ; assign file attributes
               fAtt@title         = "Time series data (year)"
               ; fAtt@source_file = f
               fAtt@creation_date = systemfunc("date")
               fAtt@month         = month_names_lev(k)
               fAtt@model         = info@fig12_obs_name
               fAtt@latitude      = lat_sel(k)
               fAtt@level         = lev_sel(k)
               ; copy attributes from the source file
               fileattdef(fout,fAtt)            ; copy file attribute
               dimNames = (/"year"/)
               dimSizes = (/dimsizes(v_wgt0&year)/)
               dimUnlim = (/False/)
               filedimdef(fout,dimNames,dimSizes,dimUnlim)
               ;=========================================================
               var_name = variable
               filevardef(fout,"year",typeof(v_wgt0&year),"year")
               filevardef(fout,var_name,typeof(v_wgt),(/"year"/))
               setfileoption(fout,"DefineMode",False)
               fout->year       = (/v_wgt0&year/)
               fout->$var_name$ = (/v_wgt0(year|:,lat|indn_lat,{plev|lev_sel})/)
               print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
               print("!")
               print("!  Files "+f_out+" are been produced in "+diri_out+"/E06FIG12/  !!")
               print("!")
               print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            else
               print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
            end if
            delete (indn_lat)
            indn_lat=ind_nearest_coord(info@fig12_lat_sel(k),v_wgt0_STD&lat,0)
            if (.not.ismissing(indn_lat)) then
               print("Note.The nearest coordinates in the standard deviation Observation file"+info@fig12_obs_stdv_file+" is:"+ v_wgt0_STD&lat(indn_lat))
               v_wgt1_STD(0,:)=v_wgt0(lat|indn_lat,{plev|lev_sel(k)},year|:)+v_wgt0_STD(lat|indn_lat,{plev|lev_sel(k)},year|:)
               v_wgt1_STD(1,:)=v_wgt0(lat|indn_lat,{plev|lev_sel(k)},year|:)-v_wgt0_STD(lat|indn_lat,{plev|lev_sel(k)},year|:)
               v_wgt1_STD!1="time"
               tempo(0,:)=v_wgt0&year
               tempo(1,:)=v_wgt0&year
               do j=0,dimsizes(v_wgt0&year)-1
                  das=0
                  add_line(wks,plot_lev(k),tempo(:,j),v_wgt1_STD(:,j),mark_color,das)
               end do
            end if

            ; OUTPUT FILE
            if (OPTIONS_FILE.eq."True") then
               fili = (project_name+"_"+info@fig12_obs_name+"_"+variable+"_std"+".nc")
               f_out = diri_out+"/E06FIG12/E06FIG12_"+month_names_lev(k)+"_"+lat_name(k)+lev_sel(k)+"hPa_"+fili

               system("/bin/rm -f "+f_out)
               fout = addfile(f_out,"c")
               setfileoption(fout,"DefineMode",True)

               ;=========================================================
               ; create global attributes of the file
               ;=========================================================
               fAtt               = True            ; assign file attributes
               fAtt@title         = "Time series data (year)"
               ; fAtt@source_file = f
               fAtt@creation_date = systemfunc("date")
               fAtt@month         = month_names_lev(k)
               fAtt@model         = info@fig12_obs_name
               fAtt@latitude      = lat_sel(k)
               fAtt@level         = lev_sel(k)
               ; copy attributes from the source file
               fileattdef(fout,fAtt)            ; copy file attributes

               dimNames = (/"year"/)
               dimSizes = (/dimsizes(v_wgt0_STD&year)/)
               dimUnlim = (/False/)
               filedimdef(fout,dimNames,dimSizes,dimUnlim)
               ;=========================================================
               var_name = variable+"@STD"

               filevardef(fout,"year",typeof(v_wgt0_STD&year),"year")
               filevardef(fout,var_name,typeof(v_wgt0_STD),(/"year"/))
               setfileoption(fout,"DefineMode",False)
               fout->year       = (/v_wgt0_STD&year/)
               fout->$var_name$ = (/v_wgt0_STD(year|:,lat|indn_lat,{plev|lev_sel})/)

               print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
               print("!")
               print("!  Files "+f_out+" are been produced in "+diri_out+"/E06FIG12/  !!")
               print("!")
               print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            else
               print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
            end if
            delete(v_wgt)
            delete(v_wgt_STD)
            delete (indn_lat)
         end if
      end do
      ;END LOOP LEVEL
      ;###############################################################

      dim_lev_sel = dimsizes(lev_sel)
   else
      dim_lev_sel = 0
   end if          ;END IF LEVEL

;*********************************************************************
;*********************************************************************

   plot1=new((dim_lat_sel+dim_lev_sel),"graphic")
   Pres               = True
   Pres@gsnFrame      = False
   Pres@gsnPanelBottom= 0.33
   Pres@gsnPanelXWhiteSpacePercent=10
   do j=0,dim_lat_sel-1
      plot1(2*j)=plot_lat(j)
      plot1(2*j+1)=plot_lev(j)
   end do
   gsn_panel(wks,plot1,(/dim_lat_sel,2/),Pres)
   add_legenda(wks,models@name,color,dash,"bottom")
   print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   print("!")
   print("!  Figure E06FIG12_"+variable+".ps is produced in "+diri_out)
   print("!")
   print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

   delete(month_n_lat)
   delete(dim_lat_sel)
   delete(lev_sel)
   delete(color)
   delete(dim_MOD)
   delete(month_names_lat)
   delete(plot_num)
   delete(yMin_lev)
   delete(dim_lev_sel)
   delete(lat_sel)
   delete(k)
   delete(Pres)
   delete(j)
   delete(i)
   delete(imod)
   delete(xaxis_year_flt)
   delete(month_n_lev)
   delete(models)
   delete(xaxis_year)
   delete(month_names_lev)
   delete(yMax_lev)
   delete(dash)
   delete(zlev_sel)
   delete(lat_name)

   print ("=================== FINISH (E06FIG12) ====================")
end
;#################################################################################
