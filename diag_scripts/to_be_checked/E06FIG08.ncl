;##################################################################################
;# This NCL script reads in CF compliant netCDF CCMVal model output as specified  #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and produces FIGURE 8 of                                                       #
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
;# Description: This script is designed to plot time-height sections of water     #
;# vapor mixing ratio shown as the deviation (in parts per million by volume)     #
;# from the time mean profile from CCMs and reanalysis. Used to diagnose the      #
;# water vapor stratospheric tape recorder signal in CCMs. Variables plotted can  #
;# also be output to netcdf files.                                                #
;##################################################################################
;# INPUT/OUTPUT:                                                                  #
;# 1. Reads in                                                                    #
;#    (a) '${WKDIR}/var_att/H2O_att.ncl'  (${WKDIR} specified in main namelist)   #
;#    (b) '${WKDIR}/namelist_EO6'         (specification of model names etc.)     #
;#    (c) '${WKDIR}/diag_att/E06Diag.att' (specifies plot type)                   #
;#    (d) Model input files               (as defined in namelists)               #
;# 2. produces                                                                    #
;#    (a) '${diri_out}/E06FIG08.ps'       (plot as ps-file)                       #
;#    (b) '${diri_out}/E06Fig08/'         (plot variables as CF netCDF file)      #
;##################################################################################
; This script takes input from the /var_att/H2O_att.ncl file
; info attributes required:
; fig08_lat_min         start latitude for spatial averaging
; fig08_lat_max         end latitude for spatial averaging
; fig08_obs_name        (optional) observation name
; fig08_obs_file        (optional) file name of file containing observation
;                       already as climatological mean
;
;**********************************************************************************
load "./plot_type/CCMVal_FUNCTION/new_function.ncl"
load "./plot_type/CCMVal_FUNCTION/Climatology_CCMVal.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
load "./plot_type/CCMVal_FUNCTION/proc_E06FIG08.ncl"

begin
; fields supported T3M T3D T2Mz T2Dz
   print("Begin ./plot_type/E06FIG08.ncl")
   diri_out = getenv("ESMValTool_plot_dir")
   project_name = getenv("ESMValTool_project")
   OPTIONS_FILE = getenv("ESMValTool_write_plot_vars")

   ;******************************************************************
   ;check the input parameters & translate them
   year_min=stringtointeger(models@start_date)
   year_max=stringtointeger(models@end_date)

   if(isvar("lat_min")) then
      delete(lat_min)
   end if
   if(isvar("lat_max")) then
      delete(lat_max)
   end if
   if(isvar("lat_name_min")) then
      delete(lat_name_min)
   end if
   if(isvar("lat_name_max")) then
      delete(lat_name_max)
   end if
   if(isatt(info,"fig08_lat_min") .and. isatt(info,"fig08_lat_max")) then
      lat_min = info@fig08_lat_min
      lat_max = info@fig08_lat_max
      lat_name_max = lat_names(lat_max)
      lat_name_min = lat_names(lat_min)
   else
      print("fatal:plot requires variable attributes 'fig08_lat_min' and 'fig08_lat_max'")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if

   ; optional input parameters
   if(isatt(info,"fig08_obs_name") .and. isatt(info,"fig08_obs_file")) then
      OBS = info@fig08_obs_name
      obs_filename = info@fig08_obs_file
      NOBS = "False"
   else
      NOBS = "True"
   end if

   ;******************************************************************
   dim_MOD = dimsizes(models@name)
   if (NOBS.ne."True") then 
      dim_OBS = dimsizes(OBS)
      Obs_mod = new(dim_OBS+dim_MOD,"string")
      Obs_mod(:dim_OBS-1) = OBS
      Obs_mod(dim_OBS:) = models@name
   else
      dim_OBS = 0
      Obs_mod = new(dim_MOD,"string")
      Obs_mod = models@name
   end if

   wks = gsn_open_wks ("ps",diri_out+"/E06FIG08")

   ; open workstation

   cmap = "WhViBlGrYeOrRe"
   gsn_define_colormap(wks,cmap)

   plot = new((/dim_OBS+dim_MOD/),"graphic")

;*****************************************************************
;OBS

   if (dim_OBS.ne.0) then
      do iobs=0,dim_OBS-1
         f_obs = addfile(obs_filename(iobs),"r")
         A1 = f_obs->$variable$      ;Observation by HALOE are already processed:mean
         variable2=variable+"@STDEV"
         A2 = f_obs->$variable2$
         A1!0 = "month"
         A1!1 = "plev"
         A1!2 = "lat"
         A2!0 = "month"
         A2!1 = "plev"
         A2!2 = "lat"
         A0= new((/2,dimsizes(A1&month),dimsizes(A1&plev),dimsizes(A1&lat)/),typeof(A1))
         A0(0,:,:,:)=A1(month|:,plev|:,lat|:)
         A0(1,:,:,:)=A2(month|:,plev|:,lat|:)
         A0!0="case"
         A0!1 = "month"
         A0!2 = "plev"
         A0!3 = "lat"
         A0&month=A1&month
         A0&plev=A1&plev
         A0&lat=A1&lat

         delete(A1)
         delete(A2)
         if (A0@units.eq."m^3/m^3") then
            A0=A0*10^6
            A0@units="ppm"
         end if
         print("processing "+OBS(iobs)+" ("+(iobs+1)+" of " +dim_OBS+")")

         v_wgt = figure8(A0,"C2Mz",lat_max,lat_min,1980,1990)
         plev = v_wgt&plev

         delete(A0)
         if (OPTIONS_FILE.eq."True") then
            fili = str_get_field(obs_filename(iobs),str_fields_count(obs_filename(iobs),"/"),"/")
            f_out = diri_out+"/E06FIG08/E06FIG08_"+OBS(iobs)+"_"+lat_name_max+lat_name_min+"_"+fili
            system("/bin/rm -f "+f_out)
            fout = addfile(f_out,"c")
            setfileoption(fout,"DefineMode",True)
            ;=========================================================
            ; create global attributes of the file
            ;=========================================================
            fAtt               = True           ;assign file attributes
            fAtt@title         = "Climatological data (plev,month)"
            ; fAtt@source_file = f
            fAtt@creation_date = systemfunc("date")
            fAtt@model         = OBS(iobs)
            ; copy attributes from the source file
            fileattdef(fout,fAtt)            ; copy file attributes

            dimNames = (/"plev","month"/)
            dimSizes = (/dimsizes(plev),dimsizes(v_wgt&month)/)
            dimUnlim = (/False,False/)
            filedimdef(fout,dimNames,dimSizes,dimUnlim)

            ;=========================================================
            var_name = variable
            var_name_std =variable+"_std"

            filevardef(fout,"plev",typeof(plev),"plev")
            filevardef(fout,var_name,typeof(v_wgt),(/"plev","month"/))
            filevardef(fout,var_name_std,typeof(v_wgt),(/"plev","month"/))
            setfileoption(fout,"DefineMode",False)

            fout->plev       = (/plev/)
            fout->$var_name$ = (/v_wgt(case|0,plev|:,month|:)/)
            fout->$var_name_std$ =(/v_wgt(case|1,plev|:,month|:)/)
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print("!")
            print("!  Files are been produced in "+diri_out+"/E06FIG08/  !!")
            print("!")
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
         else
            print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
         end if

         ;============================================================
         res                      = True                  ; plot mods desired
         res@tiMainString         = OBS(iobs)     ; title
         res@gsnDraw              = False
         res@gsnFrame             = False
         res@lbLabelBarOn         = False
         res@cnLineLabelsOn       = False
         res@tiYAxisFont          = "helvetica"
         res@tiMainFont           = "helvetica"
         res@tmXBLabelFont        = "helvetica"
         res@tmYLLabelFont        = "helvetica"
         res@lbLabelFont          = "helvetica"
         res@cnLevelSelectionMode = "ManualLevels"        ; manual contour levels
         res@cnLevelSpacingF      = 0.1                   ; contour interval
         res@cnMinLevelValF       = -1.0                  ; min level
         res@cnMaxLevelValF       = 1.0                   ; max level
         res@cnFillOn             = True                  ; turn on color fill
         res@cnInfoLabelOn        = False                 ; turn off contour label
         res@gsnSpreadColors      = True                  ; use full range of colors
         res@gsnSpreadColorStart  = 10
         res@lbLabelAutoStride    = True                  ; optimal labels
         res@gsnLeftString        = ""
         res@gsnRightString       = ""
         res@gsnLeftStringFontHeightF = 0.015
         res@tmXBMode             = "Explicit"            ; explicit labels
         ; location of labels
         res@gsnAddCyclic         = True
         res@tmXBValues = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23/)
         res@tmXBLabels = (/"J","F","M","A","M","J","J","A","S","O","N","D","J","J","F","M","A","M","J","J","A","S","O","N","D","J"/)                ; labels themselves
         ; res@tmYLValues = (/100,70,50,30,20,10/)
         ; res@tmYLLabels = (/"100","70","50","30","20","10"/)
         res@trYMaxF              = 100
         res@trYMinF              = 10
         ;*** COMMENT_OUT FOR NCL6.0.0 ***  res@trYAxisType          = "LogAxis"
         res@vpHeightF            = 0.4                    ; change aspect ratio of plot
         res@vpWidthF             = 0.75
         res@lbLabelFontHeightF   = 0.011
         res@tiYAxisString        = "Pressure ["+plev@units+"]"    ; y-axis label

         plot(iobs) = gsn_csm_pres_hgt (wks,v_wgt(case|0,plev|:,month|:),res)
         delete(plev)
         delete(v_wgt)
         delete(res)
      end do
   end if

;*****************************************************************
;MODEL

   do imod=0,dim_MOD-1

      A0 = extract_data(imod,data,-1,0,0)
      print("processing "+models@name(imod)+" ("+(imod+1)+" of " +dim_MOD+")")
      v_wgt = figure8(A0,field_number,lat_max,lat_min,year_max(imod),year_min(imod))
      plev = v_wgt&plev

      delete(A0)

      if (OPTIONS_FILE.eq."True") then

         fili = "E06FIG08/E06FIG08_"+project_name+"_"+Obs_mod(imod)+"_"+lat_name_max+lat_name_min+"_"+variable+".nc"
         f_out = diri_out+"/"+fili
         print("OUTPUT FILE: "+f_out)
         system("/bin/rm -f "+f_out)
         fout = addfile(f_out,"c")
         setfileoption(fout,"DefineMode",True)

         ;============================================================
         ; create global attributes of the file
         ;============================================================
         fAtt               = True            ; assign file attributes
         fAtt@title         = "Climatological data (plev,month)"
         ; fAtt@source_file = f
         ;fAtt@creation_date = systemfunc("date")
         fAtt@model= Obs_mod(imod)

         ; copy attributes from the source file
         fileattdef(fout,fAtt)            ; copy file attributes

         dimNames = (/"plev","month"/)
         dimSizes = (/dimsizes(plev),dimsizes(v_wgt&month)/)
         dimUnlim = (/False,False/)
         filedimdef(fout,dimNames,dimSizes,dimUnlim)

         ;============================================================
         var_name = variable
         var_name_std=variable+"_std"

         filevardef(fout,"plev",typeof(plev),"plev")
         filevardef(fout,var_name,typeof(v_wgt),(/"plev","month"/))
         filevardef(fout,var_name_std,typeof(v_wgt),(/"plev","month"/))
         setfileoption(fout,"DefineMode",False)

         fout->plev       = (/plev/)
         fout->$var_name$ = (/v_wgt(case|0,plev|:,month|:)/)
         fout->$var_name_std$ = (/v_wgt(case|1,plev|:,month|:)/)
         print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
         print("!")
         print("!  Files are been produced in "+diri_out+"/E06FIG08/  !!")
         print("!")
         print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      else
         print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
      end if

      ;===============================================================
      res                      = True                  ; plot mods desired
      res@tiMainString         = models@name(imod)     ; title
      res@gsnDraw              = False
      res@gsnFrame             = False
      res@lbLabelBarOn         = False
      res@cnLineLabelsOn       = False
      res@tiYAxisFont          = "helvetica"
      res@tiMainFont           = "helvetica"
      res@tmXBLabelFont        = "helvetica"
      res@tmYLLabelFont        = "helvetica"
      res@lbLabelFont          = "helvetica"
      res@cnLevelSelectionMode = "ManualLevels"        ; manual contour levels
      res@cnLevelSpacingF      = 0.1                   ; contour interval
      res@cnMinLevelValF       = -1.0                  ; min level
      res@cnMaxLevelValF       = 1.0                   ; max level
      res@lbLabelStrings       =(/"-1","-0.8","-0.6","-0.4","-0.2","0.","0.2","0.4","0.6","0.8","1."/)
      res@cnFillOn             = True                  ; turn on color fill
      res@cnInfoLabelOn        = False                 ; turn off contour label
      res@gsnSpreadColors      = True                  ; use full range of colors
      res@gsnSpreadColorStart  = 10
      res@lbLabelAutoStride    = True                  ; optimal labels
      res@gsnLeftString        = ""
      res@gsnRightString       = ""
      res@gsnLeftStringFontHeightF = 0.015
      res@tmXBMode             = "Explicit"            ; explicit labels
      ; location of labels
      res@gsnAddCyclic         = True
      res@tmXBValues = (/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23/)
      res@tmXBLabels = (/"J","F","M","A","M","J","J","A","S","O","N","D","J","J","F","M","A","M","J","J","A","S","O","N","D","J"/)                ; labels themselves

      ; res@tmYLValues = (/100,70,50,30,20,10/)
      ; res@tmYLLabels = (/"100","70","50","30","20","10"/)
      res@trYMaxF              = 100
      res@trYMinF              = 10
      ;*** COMMENT_OUT FOR NCL6.0.0 ***  res@trYAxisType          = "LogAxis"
      res@vpHeightF            = 0.4                    ; change aspect ratio of plot
      res@vpWidthF             = 0.75
      res@lbLabelFontHeightF   = 0.011
      res@tiYAxisString        = "Pressure ["+plev@units+"]"    ; y-axis label

      plot(imod+dim_OBS) = gsn_csm_pres_hgt (wks,v_wgt(case|0,plev|:,month|:),res)
      delete(plev)
      delete(v_wgt)
      delete(res)
   end do

   Pres = True
   Pres@gsnPanelLabelBar  = True
   Pres@lbLabelAutoStride = True
   if (isatt(info,"long_name")) then
      Pres@txString      = info@long_name+" anomaly time height sections (averaged from "+lat_name_min+" to "+lat_name_max+")"
      Pres@txFontHeightF = 0.015
   end if
   n_p = (dim_MOD+dim_OBS)/3
   if (mod((dim_MOD+dim_OBS),3) .ne. 0) then
      n_p = n_p+1
   end if
   gsn_panel(wks,plot,(/n_p,3/),Pres)

   print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   print("!")
   print("!  Figure E06FIG08.ps is produced in "+diri_out)
   print("!")
   print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

   print ("=================== FINISH (E06FIG08) ====================")
   delete(lat_min)
   delete(dim_MOD)
   delete(lat_name_max)
   delete(lat_name_min)
   delete(dim_OBS)
   delete(year_max)
   delete(Pres)
   delete(i)
   delete(lat_max)
   delete(imod)
   delete(models)
   delete(year_min)
   delete(n_p)

end
;##################################################################################
