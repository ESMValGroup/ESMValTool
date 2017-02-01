;##################################################################################
;# This NCL script reads in CF compliant netCDF CCMVal model output as specified  #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and produces FIGURE 14 of                                                      #
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
;# Description: This script is designed to plot modeled zonal mean total column   #
;# ozone climatologies in comparison to observations. Variables plotted can       #
;# also be output to netcdf files.                                                #
;##################################################################################
;# INPUT/OUTPUT:                                                                  #
;# 1. Reads in                                                                    #
;#    (a) '${WKDIR}/var_att/toz_att.ncl'  (${WKDIR} specified in main namelist)   #
;#    (b) '${WKDIR}/namelist_EO6'         (specification of model names etc.)     #
;#    (c) '${WKDIR}/diag_att/E06Diag.att' (specifies plot type)                   #
;#    (d) Model input files               (as defined in namelists)               #
;# 2. produces                                                                    #
;#    (a) '${diri_out}/E06FIG14.ps'       (plot as ps-file)                       #
;#    (b) '${diri_out}/E06Fig14/'         (plot variables as CF netCDF file)      #
;##################################################################################
; This script takes input from the /var_att/toz_att.ncl file
; no info attributes required
;
;**********************************************************************************
load "./plot_type/CCMVal_FUNCTION/new_function.ncl"
load "./plot_type/CCMVal_FUNCTION/Climatology_CCMVal.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
load "./plot_type/CCMVal_FUNCTION/proc_E06FIG14.ncl"
load "./plot_type/get_field_number.ncl"

begin
; fields supported T2Ms T2Ds
   print("Begin ./plot_type/E06FIG14.ncl")
   diri_out = getenv("ESMValTool_plot_dir")
   project_name = getenv("ESMValTool_project")
   OPTIONS_FILE = getenv("ESMValTool_write_plot_vars")

   ;******************************************************************
   ;check the input parameters
   year_min=stringtointeger(models@start_date)
   year_max=stringtointeger(models@end_date)

   ; optional input parameters
   if(isatt(info,"fig14_obs_name") .and. isatt(info,"fig14_obs_file")) then
      OBS = info@fig14_obs_name
      obs_filename = info@fig14_obs_file
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

   if (project_name.eq."CCMVal2") then
      colors = CCMVal2_ColorTable(Obs_mod)
      dashs = CCMVal2_DashTable(Obs_mod)
   else
      colors = CCMVal1_ColorTable(Obs_mod)
      dashs = CCMVal1_DashTable(Obs_mod)
   end if

   wks = gsn_open_wks ("ps",diri_out+"/E06FIG14_"+variable)

   ; open workstation
   cmap = (/"White","Gray1","Gray75","maroon4","DarkOrchid4","RoyalBlue4", "RoyalBlue","CadetBlue3", "CadetBlue1","DarkOliveGreen3", "Chartreuse3", \
           "GreenYellow", "Yellow","Orange","OrangeRed","OrangeRed4",\
           "Black"/)
   gsn_define_colormap(wks,cmap)

   plot = new((/dim_OBS+dim_MOD/),"graphic")

;******************************************************************
;OBS

   if (dim_OBS.ne.0) then
      do iobs=0,dim_OBS-1
         f_obs = addfile(obs_filename(iobs),"r")
         A0 = f_obs->$variable$      ;Observation by HALOE are already processed:mean
         A0!0 = "month"
         A0!1 = "lat"
         ;A0 = extract_data(obs,data,iobs)
         print("processing "+OBS(iobs)+" ("+(iobs+1)+" of " +dim_OBS+")")

         v_wgt=A0(month|:11,lat|:)
         if (OBS(iobs).eq."NIWA") then
            v_wgt=where(v_wgt.ne.0,v_wgt,v_wgt@_FillValue)
         end if

         uniform_coordinates(v_wgt,"T1Ms")
         lat=v_wgt&lat
         lat = coord_list("lat",lat)
         v_wgt&lat=lat
         delete(A0)
         ;A0 = extract_data(obs,data,iobs)

         if (OPTIONS_FILE.eq."True") then
            fili = (project_name+"_"+OBS(iobs)+"_"+variable+".nc")
            f_out = diri_out+"/E06FIG14/E06FIG14_"+fili
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
            fAtt@model         = OBS(iobs)
            ; copy attributes from the source file
            fileattdef(fout,fAtt)            ; copy file attributes

            dimNames = (/"lat","month"/)
            dimSizes = (/dimsizes(v_wgt&lat),dimsizes(v_wgt&month)/)
            dimUnlim = (/False,False/)
            filedimdef(fout,dimNames,dimSizes,dimUnlim)

            ;============================================================
            var_name = variable

            filevardef(fout,"lat",typeof(v_wgt&lat),"lat")
            filevardef(fout,var_name,typeof(v_wgt),(/"lat","month"/))

            setfileoption(fout,"DefineMode",False)

            fout->lat        = (/v_wgt&lat/)
            fout->$var_name$ = (/v_wgt(lat|:,month|:)/)
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print("!")
            print("!  Files "+f_out+" are been produced in "+diri_out+"/E06FIG14/  !!")
            print("!")
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
         else
            print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
         end if

         ;===============================================================
         res                      = True                  ; plot mods desired
         res@gsnDraw              = False
         res@gsnFrame             = False
         res@tiMainString         = OBS(iobs)             ; title
         res@lbLabelBarOn         = False
         res@tiYAxisFont          = "helvetica"
         res@tiMainFont           = "helvetica"
         res@tmXBLabelFont        = "helvetica"
         res@tmYLLabelFont        = "helvetica"
         res@lbLabelFont          = "helvetica"
         res@cnLevelSelectionMode = "ManualLevels"        ; manual contour levels
         res@cnLevelSpacingF      = 25                    ; contour interval
         res@cnMinLevelValF       = 150                   ; min level
         res@cnMaxLevelValF       = 475                   ; max level
         ;res@cnLevelSpacingF      = 10                    ; contour interval
         ;res@cnMinLevelValF       = 0                   ; min level
         ;res@cnMaxLevelValF       = 150                   ; max level
         res@cnLineLabelsOn       = False                 ; turn off line labels
         res@cnFillOn             = True                  ; turn on color fill
         res@cnInfoLabelOn        = False                 ; turn off contour label
         res@gsnSpreadColors      = True                  ; use full range of colors
         res@lbLabelAutoStride    = True                  ; optimal labels
         res@gsnLeftString        = ""
         res@gsnRightString       = ""
         res@gsnLeftStringFontHeightF = 0.015
         res@tmXBMode             = "Explicit"            ; explicit labels
         ; location of labels
         res@gsnAddCyclic         = True
         res@tmXBValues           = (/0,1,2,3,4,5,6,7,8,9,10,11,12/)
         res@tmXBLabels           = (/"J","F","M","A","M","J","J","A","S","O","N","D","J"/)   ; labels themselves
         res@tmYLValues           = (/-90,-60,-30,0,30,60,90/)
         res@tmYLLabels           =(/"-90","-60","-30","0","30","60","90"/)
         res@vpHeightF            = 0.4                   ; change aspect ratio of plot
         res@vpWidthF             = 0.75
         res@lbLabelFontHeightF   = 0.011
         res@tiYAxisString        = "Latitude"            ; y-axis label

         plot(iobs) = gsn_csm_lat_time (wks,v_wgt(lat|:,month|:),res)
         delete(v_wgt)
         delete(lat)
         delete(res)
      end do
   end if

;******************************************************************
;MODEL

   do imod=0,dim_MOD-1

      A0 = extract_data(imod,data,-1,0,0)
      print("Model: "+models@name(imod))
      v_wgt = figure14(A0,field_number,models@name(imod),year_min(imod),year_max(imod))
      delete(A0)
      if (OPTIONS_FILE.eq."True") then

         fili = (project_name+"_"+models@case_name(imod)+"_"+models@name(imod)+"_"+models@ensemble(imod)+"_"+field_number+"_"+variable+".nc")
         f_out = diri_out+"/E06FIG14/E06FIG14_"+year_max(imod)+"_"+year_min(imod)+fili
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
         fAtt@model         = models@name(imod)
         ; copy attributes from the source file
         fileattdef(fout,fAtt)            ; copy file attributes

         dimNames = (/"lat","month"/)
         dimSizes = (/dimsizes(v_wgt&lat),dimsizes(v_wgt&month)/)
         dimUnlim = (/False,False/)
         filedimdef(fout,dimNames,dimSizes,dimUnlim)

         ;============================================================
         var_name = variable

         filevardef(fout,"lat",typeof(v_wgt&lat),"lat")
         filevardef(fout,var_name,typeof(v_wgt),(/"lat","month"/))
         setfileoption(fout,"DefineMode",False)

         fout->lat        = (/v_wgt&lat/)
         fout->$var_name$ = (/v_wgt(lat|:,month|:)/)
         print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
         print("!")
         print("!  Files "+f_out+" are been produced in "+diri_out+"/E06FIG14/  !!")
         print("!")
         print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      else
         print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
      end if

      ;===============================================================
      res                      = True                  ; plot mods desired
      res@gsnDraw              = False
      res@gsnFrame             = False
      res@lbLabelBarOn         = False
      res@tiMainString         = models@name(imod)+" "+year_min(imod)+"_"+year_max(imod)     ; title
      ;res@tiMainString         = "ACC&SPARC "+year_min(imod)+"_"+year_max(imod)
      res@tiYAxisFont          = "helvetica"
      res@tiMainFont           = "helvetica"
      res@tmXBLabelFont        = "helvetica"
      res@tmYLLabelFont        = "helvetica"
      res@lbLabelFont          = "helvetica"
      res@cnLevelSelectionMode = "ManualLevels"        ; manual contour levels
      res@cnLevelSpacingF      = 25                    ; contour interval
      res@cnMinLevelValF       = 150                   ; min level
      res@cnMaxLevelValF       = 475                   ; max level
      ;res@cnLevelSpacingF      = 3                  ; contour interval
      ;res@cnMinLevelValF       = 10                   ; min level
      ;res@cnMaxLevelValF       = 46                   ; max level
      res@cnLineLabelsOn       = False                 ; turn off line labels
      res@cnFillOn             = True                  ; turn on color fill
      res@cnInfoLabelOn        = False                 ; turn off contour label
      res@gsnSpreadColors      = True                  ; use full range of colors
      res@lbLabelAutoStride    = True                  ; optimal labels
      res@gsnLeftString        = ""
      res@gsnRightString       = ""
      res@gsnLeftStringFontHeightF = 0.015
      res@tmXBMode             = "Explicit"            ; explicit labels
      ; location of labels
      res@gsnAddCyclic=True
      res@tmXBValues           = (/0,1,2,3,4,5,6,7,8,9,10,11,12/)
      res@tmXBLabels           = (/"J","F","M","A","M","J","J","A","S","O","N","D","J"/)    ; labels themselves
      res@tmYLValues           = (/-90,-60,-30,0,30,60,90/)
      res@tmYLLabels           = (/"-90","-60","-30","0","30","60","90"/)
      res@vpHeightF            = 0.4                   ; change aspect ratio of plot
      res@vpWidthF             = 0.75
      res@lbLabelFontHeightF   = 0.011
      res@tiYAxisString        = "Latitude"            ; y-axis label

      plot(imod+dim_OBS) = gsn_csm_lat_time (wks,v_wgt(lat|:,month|:),res)
      delete(v_wgt)
      delete(res)
   end do

   Pres = True
   if (isatt(info,"long_name")) then
      Pres@txString = info@long_name
   end if
   Pres@gsnPanelLabelBar = True
   Pres@lbLabelAutoStride    = True
   n_p = (dim_MOD+dim_OBS)/3
   if (mod((dim_MOD+dim_OBS),3) .ne. 0) then
      n_p = n_p+1
   end if
   gsn_panel(wks,plot,(/n_p,3/),Pres)
   print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   print("!")
   print("!  Figure E06FIG14.ps is produced in "+diri_out)
   print("!")
   print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   delete(cmap)
   delete(dim_MOD)
   delete(year_max)
   delete(dashs)
   delete(colors)
   delete(Pres)
   delete(i)
   delete(imod)
   delete(models)
   delete(year_min)
   delete(n_p)
   delete(Obs_mod)

   print ("=================== FINISH (E06FIG14) ====================")

end
;##################################################################################
