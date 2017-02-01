;##################################################################################
;# This NCL script reads in CF compliant netCDF CCMVal model output as specified  #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and produces FIGURE 9 of                                                       #
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
;# Description: This script is designed to plot the vertical variation of         #
;# (a) amplitude and (b) phase lag of the annual cycle of water vapor averaged    #
;# between lat_min and lat_max. The amplitude is normalized to unity and the      #
;# phase lag is set to zero at the level where the amplitude is maximum (between  #
;# 16 and 20 km), which varies between the CCMs. The vertical axis in both plots  #
;# is the distance from the level of maximum amplitude. Variables plotted can     #
;# also be output to netcdf files.                                                #
;##################################################################################
; This script takes input from the /var_att/H2O_att.ncl file
; info attributes required (are the same as in E06FIG08):
; fig09_yr_min          start year of plot
; fig09_yr_max          end year of plot
; fig09_lat_min         start latitude for spatial averaging
; fig09_lat_max         end latitude for spatial averaging
; fig09_obs_name        (optional) observation name
; fig09_obs_file        (optional) file name of file containing observation
;                       data for plot
;
;**********************************************************************************
load "./plot_type/CCMVal_FUNCTION/new_function.ncl"
load "./plot_type/CCMVal_FUNCTION/Climatology_CCMVal.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
load "./plot_type/CCMVal_FUNCTION/proc_E06FIG08.ncl"

begin
;  fields supported T3M T3D T2Mz T2Dz
   print("Begin ./plot_type/E06FIG09.ncl")
   diri_out = getenv("ESMValTool_plot_dir")
   project_name = getenv("ESMValTool_project")
   OPTIONS_FILE = getenv("ESMValTool_write_plot_vars")

   ;******************************************************************
   ;check the input parameters & translate them
   year_min=stringtointeger(models@start_date)
   year_max=stringtointeger(models@end_date)

   if(isatt(info,"fig08_lat_min") .and. isatt(info,"fig08_lat_max")) then
      lat_min = info@fig08_lat_min
      lat_max = info@fig08_lat_max
      lat_name_max = lat_names(lat_max)
      lat_name_min = lat_names(lat_min)
   else
      print("fatal:plot requires variable attributes 'fig08_lat_min' and 'fig08_lat_max'")
      status_exit(1)
   end if

   ; optional input parameters
   if(isatt(info,"fig08_obs_name") .and. isatt(info,"fig08_obs_file")) then
      OBS = info@fig08_obs_name
      obs_filename = info@fig08_obs_file
   else
      OBS = (/"NaN"/)
   end if

   ;******************************************************************
   dim_OBS = dimsizes(OBS)
   dim_MOD = dimsizes(models@name)
   if (OBS.ne."NaN") then
      Obs_mod = new(dim_OBS+dim_MOD,typeof(OBS))
      Obs_mod(:dim_OBS-1) = OBS
      Obs_mod(dim_OBS:) = models@name
   else
      dim_OBS = 0
      Obs_mod = new(dim_MOD,typeof(OBS))
      Obs_mod = models@name
   end if

   if (project_name.eq."CCMVal2") then
      line_color = CCMVal2_ColorTable(models@name)
      line_dash = CCMVal2_DashTable(models@name)
   else
      line_color = CCMVal1_ColorTable(models@name)
      line_dash = CCMVal1_DashTable(models@name)  ; line_dash = CCMVal_Dash(models@name)
   end if

   if (project_name.eq."CCMVal2") then
      line_colors = CCMVal2_ColorTable(Obs_mod)
      line_dashs = CCMVal2_DashTable(Obs_mod)
   else
      line_colors = CCMVal1_ColorTable(Obs_mod)
      line_dashs = CCMVal1_DashTable(Obs_mod)    ; line_dashs = CCMVal_Dash( Obs_mod)
   end if

   wks = gsn_open_wks ("ps",diri_out+"/E06FIG09")

   plot = new((/2/),"graphic")

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
         delete(A0)
      end do
   else
      A0 = extract_data(imod,data,-1,0,0)
      print("processing "+models@name(imod)+" ("+(imod+1)+" of " +dim_MOD+")")
      v_wgt = figure8(A0,field_number,lat_max,lat_min,year_max(imod),year_min(imod))

      plev = v_wgt&plev
      f_obs = addfile(obs_filename(iobs),"r")
      delete(A0)
      v_wgt = figure8(A0,field_number,lat_max,lat_min,year_max(imod),year_min(imod))

      delete(A0)
   end if
   plev = v_wgt&plev
   plev!0="plev"
   plev&plev=plev

   ;============================================================
   vvgt=new((/dimsizes(plev),24/),typeof(v_wgt))

   vvgt(:,0:23)=v_wgt(case|0,{plev|:},month|:)
   vvgt!0="plev"
   plev=vvgt&plev
   vvgt!1="month"

   vvgt=where(ismissing(vvgt),0.,vvgt)
   ffa=fourier_info(vvgt(plev|:,month|0:11),1,1.0)   ;fourier analysis max

   ffa!1="plev"
   ffa&plev=plev

   ind_pl=ind((ffa&plev.gt.96).and.(ffa&plev.le.100))
   ffa&plev=plev
   fase=ffa(1,:,0)
   zlevel =(16*log10(10^3/plev))
   ;max_ind=ind(ffa(0,:,0).eq.(max(ffa(0,:,0))))
   max_ind=ind(ffa(0,:,0).eq.(max(ffa(0,ind_pl,0))))
   ;min_ind=ind(ffa(0,:,0).eq.(min(ffa(0,ind_pl,0))))

   ffa_m=fourier_info(-vvgt(plev|:,month|0:11),1,1.0) ;fourier analysis min
   ffa_m!1="plev"
   ffa_m&plev=plev
   fase_m=ffa_m(1,:,0)
   zlev=  (16*log10(10^3/plev))
   max_ind_m=ind(ffa(0,:,0).eq.(max(ffa_m(0,:,0))))
   ;min_ind_m=ind(ffa(0,:,0).eq.(min(ffa_m(0,ind_pl,0))))
   delete(ind_pl)

   amp_rel=ffa(0,max_ind:,0)/(ffa(0,max_ind,0))

   z_rel=zlev(max_ind:)-zlev(max_ind)

   ph=fase(max_ind:)-fase(max_ind)   ;phase max
   ph_m=fase_m(max_ind_m:)-fase_m(max_ind_m) ;phase min

   a=(/dimsizes(ph),dimsizes(ph_m)/)

   dime=min(a)
    do i=1,dime-1
   ;do i=1,dimsizes(ph)-1
      if ((ph(i).le.ph(i-1))).and.(abs(ph(i)-ph(i-1)).ge.abs(ph(i)-ph(i-1)+12)) then
         ph(i)=12+ph(i)
      end if

      if ((ph_m(i).le.ph_m(i-1))).and.(abs(ph_m(i)-ph_m(i-1)).ge.abs(ph_m(i)-ph_m(i-1)+12)) then
         ph_m(i)=12+ph_m(i)
      end if

   end do

   phase=(ph(:dime-1)+ph_m(:dime-1))/2

   res                      = True                  ; plot mods desired
   ;   res@tiMainString         = models@name(imod)     ; title
   res@gsnDraw              = False
   res@gsnFrame             = False

   res@trYMaxF              = 12
   res@trYMinF              = 0

   res@vpHeightF            = 0.7                    ; change aspect ratio of plot
   res@vpWidthF             = 0.4

   res@tiYAxisString        = "z-z0 [Km]"    ; y-axis label
   if (dim_OBS.ne.0) then
      res@xyMarkLineModes   = "Markers"                ; choose which have markers
      res@xyMarkers         =  16                      ; choose type of marker
      res@xyMarkerColor     = "Gray0"                  ; Marker color
      res@xyMarkerSizeF     = 0.01                     ; Marker size (default 0.01)
   else
      res@xyLineColors     = line_color(0)              ; line color
      res@xyDashPatterns   = line_dash(0)
   end if
   res@trXMaxF              = 1.2
   res@trXMinF              =0.
   res@tiXAxisString        = "Relative Amplitude"
   plot(0) = gsn_csm_xy (wks,amp_rel,z_rel,res)
   res@trXMaxF              = 18
   res@trXMinF              = -2
   res@tiXAxisString        = "Phase Lag (months)"
   plot(1) =gsn_csm_xy (wks,phase,z_rel,res)
   if (OPTIONS_FILE.eq."True") then
      fili = "/E06FIG09/E06FIG09_"+project_name+"_"+OBS+"_"+lat_name_max+lat_name_min+"_"+variable+".nc"
      f_out = diri_out+"/"+fili
      print("output file: "+f_out)
      system("/bin/rm -f "+f_out)
      fout = addfile(f_out,"c")
      setfileoption(fout,"DefineMode",True)

      ;============================================================
      ; create global attributes of the file
      ;============================================================
      fAtt               = True            ; assign file attributes
      fAtt@title         = "Climatological data (plev,month)"
      ; fAtt@source_file = f
      fAtt@creation_date = systemfunc("date")
      fAtt@model= OBS

      ; copy attributes from the source file
      fileattdef(fout,fAtt)            ; copy file attributes

      dimNames = (/"z-z0"/)
      dimSizes = (/dimsizes(z_rel)/)
      dimUnlim = (/False/)
      filedimdef(fout,dimNames,dimSizes,dimUnlim)

      ;============================================================
      var_name = variable

      filevardef(fout,"RelativeAmplitude",typeof(z_rel),"z-z0")
      filevardef(fout,"Phase_Lag",typeof(phase),(/"z-z0"/))
      filevardef(fout,"z_z0",typeof(z_rel),(/"z-z0"/))
      setfileoption(fout,"DefineMode",False)

      fout->RelativeAmplitude=(/amp_rel/)
      fout->Phase_Lag = (/phase/)
      fout->z_z0 = (/z_rel/)
      print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("!")
      print("!  Files are been produced in "+diri_out+"/E06FIG09/  !!")
      print("!")
      delete(fAtt)
      delete(fout)
      delete(fili)
      print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   else
      print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
   end if

   delete(plev)
   delete(v_wgt)
   delete(vvgt)
   delete(res)
   delete(amp_rel)
   delete(z_rel)
   delete(ffa)
   delete(ffa_m)
   delete(max_ind)

   delete(ph)
   delete(ph_m)
   delete(zlev)
   delete(phase)
   delete(fase)
   delete(fase_m)

;*****************************************************************
;MODEL

   do imod=0,dim_MOD-1

      A0 = extract_data(imod,data,-1,0,0)
      print("processing "+models@name(imod)+" ("+(imod+1)+" of " +dim_MOD+")")
      v_wgt = figure8(A0,field_number,lat_max,lat_min,year_max(imod),year_min(imod))
      plev = v_wgt&plev

      plev!0="plev"
      plev&plev=plev
      vvgt=new((/dimsizes(plev),24/),typeof(v_wgt))
      vvgt(:,0:23)=v_wgt(case|0,{plev|:},month|:)
      vvgt!0="plev"
      plev=(/vvgt&plev/)
      vvgt!1="month"
      vvgt=where(ismissing(vvgt),0.,vvgt)

      ffa=fourier_info(vvgt(plev|:,month|0:11),1,1.0)   ;fourier analysis max
      ffa!1="plev"
      ffa&plev=plev
      ind_pl=ind((ffa&plev.gt.96).and.(ffa&plev.le.101))
      ffa&plev=plev
      fase=ffa(1,:,0)

      zlev =(16*log10(10^3/plev))
      ;max_ind=ind(ffa(0,:,0).eq.(max(ffa(0,ind_pl,0))))
      ; min_ind=ind(ffa(0,:,0).eq.(min(ffa(0,ind_pl,0))))
      max_ind=ind(ffa(0,:,0).eq.(max(ffa(0,ind_pl(0):,0))))

      min_ind=ind(ffa(0,:,0).eq.(min(ffa(0,ind_pl(0):,0))))
      ;max_ind=maxind(ffa(0,ind_pl(::-1),0))
      ;min_ind=minind(ffa(0,ind_pl(::-1),0))

      ffa_m=fourier_info(-vvgt(plev|:,month|0:11),1,1.0) ;fourier analysis min
      ffa_m!1="plev"
      ffa_m&plev=plev
      fase_m=ffa_m(1,:,0)

      ;max_ind_m=maxind(ffa_m(0,ind_pl(::-1),0))
      ;min_ind_m=minind(ffa_m(0,ind_pl(::-1),0))

      max_ind_m=ind(ffa(0,:,0).eq.(max(ffa_m(0,ind_pl(0):,0))))
      min_ind_m=ind(ffa(0,:,0).eq.(min(ffa_m(0,ind_pl(0):,0))))
      delete(ind_pl)
      amp_rel=ffa(0,max_ind_m:,0)/(ffa(0,max_ind_m,0))
      z_rel=zlev(max_ind_m:)-zlev(max_ind_m)

      ph=fase(max_ind:)-fase(max_ind)
      ph_m=fase_m(max_ind_m:)-fase_m(max_ind_m)

      a=(/dimsizes(ph),dimsizes(ph_m)/)

      dime=min(a)

      do i=1,dime-1
      ;do i=1,dimsizes(ph)-1
         if ((ph(i).lt.ph(i-1))).and.(abs(ph(i)-ph(i-1)).gt.abs(ph(i)-ph(i-1)+12)) then
            ph(i)=12+ph(i)
         end if
         if ((ph(i).lt.ph(i-1))).and.(abs(ph(i)-ph(i-1)).gt.abs(ph(i)-ph(i-1)+12)) then
            ph(i)=12+ph(i)
         end if
         if ((ph_m(i).lt.ph_m(i-1))).and.(abs(ph_m(i)-ph_m(i-1)).gt.abs(ph_m(i)-ph_m(i-1)+12)) then
            ph_m(i)=12+ph_m(i)
         end if
         if ((ph_m(i).lt.ph_m(i-1))).and.(abs(ph_m(i)-ph_m(i-1)).gt.abs(ph_m(i)-ph_m(i-1)+12)) then
            ph_m(i)=12+ph_m(i)
         end if
      end do

      phase=(ph(:dime-1)+ph_m(:dime-1))/2
      ;===============================================================

      add_line(wks,plot(0),amp_rel,z_rel,line_color(imod),line_dash(imod))
      add_line(wks,plot(1),phase,z_rel,line_color(imod),line_dash(imod))

      if (OPTIONS_FILE.eq."True") then
         fili = "/E06FIG09/E06FIG09_"+project_name+"_"+Obs_mod(imod)+"_"+lat_name_max+lat_name_min+"_"+variable+".nc"
         f_out = diri_out+"/"+fili
         print("Model ["+Obs_mod(imod)+"]: "+f_out)
         system("/bin/rm -f "+f_out)
         fout = addfile(f_out,"c")
         setfileoption(fout,"DefineMode",True)

         ;============================================================
         ; create global attributes of the file
         ;============================================================
         fAtt               = True            ; assign file attributes
         fAtt@title         = "Climatological data (plev,month)"
         ; fAtt@source_file = f
         fAtt@creation_date = systemfunc("date")
         fAtt@model= Obs_mod(imod)
         setfileoption(fout,"DefineMode",True)

         ; copy attributes from the source file
         fileattdef(fout,fAtt)            ; copy file attributes

         dimNames = (/"z-z0"/)
         dimSizes = (/dimsizes(z_rel)/)
         dimUnlim = (/False/)
         filedimdef(fout,dimNames,dimSizes,dimUnlim)

         ;============================================================
         var_name = variable

         filevardef(fout,"RelativeAmplitude",typeof(z_rel),"z-z0")
         filevardef(fout,"Phase_Lag",typeof(phase),(/"z-z0"/))
         filevardef(fout,"z_z0",typeof(z_rel),(/"z-z0"/))
         setfileoption(fout,"DefineMode",False)

         fout->RelativeAmplitude=(/amp_rel/)
         fout->Phase_Lag = (/phase/)
         fout->z_z0 = (/z_rel/)
         print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
         print("!")
         print("!  Files are been produced in "+diri_out+"/E06FIG09/  !!")
         print("!")
         print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
         delete(fAtt)
      else
         print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
      end if

      delete(plev)
      delete(v_wgt)
      delete(vvgt)
      delete(amp_rel)
      delete(z_rel)
      delete(ffa)
      delete(ffa_m)
      delete(max_ind)
      delete(ph)
      delete(ph_m)
      delete(zlev)
      delete(phase)
      delete(fase)
      delete(fase_m)
      delete(A0)
   end do

   Pres                 = True
   Pres@gsnFrame        = False
   Pres@gsnPanelBottom  = 0.33
   Pres@gsnPanelXWhiteSpacePercent = 10

   gsn_panel(wks,plot,(/1,2/),Pres)
   add_legenda(wks,Obs_mod,line_colors,line_dashs,"bottom")

   print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   print("!")
   print("!  Figure E06FIG09.ps is produced in "+diri_out)
   print("!")
   print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   print ("=================== FINISH (E06FIG09) ====================")
end
;##################################################################################
