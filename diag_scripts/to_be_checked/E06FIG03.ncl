;##################################################################################
;# This NCL script reads in CF compliant netCDF CCMVal model output as specified  #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and produces FIGURE 3 of                                                       #
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
;# Description: This script is designed to plot the annual mean of heat fluxes at #
;# 100 hPa versus temperatures. Shown are all the years in the selected interval  #
;# for each model simulation compared to observations.                            #
;# Variables plotted can also be output to netcdf files.                          #
;#  The number of panels in the plot is determined by the number of latitude      #
;#  ranges selected in the info array (info@fig03_lat_min ...)                    #
;##################################################################################
;# INPUT/OUTPUT:                                                                  #
;# 1. Reads in                                                                    #
;#    (a) '${WKDIR}/var_att/ta_att.ncl'   (${WKDIR} specified in main namelist)   #
;#    (b) '${WKDIR}/namelist_EO6'         (specification of model names etc.)     #
;#    (c) '${WKDIR}/diag_att/E06Diag.att' (specifies plot type)                   #
;#    (d) Model input files               (as defined in namelists)               #
;# 2. produces                                                                    #
;#    (a) '${diri_out}/E06FIG03.ps'       (plot as ps-file)                       #
;#    (b) '${diri_out}/E06FIG03/E06FIG03_*.nc (plot variables as CF netCDF file)  #
;##################################################################################
; This script takes input from the ta_att.ncl file
; info attributes required:
; fig03_lev_sel         selected pressure level
; fig03_lat_min         array of latitude (start) values for area averaging
; fig03_lat_max         array of latitude (end) values for area averaging
; fig03_month_sel       array of month for averaging.  Note: this array must have the same
;                       dimensions as fig03_lat_min and fig03_lat_max
; fig03_vs_var ="vt100"
; fig03_vs_fieldModels="T1ms"
; fig03_vs_filedir
; fig03_vs_lev_sel       selected pressure level
; fig03_vs_lat_min       array of latitude (start) values for area averaging
; fig03_vs_lat_max       array of latitude (end) values for area averaging
; fig03_vs_month_sel     array of month for averaging.  Note: this array must have the same
;                       dimensions as fig03_lat_min and fig03_lat_max

;**********************************************************************************
load "./plot_type/CCMVal_FUNCTION/new_function.ncl"
load "./plot_type/CCMVal_FUNCTION/Climatology_CCMVal.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
load "./plot_type/CCMVal_FUNCTION/proc_E06FIG03.ncl"

begin
; fields supported T3M T3D T2Mz T2Dz
   print("Begin ./plot_type/E06FIG03.ncl")
   diri_out = getenv("ESMValTool_plot_dir")
   project_name = getenv("ESMValTool_project")
   OPTIONS_FILE = getenv("ESMValTool_write_plot_vars")

   MODEL=models@name
   dim_MOD = dimsizes(models@name)-1

   year_min=stringtointeger(models@start_date)
   year_max=stringtointeger(models@end_date)

   ;******************************************************************

   if (isatt(info,"fig03_vs_var")) then
      var1=info@fig03_vs_var
   else
      print("fatal:plot requires variable attribute 'vs_var'.")
      print("fatal:In /var_att/ta_att.ncl info_03@vs_var must point to a variable in the same folders as indicated in the namelist_E06")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if
   if(isvar("month_sele")) then
      delete(month_sele)
   end if
   if(isvar("lev_sel")) then
      delete(lat_min)
   end if
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

   if(isatt(info,"fig03_lat_max") .and. isatt(info,"fig03_lat_min").and.isatt(info,"fig03_lev_sel")) then
      lat_max = info@fig03_lat_max
      lat_min = info@fig03_lat_min
      lat_name_max = lat_names(lat_max)
      lat_name_min = lat_names(lat_min)
      lev_sel= info@fig03_lev_sel
   else
      print("fatal:plot requires variable attributes 'fig03_lat_max' and 'fig03_lat_min' and 'fig03_lev_sel'")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if

   if(isatt(info,"fig03_month_sel")) then
      month_sele = info@fig03_month_sel
   else
      print("fatal:plot requires variable attribute 'fig03_month_sel'")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if

   if ((dimsizes(lat_max).ne.dimsizes(lat_min)) .or. (dimsizes(lat_max).ne.dimsizes(month_sele(:,0)))) then
      print("fatal:Error!!!!The number of the month must to be the same of the number of lat_max & lat_min!!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if

   if(isatt(info,"fig03_vs_lat_max") .and. isatt(info,"fig03_vs_lat_min").and.isatt(info,"fig03_vs_lev_sel").and.isatt(info,"fig03_vs_month_sel")) then
      lat_max_vs = info@fig03_vs_lat_max
      lat_min_vs = info@fig03_vs_lat_min
      lat_name_max_vs = lat_names(lat_max_vs)
      lat_name_min_vs = lat_names(lat_min_vs)
      lev_sel_vs= info@fig03_vs_lev_sel
      month_sel_vs=info@fig03_vs_month_sel
  else
      print("fatal:plot requires variable attributes 'fig03_vs_lat_max' and 'fig03_vs_lat_min' and 'fig03_vs_lev_sel'")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if

   if ((dimsizes(lat_max_vs).ne.dimsizes(lat_min_vs)) .or. (dimsizes(lat_max_vs).ne.dimsizes(month_sel_vs(:,0)))) then
      print("fatal:Error!!!!The number of the month must to be the same of the number of lat_max & lat_min!!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if

   if(isatt(info,"fig03_vs_var") .and. isatt(info,"fig03_vs_field")) then
      var_vs = info@fig03_vs_var
      field_vs = info@fig03_vs_field
   else
      print("fatal:plot requires variable attributes 'fig03_vs_field' and 'fig03_vs_var'")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if

   if (isatt(info,"fig03_obs_file").and.isatt(info,"fig03_obs_file_vs").and.isatt(info,"fig03_obs")) then
       dim_OBS=dimsizes(info@fig03_obs)
   else
       dim_OBS=0
   end if

;***************************************begin*******************************************
   nyear=min(year_max)-max(year_min)+1
   MODELOBS=new(dimsizes(MODEL)+dim_OBS,typeof(MODEL))

   MODELOBS(0:dimsizes(MODEL)-1)=MODEL
   MODELOBS(dimsizes(MODEL):)=info@fig03_obs
   if (project_name.eq."CCMVal2") then
      colors = CCMVal2_ColorTable(MODELOBS)
      dash = CCMVal2_DashTable(MODELOBS)
   else
      colors = CCMVal1_ColorTable(MODELOBS)
      dash = CCMVal1_DashTable(MODELOBS)
   end if
;*****************************************************************************************

   v_wgt_sel1=new((/dimsizes(lat_min_vs),dimsizes(MODEL)+dim_OBS,nyear/),"double")
   v_wgt_sel2=new((/dimsizes(lat_min),dimsizes(MODEL)+dim_OBS,nyear/),"double")

   do k=0,dimsizes(lat_min)-1
      do imod=0,dimsizes(MODEL)-1
         fili1=models@dir(imod)+"/"+project_name+"_"+models@case_name(imod)+"_"+models@name(imod)+"_"+models@ensemble(imod)+"_"+field_vs+"_"+var_vs+".nc"
         f1= addfile(fili1,"r")
         A0 = extract_data(imod,data,-1,0,0)

         Ae=f1->$var_vs$
         v_wgt_sel1(k,imod,:)=figure3(Ae,field_vs,lev_sel_vs,lat_min_vs(k),lat_max_vs(k),month_sel_vs(k,:),models@name(imod),max(year_min),min(year_max))
				 v_wgt_sel2(k,imod,:)=figure3(A0,field_number,lev_sel,lat_min(k),lat_max(k),month_sele(k,:),models@name(imod),max(year_min),min(year_max))
         year=v_wgt_sel1&time
         delete(A0)
         delete(Ae)
      end do
   end do
   v_wgt_sel1!0="case"
   v_wgt_sel2!0="case"

   v_wgt_sel1!1="model"
   v_wgt_sel2!1="model"

   v_wgt_sel1!2="year"
   v_wgt_sel2!2="year"

   v_wgt_sel1&model=MODELOBS
   v_wgt_sel2&model=MODELOBS

   v_wgt_sel1&year=year
   v_wgt_sel2&year=year

   if (var_vs.eq."vt100") then
      v_wgt_sel1=abs(v_wgt_sel1)
   end if
   do k=0,dimsizes(lat_min)-1
      if (dim_OBS.ne.0) then
         if (isatt(info,"fig03_obs_file").and.isatt(info,"fig03_obs_file_vs").and.isatt(info,"fig03_obs")) then
            do iob=0,dim_OBS-1
               f_obs = addfile(info@fig03_obs_file(iob),"r")
               v_wgt = f_obs->$variable$      ;Observation by ERA are already processed:mean
               tim=doubletointeger(cd_calendar(v_wgt&time,1)/100)
               v_wgt&time = tim
               delete(tim)
               delete(f_obs)
               f_obs = addfile(info@fig03_obs_file_vs(iob),"r")
               v_wgt_vs = f_obs->$var_vs$    ;Observation by ERA  are already processed
               if (var_vs.eq."vt100") then
                  v_wgt_vs=abs(v_wgt_vs)
               end if
               delete(f_obs)
               tim=doubletointeger(cd_calendar(v_wgt_vs&time,1)/100)
               v_wgt_vs&time = tim
               delete(tim)
               if (lat_min_vs(k).gt.0) then
                  ind_lat=ind(v_wgt&lat.gt.0.)
               else
                  ind_lat=ind(v_wgt&lat.le.0.)
               end if
               ind_time=ind(v_wgt_vs&time.ge.min(year).and.v_wgt_vs&time.le.max(year))
               if (.not.all(ismissing(ind_time))) then
                  if (dimsizes(ind_lat).gt.1) then
                     v_wgt_sel1(k,dimsizes(MODEL)+iob,:)=(/dim_avg(v_wgt_vs(time|ind_time,lat|ind_lat))/)
                     v_wgt_sel2(k,dimsizes(MODEL)+iob,:)=(/dim_avg(v_wgt(time|ind_time,lat|ind_lat))/)
                  else
                     v_wgt_sel1(k,dimsizes(MODEL)+iob,:)=(/v_wgt_vs(time|ind_time,lat|ind_lat)/)
                     v_wgt_sel2(k,dimsizes(MODEL)+iob,:)=(/v_wgt(time|ind_time,lat|ind_lat)/)
                  end if
               else
                  print("fatal:Observation time outside the models range time ")
                  status_exit(1)
               end if
               delete(ind_lat)
               delete(ind_time)
               delete(v_wgt_vs)
               delete(v_wgt)
            end do
         end if
      end if
   end do

;**************************************************************************

   y=new((/dimsizes(lat_min_vs),dimsizes(MODELOBS),dimsizes(year)/),"double")
   do k=0,dimsizes(lat_min)-1
      do i=0,dimsizes(MODELOBS)-1

         rc=regline(v_wgt_sel1(k,i,:),v_wgt_sel2(k,i,:))
         y(k,i,:)=rc*(v_wgt_sel1(k,i,:)-rc@xave)+rc@yave
         delete(rc)
      end do
   end do

;******************plot************************************************************

   wks=gsn_open_wks ("ps",diri_out+"/E06FIG03")

   plot1=new(dimsizes(lat_min),"graphic")
   plot=new(dimsizes(lat_min),"graphic")

;***************************************************************
   indexe=new(dimsizes(MODEL)+dim_OBS,"integer")
   indexe(:) = 0  ;(*)
   indexe(:dimsizes(MODEL)-1)=2
   if (dim_OBS.gt.0) then
      indexe(dimsizes(MODEL))=16     ;(* this line will raise error if dim_OBS < 1 )
   end if
   ;indexe(dimsizes(MODEL)+1:)=0  (* this line will raise error if dim_OBS < 2 )

   do k=0,dimsizes(lat_min)-1

      res=True
      res@gsnDraw=False
      res@gsnFrame=False

      res@trXMaxF=info@fig03_xmax_val(k)
      res@trXMinF=info@fig03_xmin_val(k)
      res@trYMaxF=info@fig03_ymax_val(k)
      res@trYMinF=info@fig03_ymin_val(k)
      ;res@tiMainString=""
      MarkLinesMod=new((/dimsizes(MODEL)+dim_OBS/),"string")
      do j=0,dimsizes(MODEL)+dim_OBS-1
         MarkLinesMod(j)="Markers"
         ; MarkLineMod(dimsizes(MODEL)+i)="Lines"
      end do
      res@xyMarkLineModes=MarkLinesMod
      Markersty=new(dimsizes(MODEL)+dim_OBS,"integer")
      ;res@xyMarkers = 2
      MarkLineColors=new((/dimsizes(MODELOBS)/),"string")
      MarkLineColors(0:dimsizes(MODELOBS)-1)=colors(:)
      ;MarkLineColors(dimsizes(MODEL)-1:dimsizes(MarkLineColors)-1)=colors(:)
      res@xyMarkerColors=MarkLineColors
      res@xyMarkerSizeF=0.015
      res@xyMarkers=indexe
      res@tiYAxisString    = "Mean "+v_wgt_sel2@long_name+" at "+lev_sel+"hPa, "+month_sele(k,0)+" and "+month_sele(k,1)     ; y-axis label
      res@tiXAxisString    = v_wgt_sel1@long_name+", "+month_sel_vs(k,0)+" and "+month_sel_vs(k,1)
      plot(k)=gsn_csm_xy (wks,v_wgt_sel1(k,:,:),v_wgt_sel2(k,:,:),res)
      res1=True
      res1=True
      res1@gsnDraw=False
      res1@gsnFrame=False
      res1@xyLineColors=MarkLineColors
      MarkLineMod2=new((/dimsizes(MODELOBS)/),"string")
      do j=0,dimsizes(MODELOBS)-1
         MarkLineMod2(j)="Lines"
      end do
      res1@xyMarkLineModes=MarkLineMod2
      plot1(k)=gsn_csm_xy(wks,v_wgt_sel1(k,:,:),y(k,:,:),res1)
      overlay(plot(k),plot1(k))
      ;draw(plot(k))
      ;frame(wks)
   end do

   Pres=True
   ;Pres@gsnPanelRowSpec=True
   ;num_p=floattointeger((dimsizes(season))/2)
   add_legenda_mark(wks,MODELOBS,colors,indexe)

   gsn_panel(wks,plot,(/2,2/),Pres)
   v_wgt_sel1!0="case"
   v_wgt_sel2!0="case"

   v_wgt_sel1!1="model"
   v_wgt_sel2!1="model"

   v_wgt_sel1!2="year"
   v_wgt_sel2!2="year"

   v_wgt_sel1&model=MODELOBS
   v_wgt_sel2&model=MODELOBS

   v_wgt_sel1&year=year
   v_wgt_sel2&year=year

   ;res@MarkLinesModes="Lines"
   ;plot1=gsn_csm_xy (wks,v_wgt_sel1,y,res)
   l_min=lat_names(lat_min)
   l_max=lat_names(lat_max)
   l_min_vs=lat_names(lat_min_vs)
   l_max_vs=lat_names(lat_max_vs)

;#################################################################################
   if (OPTIONS_FILE.eq."True") then
      ;do i=0,dimsizes(MODEL)-1
      do k=0,dimsizes(lat_min)-1
         var_name1=variable
         var_name2=var_vs

         filo1= (project_name+var_name1+"_"+l_min(k)+"_"+l_max(k)+".nc")
         filo2= (project_name+var_name2+"_"+l_min_vs(k)+"_"+l_max_vs(k)+".nc")
         f_out1=diri_out+"/E06FIG03/"+"E06FIG03_"+filo1
         f_out2=diri_out+"/E06FIG03/"+"E06FIG03_"+filo2

         ;########################################################################
         ; OUTPUT FILE
         system("/bin/rm -f "+f_out1)
         fout1=addfile(f_out1,"c")
         setfileoption(fout1,"DefineMode",True)

         ;===================================================================
         ; create global attributes of the file
         ;===================================================================
         fAtt               = True            ; assign file attributes
         fAtt@title="Selected Months latitude averaged data (model,year)"
         fAtt@creation_date=systemfunc("date")
         do im=0,dimsizes(MODELOBS)-1
            str=unique_string("model_")
            fAtt@$str$= MODELOBS(im)
         end do

         ;copiare attributi del file sorgente
         fileattdef(fout1,fAtt)            ; copy file attributes

         dimNames = (/"year","model"/)
         dimSizes = (/-1,dimsizes(MODELOBS)/)
         dimUnlim = (/True,False/)
         filedimdef(fout1,dimNames,dimSizes,dimUnlim)
         ;===================================================================

         filevardef(fout1, "year"  ,typeof(year),"year")
         filevardef(fout1, variable  ,typeof(v_wgt_sel1)  ,(/"year","model"/))
         ;===================================================================

         setfileoption(fout1,"DefineMode",False)

         fout1->year    = (/year/)
         fout1->$variable$      = (/v_wgt_sel1(year|:,case|k,model|:)/)

         ; OUTPUT FILE
         system("/bin/rm -f "+f_out2)
         fout2=addfile(f_out2,"c")
         setfileoption(fout2,"DefineMode",True)

         ;===================================================================
         ; create global attributes of the file
         ;===================================================================
         delete(fAtt)
         fAtt               = True            ; assign file attributes

         fAtt@title="Selected Months latitude averaged data (model,year)"
         fAtt@creation_date=systemfunc("date")
         do im=0,dimsizes(MODELOBS)-1
            str=unique_string("model_")
            fAtt@$str$= MODELOBS(im)
         end do
         ;copiare attributi del file sorgente
         fileattdef(fout2,fAtt)            ; copy file attributes
         dimNames = (/"year","model"/)
         dimSizes = (/-1,dimsizes(MODELOBS)/)
         dimUnlim = (/True,False/)
         filedimdef(fout2,dimNames,dimSizes,dimUnlim)

         ;===================================================================

         filevardef(fout2, "year"  ,typeof(year),"year")
         filevardef(fout2, var_name2  ,typeof(v_wgt_sel1)  ,(/"year","model"/))
         ;===================================================================
         setfileoption(fout2,"DefineMode",False)
         fout2->year    = (/year/)
         fout2->$var_name2$      = (/v_wgt_sel2(year|:,case|k,model|:)/)
         delete(fAtt)
      end do
   else
      print("OPTION_FILE is"+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
   end if

   print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   print("!  Figure E06FIG03.ps is produced in "+diri_out)
   print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   print ("=================== FINISH (E06FIG03) ====================")

end
;##################################################################################
