;##################################################################################
;# This NCL script reads in CF compliant netCDF CCMVal model output as specified  #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and produces FIGURE 1 of                                                       #
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
;#  Description: This script is designed to plot profiles of climatological mean  #
;#  temperature biases at polar latitudes. The biases are relative to ERA-40      #
;#  reanalyses. Also plotted is the multi-model mean and the standard deviation   #
;#  of the ERA-40 data about the climatological mean. Variables plotted can also  #
;#  be output to netcdf files.                                                    #
;#  The number of panels in the plot is determined by the number of latitude      #
;#  ranges selected in the info array (info@fig01_lat_min ...)                    #
;##################################################################################
;# INPUT/OUTPUT:                                                                  #
;# 1. Reads in                                                                    #
;#    (a) '${WKDIR}/var_att/ta_att.ncl'   (${WKDIR} specified in main namelist)   #
;#    (b) '${WKDIR}/namelist_EO6'         (specification of model names etc.)     #
;#    (c) '${WKDIR}/diag_att/E06Diag.att' (specifies plot type)                   #
;#    (d) Model input files               (as defined in namelists)               #
;# 2. produces                                                                    #
;#    (a) '${diri_out}/E06FIG01.ps'       (plot as ps-file)                       #
;#    (b) '${diri_out}/E06FIG01/E06FIG01_*.nc (plot variables as CF netCDF file)  #
;##################################################################################
; This script takes input from the ta_att.ncl file
; info attributes required:
; fig01_lat_min         array of latitude (start) values for area averaging
; fig01_lat_max         array of latitude (end) values for area averaging
; fig01_season          array of for averaging. (DJF, MAM, JJA, SON or ANN for
;                       annual mean) Note: this array must have the same
;                       dimensions as fig01_lat_min and fig01_lat_max
; fig01_refModel        name of the model or observation file as in the namelist_E06
;                       used to calculate bias and the standard deviation.
;                       The standard deviation is plotted as shadow line.

;**********************************************************************************
load "./plot_type/CCMVal_FUNCTION/new_function.ncl"
load "./plot_type/CCMVal_FUNCTION/Climatology_CCMVal.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
load "./plot_type/CCMVal_FUNCTION/proc_E06FIG01.ncl"
load "./plot_type/get_field_number.ncl"
load "./plot_type/uniform_coordinates.ncl"

begin

; fields supported T3M T3D T2Mz T2Dz
   print("Begin ./plot_type/E06FIG01.ncl")
   diri_out = getenv("ESMValTool_plot_dir")
   project_name = getenv("ESMValTool_project")
   OPTIONS_FILE = getenv("ESMValTool_write_plot_vars")

   dim_MOD = dimsizes(models@name)-1

   year_min=stringtointeger(models@start_date)
   year_max=stringtointeger(models@end_date)


   ;******************************************************************

   if (isatt(info,"fig01_refModel")) then
      idx_ref = ind(models@name .eq. info@fig01_refModel)
      if (all(ismissing(idx_ref))) then
         print("fatal:refModel indicated in /var_att/ta_att.ncl isn't included in model list (namelist).If refModel is not included bias cannot be calculated.")
         print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
         print("fatal:!No files or plots are been produced!")
         print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
         status_exit(1)
      end if
   else
      print("fatal:plot requires variable attribute 'refModel'.")
      print("fatal:In /var_att/ta_att.ncl info_01@refModel must point to a Model or Observation Name as indicated in the namelist_E06")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if
   if (isatt(info,"fig01_climObs")) then
      dim_clim = dimsizes(info@fig01_climObs)
      if (dim_clim .ne. dimsizes(info@fig01_climObs_file)) then
          print("fatal:dimensions of 'fig01_climObs' and 'fig01_climObs_file' do not match.")
          print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
          print("fatal:!No files or plots are been produced!")
          print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
          status_exit(1)
      end if
   else
      dim_clim = 0
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
   if(isatt(info,"fig01_lat_max") .and. isatt(info,"fig01_lat_min")) then
      lat_max = info@fig01_lat_max
      lat_min = info@fig01_lat_min
      lat_name_max = lat_names(lat_max)
      lat_name_min = lat_names(lat_min)
   else
      print("fatal:plot requires variable attributes 'fig01_lat_max' and 'fig01_lat_min'")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if

   if(isvar("season")) then
      delete(season)
   end if
   if(isatt(info,"fig01_season")) then
      season = info@fig01_season
   else
      print("fatal:plot requires variable attribute 'fig01_season'")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if

   if (((dimsizes("lat_max")).ne.(dimsizes("lat_min"))) .or. ((dimsizes("lat_max")).ne.(dimsizes("season")))) then
      print("fatal:Error!!!!The number of season must to be the same of the number of lat_max & lat_min!!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if

   ; optional input parameters
   if(isatt(info,"fig01_xmax")) then
      x_Max = info@fig01_xmax
   else
      x_Max = 20.
   end if
   if(isatt(info,"fig01_xmin")) then
      x_Min = info@fig01_xmin
   else
      x_Min = -20.
   end if

   ;*************************names and colors**************************
   colors=new(dim_MOD+dim_clim,string)
   dashs=new(dim_MOD+dim_clim,integer)
   MarkLin=new(dim_MOD+dim_clim,string)
   modelli=new(dim_MOD+dim_clim,string)

   jmod=0

   do imod=0,dim_MOD
      if (imod.ne.idx_ref) then
         if (project_name.eq."CCMVal2") then
            colors(jmod)= CCMVal2_ColorTable(models@name(imod))
            dashs(jmod) = CCMVal2_DashTable(models@name(imod))
         else
            colors(jmod)= CCMVal1_ColorTable(models@name(imod))
            dashs(jmod) = CCMVal1_DashTable(models@name(imod))
         end if

         modelli(jmod) = models@name(imod)
         MarkLin(jmod) = "Lines"
         jmod=jmod+1
      end if
   end do
   do imod=0,dim_clim-1
      if (project_name.eq."CCMVal2") then
         colors(jmod) = CCMVal2_ColorTable(info@fig01_climObs(imod))
         dashs(jmod) = CCMVal2_DashTable(info@fig01_climObs(imod))
      else
         colors(jmod) = CCMVal1_ColorTable(info@fig01_climObs(imod))
         dashs(jmod) = CCMVal1_DashTable(info@fig01_climObs(imod))
      end if

      modelli(jmod) = info@fig01_climObs(imod)
      MarkLin(jmod) = "Markers"
      jmod=jmod+1
   end do


   ;******************************************************************
   ; read obs file and extract data

   E40 = extract_data(idx_ref,data,-1,0,0)
   field = get_field_number(E40)
   pout1 = E40&plev

   ;*****************************************************************
   v_wgt_E40 = new((/2,dimsizes(season),dimsizes(pout1)/),"double")
   v_wgt_diff = new((/dim_MOD+dim_clim,dimsizes(season),dimsizes(pout1)/),"double")
   v_wgt_val = new((/2,dim_MOD+dim_clim,dimsizes(season),dimsizes(pout1)/),"double")

   v_wgt_diff!0 = "model"
   v_wgt_diff!1 = "season"
   v_wgt_diff!2 = "plev"
   v_wgt_diff&model = modelli
   v_wgt_diff&season = season
  ; v_wgt_diff&plev = pout

   v_wgt_val!0 = "case"
   v_wgt_val!1 = "model"
   v_wgt_val!2 = "season"
   v_wgt_val!3 = "plev"
   v_wgt_val&model = modelli
   v_wgt_val&season = season
  ; v_wgt_val&plev = pout

   do k=0,dimsizes(season)-1
      print("*** processing season: "+season(k)+" (set "+(k+1)+" of " +dimsizes(season)+") ***")
      print("processing "+models@name(idx_ref)+" (refModel)")
      v_wgt_E40(:,k,:) = figure1(E40,field,season(k),models@name(idx_ref),lat_max(k),lat_min(k),year_min(idx_ref),year_max(idx_ref))
      pout = v_wgt_E40&plev

      jmod=0
      do imod=0,dim_MOD
         if (imod.ne.idx_ref) then
            A0 = extract_data(imod,data,-1,0,0)
            print("processing "+models@name(imod)+" (MODEL: "+(jmod+1)+" of " +dim_MOD+")")
            v_wgt = figure1(A0,field_number,season(k),models@name(imod),lat_max(k),lat_min(k),year_min(imod),year_max(imod))
            pin = v_wgt&plev
            v_wgt_c = int2p(pin,v_wgt,pout,2)
            v_wgt_diff(jmod,k,:) = v_wgt_c(0,0,:)-v_wgt_E40(0,k,:)
            v_wgt_val(:,jmod,k,:)=(/v_wgt_c(:,0,:)/)
            jmod=jmod+1
            delete(v_wgt)
            delete(v_wgt_c)
            delete(A0)
            delete(pin)
         end if
      end do
      if (dim_clim .gt. 0) then
         obs_filename=info@fig01_climObs_file
         do imod=0,dim_clim-1
            f_obs = addfile(obs_filename(imod),"r")
            A0 = f_obs->$variable$      ;Observation by NCEP are already processed:mean
            print("processing "+info@fig01_climObs(imod)+" (OBS: "+(imod+1)+" of " +dim_clim+")")
            v_wgt = figure1(A0,"C2Mz",season(k),info@fig01_climObs(imod),lat_max(k),lat_min(k),0.,0.)
            pin = v_wgt&plev
            v_wgt_c = int2p(pin,v_wgt,pout,2)
            v_wgt_diff(jmod,k,:)= v_wgt_c(0,0,:)-v_wgt_E40(0,k,:)
            v_wgt_val(:,jmod,k,:)=(/v_wgt_c(:,0,:)/)
            jmod=jmod+1
            delete(A0)
            delete(f_obs)
            delete(v_wgt)
            delete(v_wgt_c)
            delete(pin)
         end do
      end if
      ; v_wgt_diff&plev = pout


      ;**************  OUTPUT FILES  *********************************
      if (OPTIONS_FILE.eq."True") then
         do imod=0,dim_MOD+dim_clim-1
            jmod=ind(models@name .eq. modelli(imod))
            if (all(ismissing(jmod))) then
               fili = (project_name+"_"+info@fig01_climObs(imod-dim_MOD)+"_"+field_number+"_"+variable+".nc")
               f_out = diri_out+"/E06FIG01/E06FIG01_"+season(k)+"_"+lat_name_max(k)+lat_name_min(k)+"_"+fili
            else
               fili = (project_name+"_"+models@case_name(jmod)+"_"+models@name(jmod)+"_"+models@ensemble(jmod)+"_"+field_number+"_"+variable+".nc")
               f_out = diri_out+"/E06FIG01/E06FIG01_"+season(k)+"_"+lat_name_max(k)+lat_name_min(k)+"_"+fili
            end if
            system("/bin/rm -f "+f_out)
            fout = addfile(f_out,"c")
            setfileoption(fout,"DefineMode",True)

            ;=========================================================
            ; create global attributes of the file
            ;=========================================================
            fAtt               = True            ; assign file attributes
            fAtt@title         = "Seasonal-mean data (season,plev)"
            ; fAtt@source_file = f
            ;fAtt@creation_date = systemfunc("date")
            fAtt@season        = season(k)
            fAtt@model         = modelli(imod) ; models@name(jmod)
            fAtt@comment       = "This file contains the output of FIG01.The standard deviation and the values ofReference model/Observation are also included"
            var_name1=0.
            var_name2=0.
            var_name3=0.
            var_name1@long_name= "bias vs "+info@fig01_refModel
            var_name2@long_name= " nv=0, climatology of "+info@fig01_refModel +" nv=1, standard deviation of the climatological mean"
            var_name3@long_name= " nv=0, climatology"  +" nv=1, standard deviation of the climatological mean"

           ; copy attributes from the source file
            fileattdef(fout,fAtt)            ; copy file attributes
            dimNames = (/"nv","plev"/)
            dimSizes = (/2,dimsizes(pout)/)
            dimUnlim = (/False,False/)
            filedimdef(fout,dimNames,dimSizes,dimUnlim)

            ;=========================================================
            var_name = variable+"_diff"
            var_namev = variable

            filevardef(fout,"plev",typeof(pout),"plev")
            filevardef(fout,var_name,typeof(v_wgt_diff),(/"plev"/))
            filevardef(fout,info@fig01_refModel,typeof(v_wgt_E40),(/"nv","plev"/))
            filevardef(fout,var_namev,typeof(v_wgt_val),(/"nv","plev"/))
            filevarattdef(fout,var_name,var_name1)
            filevarattdef(fout,info@fig01_refModel,var_name2)
            filevarattdef(fout,var_namev,var_name3)
            setfileoption(fout,"DefineMode",False)
            fout->plev       = (/pout/)
            fout->$var_name$ = (/v_wgt_diff(model|imod,season|k,plev|:)/)
            fout->$info@fig01_refModel$ = (/v_wgt_E40(case|:,season|k,plev|:)/)
            fout->$var_namev$ = (/v_wgt_val(case|:,model|imod,season|k,plev|:)/)
         end do
         print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
         print("!")
         print("!  Files are been produced in "+diri_out+"/E06FIG01/  !!")
         print("!")
         print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      else
         print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
      end if
   end do

   ;******************************************************************

   ;****************  plot  ******************************************
   plot_diff = new((/dimsizes(season)/),"graphic")
   plot_ERA = new((/dimsizes(season)/),"graphic")

   wks = gsn_open_wks ("ps",diri_out+"/E06FIG01")

   do i=0,dimsizes(season)-1

      ;***************************************************************
      xp = new((/2*dimsizes(pout),dimsizes(season)/),typeof(v_wgt_E40))
      yp = new((/2*dimsizes(pout),dimsizes(season)/),typeof(v_wgt_E40))

      do j=0,dimsizes(pout)-1
       ;  xp(j,i) = -((v_wgt_E40(1,i,j)))*0.468
         xp(j,i) = -((v_wgt_E40(1,i,j)))
         yp(j,i) = v_wgt_E40&plev(j)
         yp(dimsizes(pout)+j,i) = v_wgt_E40&plev(dimsizes(pout)-1-j)
       ;  xp(dimsizes(pout)+j,i) = (v_wgt_E40(1,i,(dimsizes(pout)-1-j)))*0.468
         xp(dimsizes(pout)+j,i) = (v_wgt_E40(1,i,(dimsizes(pout)-1-j)))
      end do

      res                   = True
      res@gsnDraw           = False
      res@gsnFrame          = False
      res@tfPolyDrawOrder   = "Predraw"
      res@trXMaxF           = x_Max
      res@trXMinF           = x_Min
      res@trYMaxF           = 500
      res@trYMinF           = 1
      res@tmYLMode          = "Explicit"
      res@tmYLValues        = (/1,3,10,30,100,300,1000/)
      res@tmYLLabels        = (/1,3,10,30,100,300,1000/)
      res@tmLabelAutoStride = True

      res@vpHeightF         = 0.2                    ; change aspect ratio of plot
      res@vpWidthF          = 1./dimsizes(season)

      res@tiYAxisString     = "Pressure (hPa)"       ; y-axis label
      res@tiXAxisString     = "Model "+E40@long_name+" bias "+"("+E40@units+")"     ; y-axis label
      res@tiMainString      = season(i)+" "+lat_name_max(i)+"-"+lat_name_min(i)

      res@trYAxisType       = "LogAxis"
      res@trYReverse        = True
      ;res@xyMonoMarker      = True

      ;res@xyMarkerSizeF=0.1
      res@xyLineColors      = colors                 ; line color
      res@xyDashPatterns    = dashs
      res@tmXMajorGrid      = True                   ; implement x grid
      res@tmXMajorGridThicknessF = 1.0               ; 2.0 is default
      res@tmXMajorGridLineDashPattern = 2            ; select short dash lines
      ;res@xyMarkerSizeF=0.001
      ;res@xyMonoMarkLineMode = True
      res@xyMarkLineModes    = MarkLin
      ;res@xyMonoMarker       =True
      res@xyMarkers          = dashs
      plot_diff(i) = gsn_csm_xy (wks,v_wgt_diff(model|:,season|i,plev|:),pout,res)

      gsres = True
      gsres@gsFillColor = "Snow3"
      plot_ERA(i) = gsn_add_polygon(wks,plot_diff(i),xp(:,i),yp(:,i),gsres)
   end do

   delete(E40)

   Pres                = True
   Pres@gsnFrame       = False
   Pres@gsnPanelBottom = 0.33
   m_d = mod(dimsizes(season),2)
   if (m_d.eq.0) then
      n_p = dimsizes(season)/2
   else
      n_p = dimsizes(season)/2+1
   end if
   gsn_panel(wks,plot_diff,(/n_p,2/),Pres)

   add_legenda(wks,modelli,colors,dashs,"bottom")
   print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   print("!")
   print("!  Figure E06FIG01.ps is produced in "+diri_out)
   print("!")
   print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   print ("=================== FINISH (E06FIG01) ====================")

   delete(year_min)
   delete(year_max)
   delete(lat_min)
   delete(lat_max)
   delete(lat_name_max)
   delete(lat_name_min)
   delete(season)
   delete(v_wgt_E40)
   delete(v_wgt_diff)
   delete(imod)
   delete(colors)
   delete(dashs)
   delete(m_d)
   delete(plot_diff)
   delete(plot_ERA)
   delete(xp)
   delete(yp)
   delete(Pres)
   delete(res)
   delete(k)
   delete(j)
   delete(i)
   delete(pout)
   delete(n_p)
   delete(dim_MOD)
   delete(dim_clim)
   delete(jmod)
   delete(idx_ref)
   delete(x_Min)
   delete(gsres)
   delete(var_name1)
   delete(var_name2)
   delete(var_name3)
   delete(fAtt)
   delete(modelli)
   delete(f_out)
   delete(fili)

end
;##################################################################################
