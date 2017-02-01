;##################################################################################
;# This NCL script reads in CF compliant netCDF CCMVal model output as specified  #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and produces FIGURE 7 of                                                       #
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
; Description: This script plots the seasonal cycle of a particular variable      #
; (e.g., temperature or H2O for figure 7 of E06) on a chosen pressure level (area #
; averaged). Variables plotted can also be output to netcdf files.                #
; The number of panels in the plot is determined by the number of latitude        #
; ranges selected in the info array (info@fig07_lat_min ...)                      #
;;#################################################################################
;# INPUT/OUTPUT:                                                                  #
;# 1. Reads in                                                                    #
;#    (a) '${WKDIR}/var_att/ta_att.ncl'   (${WKDIR} specified in main namelist)   #
;#    (b) '${WKDIR}/namelist_EO6'         (specification of model names etc.)     #
;#    (c) '${WKDIR}/diag_att/E06Diag.att' (specifies plot type)                   #
;#    (d) Model input files               (as defined in namelists)               #
;# 2. produces                                                                    #
;#    (a) '${diri_out}/E06FIG07.ps'       (plot as ps-file)                       #
;#    (b) '${diri_out}/E06FIG07/'         (plot variables as CF netCDF file)      #
;##################################################################################
; This script takes input from the ta_att.ncl or H2O_att.ncl files (other variables
; can also be chosen)
; info attributes required:
; fig07_lat_min         array of latitude (start) values for area averaging
; fig07_lat_max         array of latitude (end) values for area averaging
; fig07_lev_sel         array of selected pressure level (hPa). Note: this array
;                       must have the same dimensions as fig07_lat_min and
;                       fig07_lat_max

; fig07_refModel        name of the model or observation file as in the namelist_E06
;                       used to calculate the standard deviation.
;                       The standard deviation is plotted as shadow line.

;**********************************************************************************
load "./plot_type/CCMVal_FUNCTION/new_function.ncl"
load "./plot_type/CCMVal_FUNCTION/Climatology_CCMVal.ncl"
load "./plot_type/CCMVal_FUNCTION/proc_E06FIG07.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
load "./plot_type/get_field_number.ncl"

begin

   ;fields supported T3M T3D T2Mz T2Dz T2Ms T2Ds
   print("Begin ./plot_type/E06FIG07.ncl")
   diri_out = getenv("ESMValTool_plot_dir")
   project_name = getenv("ESMValTool_project")
   OPTIONS_FILE = getenv("ESMValTool_write_plot_vars")

   ;******************************************************************

   ;year_min and year_max = range of years to use for climatology
   year_min=stringtointeger(models@start_date)
   year_max=stringtointeger(models@end_date)

   if(isvar("idx_ref")) then
      delete(idx_ref)
   end if
   if (isatt(info,"fig07_refModel")) then
      OBS=info@fig07_refModel

      idx_ref = ind(models@name .eq. info@fig07_refModel)
      if (all(ismissing(idx_ref))) then
         print("refModel isn't included in model list (namelist)")
      end if
   else
      print("no 'refModel'")
      idx_ref=new(1,"integer")
   end if
   if(isatt(info,"fig07_obs_name") .and. isatt(info,"fig07_obs_file") \
                                   .and. isatt(info,"fig07_obs_stdv_file")) then
      OBS = info@fig07_obs_name
      dim_OBS=dimsizes(OBS)
      obs_filename = info@fig07_obs_file
      obs_stdv_file = info@fig07_obs_stdv_file
      idx_ref=new(1,"integer")
   else
      print("no 'OBS'")
   end if


   ;check the input parameters & translate them
   ;lat_min, lat_max = latitude ranges for spatial averaging
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
   if(isatt(info,"fig07_lat_min") .and. isatt(info,"fig07_lat_max")) then
      lat_min = info@fig07_lat_min
      lat_max = info@fig07_lat_max
      lat_name_max = lat_names(lat_max)
      lat_name_min = lat_names(lat_min)
   else
      print("fatal:plot requires variable attributes 'fig07_lat_min' and 'fig07_lat_max'")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("fatal:!No files or plots are been produced!")
      print("fatal:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      status_exit(1)
   end if

   ; lev_sel = pressure level selected
   if(isvar("lev_sel")) then
      delete(lev_sel)
   end if
   if(isatt(info,"fig07_lev_sel")) then
      lev_sel = info@fig07_lev_sel

      if (lev_sel.eq."cold_point") then
         lev_sel1=-999
      else
         lev_sel1=lev_sel
      end if
   else
      lev_sel=""
      lev_sel1=-999
   end if
   lev_sel1@_FillValue=-999
   ;******************************************************************
   if(isatt(info,"fig07_obs_name") .and. isatt(info,"fig07_obs_file") \
                                   .and. isatt(info,"fig07_obs_stdv_file")) then
      Obs_mod=new(dimsizes(models@name)+1,typeof(models@name))
      Obs_mod(dimsizes(Obs_mod)-1) = info@fig07_obs_name
      Obs_mod(:dimsizes(Obs_mod)-2) = models@name
   else
      Obs_mod = models@name
   end if
   dim_MOD=dimsizes(Obs_mod)

   colors= new(dim_MOD, "string")
   dashs = new(dim_MOD, "integer")
   modelli=new(dim_MOD, "string")

   v_wgt = new((/dimsizes(lat_max), dim_MOD,12/), "double")
   v_wgt_sdt = new((/dimsizes(lat_max), dim_MOD,12/), "double")

   ;*****************************************************************
   ; loop over latitude ranges
   do k = 0,dimsizes(lat_max)-1
      jmod=0
      ;OBS (ERA40)
      if (isatt(info,"fig07_obs_name") .and. isatt(info,"fig07_obs_file") \
                                   .and. isatt(info,"fig07_obs_stdv_file")) then

         f_obs = addfile(obs_filename,"r")
         v = f_obs->$variable$      ;Observation by HALOE are already processed:mean
         delete(f_obs)
         f_obs = addfile(obs_stdv_file,"r")
         vName = variable + "@STDEV"
         v_STD = f_obs->$vName$     ;Observation by HALOE  are already processed:standard deviation
         delete(f_obs)
         v!0 = "month"
         v!1 = "plev"
         v!2 = "lat"
         v_STD!0 = "month"
         v_STD!1 = "plev"
         v_STD!2 = "lat"
         if ((v@units.eq."m^3/m^3").and.(variable.eq."HCl")) then
            v=v*10^9
            v_STD=v_STD*10^9
            v@units="ppbv"
            v_STD@units="ppbv"
         end if
         if (v@units.eq."m^3/m^3") then
            v=v*10^6
            v_STD=v_STD*10^6
            v@units="ppm"
            v_STD@units="ppm"
         end if
         if ((lat_max.eq.0.).and.(lat_min.eq.0.)) then
            ind_lat = ind((v&lat.ge.(-10.)).and.(v&lat.le.(10.)))
         else
            ind_lat = ind((v&lat.ge.(lat_min)).and.(v&lat.le.lat_max))
         end if
         if (lev_sel.eq."cold_point") then

            v_m=dim_avg_Wrap(v(month|:,lat|ind_lat))
            v_STD_m=dim_avg_Wrap(v_STD(month|:,lat|ind_lat))
            v_wgt(k,dimsizes(Obs_mod)-1,:)=(/v_m(month|:)/)
            v_wgt_sdt(k,dimsizes(Obs_mod)-1,:)=(/v_STD_m(month|:)/)
            delete(v_m)
            delete(v_STD_m)
            if (project_name.eq."CCMVal2") then
               colors(dimsizes(Obs_mod)-1) = CCMVal2_ColorTable(Obs_mod(dimsizes(Obs_mod)-1))
               dashs(dimsizes(Obs_mod)-1) = CCMVal2_DashTable(Obs_mod(dimsizes(Obs_mod)-1))
            else
               colors(dimsizes(Obs_mod)-1) = CCMVal1_ColorTable(Obs_mod(dimsizes(Obs_mod)-1))
               dashs(dimsizes(Obs_mod)-1) = CCMVal1_DashTable(Obs_mod(dimsizes(Obs_mod)-1))
            end if
            modelli(dimsizes(Obs_mod)-1)=Obs_mod(dimsizes(Obs_mod)-1)
            delete(ind_lat)
            delete(lev_sel)
            lev_sel1=-999

         else
            v_m=dim_avg_Wrap(v(month|:,{plev|lev_sel},lat|ind_lat))
            v_STD_m=dim_avg_Wrap(v_STD(month|:,{plev|lev_sel},lat|ind_lat))
            v_wgt(k,dimsizes(Obs_mod)-1,:)=(/v_m(month|:)/)
            v_wgt_sdt(k,dimsizes(Obs_mod)-1,:)=(/v_STD_m(month|:)/)
            delete(v_m)
            delete(v_STD_m)
            if (project_name.eq."CCMVal2") then
               colors(dimsizes(Obs_mod)-1) = CCMVal2_ColorTable(Obs_mod(dimsizes(Obs_mod)-1))
               dashs(dimsizes(Obs_mod)-1) = CCMVal2_DashTable(Obs_mod(dimsizes(Obs_mod)-1))
            else
               colors(dimsizes(Obs_mod)-1) = CCMVal1_ColorTable(Obs_mod(dimsizes(Obs_mod)-1))
               dashs(dimsizes(Obs_mod)-1) = CCMVal1_DashTable(Obs_mod(dimsizes(Obs_mod)-1))
            end if
            modelli(dimsizes(Obs_mod)-1)=Obs_mod(dimsizes(Obs_mod)-1)
            delete(ind_lat)
         end if
      else
         if (isatt(info,"fig07_obs_name") .and. isatt(info,"fig07_obs_file") ) then

            f_obs = addfile(obs_filename,"r")
            v = f_obs->$variable$      ;Observation by HALOE are already processed:mean

            v!0 = "month"
            v!1 = "plev"
            v!2 = "lat"

            if ((lat_max.eq.0.).and.(lat_min.eq.0.)) then
               ind_lat = ind((v&lat.ge.(-5.)).and.(v&lat.le.(5.)))
            else
               ind_lat = ind((v&lat.ge.(lat_min)).and.(v&lat.le.lat_max))
            end if
            v_m=dim_avg_Wrap(v(month|:,{plev|lev_sel},lat|ind_lat))
            v_wgt(k,dimsizes(Obs_mod)-1,:)=(/v_m(month|:)/)
            delete(v_m)
            if (project_name.eq."CCMVal2") then
               colors(dimsizes(Obs_mod)-1) = CCMVal2_ColorTable(Obs_mod(dimsizes(Obs_mod)-1))
               dashs(dimsizes(Obs_mod)-1) = CCMVal2_DashTable(Obs_mod(dimsizes(Obs_mod)-1))
            else
               colors(dimsizes(Obs_mod)-1) = CCMVal1_ColorTable(Obs_mod(dimsizes(Obs_mod)-1))
               dashs(dimsizes(Obs_mod)-1) = CCMVal1_DashTable(Obs_mod(dimsizes(Obs_mod)-1))
            end if
            modelli(dimsizes(Obs_mod)-1)=Obs_mod(dimsizes(Obs_mod)-1)
            delete(ind_lat)
         end if
         if (.not.all(ismissing(idx_ref))) then
            E40 = extract_data(idx_ref,data,-1,0,0)
            field = get_field_number(E40)
            print("Model: "+models@name(idx_ref))
            v_wgt(k,dimsizes(Obs_mod)-1,:) = figure7(E40,field,lat_max(k),lat_min(k),lev_sel(k),year_min(idx_ref),year_max(idx_ref),"mean")
            v_wgt_sdt(k,dimsizes(Obs_mod)-1,:) = figure7(E40,field,lat_max(k),lat_min(k),lev_sel(k),year_min(idx_ref),year_max(idx_ref),"stdv")
            if (project_name.eq."CCMVal2") then
               colors(dimsizes(Obs_mod)-1) = CCMVal2_ColorTable(Obs_mod(idx_ref))
               dashs(dimsizes(Obs_mod)-1) = CCMVal2_DashTable(Obs_mod(idx_ref))
            else
               colors(dimsizes(Obs_mod)-1) = CCMVal1_ColorTable(Obs_mod(idx_ref))
               dashs(dimsizes(Obs_mod)-1) = CCMVal1_DashTable(Obs_mod(idx_ref))
            end if
            modelli(dimsizes(Obs_mod)-1)=Obs_mod(idx_ref)
            delete(E40)
            jmod=0
         end if

      end if
      ; MODELS
      do imod = 0,dimsizes(models@name)-1
         if (.not.all(ismissing(idx_ref))) then
            if (imod.ne.idx_ref) then
               A0 = extract_data(imod,data,-1,0,0)
               field=get_field_number(A0)

               print("Model: "+models@name(imod))
               v_wgt(k,jmod,:) = figure7(A0,field,lat_max(k),lat_min(k),lev_sel1(k),year_min(imod),year_max(imod),"mean")
               v_wgt_sdt(k,jmod,:) = figure7(A0,field,lat_max(k),lat_min(k),lev_sel1(k),year_min(imod),year_max(imod),"stdv")
               if (project_name.eq."CCMVal2") then
                  colors(jmod) = CCMVal2_ColorTable(Obs_mod(imod))
                  dashs(jmod) = CCMVal2_DashTable(Obs_mod(imod))
               else
                  colors(jmod) = CCMVal1_ColorTable(Obs_mod(imod))
                  dashs(jmod) = CCMVal1_DashTable(Obs_mod(imod))
               end if
               modelli(jmod)=Obs_mod(imod)
               jmod=jmod+1
               delete(A0)

            end if
         else
            A0 = extract_data(imod,data,-1,0,0)
            field=get_field_number(A0)

            print("Model: "+models@name(imod))
            v_wgt(k,imod,:) = figure7(A0,field,lat_max(k),lat_min(k),lev_sel1(k),year_min(imod),year_max(imod),"mean")
            v_wgt_sdt(k,imod,:) = figure7(A0,field,lat_max(k),lat_min(k),lev_sel1(k),year_min(imod),year_max(imod),"stdv")
            if (project_name.eq."CCMVal2") then
               colors(imod) = CCMVal2_ColorTable(Obs_mod(imod))
               dashs(imod) = CCMVal2_DashTable(Obs_mod(imod))
            else
               colors(imod) = CCMVal1_ColorTable(Obs_mod(imod))
               dashs(imod) = CCMVal1_DashTable(Obs_mod(imod))
            end if
            modelli(imod)=Obs_mod(imod)
            delete(A0)
         end if
      end do
   end do

   sel = new(dimsizes(lat_max), typeof(lat_max))
   nv = new(1,"integer")
   nv = 3
   sel_point = new((/dimsizes(lat_max),nv/),typeof(lat_max))
   sel_point(:,0) = lat_min
   sel_point(:,1) = lat_max
   sel_point(:,2) = lev_sel1
   sel_point!0 = "sel"
   sel_point!1 = "nv"

   v_wgt!0 = "sel"
   v_wgt!1 = "model"
   v_wgt!2 = "month"

   v_wgt_sdt!0 = "sel"
   v_wgt_sdt!1 = "model"
   v_wgt_sdt!2 = "month"

   month = (/1,2,3,4,5,6,7,8,9,10,11,12/)
   ;month = (/"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"/)

   v_wgt&model = modelli
   v_wgt&month = month

   ;*****************************************************************
   wks = gsn_open_wks ("ps",diri_out+"/E06FIG07_"+variable)


   ;*****************************************************************
   plot = new((/dimsizes(lat_max)/),"graphic")
   plot_b = new((/dimsizes(lat_max)/),"graphic")
   plot_ERA = new((/dimsizes(lat_max)/),"graphic")

   ;*****************************************************************
   ; shading for observation data (standard deviation)
   if (.not.all(ismissing(idx_ref))).or.(isatt(info,"fig07_obs_stdv_file"))then
      do il = 0,dimsizes(lev_sel)-1
         xp = new((/2*dimsizes(month),dimsizes(lev_sel)/),typeof(v_wgt))
         yp = new((/2*dimsizes(month),dimsizes(lev_sel)/),typeof(v_wgt))

         do j = 0,dimsizes(month)-1
            xp(j,il) = ((v_wgt(il,dimsizes(Obs_mod)-1,j)+v_wgt_sdt(il,dimsizes(Obs_mod)-1,j)))
            yp(j,il) = v_wgt&month(j)
            yp(dimsizes(month)+j,il) = v_wgt&month(dimsizes(month)-1-j)
            xp(dimsizes(month)+j,il) = (v_wgt(il,dimsizes(Obs_mod)-1,(dimsizes(month)-1-j))-v_wgt_sdt(il,dimsizes(Obs_mod)-1,(dimsizes(month)-1-j)))
         end do
      end do
   end if

   ; plots
   do il = 0,dimsizes(lev_sel)-1
      res                 = True
      res@gsnDraw         = False
      res@gsnFrame        = False

      res@tfPolyDrawOrder = "Predraw"

      res@trXMaxF         = 13
      res@trXMinF         = 0

      res@tmXBMode        = "Explicit"
      res@tmXBValues      = (/1,2,3,4,5,6,7,8,9,10,11,12/)
      res@tmXBLabels      = (/"J","F","M","A","M","J","J","A","S","O","N","D"/)

      res@vpHeightF       = 0.4                    ; change aspect ratio of plot
      res@vpWidthF        = 0.8

      res@tiYAxisString   = v_wgt@long_name+"("+info@units+")"    ; y-axis label
      res@tiXAxisString   = "Month"     ; y-axis label
      if (isatt(info,"long_name")) then
         varLongName = info@long_name
      else
         varLongName = v_wgt@long_name
      end if
      if (.not.ismissing(lev_sel1(il))) then
         if (lat_name_min(il).ne.lat_name_max(il)) then
            res@tiMainString = lev_sel(il)+" hPa "+varLongName+" at "+lat_name_min(il)+" "+lat_name_max(il)
         else
            res@tiMainString = lev_sel(il)+" hPa "+varLongName+" at "+lat_name_min(il)
         end if
      else
         if (lat_name_min(il).ne.lat_name_max(il)) then
            res@tiMainString = varLongName+" at "+lat_name_min(il)+" "+lat_name_max(il)
         else
            res@tiMainString = varLongName+" at "+lat_name_min(il)
         end if
      end if

      res@xyLineColors    = colors                    ; line color
      res@xyDashPatterns  = dashs

      plot_b(il) = gsn_csm_xy(wks,month,v_wgt(il,:,:),res)

      gsres = True
      gsres@gsFillColor = "Snow3"

      if (.not.all(ismissing(idx_ref))).or.(isatt(info,"fig07_obs_stdv_file"))then
         plot_ERA(il) = gsn_add_polygon(wks,plot_b(il),yp(:,il),xp(:,il),gsres)
      end if

   end do

   Pres = True
   Pres@gsnFrame=False
   Pres@gsnPanelBottom=0.33

   add_legenda(wks,modelli,colors,dashs,"bottom")
   gsn_panel(wks,plot_b,(/dimsizes(lev_sel),1/),Pres)
   print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   print("!")
   print("!  Figure E06FIG07_"+variable+".ps is produced in "+diri_out)
   print("!")
   print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

   ;*****************************************************************
   model1 = new((/dimsizes(Obs_mod),14/),"character")
   do j = 0,dimsizes(Obs_mod)-1
      modil = stringtochar(modelli(j))
      model1(j,:dimsizes(modil)-1) = modil
      delete(modil)
   end do

   ;*****************************************************************
   if (OPTIONS_FILE.eq."True") then

      ; OUTPUT FILE
      do k = 0,dimsizes(lat_max)-1
         fol = "/E06FIG07/E06FIG07_"+lev_sel(k)+"_"+lat_name_max(k)+lat_name_min(k)+project_name+"_"+variable+".nc"
         f_out = diri_out+fol
         system("/bin/rm -f "+f_out)
         fout = addfile(f_out,"c")
         setfileoption(fout,"DefineMode",True)

         ;*****************************************************************
         ; create global attributes of the file
         ;*****************************************************************
         fAtt              = True            ; assign file attributes
         fAtt@title        = "Annual cycle (sel,model,month)"
         fAtt@creation_date= systemfunc("date")
         do im=0,dimsizes(Obs_mod)-1
            str=unique_string("model_")
            fAtt@$str$= Obs_mod(im)
         end do
         ;fAtt@model=obs_model
         ; copy attributes from the source file
         fileattdef(fout,fAtt)            ; copy file attributes
         dimNames = (/"model","month","nv","n_ch"/)
         dimSizes = (/dimsizes(Obs_mod),dimsizes(month),3,14/)
         dimUnlim = (/False,False,False,False/)
         filedimdef(fout,dimNames,dimSizes,dimUnlim)

         ;*****************************************************************
         var_name = variable+"_annual_cycle"
         sdt_var_name = variable+"_annual_cycle_sdt"

         filevardef(fout,"sel_pnt",typeof(sel_point),(/"nv"/))
         filevardef(fout,"model","character", (/"model","n_ch"/))
         filevardef(fout,sdt_var_name,typeof(v_wgt_sdt),(/"model","month"/))
         filevardef(fout,var_name,typeof(v_wgt),(/"model","month"/))

         ;*****************************************************************
         filevarattdef(fout,sdt_var_name,v_wgt_sdt)                    ; copy var attributes
         filevarattdef(fout,var_name,v_wgt)                    ; copy var attributes
         filevarattdef(fout,"sel_pnt",sel_point)                     ; copy time attributes
         filevarattdef(fout,"model",model1)                       ; copy lat attributes

         setfileoption(fout,"DefineMode",False)

         fout->sel_pnt = (/sel_point(k,:)/)
         fout->model = (/model1/)

         fout->$var_name$     = (/v_wgt(k,:,:)/)
         fout->$sdt_var_name$ = (/v_wgt_sdt(k,:,:)/)
      end do
      delete(fol)
      delete(f_out)
      delete(sdt_var_name)
      delete(dimUnlim)
      delete(dimNames)
      delete(dimSizes)
      delete(var_name)
      delete(fAtt)
      delete(str)
      delete(im)
      print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
      print("!")
      print("!  Files are been produced in "+diri_out+"/E06FIG07/  !!")
      print("!")
      print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
   else
      print("OPTION_FILE is "+OPTIONS_FILE+"  NO OUTPUT FILES ARE CREATED ")
   end if
   delete(lat_min)
   delete(lev_sel)
   delete(model1)
   delete(sel)
   delete(dim_MOD)
   delete(lat_name_max)
   delete(lat_name_min)
   delete(year_max)
   delete(dashs)
   delete(month)
   delete(colors)
   delete(v_wgt_sdt)
   delete(res)
   delete(k)
   delete(j)
   delete(Pres)
   delete(lat_max)
   delete(gsres)
   delete(imod)
   delete(v_wgt)
   delete(sel_point)
   delete(year_min)
   delete(Obs_mod)
   delete(modelli)
   ;delete(xp)
   ;delete(yp)
   delete(nv)
   ;delete(OBS)
   delete(idx_ref)
   delete(jmod)
   delete(il)
   delete(varLongName)
   delete(diri_out)
   delete(project_name)
   delete(lev_sel1)

   print ("=================== FINISH (E06FIG07) ====================")

end
;*****************************************************************
