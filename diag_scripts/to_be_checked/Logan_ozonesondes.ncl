;##################################################################################
;# This NCL script reads in CF compliant netCDF CCMVal model output as specified  #
;# in the CCMVal-2 data request available at                                      #
;# http://www.pa.op.dlr.de/CCMVal/DataRequests/CCMVal-2_Datarequest_FINAL.pdf     #
;# and produces FIGURE 11 of                                                      #
;# Cionni, I., V. Eyring, J. F. Lamarque, W. J. Randel, D. S. Stevenson, F. Wu,   #
;# G. E. Bodeker, T. G. Shepherd, D. T. Shindell, and D. W. Waugh, Ozone database #
;# in support of CMIP5 simulations: results and corresponding radiative forcing,  #
;# Atmos. Chem. Phys. Discuss., 11, 10875-10933, doi:10.5194/acpd-11-10875-2011,  #
;# 2011.                                                                          #
;# (hereinafter "Logan_ozone")                                                    #
;##################################################################################
;# Authors and Contact:                                                           #
;#                                                                                #
;# Irene Cionni (irene.cionni@enea.it) & Veronika Eyring (veronika.eyring@dlr.de) #
;##################################################################################
;# The NCL script has been developed for use with the CCMVal diagnostic tool      #
;# (PIs: Andrew Gettelman & Veronika Eyring, co-I: Hisako Shiona)                 #
;##################################################################################
;# The NCL script is called by the Python script main.py                          #
;##################################################################################
;# Description: Comparison of tropospheric {tropo} or stratospheric {strato} ozone#
;# from the models with ozonesonde observations at 40 different sites.            #
;# For each of the stations, a series of plots is shown.                          #
;# (A) First plot: location of the individual station                             #
;# (B) Second plot(s): 1 for one model, 2 for two models, ..., N for N models     #
;# The first column shows the annual cycle in observations from                   #
;# Logan et al. (1999), the second column is the equivalent plot from the model.  #
;# The final column compares all points, with bars indicating the standard        #
;# deviation in the observations. Where the model overpredicts (underpredicts)    #
;# observations by more than one standard deviation, the point is plotted in red  #
;# (blue); these points are shown in the second column by the solid (dashed)      #
;# contours.                                                                      #
;# (C) Third plot: Annual mean vertical profile for all indvidual models and the  #
;# observations, as well as for the multi-model mean                              #
;#                                                                                #
;# The number of panels in the plot is determined by the number of models         #
;# in the namelist file                                                           #
;##################################################################################
;# INPUT/OUTPUT:                                                                  #
;# 1. Reads in                                                                    #
;#    (a) '${WKDIR}/var_att/o3_att.ncl'   (${WKDIR} specified in main namelist)   #
;#    (b) '${WKDIR}/namelist_EO6'         (specification of model names etc.)     #
;#    (c) '${WKDIR}/diag_att/E06Diag.att' (specifies plot type)                   #
;#    (d) Model input files               (as defined in namelists, T3M only!     #
;#    (e) '${WKDIR}/input_data/sondes/'   (Logan files)                           #
;# 2. produces                                                                    #
;#    (a) '${diri_out}/{tropo}{strato}_ozonesonde.ps'       (plot as ps-file)     #
;##################################################################################
; This script takes input from the o3_att.ncl file                                #
; info attributes required:                                                       #
; info@zone             {tropo} for tropospheric comparison                       #
;                       {strato} for stratospheric comparison                     #
; info@period           {ALL}: the program will use all years that are available  #
;                       in the observational database (note: varies from station  #
;                       to station)                                               #
;                       {specific}: all years that are specified (e.g. 1986-1987) #
;                       are used.                                                 #
; *********************************************************************************
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

load "./plot_type/CCMVal_FUNCTION/new_function.ncl"
load "./plot_type/CCMVal_FUNCTION/Climatology_CCMVal.ncl"
load "./plot_type/CCMVal_FUNCTION/misc_function.ncl"
load "./plot_type/get_field_number.ncl"
load "./plot_type/uniform_coordinates.ncl"
; ********************************************************************

begin
   
;************************************************
; Generate some random data
;************************************************
   ; fields supported T3M
   print("Begin ./plot_type/Logan_ozonesondes.ncl")
   diri_out = getenv("ESMValTool_plot_dir")

   ;; The 'models' array has been loaded in main.ncl
   dim_MOD = dimsizes(models@name)

   year_min=stringtointeger(models@start_date)
   year_max=stringtointeger(models@end_date)

   colors=new(dim_MOD+2,string)
   dashs=new(dim_MOD+2,integer)

   MarkLin=new(dim_MOD+2,string)
   modelli=new(dim_MOD+2,string)

   jmod=0
   ;; The 'models' array has been loaded in main.ncl
   do imod=0,dim_MOD-1
      project_name = models@project(imod)
      if (project_name.eq."CCMVal2") then
         colors(jmod)= CCMVal2_ColorTable(models@name(imod))
         dashs(jmod) = CCMVal2_DashTable(models@name(imod))
      else
         colors(jmod)= CCMVal1_ColorTable(models@name(imod))
         dashs(jmod) = CCMVal1_DashTable(models@name(imod))
      end if
      modelli(jmod)=models@name(imod)
      MarkLin(jmod)="Lines"
      jmod=jmod+1
   end do

   MarkLin(dim_MOD)="Lines"
   modelli(dim_MOD)="Multimodel Mean"
   colors(dim_MOD)="Black"
   dashs(dim_MOD)= 0

   MarkLin(dim_MOD+1)="Markers"
   colors(dim_MOD+1) = "Black"
   dashs(dim_MOD+1) = 0

   modelli(dim_MOD+1) ="OzoneSondes"

   filedir  = "plot_type/input_data/SONDES/"
   filenames=systemfunc("ls " + filedir + "sonde*")
   print(filenames)

   fileatt="codes"

   attrib=info@zone

   modelname=models
   modal=info@period
   year_s=1980
   year_e=1995

   wks   = gsn_open_wks ("ps",diri_out+attrib+"_OzoneSondes")               ; open workstation
   map=new(dimsizes(filenames),"graphic")
   plot1=new(2*dimsizes(filenames),"graphic")
   plot2=new(dimsizes(filenames)*dimsizes(modelli),"graphic")
   pol0=new(dimsizes(filenames),"graphic")
   plot2st=new(dimsizes(filenames)*dimsizes(modelli),"graphic")
   plot4=new(dimsizes(filenames)*dimsizes(modelli),"graphic")
   plot4b=new(dimsizes(filenames)*dimsizes(modelli),"graphic")
   plot_diff=new(dimsizes(filenames)*dimsizes(modelli),"graphic")
   plota=new(3*dimsizes(filenames)*dimsizes(modelli),"graphic")
   plota_sel=new(3*3,"graphic")
   plot_xy=new(dimsizes(filenames)+8,"graphic")
   plotu=new(dimsizes(filenames)*dimsizes(modelli),"graphic")
   plotn=new(dimsizes(filenames)*dimsizes(modelli),"graphic")

   lat_modello=new(dimsizes(filenames),"float")
;*************************************************************

   ic=0
   do ik=0,dimsizes(filenames)-3
      print(filenames(ik))
      station_number=asciiread(filenames(ik),1,"integer")
      sonde_att=asciiread(filedir+fileatt,-1,"float")
      sta_sonde_att = onedtond(sonde_att,(/36,4/)) ; npts x ncol

      sonde_number=sta_sonde_att(:,0)
      sonde_lat=sta_sonde_att(:,1)
      sonde_lon=sta_sonde_att(:,2)
      sonde_lon=where(sonde_lon.lt.0,360+sonde_lon,sonde_lon)

      sonde_alt=sta_sonde_att(:,3)

      data1   = asciiread(filedir+fileatt,-1,"string")
      sonde_names_a=(str_get_cols(data1, 27,45))
      sonde_names=sonde_names_a(1:)

      delete(sonde_names_a)
      ind_stat=ind(sonde_number.eq.station_number)

      if ismissing(ind_stat) then
         delete(sta_sonde_att)
         delete(sonde_number)
         delete(sonde_lat)
         delete(sonde_lon)
         delete(sonde_alt)
         delete(data1)
         delete(sonde_names)

         ;fileatt="sonde.codes"
         fileatt="codes"

	 print(filedir+fileatt)
	 sonde_att=asciiread(filedir+fileatt,-1,"float")
         sta_sonde_att = onedtond(sonde_att,(/36,4/)) ; npts x ncol

         sonde_number=sta_sonde_att(:,0)
         sonde_lat=sta_sonde_att(:,1)
         sonde_lon=sta_sonde_att(:,2)
         sonde_lon=where(sonde_lon.lt.0,360+sonde_lon,sonde_lon)

         sonde_alt=sta_sonde_att(:,3)

         data1   = asciiread(filedir+fileatt,-1,"string")
         sonde_names_a=(str_get_cols(data1, 27,45))
         sonde_names=sonde_names_a(1:)
         delete(sonde_names_a)
      end if

      ind_stat=ind(sonde_number.eq.station_number)
      if ismissing(ind_stat) then
         exit
      end if

;//////////////////////////////////////////////////////////
      head = asciiread(filenames(ik),1,"string")
      headcar = stringtochar(head)
      dim_h = dimsizes(headcar)

      month_start= stringtoint(chartostring(headcar(dim_h-15:dim_h-14)))
      year_start= stringtoint(charactertostring(headcar(dim_h-12:dim_h-11)))

      month_end=stringtoint(charactertostring(headcar(dim_h-7:dim_h-6)))
      year_end=stringtoint(charactertostring(headcar(dim_h-4:dim_h-3)))

      if (year_start.le.50)then
         year_start=year_start+2000
      else
         year_start=year_start+1900
      end if

      if (year_end.le.50)then
         year_end=year_end+2000
      else
         year_end=year_end+1900
      end if

      data_O3=asciiread(filenames(ik), -1,"string")

      pressure    = stringtofloat(str_get_field(data_O3,1," "))
      ind_val=ind(.not.ismissing(pressure))
      ind_start=maxind(pressure)
      ind_stop1=minind(pressure)
      ind_stop2=minind(pressure(ind_stop1+1:))
      ind_stop3=minind(pressure(ind_stop1+1+ind_stop2+1:))

      da_O3=new((/dimsizes(pressure),12/),"float")
      do k=0,11
          da_O3(:,k)=stringtofloat(str_get_field(data_O3,2+k," "))
      end do

      plev=pressure(ind_start:ind_stop1)
      plev1=pressure(ind_stop1+2:ind_stop1+1+ind_stop2)
      plev2=pressure(ind_stop1+2+ind_stop2+1:ind_stop1+ind_stop2+ind_stop3+2)

      O3=da_O3(ind_start:ind_stop1,:)
      nO3=da_O3(ind_stop1+2:ind_stop1+1+ind_stop2,:)
      st_O3=da_O3(ind_stop1+2+ind_stop2+1:ind_stop1+ind_stop2+ind_stop3+2,:)
      do il=0,dimsizes(plev)-1
         O3(il,:)=1000*O3(il,:)/plev(il)
         st_O3(il,:)=1000*st_O3(il,:)/plev(il)
      end do

      plev@units="hPa"
      O3!0="plev"
      O3&plev=plev
      O3!1="month"
      O3&month=(/1,2,3,4,5,6,7,8,9,10,11,12/)

      st_O3!0="plev"
      st_O3&plev=plev
      st_O3!1="month"
      st_O3&month=(/1,2,3,4,5,6,7,8,9,10,11,12/)

      ; There are some stations with lat/lon values outside range, so
      ; set these to missing.
      ;
      res             = True                              ; plot mods desired
      res@gsnFrame    = False                             ;  advance yet
      res@gsnMaximize = True
      ;res@tfPolyDrawOrder   = "PostDraw"
      nc1 = NhlNewColor(wks,.8,.8,.8)
      res@tiMainString = "Location of Ozone Sonde station"
      res@mpFillColors = (/-1,-1,nc1,-1/)    ;assign light gray to land masses

      map(ik) = gsn_csm_map(wks,res)                ; Create the map plot, don't draw

      mkres               = True
      mkres@gsMarkerIndex = 2
      mkres@gsMarkerSizeF=0.015
      mkres@gsMarkerThicknessF=2.
      mkres@gsMarkerColor = "green"
      mkres@tfPolyDrawOrder   = "Predraw"
      gsn_polymarker(wks,map,sonde_lon(ind_stat),sonde_lat(ind_stat),mkres)
      txres=True
      txres@txFontHeightF=0.01

      gsn_text(wks,map(ik),sonde_names(ind_stat),sonde_lon(ind_stat)+15,sonde_lat(ind_stat)+6,txres)
      txres@txFontHeightF=0.008

      gsn_text(wks,map(ik),month_start+" "+year_start,sonde_lon(ind_stat)+15,sonde_lat(ind_stat)+1,txres)
      gsn_text(wks,map(ik),month_end+" "+year_end,sonde_lon(ind_stat)+15,sonde_lat(ind_stat)-3,txres)
      ;draw(map8ik))
      frame(wks)
      delete(txres)
      delete(mkres)
      delete(res)
      gsn_define_colormap(wks,"testcmap")                ; choose colormap
      res1                      = True                  ; plot mods desired
      res1@gsnFrame=False
      res1@gsnDraw=False
      res1@tiYAxisFontHeightF          = 0.034
      res1@tmXBLabelFontHeightF        = 0.04
      res1@tmYLLabelFontHeightF        = 0.043
      ;res1@txFontHeightF     = 0.035
      ;res1@tiMainFontHeightF= 0.022
      res1@cnFillOn             = True                  ; turn on color
      res1@cnLinesOn            = False                 ; no contour lines
      res1@gsnSpreadColors      = True                  ; use full range of colors
      res1@lbLabelAutoStride    = True                  ; optimal labels
      res1@tmXBMode             = "Explicit"            ; explicit labels
      res1@gsnLeftString        = "Ozonesonde   lat= "+sonde_lat(ind_stat)+" lon= "+sonde_lon(ind_stat)
      res1@gsnLeftStringOrthogonalPosF=0.05
      res1@gsnRightStringOrthogonalPosF=0.05
      res1@gsnLeftStringFontHeightF=0.03
      res1@gsnRightStringFontHeightF=0.03
      ; res1@gsnRightString        = " lat= "+sonde_lat(ind_stat)+" lon= "+sonde_lon(ind_stat)
      ; res1@tiMainString=sonde_names(ind_stat)+"  "+year_start+" - "+year_end
      res1@tmXTBorderOn=False
      res1@tmXBValues = (/1,2,3,4,5,6,7,8,9,10,11,12/)
      res1@tmXBLabels = (/"J","F","M","A","M","J","J","A","S","O","N","D"/)                ; labels themselves
      ;res1@trYAxisType          = "LogAxis"
      ; res1@vpHeightF            = 0.5                    ; change aspect ratio of plot
      ; res1@vpWidthF             = 0.5
      res1@vpXF             = 0.01                  ; position and sizes
      res1@vpYF             = 0.94                  ; for XY plot
      res1@vpWidthF         = 0.74
      res1@vpHeightF        = 0.74
      ; res1@lbLabelBarOn =True
      ; res1@lbOrientation="Vertical"
      ; res1@lbLabelPosition="Left"
      ; res1@lbLeftMarginF =0.35

      if(attrib.eq."tropo") then
         res1@trYMaxF           = 800.
         res1@trYMinF           = 200.
         res1@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
         res1@cnMinLevelValF       = 10.                  ; set min contour level
         res1@cnMaxLevelValF       = 200.                  ; set max contour level
         res1@cnLevelSpacingF      = 3.                  ; set contour spacing
         ; res1@cnExplicitLabelBarLabelsOn=True
         ; res1@lbLabelAlignment="ExternalEdges"
         ; res1@lbLabelStrings=(/0.,50., 100.,150.,200.,250./)
      end if
      if(attrib.eq."strato") then
         res1@trYMaxF           = 200.
         res1@trYMinF           = 10.
         res1@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
         res1@cnMinLevelValF       = 0.                  ; set min contour level
         res1@cnMaxLevelValF       = 6.                  ; set max contour level
         res1@cnLevelSpacingF      = 0.2
         ; set contour spacing
         st_O3=st_O3*10^-3
         O3=O3*10^-3
      end if
      res1@lbLabelBarOn =False
      res1@cnLinesOn           = False
      res1@cnLineLabelsOn      = False
      res1@cnInfoLabelOn       = False
      ; res@gsnMainString =sonde_names(ind_stat)+" "+month_start+"/"+year_start+" - "+month_end+"/"+year_end

      plot1(2*ik) = gsn_csm_pres_hgt (wks,O3,res1)
      res1@gsnLeftString       = "Standard Deviation Ozonesonde"

      if(attrib.eq."tropo") then
         res1@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
         res1@cnMinLevelValF       = 0.                  ; set min contour level
         res1@cnMaxLevelValF       = 50.                  ; set max contour level
         res1@cnLevelSpacingF      = 3.                  ; set contour spacing
      end if

      if(attrib.eq."strato") then
         res1@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
         res1@cnMinLevelValF       = 0.                  ; set min contour level
         res1@cnMaxLevelValF       = 1.                  ; set max contour level
         res1@cnLevelSpacingF      = 0.02                  ; set contour spacing
      end if
      plot1(2*ik+1)= gsn_csm_pres_hgt (wks,st_O3,res1)

      Pres = True
      Pres@txString ="Climatology at "+sonde_names(ind_stat)+" since "+month_start+"/"+year_start+" to "+month_end+"/"+year_end

;     gsn_panel(wks,plot1(2*ik:2*ik+1),(/1,2/),Pres)
      delete(res1)
      delete(Pres)

;***************************************************************************************
;********************************** Comparison mit model ********************************
;***************************************************************************************

      O3_sel=new((/dim_MOD+1,12,dimsizes(plev)/),"double")
      O3_sel_st=new((/dim_MOD+1,12,dimsizes(plev)/),"double")
      O3_sel!0="model"
      O3_sel!1="month"
      O3_sel!2="plev"
      O3_sel&model=modelli(0:dim_MOD)
      O3_sel&month=(/1,2,3,4,5,6,7,8,9,10,11,12/)
      O3_sel&plev=plev

      O3_sel_st!0="model"
      O3_sel_st!1="month"
      O3_sel_st!2="plev"
      O3_sel_st&model=modelli(0:dim_MOD)
      O3_sel_st&month=(/1,2,3,4,5,6,7,8,9,10,11,12/)
      O3_sel_st&plev=plev

      do im=0,dim_MOD-1
         A0 = extract_data(im,data,-1,0,0)
         print(modelli(im))
         tim=A0&time

         option=0
         option@calendar=tim@calendar
         
	 if (modal.eq."specific") then
            time_start=doubletofloat(cd_inv_calendar(year_s,1,1,0,0,0,tim@units,option))
            time_end=doubletofloat(cd_inv_calendar(year_e,12,30,23,59,59,tim@units,option))
         else
            time_start=doubletofloat(cd_inv_calendar(year_start,month_start,1,0,0,0,tim@units,option))
            time_end=doubletofloat(cd_inv_calendar(year_end,month_end,30,23,59,59,tim@units,option))
         end if
;inserire if

         O3_mod=A0({time_start:time_end},:,:,:)
         delete(A0)
         delete(tim)
         O3_mod=O3_mod


;///////////////level interpolation///////////////////////////////
         pleva=O3_mod&plev
         if (pleva(0).lt.pleva(1)) then
            O3_mod1=O3_mod(:,::-1,:,:)
            delete(O3_mod)
            O3_mod=O3_mod1
            delete(O3_mod1)
         end if
         delete(pleva)

         O3_int=int2p_n(O3_mod&plev,O3_mod,plev,2,1)
         O3_int!0="time"
         O3_int!1="plev"
         O3_int!2="lat"
         O3_int!3="lon"
         O3_int&time=O3_mod&time
         plev@units="hPa"
         O3_int&plev=plev
         O3_int&lat=O3_mod&lat
         O3_int&lon=O3_mod&lon

;//////////////////////Climatology mean and standard deviation/////////////////////////////////////

         dimiO3=dimsizes(O3_int)
         dimiO3_12=dimiO3(0)/12.-floor(dimiO3(0)/12)
         if(dimiO3_12.gt.0.1) then
            dimitime_n=12*ceil(dimiO3(0)/12)
            O3_int1=new((/floattointeger(dimitime_n)+12,dimsizes(O3_int&plev),dimsizes(O3_int&lat),dimsizes(O3_int&lon)/),typeof(O3_int))
            
	    printVarSummary(O3_int1)
            printVarSummary(O3_int1(month_start-1:floattointeger(dimitime_n)+month_end-1,:,:,:))
            O3_int1(month_start-1:floattointeger(dimitime_n)+month_end-1,:,:,:)=O3_int
	    if (month_start.gt.1) then
               O3_int1(0:month_start-1,:,:,:)=O3_int@_FillValue
	    
	    end if
            if (month_end.lt.12) then
               O3_int1(dimitime_n+month_end:,:,:,:)=O3_int@_FillValue

            end if
            O3_int!0="time"
            O3_int1!0="time"
            O3_int1!1="plev"
            O3_int1!2="lat"
            O3_int1!3="lon"
            O3_int1&plev=plev
            O3_int1&lat=O3_mod&lat
            O3_int1&lon=O3_mod&lon

            O3_clim=clmMonTLLL(O3_int1)
            O3_std=stdMonTLLL(O3_int1)
            delete(O3_int1)
         else

            O3_clim=clmMonTLLL(O3_int)
            O3_std=stdMonTLLL(O3_int)
         end if

;/////////////////////////bilinear interpolation/////////////////////////////////////

         lat=O3_clim&lat
         if (lat(0).gt.lat(1)) then
            O3_clim1=O3_clim(:,:,::-1,:)
            delete(O3_clim)
            O3_clim=O3_clim1
            delete(O3_clim1)
            O3_std1=O3_std(:,:,::-1,:)
            delete(O3_std)
            O3_std=O3_std1
            delete(O3_std1)
         end if

         O3_sel(im,:,:) = (/linint2_points(O3_clim&lon,O3_clim&lat,O3_clim, False, sonde_lon(ind_stat),sonde_lat(ind_stat), 0)/)
         O3_sel_st(im,:,:)= (/linint2_points(O3_std&lon,O3_std&lat,O3_std, False, sonde_lon(ind_stat),sonde_lat(ind_stat), 0)/)

;////////////////////////////////////////////////////////////////////////////////////
         lat=O3_clim&lat
         lon=O3_clim&lon
         ind_lat_mod=ind_nearest_coord (sonde_lat(ind_stat),lat, 0)
         ind_lon_mod=ind_nearest_coord (sonde_lon(ind_stat),lon, 0)

         delete(O3_clim)
         delete(O3_std)
         delete(O3_int)

         d_O3=  O3_sel(model|im,plev|:,month|:)-O3(:,:)
         d_O3_n=1.d*O3
         d_O3_n!0="plev"
         d_O3_n!1="month"
         d_O3_n&plev=plev

         st_O3=where(st_O3.eq.0,st_O3@_FillValue,st_O3)
         d_O3_n=where((d_O3.gt.st_O3(:,:)),d_O3/st_O3(:,:),0.)
         d_O3_n=where((d_O3.lt.-st_O3(:,:)),d_O3/st_O3(:,:),d_O3_n)

         if (all(d_O3_n).eq.0.) then
            d_O3_n(0,0)=1.
         end if

;/////////////////////////////plot//////////////////////////////////////////

         ; gsn_define_colormap(wks,"testcmap")                ; choose; colormap
         res2                      = True                  ; plot mods desired
         res2@gsnFrame=False
         res2@gsnDraw=False
         res2@cnFillOn             = True                  ; turn on color
         res2@cnLinesOn            = False                 ; no contour lines
         res2@gsnSpreadColors      = True                  ; use full range of colors
         res2@lbLabelAutoStride    = True                  ; optimal labels
         res2@tmXBMode             = "Explicit"            ; explicit labels
         res2@gsnLeftString        = modelli(im)
         res2@gsnRightString        = " lat= "+lat(ind_lat_mod)+" lon= "+lon(ind_lon_mod)
         ;res2@tiMainString=sonde_names(ind_stat)+"  "+year_start+" - "+year_end
         res2@tiMainFontHeightF=0.06
         res2@tiYAxisFontHeightF          = 0.034
         res2@tmXBLabelFontHeightF        = 0.04
         res2@tmYLLabelFontHeightF        = 0.043
         res2@gsnLeftStringOrthogonalPosF=0.05
         res2@gsnRightStringOrthogonalPosF=0.05
         res2@gsnLeftStringFontHeightF=0.02
         res2@gsnRightStringFontHeightF=0.02

         ;res2@gsnLeftStringParallelPosF=-0.1
         ;res2@gsnLeftRightParallelPosF=0.1
         ;lat_modello(ik)=lat(ind_lat_mod)
         res2@lbLabelBarOn =False

         res2@tmXBValues = (/0,1,2,3,4,5,6,7,8,9,10,11/)
         res2@tmXBLabels = (/"J","F","M","A","M","J","J","A","S","O","N","D"/)                ; labels themselves
         ;res1@trYAxisType          = "LogAxis"
         ;res2@vpHeightF            = 0.5                   ; change aspect ratio of plot
         ;res2@vpWidthF             = 0.5
         res2@vpXF             = 0.33                  ; position and sizes
         res2@vpYF             = 0.94                  ; for XY plot
         res2@vpWidthF         = 0.74
         res2@vpHeightF        = 0.74

         if(attrib.eq."tropo") then
            res2@trYMaxF           = 800.
            res2@trYMinF           = 200.
            res2@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
            res2@cnMinLevelValF       = 10.                  ; set min contour level
            res2@cnMaxLevelValF       = 200.                  ; set max contour level
            res2@cnLevelSpacingF      = 3.                  ; set contour spacing
            O3_sel_st(im,:,:)=O3_sel_st(im,:,:)*10^3
            O3_sel(im,:,:)=O3_sel(im,:,:)*10^3
         end if
         if(attrib.eq."strato") then
            res2@trYMaxF           = 200.
            res2@trYMinF           = 10.
            res2@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
            res2@cnMinLevelValF       = 0.                  ; set min contour level
            res2@cnMaxLevelValF       = 6.                  ; set max contour level
            res2@cnLevelSpacingF      = 0.2
            ; set contour spacing
         end if

         res2@cnLinesOn           = False
         res2@cnLineLabelsOn      = False
         res2@cnInfoLabelOn       = False
         res2@lbLabelBarOn =True
         res2@lbOrientation="Vertical"
         res2@lbLabelPosition="Left"
         res2@lbLeftMarginF =0.35
         plot2(ic) = gsn_csm_pres_hgt (wks,O3_sel(model|im,plev|:,month|:),res2)
         res2@gsnLeftString       = "Standard Deviation "+modelname
         res2@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
         if(attrib.eq."tropo") then

            res2@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
            res2@cnMinLevelValF       = 0.                  ; set min contour level
            res2@cnMaxLevelValF       = 50.                  ; set max contour level
            res2@cnLevelSpacingF      = 3.                  ; set contour spacing
         end if
         if(attrib.eq."strato") then

            res2@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
            res2@cnMinLevelValF       = 0.                  ; set min contour level
            res2@cnMaxLevelValF       = 1.                  ; set max contour level
            res2@cnLevelSpacingF      = 0.02                  ; set contour spacing
         end if

         plot2st(ic)= gsn_csm_pres_hgt (wks,O3_sel_st(model|im,plev|:,month|:),res2)

         ;Pres = True
         ;Pres@txString ="Climatology at "+sonde_names(ind_stat)+" since "+month_start+"/"+year_start+" to "+month_end+"/"+year_end

         ;gsn_panel(wks,plot2(ic:ic+1),(/1,2/),Pres)

         resxy                   = True                     ; plot mods desired
         resxy@gsnFrame=False
         resxy@gsnDraw=False
         ;resxy@vpXF             = 0.5                  ; position and sizes
         ;resxy@vpYF             = 0.94                  ; for XY plot
         resxy@vpWidthF         = 0.6
         resxy@vpHeightF        = 0.6

         resxy@gsnMaximize       = True                     ; maximize plot
         ;resxy@tiMainString      = sonde_names(ind_stat)+" "+month_start+"/"+year_start+" - "+month_end+"/"+year_end

         resxy@xyMarkLineMode    = "Markers"                ; choose to use markers
         resxy@xyMarkers         = (/4,4,16/)                  ; choose type of marker
         resxy@xyMarkerColors     = (/"Purple","red","ForestGreen"/)            ; Marker color
         resxy@xyMarkerSizeF     = 0.01                     ; Marker size (default 0.01)
         resxy@tiXAxisString="O3 OzoneSondes  [ ppbv ]"
         resxy@tiYAxisString="O3 "+modelli(im)+" [ ppbv ]"
         resxy@tiXAxisFontHeightF=0.035
         resxy@tiYAxisFontHeightF=0.035

         ;resxy@tiMainString=sonde_names(ind_stat)+"  "+year_start+" - "+year_end

         if (attrib.eq."tropo") then
            O3_a=ndtooned(O3(month|:,{plev|800.:200.}))
            O3_b=ndtooned(O3_sel(model|im,month|:,{plev|800.:200.}))
            O3_d=new((/3,dimsizes(O3_a)/),typeof(O3_a))
            d_O3_nd=ndtooned(d_O3_n(month|:,{plev|800:200}))
            O3_d(0,:)=where(d_O3_nd.lt.0,O3_a,O3_a@_FillValue)
            O3_d(1,:)=where(d_O3_nd.gt.0,O3_a,O3_a@_FillValue)
            O3_d(2,:)=where(d_O3_nd.eq.0,O3_a,O3_a@_FillValue)
         end if
         if (attrib.eq."strato") then
            O3_a=ndtooned(O3(month|:,{plev|200:10}))
            O3_b=ndtooned(O3_sel(model|im,month|:,{plev|200:10}))
            O3_d=new((/3,dimsizes(O3_a)/),typeof(O3_a))
            d_O3_nd=ndtooned(d_O3_n(month|:,{plev|200:10}))
            O3_d(0,:)=where(d_O3_nd.lt.0,O3_a,O3_a@_FillValue)
            O3_d(1,:)=where(d_O3_nd.gt.0,O3_a,O3_a@_FillValue)
            O3_d(2,:)=where(d_O3_nd.eq.0,O3_a,O3_a@_FillValue)
         end if

         plot_xy(ic)=gsn_csm_xy (wks,O3_d,O3_b,resxy)

         if (attrib.eq."tropo") then
            err=new((/2,dimsizes(ndtooned(O3(month|:,{plev|800.:200.})))/),typeof(O3_sel))
            err(0,:)=ndtooned(O3(month|:,{plev|800.:200.}))-ndtooned(st_O3(month|:,{plev|800.:200.}))
            err(1,:)=ndtooned(O3(month|:,{plev|800.:200.}))+ndtooned(st_O3(month|:,{plev|800.:200.}))
            aerr=new((/2,dimsizes(ndtooned(O3(month|:,{plev|800.:200.})))/),typeof(O3_sel))
            aerr(0,:)=ndtooned(O3_sel(model|im,month|:,{plev|800.:200.}))
            aerr(1,:)=ndtooned(O3_sel(model|im,month|:,{plev|800.:200.}))

            do ip=0,dimsizes(ndtooned(O3(month|:,{plev|800:200})))-1
               if (.not.any(ismissing(err(:,ip))).and.(.not.any(ismissing(aerr(:,ip))))) then
                  if (.not.ismissing(O3_d(2,ip))) then
                     varstring=unique_string("plot_xy")
                     resxy1=True
                     resxy1@gsLineColor     = "ForestGreen"            ; Marker color
                     plot_xy@$varstring$ = gsn_add_polyline(wks,plot_xy(ic),err(:,ip),aerr(:,ip),resxy1)
                  end if
                  if (.not.ismissing(O3_d(1,ip))) then
                     varstring=unique_string("plot_xy")
                     resxy1=True
                     resxy1@gsLineColor     = "red"            ; Marker color
                     plot_xy@$varstring$ = gsn_add_polyline(wks,plot_xy(ic),err(:,ip),aerr(:,ip),resxy1)
                  end if
                  if (.not.ismissing(O3_d(0,ip))) then
                     varstring=unique_string("plot_xy")
                     resxy1=True
                     resxy1@gsLineColor     = "Purple"            ; Marker color
                     plot_xy@$varstring$ = gsn_add_polyline(wks,plot_xy(ic),err(:,ip),aerr(:,ip),resxy1)
                  end if
               end if
            end do
         end if
         if (attrib.eq."strato") then
            err=new((/2,dimsizes(ndtooned(O3(month|:,{plev|200:10})))/),typeof(O3_sel))
            err(0,:)=ndtooned(O3(month|:,{plev|200.:10.}))-ndtooned(st_O3(month|:,{plev|200.:10.}))
            err(1,:)=ndtooned(O3(month|:,{plev|200.:10.}))+ndtooned(st_O3(month|:,{plev|200.:10.}))
            aerr=new((/2,dimsizes(ndtooned(O3(month|:,{plev|200.:10.})))/),typeof(O3_sel))
            aerr(0,:)=ndtooned(O3_sel(model|im,month|:,{plev|200.:10.}))
            aerr(1,:)=ndtooned(O3_sel(model|im,month|:,{plev|200.:10.}))

            do ip=0,dimsizes(ndtooned(O3(month|:,{plev|200:10})))-1
               if (.not.any(ismissing(err(:,ip))).and.(.not.any(ismissing(aerr(:,ip))))) then
                  if (.not.ismissing(O3_d(2,ip))) then
                     varstring=unique_string("plot_xy")
                     resxy1=True
                     resxy1@gsLineColor     = "ForestGreen"            ; Marker color
                     plot_xy@$varstring$ = gsn_add_polyline(wks,plot_xy(ic),err(:,ip),aerr(:,ip),resxy1)
                  end if
                  if (.not.ismissing(O3_d(1,ip))) then
                     varstring=unique_string("plot_xy")
                     resxy1=True
                     resxy1@gsLineColor     = "red"            ; Marker color
                     plot_xy@$varstring$ = gsn_add_polyline(wks,plot_xy(ic),err(:,ip),aerr(:,ip),resxy1)
                  end if
                  if (.not.ismissing(O3_d(0,ip))) then
                     varstring=unique_string("plot_xy")
                     resxy1=True
                     resxy1@gsLineColor     = "Purple"            ; Marker color
                     plot_xy@$varstring$ = gsn_add_polyline(wks,plot_xy(ic),err(:,ip),aerr(:,ip),resxy1)
                  end if
               end if
            end do
         end if

         ;draw(plot_xy(ik))

 ;////////////////////diff//////////////////////////////////////////////////
         ;frame(wks)

         reso              = True                         ; plot mods desired
         reso@gsnDraw      = False                        ; don't draw yet
         reso@gsnFrame     = False                        ; don't advance frame yet
         ;reso@tiYAxisString =" "
         ;reso@cnLinesOn           = False
         ;reso@cnLineLabelsOn      = False
         reso@cnInfoLabelOn       = False
         reso@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
         reso@cnMinLevelValF       = -10.                  ; set min contour level
         reso@cnMaxLevelValF       = 10.                  ; set max contour level
         reso@cnLevelSpacingF      = 1.                  ; set contour spacing

         reso@gsnContourZeroLineThicknessF = 2. 	  ; doubles thickness of zero contour
         reso@gsnContourNegLineDashPattern = 1 		  ; sets negative contours to dash pattern 1
         reso@cnLineColor  = "Black"                        ; color of first contour
         d_O3_n!0="plev"
         d_O3_n!1="month"
         d_O3_n&plev=plev

         plotu(ic)  = gsn_csm_pres_hgt(wks, d_O3_n(plev|:,month|:), reso )         ; draw 1st contour
;**

         getvalues plotu@contour
            "pmAnnoManagers" : am_ids
         end getvalues
         index = ind(NhlName(am_ids).eq."right_axis")
         if(.not.ismissing(index)) then
            NhlRemoveAnnotation(plotu@contour,am_ids(index))
         end if

         overlay(plot2(ic),plotu(ic))

 ;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         plota(3*ic)=plot1(ik)
         plota(3*ic+1)=plot2(ic)
         plota(3*ic+2)=plot_xy(ic)
         Pres3 = True
         ;Pres3@HeightF =0.01
         Pres3@txString =sonde_names(ind_stat)+" "+month_start+"/"+year_start+" - "+month_end+"/"+year_end
         ;titles(ik)=sonde_names(ind_stat)+" "+month_start+"/"+year_start+" - "+month_end+"/"+year_end
         Pres3@gsnPanelScalePlotIndex =2
         x_ndc=(/0.08,0.38,0.74/)
         ;Pres3@gsnPanelYWhiteSpacePercent = 10
         Pres3@gsnPanelXWhiteSpacePercent = 13
         ;Pres3@gsnPanelTop    = 0.8
         ;Pres3@gsnPanelBottom =  .2
         Pres3@gsnPanelXF=x_ndc
         Pres3@gsnPanelDebug=True
         Pres3@pmLabelBarOrthogonalPosF=0.23
         Pres3@pmLabelBarParallelPosF=-0.2
         ;gsn_panel(wks,plot_fin,(/6,6,6,6,6/),pres)
         Pres3@lbLabelAutoStride    = True
         gsn_panel(wks,plota(3*ic:3*ic+2),(/1,3/),Pres3)
         ic=ic+1
         delete(O3_a)
         delete(O3_b)
         delete(O3_d)
         delete(O3_mod)
         delete(lat)
         delete(lon)
      end do  ;fine loop sui modelli

      O3_sel(dim_MOD,:,:)=dim_avg(O3_sel(month|:,plev|:,model|:))
;/////////////////////////////COMBINATION PLOTS///////////////////////////

;//////////ANNUAL VERTICAL PROFILE/////////////////////////////////
      O3_ob_mod=new((/dim_MOD+2,dimsizes(plev)/),"double",O3_sel@_FillValue)
      O3_ob_mod@units="ppb"
      O3_ob_mod@long_name="Ozone [ppb]"

      O3_ob_m=dim_avg(O3)
      O3_ob_max=dim_max(O3)
      O3_ob_min=dim_min(O3)

      O3_mod_m=dim_avg_Wrap(O3_sel(model|:,plev|:,month|:))
      O3_mod_max=O3_mod_m
      O3_mod_min=O3_mod_m
      O3_mod_max=dim_max(O3_sel(model|:,plev|:,month|:))
      O3_mod_min=dim_min(O3_sel(model|:,plev|:,month|:))

      O3_ob_mod(dim_MOD+1,:)=(/O3_ob_m/)
      O3_ob_mod(0:dim_MOD,:)=(/O3_mod_m/)

      xp = new((/dim_MOD+1,2*dimsizes(plev)/),typeof(O3_mod_m))
      yp = new((/dim_MOD+1,2*dimsizes(plev)/),typeof(O3_mod_m))

      do j=0,dimsizes(plev)-1
         xp(:,j) = O3_mod_min(:,j)
         yp(:,j) = plev(j)
         yp(:,dimsizes(plev)+j) = plev(dimsizes(plev)-1-j)
         xp(:,dimsizes(plev)+j) = (O3_mod_max(:,(dimsizes(plev)-1-j)))
      end do
      xp1 = new((/2*dimsizes(plev)/),typeof(O3_mod_m))
      yp1 = new((/2*dimsizes(plev)/),typeof(O3_mod_m))

      do j=0,dimsizes(plev)-1
         xp1(j) = O3_ob_min(j)
         yp1(j) = plev(j)
         yp1(dimsizes(plev)+j) = plev(dimsizes(plev)-1-j)
         xp1(dimsizes(plev)+j) = (O3_ob_max((dimsizes(plev)-1-j)))
      end do
      ; gsn_define_colormap(wks,"testmap")

      res3                   = True
      res3@gsnDraw           = False
      res3@gsnFrame          = False
      res3@tfPolyDrawOrder   = "Predraw"

      if(attrib.eq."tropo") then
         res3@trYMaxF           = 1000.
         res3@trYMinF           = 200.
         res3@tmYLLabels        = (/200,100,50,30,20,10/)
         res3@tiXAxisString     = "O3 [ppm]"       ; y-axis label
         res3@trXMaxF           = 10
         res3@trXMinF           =0
      end if

      res3@tmLabelAutoStride = True
      res3@vpHeightF         = 0.6                    ; change aspect ratio of plot
      res3@vpWidthF          = 0.4
      res3@tiYAxisString     = "Pressure [hPa]"       ; y-axis label
      res3@trYAxisType       = "LogAxis"
      res3@trYReverse        = True
      res3@tiMainFontHeightF =0.024
      res3@gsnRightStringFontHeightF=0.012
      res3@gsnLeftStringFontHeightF=0.012
      ;res3@gsnLeftStringOrthogonalPosF=-0.52
      ;res3@gsnRightStringOrthogonalPosF=-0.52

      res3@tiMainString="Annual mean at "+sonde_names(ind_stat)
      res3@gsnLeftString=month_start+"/"+year_start+" to "+month_end+"/"+year_end
      res3@gsnRightString="lat= "+sonde_lat(ind_stat)+" lon= "+sonde_lon(ind_stat)

      ; res@xyMarkerSizeF=0.1
      res3@xyLineColors      = colors                 ; line color
      res3@xyDashPatterns    = dashs
      res3@xyLineThicknessF =3.
      scattlines=colors
      scattlines(0:dim_MOD)= "Lines"

      scattlines(dim_MOD+1)="Markers"
      res3@xyMarkLineModes    = scattlines

      res3@pmLegendDisplayMode    = "Always"            ; turn on legend
      res3@pmLegendSide           = "Top"               ; Change location of
      res3@pmLegendParallelPosF   = 0.7                 ; move units right
      res3@pmLegendOrthogonalPosF = -0.8               ; move units down
      res3@pmLegendWidthF         = 0.1                ; Change width and
      res3@pmLegendHeightF        = 0.1                ; height of legend.
      res3@lgPerimOn              = False               ; turn off box around
      res3@lgLabelFontHeightF     = .012                 ; label font height
      ; res3@tmXBLabelDeltaF =2
      res3@xyExplicitLegendLabels =modelli         ; create explicit labels
      plot_diff(ik) = gsn_csm_xy (wks,O3_ob_mod,plev,res3)

      gsres = True
      gsres@gsFillColor = "Snow3"
      pol0(ik)=gsn_add_polygon(wks,plot_diff(ik),xp1(:),yp1(:),gsres)
      ; gsres1 = True
      ; gsres1@gsFillColor        ="(/0.105,0.675,0.8/)"
      ; gsres1@gsFillScaleF=0.4
      ; gsres1@gsFillIndex=17
      ; pol(ik)=gsn_add_polygon(wks,plot_diff(ik),xp(:),yp(:),gsres1)

      draw(plot_diff(ik))

      delete(headcar)
      print(ik)
      frame(wks)
   end do    ;fine loop stazioni
   print("Figures has been produced in "+diri_out+attrib+"_OzoneSondes" )

end
