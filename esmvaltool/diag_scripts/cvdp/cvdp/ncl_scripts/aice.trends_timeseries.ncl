; Calculates SIC hemispheric trends and extent
;
; Variables used: sic
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: aice.trends_timeseries.ncl")
  
  SCALE_TIMESERIES = getenv("SCALE_TIMESERIES")  
  OUTPUT_DATA      = getenv("OUTPUT_DATA")  
  PNG_SCALE        = tofloat(getenv("PNG_SCALE"))
  OPT_CLIMO        = getenv("OPT_CLIMO")
  CLIMO_SYEAR      = toint(getenv("CLIMO_SYEAR"))
  CLIMO_EYEAR      = toint(getenv("CLIMO_EYEAR"))
  OUTPUT_TYPE      = getenv("OUTPUT_TYPE") 
  COLORMAP         = getenv("COLORMAP")
  
  nsim = numAsciiRow("namelist_byvar/namelist_aice_nh")
  na = asciiread("namelist_byvar/namelist_aice_nh",(/nsim/),"string")
  names = new(nsim,"string")
  paths = new(nsim,"string")
  syear = new(nsim,"integer",-999)
  eyear = new(nsim,"integer",-999)
  delim = "|"

  do gg = 0,nsim-1
     names(gg) = str_strip(str_get_field(na(gg),1,delim))
     paths(gg) = str_strip(str_get_field(na(gg),2,delim))
     syear(gg) = stringtointeger(str_strip(str_get_field(na(gg),3,delim)))
     eyear(gg) = stringtointeger(str_strip(str_get_field(na(gg),4,delim)))
  end do
  nyr = eyear-syear+1
  nyr_max = max(nyr)

  nsim_sh = numAsciiRow("namelist_byvar/namelist_aice_sh")
  na_sh = asciiread("namelist_byvar/namelist_aice_sh",(/nsim/),"string")
  names_sh = new(nsim,"string")
  paths_sh = new(nsim,"string")
  syear_sh = new(nsim,"integer",-999)
  eyear_sh = new(nsim,"integer",-999)
  do gg = 0,nsim-1
     names_sh(gg) = str_strip(str_get_field(na_sh(gg),1,delim))
     paths_sh(gg) = str_strip(str_get_field(na_sh(gg),2,delim))
     syear_sh(gg) = stringtointeger(str_strip(str_get_field(na_sh(gg),3,delim)))
     eyear_sh(gg) = stringtointeger(str_strip(str_get_field(na_sh(gg),4,delim)))
  end do
  nyr_sh = eyear_sh-syear_sh+1
  nyr_max_sh = max(nyr_sh)
 
  wks_type = OUTPUT_TYPE
  if (wks_type.eq."png") then
     wks_type@wkWidth = 1500*PNG_SCALE
     wks_type@wkHeight = 1500*PNG_SCALE
  end if
  wks_trends_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.trends.djf")
  wks_trends_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.trends.mam")
  wks_trends_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.trends.jja")
  wks_trends_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.trends.son")
  wks_trends_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.trends.ann")
  wks_trends_mon = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.trends.mon")

  wks_iceext_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.extent.djf")
  wks_iceext_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.extent.mam")
  wks_iceext_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.extent.jja")
  wks_iceext_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.extent.son")
  wks_iceext_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.extent.ann")
  wks_iceext_febmar = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.extent.febmar")
  wks_iceext_sep = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.extent.sep")
  wks_iceext_mon = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.extent.mon")
 
  if (COLORMAP.eq.0) then
     gsn_define_colormap(wks_trends_djf,"ncl_default")   
     gsn_define_colormap(wks_trends_mam,"ncl_default")  
     gsn_define_colormap(wks_trends_jja,"ncl_default") 
     gsn_define_colormap(wks_trends_son,"ncl_default") 
     gsn_define_colormap(wks_trends_ann,"ncl_default") 
     gsn_define_colormap(wks_trends_mon,"ncl_default") 
  end if
  if (COLORMAP.eq.1) then
     gsn_define_colormap(wks_trends_djf,"BlueDarkRed18")     
     gsn_define_colormap(wks_trends_mam,"BlueDarkRed18")    
     gsn_define_colormap(wks_trends_jja,"BlueDarkRed18")   
     gsn_define_colormap(wks_trends_son,"BlueDarkRed18")   
     gsn_define_colormap(wks_trends_ann,"BlueDarkRed18")   
     gsn_define_colormap(wks_trends_mon,"BlueDarkRed18")   
  end if

  plot_trends_nh_djf = new(nsim,"graphic")  
  plot_trends_nh_mam = new(nsim,"graphic")  
  plot_trends_nh_jja = new(nsim,"graphic")  
  plot_trends_nh_son = new(nsim,"graphic")   
  plot_trends_nh_ann = new(nsim,"graphic")
  plot_trends_nh_mon = new(nsim,"graphic")  
  plot_iceext_nh_djf = new(nsim,"graphic")  
  plot_iceext_nh_mam = new(nsim,"graphic")  
  plot_iceext_nh_jja = new(nsim,"graphic")  
  plot_iceext_nh_son = new(nsim,"graphic")   
  plot_iceext_nh_ann = new(nsim,"graphic")
  plot_iceext_nh_mon = new(nsim,"graphic")
  plot_iceext_nh_mon_anom = new(nsim,"graphic")
  plot_iceext_nh_feb = new(nsim,"graphic")
  plot_iceext_nh_mar = new(nsim,"graphic")
  plot_iceext_nh_sep = new(nsim,"graphic")
  plot_iceext_nh_climo = new(nsim,"graphic")

  plot_trends_sh_djf = new(nsim,"graphic")  
  plot_trends_sh_mam = new(nsim,"graphic")  
  plot_trends_sh_jja = new(nsim,"graphic")  
  plot_trends_sh_son = new(nsim,"graphic")   
  plot_trends_sh_ann = new(nsim,"graphic")
  plot_trends_sh_mon = new(nsim,"graphic")
  plot_iceext_sh_djf = new(nsim,"graphic")  
  plot_iceext_sh_mam = new(nsim,"graphic")  
  plot_iceext_sh_jja = new(nsim,"graphic")  
  plot_iceext_sh_son = new(nsim,"graphic")   
  plot_iceext_sh_ann = new(nsim,"graphic")
  plot_iceext_sh_mon = new(nsim,"graphic")
  plot_iceext_sh_mon_anom = new(nsim,"graphic")
  plot_iceext_sh_feb = new(nsim,"graphic")
  plot_iceext_sh_mar = new(nsim,"graphic")
  plot_iceext_sh_sep = new(nsim,"graphic")
  plot_iceext_sh_climo = new(nsim,"graphic")

  if (isfilepresent2("obs_aice_nh")) then
     plot_iceext_nh_obs = new((/6,nsim/),"graphic")  
  end if
  if (isfilepresent2("obs_aice_sh")) then
     plot_iceext_sh_obs = new((/6,nsim/),"graphic")  
  end if

  time_mon2 = ispan(0,11,1)
  time_mon2@units = "months since 0000-01-01 00:00:00"
  time_mon2@long_name = "Time"
  time_mon2@standard_name = "time"
  time_mon2@calendar = "standard"
  time_mon2!0 = "time_mon2"
  time_mon2&time_mon2 = time_mon2

  do ee = 0,nsim-1
     aice_nh_flag = 0
     aice_nh = data_read_in_ice(paths(ee),"aice_nh",syear(ee),eyear(ee))    ; read in data, orient lats/lons correctly, set time coordinate variable up
     if (isatt(aice_nh,"is_all_missing")) then
        delete(aice_nh)
        aice_nh_flag = 1
     end if  

     ph_s = ""  ; flag that will be used to denote if pole hole area filled in via pole_hole_area attribute
     if (aice_nh_flag.eq.0) then
        if (isatt(aice_nh,"area")) then
           area3d = aice_nh@area
           area3d_c = conform(aice_nh,area3d,(/1,2/))
           aice_nh_sum = aice_nh
;           aice_nh_sum = (/ (where((aice_nh/100.).ge.0.15,(aice_nh/100.),aice_nh@_FillValue))*area3d_c /) ;   ice area calculation (all cells > 15% kept) 
           aice_nh_sum = (/ where(aice_nh.ge.15,1.,aice_nh@_FillValue) /) ; ice extent calculation (all cells greater than 15% treated as 100% covered)
           aice_nh_sum = aice_nh_sum*area3d_c
           wgts = aice_nh_sum(0,:,:)
           wgts = 1.
           aice_nh_sum_mon = aice_nh_sum(:,0,0)   ; preallocate array to retain metadata
           aice_nh_sum_mon = (/ wgt_areasum2(aice_nh_sum,wgts,0) /)
;           do gg = 0,dimsizes(aice_nh&time)-1
;              aice_nh_sum_mon(gg) = (/ sum(aice_nh_sum(gg,:,:)) /)
;           end do

           if (isatt(aice_nh,"pole_hole_area")) then   ; special attribute set up to account for pole hole in grids. NSIDC assumes hole is 100% ice covered as of Jan 2016 
              ph_area = todouble(aice_nh_sum_mon@pole_hole_area)  ; format: start YYYYMM, end YYYYMM, area, start YYYYMM, end YYYYMM, area, etc.
              dimZ_ph = dimsizes(ph_area)/3                       ; only used for Northern Hemisphere
              temp_area_arr = aice_nh_sum_mon
              temp_area_arr&time = cd_calendar(aice_nh_sum_mon&time,1)
;              printVarSummary(temp_area_arr)
;              print(ph_area)
              do gg = 0,dimZ_ph-1
                 temp_area_arr({ph_area(gg*3):ph_area(gg*3+1)}) = temp_area_arr({ph_area(gg*3):ph_area(gg*3+1)}) + tofloat(ph_area(gg*3+2))
;                 print(ph_area(gg*3)+" "+ph_area(gg*3+1)+" "+tofloat(ph_area(gg*3+2)))
              end do
              aice_nh_sum_mon = (/ temp_area_arr /)
              delete([/ph_area,dimZ_ph,temp_area_arr/])
              ph_s = "*"
           end if                

           aice_nh_sum_mon = aice_nh_sum_mon/1.e12
           aice_nh_sum_mon@units = "10^12 m2"
           aice_nh_sum_mon@long_name = "sea_ice_extent"
           if (isatt(aice_nh_sum_mon,"coordinates")) then
              delete(aice_nh_sum_mon@coordinates)
           end if
           delete([/aice_nh_sum,area3d,area3d_c,wgts/])
          
           taice = new((/dimsizes(aice_nh_sum_mon),1,1/),typeof(aice_nh_sum_mon))
           taice!0 = "time"
           taice!1 = "lat"
           taice!2 = "lon"
           taice(:,0,0) = (/ aice_nh_sum_mon /) ; convert to 3D array so we can use clmMonTLL and calcMonAnomTLL
           if (OPT_CLIMO.eq."Full") then
              taice = rmMonAnnCycTLL(taice)
           else
              check_custom_climo(names(ee),syear(ee),eyear(ee),CLIMO_SYEAR,CLIMO_EYEAR)
              temp_arr = taice
              temp_arr&time = cd_calendar(aice_nh_sum_mon&time,-1)
              if (CLIMO_SYEAR.lt.0) then
                 climo = clmMonTLL(temp_arr({(eyear(ee)+CLIMO_SYEAR)*100+1:(eyear(ee)+CLIMO_EYEAR)*100+12},:,:))       
              else
                 climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
              end if              
              delete(temp_arr)
              taice   = calcMonAnomTLL(taice,climo) 
              delete(climo)
           end if
           aice_nh_sum_mon_anom = aice_nh_sum_mon
           aice_nh_sum_mon_anom = (/ taice(:,0,0) /)
           aice_nh_sum_mon_anom@long_name = "sea_ice_extent_anomaly"
           delete(taice)

           aice_nh_sum_feb = aice_nh_sum_mon(1::12)
           delete(aice_nh_sum_feb&time)
           aice_nh_sum_feb@long_name = "February sea_ice_extent"
           aice_nh_sum_feb!0 = "TIME"
           aice_nh_sum_feb&TIME = ispan(syear(ee),eyear(ee),1)
           aice_nh_sum_feb&TIME@units = "YYYY"
           aice_nh_sum_feb&TIME@long_name = "time"
           aice_nh_sum_mar = aice_nh_sum_mon(2::12)
           aice_nh_sum_mar@long_name = "March sea_ice_extent"
           copy_VarCoords(aice_nh_sum_feb,aice_nh_sum_mar)
           aice_nh_sum_sep = aice_nh_sum_mon(8::12)
           aice_nh_sum_sep@long_name = "September sea_ice_extent"
           copy_VarCoords(aice_nh_sum_feb,aice_nh_sum_sep)

           aice_nh_sum_climo = new(12,typeof(aice_nh_sum_mon))
           copy_VarAtts(aice_nh_sum_mon,aice_nh_sum_climo)
           aice_nh_sum_climo@long_name = "climatological_sea_ice_extent"
           aice_nh_sum_climo!0 = "time_mon2"
           aice_nh_sum_climo&time_mon2 = time_mon2
           do gg = 0,11
              aice_nh_sum_climo(gg) = (/ avg(aice_nh_sum_mon(gg::12)) /)
           end do    

           temp = runave_Wrap(aice_nh_sum_mon,3,0)
           aice_nh_sum_djf = temp(0::12)
           aice_nh_sum_mam = temp(3::12)
           aice_nh_sum_jja = temp(6::12)
           aice_nh_sum_son = temp(9::12)
           delete(temp)
           temp = runave_Wrap(aice_nh_sum_mon,12,0)
           aice_nh_sum_ann = temp(5::12)
           delete(temp)

           delete(aice_nh_sum_djf&time)
           aice_nh_sum_djf!0 = "TIME"
           aice_nh_sum_djf&TIME = ispan(syear(ee),eyear(ee),1)
           aice_nh_sum_djf&TIME@units = "YYYY"
           aice_nh_sum_djf&TIME@long_name = "time"
           copy_VarMeta(aice_nh_sum_djf,aice_nh_sum_mam)
           copy_VarMeta(aice_nh_sum_djf,aice_nh_sum_jja)
           copy_VarMeta(aice_nh_sum_djf,aice_nh_sum_son)
           copy_VarMeta(aice_nh_sum_djf,aice_nh_sum_ann)
        end if 

        if (OPT_CLIMO.eq."Full") then
           aice_nh = rmMonAnnCycTLL(aice_nh)
        else
           check_custom_climo(names(ee),syear(ee),eyear(ee),CLIMO_SYEAR,CLIMO_EYEAR)
           temp_arr = aice_nh
           delete(temp_arr&time)
           temp_arr&time = cd_calendar(aice_nh&time,-1)
           if (CLIMO_SYEAR.lt.0) then
              climo = clmMonTLL(temp_arr({(eyear(ee)+CLIMO_SYEAR)*100+1:(eyear(ee)+CLIMO_EYEAR)*100+12},:,:))       
           else
              climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
           end if             
           delete(temp_arr)
           aice_nh   = calcMonAnomTLL(aice_nh,climo) 
           delete(climo)
        end if

        dimZ = dimsizes(aice_nh)
        dim_j = dimZ(1)
        dim_i = dimZ(2)
        tttt = dtrend_msg_n(ispan(0,dimsizes(aice_nh&time)-1,1),aice_nh,False,True,0)
        aice_nh_trends_mon = aice_nh(0,:,:)
        aice_nh_trends_mon = (/ onedtond(tttt@slope, (/dim_j,dim_i/) ) /)
        aice_nh_trends_mon = aice_nh_trends_mon*dimsizes(aice_nh&time)
        aice_nh_trends_mon@units = aice_nh@units+" "+nyr(ee)+"yr~S~-1~N~"
        delete(tttt)
     
        aice_nh_seas = runave_n_Wrap(aice_nh,3,0,0)
        aice_nh_seas(0,:,:) = (/ dim_avg_n(aice_nh(:1,:,:),0) /)
        aice_nh_seas(dimsizes(aice_nh&time)-1,:,:) = (/ dim_avg_n(aice_nh(dimsizes(aice_nh&time)-2:,:,:),0) /)
        aice_nh_ann = runave_n_Wrap(aice_nh,12,0,0)
        delete(aice_nh)
     
        aice_nh_trends_seas = aice_nh_seas(:3,:,:)
        aice_nh_trends_seas = aice_nh_trends_seas@_FillValue
        aice_nh_trends_ann  = aice_nh_trends_seas(0,:,:)
        do ff = 0,4
           if (ff.le.3) then
              tarr = aice_nh_seas(ff*3::12,:,:)     
           end if  
           if (ff.eq.4) then
              tarr = aice_nh_ann(5::12,:,:)
           end if
           tttt = dtrend_msg_n(ispan(0,dimsizes(tarr&time)-1,1),tarr,False,True,0)   
           if (ff.le.3) then
              aice_nh_trends_seas(ff,:,:) = (/ onedtond(tttt@slope, (/dim_j,dim_i/) ) /)
           end if
           if (ff.eq.4) then
              aice_nh_trends_ann = (/ onedtond(tttt@slope, (/dim_j,dim_i/) ) /)
           end if
           delete([/tarr,tttt/])        
        end do
        aice_nh_trends_seas = aice_nh_trends_seas*nyr(ee)
        aice_nh_trends_seas@units = aice_nh_seas@units+" "+nyr(ee)+"yr~S~-1~N~"
        aice_nh_trends_ann = aice_nh_trends_ann*nyr(ee)
        aice_nh_trends_ann@units = aice_nh_ann@units+" "+nyr(ee)+"yr~S~-1~N~"         
        delete([/aice_nh_seas,aice_nh_ann,dim_j,dim_i,dimZ/])
     end if

     aice_sh_flag = 0
     aice_sh = data_read_in_ice(paths_sh(ee),"aice_sh",syear_sh(ee),eyear_sh(ee))    ; read in data, orient lats/lons correctly, set time coordinate variable up
     if (isatt(aice_sh,"is_all_missing")) then
        delete(aice_sh)
        aice_sh_flag = 1
     end if  
     if (aice_sh_flag.eq.0) then
        if (isatt(aice_sh,"area")) then
           area3d = aice_sh@area
           area3d_c = conform(aice_sh,area3d,(/1,2/))
           aice_sh_sum = aice_sh
;           aice_sh_sum = (/ (where((aice_sh/100.).ge.0.15,(aice_sh/100.),aice_sh@_FillValue))*area3d_c /) ;   ice area calculation (all cells > 15% kept) 
           aice_sh_sum = (/ where(aice_sh.ge.15,1.,aice_sh@_FillValue) /) ; ice extent calculation (all cells greater than 15% treated as 100% covered)
           aice_sh_sum = aice_sh_sum*area3d_c
           wgts = aice_sh_sum(0,:,:)
           wgts = 1.
           aice_sh_sum_mon = aice_sh_sum(:,0,0)   ; preallocate array to retain metadata
           aice_sh_sum_mon = (/ wgt_areasum2(aice_sh_sum,wgts,0) /)
;           do gg = 0,dimsizes(aice_sh&time)-1
;              aice_sh_sum_mon(gg) = (/ sum(aice_sh_sum(gg,:,:)) /)
;           end do
           aice_sh_sum_mon = aice_sh_sum_mon/1.e12
           aice_sh_sum_mon@units = "10^12 m2"
           aice_sh_sum_mon@long_name = "sea_ice_extent"
           if (isatt(aice_sh_sum_mon,"coordinates")) then
              delete(aice_sh_sum_mon@coordinates)
           end if
           delete([/aice_sh_sum,area3d,area3d_c,wgts/])

           taice = new((/dimsizes(aice_sh_sum_mon),1,1/),typeof(aice_sh_sum_mon))
           taice!0 = "time"
           taice!1 = "lat"
           taice!2 = "lon"
           taice(:,0,0) = (/ aice_sh_sum_mon /) ; convert to 3D array so we can use clmMonTLL and calcMonAnomTLL
           if (OPT_CLIMO.eq."Full") then
              taice = rmMonAnnCycTLL(taice)
           else
              check_custom_climo(names_sh(ee),syear_sh(ee),eyear_sh(ee),CLIMO_SYEAR,CLIMO_EYEAR)
              temp_arr = taice
              temp_arr&time = cd_calendar(aice_sh_sum_mon&time,-1)
              if (CLIMO_SYEAR.lt.0) then
                 climo = clmMonTLL(temp_arr({(eyear(ee)+CLIMO_SYEAR)*100+1:(eyear(ee)+CLIMO_EYEAR)*100+12},:,:))       
              else
                 climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
              end if               
              delete(temp_arr)
              taice   = calcMonAnomTLL(taice,climo) 
              delete(climo)
           end if
           aice_sh_sum_mon_anom = aice_sh_sum_mon
           aice_sh_sum_mon_anom = (/ taice(:,0,0) /)
           aice_sh_sum_mon_anom@long_name = "sea_ice_extent_anomaly"
           delete(taice)

           aice_sh_sum_feb = aice_sh_sum_mon(1::12)
           delete(aice_sh_sum_feb&time)
           aice_sh_sum_feb@long_name = "February sea_ice_extent"
           aice_sh_sum_feb!0 = "TIME"
           aice_sh_sum_feb&TIME = ispan(syear_sh(ee),eyear_sh(ee),1)
           aice_sh_sum_feb&TIME@units = "YYYY"
           aice_sh_sum_feb&TIME@long_name = "time"
           aice_sh_sum_mar = aice_sh_sum_mon(2::12)
           aice_sh_sum_mar@long_name = "March sea_ice_extent"
           copy_VarCoords(aice_sh_sum_feb,aice_sh_sum_mar)
           aice_sh_sum_sep = aice_sh_sum_mon(8::12)
           aice_sh_sum_sep@long_name = "September sea_ice_extent"
           copy_VarCoords(aice_sh_sum_feb,aice_sh_sum_sep)

           aice_sh_sum_climo = new(12,typeof(aice_sh_sum_mon))
           copy_VarAtts(aice_sh_sum_mon,aice_sh_sum_climo)
           aice_sh_sum_climo@long_name = "climatological_sea_ice_extent"
           aice_sh_sum_climo!0 = "time_mon2"
           aice_sh_sum_climo&time_mon2 = time_mon2
           do gg = 0,11
              aice_sh_sum_climo(gg) = (/ avg(aice_sh_sum_mon(gg::12)) /)
           end do    

           temp = runave_Wrap(aice_sh_sum_mon,3,0)
           aice_sh_sum_djf = temp(0::12)
           aice_sh_sum_mam = temp(3::12)
           aice_sh_sum_jja = temp(6::12)
           aice_sh_sum_son = temp(9::12)
           delete(temp)
           temp = runave_Wrap(aice_sh_sum_mon,12,0)
           aice_sh_sum_ann = temp(5::12)
           delete(temp)

           delete(aice_sh_sum_djf&time)
           aice_sh_sum_djf!0 = "TIME"
           aice_sh_sum_djf&TIME = ispan(syear_sh(ee),eyear_sh(ee),1)
           aice_sh_sum_djf&TIME@units = "YYYY"
           aice_sh_sum_djf&TIME@long_name = "time"
           copy_VarMeta(aice_sh_sum_djf,aice_sh_sum_mam)
           copy_VarMeta(aice_sh_sum_djf,aice_sh_sum_jja)
           copy_VarMeta(aice_sh_sum_djf,aice_sh_sum_son)
           copy_VarMeta(aice_sh_sum_djf,aice_sh_sum_ann)
        end if 

        if (OPT_CLIMO.eq."Full") then
           aice_sh = rmMonAnnCycTLL(aice_sh)
        else
           check_custom_climo(names_sh(ee),syear_sh(ee),eyear_sh(ee),CLIMO_SYEAR,CLIMO_EYEAR)
           temp_arr = aice_sh
           delete(temp_arr&time)
           temp_arr&time = cd_calendar(aice_sh&time,-1)
           if (CLIMO_SYEAR.lt.0) then
              climo = clmMonTLL(temp_arr({(eyear(ee)+CLIMO_SYEAR)*100+1:(eyear(ee)+CLIMO_EYEAR)*100+12},:,:))       
           else
              climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
           end if               
           delete(temp_arr)
           aice_sh   = calcMonAnomTLL(aice_sh,climo) 
           delete(climo)
        end if

        dimZ = dimsizes(aice_sh)
        dim_j = dimZ(1)
        dim_i = dimZ(2)
        tttt = dtrend_msg_n(ispan(0,dimsizes(aice_sh&time)-1,1),aice_sh,False,True,0)
        aice_sh_trends_mon = aice_sh(0,:,:)
        aice_sh_trends_mon = (/ onedtond(tttt@slope, (/dim_j,dim_i/) ) /)
        aice_sh_trends_mon = aice_sh_trends_mon*dimsizes(aice_sh&time)
        aice_sh_trends_mon@units = aice_sh@units+" "+nyr_sh(ee)+"yr~S~-1~N~"
        delete(tttt)
     
        aice_sh_seas = runave_n_Wrap(aice_sh,3,0,0)
        aice_sh_seas(0,:,:) = (/ dim_avg_n(aice_sh(:1,:,:),0) /)
        aice_sh_seas(dimsizes(aice_sh&time)-1,:,:) = (/ dim_avg_n(aice_sh(dimsizes(aice_sh&time)-2:,:,:),0) /)
        aice_sh_ann = runave_n_Wrap(aice_sh,12,0,0)
        delete(aice_sh)
     
        aice_sh_trends_seas = aice_sh_seas(:3,:,:)
        aice_sh_trends_seas = aice_sh_trends_seas@_FillValue
        aice_sh_trends_ann  = aice_sh_trends_seas(0,:,:)
        do ff = 0,4
           if (ff.le.3) then
              tarr = aice_sh_seas(ff*3::12,:,:)     
           end if  
           if (ff.eq.4) then
              tarr = aice_sh_ann(5::12,:,:)
           end if
           tttt = dtrend_msg_n(ispan(0,dimsizes(tarr&time)-1,1),tarr,False,True,0)   
           if (ff.le.3) then
              aice_sh_trends_seas(ff,:,:) = (/ onedtond(tttt@slope, (/dim_j,dim_i/) ) /)
           end if
           if (ff.eq.4) then
              aice_sh_trends_ann = (/ onedtond(tttt@slope, (/dim_j,dim_i/) ) /)
           end if
           delete([/tarr,tttt/])        
        end do
        aice_sh_trends_seas = aice_sh_trends_seas*nyr_sh(ee)
        aice_sh_trends_seas@units = aice_sh_seas@units+" "+nyr_sh(ee)+"yr~S~-1~N~"
        aice_sh_trends_ann = aice_sh_trends_ann*nyr_sh(ee)
        aice_sh_trends_ann@units = aice_sh_ann@units+" "+nyr_sh(ee)+"yr~S~-1~N~"         
        delete([/aice_sh_seas,aice_sh_ann,dim_j,dim_i,dimZ/])
     end if

     if (OUTPUT_DATA.eq."True".and.aice_nh_flag.eq.0) then
        modname = str_sub_str(names(ee)," ","_")
        bc = (/"/","'","(",")"/)
        do gg = 0,dimsizes(bc)-1
           modname = str_sub_str(modname,bc(gg),"_")
        end do
        fn = getenv("OUTDIR")+modname+".cvdp_data.aice.trends_timeseries.nh."+syear(ee)+"-"+eyear(ee)+".nc"
        if (.not.isfilepresent2(fn)) then
           z = addfile(fn,"c")
           z@source = "NCAR Climate Analysis Section's Climate Variability Diagnostics Package v"+getenv("VERSION")
           z@notes = "Data from "+names(ee)+" from "+syear(ee)+"-"+eyear(ee)
           if (OPT_CLIMO.eq."Full") then
              z@climatology = syear(ee)+"-"+eyear(ee)+" climatology removed prior to all calculations (other than means)"
           else
              if (CLIMO_SYEAR.lt.0) then
                 z@climatology = (eyear(ee)+CLIMO_SYEAR)+"-"+(eyear(ee)+CLIMO_EYEAR)+" climatology removed prior to all calculations (other than means)"
              else
                 z@climatology = CLIMO_SYEAR+"-"+CLIMO_EYEAR+" climatology removed prior to all calculations (other than means)"
              end if
           end if
           z@Conventions = "CF-1.6"
        else
           z = addfile(fn,"w")
        end if
        nh_trends_seas = aice_nh_trends_seas
        if (isatt(nh_trends_seas,"lat2d")) then    ; if there is a lat2d there will be a lon2d
           LAT2D = nh_trends_seas@lat2d
           LON2D = nh_trends_seas@lon2d
           delete(nh_trends_seas@lat2d)
           delete(nh_trends_seas@lon2d)
           copy_VarCoords(nh_trends_seas(0,:,:),LAT2D)
           copy_VarCoords(nh_trends_seas(0,:,:),LON2D)
           z->lat2d_ice_nh = set_varAtts(LAT2D,"Northern Hemisphere ice grid 2-dimensional latitudes","","")
           z->lon2d_ice_nh = set_varAtts(LON2D,"Northern Hemisphere ice grid 2-dimensional longitudes","","")
           delete([/LAT2D,LON2D/])
           nh_trends_seas@coordinates ="lat2d_ice_nh lon2d_ice_nh"
        end if           
        if (isatt(nh_trends_seas,"area")) then
           delete(nh_trends_seas@area)
        end if
        nh_trends_ann = (/ aice_nh_trends_ann /)
        copy_VarMeta(nh_trends_seas(0,:,:),nh_trends_ann)
        nh_trends_mon = (/ aice_nh_trends_mon /)
        copy_VarMeta(nh_trends_seas(0,:,:),nh_trends_mon)
        z->sic_nh_trends_djf = set_varAtts(nh_trends_seas(0,:,:),"Northern Hemisphere sic trends (DJF)","","")
        z->sic_nh_trends_mam = set_varAtts(nh_trends_seas(1,:,:),"Northern Hemisphere sic trends (MAM)","","")
        z->sic_nh_trends_jja = set_varAtts(nh_trends_seas(2,:,:),"Northern Hemisphere sic trends (JJA)","","")
        z->sic_nh_trends_son = set_varAtts(nh_trends_seas(3,:,:),"Northern Hemisphere sic trends (SON)","","")
        z->sic_nh_trends_ann = set_varAtts(nh_trends_ann,"Northern Hemisphere sic trends (annual)","","")
        z->sic_nh_trends_mon = set_varAtts(nh_trends_mon,"Northern Hemisphere sic trends (monthly)","","")
        
        if (isvar("aice_nh_sum_djf")) then
           nh_sum_climo = (/ aice_nh_sum_climo /)
           copy_VarAtts(aice_nh_sum_djf,nh_sum_climo)
           nh_sum_climo!0 = "time_mon2"
           nh_sum_climo&time_mon2 = time_mon2
           if (isatt(nh_sum_climo,"lat2d")) then
              delete(nh_sum_climo@lat2d)
              delete(nh_sum_climo@lon2d)
           end if
           if (isatt(nh_sum_climo,"area")) then
              delete(nh_sum_climo@area)
           end if
           z->sic_nh_extent_climo = set_varAtts(nh_sum_climo,"Northern Hemisphere sic extent climatology","","")
           nh_sum_djf = aice_nh_sum_djf
           if (isatt(nh_sum_djf,"lat2d")) then
              delete(nh_sum_djf@lat2d)
              delete(nh_sum_djf@lon2d)
           end if
           if (isatt(nh_sum_djf,"area")) then
              delete(nh_sum_djf@area)
           end if
           nh_sum_mam = (/ aice_nh_sum_mam /)
           copy_VarMeta(nh_sum_djf,nh_sum_mam)
           nh_sum_jja = (/ aice_nh_sum_jja /)
           copy_VarMeta(nh_sum_djf,nh_sum_jja)
           nh_sum_son = (/ aice_nh_sum_son /)
           copy_VarMeta(nh_sum_djf,nh_sum_son)
           nh_sum_ann = (/ aice_nh_sum_ann /)
           copy_VarMeta(nh_sum_djf,nh_sum_ann)
           nh_sum_feb = (/ aice_nh_sum_feb /)
           copy_VarMeta(nh_sum_djf,nh_sum_feb)
           nh_sum_mar = (/ aice_nh_sum_mar /)
           copy_VarMeta(nh_sum_djf,nh_sum_mar)
           nh_sum_sep = (/ aice_nh_sum_sep /)
           copy_VarMeta(nh_sum_djf,nh_sum_sep)
           nh_sum_mon = aice_nh_sum_mon
           nh_sum_mon_anom = aice_nh_sum_mon_anom
           if (isatt(nh_sum_mon,"lat2d")) then
              delete(nh_sum_mon@lat2d)
              delete(nh_sum_mon@lon2d)
              delete(nh_sum_mon_anom@lat2d)
              delete(nh_sum_mon_anom@lon2d)
           end if
           if (isatt(nh_sum_mon,"area")) then
              delete(nh_sum_mon@area)
              delete(nh_sum_mon_anom@area)
           end if
           z->sic_nh_extent_djf = set_varAtts(nh_sum_djf,"Northern Hemisphere sic extent timeseries (DJF)","","")
           z->sic_nh_extent_mam = set_varAtts(nh_sum_mam,"Northern Hemisphere sic extent timeseries (MAM)","","")
           z->sic_nh_extent_jja = set_varAtts(nh_sum_jja,"Northern Hemisphere sic extent timeseries (JJA)","","")
           z->sic_nh_extent_son = set_varAtts(nh_sum_son,"Northern Hemisphere sic extent timeseries (SON)","","")
           z->sic_nh_extent_ann = set_varAtts(nh_sum_ann,"Northern Hemisphere sic extent timeseries (annual)","","")
           z->sic_nh_extent_mon = set_varAtts(nh_sum_mon,"Northern Hemisphere sic extent timeseries (monthly)","","")
           z->sic_nh_extent_mon_anom = set_varAtts(nh_sum_mon_anom,"Northern Hemisphere sic extent anomaly timeseries (monthly)","","")
           z->sic_nh_extent_feb = set_varAtts(nh_sum_feb,"Northern Hemisphere sic extent timeseries (February)","","")
           z->sic_nh_extent_mar = set_varAtts(nh_sum_mar,"Northern Hemisphere sic extent timeseries (March)","","")
           z->sic_nh_extent_sep = set_varAtts(nh_sum_sep,"Northern Hemisphere sic extent timeseries (September)","","")
           delete([/nh_sum_djf,nh_sum_mam,nh_sum_jja,nh_sum_son,nh_sum_ann,nh_sum_feb,nh_sum_mar,nh_sum_sep,nh_sum_mon,nh_sum_mon_anom/])
        end if
        delete([/nh_trends_seas,nh_trends_ann,nh_trends_mon/])
        delete(z)
     end if
     if (OUTPUT_DATA.eq."True".and.aice_sh_flag.eq.0) then
        modname = str_sub_str(names_sh(ee)," ","_")
        bc = (/"/","'","(",")"/)
        do gg = 0,dimsizes(bc)-1
           modname = str_sub_str(modname,bc(gg),"_")
        end do
        fn = getenv("OUTDIR")+modname+".cvdp_data.aice.trends_timeseries.sh."+syear_sh(ee)+"-"+eyear_sh(ee)+".nc"
        if (.not.isfilepresent2(fn)) then
           z = addfile(fn,"c")
           z@source = "NCAR Climate Analysis Section's Climate Variability Diagnostics Package v"+getenv("VERSION")
           z@notes = "Data from "+names_sh(ee)+" from "+syear_sh(ee)+"-"+eyear_sh(ee)
           if (OPT_CLIMO.eq."Full") then
              z@climatology = syear_sh(ee)+"-"+eyear_sh(ee)+" climatology removed prior to all calculations (other than means)"
           else
              if (CLIMO_SYEAR.lt.0) then
                 z@climatology = (eyear(ee)+CLIMO_SYEAR)+"-"+(eyear(ee)+CLIMO_EYEAR)+" climatology removed prior to all calculations (other than means)"
              else
                 z@climatology = CLIMO_SYEAR+"-"+CLIMO_EYEAR+" climatology removed prior to all calculations (other than means)"
              end if
           end if
           z@Conventions = "CF-1.6"
        else
           z = addfile(fn,"w")
        end if
        sh_trends_seas = aice_sh_trends_seas
        if (isatt(sh_trends_seas,"lat2d")) then    ; if there is a lat2d there will be a lon2d
           LAT2D = sh_trends_seas@lat2d
           LON2D = sh_trends_seas@lon2d
           delete(sh_trends_seas@lat2d)
           delete(sh_trends_seas@lon2d)
           copy_VarCoords(sh_trends_seas(0,:,:),LAT2D)
           copy_VarCoords(sh_trends_seas(0,:,:),LON2D)
           z->lat2d_ice_sh = set_varAtts(LAT2D,"Northern Hemisphere ice grid 2-dimensional latitudes","","")
           z->lon2d_ice_sh = set_varAtts(LON2D,"Northern Hemisphere ice grid 2-dimensional longitudes","","")
           delete([/LAT2D,LON2D/])
           sh_trends_seas@coordinates ="lat2d_ice_sh lon2d_ice_sh"
        end if  
        if (isatt(sh_trends_seas,"area")) then
           delete(sh_trends_seas@area)
        end if
        sh_trends_seas!1 = "j2"
        sh_trends_seas!2 = "i2"
        sh_trends_seas@long_name = sh_trends_seas@long_name+" trends"
        sh_trends_ann = (/ aice_sh_trends_ann /)
        copy_VarMeta(sh_trends_seas(0,:,:),sh_trends_ann)
        sh_trends_mon = (/ aice_sh_trends_mon /)
        copy_VarMeta(sh_trends_seas(0,:,:),sh_trends_mon)
        z->sic_sh_trends_djf = set_varAtts(sh_trends_seas(0,:,:),"Southern Hemisphere sic trends (DJF)","","")
        z->sic_sh_trends_mam = set_varAtts(sh_trends_seas(1,:,:),"Southern Hemisphere sic trends (MAM)","","")
        z->sic_sh_trends_jja = set_varAtts(sh_trends_seas(2,:,:),"Southern Hemisphere sic trends (JJA)","","")
        z->sic_sh_trends_son = set_varAtts(sh_trends_seas(3,:,:),"Southern Hemisphere sic trends (SON)","","")
        z->sic_sh_trends_ann = set_varAtts(sh_trends_ann,"Southern Hemisphere sic trends (annual)","","")
        z->sic_sh_trends_mon = set_varAtts(sh_trends_mon,"Southern Hemisphere sic trends (monthly)","","")
        
        if (isvar("aice_sh_sum_djf")) then
           sh_sum_climo = (/ aice_sh_sum_climo /)
           copy_VarAtts(aice_sh_sum_djf,sh_sum_climo)
           sh_sum_climo!0 = "time_mon2"
           sh_sum_climo&time_mon2 = time_mon2
           if (isatt(sh_sum_climo,"lat2d")) then
              delete(sh_sum_climo@lat2d)
              delete(sh_sum_climo@lon2d)
           end if
           if (isatt(sh_sum_climo,"area")) then
              delete(sh_sum_climo@area)
           end if
           z->sic_sh_extent_climo = set_varAtts(sh_sum_climo,"Southern Hemisphere sic extent climatology","","")
           sh_sum_djf = aice_sh_sum_djf
           if (isatt(sh_sum_djf,"lat2d")) then
              delete(sh_sum_djf@lat2d)
              delete(sh_sum_djf@lon2d)
           end if
           if (isatt(sh_sum_djf,"area")) then
              delete(sh_sum_djf@area)
           end if
           sh_sum_mam = (/ aice_sh_sum_mam /)
           copy_VarMeta(sh_sum_djf,sh_sum_mam)
           sh_sum_jja = (/ aice_sh_sum_jja /)
           copy_VarMeta(sh_sum_djf,sh_sum_jja)
           sh_sum_son = (/ aice_sh_sum_son /)
           copy_VarMeta(sh_sum_djf,sh_sum_son)
           sh_sum_ann = (/ aice_sh_sum_ann /)
           copy_VarMeta(sh_sum_djf,sh_sum_ann)
           sh_sum_feb = (/ aice_sh_sum_feb /)
           copy_VarMeta(sh_sum_djf,sh_sum_feb)
           sh_sum_mar = (/ aice_sh_sum_mar /)
           copy_VarMeta(sh_sum_djf,sh_sum_mar)
           sh_sum_sep = (/ aice_sh_sum_sep /)
           copy_VarMeta(sh_sum_djf,sh_sum_sep)
           sh_sum_mon = aice_sh_sum_mon
           sh_sum_mon_anom = aice_sh_sum_mon_anom
           if (isatt(sh_sum_mon,"lat2d")) then
              delete(sh_sum_mon@lat2d)
              delete(sh_sum_mon@lon2d)
              delete(sh_sum_mon_anom@lat2d)
              delete(sh_sum_mon_anom@lon2d)
           end if
           if (isatt(sh_sum_mon,"area")) then
              delete(sh_sum_mon@area)
              delete(sh_sum_mon_anom@area)
           end if
           z->sic_sh_extent_djf = set_varAtts(sh_sum_djf,"Southern Hemisphere sic extent timeseries (DJF)","","")
           z->sic_sh_extent_mam = set_varAtts(sh_sum_mam,"Southern Hemisphere sic extent timeseries (MAM)","","")
           z->sic_sh_extent_jja = set_varAtts(sh_sum_jja,"Southern Hemisphere sic extent timeseries (JJA)","","")
           z->sic_sh_extent_son = set_varAtts(sh_sum_son,"Southern Hemisphere sic extent timeseries (SON)","","")
           z->sic_sh_extent_ann = set_varAtts(sh_sum_ann,"Southern Hemisphere sic extent timeseries (annual)","","")
           z->sic_sh_extent_mon = set_varAtts(sh_sum_mon,"Southern Hemisphere sic extent timeseries (monthly)","","")
           z->sic_sh_extent_mon_anom = set_varAtts(sh_sum_mon_anom,"Southern Hemisphere sic extent anomaly timeseries (monthly)","","")
           z->sic_sh_extent_feb = set_varAtts(sh_sum_feb,"Southern Hemisphere sic extent timeseries (February)","","")
           z->sic_sh_extent_mar = set_varAtts(sh_sum_mar,"Southern Hemisphere sic extent timeseries (March)","","")
           z->sic_sh_extent_sep = set_varAtts(sh_sum_sep,"Southern Hemisphere sic extent timeseries (September)","","")
           delete([/sh_sum_djf,sh_sum_mam,sh_sum_jja,sh_sum_son,sh_sum_ann,sh_sum_feb,sh_sum_mar,sh_sum_sep,sh_sum_mon,sh_sum_mon_anom/])
        end if
        delete([/sh_trends_seas,sh_trends_ann,sh_trends_mon/])
        delete(z)
     end if
     if (aice_nh_flag.eq.0) then
        aice_nh_trends_seas = where(aice_nh_trends_seas.eq.0,aice_nh_trends_seas@_FillValue,aice_nh_trends_seas)
        aice_nh_trends_ann = where(aice_nh_trends_ann.eq.0,aice_nh_trends_ann@_FillValue,aice_nh_trends_ann)
        aice_nh_trends_mon = where(aice_nh_trends_mon.eq.0,aice_nh_trends_mon@_FillValue,aice_nh_trends_mon)
     end if
     if (aice_sh_flag.eq.0) then
        aice_sh_trends_seas = where(aice_sh_trends_seas.eq.0,aice_sh_trends_seas@_FillValue,aice_sh_trends_seas)
        aice_sh_trends_ann = where(aice_sh_trends_ann.eq.0,aice_sh_trends_ann@_FillValue,aice_sh_trends_ann)
        aice_sh_trends_mon = where(aice_sh_trends_mon.eq.0,aice_sh_trends_mon@_FillValue,aice_sh_trends_mon)
     end if
;==========================================================================================
     res = True
     res@gsnPolar = "NH"  
     res@mpMinLatF    = 40.
     res@mpCenterLonF = 0.
     res@mpGeophysicalLineColor = "gray42"
     if (wks_type.eq."png") then
        res@mpGeophysicalLineThicknessF = 2.  
     else
        res@mpGeophysicalLineThicknessF = 1.  
     end if     
     res@mpGridAndLimbOn = False
     res@mpLandFillColor  = "gray75"           
     res@mpFillDrawOrder  = "PostDraw"       
     res@mpPerimDrawOrder  = "PostDraw"     

     res@mpOutlineOn = True  
     res@gsnDraw      = False
     res@gsnFrame     = False
     res@gsnAddCyclic = True
     res@cnLineLabelsOn = False
     res@cnFillOn        = True
     res@cnLinesOn       = False
     res@trGridType = "TriangularMesh"
;     res@cnFillMode = "RasterFill"
     res@lbLabelBarOn    = False
  
     res@cnLevelSelectionMode = "ExplicitLevels"
     res@cnLevels = (/-90,-70,-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30,35,40,45,50,70,90/)
     cmap = gsn_retrieve_colormap(wks_trends_djf)
     dimc = dimsizes(cmap)
     cmap2 = cmap(::-1,:)
     res@cnFillPalette = cmap2(:dimc(0)-3,:)
     delete([/cmap,cmap2,dimc/])

     res@gsnLeftStringOrthogonalPosF = -0.03
     res@gsnLeftStringParallelPosF = .005
     res@gsnRightStringOrthogonalPosF = -0.03
     res@gsnRightStringParallelPosF = 0.96
     res@gsnRightString = ""
     res@gsnLeftString = ""
     if (nsim.le.5) then
        res@gsnLeftStringFontHeightF = 0.018
        res@gsnCenterStringFontHeightF = 0.022
        res@gsnRightStringFontHeightF = 0.018
     else
        res@gsnLeftStringFontHeightF = 0.024
        res@gsnCenterStringFontHeightF = 0.028
        res@gsnRightStringFontHeightF = 0.024     
     end if

     xyres = True
     xyres@gsnDraw = False
     xyres@gsnFrame = False
     xyres@gsnYRefLine = 0.0
     xyres@gsnYRefLineColor = "gray18"     
     if (wks_type.eq."png") then
        xyres@xyLineThicknessF = 4.
     else
        xyres@xyLineThicknessF = 2.
     end if
     if (isfilepresent2("obs_aice_nh").and.ee.eq.0) then
        xyres@xyLineColor = "black"
     else
        xyres@xyLineColor = "royalblue"
     end if
     xyres@tiYAxisString = ""
     if (nsim.le.5) then
        xyres@tmXBLabelFontHeightF = 0.0125
        xyres@tmYLLabelFontHeightF = 0.0125
        xyres@gsnLeftStringFontHeightF = 0.0155 
        xyres@gsnCenterStringFontHeightF = 0.0155       
        xyres@gsnRightStringFontHeightF = 0.012     
     else
        xyres@tmXBLabelFontHeightF = 0.018
        xyres@tmYLLabelFontHeightF = 0.018
        xyres@gsnLeftStringFontHeightF = 0.024
        xyres@gsnCenterStringFontHeightF = 0.024   
        xyres@gsnRightStringFontHeightF = 0.020     
     end if
     xyres@gsnLeftStringOrthogonalPosF = 0.025
     xyres@gsnRightStringOrthogonalPosF = xyres@gsnLeftStringOrthogonalPosF
     xyres@vpXF = 0.05
     xyres@vpHeightF = 0.15
     if (SCALE_TIMESERIES.eq."True") then
        xyres@vpWidthF = 0.9*((nyr(ee)*1.)/nyr_max)
     else
        xyres@vpWidthF = 0.9
     end if
     xyres@xyMonoDashPattern = True
     xyres@gsnLeftString = ""     
     xyres@gsnCenterString = ""
     xyres@gsnRightString = ""
     
     xyres_c = xyres
     xyres_c@trXMinF = 1
     xyres_c@trXMaxF = 12
     xyres_c@vpWidthF = 0.65
     xyres_c@vpHeightF = 0.45
     if (isfilepresent2("obs_aice_nh").and.ee.eq.0) then
        xyres_c@xyLineColor = "black"
     else
        xyres_c@xyLineColors = (/"royalblue","gray60"/)
     end if
     xyres_c@tmXBMode    = "Explicit"        ; explicit labels
     xyres_c@tmXBValues     = ispan(1,12,1)
     xyres_c@tmXBLabels = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
     xyres_c@tmXTOn = False
     xyres_c@gsnLeftStringOrthogonalPosF = 0.025
     xyres_c@gsnCenterStringOrthogonalPosF = xyres_c@gsnLeftStringOrthogonalPosF
     xyres_c@gsnRightStringOrthogonalPosF = xyres_c@gsnLeftStringOrthogonalPosF

     xyres@trXMinF = syear(ee)-.5
     xyres@trXMaxF = eyear(ee)+0.5

     xyres2 = xyres
     xyres2@xyLineColor = "gray60"
     xyres2@xyCurveDrawOrder = "PreDraw"

     if (aice_nh_flag.eq.0) then 
        res@gsnLeftString = syear(ee)+"-"+eyear(ee)
        res@gsnRightString = aice_nh_trends_seas@units
        res@gsnCenterString = names(ee)

        plot_trends_nh_djf(ee) = gsn_csm_contour_map(wks_trends_djf,aice_nh_trends_seas(0,:,:),res)
        plot_trends_nh_mam(ee) = gsn_csm_contour_map(wks_trends_mam,aice_nh_trends_seas(1,:,:),res)
        plot_trends_nh_jja(ee) = gsn_csm_contour_map(wks_trends_jja,aice_nh_trends_seas(2,:,:),res)
        plot_trends_nh_son(ee) = gsn_csm_contour_map(wks_trends_son,aice_nh_trends_seas(3,:,:),res)
        plot_trends_nh_ann(ee) = gsn_csm_contour_map(wks_trends_ann,aice_nh_trends_ann,res)
        plot_trends_nh_mon(ee) = gsn_csm_contour_map(wks_trends_mon,aice_nh_trends_mon,res)
        delete([/aice_nh_trends_seas,aice_nh_trends_ann,aice_nh_trends_mon/])

        if (isvar("aice_nh_sum_djf")) then
           xyres@gsnLeftString = names(ee)+ph_s  
           if (isfilepresent2("obs_aice_nh").and.ee.eq.0) then
              aice_nh_sum_djf_obs_min = min(aice_nh_sum_djf)
              aice_nh_sum_djf_obs_max = max(aice_nh_sum_djf)
              aice_nh_sum_mam_obs_min = min(aice_nh_sum_mam)
              aice_nh_sum_mam_obs_max = max(aice_nh_sum_mam)
              aice_nh_sum_jja_obs_min = min(aice_nh_sum_jja)
              aice_nh_sum_jja_obs_max = max(aice_nh_sum_jja)
              aice_nh_sum_son_obs_min = min(aice_nh_sum_son)
              aice_nh_sum_son_obs_max = max(aice_nh_sum_son)
              aice_nh_sum_ann_obs_min = min(aice_nh_sum_ann)
              aice_nh_sum_ann_obs_max = max(aice_nh_sum_ann)
              aice_nh_sum_mon_obs_min = min(aice_nh_sum_mon)
              aice_nh_sum_mon_obs_max = max(aice_nh_sum_mon)
              aice_nh_sum_mon_anom_obs_min = min(aice_nh_sum_mon_anom)
              aice_nh_sum_mon_anom_obs_max = max(aice_nh_sum_mon_anom)
              aice_nh_sum_feb_obs_min = min(aice_nh_sum_feb)
              aice_nh_sum_feb_obs_max = max(aice_nh_sum_feb)
              aice_nh_sum_mar_obs_min = min(aice_nh_sum_mar)
              aice_nh_sum_mar_obs_max = max(aice_nh_sum_mar)
              aice_nh_sum_sep_obs_min = min(aice_nh_sum_sep)
              aice_nh_sum_sep_obs_max = max(aice_nh_sum_sep)
              aice_nh_obs_djf = aice_nh_sum_djf
              aice_nh_obs_mam = aice_nh_sum_mam
              aice_nh_obs_jja = aice_nh_sum_jja
              aice_nh_obs_son = aice_nh_sum_son
              aice_nh_obs_ann = aice_nh_sum_ann
              aice_nh_obs_mon = aice_nh_sum_mon
              aice_nh_obs_mon_anom = aice_nh_sum_mon_anom
              aice_nh_obs_feb = aice_nh_sum_feb
              aice_nh_obs_mar = aice_nh_sum_mar
              aice_nh_obs_sep = aice_nh_sum_sep
              aice_nh_obs_climo = aice_nh_sum_climo
           end if

           xyres_c@gsnLeftString = syear(ee)+"-"+eyear(ee)  
           xyres_c@gsnCenterString = names(ee)+ph_s   ; ph_s is used to denote extent timeseries that have had pole hole areas added in   
           xyres_c@gsnRightString = "("+str_sub_str(str_sub_str(aice_nh_sum_climo@units,"m2","m~S~2~N~"),"10^12","10~S~12~N~")+")"
           if (ee.eq.0) then
              plot_iceext_nh_climo(ee) = gsn_csm_xy(wks_iceext_mon,ispan(1,12,1),aice_nh_sum_climo,xyres_c)
           end if    
           if (ee.ge.1.) then
              if (isvar("aice_nh_obs_climo")) then
                 tarr = new((/2,12/),float)
                 tarr(0,:) = (/ tofloat(aice_nh_sum_climo) /)
                 tarr(1,:) = (/ tofloat(aice_nh_obs_climo) /)
                 plot_iceext_nh_climo(ee) = gsn_csm_xy(wks_iceext_mon,ispan(1,12,1),tarr,xyres_c)
                 delete(tarr)
              else
                 plot_iceext_nh_climo(ee) = gsn_csm_xy(wks_iceext_mon,ispan(1,12,1),aice_nh_sum_climo,xyres_c)
              end if
           end if
              
           if (ee.ge.1.and.isvar("aice_nh_sum_djf_obs_min")) then
              xyres@trYMinF = min((/min(aice_nh_sum_djf),aice_nh_sum_djf_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_nh_sum_djf),aice_nh_sum_djf_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_nh_sum_djf)-1,1),aice_nh_sum_djf,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr(ee),2,True)+aice_nh_sum_djf@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_nh_djf(ee)     = gsn_csm_xy(wks_iceext_djf,ispan(syear(ee),eyear(ee),1),aice_nh_sum_djf,xyres)
           if (ee.ge.1.and.isvar("aice_nh_obs_djf")) then
              plot_iceext_nh_obs_djf = gsn_csm_xy(wks_iceext_djf,ispan(syear(0),eyear(0),1),aice_nh_obs_djf,xyres2)
              overlay(plot_iceext_nh_djf(ee),plot_iceext_nh_obs_djf)
           end if

           if (ee.ge.1.and.isvar("aice_nh_sum_mam_obs_min")) then
              xyres@trYMinF = min((/min(aice_nh_sum_mam),aice_nh_sum_mam_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_nh_sum_mam),aice_nh_sum_mam_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_nh_sum_mam)-1,1),aice_nh_sum_mam,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr(ee),2,True)+aice_nh_sum_mam@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_nh_mam(ee)     = gsn_csm_xy(wks_iceext_mam,ispan(syear(ee),eyear(ee),1),aice_nh_sum_mam,xyres)
           if (ee.ge.1.and.isvar("aice_nh_obs_mam")) then
              plot_iceext_nh_obs_mam = gsn_csm_xy(wks_iceext_mam,ispan(syear(0),eyear(0),1),aice_nh_obs_mam,xyres2)
              overlay(plot_iceext_nh_mam(ee),plot_iceext_nh_obs_mam)
           end if      

           if (ee.ge.1.and.isvar("aice_nh_sum_jja_obs_min")) then
              xyres@trYMinF = min((/min(aice_nh_sum_jja),aice_nh_sum_jja_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_nh_sum_jja),aice_nh_sum_jja_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_nh_sum_jja)-1,1),aice_nh_sum_jja,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr(ee),2,True)+aice_nh_sum_jja@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_nh_jja(ee)     = gsn_csm_xy(wks_iceext_jja,ispan(syear(ee),eyear(ee),1),aice_nh_sum_jja,xyres)
           if (ee.ge.1.and.isvar("aice_nh_obs_jja")) then
              plot_iceext_nh_obs_jja = gsn_csm_xy(wks_iceext_jja,ispan(syear(0),eyear(0),1),aice_nh_obs_jja,xyres2)
              overlay(plot_iceext_nh_jja(ee),plot_iceext_nh_obs_jja)
           end if

           if (ee.ge.1.and.isvar("aice_nh_sum_son_obs_min")) then
              xyres@trYMinF = min((/min(aice_nh_sum_son),aice_nh_sum_son_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_nh_sum_son),aice_nh_sum_son_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_nh_sum_son)-1,1),aice_nh_sum_son,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr(ee),2,True)+aice_nh_sum_son@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_nh_son(ee)     = gsn_csm_xy(wks_iceext_son,ispan(syear(ee),eyear(ee),1),aice_nh_sum_son,xyres)
           if (ee.ge.1.and.isvar("aice_nh_obs_son")) then
              plot_iceext_nh_obs_son = gsn_csm_xy(wks_iceext_son,ispan(syear(0),eyear(0),1),aice_nh_obs_son,xyres2)
              overlay(plot_iceext_nh_son(ee),plot_iceext_nh_obs_son)
           end if

           if (ee.ge.1.and.isvar("aice_nh_sum_ann_obs_min")) then
              xyres@trYMinF = min((/min(aice_nh_sum_ann),aice_nh_sum_ann_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_nh_sum_ann),aice_nh_sum_ann_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_nh_sum_ann)-1,1),aice_nh_sum_ann,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr(ee),2,True)+aice_nh_sum_ann@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_nh_ann(ee)     = gsn_csm_xy(wks_iceext_ann,ispan(syear(ee),eyear(ee),1),aice_nh_sum_ann,xyres)
           if (ee.ge.1.and.isvar("aice_nh_obs_ann")) then
              plot_iceext_nh_obs_ann = gsn_csm_xy(wks_iceext_ann,ispan(syear(0),eyear(0),1),aice_nh_obs_ann,xyres2)
              overlay(plot_iceext_nh_ann(ee),plot_iceext_nh_obs_ann)
           end if
           delete(tttt)

           if (ee.ge.1.and.isvar("aice_nh_sum_feb_obs_min")) then
              xyres@trYMinF = min((/min(aice_nh_sum_feb),aice_nh_sum_feb_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_nh_sum_feb),aice_nh_sum_feb_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_nh_sum_feb)-1,1),aice_nh_sum_feb,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr(ee),2,True)+aice_nh_sum_feb@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_nh_feb(ee)     = gsn_csm_xy(wks_iceext_febmar,ispan(syear(ee),eyear(ee),1),aice_nh_sum_feb,xyres)
           if (ee.ge.1.and.isvar("aice_nh_obs_feb")) then
              plot_iceext_nh_obs_feb = gsn_csm_xy(wks_iceext_febmar,ispan(syear(0),eyear(0),1),aice_nh_obs_feb,xyres2)
              overlay(plot_iceext_nh_feb(ee),plot_iceext_nh_obs_feb)
           end if
           delete(tttt)

           if (ee.ge.1.and.isvar("aice_nh_sum_mar_obs_min")) then
              xyres@trYMinF = min((/min(aice_nh_sum_mar),aice_nh_sum_mar_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_nh_sum_mar),aice_nh_sum_mar_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_nh_sum_mar)-1,1),aice_nh_sum_mar,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr(ee),2,True)+aice_nh_sum_mar@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_nh_mar(ee)     = gsn_csm_xy(wks_iceext_febmar,ispan(syear(ee),eyear(ee),1),aice_nh_sum_mar,xyres)
           if (ee.ge.1.and.isvar("aice_nh_obs_mar")) then
              plot_iceext_nh_obs_mar = gsn_csm_xy(wks_iceext_febmar,ispan(syear(0),eyear(0),1),aice_nh_obs_mar,xyres2)
              overlay(plot_iceext_nh_mar(ee),plot_iceext_nh_obs_mar)
           end if
           delete(tttt)

           if (ee.ge.1.and.isvar("aice_nh_sum_sep_obs_min")) then
              xyres@trYMinF = min((/min(aice_nh_sum_sep),aice_nh_sum_sep_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_nh_sum_sep),aice_nh_sum_sep_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_nh_sum_sep)-1,1),aice_nh_sum_sep,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr(ee),2,True)+aice_nh_sum_sep@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_nh_sep(ee)     = gsn_csm_xy(wks_iceext_sep,ispan(syear(ee),eyear(ee),1),aice_nh_sum_sep,xyres)
           if (ee.ge.1.and.isvar("aice_nh_obs_sep")) then
              plot_iceext_nh_obs_sep = gsn_csm_xy(wks_iceext_sep,ispan(syear(0),eyear(0),1),aice_nh_obs_sep,xyres2)
              overlay(plot_iceext_nh_sep(ee),plot_iceext_nh_obs_sep)
           end if
           delete(tttt)

           if (ee.ge.1.and.isvar("aice_nh_sum_mon_obs_min")) then
              xyres@trYMinF = min((/min(aice_nh_sum_mon),aice_nh_sum_mon_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_nh_sum_mon),aice_nh_sum_mon_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_nh_sum_mon)-1,1),aice_nh_sum_mon,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr(ee),2,True)+aice_nh_sum_mon@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_nh_mon(ee)     = gsn_csm_xy(wks_iceext_mon,fspan(syear(ee),eyear(ee)+.91667,dimsizes(aice_nh_sum_mon)),aice_nh_sum_mon,xyres)
           if (ee.ge.1.and.isvar("aice_nh_obs_mon")) then
              plot_iceext_nh_obs_mon = gsn_csm_xy(wks_iceext_mon,fspan(syear(0),eyear(0)+.91667,dimsizes(aice_nh_obs_mon)),aice_nh_obs_mon,xyres2)
              overlay(plot_iceext_nh_mon(ee),plot_iceext_nh_obs_mon)
           end if
           delete(tttt)

           if (ee.ge.1.and.isvar("aice_nh_sum_mon_anom_obs_min")) then
              xyres@trYMinF = min((/min(aice_nh_sum_mon_anom),aice_nh_sum_mon_anom_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_nh_sum_mon_anom),aice_nh_sum_mon_anom_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_nh_sum_mon_anom)-1,1),aice_nh_sum_mon_anom,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr(ee),2,True)+aice_nh_sum_mon_anom@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_nh_mon_anom(ee)     = gsn_csm_xy(wks_iceext_mon,fspan(syear(ee),eyear(ee)+.91667,dimsizes(aice_nh_sum_mon_anom)),aice_nh_sum_mon_anom,xyres)
           if (ee.ge.1.and.isvar("aice_nh_obs_mon_anom")) then
              plot_iceext_nh_obs_mon_anom = gsn_csm_xy(wks_iceext_mon,fspan(syear(0),eyear(0)+.91667,dimsizes(aice_nh_obs_mon_anom)),aice_nh_obs_mon_anom,xyres2)
              overlay(plot_iceext_nh_mon_anom(ee),plot_iceext_nh_obs_mon_anom)
           end if
           delete([/aice_nh_sum_djf,aice_nh_sum_mam,aice_nh_sum_jja,aice_nh_sum_son,aice_nh_sum_ann,aice_nh_sum_feb,aice_nh_sum_mar,aice_nh_sum_sep,aice_nh_sum_mon,aice_nh_sum_mon_anom,tttt/])
        end if
     end if
     
     delete(res@mpMinLatF)
     res@mpMaxLatF = -45.
     res@gsnPolar = "SH"  
     if (SCALE_TIMESERIES.eq."True") then
        xyres@vpWidthF = 0.9*((nyr_sh(ee)*1.)/nyr_max_sh)
     else
        xyres@vpWidthF = 0.9
     end if     
     xyres2@vpWidthF = xyres@vpWidthF
     xyres@trXMinF = syear_sh(ee)-.5
     xyres@trXMaxF = eyear_sh(ee)+0.5
     if (aice_sh_flag.eq.0) then 
        res@gsnLeftString = syear_sh(ee)+"-"+eyear_sh(ee)
        res@gsnRightString = aice_sh_trends_seas@units
        res@gsnCenterString = names_sh(ee)
        plot_trends_sh_djf(ee) = gsn_csm_contour_map(wks_trends_djf,aice_sh_trends_seas(0,:,:),res)
        plot_trends_sh_mam(ee) = gsn_csm_contour_map(wks_trends_mam,aice_sh_trends_seas(1,:,:),res)
        plot_trends_sh_jja(ee) = gsn_csm_contour_map(wks_trends_jja,aice_sh_trends_seas(2,:,:),res)
        plot_trends_sh_son(ee) = gsn_csm_contour_map(wks_trends_son,aice_sh_trends_seas(3,:,:),res)
        plot_trends_sh_ann(ee) = gsn_csm_contour_map(wks_trends_ann,aice_sh_trends_ann,res)
        plot_trends_sh_mon(ee) = gsn_csm_contour_map(wks_trends_mon,aice_sh_trends_mon,res)
        delete([/aice_sh_trends_seas,aice_sh_trends_ann,aice_sh_trends_mon/])

        if (isvar("aice_sh_sum_djf")) then
           xyres@gsnLeftString = names_sh(ee)  
           if (isfilepresent2("obs_aice_sh").and.ee.eq.0) then
              aice_sh_sum_djf_obs_min = min(aice_sh_sum_djf)
              aice_sh_sum_djf_obs_max = max(aice_sh_sum_djf)
              aice_sh_sum_mam_obs_min = min(aice_sh_sum_mam)
              aice_sh_sum_mam_obs_max = max(aice_sh_sum_mam)
              aice_sh_sum_jja_obs_min = min(aice_sh_sum_jja)
              aice_sh_sum_jja_obs_max = max(aice_sh_sum_jja)
              aice_sh_sum_son_obs_min = min(aice_sh_sum_son)
              aice_sh_sum_son_obs_max = max(aice_sh_sum_son)
              aice_sh_sum_ann_obs_min = min(aice_sh_sum_ann)
              aice_sh_sum_ann_obs_max = max(aice_sh_sum_ann)
              aice_sh_sum_mon_obs_min = min(aice_sh_sum_mon)
              aice_sh_sum_mon_obs_max = max(aice_sh_sum_mon)
              aice_sh_sum_mon_anom_obs_min = min(aice_sh_sum_mon_anom)
              aice_sh_sum_mon_anom_obs_max = max(aice_sh_sum_mon_anom)
              aice_sh_sum_feb_obs_min = min(aice_sh_sum_feb)
              aice_sh_sum_feb_obs_max = max(aice_sh_sum_feb)
              aice_sh_sum_mar_obs_min = min(aice_sh_sum_mar)
              aice_sh_sum_mar_obs_max = max(aice_sh_sum_mar)
              aice_sh_sum_sep_obs_min = min(aice_sh_sum_sep)
              aice_sh_sum_sep_obs_max = max(aice_sh_sum_sep)
              aice_sh_obs_djf = aice_sh_sum_djf
              aice_sh_obs_mam = aice_sh_sum_mam
              aice_sh_obs_jja = aice_sh_sum_jja
              aice_sh_obs_son = aice_sh_sum_son
              aice_sh_obs_ann = aice_sh_sum_ann
              aice_sh_obs_mon = aice_sh_sum_mon
              aice_sh_obs_mon_anom = aice_sh_sum_mon_anom
              aice_sh_obs_feb = aice_sh_sum_feb
              aice_sh_obs_mar = aice_sh_sum_mar
              aice_sh_obs_sep = aice_sh_sum_sep
              aice_sh_obs_climo = aice_sh_sum_climo
           end if

           xyres_c@gsnLeftString = syear_sh(ee)+"-"+eyear_sh(ee)  
           xyres_c@gsnCenterString = names_sh(ee)  
           xyres_c@gsnRightString = "("+str_sub_str(str_sub_str(aice_sh_sum_climo@units,"m2","m~S~2~N~"),"10^12","10~S~12~N~")+")"
           if (ee.eq.0) then
              plot_iceext_sh_climo(ee) = gsn_csm_xy(wks_iceext_mon,ispan(1,12,1),aice_sh_sum_climo,xyres_c)
           end if    
           if (ee.ge.1.) then
              if (isvar("aice_sh_obs_climo")) then
                 tarr = new((/2,12/),float)
                 tarr(0,:) = (/ tofloat(aice_sh_sum_climo) /)
                 tarr(1,:) = (/ tofloat(aice_sh_obs_climo) /)
                 plot_iceext_sh_climo(ee) = gsn_csm_xy(wks_iceext_mon,ispan(1,12,1),tarr,xyres_c)
                 delete(tarr)
              else
                 plot_iceext_sh_climo(ee) = gsn_csm_xy(wks_iceext_mon,ispan(1,12,1),aice_sh_sum_climo,xyres_c)
              end if
           end if

           if (ee.ge.1.and.isvar("aice_sh_sum_djf_obs_min")) then
              xyres@trYMinF = min((/min(aice_sh_sum_djf),aice_sh_sum_djf_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_sh_sum_djf),aice_sh_sum_djf_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_sh_sum_djf)-1,1),aice_sh_sum_djf,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr_sh(ee),2,True)+aice_sh_sum_djf@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_sh_djf(ee)     = gsn_csm_xy(wks_iceext_djf,ispan(syear_sh(ee),eyear_sh(ee),1),aice_sh_sum_djf,xyres)
           if (ee.ge.1.and.isvar("aice_sh_obs_djf")) then
              plot_iceext_sh_obs_djf = gsn_csm_xy(wks_iceext_djf,ispan(syear_sh(0),eyear_sh(0),1),aice_sh_obs_djf,xyres2)
              overlay(plot_iceext_sh_djf(ee),plot_iceext_sh_obs_djf)
           end if

           if (ee.ge.1.and.isvar("aice_sh_sum_mam_obs_min")) then
              xyres@trYMinF = min((/min(aice_sh_sum_mam),aice_sh_sum_mam_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_sh_sum_mam),aice_sh_sum_mam_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_sh_sum_mam)-1,1),aice_sh_sum_mam,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr_sh(ee),2,True)+aice_sh_sum_mam@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_sh_mam(ee)     = gsn_csm_xy(wks_iceext_mam,ispan(syear_sh(ee),eyear_sh(ee),1),aice_sh_sum_mam,xyres)
           if (ee.ge.1.and.isvar("aice_sh_obs_mam")) then
              plot_iceext_sh_obs_mam = gsn_csm_xy(wks_iceext_mam,ispan(syear_sh(0),eyear_sh(0),1),aice_sh_obs_mam,xyres2)
              overlay(plot_iceext_sh_mam(ee),plot_iceext_sh_obs_mam)
           end if      

           if (ee.ge.1.and.isvar("aice_sh_sum_jja_obs_min")) then
              xyres@trYMinF = min((/min(aice_sh_sum_jja),aice_sh_sum_jja_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_sh_sum_jja),aice_sh_sum_jja_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_sh_sum_jja)-1,1),aice_sh_sum_jja,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr_sh(ee),2,True)+aice_sh_sum_jja@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_sh_jja(ee)     = gsn_csm_xy(wks_iceext_jja,ispan(syear_sh(ee),eyear_sh(ee),1),aice_sh_sum_jja,xyres)
           if (ee.ge.1.and.isvar("aice_sh_obs_jja")) then
              plot_iceext_sh_obs_jja = gsn_csm_xy(wks_iceext_jja,ispan(syear_sh(0),eyear_sh(0),1),aice_sh_obs_jja,xyres2)
              overlay(plot_iceext_sh_jja(ee),plot_iceext_sh_obs_jja)
           end if

           if (ee.ge.1.and.isvar("aice_sh_sum_son_obs_min")) then
              xyres@trYMinF = min((/min(aice_sh_sum_son),aice_sh_sum_son_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_sh_sum_son),aice_sh_sum_son_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_sh_sum_son)-1,1),aice_sh_sum_son,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr_sh(ee),2,True)+aice_sh_sum_son@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_sh_son(ee)     = gsn_csm_xy(wks_iceext_son,ispan(syear_sh(ee),eyear_sh(ee),1),aice_sh_sum_son,xyres)
           if (ee.ge.1.and.isvar("aice_sh_obs_son")) then
              plot_iceext_sh_obs_son = gsn_csm_xy(wks_iceext_son,ispan(syear_sh(0),eyear_sh(0),1),aice_sh_obs_son,xyres2)
              overlay(plot_iceext_sh_son(ee),plot_iceext_sh_obs_son)
           end if

           if (ee.ge.1.and.isvar("aice_sh_sum_ann_obs_min")) then
              xyres@trYMinF = min((/min(aice_sh_sum_ann),aice_sh_sum_ann_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_sh_sum_ann),aice_sh_sum_ann_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_sh_sum_ann)-1,1),aice_sh_sum_ann,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr_sh(ee),2,True)+aice_sh_sum_ann@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_sh_ann(ee)     = gsn_csm_xy(wks_iceext_ann,ispan(syear_sh(ee),eyear_sh(ee),1),aice_sh_sum_ann,xyres)
           if (ee.ge.1.and.isvar("aice_sh_obs_ann")) then
              plot_iceext_sh_obs_ann = gsn_csm_xy(wks_iceext_ann,ispan(syear_sh(0),eyear_sh(0),1),aice_sh_obs_ann,xyres2)
              overlay(plot_iceext_sh_ann(ee),plot_iceext_sh_obs_ann)
           end if
           delete(tttt)

           if (ee.ge.1.and.isvar("aice_sh_sum_feb_obs_min")) then
              xyres@trYMinF = min((/min(aice_sh_sum_feb),aice_sh_sum_feb_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_sh_sum_feb),aice_sh_sum_feb_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_sh_sum_feb)-1,1),aice_sh_sum_feb,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr_sh(ee),2,True)+aice_sh_sum_feb@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_sh_feb(ee)     = gsn_csm_xy(wks_iceext_febmar,ispan(syear_sh(ee),eyear_sh(ee),1),aice_sh_sum_feb,xyres)
           if (ee.ge.1.and.isvar("aice_sh_obs_feb")) then
              plot_iceext_sh_obs_feb = gsn_csm_xy(wks_iceext_febmar,ispan(syear_sh(0),eyear_sh(0),1),aice_sh_obs_feb,xyres2)
              overlay(plot_iceext_sh_feb(ee),plot_iceext_sh_obs_feb)
           end if
           delete(tttt)

           if (ee.ge.1.and.isvar("aice_sh_sum_mar_obs_min")) then
              xyres@trYMinF = min((/min(aice_sh_sum_mar),aice_sh_sum_mar_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_sh_sum_mar),aice_sh_sum_mar_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_sh_sum_mar)-1,1),aice_sh_sum_mar,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr_sh(ee),2,True)+aice_sh_sum_mar@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_sh_mar(ee)     = gsn_csm_xy(wks_iceext_febmar,ispan(syear_sh(ee),eyear_sh(ee),1),aice_sh_sum_mar,xyres)
           if (ee.ge.1.and.isvar("aice_sh_obs_mar")) then
              plot_iceext_sh_obs_mar = gsn_csm_xy(wks_iceext_febmar,ispan(syear_sh(0),eyear_sh(0),1),aice_sh_obs_mar,xyres2)
              overlay(plot_iceext_sh_mar(ee),plot_iceext_sh_obs_mar)
           end if
           delete(tttt)

           if (ee.ge.1.and.isvar("aice_sh_sum_sep_obs_min")) then
              xyres@trYMinF = min((/min(aice_sh_sum_sep),aice_sh_sum_sep_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_sh_sum_sep),aice_sh_sum_sep_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_sh_sum_sep)-1,1),aice_sh_sum_sep,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr_sh(ee),2,True)+aice_sh_sum_sep@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_sh_sep(ee)     = gsn_csm_xy(wks_iceext_sep,ispan(syear_sh(ee),eyear_sh(ee),1),aice_sh_sum_sep,xyres)
           if (ee.ge.1.and.isvar("aice_sh_obs_sep")) then
              plot_iceext_sh_obs_sep = gsn_csm_xy(wks_iceext_sep,ispan(syear_sh(0),eyear_sh(0),1),aice_sh_obs_sep,xyres2)
              overlay(plot_iceext_sh_sep(ee),plot_iceext_sh_obs_sep)
           end if
           delete(tttt)

           if (ee.ge.1.and.isvar("aice_sh_sum_mon_obs_min")) then
              xyres@trYMinF = min((/min(aice_sh_sum_mon),aice_sh_sum_mon_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_sh_sum_mon),aice_sh_sum_mon_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_sh_sum_mon)-1,1),aice_sh_sum_mon,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr_sh(ee),2,True)+aice_sh_sum_mon@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_sh_mon(ee)     = gsn_csm_xy(wks_iceext_mon,fspan(syear_sh(ee),eyear_sh(ee)+.91667,dimsizes(aice_sh_sum_mon)),aice_sh_sum_mon,xyres)
           if (ee.ge.1.and.isvar("aice_sh_obs_mon")) then
              plot_iceext_sh_obs_mon = gsn_csm_xy(wks_iceext_mon,fspan(syear_sh(0),eyear_sh(0)+.91667,dimsizes(aice_sh_obs_mon)),aice_sh_obs_mon,xyres2)
              overlay(plot_iceext_sh_mon(ee),plot_iceext_sh_obs_mon)
           end if
           delete(tttt)
         
           if (ee.ge.1.and.isvar("aice_sh_sum_mon_anom_obs_min")) then
              xyres@trYMinF = min((/min(aice_sh_sum_mon_anom),aice_sh_sum_mon_anom_obs_min/))-1
              xyres@trYMaxF = max((/max(aice_sh_sum_mon_anom),aice_sh_sum_mon_anom_obs_max/))+1
              xyres2@trYMinF = xyres@trYMinF
              xyres2@trYMaxF = xyres@trYMaxF
           end if
           tttt = dtrend_msg(ispan(0,dimsizes(aice_sh_sum_mon_anom)-1,1),aice_sh_sum_mon_anom,False,True) 
           xyres@gsnRightString = str_sub_str(str_sub_str(decimalPlaces(tttt@slope*nyr_sh(ee),2,True)+aice_sh_sum_mon_anom@units,"m2","m~S~2~N~"),"10^12"," 10~S~12~N~")
           plot_iceext_sh_mon_anom(ee)     = gsn_csm_xy(wks_iceext_mon,fspan(syear_sh(ee),eyear_sh(ee)+.91667,dimsizes(aice_sh_sum_mon_anom)),aice_sh_sum_mon_anom,xyres)
           if (ee.ge.1.and.isvar("aice_sh_obs_mon_anom")) then
              plot_iceext_sh_obs_mon_anom = gsn_csm_xy(wks_iceext_mon,fspan(syear_sh(0),eyear_sh(0)+.91667,dimsizes(aice_sh_obs_mon_anom)),aice_sh_obs_mon_anom,xyres2)
              overlay(plot_iceext_sh_mon_anom(ee),plot_iceext_sh_obs_mon_anom)
           end if

           delete([/aice_sh_sum_djf,aice_sh_sum_mam,aice_sh_sum_jja,aice_sh_sum_son,aice_sh_sum_ann,aice_sh_sum_feb,aice_sh_sum_mar,aice_sh_sum_sep,aice_sh_sum_mon,aice_sh_sum_mon_anom,tttt/])
        end if
     end if
     delete([/res,xyres,xyres2,xyres_c,ph_s/])
;     print("Done with ee = "+ee)
  end do     
  delete(time_mon2)
     
  panres = True
  panres@gsnMaximize = True  
  panres@gsnPaperOrientation = "portrait"
  panres@gsnPanelLabelBar = True
  panres@gsnPanelYWhiteSpacePercent = 3.0
  panres@pmLabelBarHeightF = 0.05
  panres@pmLabelBarWidthF = 0.65
  panres@lbTitleOn = False
  panres@lbBoxLineColor = "gray70"
  panres@lbLabelFontHeightF = 0.013
  panres@lbLabelStride = 1
  if (nsim.le.4) then
     if (nsim.eq.1) then
        panres@txFontHeightF = 0.0185
        panres@gsnPanelBottom = 0.50
        panres@pmLabelBarWidthF = 0.5
        panres@pmLabelBarHeightF = 0.04
        panres@lbLabelFontHeightF = 0.01
     else
        panres@txFontHeightF = 0.0125
        panres@gsnPanelBottom = 0.50
        panres@pmLabelBarWidthF = 0.5
        panres@pmLabelBarHeightF = 0.04
        panres@lbLabelFontHeightF = 0.01
     end if
  else
     panres@txFontHeightF = 0.016
     panres@gsnPanelBottom = 0.05
  end if

  ncol = floattointeger(sqrt(nsim))
  nrow = (nsim/ncol)+mod(nsim,ncol) 

  ncol_sh = floattointeger(sqrt(nsim_sh))
  nrow_sh = (nsim_sh/ncol_sh)+mod(nsim_sh,ncol_sh)   
  
  panres@txString = "SIC Trends (DJF)"  
  gsn_panel2(wks_trends_djf,plot_trends_nh_djf,(/nrow,ncol/),panres)
  gsn_panel2(wks_trends_djf,plot_trends_sh_djf,(/nrow_sh,ncol_sh/),panres)
  delete(wks_trends_djf)
  
  panres@txString = "SIC Trends (MAM)"
  gsn_panel2(wks_trends_mam,plot_trends_nh_mam,(/nrow,ncol/),panres)
  gsn_panel2(wks_trends_mam,plot_trends_sh_mam,(/nrow_sh,ncol_sh/),panres)
  delete(wks_trends_mam)
  
  panres@txString = "SIC Trends (JJA)"
  gsn_panel2(wks_trends_jja,plot_trends_nh_jja,(/nrow,ncol/),panres)
  gsn_panel2(wks_trends_jja,plot_trends_sh_jja,(/nrow_sh,ncol_sh/),panres)
  delete(wks_trends_jja)
  
  panres@txString = "SIC Trends (SON)"
  gsn_panel2(wks_trends_son,plot_trends_nh_son,(/nrow,ncol/),panres)
  gsn_panel2(wks_trends_son,plot_trends_sh_son,(/nrow_sh,ncol_sh/),panres)
  delete(wks_trends_son)
  
  panres@txString = "SIC Trends (Annual)"
  gsn_panel2(wks_trends_ann,plot_trends_nh_ann,(/nrow,ncol/),panres)
  gsn_panel2(wks_trends_ann,plot_trends_sh_ann,(/nrow_sh,ncol_sh/),panres)
  delete(wks_trends_ann)
  
  panres@txString = "SIC Trends (Monthly)"
  gsn_panel2(wks_trends_mon,plot_trends_nh_mon,(/nrow,ncol/),panres)
  gsn_panel2(wks_trends_mon,plot_trends_sh_mon,(/nrow_sh,ncol_sh/),panres)
  delete(wks_trends_mon)
  delete(panres)

  panres2 = True
  panres2@gsnMaximize = True
  panres2@gsnPaperOrientation = "portrait"
  panres2@gsnPanelYWhiteSpacePercent = 3.0  
  if (nsim.le.5) then
     panres2@txFontHeightF = 0.024
  else
     panres2@txFontHeightF = 0.016
  end if
  panres3 = panres2
  if (SCALE_TIMESERIES.eq."True") then
     tt = ind(nyr.eq.nyr_max)
     panres2@gsnPanelScalePlotIndex = tt(0)
     delete(tt)

     tt = ind(nyr_sh.eq.nyr_max_sh)
     panres3@gsnPanelScalePlotIndex = tt(0)
     delete(tt)
  end if

  if (nsim.le.12) then
     lp    = (/nsim,1/)
     lp_sh = (/nsim_sh,1/)
  else
     lp    = (/nrow,ncol/)   ;(/nsim/2+1,nsim/8+1/)
     lp_sh = (/nrow_sh,ncol_sh/)  
  end if

  panres2@txString = "SIC NH Extent (DJF)"
  gsn_panel2(wks_iceext_djf,plot_iceext_nh_djf,lp,panres2)  
  panres3@txString = "SIC SH Extent (DJF)"
  gsn_panel2(wks_iceext_djf,plot_iceext_sh_djf,lp_sh,panres3)  
  delete(wks_iceext_djf)
  
  panres2@txString = "SIC NH Extent (MAM)"
  gsn_panel2(wks_iceext_mam,plot_iceext_nh_mam,lp,panres2)  
  panres3@txString = "SIC SH Extent (MAM)"
  gsn_panel2(wks_iceext_mam,plot_iceext_sh_mam,lp_sh,panres3)    
  delete(wks_iceext_mam)
  
  panres2@txString = "SIC NH Extent (JJA)"
  gsn_panel2(wks_iceext_jja,plot_iceext_nh_jja,lp,panres2)  
  panres3@txString = "SIC SH Extent (JJA)"
  gsn_panel2(wks_iceext_jja,plot_iceext_sh_jja,lp_sh,panres3)     
  delete(wks_iceext_jja)
  
  panres2@txString = "SIC NH Extent (SON)"
  gsn_panel2(wks_iceext_son,plot_iceext_nh_son,lp,panres2)  
  panres3@txString = "SIC SH Extent (SON)"
  gsn_panel2(wks_iceext_son,plot_iceext_sh_son,lp_sh,panres3)   
  delete(wks_iceext_son)
  
  panres2@txString = "SIC NH Extent (Annual)"
  gsn_panel2(wks_iceext_ann,plot_iceext_nh_ann,lp,panres2)  
  panres3@txString = "SIC SH Extent (Annual)"
  gsn_panel2(wks_iceext_ann,plot_iceext_sh_ann,lp_sh,panres3)     
  delete(wks_iceext_ann)

  panres2@txString = "SIC NH Extent (February)"
  gsn_panel2(wks_iceext_febmar,plot_iceext_nh_mar,lp,panres2)  
  panres3@txString = "SIC SH Extent (Feburary)"
  gsn_panel2(wks_iceext_febmar,plot_iceext_sh_mar,lp_sh,panres3)     

  panres2@txString = "SIC NH Extent (March)"
  gsn_panel2(wks_iceext_febmar,plot_iceext_nh_mar,lp,panres2)  
  panres3@txString = "SIC SH Extent (March)"
  gsn_panel2(wks_iceext_febmar,plot_iceext_sh_mar,lp_sh,panres3)     
  delete(wks_iceext_febmar)

  panres2@txString = "SIC NH Extent (September)"
  gsn_panel2(wks_iceext_sep,plot_iceext_nh_sep,lp,panres2)  
  panres3@txString = "SIC SH Extent (September)"
  gsn_panel2(wks_iceext_sep,plot_iceext_sh_sep,lp_sh,panres3)     
  delete(wks_iceext_sep)
  
  panres2@txString = "SIC NH Extent (Monthly)"
  gsn_panel2(wks_iceext_mon,plot_iceext_nh_mon,lp,panres2)  
  panres3@txString = "SIC SH Extent (Monthly)"
  gsn_panel2(wks_iceext_mon,plot_iceext_sh_mon,lp_sh,panres3)  
  panres2@txString = "SIC NH Extent (Monthly Anomalies)"
  gsn_panel2(wks_iceext_mon,plot_iceext_nh_mon_anom,lp,panres2)  
  panres3@txString = "SIC SH Extent (Monthly Anomalies)"
  gsn_panel2(wks_iceext_mon,plot_iceext_sh_mon_anom,lp_sh,panres3)  

  if (nsim.le.5) then
     panres2@txFontHeightF = 0.017
  else
     panres2@txFontHeightF = 0.014
  end if
  ncol = floattointeger(sqrt(nsim))
  nrow = (nsim/ncol)+mod(nsim,ncol) 

  ncol_sh = floattointeger(sqrt(nsim_sh))
  nrow_sh = (nsim_sh/ncol_sh)+mod(nsim_sh,ncol_sh)

  panres2@txString = "SIC NH Extent Climatology"
  gsn_panel2(wks_iceext_mon,plot_iceext_nh_climo,(/nrow,ncol/),panres2)  
  panres2@txString = "SIC SH Extent Climatology"
  gsn_panel2(wks_iceext_mon,plot_iceext_sh_climo,(/nrow_sh,ncol_sh/),panres2)        
  delete(wks_iceext_mon)
;--------------------------------------------------------------------------------
  OUTDIR = getenv("OUTDIR") 
  if (wks_type.eq."png") then  
     system("mv "+OUTDIR+"aice.trends.djf.000001.png "+OUTDIR+"aice.trends.nh.djf.png") 
     system("mv "+OUTDIR+"aice.trends.djf.000002.png "+OUTDIR+"aice.trends.sh.djf.png") 
     system("mv "+OUTDIR+"aice.trends.mam.000001.png "+OUTDIR+"aice.trends.nh.mam.png") 
     system("mv "+OUTDIR+"aice.trends.mam.000002.png "+OUTDIR+"aice.trends.sh.mam.png")
     system("mv "+OUTDIR+"aice.trends.jja.000001.png "+OUTDIR+"aice.trends.nh.jja.png") 
     system("mv "+OUTDIR+"aice.trends.jja.000002.png "+OUTDIR+"aice.trends.sh.jja.png")
     system("mv "+OUTDIR+"aice.trends.son.000001.png "+OUTDIR+"aice.trends.nh.son.png") 
     system("mv "+OUTDIR+"aice.trends.son.000002.png "+OUTDIR+"aice.trends.sh.son.png") 
     system("mv "+OUTDIR+"aice.trends.ann.000001.png "+OUTDIR+"aice.trends.nh.ann.png") 
     system("mv "+OUTDIR+"aice.trends.ann.000002.png "+OUTDIR+"aice.trends.sh.ann.png")
     system("mv "+OUTDIR+"aice.trends.mon.000001.png "+OUTDIR+"aice.trends.nh.mon.png") 
     system("mv "+OUTDIR+"aice.trends.mon.000002.png "+OUTDIR+"aice.trends.sh.mon.png")

     if (isfilepresent2(OUTDIR+"aice.extent.djf.000001.png")) then
        system("mv "+OUTDIR+"aice.extent.djf.000001.png "+OUTDIR+"aice.extent.nh.djf.png") 
        system("mv "+OUTDIR+"aice.extent.djf.000002.png "+OUTDIR+"aice.extent.sh.djf.png") 
        system("mv "+OUTDIR+"aice.extent.mam.000001.png "+OUTDIR+"aice.extent.nh.mam.png") 
        system("mv "+OUTDIR+"aice.extent.mam.000002.png "+OUTDIR+"aice.extent.sh.mam.png")
        system("mv "+OUTDIR+"aice.extent.jja.000001.png "+OUTDIR+"aice.extent.nh.jja.png") 
        system("mv "+OUTDIR+"aice.extent.jja.000002.png "+OUTDIR+"aice.extent.sh.jja.png")
        system("mv "+OUTDIR+"aice.extent.son.000001.png "+OUTDIR+"aice.extent.nh.son.png") 
        system("mv "+OUTDIR+"aice.extent.son.000002.png "+OUTDIR+"aice.extent.sh.son.png") 
        system("mv "+OUTDIR+"aice.extent.ann.000001.png "+OUTDIR+"aice.extent.nh.ann.png") 
        system("mv "+OUTDIR+"aice.extent.ann.000002.png "+OUTDIR+"aice.extent.sh.ann.png")  
        system("mv "+OUTDIR+"aice.extent.febmar.000001.png "+OUTDIR+"aice.extent.nh.feb.png") 
        system("mv "+OUTDIR+"aice.extent.febmar.000002.png "+OUTDIR+"aice.extent.sh.feb.png")  
        system("mv "+OUTDIR+"aice.extent.febmar.000003.png "+OUTDIR+"aice.extent.nh.mar.png") 
        system("mv "+OUTDIR+"aice.extent.febmar.000004.png "+OUTDIR+"aice.extent.sh.mar.png")  
        system("mv "+OUTDIR+"aice.extent.sep.000001.png "+OUTDIR+"aice.extent.nh.sep.png") 
        system("mv "+OUTDIR+"aice.extent.sep.000002.png "+OUTDIR+"aice.extent.sh.sep.png")  
        system("mv "+OUTDIR+"aice.extent.mon.000001.png "+OUTDIR+"aice.extent.nh.mon.png") 
        system("mv "+OUTDIR+"aice.extent.mon.000002.png "+OUTDIR+"aice.extent.sh.mon.png") 
        system("mv "+OUTDIR+"aice.extent.mon.000003.png "+OUTDIR+"aice.extent.anom.nh.mon.png") 
        system("mv "+OUTDIR+"aice.extent.mon.000004.png "+OUTDIR+"aice.extent.anom.sh.mon.png") 
        system("mv "+OUTDIR+"aice.extent.mon.000005.png "+OUTDIR+"aice.extent.nh.climo.png") 
        system("mv "+OUTDIR+"aice.extent.mon.000006.png "+OUTDIR+"aice.extent.sh.climo.png")
     end if 
  else
     system("psplit "+OUTDIR+"aice.trends.djf.ps "+OUTDIR+"aice_tr")
     system("mv "+OUTDIR+"aice_tr0001.ps "+OUTDIR+"aice.trends.nh.djf.ps") 
     system("mv "+OUTDIR+"aice_tr0002.ps "+OUTDIR+"aice.trends.sh.djf.ps")
     system("psplit "+OUTDIR+"aice.trends.mam.ps "+OUTDIR+"aice_tr")
     system("mv "+OUTDIR+"aice_tr0001.ps "+OUTDIR+"aice.trends.nh.mam.ps") 
     system("mv "+OUTDIR+"aice_tr0002.ps "+OUTDIR+"aice.trends.sh.mam.ps")
     system("psplit "+OUTDIR+"aice.trends.jja.ps "+OUTDIR+"aice_tr")
     system("mv "+OUTDIR+"aice_tr0001.ps "+OUTDIR+"aice.trends.nh.jja.ps") 
     system("mv "+OUTDIR+"aice_tr0002.ps "+OUTDIR+"aice.trends.sh.jja.ps")
     system("psplit "+OUTDIR+"aice.trends.son.ps "+OUTDIR+"aice_tr")
     system("mv "+OUTDIR+"aice_tr0001.ps "+OUTDIR+"aice.trends.nh.son.ps") 
     system("mv "+OUTDIR+"aice_tr0002.ps "+OUTDIR+"aice.trends.sh.son.ps")
     system("psplit "+OUTDIR+"aice.trends.ann.ps "+OUTDIR+"aice_tr")
     system("mv "+OUTDIR+"aice_tr0001.ps "+OUTDIR+"aice.trends.nh.ann.ps") 
     system("mv "+OUTDIR+"aice_tr0002.ps "+OUTDIR+"aice.trends.sh.ann.ps")
     system("psplit "+OUTDIR+"aice.trends.mon.ps "+OUTDIR+"aice_tr")
     system("mv "+OUTDIR+"aice_tr0001.ps "+OUTDIR+"aice.trends.nh.mon.ps") 
     system("mv "+OUTDIR+"aice_tr0002.ps "+OUTDIR+"aice.trends.sh.mon.ps")
     system("rm "+OUTDIR+"aice.trends.???.ps")

     if (isfilepresent2(OUTDIR+"aice.extent.djf.ps")) then
        system("psplit "+OUTDIR+"aice.extent.djf.ps "+OUTDIR+"aice_tr")
        system("mv "+OUTDIR+"aice_tr0001.ps "+OUTDIR+"aice.extent.nh.djf.ps") 
        system("mv "+OUTDIR+"aice_tr0002.ps "+OUTDIR+"aice.extent.sh.djf.ps")
        system("psplit "+OUTDIR+"aice.extent.mam.ps "+OUTDIR+"aice_tr")
        system("mv "+OUTDIR+"aice_tr0001.ps "+OUTDIR+"aice.extent.nh.mam.ps") 
        system("mv "+OUTDIR+"aice_tr0002.ps "+OUTDIR+"aice.extent.sh.mam.ps")
        system("psplit "+OUTDIR+"aice.extent.jja.ps "+OUTDIR+"aice_tr")
        system("mv "+OUTDIR+"aice_tr0001.ps "+OUTDIR+"aice.extent.nh.jja.ps") 
        system("mv "+OUTDIR+"aice_tr0002.ps "+OUTDIR+"aice.extent.sh.jja.ps")
        system("psplit "+OUTDIR+"aice.extent.son.ps "+OUTDIR+"aice_tr")
        system("mv "+OUTDIR+"aice_tr0001.ps "+OUTDIR+"aice.extent.nh.son.ps") 
        system("mv "+OUTDIR+"aice_tr0002.ps "+OUTDIR+"aice.extent.sh.son.ps")
        system("psplit "+OUTDIR+"aice.extent.ann.ps "+OUTDIR+"aice_tr")
        system("mv "+OUTDIR+"aice_tr0001.ps "+OUTDIR+"aice.extent.nh.ann.ps") 
        system("mv "+OUTDIR+"aice_tr0002.ps "+OUTDIR+"aice.extent.sh.ann.ps")
        system("psplit "+OUTDIR+"aice.extent.febmar.ps "+OUTDIR+"aice_tr")
        system("mv "+OUTDIR+"aice_tr0001.ps "+OUTDIR+"aice.extent.nh.feb.ps") 
        system("mv "+OUTDIR+"aice_tr0002.ps "+OUTDIR+"aice.extent.sh.feb.ps")
        system("mv "+OUTDIR+"aice_tr0003.ps "+OUTDIR+"aice.extent.nh.mar.ps") 
        system("mv "+OUTDIR+"aice_tr0004.ps "+OUTDIR+"aice.extent.sh.mar.ps")
        system("psplit "+OUTDIR+"aice.extent.sep.ps "+OUTDIR+"aice_tr")
        system("mv "+OUTDIR+"aice_tr0001.ps "+OUTDIR+"aice.extent.nh.sep.ps") 
        system("mv "+OUTDIR+"aice_tr0002.ps "+OUTDIR+"aice.extent.sh.sep.ps")
        system("psplit "+OUTDIR+"aice.extent.mon.ps "+OUTDIR+"aice_tr")
        system("mv "+OUTDIR+"aice_tr0001.ps "+OUTDIR+"aice.extent.nh.mon.ps") 
        system("mv "+OUTDIR+"aice_tr0002.ps "+OUTDIR+"aice.extent.sh.mon.ps")
        system("mv "+OUTDIR+"aice_tr0003.ps "+OUTDIR+"aice.extent.anom.nh.mon.ps") 
        system("mv "+OUTDIR+"aice_tr0004.ps "+OUTDIR+"aice.extent.anom.sh.mon.ps")
        system("mv "+OUTDIR+"aice_tr0005.ps "+OUTDIR+"aice.extent.nh.climo.ps") 
        system("mv "+OUTDIR+"aice_tr0006.ps "+OUTDIR+"aice.extent.sh.climo.ps")
        system("rm "+OUTDIR+"aice.extent.???.ps "+OUTDIR+"aice.extent.febmar.ps")
     end if
  end if
  print("Finished: aice.trends_timeseries.ncl")
end
