; Calculates SIC hemispheric means and standard deviations
;
; Variables used: sic
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CVDP_SCRIPTS/functions.ncl"

begin
  print("Starting: aice.mean_stddev.ncl")
  
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
  wks_stddev_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.stddev.djf")
  wks_stddev_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.stddev.mam")
  wks_stddev_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.stddev.jja")
  wks_stddev_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.stddev.son")
  wks_stddev_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.stddev.ann")
  wks_mean_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.mean.djf")
  wks_mean_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.mean.mam")
  wks_mean_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.mean.jja")
  wks_mean_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.mean.son")
  wks_mean_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"aice.mean.ann")

  if (COLORMAP.eq.0) then
     gsn_define_colormap(wks_stddev_djf,"BkBlAqGrYeOrReViWh200")    
     gsn_define_colormap(wks_stddev_mam,"BkBlAqGrYeOrReViWh200")   
     gsn_define_colormap(wks_stddev_jja,"BkBlAqGrYeOrReViWh200")   
     gsn_define_colormap(wks_stddev_son,"BkBlAqGrYeOrReViWh200")   
     gsn_define_colormap(wks_stddev_ann,"BkBlAqGrYeOrReViWh200")   
     gsn_define_colormap(wks_mean_djf,"BkBlAqGrYeOrReViWh200")   
     gsn_define_colormap(wks_mean_mam,"BkBlAqGrYeOrReViWh200")  
     gsn_define_colormap(wks_mean_jja,"BkBlAqGrYeOrReViWh200") 
     gsn_define_colormap(wks_mean_son,"BkBlAqGrYeOrReViWh200") 
     gsn_define_colormap(wks_mean_ann,"BkBlAqGrYeOrReViWh200") 
  end if
  if (COLORMAP.eq.1) then
     gsn_define_colormap(wks_stddev_djf,"cb_rainbow")     
     gsn_define_colormap(wks_stddev_mam,"cb_rainbow")    
     gsn_define_colormap(wks_stddev_jja,"cb_rainbow")   
     gsn_define_colormap(wks_stddev_son,"cb_rainbow")   
     gsn_define_colormap(wks_stddev_ann,"cb_rainbow")  
     gsn_define_colormap(wks_mean_djf,"cb_rainbow")     
     gsn_define_colormap(wks_mean_mam,"cb_rainbow")    
     gsn_define_colormap(wks_mean_jja,"cb_rainbow")   
     gsn_define_colormap(wks_mean_son,"cb_rainbow")   
     gsn_define_colormap(wks_mean_ann,"cb_rainbow")   
  end if


  plot_mean_nh_djf = new(nsim,"graphic")  
  plot_mean_nh_mam = new(nsim,"graphic")  
  plot_mean_nh_jja = new(nsim,"graphic")  
  plot_mean_nh_son = new(nsim,"graphic")   
  plot_mean_nh_ann = new(nsim,"graphic")  
  plot_stddev_nh_djf = new(nsim,"graphic")  
  plot_stddev_nh_mam = new(nsim,"graphic")  
  plot_stddev_nh_jja = new(nsim,"graphic")  
  plot_stddev_nh_son = new(nsim,"graphic")   
  plot_stddev_nh_ann = new(nsim,"graphic")

  plot_mean_sh_djf = new(nsim,"graphic")  
  plot_mean_sh_mam = new(nsim,"graphic")  
  plot_mean_sh_jja = new(nsim,"graphic")  
  plot_mean_sh_son = new(nsim,"graphic")   
  plot_mean_sh_ann = new(nsim,"graphic")  
  plot_stddev_sh_djf = new(nsim,"graphic")  
  plot_stddev_sh_mam = new(nsim,"graphic")  
  plot_stddev_sh_jja = new(nsim,"graphic")  
  plot_stddev_sh_son = new(nsim,"graphic")   
  plot_stddev_sh_ann = new(nsim,"graphic")
  
  do ee = 0,nsim-1
     aice_nh_flag = 0
     aice_nh = data_read_in_ice(paths(ee),"aice_nh",syear(ee),eyear(ee))    ; read in data, orient lats/lons correctly, set time coordinate variable up
     if (isatt(aice_nh,"is_all_missing")) then
        delete(aice_nh)
        aice_nh_flag = 1
     end if  

     if (aice_nh_flag.eq.0) then
        do ff = 0,1     
           aice_nhT = aice_nh
           if (ff.eq.1) then
              if (OPT_CLIMO.eq."Full") then
                 aice_nhT = rmMonAnnCycTLL(aice_nhT)
              else
                 check_custom_climo(names(ee),syear(ee),eyear(ee),CLIMO_SYEAR,CLIMO_EYEAR)
                 temp_arr = aice_nhT
                 delete(temp_arr&time)
                 temp_arr&time = cd_calendar(aice_nhT&time,-1)
                 if (CLIMO_SYEAR.lt.0) then
                    climo = clmMonTLL(temp_arr({(eyear(ee)+CLIMO_SYEAR)*100+1:(eyear(ee)+CLIMO_EYEAR)*100+12},:,:))       
                 else
                    climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
                 end if               
                 delete(temp_arr)
                 aice_nhT   = calcMonAnomTLL(aice_nhT,climo) 
                 delete(climo)
              end if
           end if     
           aice_nh_seas = runave_n_Wrap(aice_nhT,3,0,0)
           aice_nh_seas(0,:,:) = (/ dim_avg_n(aice_nhT(:1,:,:),0) /)
           aice_nh_seas(dimsizes(aice_nhT&time)-1,:,:) = (/ dim_avg_n(aice_nhT(dimsizes(aice_nhT&time)-2:,:,:),0) /)
           aice_nh_ann = runave_n_Wrap(aice_nhT,12,0,0)
           delete(aice_nhT)
        
           if (ff.eq.0) then
              aice_nh_mean_djf = dim_avg_n_Wrap(aice_nh_seas(0::12,:,:),0)
              aice_nh_mean_mam = dim_avg_n_Wrap(aice_nh_seas(3::12,:,:),0)
              aice_nh_mean_jja = dim_avg_n_Wrap(aice_nh_seas(6::12,:,:),0)
              aice_nh_mean_son = dim_avg_n_Wrap(aice_nh_seas(9::12,:,:),0)
              aice_nh_mean_ann = dim_avg_n_Wrap(aice_nh_ann(5::12,:,:),0)
           end if
           if (ff.eq.1) then     
              aice_nh_sd_djf = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),aice_nh_seas(0::12,:,:),False,False,0),0)
              aice_nh_sd_mam = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),aice_nh_seas(3::12,:,:),False,False,0),0)
              aice_nh_sd_jja = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),aice_nh_seas(6::12,:,:),False,False,0),0)
              aice_nh_sd_son = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),aice_nh_seas(9::12,:,:),False,False,0),0)
              aice_nh_sd_ann = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr(ee)-1,1),aice_nh_ann(5::12,:,:),False,False,0),0)
           end if
           delete([/aice_nh_seas,aice_nh_ann/])
        end do
        delete(aice_nh)
        copy_VarMeta(aice_nh_mean_djf,aice_nh_sd_djf)
        copy_VarMeta(aice_nh_mean_mam,aice_nh_sd_mam)
        copy_VarMeta(aice_nh_mean_jja,aice_nh_sd_jja)
        copy_VarMeta(aice_nh_mean_son,aice_nh_sd_son)
        copy_VarMeta(aice_nh_mean_ann,aice_nh_sd_ann)
     end if

     aice_sh_flag = 0
     aice_sh = data_read_in_ice(paths_sh(ee),"aice_sh",syear_sh(ee),eyear_sh(ee))    ; read in data, orient lats/lons correctly, set time coordinate variable up
     if (isatt(aice_sh,"is_all_missing")) then
        delete(aice_sh)
        aice_sh_flag = 1
     end if  
     if (aice_sh_flag.eq.0) then
        do ff = 0,1     
           aice_shX = aice_sh
           if (ff.eq.1) then
              if (OPT_CLIMO.eq."Full") then
                 aice_shX = rmMonAnnCycTLL(aice_shX)
              else
                 check_custom_climo(names_sh(ee),syear_sh(ee),eyear_sh(ee),CLIMO_SYEAR,CLIMO_EYEAR)
                 temp_arr = aice_shX
                 delete(temp_arr&time)
                 temp_arr&time = cd_calendar(aice_shX&time,-1)
                 if (CLIMO_SYEAR.lt.0) then
                    climo = clmMonTLL(temp_arr({(eyear(ee)+CLIMO_SYEAR)*100+1:(eyear(ee)+CLIMO_EYEAR)*100+12},:,:))       
                 else
                    climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
                 end if               
                 delete(temp_arr)
                 aice_shX   = calcMonAnomTLL(aice_shX,climo) 
                 delete(climo)
              end if
           end if    
           aice_sh_seas = runave_n_Wrap(aice_shX,3,0,0)
           aice_sh_seas(0,:,:) = (/ dim_avg_n(aice_shX(:1,:,:),0) /)
           aice_sh_seas(dimsizes(aice_shX&time)-1,:,:) = (/ dim_avg_n(aice_shX(dimsizes(aice_shX&time)-2:,:,:),0) /)
           aice_sh_ann = runave_n_Wrap(aice_shX,12,0,0)
           delete(aice_shX)
        
           if (ff.eq.0) then
              aice_sh_mean_djf = dim_avg_n_Wrap(aice_sh_seas(0::12,:,:),0)
              aice_sh_mean_mam = dim_avg_n_Wrap(aice_sh_seas(3::12,:,:),0)
              aice_sh_mean_jja = dim_avg_n_Wrap(aice_sh_seas(6::12,:,:),0)
              aice_sh_mean_son = dim_avg_n_Wrap(aice_sh_seas(9::12,:,:),0)
              aice_sh_mean_ann = dim_avg_n_Wrap(aice_sh_ann(5::12,:,:),0)
           end if
           if (ff.eq.1) then     
              aice_sh_sd_djf = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr_sh(ee)-1,1),aice_sh_seas(0::12,:,:),False,False,0),0)
              aice_sh_sd_mam = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr_sh(ee)-1,1),aice_sh_seas(3::12,:,:),False,False,0),0)
              aice_sh_sd_jja = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr_sh(ee)-1,1),aice_sh_seas(6::12,:,:),False,False,0),0)
              aice_sh_sd_son = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr_sh(ee)-1,1),aice_sh_seas(9::12,:,:),False,False,0),0)
              aice_sh_sd_ann = dim_stddev_n_Wrap(dtrend_msg_n(ispan(0,nyr_sh(ee)-1,1),aice_sh_ann(5::12,:,:),False,False,0),0)
           end if
           delete([/aice_sh_seas,aice_sh_ann/])
        end do
        delete(aice_sh)
        copy_VarMeta(aice_sh_mean_djf,aice_sh_sd_djf)
        copy_VarMeta(aice_sh_mean_mam,aice_sh_sd_mam)
        copy_VarMeta(aice_sh_mean_jja,aice_sh_sd_jja)
        copy_VarMeta(aice_sh_mean_son,aice_sh_sd_son)
        copy_VarMeta(aice_sh_mean_ann,aice_sh_sd_ann)
     end if

     if (OUTPUT_DATA.eq."True".and.aice_nh_flag.eq.0) then
        modname = str_sub_str(names(ee)," ","_")
        bc = (/"/","'","(",")"/)
        do gg = 0,dimsizes(bc)-1
           modname = str_sub_str(modname,bc(gg),"_")
        end do
        fn = getenv("OUTDIR")+modname+".cvdp_data.aice.mean_stddev.nh."+syear(ee)+"-"+eyear(ee)+".nc"
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
        nh_mean_djf = aice_nh_mean_djf
        if (isatt(nh_mean_djf,"lat2d")) then    ; if there is a lat2d there will be a lon2d
           delete(nh_mean_djf@lat2d)
           delete(nh_mean_djf@lon2d)
           LAT2D = aice_nh_mean_djf@lat2d
           copy_VarCoords(nh_mean_djf,LAT2D)
           LON2D = aice_nh_mean_djf@lon2d
           copy_VarCoords(nh_mean_djf,LON2D)
           z->lat2d_ice_nh = set_varAtts(LAT2D,"Northern Hemisphere ice grid 2-dimensional latitudes","","")
           z->lon2d_ice_nh = set_varAtts(LON2D,"Northern Hemisphere ice grid 2-dimensional longitudes","","")
           delete([/LAT2D,LON2D/])
           nh_mean_djf@coordinates ="lat2d_ice_nh lon2d_ice_nh"
        end if  
        if (isatt(nh_mean_djf,"area")) then
           delete(nh_mean_djf@area)
        end if
        nh_mean_djf@long_name = nh_mean_djf@long_name+" mean"
        nh_mean_mam = (/ aice_nh_mean_mam /)
        copy_VarMeta(nh_mean_djf,nh_mean_mam)
        nh_mean_jja = (/ aice_nh_mean_jja /)
        copy_VarMeta(nh_mean_djf,nh_mean_jja)
        nh_mean_son = (/ aice_nh_mean_son /)
        copy_VarMeta(nh_mean_djf,nh_mean_son)
        nh_mean_ann = (/ aice_nh_mean_ann /)
        copy_VarMeta(nh_mean_djf,nh_mean_ann)
        nh_sd_djf = aice_nh_sd_djf
        if (isatt(nh_sd_djf,"lat2d")) then
           delete(nh_sd_djf@lat2d)
           delete(nh_sd_djf@lon2d)
           nh_sd_djf@coordinates ="lat2d_ice_nh lon2d_ice_nh"
        end if
        if (isatt(nh_sd_djf,"area")) then
           delete(nh_sd_djf@area)
        end if
        nh_sd_djf@long_name = nh_sd_djf@long_name+" standard deviation"
        nh_sd_mam = (/ aice_nh_sd_mam /)
        copy_VarMeta(nh_sd_djf,nh_sd_mam)
        nh_sd_jja = (/ aice_nh_sd_jja /)
        copy_VarMeta(nh_sd_djf,nh_sd_jja)
        nh_sd_son = (/ aice_nh_sd_son /)
        copy_VarMeta(nh_sd_djf,nh_sd_son)
        nh_sd_ann = (/ aice_nh_sd_ann /)
        copy_VarMeta(nh_sd_djf,nh_sd_ann)
        z->sic_nh_spatialmean_djf     = set_varAtts(nh_mean_djf,"Northern Hemisphere sic mean (DJF)","","")
        z->sic_nh_spatialmean_mam     = set_varAtts(nh_mean_mam,"Northern Hemisphere sic mean (MAM)","","")
        z->sic_nh_spatialmean_jja     = set_varAtts(nh_mean_jja,"Northern Hemisphere sic mean (JJA)","","")
        z->sic_nh_spatialmean_son     = set_varAtts(nh_mean_son,"Northern Hemisphere sic mean (SON)","","")
        z->sic_nh_spatialmean_ann     = set_varAtts(nh_mean_ann,"Northern Hemisphere sic mean (annual)","","")
        
        z->sic_nh_spatialstddev_djf     = set_varAtts(nh_sd_djf,"Northern Hemisphere sic standard deviation (DJF)","","")
        z->sic_nh_spatialstddev_mam     = set_varAtts(nh_sd_mam,"Northern Hemisphere sic standard deviation (MAM)","","")
        z->sic_nh_spatialstddev_jja     = set_varAtts(nh_sd_jja,"Northern Hemisphere sic standard deviation (JJA)","","")
        z->sic_nh_spatialstddev_son     = set_varAtts(nh_sd_son,"Northern Hemisphere sic standard deviation (SON)","","")
        z->sic_nh_spatialstddev_ann     = set_varAtts(nh_sd_ann,"Northern Hemisphere sic standard deviation (annual)","","")
        delete([/nh_mean_djf,nh_mean_mam,nh_mean_jja,nh_mean_son,nh_mean_ann/])
        delete([/nh_sd_djf,nh_sd_mam,nh_sd_jja,nh_sd_son,nh_sd_ann/])
        delete(z)
     end if
     if (OUTPUT_DATA.eq."True".and.aice_sh_flag.eq.0) then
        modname = str_sub_str(names_sh(ee)," ","_")
        bc = (/"/","'","(",")"/)
        do gg = 0,dimsizes(bc)-1
           modname = str_sub_str(modname,bc(gg),"_")
        end do
        fn = getenv("OUTDIR")+modname+".cvdp_data.aice.mean_stddev.sh."+syear_sh(ee)+"-"+eyear_sh(ee)+".nc"
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
        sh_mean_djf = aice_sh_mean_djf
        if (isatt(sh_mean_djf,"lat2d")) then    ; if there is a lat2d there will be a lon2d
           delete(sh_mean_djf@lat2d)
           delete(sh_mean_djf@lon2d)
           LAT2D = aice_sh_mean_djf@lat2d
           copy_VarCoords(sh_mean_djf,LAT2D)
           LON2D = aice_sh_mean_djf@lon2d
           copy_VarCoords(sh_mean_djf,LON2D)
           z->lat2d_ice_sh = set_varAtts(LAT2D,"Southern Hemisphere ice grid 2-dimensional latitudes","","")
           z->lon2d_ice_sh = set_varAtts(LON2D,"Southern Hemisphere ice grid 2-dimensional longitudes","","")
           delete([/LAT2D,LON2D/])
           sh_mean_djf@coordinates ="lat2d_ice_sh lon2d_ice_sh"
        end if  
        if (isatt(sh_mean_djf,"area")) then
           delete(sh_mean_djf@area)
        end if
        sh_mean_djf@long_name = sh_mean_djf@long_name+" mean"
        sh_mean_djf!0 = "j2"
        sh_mean_djf!1 = "i2"
        sh_mean_mam = (/ aice_sh_mean_mam /)
        copy_VarMeta(sh_mean_djf,sh_mean_mam)
        sh_mean_jja = (/ aice_sh_mean_jja /)
        copy_VarMeta(sh_mean_djf,sh_mean_jja)
        sh_mean_son = (/ aice_sh_mean_son /)
        copy_VarMeta(sh_mean_djf,sh_mean_son)
        sh_mean_ann = (/ aice_sh_mean_ann /)
        copy_VarMeta(sh_mean_djf,sh_mean_ann)
        sh_sd_djf = aice_sh_sd_djf
        if (isatt(sh_sd_djf,"lat2d")) then
           delete(sh_sd_djf@lat2d)
           delete(sh_sd_djf@lon2d)
           sh_sd_djf@coordinates ="lat2d_ice_sh lon2d_ice_sh"
        end if
        if (isatt(sh_sd_djf,"area")) then
           delete(sh_sd_djf@area)
        end if
        sh_sd_djf@long_name = sh_sd_djf@long_name+" standard deviation"
        sh_sd_djf!0 = "j2"
        sh_sd_djf!1 = "i2"
        sh_sd_mam = (/ aice_sh_sd_mam /)
        copy_VarMeta(sh_sd_djf,sh_sd_mam)
        sh_sd_jja = (/ aice_sh_sd_jja /)
        copy_VarMeta(sh_sd_djf,sh_sd_jja)
        sh_sd_son = (/ aice_sh_sd_son /)
        copy_VarMeta(sh_sd_djf,sh_sd_son)
        sh_sd_ann = (/ aice_sh_sd_ann /)
        copy_VarMeta(sh_sd_djf,sh_sd_ann)
        z->sic_sh_spatialmean_djf     = set_varAtts(sh_mean_djf,"Southern Hemisphere sic mean (DJF)","","")
        z->sic_sh_spatialmean_mam     = set_varAtts(sh_mean_mam,"Southern Hemisphere sic mean (MAM)","","")
        z->sic_sh_spatialmean_jja     = set_varAtts(sh_mean_jja,"Southern Hemisphere sic mean (JJA)","","")
        z->sic_sh_spatialmean_son     = set_varAtts(sh_mean_son,"Southern Hemisphere sic mean (SON","","")
        z->sic_sh_spatialmean_ann     = set_varAtts(sh_mean_ann,"Southern Hemisphere sic mean (annual)","","")
        
        z->sic_sh_spatialstddev_djf     = set_varAtts(sh_sd_djf,"Southern Hemisphere sic standard deviation (DJF)","","")
        z->sic_sh_spatialstddev_mam     = set_varAtts(sh_sd_mam,"Southern Hemisphere sic standard deviation (MAM)","","")
        z->sic_sh_spatialstddev_jja     = set_varAtts(sh_sd_jja,"Southern Hemisphere sic standard deviation (JJA)","","")
        z->sic_sh_spatialstddev_son     = set_varAtts(sh_sd_son,"Southern Hemisphere sic standard deviation (SON)","","")
        z->sic_sh_spatialstddev_ann     = set_varAtts(sh_sd_ann,"Southern Hemisphere sic standard deviation (annual)","","")
        delete([/sh_mean_djf,sh_mean_mam,sh_mean_jja,sh_mean_son,sh_mean_ann/])
        delete([/sh_sd_djf,sh_sd_mam,sh_sd_jja,sh_sd_son,sh_sd_ann/])
        delete(z)
     end if
     if (aice_nh_flag.eq.0) then
        aice_nh_mean_djf = where(aice_nh_mean_djf.lt.1,aice_nh_mean_djf@_FillValue,aice_nh_mean_djf)
        aice_nh_mean_mam = where(aice_nh_mean_mam.lt.1,aice_nh_mean_mam@_FillValue,aice_nh_mean_mam)
        aice_nh_mean_jja = where(aice_nh_mean_jja.lt.1,aice_nh_mean_jja@_FillValue,aice_nh_mean_jja)
        aice_nh_mean_son = where(aice_nh_mean_son.lt.1,aice_nh_mean_son@_FillValue,aice_nh_mean_son)
        aice_nh_mean_ann = where(aice_nh_mean_ann.lt.1,aice_nh_mean_ann@_FillValue,aice_nh_mean_ann)
        aice_nh_sd_djf = where(aice_nh_sd_djf.eq.0,aice_nh_sd_djf@_FillValue,aice_nh_sd_djf)
        aice_nh_sd_mam = where(aice_nh_sd_mam.eq.0,aice_nh_sd_mam@_FillValue,aice_nh_sd_mam)
        aice_nh_sd_jja = where(aice_nh_sd_jja.eq.0,aice_nh_sd_jja@_FillValue,aice_nh_sd_jja)
        aice_nh_sd_son = where(aice_nh_sd_son.eq.0,aice_nh_sd_son@_FillValue,aice_nh_sd_son)
        aice_nh_sd_ann = where(aice_nh_sd_ann.eq.0,aice_nh_sd_ann@_FillValue,aice_nh_sd_ann)
     end if
     if (aice_sh_flag.eq.0) then
        aice_sh_mean_djf = where(aice_sh_mean_djf.lt.1,aice_sh_mean_djf@_FillValue,aice_sh_mean_djf)
        aice_sh_mean_mam = where(aice_sh_mean_mam.lt.1,aice_sh_mean_mam@_FillValue,aice_sh_mean_mam)
        aice_sh_mean_jja = where(aice_sh_mean_jja.lt.1,aice_sh_mean_jja@_FillValue,aice_sh_mean_jja)
        aice_sh_mean_son = where(aice_sh_mean_son.lt.1,aice_sh_mean_son@_FillValue,aice_sh_mean_son)
        aice_sh_mean_ann = where(aice_sh_mean_ann.lt.1,aice_sh_mean_ann@_FillValue,aice_sh_mean_ann)
        aice_sh_sd_djf = where(aice_sh_sd_djf.eq.0,aice_sh_sd_djf@_FillValue,aice_sh_sd_djf)
        aice_sh_sd_mam = where(aice_sh_sd_mam.eq.0,aice_sh_sd_mam@_FillValue,aice_sh_sd_mam)
        aice_sh_sd_jja = where(aice_sh_sd_jja.eq.0,aice_sh_sd_jja@_FillValue,aice_sh_sd_jja)
        aice_sh_sd_son = where(aice_sh_sd_son.eq.0,aice_sh_sd_son@_FillValue,aice_sh_sd_son)
        aice_sh_sd_ann = where(aice_sh_sd_ann.eq.0,aice_sh_sd_ann@_FillValue,aice_sh_sd_ann)
     end if
;==========================================================================================
     res = True
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
     res@mpMinLatF    = 40.
     res@mpCenterLonF = 0.
     res@gsnPolar = "NH"     
     res@gsnDraw      = False
     res@gsnFrame     = False
     res@gsnAddCyclic = True
     res@cnLevelSelectionMode = "ExplicitLevels"
     res@cnLineLabelsOn = False
     res@cnFillOn        = True
     res@cnLinesOn       = False
     res@trGridType = "TriangularMesh"
;     res@cnFillMode = "RasterFill"
     res@lbLabelBarOn    = False

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

    
     res@cnLevelSelectionMode = "ExplicitLevels"

     sres = res
     res@cnLevels = (/4,8,12,16,20,24,28,32,36/)
     sres@cnLevels = (/5,10,15,20,30,40,50,60,70,80,85,90,95,99/)
     contour_means = sres@cnLevels   ; for use in paneling section
     contour_sd = res@cnLevels
     if (COLORMAP.eq.0) then
        res@cnFillColors =  (/42,29,80,95,105,120,140,161,170,193/) ; radar: (/5,6,7,8,9,11,12,13,14,15/)
        sres@cnFillColors = (/52,42,34,24,65,80,95,105,120,140,155,161,170,184,193/)        
     end if
     if (COLORMAP.eq.1) then
        res@cnFillColors = (/2,18,34,50,66,82,98,114,137,162/)
        sres@cnFillColors = (/8,26,38,50,62,74,86,98,110,122,134,146,158,170,182/)        
     end if
   
     if (aice_nh_flag.eq.0) then
        res@gsnLeftString = syear(ee)+"-"+eyear(ee) 
        res@gsnRightString = aice_nh_sd_djf@units
        res@gsnCenterString = names(ee)  
        plot_stddev_nh_djf(ee) = gsn_csm_contour_map(wks_stddev_djf,aice_nh_sd_djf,res)
        plot_stddev_nh_mam(ee) = gsn_csm_contour_map(wks_stddev_mam,aice_nh_sd_mam,res)
        plot_stddev_nh_jja(ee) = gsn_csm_contour_map(wks_stddev_jja,aice_nh_sd_jja,res)
        plot_stddev_nh_son(ee) = gsn_csm_contour_map(wks_stddev_son,aice_nh_sd_son,res)
        plot_stddev_nh_ann(ee) = gsn_csm_contour_map(wks_stddev_ann,aice_nh_sd_ann,res)

        sres@gsnLeftString = syear(ee)+"-"+eyear(ee) 
        sres@gsnRightString = aice_nh_mean_djf@units
        sres@gsnCenterString = names(ee)
        plot_mean_nh_djf(ee) = gsn_csm_contour_map(wks_mean_djf,aice_nh_mean_djf,sres)
        plot_mean_nh_mam(ee) = gsn_csm_contour_map(wks_mean_mam,aice_nh_mean_mam,sres)
        plot_mean_nh_jja(ee) = gsn_csm_contour_map(wks_mean_jja,aice_nh_mean_jja,sres)
        plot_mean_nh_son(ee) = gsn_csm_contour_map(wks_mean_son,aice_nh_mean_son,sres)
        plot_mean_nh_ann(ee) = gsn_csm_contour_map(wks_mean_ann,aice_nh_mean_ann,sres)
        delete([/aice_nh_sd_djf,aice_nh_sd_mam,aice_nh_sd_jja,aice_nh_sd_son,aice_nh_sd_ann/])
        delete([/aice_nh_mean_djf,aice_nh_mean_mam,aice_nh_mean_jja,aice_nh_mean_son,aice_nh_mean_ann/])
     end if
     
     delete(res@mpMinLatF)
     delete(sres@mpMinLatF)
     res@mpMaxLatF = -45.
     res@gsnPolar = "SH"  
     sres@mpMaxLatF = -45.
     sres@gsnPolar = "SH"  
  
     if (aice_sh_flag.eq.0) then
        res@gsnLeftString = syear_sh(ee)+"-"+eyear_sh(ee) 
        res@gsnRightString = aice_sh_sd_djf@units
        res@gsnCenterString = names_sh(ee)
        plot_stddev_sh_djf(ee) = gsn_csm_contour_map(wks_stddev_djf,aice_sh_sd_djf,res)
        plot_stddev_sh_mam(ee) = gsn_csm_contour_map(wks_stddev_mam,aice_sh_sd_mam,res)
        plot_stddev_sh_jja(ee) = gsn_csm_contour_map(wks_stddev_jja,aice_sh_sd_jja,res)
        plot_stddev_sh_son(ee) = gsn_csm_contour_map(wks_stddev_son,aice_sh_sd_son,res)
        plot_stddev_sh_ann(ee) = gsn_csm_contour_map(wks_stddev_ann,aice_sh_sd_ann,res)

        sres@gsnLeftString = syear_sh(ee)+"-"+eyear_sh(ee) 
        sres@gsnRightString = aice_sh_mean_djf@units
        sres@gsnCenterString = names_sh(ee)
        plot_mean_sh_djf(ee) = gsn_csm_contour_map(wks_mean_djf,aice_sh_mean_djf,sres)
        plot_mean_sh_mam(ee) = gsn_csm_contour_map(wks_mean_mam,aice_sh_mean_mam,sres)
        plot_mean_sh_jja(ee) = gsn_csm_contour_map(wks_mean_jja,aice_sh_mean_jja,sres)
        plot_mean_sh_son(ee) = gsn_csm_contour_map(wks_mean_son,aice_sh_mean_son,sres)
        plot_mean_sh_ann(ee) = gsn_csm_contour_map(wks_mean_ann,aice_sh_mean_ann,sres)
        delete([/aice_sh_sd_djf,aice_sh_sd_mam,aice_sh_sd_jja,aice_sh_sd_son,aice_sh_sd_ann/])
        delete([/aice_sh_mean_djf,aice_sh_mean_mam,aice_sh_mean_jja,aice_sh_mean_son,aice_sh_mean_ann/])
     end if
     delete([/res,sres/])
  end do     
     
  panres = True
  panres@gsnMaximize = True  
  panres@gsnPaperOrientation = "portrait"
  panres@gsnPanelLabelBar = True
  panres@gsnPanelYWhiteSpacePercent = 3.0
  panres@pmLabelBarHeightF = 0.05
  panres@pmLabelBarWidthF = 0.65
  panres@lbTitleOn = False
  panres@lbBoxLineColor = "gray70"
  if (nsim.le.4) then
     if (nsim.eq.1) then
        panres@txFontHeightF = 0.022
        panres@gsnPanelBottom = 0.50
     else
        panres@txFontHeightF = 0.0145
        panres@gsnPanelBottom = 0.50
     end if
  else
     panres@txFontHeightF = 0.016
     panres@gsnPanelBottom = 0.05
  end if
  panres@lbLabelFontHeightF = 0.013
  panres@lbLabelStride = 1
  ncol = floattointeger(sqrt(nsim))
  nrow = (nsim/ncol)+mod(nsim,ncol)  

  c_sd = (/4,8,12,16,20,24,28,32,36/)
  if (dimsizes(c_sd).eq.dimsizes(contour_sd).and.all((c_sd - contour_sd).eq.0)) then  ; needed to make sure contour intervals
     panres@lbLabelAlignment = "ExternalEdges"                               ; set above match the labels set here
     panres@lbLabelStrings = (/"0","4","8","12","16","20","24","28","32","36"," "/)
  end if
  delete([/c_sd,contour_sd/])  

  panres@txString = "SIC Standard Deviations (DJF)"
  gsn_panel2(wks_stddev_djf,plot_stddev_nh_djf,(/nrow,ncol/),panres)
  gsn_panel2(wks_stddev_djf,plot_stddev_sh_djf,(/nrow,ncol/),panres)
  delete(wks_stddev_djf)
  
  panres@txString = "SIC Standard Deviations (MAM)"
  gsn_panel2(wks_stddev_mam,plot_stddev_nh_mam,(/nrow,ncol/),panres)
  gsn_panel2(wks_stddev_mam,plot_stddev_sh_mam,(/nrow,ncol/),panres)
  delete(wks_stddev_mam)
  
  panres@txString = "SIC Standard Deviations (JJA)"
  gsn_panel2(wks_stddev_jja,plot_stddev_nh_jja,(/nrow,ncol/),panres)
  gsn_panel2(wks_stddev_jja,plot_stddev_sh_jja,(/nrow,ncol/),panres)
  delete(wks_stddev_jja)

  panres@txString = "SIC Standard Deviations (SON)"
  gsn_panel2(wks_stddev_son,plot_stddev_nh_son,(/nrow,ncol/),panres)
  gsn_panel2(wks_stddev_son,plot_stddev_sh_son,(/nrow,ncol/),panres)
  delete(wks_stddev_son)

  panres@txString = "SIC Standard Deviations (Annual)"
  gsn_panel2(wks_stddev_ann,plot_stddev_nh_ann,(/nrow,ncol/),panres)
  gsn_panel2(wks_stddev_ann,plot_stddev_sh_ann,(/nrow,ncol/),panres)
  delete(wks_stddev_ann)

  if (isatt(panres,"lbLabelAlignment")) then
     delete(panres@lbLabelAlignment)
     delete(panres@lbLabelStrings)
  end if

  c_me = (/5,10,15,20,30,40,50,60,70,80,85,90,95,99/)
  if (dimsizes(c_me).eq.dimsizes(contour_means).and.all((c_me - contour_means).eq.0)) then  ; needed to make sure contour intervals
     panres@lbLabelAlignment = "ExternalEdges"                              ; set above match the labels set here
     panres@lbLabelStrings = (/"1","5","10","15","20","30","40","50","60","70","80","85","90","95","99"," "/)
  end if
  delete([/c_me,contour_means/])

  panres@txString = "SIC Means (DJF)"
  gsn_panel2(wks_mean_djf,plot_mean_nh_djf,(/nrow,ncol/),panres)
  gsn_panel2(wks_mean_djf,plot_mean_sh_djf,(/nrow,ncol/),panres)
  delete(wks_mean_djf)
  
  panres@txString = "SIC Means (MAM)"
  gsn_panel2(wks_mean_mam,plot_mean_nh_mam,(/nrow,ncol/),panres)
  gsn_panel2(wks_mean_mam,plot_mean_sh_mam,(/nrow,ncol/),panres)
  delete(wks_mean_mam)
  
  panres@txString = "SIC Means (JJA)"
  gsn_panel2(wks_mean_jja,plot_mean_nh_jja,(/nrow,ncol/),panres)
  gsn_panel2(wks_mean_jja,plot_mean_sh_jja,(/nrow,ncol/),panres)
  delete(wks_mean_jja)

  panres@txString = "SIC Means (SON)"
  gsn_panel2(wks_mean_son,plot_mean_nh_son,(/nrow,ncol/),panres)
  gsn_panel2(wks_mean_son,plot_mean_sh_son,(/nrow,ncol/),panres)
  delete(wks_mean_son)

  panres@txString = "SIC Means (Annual)"
  gsn_panel2(wks_mean_ann,plot_mean_nh_ann,(/nrow,ncol/),panres)
  gsn_panel2(wks_mean_ann,plot_mean_sh_ann,(/nrow,ncol/),panres)
  delete(wks_mean_ann)
  delete(panres)
;--------------------------------------------------------------------------------
  OUTDIR = getenv("OUTDIR") 
  if (wks_type.eq."png") then  
     system("mv "+OUTDIR+"aice.stddev.djf.000001.png "+OUTDIR+"aice.stddev.nh.djf.png") 
     system("mv "+OUTDIR+"aice.stddev.djf.000002.png "+OUTDIR+"aice.stddev.sh.djf.png") 
     system("mv "+OUTDIR+"aice.stddev.mam.000001.png "+OUTDIR+"aice.stddev.nh.mam.png") 
     system("mv "+OUTDIR+"aice.stddev.mam.000002.png "+OUTDIR+"aice.stddev.sh.mam.png")
     system("mv "+OUTDIR+"aice.stddev.jja.000001.png "+OUTDIR+"aice.stddev.nh.jja.png") 
     system("mv "+OUTDIR+"aice.stddev.jja.000002.png "+OUTDIR+"aice.stddev.sh.jja.png")
     system("mv "+OUTDIR+"aice.stddev.son.000001.png "+OUTDIR+"aice.stddev.nh.son.png") 
     system("mv "+OUTDIR+"aice.stddev.son.000002.png "+OUTDIR+"aice.stddev.sh.son.png") 
     system("mv "+OUTDIR+"aice.stddev.ann.000001.png "+OUTDIR+"aice.stddev.nh.ann.png") 
     system("mv "+OUTDIR+"aice.stddev.ann.000002.png "+OUTDIR+"aice.stddev.sh.ann.png")

     system("mv "+OUTDIR+"aice.mean.djf.000001.png "+OUTDIR+"aice.mean.nh.djf.png") 
     system("mv "+OUTDIR+"aice.mean.djf.000002.png "+OUTDIR+"aice.mean.sh.djf.png") 
     system("mv "+OUTDIR+"aice.mean.mam.000001.png "+OUTDIR+"aice.mean.nh.mam.png") 
     system("mv "+OUTDIR+"aice.mean.mam.000002.png "+OUTDIR+"aice.mean.sh.mam.png")
     system("mv "+OUTDIR+"aice.mean.jja.000001.png "+OUTDIR+"aice.mean.nh.jja.png") 
     system("mv "+OUTDIR+"aice.mean.jja.000002.png "+OUTDIR+"aice.mean.sh.jja.png")
     system("mv "+OUTDIR+"aice.mean.son.000001.png "+OUTDIR+"aice.mean.nh.son.png") 
     system("mv "+OUTDIR+"aice.mean.son.000002.png "+OUTDIR+"aice.mean.sh.son.png") 
     system("mv "+OUTDIR+"aice.mean.ann.000001.png "+OUTDIR+"aice.mean.nh.ann.png") 
     system("mv "+OUTDIR+"aice.mean.ann.000002.png "+OUTDIR+"aice.mean.sh.ann.png")  
  else
     system("psplit "+OUTDIR+"aice.stddev.djf.ps "+OUTDIR+"aice_sd")
     system("mv "+OUTDIR+"aice_sd0001.ps "+OUTDIR+"aice.stddev.nh.djf.ps") 
     system("mv "+OUTDIR+"aice_sd0002.ps "+OUTDIR+"aice.stddev.sh.djf.ps")
     system("psplit "+OUTDIR+"aice.stddev.mam.ps "+OUTDIR+"aice_sd")
     system("mv "+OUTDIR+"aice_sd0001.ps "+OUTDIR+"aice.stddev.nh.mam.ps") 
     system("mv "+OUTDIR+"aice_sd0002.ps "+OUTDIR+"aice.stddev.sh.mam.ps")
     system("psplit "+OUTDIR+"aice.stddev.jja.ps "+OUTDIR+"aice_sd")
     system("mv "+OUTDIR+"aice_sd0001.ps "+OUTDIR+"aice.stddev.nh.jja.ps") 
     system("mv "+OUTDIR+"aice_sd0002.ps "+OUTDIR+"aice.stddev.sh.jja.ps")
     system("psplit "+OUTDIR+"aice.stddev.son.ps "+OUTDIR+"aice_sd")
     system("mv "+OUTDIR+"aice_sd0001.ps "+OUTDIR+"aice.stddev.nh.son.ps") 
     system("mv "+OUTDIR+"aice_sd0002.ps "+OUTDIR+"aice.stddev.sh.son.ps")
     system("psplit "+OUTDIR+"aice.stddev.ann.ps "+OUTDIR+"aice_sd")
     system("mv "+OUTDIR+"aice_sd0001.ps "+OUTDIR+"aice.stddev.nh.ann.ps") 
     system("mv "+OUTDIR+"aice_sd0002.ps "+OUTDIR+"aice.stddev.sh.ann.ps")
     system("rm "+OUTDIR+"aice.stddev.???.ps")

     system("psplit "+OUTDIR+"aice.mean.djf.ps "+OUTDIR+"aice_m")
     system("mv "+OUTDIR+"aice_m0001.ps "+OUTDIR+"aice.mean.nh.djf.ps") 
     system("mv "+OUTDIR+"aice_m0002.ps "+OUTDIR+"aice.mean.sh.djf.ps")
     system("psplit "+OUTDIR+"aice.mean.mam.ps "+OUTDIR+"aice_m")
     system("mv "+OUTDIR+"aice_m0001.ps "+OUTDIR+"aice.mean.nh.mam.ps") 
     system("mv "+OUTDIR+"aice_m0002.ps "+OUTDIR+"aice.mean.sh.mam.ps")
     system("psplit "+OUTDIR+"aice.mean.jja.ps "+OUTDIR+"aice_m")
     system("mv "+OUTDIR+"aice_m0001.ps "+OUTDIR+"aice.mean.nh.jja.ps") 
     system("mv "+OUTDIR+"aice_m0002.ps "+OUTDIR+"aice.mean.sh.jja.ps")
     system("psplit "+OUTDIR+"aice.mean.son.ps "+OUTDIR+"aice_m")
     system("mv "+OUTDIR+"aice_m0001.ps "+OUTDIR+"aice.mean.nh.son.ps") 
     system("mv "+OUTDIR+"aice_m0002.ps "+OUTDIR+"aice.mean.sh.son.ps")
     system("psplit "+OUTDIR+"aice.mean.ann.ps "+OUTDIR+"aice_m")
     system("mv "+OUTDIR+"aice_m0001.ps "+OUTDIR+"aice.mean.nh.ann.ps") 
     system("mv "+OUTDIR+"aice_m0002.ps "+OUTDIR+"aice.mean.sh.ann.ps")
     system("rm "+OUTDIR+"aice.mean.???.ps")
  end if
  print("Finished: aice.mean_stddev.ncl")
end
