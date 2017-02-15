load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messagging.ncl"

load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/calc_seasonal.ncl"
load "./diag_scripts/uniform_coordinates.ncl"

begin
 print("1") 
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output("plot SET8 - " + diag_script + ".ncl (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    var1 = variables(1)
    var2=  variables(2)
    print(var0)
    print(var1)
    print(var2)

    field_type0 = field_types(0)
  
 
    do idx_mod = 0, dim_MOD - 1
        data0 = read_data(idx_mod, var0, field_type0)
        data1 = read_data(idx_mod, var1, field_type0)
        data2 = read_data(idx_mod, var2, field_type0)
    end do
    printVarSummary(data0)
    printVarSummary(data1)
end

begin
    vardeffile = "interface_data/" + "ps_info.tmp"
    
    loadscript(vardeffile)
end
begin


    plot_dir = getenv("ESMValTool_plot_dir")
    diag_script_base = basename(diag_script)
    output_dir = get_output_dir(plot_dir, diag_script_base)

    file_type = getenv("ESMValTool_output_file_type")
    if(ismissing(file_type)) then
        file_type = "ps"
    end if
    file_type = str_lower(file_type)

    color_type = "COLOR"

      
    ;diag_script_base = basename(diag_script)
    ;aux_title_info=str_join(,"_")
     
    ;output_filename = interface_get_figure_filename(diag_script_base,\
    ;                                            var0,\
    ;                                            "T2MS",\
    ;                                            aux_title_info,\
    ;                                            0)
    ; print(output_filename)
     
    MODELS=models@name
    dim_MOD = dimsizes(models@name)
    do imod = 0, dim_MOD-1
        ;print(MODELS(imod))
        output_filename= "toz_Amon_"+MODELS(imod)+"_historical_r1i1p1.nc"
        tro3_t = extract_data(imod, data0, -1, 0, 0)     
;        ps_t = extract_data(imod, data1, -1, 0, 0)
;        ta_t = extract_data(imod, data2, -1, 0, 0)
        ;uniform_coordinates(tro3_t, field_type0)
        ;uniform_coordinates(ps_t, field_type0)
        printVarSummary(tro3_t)
        time=tro3_t&time
        lat=tro3_t&lat
        lon=tro3_t&lon
       ; plev=tro3_t&plev
        tropoz=new((/dimsizes(time),dimsizes(lat),dimsizes(lon)/),typeof(tro3_t))
 delete(tro3_t)       
        
        do itim=0,dimsizes(time)-1
     
           tro3_t = extract_data(imod, data0, -1, 0, 0) 
           tro3=tro3_t(itim,:,:,:)
           delete(tro3_t)
           plev=tro3&plev    
          
            if (plev@units.ne."Pa") then
                plev=plev*100.
                plev@units="Pa"
            end if
         ;   print(plev)
            ps_t = extract_data(imod, data1, -1, 0, 0)    
            ps=ps_t(itim,:,:)
            delete(ps_t)
          ;   printVarSummary(ps) 
            
            ta_t= extract_data(imod, data2, -1, 0, 0)
            ta=ta_t(itim,:,:,:)
           delete(ta_t) 
           ;printVarSummary(ta)
            
            if (plev(0).gt.plev(1)) then
                plev1=plev(::-1)
                tro31=tro3(::-1,:,:)
                ta1=ta(::-1,:,:)
                delete(plev)
                delete(tro3)
                delete(ta)
                plev=plev1
                tro3=tro31
                ta=ta1
                delete(plev1)
                delete(tro31)
                delete(ta1)
            end if
            ptop=10.
         
;*********************************WMO tropopause***********************
         opt=True
         
         opt@lapsec=2.0
         ;print(plev@units)
         
         if(plev@units.eq."Pa") then
            unit_plev=1
         else
            unit_plev=0
         end if
         
         trop_h=trop_wmo(plev,ta(lat|:,lon|:,plev|:),unit_plev,opt)
         
        ; printVarSummary(trop_h)
        
         trop_h!0="lat"
         trop_h!1="lon"


         trop_h&lat=ta&lat
         trop_h&lon=ta&lon
         trop_h@_FillValue=10e20
         trop_h=where(trop_h.le.0,trop_h@_FillValue,trop_h)

         trop_hD=conform_dims((/dimsizes(plev),dimsizes(lat),dimsizes(lon)/),trop_h,(/1,2/))
         delete(trop_h)
         ;printVarSummary(trop_hD) 
         pm=conform_dims((/dimsizes(plev),dimsizes(lat),dimsizes(lon)/),plev,(/0/))
         delete(ta)
         
          
;********************************************************************
         dp   = dpres_plevel(plev, ps, ptop, 0)
         pm=conform_dims((/dimsizes(plev),dimsizes(lat),dimsizes(lon)/),plev,(/0/))
         delete(ps)
;*******************************************************
            gravit=9.81
            mw=  28.966/48.         ;peso mol air/peso mol ozone
      
            constant=1./(gravit*mw)
           
            dims=dimsizes(tro3)
        
            o3nd=tro3*(10^-9)*tofloat(dp) *  constant
          delete(tro3) 
          delete(dp)
            trop_o3col_du= where(pm.gt.trop_hD,o3nd,o3nd@_FillValue)
          delete(trop_hD)
          delete(o3nd)     
          ; printVarSummary(trop_o3col_du)     
          
       
            tropoz(itim,:,:)=(/tofloat(dim_sum_n(trop_o3col_du,0)/(2.1415*10^(-5)))/)
        
                   
            delete(trop_o3col_du)

        
            delete(pm)
       
         end do

        tropoz!0="time"
        tropoz!1="lat"
        tropoz!2="lon"
        tropoz&time=time
        tropoz&lat=lat
        tropoz&lon=lon
        tropoz@units="DU"
      ;----------------create file--------------------------------------
        system("/bin/rm -f " +output_dir+output_filename )          
        fout = addfile ( output_dir+output_filename,"c")
        setfileoption(fout,"DefineMode",True)

     ;===================================================================
    ; create global attributes of the file
    ;===================================================================
        fAtt               = True            ; assign file attributes

       ; fAtt=f_in@global

        ;global_attribute ( MODEL,table_id, project_id, frequency, experiment, experiment_id,fAtt)
        ;fileattdef( fout2, fAtt )            ; copy file attributes
    ;===================================================================
        dimNames= (/"time","lat","lon"/)
        dimSizes = (/ -1   ,dimsizes(lat), dimsizes(lon) /)
        dimUnlim = (/ True , False, False/)
        filedimdef(fout,dimNames,dimSizes,dimUnlim)

    ;===================================================================
        filevardef(fout, "lon"  ,typeof(lon),"lon")
        filevardef(fout, "lat"  ,typeof(lat),"lat")
        filevardef(fout, "time" ,typeof(time),"time")
        filevardef(fout, "tropoz",typeof(tropoz)  ,(/"time","lat","lon"/))

    ;===================================================================
    ; Copy attributes associated with each variable to the file
    ; All attributes associated with each variable will be copied.
    ;====================================================================
       filevarattdef(fout,"tropoz",tropoz)                           ; copy var attributes

       filevarattdef(fout,"time",time)                     ; copy time attributes
       filevarattdef(fout,"lat",lat)                       ; copy lon attributes
       filevarattdef(fout,"lon",lon)

       fout->lon    = (/lon/)
       fout->lat    = (/lat/)
       fout->time   = (/time/)

       fout->tropoz      = (/tropoz/)
;********************************3*****************************************

       delete(time)
       delete(lat)
       delete(lon)
       delete(plev)
    end do
end
