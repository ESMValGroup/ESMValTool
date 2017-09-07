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
;    print(var0)
;    print(var1)

    field_type0 = field_types(0)
  
 
    do idx_mod = 0, dim_MOD - 1
        data0 = read_data(idx_mod, var0, field_type0)
        data1 = read_data(idx_mod, var1, field_type0)
    end do
    ;printVarSummary(data0)
    ;printVarSummary(data1)
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
        print(MODELS(imod))
        output_filename= "toz_Amon_"+MODELS(imod)+"_historical_r1i1p1.nc"
        tro3_t = extract_data(imod, data0, -1, 0, 0)
        ps_t = extract_data(imod, data1, -1, 0, 0)
        ;uniform_coordinates(tro3_t, field_type0)
        ;uniform_coordinates(ps_t, field_type0)
        printVarSummary(tro3_t)
        time=tro3_t&time
        lat=tro3_t&lat
        lon=tro3_t&lon
       ; plev=tro3_t&plev
        toz=new((/dimsizes(time),dimsizes(lat),dimsizes(lon)/),typeof(tro3_t))
        
        
        do itim=0,dimsizes(time)-1
        ;    print(itim)
            tro3=tro3_t(itim,:,:,:)
         plev=tro3_t&plev    
   
            if (plev@units.ne."Pa") then
               plev=plev*100.
                plev@units="Pa"
            end if
         
            ps=ps_t(itim,:,:)
            
            if (plev(0).gt.plev(1)) then
                plev1=plev(::-1)
                tro31=tro3(::-1,:,:)
                delete(plev)
                delete(tro3)
                plev=plev1
                tro3=tro31
                delete(plev1)
                delete(tro31)
            end if
            ptop=10.
            dp   = dpres_plevel(plev, ps, ptop, 0)
         pm=conform_dims((/dimsizes(plev),dimsizes(lat),dimsizes(lon)/),plev,(/0/))

;*******************************************************
            gravit=9.81
            mw=  28.966/48.         ;peso mol air/peso mol ozone
      
            constant=1./(gravit*mw)
           
            dims=dimsizes(tro3)
        
            o3nd=tro3*(10^-9)*tofloat(dp) *  constant
     

             
           
       
            toz(itim,:,:)=(/tofloat(dim_sum_n(o3nd,0)/(2.1415*10^(-5)))/)
                    
          
            delete(o3nd)
            delete(tro3)
            delete(ps)
            delete(dp)
            delete(pm)
         end do
        
        toz!0="time"
        toz!1="lat"
        toz!2="lon"
        toz&time=time
        toz&lat=lat
        toz&lon=lon
        toz@units="DU"
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
        filevardef(fout, "toz",typeof(toz)  ,(/"time","lat","lon"/))

    ;===================================================================
    ; Copy attributes associated with each variable to the file
    ; All attributes associated with each variable will be copied.
    ;====================================================================
       filevarattdef(fout,"toz",toz)                           ; copy var attributes

       filevarattdef(fout,"time",time)                     ; copy time attributes
       filevarattdef(fout,"lat",lat)                       ; copy lon attributes
       filevarattdef(fout,"lon",lon)

       fout->lon    = (/lon/)
       fout->lat    = (/lat/)
       fout->time   = (/time/)

       fout->toz      = (/toz/)
;********************************3*****************************************

       delete(time)
       delete(lat)
       delete(lon)
       delete(plev)
       delete(ps_t)
       delete(tro3_t)
    end do
end
