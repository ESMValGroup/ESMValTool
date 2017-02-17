;;############################################################################
;; crbon_corr_2var.ncl
;; Author: Sabrina Wenzel (DLR, Germany)                                      
;; ESMVal project PI: Veronika Eyring (DLR, Germany)                          
;;############################################################################
;; Description
;;    Reads and correlates in time spatial averages of CO2 fluxes (derived 
;;    from land and ocean carbon fluxes) and temperature
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    * ts_latrange2 min and max latitudes for spatial average
;;    * ts_lonrange2 min and max longitudes for spatial average
;;    * ec_units     units of first var
;;    * ec_v2_units  units of second var 
;;    * ec_anom      (optional) plot anomalies. 
;;    * ec_detr      (optional) plot detrendet variables 
;;
;; Caveats
;;
;; Modification history
;;    * 20151027-A_laue_ax: added call to 'write_references'
;;    * 201404-A_Sabrina Wenzel: updated according to changes in ESMValTool
;;    * 201303-A_Sabrina Wenzel: written       
;;                                                                            
;;############################################################################
load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"
load "./interface_scripts/constants.ncl"

load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/uniform_coordinates.ncl"
load "./diag_scripts/lib/ncl/scaling.ncl"
load "./diag_scripts/aux/SAMonsoon/SAMonsoon.ncl"

load "./plot_scripts/ncl/scatterplot.ncl"
load "./plot_scripts/ncl/aux_plotting.ncl"

begin
  plot_script = "corr_2var"
  
  verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
  info_output("<<<<<<<< Entering " + plot_script + ".ncl", verbosity, 4)
  info_output( plot_script + ".ncl", verbosity, 1)

; A-laue_ax+
    ; ###########################################
    ; # references                              #
    ; ###########################################

    write_references(diag_script,      \  ; script name
                     "A_wenz_sa",      \  ; authors
                     "",               \  ; contributors
                     "D_wenzel14jgr",  \  ; diag_references
                     "",               \  ; obs_references
                     (/""/))              ; proj_references
; A-laue_ax-

  dim_MOD = dimsizes(models@name)

  var0 = variables(0)
  var1 = variables(1)

  data_1st = True
  data_2nd = True
  do idx_mod = 0, dim_MOD - 1
        ;; Skip models<->parameter combinations that are explicitly
        ;; excluded in the xml-file
      curr_var_idx = ind(var0 .eq. variables)
      if (var_attr_exclude(curr_var_idx) .ne. model_attr_id(idx_mod)) then
          data_1st = read_data(idx_mod, var0, field_types(0))
      end if
      curr_var_idx = ind(var1 .eq. variables)
      if (var_attr_exclude(curr_var_idx) .ne. model_attr_id(idx_mod)) then
          data_2nd = read_data(idx_mod, var1, field_types(1))
      end if
  end do
end

;; Load var_info
begin
    var0deffile = "interface_data/" + var0 + "_info.tmp"
    var1deffile = "interface_data/" + var1 + "_info.tmp"
    loadscript(var0deffile)
    loadscript(var1deffile)
end

;; Main
begin 
  ;;****************************************************************** 
  ;; Check required diag_script_info attributes
    req_atts = (/"ts_latrange2", "ts_lonrange2", "ec_timerange"/)
    if (.not.all(isatt(diag_script_info, req_atts))) then
        error_msg("f", diag_script, "", "missing required " + \
                  "diag_script_info attribute:" + \
                  req_atts(ind(.not.isatt(diag_script_info, req_atts))))
    end if

  if(isvar("latrange_v2")) then
    delete(latrange_v2)
  end if
  if(isatt(diag_script_info, "ts_latrange2")) then
    latrange_v2 = diag_script_info@ts_latrange2
  end if
  
  if(isvar("lonrange_v2")) then
    delete(lonrange_v2)
  end if
  if(isatt(diag_script_info, "ts_lonrange2")) then
    lonrange_v2 = diag_script_info@ts_lonrange2
  end if

  if(isvar("ec_timerange")) then
    delete(ec_timerange)
  end if
  if(isatt(diag_script_info, "ec_timerange")) then
    timerange = diag_script_info@ec_timerange
  end if

  if(isatt(diag_script_info, "ec_units")) then
    UNITS = diag_script_info@ec_units
  end if
  if(isatt(diag_script_info, "ec_v2_units")) then
    UNITS_v2 = diag_script_info@ec_v2_units
  end if
    
  anom = False 
  if(isatt(diag_script_info, "ec_anom")) then
    anom = diag_script_info@ts_anom
  end if
  
  detr = False
  if(isatt(diag_script_info,"ec_detr")) then 
    detr=diag_script_info@ec_detr 
  end if
  
  ;;================================================================== 
  ;; Reference model index
  iref_co2 = get_ref_model_idx(models, var_attr_ref(0))
  iref_tas = get_ref_model_idx(models, var_attr_ref(1))
  iref = iref_tas - iref_co2
  min_idx = min((/iref_co2, iref_tas/))
  max_idx = max((/iref_co2, iref_tas/))
  iref = where(iref_co2 .eq. iref_tas, 0, 1)

  ;;time span
  xMax_year = max(toint(models@end_year))
  xMin_year = min(toint(models@start_year))
  nMax = xMax_year - xMin_year + 1

  ;;setup data holders for time lines 
  aY1reg = new((/dim_MOD-iref,nMax/),"float") 
  aY2reg = new((/dim_MOD-iref,nMax/),"float") 
  aY1reg!1 = "year"
  aY1reg&year = ispan(xMin_year,xMax_year,1)  
  aY2reg!1 = "year"
  aY2reg&year = ispan(xMin_year,xMax_year,1) 

  ;;setup data holders for slope 
  aYreg  = new((/dim_MOD-iref,nMax/),"float") 
  arc    = new((/2,dim_MOD-iref/),"float") 
  ccr    = new((/dim_MOD-iref/),"float")
  var    = new((/2,dim_MOD-iref/),"float")
  ;;================================================================== 
  ;;loop for models
  do imod = 0, dim_MOD-1 
      info_output(models@name(imod), verbosity, 1) 
    
      ;; Check for reference data sets
      if (any(imod .eq. (/iref_co2, iref_tas/))) then

          ;; Check if there's a single reference data set
          if (iref_co2 .eq. iref_tas) then
              A0 := extract_data_extended(imod, var0, data_1st, -1, 0, 0)
              A2 := extract_data_extended(imod, var1, data_2nd, -1, 0, 0)  

          ;; If there are two separate data sets, only use them once
          ;; Skip the second time (max_idx)
          else if (imod .eq. max_idx) then
              continue
          else  ; Use them the first time (min_idx)
              A0 := extract_data_extended(iref_co2, var0, data_1st, -1, 0, 0)
              A2 := extract_data_extended(iref_tas, var1, data_2nd, -1, 0, 0) 
              idx_mod := imod
          end if
          end if 

      else  ; not a reference data set
          A0 := extract_data_extended(imod, var0, data_1st, -1, 0, 0)
          A2 := extract_data_extended(imod, var1, data_2nd, -1, 0, 0)
          idx_mod := imod - 1
      end if

      if (dimsizes(dimsizes(A0)).gt.1) then
          tmp := area_operations(A0,latrange(0),latrange(1), \
                                    lonrange(0),lonrange(1),"average",True) 
      else
          tmp := A0
      end if 
      tmp_A0 := time_operations(tmp,-1,-1,"average","yearly",False);True)
      delete(tmp)

      ;; Convert units if appropriate
      if (isatt(variable_info, "plot_units")) then
          tmp_A0 := convert_units(tmp_A0,variable_info@plot_units)
      end if

      if (dimsizes(dimsizes(A2)).gt.1) then
          tmp := area_operations(A2,latrange_v2(0),latrange_v2(1), \
                                   lonrange_v2(0),lonrange_v2(1),"average",True)
      else
          tmp := A2
      end if
      tmp_A2 := time_operations(tmp,-1,-1,"average","yearly",True)
      delete([/A0,A2,tmp/])

      ;;set years with big vulkanic eruptions to _FillValues
      if (isatt(diag_script_info,"ec_volc")) then
          volc_years = diag_script_info@ec_volc
          do y = 0,dimsizes(tmp_A0&year)-1
              aY1reg(idx_mod,y) = where(any(tmp_A0&year(y).eq.volc_years), \
                                            tmp_A0@_FillValue,(/tmp_A0(y)/))
              aY2reg(idx_mod,y) = where(any(tmp_A2&year(y).eq.volc_years), \
                                            tmp_A2@_FillValue,(/tmp_A2(y)/))
          end do
      else
          aY1reg(idx_mod,:) = (/tmp_A0/)
          aY2reg(idx_mod,:) = (/tmp_A2/)
      end if
      delete([/tmp_A0,tmp_A2/])

      ;;remove 11-yr running mean for anomalies 
      if (anom .eq. True) then
          nmeanh = 5
          tmp1 = new((/nMax/),"float")
          tmp2 = new((/nMax/),"float")
          do n = nmeanh,nMax-nmeanh-1
              tmp1(n) = (/avg(aY1reg(idx_mod,n-nmeanh:n+nmeanh))/)
              tmp2(n) = (/avg(aY2reg(idx_mod,n-nmeanh:n+nmeanh))/)
          end do

          tmp1(:nmeanh-1) = tmp1(nmeanh)
          tmp1(nMax-nmeanh:) = (/tmp1(nMax-nmeanh-1)/)
          aY1reg(idx_mod,:)     = aY1reg(idx_mod,:) - tmp1
          delete(tmp1)
          
          tmp2(:nmeanh-1) = tmp2(nmeanh)
          tmp2(nMax-nmeanh:) = (/tmp2(nMax-nmeanh-1)/)
          aY2reg(idx_mod,:)     = aY2reg(idx_mod,:) - tmp2
          delete(tmp2) 
      end if 

      ;;remove trend
      if (detr .eq. True) then  
          tmpd := runave(aY1reg(idx_mod,:),2,0)
          tmpd(0)           = (/tmpd(1)/)
          aY1reg(idx_mod,:) = aY1reg(idx_mod,:)-tmpd
          delete(tmpd)
          
          tmpd := runave(aY2reg(idx_mod,:),2,0)
          tmpd(0)           = (/tmpd(1)/)
          aY2reg(idx_mod,:) = aY2reg(idx_mod,:)-tmpd
          delete(tmpd) 
      end if
      
      ;;check time range
      tmin := ind(aY2reg&year.eq.max((/timerange(0),toint(models@start_year(imod))/)))
      tmax := ind(aY2reg&year.eq.min((/timerange(1),toint(models@end_year(imod))/)))

      ;;calculate regression line:
      rc := regline(aY2reg(idx_mod,tmin:tmax), aY1reg(idx_mod,tmin:tmax))
      df := (/rc@nptxy-2 /)
      
      aYreg(idx_mod,tmin:tmax) = (/rc * aY2reg(idx_mod,tmin:tmax) + rc@yintercept/)
      arc(0,idx_mod) = (/rc/)
      arc(1,idx_mod) = (/rc@rstd/)
      ccr(idx_mod)   = (/escorc(aY2reg(idx_mod,tmin:tmax), aY1reg(idx_mod,tmin:tmax))/)
      var(0,idx_mod) = (/stddev(aY2reg(idx_mod,tmin:tmax))/)
      var(1,idx_mod) = (/stddev(aY1reg(idx_mod,tmin:tmax))/)

      delete([/rc, df, tmin, tmax/])
  end do     ;model loop
   
  ;;counting ref_models only as one model
  tmp_mod = models@name
  delete(models@name)
  models@name = new((/dim_MOD-iref/),typeof(tmp_mod))
  models@name(iref_co2:iref_tas-1) = "OBS"
  models@name(iref_tas:dim_MOD-iref-1) = tmp_mod(iref_tas+1:dim_MOD-1)
  delete(tmp_mod)
  dim_MOD = dimsizes(models@name)

  arc!1     = "model"
  arc&model = models@name 
  model     = arc&model

  ;;Optional output to NetCDF
  if (getenv("ESMValTool_write_netcdf").eq."True") then
      new_path = get_output_dir(get_work_dir(), basename(diag_script))
      new_path = new_path + basename(diag_script) + "_gIAV_" \
                          + timerange(0) + "-" + timerange(1) + ".nc"

      ;;Attach attributes to the results
      arc!0        = "case"
      arc&case     = (/"mean","stddev"/)
      arc!1        = "model"
      arc&model    = models@name
      arc@ncdf     = new_path
      arc@var      = "gIAV"
      arc@diag_script = (/diag_script/)
      arc@varIAV   = (/var/)
      arc@rIAV     = (/ccr/)

      ;;Write NetCDF output
      ncdf_outfile = ncdf_write(arc, new_path)

      ;;Write results of temporary list
      xml_name = getenv("ESMValTool_xml_name")
      suffix = get_file_suffix(xml_name, 0)
      name = str_sub_str(xml_name, suffix, "")
      temp_dir = get_output_dir(getenv("ESMValTool_wrk_dir"), "temp") + \
      name + "_gIAV.nc"

      ;;Remove previous file if already exists
      if (isfilepresent_esmval(temp_dir)) then
          system("rm -f " + temp_dir)
      end if

      ncdf_char = tochar(new_path)
      temp_list_new = new((/1, dimsizes(ncdf_char)/), character)
      temp_list_new(0, :) = ncdf_char

      ;; Create new file and add list
      temp = addfile(temp_dir, "c")
      temp->temp_list = temp_list_new
  end if

  ;; Get environment variables for plotting
  plot_dir = getenv("ESMValTool_plot_dir")
  diag_script_base = basename(diag_script)
  output_dir = get_output_dir(plot_dir, diag_script_base)

  ;;Define path and names of plot files  
  file_type = getenv("ESMValTool_output_file_type")
  if(ismissing(file_type)) then
      file_type = "PS"
  end if

  ;;define arrays for plotting
  f  = sqrt(dim_MOD)
  mm = floattointeger(ceil(f))
  colors   = project_style(diag_script_info, "colors")
  thicks   = project_style(diag_script_info, "thicks")
  plot_num = integertochar(ispan(97, 122, 1))
  plot     = new(dim_MOD, "graphic")

  ;;Call plot scripts
  plot_file = "corr_" + var0 + "-" + var1 + "_anom_" + timerange(0) + "-" + timerange(1) 
  wks       = gsn_open_wks(file_type, output_dir + plot_file)

  ;;Create panel for each model
  do imod=0,dim_MOD-1

      MainStg  = plot_num(imod) + ") " + models@name(imod) \
                 + "  ~F8~g~F21~~B~IAV~N~= " \
                 + sprintf("%4.2f", arc(0,imod)) + " GtC/K/yr;  r= " \
                 + sprintf("%4.2f", ccr(imod))
      XStg     = "Anomaly in Trop. Temp. (" + UNITS_v2 + ")"
      YStg     = "Anomaly in d(land+ocean CO~B~2~N~ flux)/dt (" + UNITS + ")"

      data_arr = (/aY2reg(imod,:), aY1reg(imod,:)/)
      data_arr!0      = "models"
      data_arr&models = (/XStg,YStg/)

      data_arr@diag_script         = (/diag_script/)
      data_arr@res_tiMainString    = MainStg
      data_arr@res_gsnCenterString = ""
      data_arr@res_tiMainFontHeightF = 0.025
      data_arr@res_xyMarkerColors  = "white"
      data_arr@res_trXMinF         = min(aY2reg(imod,:));-0.4
      data_arr@res_trXMaxF         = max(aY2reg(imod,:));0.6
      data_arr@res_trYMinF         = min(aY1reg(imod,:))
      data_arr@res_trYMaxF         = max(aY1reg(imod,:))

      if (.not.all(ismissing(data_arr))) then
          plot(imod) = scatterplot(wks, data_arr, var0, False, False)

       ;;Add year strings as markers
          add_num_markers(wks, plot(imod), False, aY2reg(imod,:), aY1reg(imod,:))
          lineres    = True
          lineres@gsLineColor = colors(imod)
          str        = unique_string("line")
          plot@$str$ = gsn_add_polyline(wks, plot(imod), aY2reg(imod,:), \
                                                         aYreg(imod,:), lineres)
      end if
      delete([/data_arr,MainStg,XStg,YStg/])
  end do
  
  ;;Common resources for panel plots
  pan             = True 
  pan@gsnFrame    = False 
  pan@gsnMaximize = True 
  pan@gsnPanelXWhiteSpacePercent = 10
  pan@gsnPanelYWhiteSpacePercent = 5
  gsn_panel(wks,plot,(/mm,mm/),pan)
  frame (wks)

  info_output("gv " + output_dir + plot_file + "." + file_type, verbosity, 1)
  delete([/plot,pan,y/])

end
;#################################################################################
