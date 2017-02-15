;;############################################################################
;; create_co2flux.ncl
;; Author: Sabrina Wenzel (DLR, Germany)                                      
;; ESMVal project PI: Veronika Eyring (DLR, Germany)                          
;;############################################################################
;; Description
;;    Reads land and ocean carbon fluxes and creates a new variable to be saved
;;    as T0M
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    * ts_latrange2 min and max latitudes for spatial average
;;    * ts_lonrange2 min and max longitudes for spatial average
;;
;; Caveats
;;    The CO2 flux can not be calculated using derive var, it therfore needs
;;    to be calculated within the diagnostic
;;
;; Modification history
;;    * 20151027-A_laue_ax: added call to 'write_references'
;;    * 201509-A_Sabrina Wenzel:  written       
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

begin
  plot_script = "create_co2flux"
  
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

  ;; fields supported
  field_number = field_types
  field_0 = field_number(0)

  var0 = variables(0)
  var1 = variables(1)
  var2 = variables(2)

  data_1st = True
  data_2nd = True
  data_3nd = True
  do idx_mod = 0, dim_MOD - 1
    data_1st = read_data(idx_mod, var0, field_0)
    data_2nd = read_data(idx_mod, var1, field_0)
    data_3nd = read_data(idx_mod, var2, field_0)
  end do
end

;; Load var_info
begin
    var0deffile = "interface_data/" + var0 + "_info.tmp"
    var1deffile = "interface_data/" + var1 + "_info.tmp"
    var2deffile = "interface_data/" + var2 + "_info.tmp"
    loadscript(var0deffile)
    loadscript(var1deffile)
    loadscript(var2deffile)
end

;; Main
begin 
  ;;****************************************************************** 
  ;;check the input parameters 
    
  if(isvar("latrange")) then
    delete(latrange)
  end if
  if(isatt(diag_script_info, "ts_latrange")) then
    latrange = diag_script_info@ts_latrange
  else
    print("fatal:plot requires variable attribute 'ts_latrange'")
    status_exit(1)
  end if
  
  if(isvar("lonrange")) then
    delete(lonrange)
  end if
  if(isatt(diag_script_info, "ts_lonrange")) then
    lonrange = diag_script_info@ts_lonrange
  else
    print("fatal:plot requires variable attribute 'ts_lonrange'")
    status_exit(1)
  end if

  if(isvar("latrange_v2")) then
    delete(latrange_v2)
  end if
  if(isatt(diag_script_info, "ts_latrange2")) then
    latrange_v2 = diag_script_info@ts_latrange2
  else
    print("fatal:plot requires variable attribute 'ts_latrange2'")
    status_exit(1)
  end if
  
  if(isvar("lonrange_v2")) then
    delete(lonrange_v2)
  end if
  if(isatt(diag_script_info, "ts_lonrange2")) then
    lonrange_v2 = diag_script_info@ts_lonrange2
  else
    print("fatal:plot requires variable attribute 'ts_lonrange2'")
    status_exit(1)
  end if

  ;;================================================================== 
  ;;loop for models
  do imod = 0,dim_MOD-1 
    
      ;;get out file path
      new_path = getenv("ESMValTool_climo_dir")
      new_path = new_path + "/CMIP5_ETHZ/"

      if (.not.isfilepresent(new_path)) then
          system("mkdir " + new_path)
      end if
      delete(new_path)
  end do

  do imod = 0,dim_MOD-1  
      info_output(models@name(imod), verbosity, 1) 

      ;;get out file name
      new_path = getenv("ESMValTool_climo_dir")
      new_path = new_path + "/" + models@project(imod) + "/"

      new_path = new_path + "co2flux_Amon_" \
                          + models@name(imod) + "_" \
                          + models@experiment(imod) + "_" \
                          + models@ensemble(imod) + "_" \
                          + models@start_year(imod) + "_" \
                          + models@end_year(imod) + ".nc"

      ;;Only process if file is not existing
      if (.not.isfilepresent(new_path).or.getenv("ESMValTool_force_processing").eq."True") then  

          ;;get data 
          A0 = extract_data_extended(imod,var0,data_1st,-1,0,0)
          A1 = extract_data_extended(imod,var1,data_2nd,-1,0,0) 

          ;;extract specific region
          ;;read coastal files for land
          filsl = systemfunc ("ls "+ models@dir(imod) + "/historical/fx/sftlf/" \
                                   + models@name(imod) + "/r0i0p0/sftlf_fx_" \
                                   + models@name(imod) + "_historical_r0i0p0.nc")
          ;;read area files for land
          fill = systemfunc ("ls "+models@dir(imod)+"/historical/fx/areacella/" \
                                  +models@name(imod)+"/r0i0p0/areacella_fx_" \
                                  +models@name(imod)+"_historical_r0i0p0.nc")

          req_files := (/filsl, fill/)
          req_vars  := (/"sftlf", "areacella"/)
          req_lines := (/"154", "158"/)
          if (.not.all(isfilepresent(req_files))) then
              error_msg("f", diag_script, "", "missing, could not find fx-file " + \
              "for " + req_vars(ind(.not.isfilepresent(req_files))) + ": " + \
              " please adjust lokal file path in line " + \
              req_lines(ind(.not.isfilepresent(req_files))))
          else 
              dumsl = addfile(filsl, "r")
              sftlf = dumsl->sftlf(:,:)

              duml  = addfile(fill, "r")
              areacella = duml->areacella(:,:)
          end if

          A0_latlon = A0 * conform(A0,areacella,(/1,2/)) * (conform(A0,sftlf,(/1,2/)))/100

          tmp_A0 = dim_sum_Wrap(dim_sum_Wrap(A0_latlon))
          tmp_A0!0 = "time"
          tmp_A0&time = (/A0&time/)
          
          delete([/req_files,req_vars,req_lines/])
          delete([/sftlf,filsl,dumsl/])
          delete([/areacella,fill,duml,A0_latlon/])

          ;;read area files for ocean
          filo = systemfunc ("ls "+models@dir(imod)+"/historical/fx/areacello/" \
                                  +models@name(imod)+"/r0i0p0/areacello_fx_" \
                                  +models@name(imod)+"_historical_r0i0p0.nc")

          ;;read coastal files for ocean
          filso = systemfunc ("ls "+ models@dir(imod) + "/historical/fx/sftof/" \
                                   + models@name(imod) + "/r0i0p0/sftof_fx_" \
                                   + models@name(imod) + "_historical_r0i0p0.nc")

          if (.not.all(isfilepresent(filo))) then
              error_msg("f", diag_script, "", "missing, could not find fx-file " + \
              "for areacello:  please adjust lokal file path in line 188")
          else 
              dumo  = addfile(filo, "r")
              areacello = dumo->areacello(:,:)

              if (models@name(imod).ne."MIROC-ESM") then
                  if (ismissing(filso)) then
                      A1_latlon = A1*conform(A1,areacello,(/1,2/))
                  else
                      if (.not.all(isfilepresent(filso))) then
                          error_msg("f", diag_script, "", "missing, could not find fx-file " + \
                          "for sftof: please adjust lokal file path in line 193")
                      else
                          dumso = addfile(filso, "r")
                          sftof = dumso->sftof(:,:)
                      end if
                      
                      A1_latlon = A1*conform(A1,areacello,(/1,2/))*(conform(A1,sftof,(/1,2/)))/100
                      delete([/sftof,filso,dumso/])
                  end if
              else
                  A1_latlon = A1
              end if
          end if
          
          ;;calculate spatial average
          tmp_A1 = dim_sum_Wrap(dim_sum_Wrap(A1_latlon))
          tmp_A1!0 = "time"
          tmp_A1&time = (/A1&time/)
          
          delete([/areacello,filo,dumo,A1_latlon/])

          ;;add land and ocean carbon fluxes
          tmp_sum = tmp_A0 + tmp_A1
          copy_VarAtts(A0,tmp_sum)
          copy_VarCoords_2(A0,tmp_sum)
          delete([/tmp_A0,tmp_A1/])

          ;;Convert units if appropriate
          if (isatt(variable_info, "plot_units")) then
              tmp_sum = convert_units(tmp_sum,variable_info@plot_units)
          end if

          ;;Output to NetCDF
          tmp_sum!0           = "time"
          tmp_sum&time        = A0&time
          tmp_sum@ncdf        = new_path
          tmp_sum@var         = "co2flux"
          tmp_sum@diag_script = diag_script
          ncdf_outfile = ncdf_write(tmp_sum, new_path)
          delete([/ncdf_outfile,A0,A1,tmp_sum/])
      else
          info_output("No recalculation needed for co2flux!", verbosity, 1)
      end if
      delete(new_path)

      ;;get out file name
      new_path = getenv("ESMValTool_climo_dir")
      new_path = new_path + "/CMIP5_ETHZ/"

      ;;Output to NetCDF   
      new_path = new_path + var2 + "_Amon_" \
                          + models@name(imod) + "_" \
                          + models@experiment(imod) + "_" \
                          + models@ensemble(imod) + "_" \
                          + models@start_year(imod) + "_" \
                          + models@end_year(imod) + ".nc" 

      ;;Only process if file is not existing
      if (.not.isfilepresent(new_path).or.getenv("ESMValTool_force_processing").eq."True") then 
          A2 = extract_data_extended(imod,var2,data_3nd,-1,0,0) 

          A2@ncdf        = new_path
          A2@var         = var2
          A2@diag_script = diag_script
          ncdf_outfile = ncdf_write(A2, new_path)
          delete([/A2,ncdf_outfile/])
      else
          info_output("No recalculation needed for " + var2 + "!", verbosity, 1)
      end if
      delete(new_path)
  end do
  
  info_output(">>>>>>>> Leaving " + plot_script + ".ncl", verbosity, 4)
end
;#################################################################################
