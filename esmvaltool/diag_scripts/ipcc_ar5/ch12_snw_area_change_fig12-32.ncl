; #############################################################################
; SNOW AREA CHANGE
; Author: Ruth Lorenz (ETH Zurich, Switzerland)
; CRESCENDO project
; #############################################################################
;
; Description
;    Calculate Snow area in a region (e.g Northern Hemisphere)
;    and season (e.g. NH spring March & April)
;    relative to reference period (e.g 1986-2005) as in Fig. 12.32 of IPCC AR5
;
; Modification history
;    20180718-A_lore_ru: ported to v2 and bugs in area calculation fixed
;    20170120-A_lore_ru: written.
;
; #############################################################################

; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
load "./interface_scripts/interface.ncl"

load "./diag_scripts/shared/latlon.ncl"
load "./diag_scripts/shared/statistics.ncl"

begin
  ; ##############################################################
  ; # Fetch general parameters, set in namelist_collins13ipcc.xml#
  ; # passed via environment variables by python code            #
  ; ##############################################################
  enter_msg(diag_script, "")

  ; 'diag_script' and 'variables' are fetched from the
  ; above 'ncl.interface' file
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (var: " + variables(0) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; 'models', 'variables' and 'field_types' are fetched from the
  ; above 'ncl.interface' file
  dim_MOD = dimsizes(dataset_info@dataset)
  dim_VAR = dimsizes(variables)
  var0 = variables(0)
  field_type0 = field_types(0)

  ; ###########################################
  ; # References Information                  #
  ; ###########################################
;  write_references(diag_script, \  ; script name
;                   "A_lore_ru", \  ; authors
;                   "", \  ; contributors
;                   "D_collins13ipcc",    \  ; diag_references
;                   "",          \  ; obs_references
;                   "P_crescendo")   ; proj_references
  ; Check required diag_script_info attributes
  req_atts = (/"scenarios", "syears", "eyears", "begin_ref_year", \
               "end_ref_year", "months"/)
  exit_if_missing_atts(diag_script_info, req_atts)
end

begin
  ; Get environment variables
  if (isatt(diag_script_info, "plot_type")) then
    file_type = diag_script_info@plot_type
  else
    file_type = "pdf"
  end if

  ; Basename of diag_script
  diag_script_base = basename(diag_script)

  ; Output netcdf directory
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)

  ; Output plot directory
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)
end
begin
  ; #############################################
  ; # Get parameters from ./variable_defs/*.ncl #
  ; # passed via the 'info' attribute           #
  ; #############################################
  if(isvar("MyParam")) then
    delete(MyParam)
  end if
  if(isatt(variable_info, "long_name")) then
    MyParam = variable_info@long_name
    log_info(" MyParam = " + MyParam)
  else
    error_msg("f", diag_script, "", "missing variable attribute " + \
              "'long_name'")
  end if
  ; ###########################################
  ; # Get data and average annual and globally#
  ; ###########################################
  ; get data from first model for historical period
  imod = 0  ; NCL array indices start from zero
  log_info("processing " + dataset_info@dataset(imod) + "_" \
           + dataset_info@exp(imod) + "_" \
           + dataset_info@ensemble(imod))

  ; See ./interface_scripts/data_handling.ncl
  A0 = read_data(imod, var0, field_type0)

  ; Check dimensions
  dims = getvardims(A0)
  if (dimsizes(dims) .lt. 2) then
    error_msg("f", diag_script, "", dimsizes(dims) + " dimensions, " + \
              "need 2 or 3")
  end if
  idx = ind(dims .eq. "lat")
  if (ismissing(idx)) then
    error_msg("f", diag_script, "", "no lat dimension")
  end if
  idx = ind(dims .eq. "lon")
  if (ismissing(idx)) then
    error_msg("f", diag_script, "", "no lon dimension")
  end if

  ; how many months included in analysis?
  months = diag_script_info@months
  dim_SEAS = dimsizes(months)

  ; we need to calculate the snow area using sftgif and sftlf
  list_sftgif = read_fx_data(imod, "sftgif")
  var_idx = ListIndexFromName(list_sftgif, "sftgif")
  land_ice_sftgif = list_sftgif[var_idx]

  ; read sftlf (land sea mask)
  list_sftlf = read_fx_data(imod, "sftlf")
  var_idx2 = ListIndexFromName(list_sftlf, "sftlf")
  land_sea_sftlf = list_sftlf[var_idx2]
  sftifl = land_sea_sftlf - land_ice_sftgif
  delete([/var_idx, var_idx2/])

  sftifl_conf = conform(A0, sftifl, (/1, 2/))
  tmp2 = (A0 * sftifl_conf) / 100  ; sftifl in percent -> / 100
  ; where > 5
  tmp3 = where(tmp2 .gt. 5., 1.0, 0)
  ; Add coordinates
  copy_VarCoords(A0, tmp3)
  ; Average over globe/region
  ; (function in ./diag_scripts/lib/ncl/latlon.ncl)
  if (isatt(diag_script_info, "ts_minlat")) then
    tmp4 = area_operations(tmp3, diag_script_info@ts_minlat, \
                           diag_script_info@ts_maxlat, \
                           diag_script_info@ts_minlon, \
                           diag_script_info@ts_maxlon, "average", True)
  else
    tmp4 = coswgt_areaave(tmp3)
  end if
  tmp5 = (tmp4 * 5.10067e14) / 2.0
  ; Add coordinates
  tmp5!0 = "time"
  tmp5&time = A0&time
  tmp6 = extract_season(tmp5, months)
  data2_glob = time_operations(tmp6, -1, -1, "average", "yearly", False)
  delete([/tmp2, tmp3, tmp4, tmp5, tmp6, sftifl_conf/])

  dim_file = dimsizes(data2_glob)  ; file size of historical file

  ; How many historical model runs?
  ind_hist = ind(dataset_info@exp .eq. "historical")
  dim_MOD_hist = dimsizes(dataset_info@dataset(ind_hist))
  ref_data2 = new((/dim_MOD_hist/), float, A0@_FillValue)
  data2_rel = new((/dim_MOD_hist, dim_file(0)/), float, \
                  A0@_FillValue)

  ; How many scenarios?
  dim_scen = dimsizes(diag_script_info@scenarios)
  ; start and end years of timeperiods
  syears = diag_script_info@syears
  eyears = diag_script_info@eyears
  nperiods = dimsizes(syears)
  nyears_full = toint(eyears(nperiods - 1)) - toint(syears(0)) + 1
  nyearsp1 =  toint(eyears(nperiods - 1)) - toint(eyears(0))

  ; How many model runs in 1st scenario?
  ind_26 = ind(dataset_info@exp .eq. diag_script_info@scenarios(0))
  dim_MOD_26 = dimsizes(dataset_info@dataset(ind_26))
  ; indicate array for scenarios
  data26_rel = new((/dim_MOD_26, nyearsp1/), float, \
                   A0@_FillValue)

  if (dim_scen .gt. 1) then
    ; How many model runs in 2nd sceanrio?
    ind_45 = ind(dataset_info@exp .eq. diag_script_info@scenarios(1))
    dim_MOD_45 = dimsizes(dataset_info@dataset(ind_45))
    data45_rel = new((/dim_MOD_45, nyearsp1/), float, \
                     A0@_FillValue)
  end if
  if (dim_scen .gt. 2) then
    ; How many model runs in 3rd scenario?
    ind_60 = ind(dataset_info@exp .eq. diag_script_info@scenarios(2))
    dim_MOD_60 = dimsizes(dataset_info@dataset(ind_60))
    data60_rel = new((/dim_MOD_60, nyearsp1/), float, \
                     A0@_FillValue)
  end if
  if (dim_scen .gt. 3) then
    ; How many model runs in 4th scenario?
    ind_85 = ind(dataset_info@exp .eq. diag_script_info@scenarios(3))
    dim_MOD_85 = dimsizes(dataset_info@dataset(ind_85))
    data85_rel = new((/dim_MOD_85, nyearsp1/), float, \
                     A0@_FillValue)
  end if
  i = 0
  i26 = 0
  i45 = 0
  i60 = 0
  i85 = 0
  do imod = 0, dim_MOD_hist - 1
    if (imod .eq. 0) then
      delete(A0)
    else
      A0 = read_data(imod, var0, field_type0)
      ; we need to calculate the snow area using sftgif (glacier) and sftlf
      list_sftgif = read_fx_data(imod, "sftgif")
      var_idx = ListIndexFromName(list_sftgif, "sftgif")
      land_ice_sftgif = list_sftgif[var_idx]

      ; read sftlf (land sea mask)
      list_sftlf = read_fx_data(imod, "sftlf")
      var_idx2 = ListIndexFromName(list_sftlf, "sftlf")
      land_sea_sftlf = list_sftlf[var_idx2]
      sftifl = land_sea_sftlf - land_ice_sftgif
      delete([/var_idx, var_idx2/])

      sftifl_conf = conform(A0, sftifl, (/1, 2/))
      tmp2 = (A0 * sftifl_conf) / 100
      ; where > 5 kg m**-2 considered snow covered
      tmp3 = where(tmp2 .gt. 5., 1., 0)
      ; Add coordinates
      copy_VarCoords(A0, tmp3)
      ; Average over globe/region
      ; (function in ./diag_scripts/lib/ncl/latlon.ncl)
      if (isatt(diag_script_info, "ts_minlat")) then
        tmp4 = area_operations(tmp3, diag_script_info@ts_minlat, \
                               diag_script_info@ts_maxlat, \
                               diag_script_info@ts_minlon, \
                               diag_script_info@ts_maxlon, "average", \
                               True)
      else
        tmp4 = coswgt_areaave(tmp3)
      end if
      tmp5 = (tmp4 * 5.10067e14) / 2.0
      ; Add coordinates
      tmp5!0 = "time"
      tmp5&time = A0&time
      tmp6 = extract_season(tmp5, months)
      data2_glob = time_operations(tmp6, -1, -1, "average", "yearly", False)
      delete([/tmp2, tmp3, tmp4, tmp5, tmp6, sftifl_conf, A0/])
    end if

    ind_start = ind(data2_glob&year .eq. diag_script_info@begin_ref_year)
    ind_end = ind(data2_glob&year .eq. diag_script_info@end_ref_year)
    ref_data2(imod) = dim_avg(data2_glob(ind_start : ind_end))
    print(ref_data2(imod))
    dim_data2g = dimsizes(data2_glob)
    if (dim_data2g .ne. dim_file(0)) then
      tdiff = dim_file(0) - dim_data2g
      data2_rel(imod, 0 : tdiff - 1) = data2_rel@_FillValue
      data2_rel(imod, tdiff :) =  (data2_glob / ref_data2(imod)) - 1
    else
      data2_rel(imod, :) = (data2_glob / ref_data2(imod)) - 1
    end if
    ; clean up
    delete([/data2_glob/])

    do s = 0, dim_scen - 1
      ; find all other runs from this model and calculate relative ts
      imodrcp = ind(dataset_info@dataset .eq. dataset_info@dataset(imod) \
                    .and. dataset_info@exp .eq. diag_script_info@scenarios(s) \
                    .and. dataset_info@ensemble .eq. \
                    dataset_info@ensemble(imod))
      if (.not.all(ismissing(imodrcp))) then
        A1 = read_data(imodrcp, var0, field_type0)
        ; we need to calculate the snow area using sftgif (glacier) and sftlf
        list_sftgif = read_fx_data(imodrcp, "sftgif")
        var_idx = ListIndexFromName(list_sftgif, "sftgif")
        land_ice_sftgif = list_sftgif[var_idx]

        ; read sftlf (land sea mask)
        list_sftlf = read_fx_data(imodrcp, "sftlf")
        var_idx2 = ListIndexFromName(list_sftlf, "sftlf")
        land_sea_sftlf = list_sftlf[var_idx2]
        sftifl = land_sea_sftlf - land_ice_sftgif  ; land areas sans glaciers
        delete([/var_idx, var_idx2/])
        sftifl_conf = conform(A1, sftifl, (/1, 2/))
        tmp2 = (A1 * sftifl_conf) / 100  ; sftifl in percent -> / 100
        ; where > 5 kg m**-2 considered snow covered
        tmp3 = where(tmp2 .gt. 5., 1., 0)
        ; Add coordinates
        copy_VarCoords(A1, tmp3)
        ; Average over globe/region
        ; (function in ./diag_scripts/lib/ncl/latlon.ncl)
        if (isatt(diag_script_info, "ts_minlat")) then
          tmp4 = area_operations(tmp3, diag_script_info@ts_minlat, \
                                 diag_script_info@ts_maxlat, \
                                 diag_script_info@ts_minlon, \
                                 diag_script_info@ts_maxlon, "average", \
                                 True)
        else
          tmp4 = coswgt_areaave(tmp3)
        end if
        tmp5 = (tmp4 * 5.10067e14) / 2.0
        ; Add coordinates
        tmp5!0 = "time"
        tmp5&time = A1&time
        tmp6 = extract_season(tmp5, months)
        data2_glob = time_operations(tmp6, -1, -1, "average", "yearly", False)
        printMinMax(data2_glob, True)
        delete([/tmp2, tmp3, tmp4, tmp5, tmp6, sftifl_conf/])
        dim_data2g = dimsizes(data2_glob)
        if (s .eq. 0) then
          if (dim_data2g .ne. nyearsp1) then
            if (dim_data2g .gt. nyearsp1) then
              error_msg("w", diag_script, "", "Length of dataset " + \
                        dataset_info@dataset(imod) + " in recipe is " + \
                        "longer than period given in recipe (" + nyearsp1 + \
                        " years)")
            else
              data26_rel(i26, 0 : dim_data2g - 1) = \
                (data2_glob / ref_data2(imod)) - 1
                data26_rel(i26, dim_data2g :) = data26_rel@_FillValue
            end if
          else
            data26_rel(i26, :) = (data2_glob / ref_data2(imod)) - 1
          end if
          i26 = i26 + 1
        end if
        if (s .eq. 1) then
          if (dim_data2g .ne. nyearsp1) then
            if (dim_data2g .gt. nyearsp1) then
              error_msg("w", diag_script, "", "Length of dataset " + \
                        dataset_info@dataset(imod) + " in recipe is " + \
                        "longer than period given in recipe (" + nyearsp1 + \
                        " years)")
            else
              data45_rel(i45, 0 : dim_data2g - 1) = \
                (data2_glob / ref_data2(imod)) - 1
              data45_rel(i45, dim_data2g:) =  data45_rel@_FillValue
            end if
          else
            data45_rel(i45, :) = (data2_glob / ref_data2(imod)) - 1
          end if
          i45 = i45 + 1
        end if
        if (s .eq. 2) then
          if (dim_data2g.ne.nyearsp1) then
            if (dim_data2g .gt. nyearsp1) then
              error_msg("w", diag_script, "", "Length of dataset " + \
                        dataset_info@dataset(imod) + " in recipe is " + \
                        "longer than period given in recipe (" + nyearsp1 + \
                        " years)")
            else
              data60_rel(i60, 0 : dim_data2g - 1) = \
                (data2_glob / ref_data2(imod)) - 1
              data60_rel(i60, dim_data2g:) =  data60_rel@_FillValue
            end if
          else
            data60_rel(i60, :) = (data2_glob / ref_data2(imod)) - 1
          end if
          i60 = i60 + 1
        end if
        if (s .eq. 3) then
          if (dim_data2g.ne.nyearsp1) then
            if (dim_data2g .gt. nyearsp1) then
              error_msg("w", diag_script, "", "Length of dataset " + \
                        dataset_info@dataset(imod) + " in recipe is " + \
                        "longer than period given in recipe (" + nyearsp1 + \
                        " years)")
            else
              data85_rel(i85, 0 : dim_data2g - 1) = \
                (data2_glob / ref_data2(imod)) - 1
              data85_rel(i85, dim_data2g :) = data85_rel@_FillValue
            end if
          else
            data85_rel(i85, :) = (data2_glob / ref_data2(imod)) - 1
          end if
          i85 = i85 + 1
        end if
        ; clean up
        delete([/A1, data2_glob/])
      end if
      delete(imodrcp)
    end do  ; scenario
    delete([/land_ice_sftgif, land_sea_sftlf, sftifl/])
  end do  ; imod

  ; average and standard deviation over all datasets
  datahist_avg = dim_avg_n_Wrap(data2_rel, 0)
  datahist_std = dim_stddev_n_Wrap(data2_rel, 0)

  ; put all data from 1850-2100 into mean_val and std_val
  mean_val = new((/dim_scen + 1, nyears_full/), typeof(data2_rel),\
                 data2_rel@_FillValue)
  std_val = new((/dim_scen + 1, nyears_full/), typeof(data2_rel),\
                data2_rel@_FillValue)
  p0 = toint(eyears(0)) - toint(syears(0))
  p1 = toint(eyears(0)) - toint(syears(0)) + 1
  mean_val(0, 0 : p0) = datahist_avg
  std_val(0, 0 : p0) = datahist_std

  data26_avg = dim_avg_n_Wrap(data26_rel, 0)
  data26_std = dim_stddev_n_Wrap(data26_rel, 0)

  mean_val(1, p1 :) = data26_avg
  std_val(1, p1 :) = data26_std

  if (dim_scen .gt. 1) then
    data45_avg = dim_avg_n_Wrap(data45_rel, 0)
    data45_std = dim_stddev_n_Wrap(data45_rel, 0)
    mean_val(2, p1 :) = data45_avg
    std_val(2, p1 :) = data45_std
  end if
  if (dim_scen .gt. 2) then
    data60_avg = dim_avg_n_Wrap(data60_rel, 0)
    data60_std = dim_stddev_n_Wrap(data60_rel, 0)
    mean_val(3, p1 :) = data60_avg
    std_val(3, p1 :) = data60_std
  end if
  if (dim_scen .gt. 3) then
    data85_avg = dim_avg_n_Wrap(data85_rel, 0)
    data85_std = dim_stddev_n_Wrap(data85_rel, 0)
    mean_val(4, p1 :) = data85_avg
    std_val(4, p1 :) = data85_std
  end if
  ; check number of datasets for timeperiods in scenarios
  ; if there are less than 5 datasets set to missing
  runs = new((/1 + dim_scen, nperiods/), integer)
  runs(0, 0) = dimsizes(ind(dataset_info@exp .eq. "historical" .and. \
                            dataset_info@end_year .eq. \
                            diag_script_info@eyears(0)))
  do scen = 0, dim_scen - 1
    do p = 0, nperiods - 1
      runs(scen + 1, p) = dimsizes(ind(dataset_info@exp .eq. \
                                       diag_script_info@scenarios(scen) \
                                       .and. dataset_info@end_year .eq. \
                                       diag_script_info@eyears(p)))
      if (runs(scen + 1, p) .lt. 5) .and. (p .ne. 0) then
        p2 = toint(eyears(p - 1)) - toint(syears(0))
        p3 = toint(eyears(p)) - toint(syears(0))
        mean_val(scen + 1, p2 : p3) = mean_val@_FillValue
        std_val(scen + 1, p2 : p3) = std_val@_FillValue
      end if
    end do
  end do

  dims_data = dimsizes(mean_val)
  to_plot = new((/3 * dims_data(0), dims_data(1)/), float)

  do j = 0, dims_data(0) - 1
    ; if variable sw or lw TOA flux: change direction
    if (var0 .eq. "rlut" .or. var0 .eq. "rsut") then
      to_plot(j * 3, :) = - mean_val(j, :)
    else
      to_plot(j * 3, :) = mean_val(j, :)
    end if
    to_plot((j * 3) + 1, :) = to_plot(j * 3, :) + 1. * std_val(j, :)
    to_plot((j * 3) + 2, :) = to_plot(j * 3, :) - 1. * std_val(j, :)
  end do
  copy_VarMeta(mean_val, to_plot)
  to_plot!0 = "metric"
  if (dim_scen .eq. 1) then
    to_plot&metric = (/"mean", "+1.std", "-1.std", \
                       "mean", "+1.std", "-1.std"/)
  else if (dim_scen .eq. 2) then
    to_plot&metric = (/"mean", "+1.std", "-1.std", \
                       "mean", "+1.std", "-1.std", \
                       "mean", "+1.std", "-1.std"/)
  else if (dim_scen .eq. 3) then
    to_plot&metric = (/"mean", "+1.std", "-1.std", \
                       "mean", "+1.std", "-1.std", \
                       "mean", "+1.std", "-1.std", \
                       "mean", "+1.std", "-1.std"/)
  else if (dim_scen .eq. 4) then
    to_plot&metric = (/"mean", "+1.std", "-1.std", \
                       "mean", "+1.std", "-1.std", \
                       "mean", "+1.std", "-1.std", \
                       "mean", "+1.std", "-1.std", \
                       "mean", "+1.std", "-1.std"/)
  end if
  end if
  end if
  end if
  to_plot!1 = "year"
  to_plot&year = ispan(toint(syears(0)), toint(eyears(nperiods - 1)), 1)
  ; ###########################################
  ; # Other Metadata: diag_script, var        #
  ; ###########################################
  ; Add to to_plot, as attributes without prefix
  if(isatt(to_plot, "diag_script")) then  ; Add to existing entries
    temp = to_plot@diag_script
    delete(to_plot@diag_script)
    to_plot@diag_script = array_append_record(temp, (/diag_script/), 0)
    delete(temp)
  else  ; Add as new attribute
    to_plot@diag_script = (/diag_script/)
  end if
  to_plot@var = var0  ; Overwrite existing entry
  if(isatt(variable_info, "long_name")) then
    to_plot@var_long_name = variable_info@long_name
  end if
  if(isatt(variable_info, "units")) then
    to_plot@var_units     = variable_info@units
  end if

  runs@var = "nr"
  runs!0 = "scenario"
  runs&scenario = diag_script_info@label
  runs!1 = "period"
  runs&period = diag_script_info@syears
  runs@diag_script = (/diag_script/)
  runs@var_long_name = "number of model runs per scenario and period"

  ; ###########################################
  ; # Output to netCDF               #
  ; ###########################################
  ; Output (controlled by diag_script_info)
  to_plot@ncdf = variable_info@diagnostic + ".nc"
  ncdf_file = work_dir + "/" + to_plot@ncdf
  ncdf_file@existing = "overwrite"
  ; Function in ~/interface_scripts/auxiliary.ncl
  ncdf_outfile = ncdf_write(to_plot, ncdf_file)
  if(isatt(diag_script_info, "model_nr") .and. \
     diag_script_info@model_nr .eq. True) then
    ncdf_file@existing = "append"
    ncdf_outfile = ncdf_write(runs, ncdf_file)
  end if
end
