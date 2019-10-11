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
;    20180718-A_lorenz_ruth: ported to v2 and bugs in area calculation fixed
;    20170120-A_lorenz_ruth: written.
;
; #############################################################################

; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/statistics.ncl"

begin
  ; ##############################################################
  ; # Fetch general parameters, set in namelist_collins13ipcc.xml#
  ; # passed via environment variables by python code            #
  ; ##############################################################
  enter_msg(DIAG_SCRIPT, "")

  ; 'datasets', 'variables' are fetched from the above 'interface.ncl' file
  dim_VAR = ListCount(variable_info)
  if (dim_VAR .eq. 3) then
    var0 = variable_info[0]@short_name
    info_items0 = select_metadata_by_name(input_file_info, var0)
    var1 = variable_info[1]@short_name
    info_items1 = select_metadata_by_name(input_file_info, var1)
    var2 = variable_info[2]@short_name
    info_items2 = select_metadata_by_name(input_file_info, var2)
  else
    error_msg("f", DIAG_SCRIPT, "", "3 variables are needed, " + \
              "'snw', 'sftlf' and 'sftgif'")
  end if

  ; Match snw variable to correct variable
  if (var0 .eq. "snw") then
    var_snw = var0
    info_items = info_items0
    idx_snw = 0
  elseif (var1 .eq. "snw") then
    var_snw = var1
    info_items = info_items1
    idx_snw = 1
  elseif (var2 .eq. "snw") then
    var_snw = var2
    info_items = info_items2
    idx_snw = 2
  end if
  dim_MOD = ListCount(info_items)

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

  ; Output netcdf directory
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)

end
begin
  ; #############################################
  ; # Get parameters from ./variable_defs/*.ncl #
  ; # passed via the 'info' attribute           #
  ; #############################################
  if(isvar("MyParam")) then
    delete(MyParam)
  end if
  if(isatt(variable_info[0], "long_name")) then
    MyParam = variable_info[0]@long_name
    log_info(" MyParam = " + MyParam)
  else
    error_msg("f", DIAG_SCRIPT, "", "missing variable attribute " + \
              "'long_name'")
  end if
  ; ###########################################
  ; # Get data and average annual and globally#
  ; ###########################################
  ; get data from first model for historical period
  imod = 0  ; NCL array indices start from zero
  log_info("processing " + info_items[imod]@dataset + "_" \
           + info_items[imod]@exp + "_" \
           + info_items[imod]@ensemble)

  ; See ./interface_scripts/data_handling.ncl
  A0 = read_data(info_items[imod])

  ; Check dimensions
  dims = getvardims(A0)
  if (dimsizes(dims) .lt. 2) then
    error_msg("f", DIAG_SCRIPT, "", dimsizes(dims) + " dimensions, " + \
              "need 2 or 3")
  end if
  idx = ind(dims .eq. "lat")
  if (ismissing(idx)) then
    error_msg("f", DIAG_SCRIPT, "", "no lat dimension")
  end if
  idx = ind(dims .eq. "lon")
  if (ismissing(idx)) then
    error_msg("f", DIAG_SCRIPT, "", "no lon dimension")
  end if

  ; how many months included in analysis?
  months = diag_script_info@months
  dim_SEAS = dimsizes(months)

  ; we need to calculate the snow area using sftgif and sftlf
  ; read sftgif and sftlf (land sea mask)
  if (var0 .eq. "sftlf") then
    land_sea_sftlf = read_data(info_items0[imod])
    var_sftlf = var0
    info_items_sftlf = info_items0
  elseif (var1 .eq. "sftlf") then
    land_sea_sftlf = read_data(info_items1[imod])
    var_sftlf = var1
    info_items_sftlf = info_items1
  elseif (var2 .eq. "sftlf") then
    land_sea_sftlf = read_data(info_items2[imod])
    var_sftlf = var2
    info_items_sftlf = info_items2
  end if

  if (var0 .eq. "sftgif") then
    land_ice_sftgif = read_data(info_items0[imod])
    var_sftgif = var0
    info_items_sftgif = info_items0
  elseif (var1 .eq. "sftgif") then
    land_ice_sftgif = read_data(info_items1[imod])
    var_sftgif = var1
    info_items_sftgif = info_items1
  elseif (var2 .eq. "sftgif") then
    land_ice_sftgif = read_data(info_items2[imod])
    var_sftgif = var2
    info_items_sftgif = info_items2
  end if
  sftifl = land_sea_sftlf - land_ice_sftgif

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
  if (isatt(diag_script_info, "reference_run")) then
    reference_run_name = diag_script_info@reference_run
  else
    log_info("No reference run name given, set to default 'historical'")
    reference_run_name = "historical"
  end if

  attsh = True
  attsh@exp = reference_run_name
  item_hist = select_metadata_by_atts(info_items, attsh)
  dim_MOD_hist = ListCount(item_hist)

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
  atts1 = True
  atts1@exp = diag_script_info@scenarios(0)
  item_rcp26 = select_metadata_by_atts(info_items, atts1)
  dim_MOD_26 = ListCount(item_rcp26)
  ; indicate array for scenarios
  data26_rel = new((/dim_MOD_26, nyearsp1/), float, \
                   A0@_FillValue)

  if (dim_scen .gt. 1) then
    ; How many model runs in 2nd sceanrio?
    atts2 = True
    atts2@exp = diag_script_info@scenarios(1)
    item_rcp45 = select_metadata_by_atts(info_items, atts2)
    dim_MOD_45 = ListCount(item_rcp45)
    data45_rel = new((/dim_MOD_45, nyearsp1/), float, \
                     A0@_FillValue)
  end if
  if (dim_scen .gt. 2) then
    ; How many model runs in 3rd scenario?
    atts3 = True
    atts3@exp = diag_script_info@scenarios(2)
    item_rcp60 = select_metadata_by_atts(info_items, atts3)
    dim_MOD_60 = ListCount(item_rcp60)
    data60_rel = new((/dim_MOD_60, nyearsp1/), float, \
                     A0@_FillValue)
  end if
  if (dim_scen .gt. 3) then
    ; How many model runs in 4th scenario?
    atts4 = True
    atts4@exp = diag_script_info@scenarios(3)
    item_rcp85 = select_metadata_by_atts(info_items, atts4)
    dim_MOD_85 = ListCount(item_rcp85)
    data85_rel = new((/dim_MOD_85, nyearsp1/), float, \
                     A0@_FillValue)
  end if
  if (dim_scen .gt. 4) then
    error_msg("f", DIAG_SCRIPT, "", "Too many different scenarios, " + \
              "code not set up for more than 4.")
  end if
  delete([/attsh, atts1, atts2, atts3, atts4/])
  delete([/item_rcp26, item_rcp45, item_rcp60, item_rcp85/])

  i = 0
  i26 = 0
  i45 = 0
  i60 = 0
  i85 = 0
  do imod = 0, dim_MOD_hist - 1
    if (imod .eq. 0) then
      delete(A0)
    else
      A0 = read_data(info_items[imod])
      ; we need to calculate the snow area using sftgif (glacier) and sftlf
      land_ice_sftgif = read_data(info_items_sftgif[imod])
      ; read sftlf (land sea mask)
      land_sea_sftlf = read_data(info_items_sftlf[imod])
      sftifl = land_sea_sftlf - land_ice_sftgif

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
      atts = True
      atts@dataset = info_items[imod]@dataset
      atts@exp = diag_script_info@scenarios(s)
      atts@ensemble = info_items[imod]@ensemble
      item_rcp = select_metadata_by_atts(info_items, atts)  ; matching list
      if (ListCount(item_rcp) .ne. 0) then
        A1 = read_data(item_rcp[0])
        ; we need to calculate the snow area using sftgif (glacier) and sftlf
        ; these are fx files and should not change between runs, use same
        ; sftifl as for historical
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
        delete([/tmp2, tmp3, tmp4, tmp5, tmp6, sftifl_conf/])
        dim_data2g = dimsizes(data2_glob)
        if (s .eq. 0) then
          if (dim_data2g .ne. nyearsp1) then
            if (dim_data2g .gt. nyearsp1) then
              error_msg("w", DIAG_SCRIPT, "", "Length of dataset " + \
                        info_items[imod]@dataset + " in recipe is " + \
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
              error_msg("w", DIAG_SCRIPT, "", "Length of dataset " + \
                        info_items[imod]@dataset + " in recipe is " + \
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
              error_msg("w", DIAG_SCRIPT, "", "Length of dataset " + \
                        info_items[imod]@dataset + " in recipe is " + \
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
              error_msg("w", DIAG_SCRIPT, "", "Length of dataset " + \
                        info_items[imod]@dataset + " in recipe is " + \
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
      delete(item_rcp)
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
  atts_runs = True
  atts_runs@exp = reference_run_name
  item_runs = select_metadata_by_atts(info_items, atts_runs)
  runs(0, 0) = ListCount(item_runs)
  do scen = 0, dim_scen - 1
    do p = 1, nperiods - 1
      atts_runs@exp = diag_script_info@scenarios(scen)
      item_runs = select_metadata_by_atts(info_items, atts_runs)
      nrp1 = NewList("fifo")
      nrp2 = NewList("fifo")
      do r = 0, ListCount(item_runs) - 1
        if (item_runs[r]@end_year .gt. toint(diag_script_info@eyears(p))) then
          ListAppend(nrp1, item_runs[r])
        elseif (item_runs[r]@end_year .eq. \
                toint(diag_script_info@eyears(p))) then
          ListAppend(nrp2, item_runs[r])
        end if
      end do
      print(nrp1)
      tmp = ListCount(nrp1)
      ; check if only found missing, in that case set tmp to zero
      if tmp .eq. 0 then
        tmp = 0
      end if
      print(nrp2)
      tmp2 = ListCount(nrp2)
      ; check if only found missing, in that case set tmp2 to zero
      if tmp2 .eq. 0 then
        tmp2 = 0
      end if
      runs(scen + 1, p) = tmp + tmp2
      delete([/tmp, tmp2, nrp1, nrp2/])
      if (runs(scen + 1, p) .le. 4) .and. (p .ne. 0) then
        p2 = toint(eyears(p - 1)) - toint(syears(0))
        p3 = toint(eyears(p)) - toint(syears(0))
        mean_val(scen + 1, p2 : p3) = mean_val@_FillValue
        std_val(scen + 1, p2 : p3) = std_val@_FillValue
        log_info("Scenario " + diag_script_info@scenarios(scen) + \
                 " in period ending " + diag_script_info@eyears(p) + \
                 " has less than 5 models, set to missing.")
      end if
    end do
  end do

  dims_data = dimsizes(mean_val)
  to_plot = new((/3 * dims_data(0), dims_data(1)/), float)

  if (isatt(diag_script_info, "spread")) then
    spread = diag_script_info@spread
  else
    spread = 1.0
  end if

  do j = 0, dims_data(0) - 1
    to_plot(j * 3, :) = mean_val(j, :)
    to_plot((j * 3) + 1, :) = to_plot(j * 3, :) + spread * std_val(j, :)
    to_plot((j * 3) + 2, :) = to_plot(j * 3, :) - spread * std_val(j, :)
  end do
  copy_VarMeta(mean_val, to_plot)
  to_plot!0 = "metric"
  str_spread = sprintf("%5.2f", spread)
  if (dim_scen .eq. 1) then
    to_plot&metric = (/"mean", "+" + str_spread + "std", \
                       "-" + str_spread + "std", \
                       "mean", "+" + str_spread + "std", \
                       "-" + str_spread + "std"/)
  elseif (dim_scen .eq. 2) then
    to_plot&metric = (/"mean", "+" + str_spread + "std", \
                       "-" + str_spread + "std", \
                       "mean", "+" + str_spread + "std", \
                       "-" + str_spread + "std", \
                       "mean", "+" + str_spread + "std", \
                       "-" + str_spread + "std"/)
  elseif (dim_scen .eq. 3) then
    to_plot&metric = (/"mean", "+" + str_spread + "std", \
                       "-" + str_spread + "std", \
                       "mean", "+" + str_spread + "std", \
                       "-" + str_spread + "std", \
                       "mean", "+" + str_spread + "std", \
                       "-" + str_spread + "std", \
                       "mean", "+" + str_spread + "std", \
                       "-" + str_spread + "std"/)
  elseif (dim_scen .eq. 4) then
    to_plot&metric = (/"mean", "+" + str_spread + "std", \
                       "-" + str_spread + "std", \
                       "mean", "+" + str_spread + "std", \
                       "-" + str_spread + "std", \
                       "mean", "+" + str_spread + "std", \
                       "-" + str_spread + "std", \
                       "mean", "+" + str_spread + "std", \
                       "-" + str_spread + "std", \
                       "mean", "+" + str_spread + "std", \
                       "-" + str_spread + "std"/)
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
    to_plot@diag_script = array_append_record(temp, (/DIAG_SCRIPT/), 0)
    delete(temp)
  else  ; Add as new attribute
    to_plot@diag_script = (/DIAG_SCRIPT/)
  end if
  to_plot@var = var_snw  ; Overwrite existing entry
  if(isatt(variable_info[idx_snw], "long_name")) then
    to_plot@var_long_name = variable_info[idx_snw]@long_name
  end if
  if(isatt(variable_info[0], "units")) then
    to_plot@units     = variable_info[idx_snw]@units
  end if

  runs@var = "nr"
  runs!0 = "scenario"
  runs&scenario = diag_script_info@label
  runs!1 = "period"
  runs&period = diag_script_info@syears
  runs@diag_script = (/DIAG_SCRIPT/)
  runs@var_long_name = "number of model runs per scenario and period"

  ; ###########################################
  ; # Output to netCDF               #
  ; ###########################################
  ; Output (controlled by diag_script_info)
  to_plot@ncdf = variable_info[0]@diagnostic + ".nc"
  ncdf_file = work_dir + "/" + to_plot@ncdf
  ncdf_file@existing = "overwrite"
  ; Function in ~/interface_scripts/auxiliary.ncl
  ncdf_outfile = ncdf_write(to_plot, ncdf_file)
  if(isatt(diag_script_info, "model_nr") .and. \
     diag_script_info@model_nr .eq. True) then
    ncdf_file2 = work_dir + "/nr_runs_" + variable_info[0]@diagnostic + ".nc"
    ncdf_outfile2 = ncdf_write(runs, ncdf_file2)
  end if

  leave_msg(DIAG_SCRIPT, "")
end
