; ############################################################################
; carbon_corr_2var.ncl
; Author: Sabrina Wenzel (DLR, Germany)
; ############################################################################
; Description
;    Reads and correlates in time spatial averages of CO2 fluxes (derived
;    from land and ocean carbon fluxes) and temperature
;
; Required diag_script_info attributes (diagnostics specific)
;    * ec_anom      (optional) plot anomalies.
;
; Caveats
;
; Modification history
;     20200115-zechlau_sabrina: code rewritten for ESMValTool v2.0
;     20151027-lauer_axel: added call to 'write_references'
;     201404??-zechlau_sabrina: updated according to changes in ESMValTool
;     201303??-zechlau_sabrina: written
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/carbon_ec/carbon_aux.ncl"

load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/scaling.ncl"
load "$diag_scripts/shared/ensemble.ncl"

load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/scatterplot.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  tmp = metadata_att_as_array(input_file_info, "short_name")
  variables = get_unique_values(tmp)
  delete(tmp)

  ; get indices associated with specific variables
  tas_idx = ind(variables.eq."tas")
  nbp_idx = ind(variables.eq."nbp")
  fgco2_idx = ind(variables.eq."fgco2")

  ; Load var_info
  var0 = variables(tas_idx)
  info0 = select_metadata_by_name(input_file_info, var0)

  var1 = variables(nbp_idx)
  info1 = select_metadata_by_name(input_file_info, var1)

  var2 = variables(fgco2_idx)
  info2 = select_metadata_by_name(input_file_info, var2)

  log_info("+++++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + " + " + var1 + ")")
  log_info("+++++++++++++++++++++++++++++++++++++++++++++")

  dim_MOD = ListCount(info0)
  datasetnames0 = metadata_att_as_array(info0, "dataset")
  datasetnames1 = metadata_att_as_array(info1, "dataset")
  datasetnames2 = metadata_att_as_array(info2, "dataset")

  ; *****************************************************************
  ; Create output plot directory
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)
  system("mkdir -p " + config_user_info@work_dir)

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

  ; Check for required settings
  exit_if_missing_atts(diag_script_info, "styleset")

  ; Set default values for non-required diag_script_info attributes
  set_default_att(diag_script_info, "ts_anom", False)
  set_default_att(diag_script_info, "volcanoes", False)

  if(isatt(diag_script_info, "ec_anom")) then
    anom = diag_script_info@ec_anom
  end if

  ; ==================================================================
  ; time span
  xMax_year = max(metadata_att_as_array(info0, "end_year"))
  xMin_year = min(metadata_att_as_array(info0, "start_year"))
  nMax = xMax_year - xMin_year + 1

  ; setup data holders for time lines
  aY1reg = new((/dim_MOD, nMax/), "float")
  aY2reg = new((/dim_MOD, nMax/), "float")
  aY1reg!1 = "year"
  aY1reg&year = ispan(xMin_year, xMax_year, 1)
  aY2reg!1 = "year"
  aY2reg&year = ispan(xMin_year, xMax_year, 1)

  ; setup data holders for slope
  aYreg  = new((/dim_MOD, nMax/), "float")
  arc    = new((/2, dim_MOD/), "float")
  ccr    = new((/dim_MOD/), "float")
  var    = new((/2, dim_MOD/), "float")

  ; ==================================================================
  ; loop for models
  do imod = 0, dim_MOD-1

    ; Read data
    A0 = read_data(info0[imod])
    A1 = read_data(info1[imod])
    A2 = read_data(info2[imod])

    idx_mod := imod

    ; special case for GCP data
    if (datasetnames1(imod).eq."GCP") then
      A1 = A1*148300000000000.  ; A1@area
    end if
    if (datasetnames2(imod).eq."GCP") then
      A2 = A2*360000000000000.  ; A2@area
    end if

    ; calculate spatial average
    if (dimsizes(dimsizes(A1)).gt.1) then
      A1_tmp = dim_sum_Wrap(dim_sum_Wrap(A1))
      A1_tmp!0 = "time"
      A1_tmp&time = (/A1&time/)
    else
      A1_tmp = A1
    end if

    if (dimsizes(dimsizes(A2)).gt.1) then
      A2_tmp = dim_sum_Wrap(dim_sum_Wrap(A2))
      A2_tmp!0 = "time"
      A2_tmp&time = (/A2&time/)
    else
      A2_tmp = A2
    end if

    ; align timeseries first (there might be better ways to do it)
    if (dimsizes(A1_tmp).ne.dimsizes(A2_tmp))then
      maxt = max((/dimsizes(A1_tmp), dimsizes(A2_tmp)/))
      tmpA = new((/2, maxt/), typeof(A1_tmp&time))
      tmpA(0, :dimsizes(A1_tmp)-1) = A1_tmp&time
      tmpA(1, :dimsizes(A2_tmp)-1) = A2_tmp&time
      dum = tmpA(1, :) - tmpA(0, :)
      ind_t = ind(dum.eq.0)
      delete([/maxt, dum, tmpA/])
    else
      ind_t = ispan(0, dimsizes(A1_tmp)-1, 1)
    end if

    ; add land and ocean carbon fluxes
    A1_co2flux = A1_tmp(ind_t) + A2_tmp(ind_t)
    copy_VarAtts(A2, A1_co2flux)
    copy_VarCoords(A2, A1_co2flux)
    delete([/A1, A1_tmp, A2_tmp, ind_t/])

    ; Convert units if appropriate
    if (isatt(diag_script_info, "plot_units_y")) then
      A1 := convert_units(A1_co2flux, diag_script_info@plot_units_y)
    end if

    ; Copute time average
    tmp_A0 := time_operations(A0, -1, -1, "average", "yearly", True)
    tmp_A1 := time_operations(A1, -1, -1, "average", "yearly", False)
    delete([/A0, A1, A2, A1_co2flux/])

    ; set years with big vulcanic eruptions to _FillValues
    if (diag_script_info@ec_volc .eq. True) then

      ; remove 2 yrs after Krakatoa, Santa Maria, Agung, El Chichon, Pinatubo
      volc_years = (/1883, 1884, 1902, 1903, 1963, \
                    1964, 1982, 1983, 1991, 1992/)

      do y = 0, dimsizes(tmp_A0&year)-1
        tmp0 = where(any(tmp_A0&year(y).eq.volc_years), \
                     tmp_A0@_FillValue, tofloat(tmp_A0(y)))
        aY1reg(idx_mod, y) = tofloat(tmp0)

        tmp1 = where(any(tmp_A1&year(y).eq.volc_years), \
                     tmp_A1@_FillValue, tofloat(tmp_A1(y)))
        aY2reg(idx_mod, y) = tofloat(tmp1)
        delete([/tmp0, tmp1/])
      end do
    else
      aY1reg(idx_mod, 0:dimsizes(tmp_A0&year)-1) = tofloat(tmp_A0)
      aY2reg(idx_mod, 0:dimsizes(tmp_A1&year)-1) = tofloat(tmp_A1)
    end if
    delete([/tmp_A0, tmp_A1/])

    ; remove 11-yr running mean for anomalies
    if (anom .eq. True) then
      nmeanh = 5
      tmp1 = new((/nMax/), "float")
      tmp2 = new((/nMax/), "float")
      do n = nmeanh, nMax-nmeanh-1
        tmp1(n) = (/avg(aY1reg(idx_mod, n-nmeanh:n+nmeanh))/)
        tmp2(n) = (/avg(aY2reg(idx_mod, n-nmeanh:n+nmeanh))/)
      end do

      tmp1(:nmeanh-1) = tmp1(nmeanh)
      tmp1(nMax-nmeanh:) = (/tmp1(nMax-nmeanh-1)/)
      aY1reg(idx_mod, :) = aY1reg(idx_mod, :) - tmp1
      delete(tmp1)

      tmp2(:nmeanh-1) = tmp2(nmeanh)
      tmp2(nMax-nmeanh:) = (/tmp2(nMax-nmeanh-1)/)
      aY2reg(idx_mod, :) = aY2reg(idx_mod, :) - tmp2
      delete(tmp2)
    end if

    ; check time range
    tmin := ind(aY1reg&year .eq. max((/diag_script_info@start_year, \
                                     toint(info0[imod]@start_year)/)))
    tmax := ind(aY1reg&year .eq. min((/diag_script_info@end_year, \
                                     toint(info0[imod]@end_year)/)))

    ; calculate regression line:
    rc := regline(aY1reg(idx_mod, tmin:tmax), aY2reg(idx_mod, tmin:tmax))
    df := (/rc@nptxy-2/)

    aYreg(idx_mod, tmin:tmax) = (/rc * aY1reg(idx_mod, tmin:tmax) \
                                 + rc@yintercept/)
    arc(0, idx_mod) = (/rc/)
    arc(1, idx_mod) = (/rc@rstd/)
    ccr(idx_mod)    = (/escorc(aY1reg(idx_mod, tmin:tmax), \
                               aY2reg(idx_mod, tmin:tmax))/)
    var(0, idx_mod) = (/stddev(aY1reg(idx_mod, tmin:tmax))/)
    var(1, idx_mod) = (/stddev(aY2reg(idx_mod, tmin:tmax))/)

    delete([/rc, df, tmin, tmax/])
  end do     ; model loop

  ; counting ref_models only as one model
  datasetnames0 = where(datasetnames0.eq.variable_info[0]@reference_dataset, \
                        "OBS", datasetnames0)

  ; Get environment variables for plotting
  plot_dir = getenv("ESMValTool_plot_dir")
  diag_script_base = DIAG_SCRIPT
  output_dir = config_user_info@plot_dir

  ; define arrays for plotting
  f  = sqrt(dim_MOD)
  mm = floattointeger(ceil(f))
  colors   = project_style(info0, diag_script_info, "colors")
  thicks   = project_style(info0, diag_script_info, "thicks")
  plot_num = integertochar(ispan(97, 122, 1))
  plot     = new(dim_MOD, "graphic")

  ; Call plot scripts
  plot_file = "corr_" + var0 + "-" + var1 + "_anom_" + \
              diag_script_info@start_year + "-" + diag_script_info@end_year
  wks       = gsn_open_wks(file_type, output_dir + plot_file)

  ; Create panel for each model
  do imod = 0, dim_MOD-1

      MainStg  = plot_num(imod) + ") " + datasetnames0(imod) \
                 + "  ~F8~g~F21~~B~IAV~N~= " \
                 + sprintf("%4.2f", arc(0, imod)) + " GtC/K/yr;  r= " \
                 + sprintf("%4.2f", ccr(imod))
      XStg     = "Anomaly in Trop. Temp. (" \
                 + diag_script_info@plot_units_x + ")"
      YStg     = "Anomaly in d(land+ocean CO~B~2~N~ flux)/dt (" \
                 + diag_script_info@plot_units_y + ")"

      data_arr = new((/2, dimsizes(aY1reg(imod, :))/), float)
      data_arr!0      = "datasets"
      data_arr!1      = "time"
      data_arr&datasets = (/XStg, YStg/)
      data_arr(0, :) = (/aY1reg(imod, :)/)
      data_arr(1, :) = (/aY2reg(imod, :)/)
      data_arr@diag_script         = DIAG_SCRIPT
      data_arr@res_tiMainString    = MainStg
      data_arr@res_gsnCenterString = ""
      data_arr@res_tiMainFontHeightF = 0.025
      data_arr@res_xyMarkerColors  = "white"
      data_arr@res_trXMinF         = min(aY1reg(imod, :))
      data_arr@res_trXMaxF         = max(aY1reg(imod, :))
      data_arr@res_trYMinF         = min(aY2reg(imod, :))
      data_arr@res_trYMaxF         = max(aY2reg(imod, :))

      if (.not.all(ismissing(data_arr))) then
          plot(imod) = scatterplot(wks, data_arr, var0, False, \
                                   False, input_file_info)

          ; Add year strings as markers
          add_num_markers(wks, plot(imod), False, aY1reg(imod, :), \
                          aY2reg(imod, :))
          lineres    = True
          lineres@gsLineColor = colors(imod)
          str        = unique_string("line")
          plot@$str$ = gsn_add_polyline(wks, plot(imod), aY1reg(imod, :), \
                                        aYreg(imod, :), lineres)
      end if
      delete([/data_arr, MainStg/])
  end do

  ; Common resources for panel plots
  pan             = True
  pan@gsnFrame    = False
  pan@gsnMaximize = True
  pan@gsnPanelXWhiteSpacePercent = 10
  pan@gsnPanelYWhiteSpacePercent = 5
  gsn_panel(wks, plot, (/mm, mm/), pan)
  frame(wks)

  delete([/plot, pan/])

  ; Save data for later use
  if (diag_script_info@write_netcdf) then
    new_path = config_user_info@run_dir
    new_path = new_path + "gIAV_" \
                        + diag_script_info@start_year + "-" \
                        + diag_script_info@end_year + ".nc"

    ; Attach attributes to the results
    arc!0        = "case"
    arc&case     = (/"mean", "stddev"/)
    arc!1        = "model"
    arc&model    = datasetnames0
    arc@ncdf     = new_path
    arc@var      = "gIAV"
    arc@diag_script = (/DIAG_SCRIPT/)
    arc@varIAV   = (/var/)
    arc@rIAV     = (/ccr/)

    ; Write NetCDF output
    ncdf_outfile = ncdf_write(arc, new_path)
  end if

  ; Call provenance logger
  log_provenance(ncdf_outfile, \
                 output_dir + "." + file_type, \
                 XStg + " vs " + YStg, \
                 (/"anomaly", "corr", "stddev"/), \
                 (/"trop", "global"/),\
                 (/""/), \
                 (/"wenzel_sabrina"/), \
                 (/"wenzel14jgr"/), \
                 metadata_att_as_array(info0, "filename"))

end
; #############################################################
