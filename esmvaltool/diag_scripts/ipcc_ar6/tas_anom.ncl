; #############################################################################
; PLOTS TIME SERIES OF NEAR-SURFACE AIR TEMPERATURE ANOMALY
; Authors: Lisa Bock (DLR, Germany)
; #############################################################################
;
; Description
;    This script is for plotting a time series of the global mean (anomaly)
;    surface temperature (GMST)
;
; Required diag_script_info attributes (diagnostic specific)
;    styleset: as in diag_scripts/shared/plot/style.ncl functions
;
; Optional diag_script_info attributes (diagnostic specific)
;    blending: if true, calculates blended surface temperature
;    ref_start: start year of reference period for anomalies
;    ref_end: end year of reference period for anomalies
;    ref_value: if true, right panel with mean values is attached
;    ref_mask: if true, model fields will be masked by reference fields
;    region: name of domain
;    plot_units: variable unit for plotting
;    y-min: set min of y-axis
;    y-max: set max of y-axis
;    header: if true, region name as header
;    volcanoes: if true, adds volcanoes to the plot
;    write_stat: if true, write multi model statistics in nc-file
;
; Caveats
;
; Modification history
;    20191030-bock_lisa: written.
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/scaling.ncl"
load "$diag_scripts/shared/ensemble.ncl"

load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/xy_line.ncl"

load "$diag_scripts/ipcc_ar6/percentiles.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  ; Check for required settings
  exit_if_missing_atts(diag_script_info, "styleset")

  ; Set default values for non-required diag_script_info attributes
  set_default_att(diag_script_info, "blending", False)
  set_default_att(diag_script_info, "ref_value", True)
  set_default_att(diag_script_info, "ref_mask", False)
  set_default_att(diag_script_info, "region", "Global")
  set_default_att(diag_script_info, "header", False)
  set_default_att(diag_script_info, "volcanoes", False)
  set_default_att(diag_script_info, "write_stat", False)

  var = metadata_att_as_array(variable_info, "short_name")

  info_tas = select_metadata_by_name(input_file_info, "tas")
  project0 = info_tas[0]@project
  exp0 = info_tas[0]@exp

  datasetnames = metadata_att_as_array(info_tas, "dataset")
  ensembles = metadata_att_as_array(info_tas, "ensemble")
  dim_MOD = ListCount(info_tas)

  if (any(var .eq. "tasa")) then
    info_tasa = select_metadata_by_name(input_file_info, "tasa")
  end if

  if (diag_script_info@blending) then
    info_tos = select_metadata_by_name(input_file_info, "tos")
    if (project0 .eq. "CMIP6") then
      info_sic = select_metadata_by_name(input_file_info, "siconc")
    else
      info_sic = select_metadata_by_name(input_file_info, "sic")
    end if
    info_sftlf = select_metadata_by_name(input_file_info, "sftlf")
  end if

  ; Create output directories
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)
  out_path = config_user_info@work_dir
  system("mkdir -p " + out_path)

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

  ; Determine time range
  start_year = min(metadata_att_as_array(info_tas, "start_year"))
  end_year = max(metadata_att_as_array(info_tas, "end_year"))
  all_years = ispan(start_year, end_year, 1)

  ; Create time coordinate
  ntime = end_year - start_year + 1
  time = new(ntime, integer)
  time = ispan(start_year, end_year, 1)

  anom_ref = new((/dim_MOD/), double)
  anom_ref!0 = "model"
  anom_ref&model = datasetnames

  ref_start = diag_script_info@ref_start
  ref_end = diag_script_info@ref_end

  ; Set index of the reference dataset and read it
  if (isvar("info_tasa")) then
    ref_ind = -1
    ref_dataset = metadata_att_as_array(info_tasa, "dataset")
    datasetnames := array_append_record(datasetnames, ref_dataset, 0)
    ensembles := array_append_record(ensembles, "1", 0)
    anom_ref := array_append_record(anom_ref, anom_ref@_FillValue, 0)
    anom_ref&model = datasetnames
    ListAppend(info_tas, info_tasa[0])
    tas_ref = read_data(info_tasa[0])
    ndim = dimsizes(tas_ref)
    if (ref_start .ne. 1961 .or. ref_end .ne. 1990) then
      tmp = time_operations(tas_ref, ref_start, ref_end, "average", \
                            "monthlyclim", True)
      do i = 0, ndim(0) - 1
        tas_ref(i, :, :) = tas_ref(i, :, :) - tmp(mod(i, 12), :, :)
      end do
      delete(tmp)
    end if
    delete(ndim)
    ; calculate field mean
    gmst_mon := area_operations(tas_ref, -90., 90., \
                                0., 360., "average", True)
    ; yearly average
    diag_ref := time_operations(gmst_mon, -1, -1, "average", \
                                "yearly", True)
  elseif (isatt(info_tas[0], "reference_dataset"))
    ref_dataset = metadata_att_as_array(info_tas, "reference_dataset")
    ref_ind = ind(datasetnames.eq.ref_dataset(0))
    tas_ref = read_data(info_tas[ref_ind])
  else
    error_msg("f", DIAG_SCRIPT, "", \
              "no reference dataset defined")
  end if

  ; Create model array
  if (ref_ind .ge. 0) then
    model_arr = new((/dim_MOD, ntime/), double)
  else
    model_arr = new((/dim_MOD+1, ntime/), double)
  end if
  model_arr!0 = "model"
  model_arr!1 = "time"
  model_arr&model = datasetnames
  model_arr&time = time
  model_arr@_FillValue = 1e+20
  model_arr = model_arr@_FillValue
  model_arr@units = "degC"
  model_arr@ref = ref_dataset
  delete(time)

  if (ref_ind .ge. 0) then
    model_fullspan = new(dim_MOD, logical)
  else
    model_fullspan = new(dim_MOD+1, logical)
  end if
  model_fullspan = 1

  ; get multi-model mean index
  mmm_ind = ind(datasetnames .eq. "MultiModelMean")

  if (ref_ind .eq. -1) then
    model_arr(dim_MOD, :) = (/diag_ref/)
  end if

  model_arr_tas = model_arr

  ; Loop over models
  do imod = 0, dim_MOD - 1

    log_info("Process dataset: " + datasetnames(imod))

    ; Read data
    tas0 := read_data(info_tas[imod])
    dnames = getVarDimNames(tas0)
    ndim = dimsizes(tas0)

    ; Convert units from K to degC
    if (info_tas[imod]@units .eq. "K") then
      tas0 = convert_units(tas0, "degC")
    end if

    if (diag_script_info@blending .and. \
        all(imod .ne. (/ref_ind, mmm_ind/))) then

      ; Read data
      tos0 := read_data(info_tos[imod])
      sic0 := read_data(info_sic[imod])
      sftlf0 := read_data(info_sftlf[imod])
      ndim_sftlf = dimsizes(sftlf0)

      if (info_tos[imod]@units .eq. "K") then
        tos0 = convert_units(tos0, "degC")
      end if

      ; force missing cells to open water/land and scale if stored as
      ; percentage
      sic0 = where(sic0.lt.0.0, 0.0, sic0)
      sic0 = where(sic0.gt.100.0, 0.0, sic0)
      if (max(sic0) .gt. 90.0) then
        sic0 = 0.01 * sic0
      end if

      sftof0 = 1. - sftlf0
      sftof0 = where(sftof0.lt.0.0, 0.0, sftof0)
      sftof0 = where(sftof0.gt.100.0, 0.0, sftof0)
      if (max(sftof0) .gt. 90.0) then
        sftof0 = 0.01 * sftof0
      end if

      ; combine land/ice masks
      ndim_sic = dimsizes(sic0)
      do m = 0, ndim_sic(0) - 1
        sic0(m, :, :) = (1.0 - sic0(m, :, :)) * sftof0(:, :)
      end do

      ; deal with missing tos through sic
      sic0 = where(tos0 .lt. -500.00, 0.0, sic0)
      sic0 = where(tos0 .gt. 500.00, 0.0, sic0)
      sic0 = where(tos0 .eq. tos0@_FillValue, 0.0, sic0)

      ; prepare missing
      tos0 = where(abs(tos0) .gt. 500.0, tos0@_FillValue, tos0)

    end if

    ; calculate tas anomaly for each grid point
    ; Determine start/end year
    start_year = info_tas[imod]@start_year
    end_year = info_tas[imod]@end_year
    if (start_year .gt. diag_script_info@ref_start) then
      ref_start = start_year
    else
      ref_start = diag_script_info@ref_start
    end if

    base_mon = time_operations(tas0, 1961, 1990, "average", \
                               "monthlyclim", True)

    tas0_anom := tas0
    do i = 0, ndim(0) - 1
      tas0_anom(i, :, :) = tas0(i, :, :) - base_mon(mod(i, 12), :, :)
    end do

    if (diag_script_info@blending .and. \
        all(imod .ne. (/ref_ind, mmm_ind/))) then

      ; calculate tos anomaly for each grid point
      base_mon = time_operations(tos0, 1961, 1990, "average", \
                                 "monthlyclim", True)

      tos0_anom := tos0
      do i = 0, ndim(0) - 1
        tos0_anom(i, :, :) = tos0(i, :, :) - base_mon(mod(i, 12), :, :)
      end do

      ; blend
      tas0 = tas0 * (1.0 - sic0) + tos0 * sic0
      tos0_anom = tas0_anom * (1.0 - sic0) + tos0_anom * sic0

      ; deal with any remaining NaNs
      tos0 = where(tos0.eq.tos0@_FillValue, tas0, tos0)
      tos0_anom = where(tos0_anom.eq.tos0@_FillValue, tas0_anom, tos0_anom)

    end if

    ; Determine start/end year
    start_year = info_tas[imod]@start_year
    end_year = info_tas[imod]@end_year

    ; masking anomaly with reference dataset
    if (diag_script_info@ref_mask) then
      tas_ref_mask := time_operations(tas_ref, start_year, \
                                      end_year, "extract", "", 0)
      tas0_anom = where(tas_ref_mask.eq.tas_ref_mask@_FillValue, \
                        tas_ref_mask@_FillValue, tas0_anom)
    end if

    ; calculate field mean
    tas_mon := area_operations(tas0, -90., 90., \
                               0., 360., "average", True)
    tas_mon_anom := area_operations(tas0_anom, -90., 90., \
                                    0., 360., "average", True)

    ; yearly average
    diag_tas := time_operations(tas_mon_anom, -1, -1, "average", \
                                "yearly", True)

    if (ref_start .ne. 1961 .or. ref_end .ne. 1990) then
      ref_diag_tas = avg(diag_tas(ref_start-start_year:ref_end-start_year))
      diag_tas = diag_tas - ref_diag_tas
      ref_diag_tas := time_operations(tas_mon, ref_start, ref_end, "average", \
                                      "annualclim", True)
    else
      ref_diag_tas = time_operations(tas_mon, 1961, 1990, "average", \
                                     "annualclim", True)
    end if

    ; save mean ref value
    anom_ref(imod) = ref_diag_tas

    procmod := diag_tas

    if (diag_script_info@blending .and. \
        all(imod .ne. (/ref_ind, mmm_ind/))) then

      if (diag_script_info@ref_mask) then
        ; masking with reference dataset
        tos0_anom = where(tas_ref_mask.eq.tas_ref_mask@_FillValue, \
                          tas_ref_mask@_FillValue, tos0_anom)
      end if

      ; calculate field mean
      gmst_mon := area_operations(tos0, -90., 90., \
                                  0., 360., "average", True)
      gmst_mon_anom := area_operations(tos0_anom, -90., 90., \
                                       0., 360., "average", True)

      ; yearly average
      diag := time_operations(gmst_mon_anom, -1, -1, "average", \
                              "yearly", True)

      if (ref_start .ne. 1961 .or. ref_end .ne. 1990) then
        ref_diag = avg(diag(ref_start-start_year:ref_end-start_year))
        diag = diag - ref_diag
        ref_diag := time_operations(gmst_mon, ref_start, ref_end, "average", \
                                    "annualclim", True)
      else
        ref_diag = time_operations(gmst_mon, 1961, 1990, "average", \
                                   "annualclim", True)
      end if

      ; save mean ref value
      anom_ref(imod) = ref_diag

      procmod := diag

    end if

    ; Match time coordinate
    date = procmod&year
    idx1 = ind(date(0).eq.model_arr&time)
    idx2 = ind(date(dimsizes(date) - 1).eq.model_arr&time)
    if (idx1 .ne. 0) then
      model_fullspan(imod) = 0
    end if
    if (idx2 .ne. (ntime - 1)) then
      model_fullspan(imod) = 0
    end if
    model_arr(imod, idx1:idx2) = (/procmod/)

    delete(procmod)
    delete(date)

  end do

  ; *************************************
  ; calculate means over ensemble members
  ; *************************************

  ndim := dimsizes(datasetnames)
  ensemble_assign = new(ndim, integer)
  ensemble_size = new(ndim, integer)

  enum = 1
  temp = datasetnames(0)
  do gg = 0, ndim - 1
    ; does the model name match what's in temp?
    if (temp.eq.datasetnames(gg)) then
      ensemble_assign(gg) = enum  ; if so, assign it the same number
    else
      enum = enum + 1              ; if not, assign it the next number
      ensemble_assign(gg) = enum
      temp = datasetnames(gg)
    end if
  end do
  ensemble_assign@model = datasetnames
  delete(temp)
  ensemble_assign@models = str_join(datasetnames + "/" + ensembles, ", ")

  model_arr_em = model_arr(:max(ensemble_assign)-1, :)
  model_arr_em = model_arr_em@_FillValue
  anom_ref_em = anom_ref(:max(ensemble_assign)-1)
  anom_ref_em = anom_ref_em@_FillValue
  datasetnames_em = datasetnames(:max(ensemble_assign)-1)
  model_fullspan_em =  model_fullspan(:max(ensemble_assign)-1)

  do gg = 1, max(ensemble_assign)         ; calculate ensemble means
    wind := ind(ensemble_assign.eq.gg)
    ensemble_size(wind) = dimsizes(wind)
    if (dimsizes(wind).eq.1) then
      model_arr_em(gg-1, :) = (/ model_arr(wind, :) /)
      anom_ref_em(gg-1) = (/ anom_ref(wind) /)
      datasetnames_em(gg-1) = datasetnames(wind)
      model_fullspan_em(gg-1) = model_fullspan(wind)
    else
      model_arr_em(gg-1, :) = (/ dim_avg_n(model_arr(wind, :), 0) /)
      anom_ref_em(gg-1) = (/ dim_avg_n(anom_ref(wind), 0) /)
      datasetnames_em(gg-1) = datasetnames(wind(0))
      model_fullspan_em(gg-1) = model_fullspan(wind(0))
    end if
  end do

  ; *************************************
  ; calculate multi-model mean
  ; *************************************

  if (ref_ind .ge. 0) then
    ind_mod := ind(datasetnames_em.ne."MultiModelMean" .and. \
                   model_fullspan_em .and. \
                   datasetnames_em.ne.variable_info[0]@reference_dataset)
  elseif (ref_ind .eq. -1)
    ind_mod := ind(datasetnames_em.ne."MultiModelMean" .and. \
                   model_fullspan_em .and. \
                   datasetnames_em.ne.ref_dataset)
  else
    ind_mod := ind(datasetnames_em.ne."MultiModelMean" .and. \
                   model_fullspan_em)
  end if
  model_arr(mmm_ind, :) = dim_avg_n(model_arr_em(ind_mod, :), 0)
  anom_ref(mmm_ind) = avg(anom_ref_em(ind_mod))

  ; *************************************
  ; calculate 5% and 95% quantiles
  ; *************************************

  q_05 = new(ntime, double)
  q_95 = new(ntime, double)

  do itime = 0, ntime - 1
    model_num = num(.not.ismissing(model_arr_em(ind_mod, itime)))
    if (ref_ind .ge. 0) then
      ind_m := ind(.not.ismissing(model_arr(:, itime)) .and. \
                   model_fullspan .and. \
                   datasetnames.ne."MultiModelMean" .and. \
                   datasetnames.ne.variable_info[0]@reference_dataset)
    elseif (ref_ind .eq. -1)
      ind_m := ind(.not.ismissing(model_arr(:, itime)) .and. \
                   model_fullspan .and. \
                   datasetnames.ne."MultiModelMean" .and. \
                   datasetnames.ne.ref_dataset)
    else
      ind_m := ind(.not.ismissing(model_arr(:, itime)) .and. \
                   model_fullspan .and. \
                   datasetnames.ne."MultiModelMean")
    end if
    ind_sort := dim_pqsort(model_arr(ind_m, itime), 1)
    count = 0.
    do imodel = 0, dimsizes(ind_m) - 1
      count = count + 1./(model_num*ensemble_size(ind_m(ind_sort(imodel))))
      if(count .gt. 0.95) then
        break
      end if
      if(ismissing(q_05(itime)) .and. count .gt. 0.05) then
        q_05(itime) =  model_arr(ind_m(ind_sort(imodel)), itime)
      end if
      q_95(itime) =  model_arr(ind_m(ind_sort(imodel)), itime)
    end do
  end do

  ; *************************************
  ; Convert time coordinate to years (required by the plot routine)
  ; *************************************

  xmin = min(model_arr&time)
  xmax = max(model_arr&time)
  tmp = todouble(model_arr&time)
  delete(model_arr&time)
  model_arr&time = tmp
  delete(tmp)

  ; *************************************
  ; output to NetCDF
  ; *************************************

  if (diag_script_info@blending) then
    var0 = "gmst"
  else
    var0 = "gsat"
  end if

  out_path1 = out_path + "tsline_gmst_anom_" + project0 + "_" \
                       + exp0 + ".nc"
  out_path1@existing = "overwrite"
  model_arr_em@ncdf = out_path
  model_arr_em@experiment = project0 + "_" + exp0
  model_arr_em@diag_script = DIAG_SCRIPT
  model_arr_em@var = var0
  model_arr_em@long_name = "Global Mean Surface Temperature"
  ncdf_outfile = ncdf_write(model_arr_em, out_path1)

  ; *************************************
  ; optional output with statistics
  ; *************************************

  if (diag_script_info@write_stat) then

    model_arr_stat = new((/ntime, 10/), double)
    model_arr_stat!0 = "time"
    model_arr_stat&time = model_arr&time
    do i = 0, ntime-1
      statb = stat_dispersion(model_arr_em(ind_mod, i), False)
      ; mean
      model_arr_stat(i, 0) = (/statb(0)/)
      ; standard deviation
      model_arr_stat(i, 1) = (/statb(1)/)
      ; spread
      model_arr_stat(i, 2) = (/statb(14)-statb(2)/)
      ; min
      model_arr_stat(i, 3) = (/statb(2)/)
      ; max
      model_arr_stat(i, 8) = (/statb(14)/)
      ; 5% quantile
      model_arr_stat(i, 4) = (/statb(24)/)
      ; 95% quantile
      model_arr_stat(i, 7) = (/statb(25)/)
      ; 10% quantile
      model_arr_stat(i, 5) = (/statb(3)/)
      ; 90% quantile
      model_arr_stat(i, 6) = (/statb(13)/)
    end do
    if (isatt(variable_info[0], "reference_dataset")) then
      model_arr_stat(:, 9) = (/model_arr(ref_ind, :)/)
    end if

    model_arr_stat(:, 0) = (/model_arr(mmm_ind, :)/)
    model_arr_stat(:, 4) = (/q_05/)
    model_arr_stat(:, 7) = (/q_95/)

    out_path2 = out_path + "tsline_gmst_anom_stat.nc"

    out_path2@existing = "overwrite"

    model_arr_stat@ncdf = out_path
    model_arr_stat@experiment = project0 + "_" + exp0
    model_arr_stat@diag_script = DIAG_SCRIPT
    model_arr_stat@var = var0 + "_mean"
    model_arr_stat@long_name = var0 + " mulit-model mean"
    ncdf_outfile = ncdf_write(model_arr_stat(:, 0), out_path2)

    out_path2@existing = "append"

    model_arr_stat@var = var0 + "_stddev"
    model_arr_stat@long_name = var0 + " standard deviation"
    ncdf_outfile = ncdf_write(model_arr_stat(:, 1), out_path2)

    model_arr_stat@var = var0 + "_spread"
    model_arr_stat@long_name = var0 + " spread"
    ncdf_outfile = ncdf_write(model_arr_stat(:, 2), out_path2)

    model_arr_stat@var = var0 + "_min"
    model_arr_stat@long_name = var0 + " minimum"
    ncdf_outfile = ncdf_write(model_arr_stat(:, 3), out_path2)

    model_arr_stat@var = var0 + "_5quantile"
    model_arr_stat@long_name = var0 + " 5% quantile"
    ncdf_outfile = ncdf_write(model_arr_stat(:, 4), out_path2)

    model_arr_stat@var = var0 + "_10quantile"
    model_arr_stat@long_name = var0 + " 10% quantile"
    ncdf_outfile = ncdf_write(model_arr_stat(:, 5), out_path2)

    model_arr_stat@var = var0 + "_90quantile"
    model_arr_stat@long_name = var0 + " 90% quantile"
    ncdf_outfile = ncdf_write(model_arr_stat(:, 6), out_path2)

    model_arr_stat@var = var0 + "_95quantile"
    model_arr_stat@long_name = var0 + " 95% quantile"
    ncdf_outfile = ncdf_write(model_arr_stat(:, 7), out_path2)

    model_arr_stat@var = var0 + "_max"
    model_arr_stat@long_name = var0 + " maximum"
    ncdf_outfile = ncdf_write(model_arr_stat(:, 8), out_path2)

    if (isatt(variable_info[0], "reference_dataset")) then
      model_arr_stat@var = var0 + "_ref"
      model_arr_stat@long_name = var0 + " reference dataset"
      ncdf_outfile = ncdf_write(model_arr_stat(:, 9), out_path2)
    end if

  end if

  ; *************************************
  ; plotting
  ; *************************************

  ; Define workstation
  outfile = config_user_info@plot_dir + var0 + "_" \
    + str_sub_str(diag_script_info@region, " ", "_") \
    + "_" + project0 + "_" + exp0 + "_anom_" \
    + start_year + "-" + end_year
  wks = gsn_open_wks(file_type, outfile)

  ; Set resources
  res = True
  res@trXMinF = xmin
  res@trXMaxF = xmax
  ; res@trXMaxF = xmax + 0.25 * (xmax - xmin)
  res@tmXBMode = "Explicit"
  if (xmax - xmin.gt.20) then
    res@tmXBValues = ispan(xmin, xmax, 10)
    res@tmXBLabels = ispan(xmin, xmax, 10)
    res@tmXBMinorValues =  ispan(xmin, xmax, 5)
  else
    res@tmXBValues = ispan(xmin, xmax, 5)
    res@tmXBLabels = ispan(xmin, xmax, 5)
    res@tmXBMinorValues =  ispan(xmin, xmax, 1)
  end if
  res@tmXBLabelAngleF = 45
  res@tmXBLabelJust     = "CenterRight"
  if (isatt(diag_script_info, "y_min")) then
    res@trYMinF = diag_script_info@y_min
  end if
  if (isatt(diag_script_info, "y_max")) then
    res@trYMaxF = diag_script_info@y_max
  end if

  res@tmXBMode         = "Manual"
  res@tmXBTickSpacingF = 20

  if (diag_script_info@header) then
    res@tiMainString = diag_script_info@region
  else
    res@tiMainString = "Anomaly of Near-Surface Air Temperature"
  end if

  res@tiYAxisString = "Change in Temperature (~F34~0~F~ C)"

  res0 = True
  ref_start = diag_script_info@ref_start
  ref_end = diag_script_info@ref_end
  res0@tiYAxisString = tostring(ref_start) + "-" + tostring(ref_end) \
    + " Mean Temperature (~F34~0~F~ C)"

  xy_line_anom(wks, anom_ref, model_arr, model_arr&time, \
               model_arr_stat, ref_start, ref_end, res, res0, info_tas)

  log_info(" wrote " + outfile + "." + file_type)

  ; ***************************************
  ; add meta data to plot (for reporting)
  ; ***************************************

  if (diag_script_info@region .eq. "Global") then
    domain = "global"
  else
    domain = "reg"
  end if

  statistics = "anomaly"

  caption = "Time series of the " + statistics + " for Surface Temperature " \
            + "with blending and masking as Cowtan et al., 2015."
  contrib_authors = (/"bock_lisa"/)

  ; Call provenance logger
  log_provenance(ncdf_outfile, \
                 outfile + "." + file_type, \
                 caption, \
                 statistics, \
                 domain, \
                 "times", \
                 contrib_authors, \
                 (/"eyring21ipcc"/), \
                 metadata_att_as_array(input_file_info, "filename"))

  leave_msg(DIAG_SCRIPT, "")

end
