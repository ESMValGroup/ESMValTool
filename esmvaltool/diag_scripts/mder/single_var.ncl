; #############################################################################
; mder/single_var.ncl
; #############################################################################
; DIAGNOSTIC SCRIPT FOR SIMPLE MDER OUTPUT FOR ONE VARIABLE
; Author: Manuel Schlund (DLR, Germany)
; CRESCENDO project
; #############################################################################
;
; Description:
;     Write quantities of a single variable to netcdf for MDER analysis.
;
; Required diag_script_info attributes:
;     calc_type: calculation type ('mean', 'trend').
;     region: considered region.
;     time_oper: operation for time aggregation.
;     time_opt: option for time aggregation.
;     area_oper: operation for spacial aggregation.
;     wdiag: diagnostic name for MDER output.
;
; Optional diag_script_info attributes:
;     wdiag_title: this name will apper in the figure title.
;
; Modification history:
;     20180906_A-schl_ma: written.
;
; #############################################################################

load "./interface_scripts/interface.ncl"

load "./diag_scripts/shared/latlon.ncl"
load "./diag_scripts/shared/scaling.ncl"
load "./diag_scripts/shared/statistics.ncl"

load "./diag_scripts/shared/mder.ncl"


begin

  enter_msg(diag_script, "")
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (var: " + variables(0) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Variable
  var0 = variables(0)
  field_type0 = field_types(0)

  ; Check field type
  plot_fields = (/"T2Ms"/)
  ind_f = any(field_type0 .eq. plot_fields)
  if (all(ismissing(ind_f))) then
    error_msg("f", diag_script, "", "cannot use field " + field_type0)
  end if
  delete([/plot_fields, ind_f/])

  ; Write references (TODO)
  ; write_references(diag_script, \
  ;                  (/"A-schl_ma"/), \
  ;                  (/""/), \
  ;                  (/""/), \
  ;                  (/""/), \
  ;                  (/"P_crescendo"/))

  ; ---------------------------------------------------------------------------
  ; Read recipe and config data
  ; ---------------------------------------------------------------------------

  ; Check required diag_script_info attributes
  req_atts = (/"calc_type", "region", "time_oper", "time_opt", "area_oper", \
              "wdiag"/)
  exit_if_missing_atts(diag_script_info, req_atts)
  calc_type = diag_script_info@calc_type
  wdiag = diag_script_info@wdiag
  delete(req_atts)

  ; Check calc_type
  if (.not. any(calc_type .eq. (/"mean", "trend"/))) then
    error_msg("f", diag_script, "", "calc_type " + calc_type + " not " + \
              "supported yet")
  end if

  ; Check size of wdiag
  if (dimsizes(wdiag) .gt. 1) then
    error_msg("f", diag_script, "", "only one diagnostic (wdiag) supported")
  end if

  ; Set region
  region = select_region(diag_script_info@region)

  ; Array sizes
  DIM_DAT = dimsizes(input_file_info@dataset)

  ; ---------------------------------------------------------------------------
  ; Read and preprocess data
  ; ---------------------------------------------------------------------------

  ; Loop over all datasets
  do idat = 0, DIM_DAT - 1
    dataset_name = input_file_info@dataset(idat)
    log_debug("Processing " + dataset_name)

    ; Read data
    var = read_data(idat, var0, field_type0)
    if (typeof(var) .eq. "double") then
      var := dble2flt(var)
    end if
    var@_FillValue = 1.0e20

    ; Temporal aggregation
    tmp = time_operations(var, \
                          toint(input_file_info@start_year(idat)), \
                          toint(input_file_info@end_year(idat)), \
                          diag_script_info@time_oper, \
                          diag_script_info@time_opt, \
                          True)

    ; Special case for lai:
    ; Since datasets have different land covers, all sparse vegetated points
    ; must be masked
    if (isStrSubset(var0, "lai")) then
      tmp = where(tmp .ge. 0.025, tmp, tmp@_FillValue)
    end if

    ; Calculate trend if desired
    if (calc_type .eq. "trend") then
      time = tmp&$tmp!0$
      rc = regCoef_n(1.0 * time, tmp, 0, 0)
      copy_VarAtts(tmp, rc)
      rc!0 = "lat"
      rc!1 = "lon"
      rc&lat = tmp&lat
      rc&lon = tmp&lon
      tmp := rc
      delete(rc)
    end if

    ; Spatial aggregation
    area_op = area_operations(tmp, \
                              region(0), \
                              region(1), \
                              region(2), \
                              region(3), \
                              diag_script_info@area_oper, \
                              True)

    ; Convert units if necessary
    if (isatt(variable_info, "plot_units")) then
      new_units = variable_info@plot_units
      if (calc_type .eq. "trend" .and. new_units .eq. "degC") then
        area_op@units = "degC"
      else
        area_op = convert_units(area_op, new_units)
      end if
    end if

    ; Units conversion for trend
    if (calc_type .eq. "trend") then
      if (input_file_info@mip(0) .eq. "Lmon") then
        temp_unit = "y-1"
      else
        area_op = area_op * 10.0
        temp_unit = "decade-1"
      end if
      if (area_op@units .ne. "1") then
        area_op@units = area_op@units + " " + temp_unit
      else
        area_op@units = temp_unit
      end if
      delete(temp_unit)
    end if

    ; Assign to global array
    if (.not. isvar("all_data")) then
      if (dimsizes(dimsizes(area_op)) .eq. 1 .and. \
          all(dimsizes(area_op) .eq. 1)) then
        dim_data = (/DIM_DAT/)
      else
        dim_data = array_append_record(DIM_DAT, dimsizes(area_op), 0)
      end if
      RANK = dimsizes(dim_data)
      all_data = new(dim_data, float)
      all_data!0 = "dataset"
      all_data&dataset = input_file_info@dataset
      copy_VarCoords_1(area_op, all_data)
    end if
    if (RANK .eq. 1) then
      all_data(idat) = area_op
    else if (RANK .eq. 2) then
      all_data(idat, :) = area_op
    else if (RANK .eq. 3) then
      all_data(idat, :, :) = area_op
    else
      all_data(idat, :, :, :) = area_op
    end if
    end if
    end if
    delete([/var, tmp, area_op/])
  end do

  ; ---------------------------------------------------------------------------
  ; Write NETCDF file for MDER
  ; ---------------------------------------------------------------------------

  output = mder_output_desired()
  if (output) then
    work_dir = output@work_dir
    opt = diag_script_info

    ; Iterate over all datasets
    do idat = 0, DIM_DAT - 1
      log_debug("Processing " + input_file_info@dataset(idat))
      new_path = work_dir + diag_file_name(idat, opt)

      ; Write files
      if (RANK .eq. 1) then
        val = (/all_data(idat), 0.0/)
      else if (RANK .eq. 2) then
        zeros = new(dimsizes(all_data(idat, :)), float)
        zeros = 0.0
        val = (/all_data(idat, :), zeros/)
      else if (RANK .eq. 3) then
        zeros = new(dimsizes(all_data(idat, :, :)), float)
        zeros = 0.0
        val = (/all_data(idat, :, :), zeros/)
      else
        zeros = new(dimsizes(all_data(idat, :, :, :)), float)
        zeros = 0.0
        val = (/all_data(idat, :, :, :), zeros/)
      end if
      end if
      end if
      val!0 = "case"
      val&case = (/"val", "stddev"/)
      dims = getvardims(all_data)
      do idim = 1, dimsizes(dims) - 1
        val!idim = dims(idim)
        val&$dims(idim)$ = all_data&$dims(idim)$
      end do
      val@units = all_data@units
      val@var = var0
      val@diag_script = (/diag_script/)
      copy_VarAtts(diag_script_info, val)
      val@wdiag := opt@wdiag
      val@dataset = input_file_info@dataset(idat)
      val@start_year = input_file_info@start_year(idat)
      val@end_year = input_file_info@end_year(idat)
      val@ncdf = new_path
      ncdf_outfile = ncdf_write(val, new_path)
      delete([/val, new_path/])
    end do
  end if
  delete(all_data)

  leave_msg(diag_script, "")

end
