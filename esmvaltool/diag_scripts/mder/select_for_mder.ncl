; #############################################################################
; mder/select_for_mder.ncl
; #############################################################################
; DIAGNOSTIC SCRIPT SELECTING MDER DIAGNOSTICS
; Author: Sabrina Wenzel (DLR, Germany)
; EMBRACE project PIs: Veronika Eyring (DLR, Germany)
;                      Klaus Wyser (SMHI, Sweden)
; #############################################################################
; Description
;     Collects output from other diagnostics and saves them in a temporary
;     file. It requires precalculation by other diagnostis scripts.
;
; Required recipe attributes:
;     region: selected region
;     time_avg: type of time average
;     wdiag: Diagnostic names
;
; Optional recipe attributes:
;
; Caveats
;
; Modification history
;   * 20180718-A_schl_ma: ported to v2.0
;   * 20131203-A_wenz_sa: written based on perfmetrics_grading.ncl
;
; #############################################################################

load "interface_scripts/interface.ncl"

load "./diag_scripts/shared/latlon.ncl"
load "./diag_scripts/shared/statistics.ncl"
load "./diag_scripts/shared/set_operators.ncl"

load "./diag_scripts/shared/mder.ncl"


begin

  if (.not. isdefined("dataset_info")) then
    error_msg("f", diag_script, "", "no dataset_info available, " + \
              "variable has to be selected as ancestor in the recipe")
  end if

  ; Get correct variable and field type
  if (isatt(diag_script_info, "derive_var")) then
    var0 = diag_script_info@derive_var
  else
    var0 = variables(0)
  end if
  field_type0 = field_types(0)

  enter_msg(diag_script, "")
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Get datasets (consider ensemble averaging if desired)
  avgens = False
  if (isatt(diag_script_info, "average_ens")) then
    avgens = diag_script_info@average_ens
  end if
  if (avgens .eq. "True") then
    idx_datasets = UNIQ(dataset_info@dataset)
    dataset_names = dataset_info@dataset(idx_datasets)
  else
    idx_datasets = ispan(0, dimsizes(dataset_info@dataset) - 1, 1)
    dataset_names = unique_labels_min((/"dataset", "ensemble"/))
  end if

  ; Array sizes
  DIM_DAT = dimsizes(idx_datasets)
  DIM_VAR = dimsizes(variables)
  DIM_DIAG = dimsizes(diag_script_info@wdiag)

end

begin

  ; Check field type
  if (.not. any(field_type0 .eq. (/"T3M", "T2Ms", "T2Mz", "T0D"/))) then
    error_msg("f", diag_script, "", "cannot use field " + field_type0)
  end if

  ; Get reference dataset
  ref_dataset = diag_script_info@ref_dataset(0)
  if (dimsizes(diag_script_info@ref_dataset) .gt. 1) then
      alt_ref_dataset = diag_script_info@ref_dataset(1)
  end if

end

begin

  opt = True

  ; Iterate over all diagnostics
  do dd = 0, DIM_DIAG - 1
    opt@wdiag = diag_script_info@wdiag(dd)
    var0_new = opt@wdiag

    ; Output file
    work_dir = config_user_info@work_dir
    system("mkdir -p " + work_dir)
    ncdf_file = work_dir + var0_new + ".nc"

    ; Read data
    do i = 0, DIM_DAT - 1
      idat = idx_datasets(i)

      ; Get input directory (using ancestor diagnostic)
      input_files = diag_script_info@input_files
      input_dir = input_files(ind(input_files .ne. \
                                  str_match(input_files, "_info.ncl")))
      if (dimsizes(input_dir) .gt. 1) then
        error_msg("f", diag_script, "", "multiple input directories " + \
                  "(ancestors in recipe) given, only one supported")
      end if

      ; Check if dataset is available
      new_path = input_dir + "/" + diag_file_name(idat, opt)
      if (fileexists(new_path)) then
        log_info("Reading " + new_path)
        data_temp = ncdf_read(new_path, var0)

        ; Define data array in the first iteration
        if (.not. isdefined("data")) then
          dim_data = array_append_record((/DIM_DAT/), dimsizes(data_temp), 0)
          data = new(dim_data, float)
          data!0 = "dataset"
          data&dataset = dataset_names
        end if
        rank = dimsizes(dim_data)

        ; Add data
        if (rank .eq. 2) then
          data(i, :) = data_temp
        else if (rank .eq. 3) then
          data(i, :, :) = data_temp
        else if (rank .eq. 4) then
          data(i, :, :, :) = data_temp
        else if (rank .eq. 5) then
          data(i, :, :, :, :) = data_temp
        else
          error_msg("f", diag_script, "", "rank " + rank + \
                    " data not supported")
        end if
        end if
        end if
        end if
        delete(data_temp)

      ; Skip dataset if no data is available
      else
        log_info("No data available for dataset: " + \
                 dataset_info@dataset(idat))
        log_info("Skipping")
      end if
    end do

    ; Define result variable
    ; FIXME: do not need for DIM_VAR dimension here, as
    ;        we always use 1 variable
    ; DIM_REF = dimsizes(diag_script_info@ref_dataset)
    ; val = new((/DIM_VAR, DIM_DAT , DIM_REF, 2/), float)
    val = new((/DIM_VAR, DIM_DAT, 2/), float)
    val!0 = "diagnostics"
    val!1 = "datasets"
    val!2 = "stat"
    val&diagnostics = var0_new
    val&datasets = dataset_names
    if (isdefined("alt_ref_idx")) then
      val@reference = (/ref_dataset, alt_ref_dataset/)
    else
      val@reference = (/ref_dataset/)
    end if
    val&stat = (/"mean", "stddev"/)
    val@ncdf_file = ncdf_file

    ; Loop over all datasets
    do idat = 0, DIM_DAT - 1

      ; Extract data of given dataset
      if (rank .eq. 2) then
        var = data(idat, 0)
        var_std = data(idat, 1)
      else if (rank .eq. 3) then
        var = data(idat, :, 0)
        var_std = data(idat, :, 1)
      else if (rank .eq. 4) then
        var = data(idat, :, :, 0)
        var_std = data(idat, :, :, 1)
      else if (rank .eq. 5) then
        var = data(idat, :, :, :, 0)
        var_std = data(idat, :, :, :, 1)
      else
        error_msg("f", diag_script, "", "rank " + rank + " data not supported")
      end if
      end if
      end if
      end if
      val(0, idat,  0) = var
      val(0, idat,  1) = var_std
      delete(var)
    end do
    delete(data)

    ; Attach attributes to the results
    val@title = "selected diagnostic"
    val@long_name = "Diagnostic for Multiple Diagnostic Ensemble Regression"
    val@diag_script = (/diag_script/)
    val@var = "diagn"

    ; Write NetCDF output
    ncdf_outfile = ncdf_write(val, val@ncdf_file)
    delete(val)
  end do

  leave_msg(diag_script, "")

end
