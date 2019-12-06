; #############################################################################
; mder/select_for_mder.ncl
; #############################################################################
;
; Description:
;     Collects output from other diagnostics and saves them in a temporary
;     file. It requires precalculation by other diagnostic scripts.
;
; Required diag_script_info attributes:
;     wdiag: diagnostic names.
;     domain: domain of the diagnostic (used for provenance tracking).
;     ref_dataset: reference dataset.
;
; Optional diag_script_info attributes:
;     average_ens: average over ensemble members of the same model (default:
;                  False).
;     derive_var: derive variable.
;
; Caveats:
;     "warning: in unique_labels_min (diag_scripts/shared/plot/style.ncl), Add
;     more attributes to prio to make labels unique! Continuing with non-unique
;     labels" is normal when using "average_ens = true".
;
; Modification history:
;     20191121-schlund_manuel: added new provenance tracking.
;     20180718-schlund_manuel: ported to v2.0.
;     20131203-wenzel_sabrina: written based on perfmetrics_grading.ncl.
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/set_operators.ncl"

load "$diag_scripts/shared/plot/style.ncl"

load "$diag_scripts/shared/mder.ncl"


begin

  enter_msg(DIAG_SCRIPT, "")

  if (.not. isdefined("input_file_info")) then
    error_msg("f", DIAG_SCRIPT, "", "no input_file_info available, " + \
              "variable has to be selected as ancestor in the recipe")
  end if

  ; Variable
  VAR0 = variable_info[0]
  DIM_VAR = ListCount(variable_info)
  if (DIM_VAR .gt. 2) then
    error_msg("w", DIAG_SCRIPT, "", "this diagnostic supports at most two " + \
              "variables, got " + DIM_VAR)
  end if

  ; Input data
  INFO0 = select_metadata_by_name(input_file_info, VAR0@short_name)
  DATASETS = metadata_att_as_array(INFO0, "dataset")

  ; Get correct variable
  if (isatt(diag_script_info, "derive_var")) then
    var0 = diag_script_info@derive_var
  else
    var0 = VAR0@short_name
  end if

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Ensemble averaging (if desired)
  avgens = False
  if (isatt(diag_script_info, "average_ens")) then
    avgens = diag_script_info@average_ens
  end if
  if (avgens) then
    DATASETS := get_unique_values(DATASETS)
    INFO0 := get_unique_items(INFO0, "dataset")
  else
    DATASETS := unique_labels_min(INFO0, (/"dataset", "ensemble"/))
  end if
  DIM_DAT = dimsizes(DATASETS)

end

begin

  ; Required attributes
  req_atts = (/"wdiag", "domain", "ref_dataset"/)
  exit_if_missing_atts(diag_script_info, req_atts)
  DIM_DIAG = dimsizes(diag_script_info@wdiag)

  ; Provenance
  AUTHORS = (/"wenzel_sabrina", "schlund_manuel"/)
  DOMAIN = diag_script_info@domain
  PLOT_TYPE = ""
  PLOT_FILE = "missing"
  REFERENCES = (/"wenzel16jclim"/)

  ; Get reference dataset(s)
  ref_dataset = diag_script_info@ref_dataset

end

begin

  opt = diag_script_info

  ; Iterate over all diagnostics
  do dd = 0, DIM_DIAG - 1
    ANCESTORS := new(DIM_DAT, string)
    opt@wdiag := diag_script_info@wdiag(dd)
    var0_new = opt@wdiag

    ; Caption and statistics
    if (isStrSubset(opt@wdiag, "_t")) then
      statistics := (/"trend", "stddev"/)
      caption = "Trend in " + var0 + "."
    else if(isStrSubset(opt@wdiag, "_c")) then
      statistics := (/"mean", "clim", "stddev"/)
      caption = "Climatological mean in " + var0 + "."
    else
      statistics := (/"mean", "stddev"/)
      caption = "Mean in " + var0 + "."
    end if
    end if

    ; Output file
    work_dir = config_user_info@work_dir
    system("mkdir -p " + work_dir)
    ncdf_file = work_dir + var0_new + ".nc"

    ; Read data
    do idat = 0, DIM_DAT - 1
      log_info(DATASETS(idat))
      info = INFO0[idat]

      ; Get input directory (using ancestor diagnostic)
      input_files := diag_script_info@input_files
      var_files = str_match(input_files, "_info.ncl")
      do ivar = 0, dimsizes(var_files) - 1
        input_files := input_files(ind(input_files .ne. var_files(ivar)))
      end do
      input_dir = input_files
      if (dimsizes(input_dir) .gt. 1) then
        error_msg("f", DIAG_SCRIPT, "", "multiple input directories " + \
                  "(ancestors in recipe) given, only one supported")
      end if

      ; Check if dataset is available
      new_path = input_dir + "/" + diag_file_name(info@filename, opt)
      if (fileexists(new_path)) then
        log_info("Reading " + new_path)
        data_temp = ncdf_read(new_path, var0)
        ANCESTORS(idat) = new_path

        ; Define data array in the first iteration
        if (.not. isdefined("data")) then
          dim_data = array_append_record((/DIM_DAT/), dimsizes(data_temp), 0)
          data = new(dim_data, float)
          data!0 = "dataset"
          data&dataset = DATASETS
        end if
        rank = dimsizes(dim_data)

        ; Add data
        if (rank .eq. 2) then
          data(idat, :) = data_temp
        else if (rank .eq. 3) then
          data(idat, :, :) = data_temp
        else if (rank .eq. 4) then
          data(idat, :, :, :) = data_temp
        else if (rank .eq. 5) then
          data(idat, :, :, :, :) = data_temp
        else
          error_msg("f", DIAG_SCRIPT, "", "rank " + rank + \
                    " data not supported")
        end if
        end if
        end if
        end if
        delete(data_temp)

      ; Skip dataset if no data is available
      else
        log_info("No data available for dataset: " + DATASETS(idat))
        log_info("Skipping")
      end if
    end do

    ; Define result variable (always one diagnostic)
    val = new((/1, DIM_DAT, 2/), float)
    val!0 = "diagnostics"
    val!1 = "datasets"
    val!2 = "stat"
    val&diagnostics = var0_new
    val&datasets = DATASETS
    val&stat = (/"mean", "stddev"/)
    val@ncdf_file = ncdf_file

    if (dimsizes(ref_dataset) .gt. 1) then
      val@reference = (/ref_dataset(dd)/)
    else
      val@reference = (/ref_dataset/)
    end if

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
        error_msg("f", DIAG_SCRIPT, "", "rank " + rank + " data not supported")
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
    val@diag_script = (/DIAG_SCRIPT/)
    val@var = "diag"

    ; Write NetCDF output
    ncdf_outfile = ncdf_write(val, val@ncdf_file)
    delete(val)

    ; Provenance tracking
    log_provenance(ncdf_outfile, PLOT_FILE, caption, statistics, DOMAIN, \
                   PLOT_TYPE, AUTHORS, REFERENCES, ANCESTORS)
  end do

  leave_msg(DIAG_SCRIPT, "")

end
