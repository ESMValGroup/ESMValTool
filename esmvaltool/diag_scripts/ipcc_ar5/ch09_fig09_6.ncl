; #############################################################################
; DIAGNOSTIC SCRIPT for reproducing IPCC ch. 9 fig. 9.6
; Author: Bettina Gier (DLR, Germany)
; CRESCENDO project
; #############################################################################
;
; Description
;    Calculated centred pattern correlations for annual mean climatologies
;    and plots them. Like IPCC ch. 9 fig 9.6
;
; Required diag_script_info attributes (diagnostics specific)
;
; Optional diag_script_info attributes (diagnostic specific)
;
; Required variable_info attributes (variable specific)
;     none
;
; Optional variable_info attributes (variable specific)
;     none
;
; Caveats
;
;
; Modification history
;    20181012-A_gier_be: Ported to v2
;    20171020-A_laue_ax: added tags for reporting
;    20170404-A_gier_be: written.
;
; #############################################################################

load "interface_scripts/interface.ncl"
load "./diag_scripts/shared/latlon.ncl"
load "./diag_scripts/shared/statistics.ncl"

begin
  enter_msg(diag_script, "")
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (var: " + variables(0) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Check and set reference dataset
  ref_model = variable_info@reference_dataset
  if (variable_info@reference_dataset.eq."None") then
    error_msg("f", diag_script, "", "no reference dataset is specified")
  end if
  log_info("reference model = " + ref_model)
  ref_ind = get_ref_dataset_idx(input_file_info, \
                                variable_info@reference_dataset)
  if isatt(variable_info, "alternative_dataset") then
    alt_ref = True
    aref_ind = get_ref_dataset_idx(input_file_info, \
                                   variable_info@reference_dataset)
  else
    alt_ref = False
  end if

  dim_MOD = dimsizes(input_file_info@dataset)
  dim_VAR = dimsizes(variables)
  var0 = variables(0)
  field_type0 = field_types(0)

  ; alltags = array_append_record(tags, \
  ; (/"PT_other", "ST_clim", "ST_corr", "DM_global"/), 0)

end

begin
  authors = "A_gier_be"

  write_references(diag_script,                  \
                   authors,                      \  ; authors
                   (/""/),                       \  ; contributors
                   (/""/),                       \  ; diag_references
                   (/""/),                       \  ; obs_references
                   (/"P_crescendo"/))               ; proj_references

  ; Check field type
  plot_fields = (/"T2Ms"/)
  if (field_type0.ne.plot_fields) then
    error_msg("f", diag_script, "", "can't use field " + field_type0 + \
              " for this diag_script, currently only (" + \
              str_join(plot_fields, ", ") + ") supported.")
  end if
  delete(plot_fields)

  ; Output directories
  ncdf_dir = config_user_info@work_dir + "pattern_cor.nc"
  modproj_dir = config_user_info@work_dir + "modprojnames.txt"
  system("mkdir -p " + config_user_info@work_dir)
end


begin
; -----------------------------------------------------------------------------
; -------------------- Compute correlation for one var ------------------------
; -----------------------------------------------------------------------------

  ; Reference model
  mod_idx = ispan(0, dim_MOD - 1, 1)
  mod_ind_woref = mod_idx(ind(mod_idx.ne.ref_ind))
  delete(mod_idx)

  ; Make output array
  all_cor = new((/dim_MOD-1/), float)
  all_cor!0 = "models"
  all_cor&models = dataset_info@dataset(mod_ind_woref)
  ; Pass on alt models
  if isatt(variable_info, "alternative_dataset").and. \
    any(dataset_info@dataset.eq.variable_info@alternative_dataset) then
    all_cor@alt_obs =  variable_info@alternative_dataset
  else
    all_cor@alt_obs = "none"
  end if

  ; Loop over models, with ref model processed first
  model_ind = array_append_record(ref_ind, mod_ind_woref, 0)
  do iloop = 0, dim_MOD - 1
    imod = model_ind(iloop)
    log_info("Processing " + dataset_info@dataset(imod))

    ; Extract model corresponding data
    var = read_data(imod, var0, field_type0)

    ; Calculate annual mean
    data_yearly = time_operations(var, \
                                  toint(input_file_info@start_year(imod)), \
                                  toint(input_file_info@end_year(imod)), \
                                  "average", "yearly", True)

    ; Mean over the years
    mean_years = dim_avg_n_Wrap(data_yearly, 0)
    delete(data_yearly)

    ; Compute centred pattern correlation (ref_model will be first so no error)
    if dataset_info@dataset(imod).eq.ref_model then
      data_ref = mean_years
    else
      ; Since ref model processed first, move all models up one spot
      all_cor(iloop-1) = pattern_cor(data_ref, mean_years, 1.0, 0)
    end if
    delete(mean_years)
    delete(var)
  end do

  ; Write data
  all_cor@corvar = var0
  all_cor@var = "cor"
  all_cor@diag_script = diag_script
  all_cor@ncdf = ncdf_dir
  ncdf_outfile = ncdf_write(all_cor, ncdf_dir)

  ; Write list of models with project (arrays get squished as attributes)
  modnames = dataset_info@dataset(mod_ind_woref)
  projnames = dataset_info@project(mod_ind_woref)
  modproj = (/modnames, projnames/)
  asciiwrite(modproj_dir, modproj)

  leave_msg(diag_script, "")

end
