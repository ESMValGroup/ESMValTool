; #############################################################################
; DIAGNOSTIC SCRIPT for reproducing IPCC ch. 9 fig. 9.6
; Author: Bettina Gier (Uni Bremen & DLR, Germany)
; CRESCENDO project
; #############################################################################
;
; Description
;    Calculates centred pattern correlations for annual mean climatologies
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
;    20190205-gier_bettina: Adapted to new ncl structure
;    20181012-gier_bettina: Ported to v2
;    20171020-lauer_axel: added tags for reporting
;    20170404-gier_bettina: written.
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"
load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/statistics.ncl"

begin
  enter_msg(DIAG_SCRIPT, "")

  var0 = variable_info[0]@short_name
  info_items = select_metadata_by_name(input_file_info, var0)
  infiles = metadata_att_as_array(info_items, "filename")
  datasetnames = metadata_att_as_array(info_items, "dataset")
  projectnames = metadata_att_as_array(info_items, "project")
  dim_MOD = ListCount(info_items)
  dim_VAR = ListCount(variable_info)

  ; Check and set reference dataset
  ref_model = variable_info[0]@reference_dataset
  if (variable_info[0]@reference_dataset.eq."None") then
    error_msg("f", DIAG_SCRIPT, "", "no reference dataset is specified")
  end if
  log_info("reference model = " + ref_model)
  ref_ind = ind(datasetnames.eq.ref_model)
  if isatt(variable_info[0], "alternative_dataset") then
    alt_ref = True
    aref_ind = ind(datasetnames.eq.variable_info[0]@alternative_dataset)
  else
    alt_ref = False
  end if

end

begin
  ; Output directories
  ncdf_dir = config_user_info@work_dir + "pattern_cor.nc"
  modproj_dir = config_user_info@work_dir + "modprojnames.txt"
  system("mkdir -p " + config_user_info@work_dir)
end

begin
; -----------------------------------------------------------------------------
; -------------------- Compute correlation for one var ------------------------
; -----------------------------------------------------------------------------

  ; Save list of preproc files for provenance in collect.ncl
  preproc_files = metadata_att_as_array(info_items, "filename")

  ; Reference model
  mod_idx = ispan(0, dim_MOD - 1, 1)
  mod_ind_woref = mod_idx(ind(mod_idx.ne.ref_ind))
  delete(mod_idx)

  ; Make output array
  all_cor = new((/dim_MOD-1/), float)
  all_cor!0 = "models"
  all_cor&models = datasetnames(mod_ind_woref)
  ; Pass on alt models
  if alt_ref then
    all_cor@alt_obs =  variable_info[0]@alternative_dataset
  else
    all_cor@alt_obs = "none"
  end if

  ; Loop over models, with ref model processed first
  model_ind = array_append_record(ref_ind, mod_ind_woref, 0)
  do iloop = 0, dim_MOD - 1
    imod = model_ind(iloop)
    log_info("Processing " + datasetnames(imod))

    ; Extract model corresponding data
    var = read_data(info_items[imod])

    ; Calculate annual mean
    data_yearly = time_operations(var, \
                                  toint(info_items[imod]@start_year), \
                                  toint(info_items[imod]@end_year), \
                                  "average", "yearly", True)

    ; Mean over the years
    mean_years = dim_avg_n_Wrap(data_yearly, 0)
    delete(data_yearly)

    ; Compute centred pattern correlation (ref_model will be first so no error)
    if datasetnames(imod).eq.ref_model then
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
  all_cor@corvar_long = variable_info[0]@long_name
  all_cor@var = "cor"
  all_cor@diag_script = DIAG_SCRIPT
  all_cor@diagnostics = variable_info[0]@diagnostic
  all_cor@ncdf = ncdf_dir
  all_cor@input = str_join(infiles, ",")
  ncdf_outfile = ncdf_write(all_cor, ncdf_dir)

  ; Write provenance
  statistics = (/"corr", "clim"/)
  domains = (/"global"/)
  plottype = "other"
  authors = (/"gier_bettina", "bock_lisa"/)
  references = (/"flato13ipcc"/)
  log_provenance(ncdf_outfile, "n/a", "n/a", statistics, domains, "other", \
                 authors, references, preproc_files)

  ; Write list of models with project (arrays get squished as attributes)
  modnames = datasetnames(mod_ind_woref)
  projnames = projectnames(mod_ind_woref)
  modproj = (/modnames, projnames/)
  asciiwrite(modproj_dir, modproj)

  leave_msg(DIAG_SCRIPT, "")

end
