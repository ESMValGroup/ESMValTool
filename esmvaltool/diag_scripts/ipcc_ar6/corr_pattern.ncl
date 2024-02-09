; #############################################################################
; DIAGNOSTIC SCRIPT for pattern correlation figure
; Author: Bettina Gier (Uni Bremen & DLR, Germany)
; IPCC-AR6 project
; #############################################################################
;
; Description
;    Calculates centred pattern correlations for annual mean climatologies
;    and plots them. Like IPCC AR5 Ch. 9 Fig 9.6
;
; Required diag_script_info attributes (diagnostics specific)
;     none
;
; Optional diag_script_info attributes (diagnostic specific)
;     none
;
; Required variable_info attributes (variable specific)
;     reference_dataset: name of reference observation
;
; Optional variable_info attributes (variable specific)
;     alternative_dataset: name of alternative observations
;
; Caveats
;
;
; Modification history
;    20210226-A_bock_lisa: Modified diagnostic based on
;                          ipcc_ar5/ch09_fig09_6.ncl
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
  all_cor@_FillValue = 1.e20
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

  ; *************************************
  ; calculate means over ensemble members
  ; *************************************

  datasetnames := datasetnames(mod_ind_woref)
  projectnames := projectnames(mod_ind_woref)

  ndim := dimsizes(datasetnames)
  ensemble_assign = new(ndim, integer)

  unique_models = get_unique_values(datasetnames)
  do gg = 0, dimsizes(unique_models) - 1
    ensemble_assign = where(datasetnames.eq.unique_models(gg), gg + 1, \
                            ensemble_assign)
  end do
  ensemble_assign@model = datasetnames

  all_cor_em = all_cor(:max(ensemble_assign)-1)
  all_cor_em = all_cor_em@_FillValue
  datasetnames_em = datasetnames(:max(ensemble_assign)-1)
  projectnames_em = projectnames(:max(ensemble_assign)-1)

  do gg = 1, max(ensemble_assign)     ; calculate ensemble means
    wind := ind(ensemble_assign.eq.gg)
    if (dimsizes(wind).eq.1) then
      all_cor_em(gg-1) = (/ all_cor(wind) /)
      all_cor_em&models(gg-1) =  datasetnames(wind)
      datasetnames_em(gg-1) = datasetnames(wind)
      projectnames_em(gg-1) = projectnames(wind)
    else
      all_cor_em(gg-1) = (/ dim_avg_n(all_cor(wind), 0) /)
      all_cor_em&models(gg-1) =  datasetnames(wind(0))
      datasetnames_em(gg-1) = datasetnames(wind(0))
      projectnames_em(gg-1) = projectnames(wind(0))
    end if
  end do

  delete(all_cor)
  all_cor = all_cor_em
  delete(all_cor_em)

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
  modnames = datasetnames_em
  projnames = projectnames_em
  modproj = (/modnames, projnames/)
  asciiwrite(modproj_dir, modproj)

  leave_msg(DIAG_SCRIPT, "")

end
