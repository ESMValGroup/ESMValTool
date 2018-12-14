; #############################################################################
; GENERAL ROUTINES FOR MODEL ENSEMBLE OPERATIONS
; #############################################################################
; Please consider using or extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
; Contents:
;    function get_start_year
;    function get_end_year
;    function find_model_index
;
; #############################################################################

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/logging.ncl"

; #############################################################################
undef("get_start_year")
function get_start_year(dataset_index[1]:numeric)
;
; Arguments:
;    dataset_index: index of array "datasets"; this function is supposed to be
;                 called inside a dataset loop
;
; Return value
;    An integer containing the first year
;
; Description
;    Chooses first year of dataset input with respect to
;    diag_script_info@range_option
;
; Caveats
;
; References
;
; Modification history
;    20140128-A_senf_da: written.
;
local funcname, scriptname, dataset_index
begin

  funcname = "get_start_year"
  scriptname = "diag_scripts/shared/ensemble.ncl"
  enter_msg(scriptname, funcname)

  if (.not.isatt(diag_script_info, "range_option")) then
    error_msg("f", "ensemble.ncl", funcname, "diag_script_info must " + \
              "contain the range_option attribute")
  end if
  if (diag_script_info@range_option.eq.0) then
    start_year = toint(input_file_info@start_year(dataset_index))
  else if (diag_script_info@range_option.eq.1) then
    start_year = max(toint(input_file_info@start_year))
  else
    error_msg("f", "ensemble.ncl", funcname, "range_option not " + \
              "properly set, must be either 0 or 1")
  end if
  end if

  leave_msg(scriptname, funcname)
  return(start_year)

end

; #############################################################################
undef("get_end_year")
function get_end_year(dataset_index:numeric)
;
; Arguments
;    dataset_index: index of array "datasets"; this function is supposed to be
;                 called inside a dataset loop
;
; Return value
;    An integer containing the last year
;
; Description
;    Chooses last year of dataset input with respect to
;    diag_script_info@range_option
;
; Caveats
;
; References
;
; Modification history
;    20140128-A_senf_da: written.
;
local funcname, scriptname, dataset_index
begin

  funcname = "get_end_year"
  scriptname = "diag_scripts/shared/ensemble.ncl"
  enter_msg(scriptname, funcname)

  if (.not.isatt(diag_script_info, "range_option")) then
    error_msg("f", "ensemble.ncl", funcname, "diag_script_info must " + \
              "contain the range_option attribute")
  end if
  if (diag_script_info@range_option.eq.0) then
    end_year = toint(input_file_info@end_year(dataset_index))
  else if (diag_script_info@range_option.eq.1) then
    end_year = min(toint(input_file_info@end_year))
  else
    error_msg("f", "ensemble.ncl", funcname, "range_option not " + \
              "properly set, must be either 0 or 1")
  end if
  end if

  leave_msg(scriptname, funcname)
  return(end_year)

end

; #############################################################################
undef("find_model_index")
function find_model_index(model_name :string,
                          models_tot :string,
                          exp_name : string,
                          experiment :string,
                          ens_name : string,
                          ensemble :string,
                          opt[1])
;
; Arguments
;    model_name: model names to match up.
;    models_tot: all model names to look in, dataset_info@dataset.
;    exp_name : experiment to match up
;    experiment: experiments to look for matching models in,
;                dataset_info@exp.
;    ens_name : ensemble to match up,
;    ensemble : ensembles to look for match in, dataset_info@ensemble
;    opt: not used yet, set to 0
;
; Return value
;    A list of indices which match the experiment, ensemble and model name
;
; Description
;    Finds all indices which match a particular experiment name, model
;    name and ensemble. For instance to calculate changes in different RCPs
;    we only consider historical runs which are also available in the
;    particular RCP, so we need to find all historical runs for the models
;    in the RCP in the namelist.
;
; Caveats
;    Fails if we look for a model name that does not exist in the namelist.
;    Is there a more elegant and simpler way of doing this???
;    e.g. without loop?
;
; References
;
; Modification history:
;    20161031-A_lore_ru: added ensemble as matching criteria
;    20160908-A_lore_ru: written.
;
local funcname, scriptname, dim_rcp, hist, tmp_hist0, tmp_ind, \
    idx_hist
begin

    funcname = "find_model_index"
    scriptname = "diag_scripts/shared/ensemble.ncl"
    enter_msg(scriptname, funcname)

    ; how many model names do we need to match up?
    dim_rcp = dimsizes(model_name)

    ; loop over all model names in RCP and find corresponding historical run
    do hist = 0, dim_rcp - 1
        if (hist .eq. 0) then
            tmp_hist0 = 0
            idx_hist = ind(experiment .eq. exp_name .and. \
                           models_tot .eq. model_name(hist) .and. \
                           ensemble .eq. ens_name(hist))
        else
            delete(idx_hist)
            tmp_ind = ind(experiment .eq. exp_name .and. \
                          models_tot .eq. model_name(hist) .and. \
                          ensemble .eq. ens_name(hist))
            idx_hist = array_append_record(tmp_hist0, tmp_ind, 0)
            delete(tmp_ind)
        end if
        delete(tmp_hist0)
        tmp_hist0 = idx_hist
    end do
    if (any(ismissing(idx_hist))) then
        error_msg("f", "ensemble.ncl", funcname, "missing index, model " + \
                  model_name(ind(ismissing(idx_hist))) + "ensemble " + \
                  ens_name(ind(ismissing(idx_hist))) + \
                  " not found in historical runs.")
    end if
    return(idx_hist)
end
