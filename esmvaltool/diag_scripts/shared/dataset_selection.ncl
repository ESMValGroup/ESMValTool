; #############################################################################
; GENERAL ROUTINES FOR SELECTING DATASETS
; #############################################################################
; Please consider using of extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
; Contents:
;    function get_obs
;    function get_mod
; #############################################################################


; #############################################################################
undef("get_obs")
function get_obs(names[*]:string, projects[*]:string, exclude[*]:string)
;
; Arguments
;    names: names of datasets
;    projects: project class of datasets
;    exclude: optional list of observational datasets to explicitely exclude
;
; Return value
;    An integer array containing all indices of observational datasets
;    (project classes OBS, OBS6, obs4mips).
;    ERA5 data are handled separately because they can be of project
;    OBS6 or nativ6. ERA5 (if present) are also returned as observational
;    dataset.
;
; Description
;    Checks the project class of all datasets and returns all indices
;    of datasets that belong to one of the observational classes.
;
; Caveats
;    ERA5 is the only "native6" observational dataset that can be found
;    at the moment. An extension to include other native6 observations
;    is needed once more of such datasets are available.
;
; References
;
; Modification history
;    20211018-lauer_axel: added option to exclude observational datasets
;    20210415-lauer_axel: written.
;
local funcname, scriptname, i, idx, idxobs
begin

  funcname = "get_obs"
  scriptname = "diag_scripts/shared/dataset_selection.ncl"
  enter_msg(scriptname, funcname)

  ; find indices of all OBS and obs4mips datasets
  ; (treat ERA5 as special case)

  idx = new(dimsizes(projects), integer)

  do i = 0, dimsizes(projects) - 1
    if (isStrSubset(str_lower(projects(i)), "obs"))
      idx(i) = 1
    else
      idx(i) = 0
    end if
    ; ERA5 might be native6 and thus not yet included in the list
    ; of observations
    if (isStrSubset(str_lower(names(i)), "era5") .and. \
        isStrSubset(str_lower(projects(i)), "native")) then
      idx(i) = 1
    end if

    ; check for observational datasets to exclude (optional)

    do j = 0, dimsizes(exclude) - 1
      if (names(i) .eq. exclude(j)) then
        idx(i) = 2
      end if
    end do
  end do

  idxobs = ind(idx .eq. 1)

  if (all(ismissing(idxobs))) then
    idxobs = -1
  else
    log_info("The following observational datasets are used for calculating" \
             + " multiobs products: " + str_join(names(idxobs), ", "))
    idxexcl = ind(idx .eq. 2)
    if (.not.all(ismissing(idxexcl))) then
      log_info("Observational datasets excluded by user: " \
               + str_join(names(idxexcl), ", "))
    end if
  end if

  leave_msg(scriptname, funcname)
  return(idxobs)

end

; #############################################################################
undef("get_mod")
function get_mod(names[*]:string, projects[*]:string)
;
; Arguments
;    names: names of datasets
;    projects: project class of datasets
;
; Return value
;    An integer array containing all indices of model datasets
;    excluding MultiModelMean / MultiModelMedian (if present).
;
; Description
;    Checks the project class of all datasets and returns all indices
;    of datasets that do not belong to one of the observational classes and
;    that are not "MultiModelMean" or "MultiModelMedian".
;
; Caveats
;    The dataset names MultiModelMean and MultiModelMedian are hardcoded and
;    need to be adjusted if the corresponding preprocessor creating these
;    datasets should change.
;
; References
;
; Modification history
;    20210415-lauer_axel: written.
;
local funcname, scriptname, i, idx, idxmod, idxobs, mm_ind1, mm_ind2
begin

  funcname = "get_mod"
  scriptname = "diag_scripts/shared/dataset_selection.ncl"
  enter_msg(scriptname, funcname)

  idx = new(dimsizes(projects), integer)
  idx = 1

  ; flag all observations (if present)

  idxobs = get_obs(names, projects, "")

  if (idxobs(0) .ne. -1) then
    idx(idxobs) = 0
  end if

  ; flag MultiModelMean and MultiModelMedian (if preset)

  mm_ind1 = ind(names .eq. "MultiModelMean")
  mm_ind2 = ind(names .eq. "MultiModelMedian")

  if (.not. ismissing(mm_ind1)) then
    idx(mm_ind1) = 0
  end if

  if (.not. ismissing(mm_ind2)) then
    idx(mm_ind2) = 0
  end if

  idxmod = ind(idx .eq. 1)

  if (all(ismissing(idxmod))) then
    idxmod = -1
  else
    log_info("The following model datasets have been found (exluding " \
             + "MultiModelMean/MultiModelMedian): " \
             + str_join(names(idxmod), ", "))
  end if

  leave_msg(scriptname, funcname)
  return(idxmod)

end
