; #############################################################################
; GENERAL SET OPERATORS
; #############################################################################
; Please consider using of extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
; Contents:
;    function UNIQ
;    function union
;    function set_inclusive_OR
;    function intersection
;    function is_array_subset
;    function relative_complement
;    function set_symmetric_difference
;
; #############################################################################

load "$diag_scripts/../interface_scripts/logging.ncl"

; #############################################################################
undef("UNIQ")
function UNIQ(a)
;
; Arguments
;    a: array to analyse, could be any type.
;
; Return value
;    Index array of first occurence of each new element.
;
; Description
;    Useful to get an array of unique elements.
;
; Caveats
;
; Reference
;
; Modification history
;    20130419-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, a, i, u, doubles, idx
begin

  funcname = "UNIQ"
  scriptname = "diag_scripts/shared/set_operators.ncl"
  enter_msg(scriptname, funcname)

  u = new(dimsizes(a), "integer", 0)
  do i = 0, dimsizes(a) - 1
    doubles = ind(a .eq. a(i))
    u(doubles(0)) = 1
    delete(doubles)
  end do
  idx = ind(u .eq. 1)

  leave_msg(scriptname, funcname)
  return(idx)

end

; #############################################################################
undef("union")
function union(set1[*],
               set2[*])
;
; Arguments
;    set1: a one-dimensional array.
;    set2: a one-dimensional array.
;
; Return value
;    An array.
;
; Description
;    Returns  the union of two sets
;
; Caveats
;
; References
;
; Modification history
;
local funcname, scriptname, merged_set, union_set, unique_entries
begin

  funcname = "union"
  scriptname = "diag_scripts/shared/set_operators.ncl"
  enter_msg(scriptname, funcname)

  merged_set = array_append_record(set1, set2, 0)
  unique_entries = UNIQ(merged_set)
  union_set = merged_set(unique_entries)

  leave_msg(scriptname, funcname)
  return(union_set)

end

; #############################################################################
undef("set_inclusive_OR")
function set_inclusive_OR(set1[*]:float,
                          set2[*]:float)
;
; Arguments
;    set1: a one-dimensional array.
;    set2: a one-dimensional array.
;
; Return value
;    An array.
;
; Description:
;   Returns sorted inclusive 'or' of two sets.
;
; Caveats
;
; References
;
; Modification history
;
local funcname, scriptname, set_tmp, dummy_year, dummy_year_idx, \
  keep_set_idx, keep_set
begin

  funcname = "set_inclusive_OR"
  scriptname = "diag_scripts/shared/set_operators.ncl"
  enter_msg(scriptname, funcname)

  dummy_year  = -999.

  set_tmp = union(set1, set2)

  ; Make sure the dummy year -999. is removed (if present)
  keep_set_idx = ind(abs(set_tmp + abs(dummy_year)) .gt. 1)
  keep_set = set_tmp(keep_set_idx)
  qsort(keep_set)

  leave_msg(scriptname, funcname)
  return(keep_set)

end

; #############################################################################
undef("intersection")
function intersection(array1[*],
                      array2[*])
;
; Arguments
;    array1: a one-dimensional array.
;    array2: a one-dimensional array.
;
; Return value
;    array: Intersection array.
;
; Description
;    Returns the intersection of array1 and array2 or 'False' if no
;    intersection is found.
;
; Caveats
;
; References
;
; Modification history
;
local funcname, scriptname, idx, intersection_array, tmp_array
begin

  funcname = "intersection"
  scriptname = "diag_scripts/shared/set_operators.ncl"
  enter_msg(scriptname, funcname)

  do idx = 0, dimsizes(array1) - 1
    curr_check = ind(array1(idx) .eq. array2)
    if (.not. ismissing(curr_check)) then
      if (isdefined("intersection_array")) then
        tmp_array2 = array_append_record(intersection_array, array1(idx), 0)
        delete(intersection_array)
        intersection_array = tmp_array2
        delete(tmp_array2)
      else
        intersection_array = (/array1(idx)/)
      end if
    end if
  end do

  if (.not. isdefined("intersection_array")) then
    intersection_array = False
  end if

  leave_msg(scriptname, funcname)
  return(intersection_array)

end

; #############################################################################
undef("is_array_subset")
function is_array_subset(subset_array[*],
                         full_array[*])
;
; Arguments
;    subset_array:  an array of dimension N.
;    full_array: an array of dimension >= N.
;
; Return value
;    A logical: True if "subset_array" a true subset of "full_array", False
;    otherwise.
;
; Description
;    Checks if an array is a subset of another array
;
; Caveats
;
; Reference
;
; Modification history
;
local funcname, scriptname, subset_size, idx, full_idx, is_subset_true
begin

  funcname = "is_array_subset"
  scriptname = "diag_scripts/shared/set_operators.ncl"
  enter_msg(scriptname, funcname)

  subset_size = dimsizes(subset_array)
  is_subset_true = True
  do idx = 0, subset_size - 1
    full_idx = ind(subset_array(idx) .eq. full_array)
    if (ismissing(full_idx)) then
      is_subset_true = False
    end if
  end do

  leave_msg(scriptname, funcname)
  return(is_subset_true)

end

; #############################################################################
undef("relative_complement")
function relative_complement(array1[*],
                             array2[*])
;
; Arguments
;    array1: an array
;    array2: another array
;
; Return value
;    logical: True if there is a complemnet, in this case the complement
;                  itself is attached as '@array'
;             False  if there isn't a complement
;
; Description
;    Substracts all elements in array1 from array2 (array2 - array1). See
;    https://en.wikipedia.org/wiki/Complement_(set_theory)#Relative_complement
;    for further details and expected behaviour
;
; Caveats
;
; Reference
;
; Modification history
;
local funcname, scriptname, array2_size, i, rel_complement, idx, remove_idx
begin

  funcname = "relative_complement"
  scriptname = "diag_scripts/shared/set_operators.ncl"
  enter_msg(scriptname, funcname)

  array2_size = dimsizes(array2)
  rel_complement = False
  do idx = 0, array2_size - 1
    remove_idx = ind(array2(idx) .eq. array1)
    if (ismissing(remove_idx)) then
      if (rel_complement .eq. False) then
        rel_complement@array = (/array2(idx)/)
        rel_complement = True
      else
        tmp = array_append_record(rel_complement@array, array2(idx), 0)
        delete_VarAtts(rel_complement, "array")
        rel_complement@array = tmp
        delete(tmp)
      end if
    end if
  end do
  if (isatt(rel_complement, "array")) then
    return_array = rel_complement@array
  else
    return_array = rel_complement
  end if

  leave_msg(scriptname, funcname)
  return(return_array)

end

; #############################################################################
undef("set_symmetric_difference")
function set_symmetric_difference(array1[*],
                                  array2[*])
;
; Arguments
;    array1: an array
;    array2: another array
;
; Return value
;    array: The symmetric difference of array1 and array2
;
; Description
;    Returns all elements only in array1 or only in array2
;    See "https://en.wikipedia.org/wiki/Symmetric_difference"
;    for furhter details and expected behaviour
;
; Caveats
;
; Reference
;
; Modification history
;
local funcname, scriptname, rel_compl1, rel_compl2
begin

  funcname = "set_symmetric_difference"
  scriptname = "diag_scripts/shared/set_operators.ncl"
  enter_msg(scriptname, funcname)

  rel_compl1 = relative_complement(array1, array2)
  rel_compl2 = relative_complement(array2, array1)

  leave_msg(scriptname, funcname)
  return(union(rel_compl1, rel_compl2))

end
