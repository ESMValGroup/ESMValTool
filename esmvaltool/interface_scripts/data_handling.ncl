; #############################################################################
; GENERAL ROUTINES FOR DATA HANDLING
; #############################################################################
; Please consider using of extend
; Check the header of each routine for documentation.
;
; Contents:
;    function read_data
;    function read_fx_data
;
; #############################################################################

undef("read_data")
function read_data(index[1]:integer,
                   variable[1]:string,
                   field_type:string)
;
; Arguments
;    index: index to current input file as defined in the <variable>_info.ncl
;           file in the preproc directory of the current variable
;    variable: name of the variable to read
;    field_type: the variable's field type
;
; Return value
;    The actual data array with all metadata.
;
; Description
;    Retrieve the variable data and metadata based on the input_file_info
;    logical defined in the <variable>_info.ncl file automatically generated
;    by the preprocessor and stored in the preproc/<preproc_id>/ subdirectory.
;
; Caveats
;
; References
;
; Modification history
;   20180212_A-righ_ma: written.
;
local funcname, scriptname, infile, fin
begin

  funcname = "read_data"
  scriptname = "interface_scripts/data_handling.ncl"
  enter_msg(scriptname, funcname)

  ; Check array boundaries first
  if (index.lt.0 .or. index.gt.(dimsizes(input_file_info@filename) - 1)) then
    error_msg("f", scriptname, funcname, "dataset index out of range")
  end if

  infile = input_file_info@filename(index)

  ; Check existance and supported format
  if (.not.isfilepresent(infile)) then
    error_msg("f", scriptname, funcname, "input file " + infile + \
              " not found or not a supported format")
  end if

  ; Check field type correspondance
  if (input_file_info@field(index) .ne. field_type) then
    error_msg("f", scriptname, funcname, "requested field type " + \
              field_type + " does not match")
  end if

  ; Open file
  log_info("Reading in file = " + infile)
  fin = addfile(infile, "r")

  ; Check that variable is in the file
  if (.not. isfilevar(fin, variable)) then
    error_msg("f", scriptname, funcname, "variable " + variable + \
              " not found in file " + infile)
  end if

  ; Read data
  data = fin->$variable$

  leave_msg(scriptname, funcname)
  return(data)

end


; #############################################################################
undef("read_fx_data")
function read_fx_data(index[1]:integer,
                      variable[1]:string)
;
; Arguments
;    index: index to current input file as defined in the <variable>_info.ncl
;           file in the preproc directory of the current variable
;    variable: name of the fx variable to read
;
; Return value
;    The a list of data array(s) with all metadata or a missing value if the
;    fx variable is not found.
;
; Description
;    Retrieve the fx variable data and metadata based on the input_file_info
;    logical defined in the <variable>_info.ncl file automatically generated
;    by the preprocessor and stored in the preproc/<preproc_id>/ subdirectory.
;    Unlike read_data, this function does not exit on error if the input file
;    is not found. This is to provide more flexibility to the calling
;    diagnostic.
;
; Caveats
;    If auxiliary coordinates are associated to the fx-variables (typically for
;    areacello), these are also read and pushed to the return list.
;
; References
;
; Modification history
;   20180526_A-righ_ma: written.
;
local funcname, scriptname, infile, fin
begin

  funcname = "read_fx_data"
  scriptname = "interface_scripts/data_handling.ncl"
  enter_msg(scriptname, funcname)

  ; Check array boundaries first
  if (index.lt.0 .or. index.gt.(dimsizes(input_file_info@filename) - 1)) then
    error_msg("f", scriptname, funcname, "dataset index out of range")
  end if

  if (.not.isatt(input_file_info, variable)) then
    error_msg("f", scriptname, funcname, "input path for fx variable " + \
              variable + " not found in the interface file, please add " + \
              "the required fx variable to the recipe")
  end if

  infile = input_file_info@$variable$(index)

  ; Return missing if variable not available
  if (infile.eq."None") then
    error_msg("w", scriptname, funcname, "fx variable " + variable + \
              " not found for dataset " + dataset_info@dataset(index))
    out = new(1, float)
    out = default_fillvalue(typeof(out))
    return(out)
  end if

  ; Open file
  log_info("Reading in file = " + infile)
  fin = addfile(infile, "r")

  ; Check that variable is in the file
  if (.not. isfilevar(fin, variable)) then
    error_msg("f", scriptname, funcname, "variable " + variable + \
              " not found in " + infile)
  end if

  ; Initialize list
  outlist = NewList("lifo")

  ; Check if auxiliary coordinates are defined
  atts = getfilevaratts(fin, variable)
  if (any("coordinates".eq.atts)) then
    ; Read auxiliary coordinates in reverse order (for correct order in list)
    coords = fin->$variable$@coordinates
    ncoords = str_fields_count(coords, " ")
    do ii = 1, ncoords
      curr_coord = str_squeeze(str_get_field(coords, ncoords + 1 - ii, " "))
      if (.not. isfilevar(fin, curr_coord)) then
        error_msg("w", scriptname, funcname, "auxiliary coordinate " + \
                  curr_coord + " not found in " + infile + ", skipping")
      else
        ListPush(outlist, fin->$curr_coord$)
      end if
    end do
  end if

  ; Push variable to the list
  ListPush(outlist, fin->$variable$)

  leave_msg(scriptname, funcname)
  return(outlist)

end
