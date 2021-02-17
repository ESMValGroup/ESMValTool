; #############################################################################
; GENERAL-PURPOSE AUXILIARY NCL ROUTINES
; #############################################################################
; Please consider using of extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
; Contents:
;    function bname
;    function basename
;    function att2var
;    function att2var_default
;    function get_ncdf_name
;    function get_ncdf_dir
;    function ncdf_read
;    function ncdf_define
;    function ncdf_write
;    function ncdf_att
;    function copy_CoordNames_n
;    function extend_var_at
;    function remove_index
;    procedure set_default_att
;    function empty_str
;
; #############################################################################

undef("bname")
function bname(a[*]:string)
;
; Arguments
;    a: (array of) file name strings to analyse.
;
; Return value
;    An array of strings, without paths and suffixes.
;
; Description
;    Strips paths and suffixes from file names.
;
; Caveats
;
; References
;
; Modification history
;    20131106-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, a, i, b, rm_prefix, rm_suffix
begin

  funcname = "bname"
  scriptname = "interface_scripts/auxiliary.ncl"
  enter_msg(scriptname, funcname)

  b = a
  do i = 0, dimsizes(a) - 1
    rm_prefix = systemfunc("basename " + a(i))
    rm_suffix = str_split(rm_prefix, ".")
    b(i) = rm_suffix(0)
    delete(rm_suffix)
  end do

  leave_msg(scriptname, funcname)
  return(b)

end

; #############################################################################
undef("basename")
function basename(name[1]:string)
;
; Arguments
;    name: input string.
;
; Return value
;    A string.
;
; Description
;    Removes the suffix from the string 'name' and returns the string.
;
; Caveats
;
; References
;
; Modification history
;    20130419-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, suffix
begin

  funcname = "basename"
  scriptname = "interface_scripts/auxiliary.ncl"
  enter_msg(scriptname, funcname)

  suffix = get_file_suffix(name, 0)
  bsn = str_sub_str(name, suffix, "")

  leave_msg(scriptname, funcname)
  return(bsn)

end

; #############################################################################
undef("att2var")
function att2var(data,
                 meta[1]:string)
;
; Arguments
;    data: variable of any type with plot specific metadata.
;    meta: string specifying the requested attribute, "var", "res_" and
;          "diag_file" get special treatment.
;
; Return value
;    Requested attribute of data is returned as variable:
;      var (string): name string, with "var_" attached as attributes
;                    (without prefix).
;      diag_file (string): name bases of registred scripts, joined into one
;                          string.
;      res (logical): attributes needed for plotting.
;      other attributes (any type): as requested by meta.
;
; Description
;    Puts attributes of data into individual variables.
;
; Caveats
;
; References
;
; Modification history
;    20131107-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, data, meta, result, temp, atts, att, iatt
begin

  funcname = "att2var"
  scriptname = "interface_scripts/auxiliary.ncl"
  enter_msg(scriptname, funcname)

  ; Attribute "var"
  if (meta .eq. "var")
    if (isatt(data, meta)) then
      result = data@$meta$
      atts = getvaratts(data)
      do iatt = 0, dimsizes(atts) - 1
        att = atts(iatt)
        if (isStrSubset(att, "var_")) then
          ; Remove prefix "var_"
          temp = str_get_cols(att, 4, strlen(att))
          ; Assign contents to variable "var"
          result@$temp$ = data@$att$
        end if
      end do
    else
      error_msg("f", scriptname, funcname, get_script_name() + \
                " requires variable attribute " + meta)
    end if

  ; Attribute "diag_script"
  elseif (meta .eq. "diag_script")
    if (isatt(data, meta)) then
      base = bname(data@$meta$)
      ; Join names of all contributing scripts
      result = str_join(base, "__")
    else
      error_msg("f", scriptname, funcname, get_script_name() + \
                " requires variable attribute " + meta)
    end if

  ; Attribute "res_"
  ; DEBUG info: different attribute types are not anticipated
  ; -> check this if it crashes
  elseif (meta .eq. "res_")
    result = True
    atts = getvaratts(data)
    do iatt = 0, dimsizes(atts) - 1
      att = atts(iatt)
      if (isStrSubset(att, "res_")) then
        ; Remove prefix "res_"
        temp = str_get_cols(att, 4, strlen(att))
        ; Assign contents to variable "res"
        result@$temp$ = data@$att$
      end if
    end do

  ; Optional/generic attributes
  else
    if (isatt(data, meta)) then
      result = data@$meta$
    else
      error_msg("w", scriptname, funcname, "attribute " + meta + \
                " not found")
      result = default_fillvalue("string")
    end if
  end if

  leave_msg(scriptname, funcname)
  return(result)

end

; #############################################################################
undef("att2var_default")
function att2var_default(data,
                         meta[1]:string,
                         default,
                         option)
;
; Arguments
;    data: variable of any type with plot specific metadata.
;    meta: string specifying the requested attribute
;          "var", "res_" and "diag_file" get special treatment.
;    default: variable of any type.
;    option: not used yet.
;
; Return value
;    Requested attribute of data is returned as variable.
;    See function  att2var for details on 'magic' strings.
;    Those get special treatment and may trigger a crash, if missing.
;    Default is returned if (non magic) data@meta is missing.
;
; Description
;    Puts attributes of data into individual variables,
;    and returns default if the requested attribute is missing.
;
; Caveats
;
; References
;
; Modification history
;    20140403-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, data, meta, result, default, option
begin

  funcname = "att2var_default"
  scriptname = "interface_scripts/auxiliary.ncl"
  enter_msg(scriptname, funcname)

  result = att2var(data, meta)
  if(dimsizes(result).eq.1) then
    if(result.eq.default_fillvalue("string")) then
      log_debug("Using default: " + meta + " = " + tostring(default))
      result = default
      copy_VarMeta(default, result)
    end if
  end if

  leave_msg(scriptname, funcname)
  return(result)

end

; #############################################################################
undef("get_ncdf_name")
function get_ncdf_name(var[1]:string)
;
; Arguments
;    var: variable name.
;
; Return value
;    A string with the full outfile name (incl. path)
;    work_dir + var + ".nc".
;
; Description
;    Fetches work_dir via other scripts.
;
; References
;
; Modification history
;    20131108-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, work_dir, diag_script, output_dir, ncdf_dir
begin

  funcname = "get_ncdf_name"
  scriptname = "interface_scripts/auxiliary.ncl"
  enter_msg(scriptname, funcname)

  ; Collect parts of outfile name
  work_dir = config_user_info@work_dir + "/"

  ; Create work dir
  system("mkdir -p " + work_dir)

  ; Assemble final name
  outfile = work_dir + var + ".nc"

  leave_msg(scriptname, funcname)
  return(outfile)

end

; #############################################################################
undef("get_ncdf_dir")
function get_ncdf_dir(work_dir[1]:string,
                      subfolder[1]:string)
;
; Arguments
;    work_dir: root path for intermediate files: expected to end with a slash.
;    subfolder: subfolder for diag_script.
;
; Return value
;    A string with the full path for netCDF output.
;
; Description
;    Concatenates/generates work_dir + subfolder.
;
; Caveats
;    diag_script may need to be extended by other contributing scripts.
;
; References
;
; Modification history
;    20131108-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, work_dir, subfolder, ncdf_dir
begin

  funcname = "get_ncdf_dir"
  scriptname = "interface_scripts/auxiliary.ncl"
  enter_msg(scriptname, funcname)

  ncdf_dir = work_dir + subfolder + "/"
  if (.not. fileexists(ncdf_dir)) then
    system("mkdir -p " + ncdf_dir)
  end if

  leave_msg(scriptname, funcname)
  return(ncdf_dir)

end

; #############################################################################
undef("ncdf_read")
function ncdf_read(infile, varname)
;
; Arguments
;    infile: full path & name of the netCDF file to be read.
;
; Return value
;    Data with metadata:
;        if varname="dummy" AND only one var in file: return that var
;        if varname="default" AND only one var in file: return that var
;        if varname="dummy" AND more than one var in file: exit
;        if varname="default" AND more than one var in file: exit
;        if varname=var AND varname exists in file: return var varname
;        if varname=var AND varname does not exist: exit
;
; Description
;    http://www.ncl.ucar.edu/Applications/r-support.shtml
;
; Caveats
;    A variable is considered a coordinate variable, if it has the same name
;    as first dimension (not sure if this is a satisfactory criterion).
;    The netCDF format does not support data type "logical" and converts True
;    (logical) to 1 (integer), False (logical) to 0 (integer), see
;    http://www.ncl.ucar.edu/Support/talk_archives/2009/0367.html
;    A workaround could be to convert logical data to string before writing to
;    netCDF, and adjust reading routine accordingly.
;    Currently ncdf_read interprets variable attributes of type "integer" and
;    dimsizes = 1 & value 0 or 1 as logicals.
;    If a named coordinate was converted to a character array in ncdf_write it
;    is converted back to string and appended to the main variable. It only
;    works if the variable has an attribute named_coordinate equal 1.
;
; References
;
; Modification history
;    20131111-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, data, infile, ncdf, varname, i, ivar, nvar, \
  vdims, var, atts, temp, defaults
begin

  funcname = "ncdf_read"
  scriptname = "interface_scripts/auxiliary.ncl"
  enter_msg(scriptname, funcname)

  ; Open file for reading
  if (fileexists(infile)) then
    ncdf = addfile(infile, "r")
  else
    error_msg("f", scriptname, funcname, infile + " does not exist")
  end if

  ; Determine (number of) variables, exclusive coordinates
  vNames = getfilevarnames(ncdf)
  nvar = new(dimsizes(vNames), integer)
  nvar(:) = 1  ; used for flagging non-coordinate vars
  do ivar = 0, dimsizes(vNames) - 1
    vdims = getfilevardims(ncdf, vNames(ivar))
    if (vdims(0) .eq. vNames(ivar)) then
      nvar(ivar) = 0  ; coordinate variable
    end if
    delete(vdims)
  end do
  var = vNames(ind(nvar .eq. 1))  ; non-coordinate vars

  ; Cases for reading
  defaults = (/"default", "dummy", "dummy_for_var", "Default", "Dummy"/)
  if (any(varname .eq. defaults)) then  ; allow synonyms
    varname = "dummy"  ; unify
  end if
  if (varname .eq. "dummy" .and. dimsizes(var) .eq. 1) then
    data = ncdf->$var$
  elseif (varname .eq. "dummy" .and. dimsizes(var) .ne. 1) then
    error_msg("f", scriptname, funcname, "please select variable " + \
              "from " + infile)
  elseif (varname .ne. "dummy" .and. any(varname .eq. vNames)) then
    i = ind(vNames .eq. varname)
    data = ncdf->$vNames(i)$
  elseif (varname .ne. "dummy" .and. .not. any(varname .eq. vNames)) then
    error_msg("f", scriptname, funcname, "variable " + varname + \
              " not in " + infile)
  else
    error_msg("f", scriptname, funcname, "case not covered")
  end if

  ; Convert attributes to type logical, if they meet above conditions
  atts = getvaratts(data)
  do i = 0, dimsizes(atts) - 1
    if (typeof(data@$atts(i)$) .eq. "integer" .and. \
        dimsizes(data@$atts(i)$) .eq. 1) then
      if (data@$atts(i)$ .eq. 1 .or. data@$atts(i)$ .eq. 0) then
        temp = new(1, logical)
        temp = data@$atts(i)$
        delete(data@$atts(i)$)
        data@$atts(i)$ = temp
        delete(temp)
      end if
    end if
  end do

  ; Converts named coordinates back to strings and append them to the data
  coord_ind = ind(nvar.eq.0)
  do j = 0, dimsizes(vNames(coord_ind)) - 1
    new_coord = ncdf->$vNames(coord_ind(j))$
    if ((isatt(new_coord, "named_coordinate")).and. \
        (new_coord@named_coordinate.eq.1)) then
      data&$vNames(coord_ind(j))$ = chartostring(new_coord)
    end if
    delete(new_coord)
  end do

  leave_msg(scriptname, funcname)
  return(data)

end


; #############################################################################
undef("ncdf_att")
function ncdf_att(infile[1],
                  att[1]:string)
;
; Arguments
;    infile: full path & name of the netCDF file to be read OR a file
;            variable already created by addfile(...).
;    att: attribute to be read.
;
; Return value
;    Contents of attribute att (could be any type supported by netCDF).
;
; Description
;    Checks existence of att in infile, then retrieves contents.
;
; Caveats
;    str_infile does not resolve the file name for a file variable.
;
; References
;
; Modification history
;    20140131-A_gott_kl: written.
;
local funcname, scriptname, f, file_atts, i, result
begin

  funcname = "ncdf_att"
  scriptname = "interface_scripts/auxiliary.ncl"
  enter_msg(scriptname, funcname)

  if (typeof(infile).eq."string") then
    f = addfile(infile, "r")
    str_infile = infile
  else if (typeof(infile).eq."file") then
    f = infile
    str_infile = "file associated with file variable argument."
  else
    error_msg("f", scriptname, funcname, "wrong type of infile argument")
  end if
  end if
  file_atts = getvaratts(f)  ; retrieve global attributes
  if (any(file_atts.eq.att)) then
    i = ind(file_atts.eq.att)
    result = f@$file_atts(i)$
  else
    error_msg("w", scriptname, funcname, att + " not found in " + str_infile)
    result = "missing"
  end if

  leave_msg(scriptname, funcname)
  return(result)

end

; #############################################################################
undef("ncdf_define")
function ncdf_define(ncdf, data)
;
; Arguments
;    ncdf: file specified.
;    data: data with metadata (attributes).
;
; Return value
;    Modifies the file specifier "ncdf".
;
; Description
;
; Caveats
;    No unlimited dimensions supported (yet).
;
; References
;    http://www.ncl.ucar.edu/Applications/method_2.shtml
;
; Modification history
;    20131112-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, data, diag_script, var, fAtt, dimNames, dimSzs, \
  dimUnlim, atts
begin

  funcname = "ncdf_define"
  scriptname = "interface_scripts/auxiliary.ncl"
  enter_msg(scriptname, funcname)

  diag_script = att2var(data, "diag_script")
  var = att2var(data, "var")

  ; Explicitly declare file definition mode. Improve efficiency.
  setfileoption(ncdf, "DefineMode", True)

  ; Create global attributes of the file
  fAtt               = True  ; assign file attributes
  fAtt@creation_date = systemfunc("date")
  fAtt@created_by    = "ESMValTool, " + diag_script

  fileattdef(ncdf, fAtt)   ; copy file attributes

  ; Predefine the coordinate variables and their dimensionality
  ; Note: to get an UNLIMITED record dimension, we set the dimensionality
  ; to -1 (or the actual size) and set the dimension name to True.
  dimNames = getvardims(data)
  dim_dim = dimsizes(dimNames)
  dimSzs = new(dimsizes(dimNames), long)
  dimUnlim = new(dimsizes(dimNames), logical)
  dimUnlim(:) = False

  do i = 0, dim_dim - 1
    dimSzs(i) = dimsizes(data&$dimNames(i)$)
  end do

  do ii = 0, dim_dim - 1
    if (typeof(data&$dimNames(ii)$).eq."string") then
      ; Add necessary infomation to the dim* arrays
      dim_coord = dimsizes(stringtochar(data&$dimNames(ii)$))
      ; If coordinate has length 1
      if (dimsizes(dim_coord).eq.1) then
        dimc = dim_coord(0)
      else
        dimc = dim_coord(1)
      end if
      new_dimSzs = array_append_record(dimSzs, tolong(dimc), 0)
      new_dimUnlim = array_append_record(dimUnlim, False, 0)
      new_dimNames = array_append_record(dimNames, dimNames(ii) + "_char", 0)

      ; Delete and reattach to previous variable
      delete([/dimSzs, dimUnlim, dimNames/])
      dimSzs = new_dimSzs
      dimUnlim = new_dimUnlim
      dimNames = new_dimNames
      delete([/new_dimSzs, new_dimUnlim, new_dimNames/])
      delete(dim_coord)
      jj = 1  ; counter initialization for later use
    end if
  end do

  filedimdef(ncdf, dimNames, dimSzs, dimUnlim)
  ; Predefine the the dimensionality of the variables to be written out
  do i = 0, dim_dim - 1
    if (typeof(data&$dimNames(i)$).eq."string") then
      ; Attach second dimension (character array) to variable
      filevardef(ncdf, dimNames(i), "character", \
                 (/getvardims(data&$dimNames(i)$), \
                 dimNames(dim_dim + jj - 1)/))
      data&$dimNames(i)$@named_coordinate = True  ; set flag attribute
      jj = jj + 1
    else
      filevardef(ncdf, dimNames(i), typeof(data&$dimNames(i)$), \
                 getvardims(data&$dimNames(i)$))
    end if
  end do
  filevardef(ncdf, var, typeof(data), getvardims(data))

  ; Copy attributes associated with each variable to the file
  do i = 0, dim_dim - 1
    ; No attribute assignment to variable if no attributes are available
    if (.not.(all(ismissing(getvaratts(data&$dimNames(i)$))))) then
        filevarattdef(ncdf, dimNames(i), data&$dimNames(i)$)
    end if
  end do
  filevarattdef(ncdf, var, data)

  ; Explicitly exit file definition mode. **NOT REQUIRED**
  setfileoption(ncdf, "DefineMode", False)

  leave_msg(scriptname, funcname)
  return(ncdf)

end

; #############################################################################
undef("ncdf_write")
function ncdf_write(data,
                    outfile_in)
;
; Arguments
;    data: data with metadata (attributes) to be written to a netCDF file.
;    outfile_in: full path & name of the netCDF file to be written. It may
;                may contain an attribute "existing" with the values
;                "overwrite" (default) / "append" / "rename"
;                default is used for "default", "dummy" or non-string values.
; Return value
;    A string with the outfile.
;
; Description
;    Creates or modifies outfile, adding a single variable. All attributes of
;    "data" are added to the variable in the netCDF file. It may rename an
;    existing file to keep it.
;    Two options for netCDF output, one disabled (see "Option"), change if
;    needed.
;    If data has coordinates which are named with strings, these names are
;    converted to characeter arrays and the coordinates are appended as two
;    dimensional variables.
;
; Caveats
;    Appending not tested.
;    netCDF does not support data type "logical" and converts True (logical)
;    to 1 (integer), False (logical) to 0 (integer).
;    A workaround could be to convert logical data to string before writing to
;    netCDF, and adjust reading routine accordingly.
;    Currently ncdf_read interprets variable attributes of type "integer" and
;    dimsizes = 1 & value 0 or 1 as logicals.
;
; References
;    http://www.ncl.ucar.edu/Document/Manuals/Ref_Manual/
;      NclFormatSupport.shtml#NetCDF
;    http://www.ncl.ucar.edu/Applications/o-netcdf.shtml
;    http://www.ncl.ucar.edu/Support/talk_archives/2009/0367.html
;
; Modification history
;    20131107-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, data, outfile, outfile_in, ncdf, varname
begin

  funcname = "ncdf_write"
  scriptname = "interface_scripts/auxiliary.ncl"
  enter_msg(scriptname, funcname)

  ; Retrieve special attributes
  varname = att2var(data, "var")

  ; Determine outfile
  defaulting = (/"default", "dummy"/)
  if (typeof(outfile_in) .ne. "string" .or. \
      .not.ismissing(ind(defaulting .eq. outfile_in))) then
    outfile = get_ncdf_name(varname)
  else
    outfile = outfile_in
  end if

  ; Determine mode for dealing with existing files
  if (isatt(outfile, "existing")) then
    existing = outfile@existing
  else
    existing = "overwrite"   ; Default
  end if

  ; Open outfile (see NCL documentation for "addfile")
  if (fileexists(outfile)) then  ; Outfile already exists
    if (existing .eq. "rename") then
      sysdate = systemfunc("date +%Y%m%d%H%M%S")
      system("mv " + outfile + " " + outfile + "." + sysdate)
      ncdf = addfile(outfile, "c")
      log_debug(" Keeping existing " + outfile + " with suffix ." + sysdate)
    elseif (existing .eq. "append") then
      ncdf = addfile(outfile, "w")
      dimNames = getvardims(data)
      dim_dim = dimsizes(dimNames)
      ; define dimensions in case not defined
      do idim = 0, dim_dim-1
        if all(dimNames(idim) .ne. getvardims(ncdf)) then
          if .not. isdefined("dims_new") then
            dims_new = dimNames(idim)
          else
            dims_new := array_append_index(dims_new, dimNames(idim), 0)
          end if
        end if
      end do
      if isdefined("dims_new") then
        dimSzs = new(dimsizes(dims_new), long)
        dimUnlim = new(dimsizes(dims_new), logical)
        dimUnlim(:) = False
        do i = 0, dim_dim - 1
          dimSzs(i) = dimsizes(data&$dimsnew(i)$)
        end do
        filedimdef(ncdf, dimNames, dimSzs, dimUnlim)
      end if
      filevardef(ncdf, varname, typeof(data), dimNames)
      ; Copy attributes associated with each variable to the file
      do i = 0, dim_dim - 1
        ; No attribute assignment to variable if no attributes are available
        ; Also skip existing dimensions
        if ((.not.(all(ismissing(getvaratts(data&$dimNames(i)$))))) .and. \
            (all(dimNames(i) .ne. getvardims(ncdf)))) then
          print("in if I should not be in")
          filevarattdef(ncdf, dimNames(i), data&$dimNames(i)$)
        end if
      end do
      filevarattdef(ncdf, varname, data)
    else
      system("rm -f " + outfile)
      ncdf = addfile(outfile, "c")
      ncdf = ncdf_define(ncdf, data)
    end if
  else  ; Outfile does not exist yet
    ncdf = addfile(outfile, "c")
    ncdf = ncdf_define(ncdf, data)
  end if

  ; Write to netCDF
  ; Option 1: q&d
;    ncdf->$varname$ = data
  ; Option 2: write with explicit file definition
  dimNames = getvardims(data)
  do i = 0, dimsizes(dimNames) - 1
    if (typeof(data&$dimNames(i)$).eq."string") then
      ; Convert string array to two dimensional character array
      char_array = stringtochar(data&$dimNames(i)$)
      ncdf->$dimNames(i)$ = (/char_array/)
      delete(char_array)
    else
      ncdf->$dimNames(i)$ = (/data&$dimNames(i)$/)
    end if
  end do
  ncdf->$varname$ = (/data/)

  log_info("Wrote variable " + varname + " to " + outfile)

  leave_msg(scriptname, funcname)
  return(outfile)

end

; #############################################################################
undef("ncdf_att")
function ncdf_att(infile[1],
                  att[1]:string)
;
; Arguments
;    infile: full path & name of the netCDF file to be read OR a file
;            variable already created by addfile(...).
;    att: attribute to be read.
;
; Return value
;    Contents of attribute att (could be any type supported by netCDF).
;
; Description
;    Checks existence of att in infile, then retrieves contents.
;
; Caveats
;    str_infile does not resolve the file name for a file variable.
;
; References
;
; Modification history
;    20140131-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, f, file_atts, i, result
begin

  funcname = "ncdf_att"
  scriptname = "interface_scripts/auxiliary.ncl"
  enter_msg(scriptname, funcname)

  if (typeof(infile).eq."string") then
    f = addfile(infile, "r")
    str_infile = infile
  elseif (typeof(infile).eq."file") then
    f = infile
    str_infile = "file associated with file variable argument."
  else
    error_msg("f", scriptname, funcname, "wrong type of infile argument")
  end if
  file_atts = getvaratts(f)  ; retrieve global attributes
  if (any(file_atts.eq.att)) then
    i = ind(file_atts.eq.att)
    result = f@$file_atts(i)$
  else
    error_msg("w", scriptname, funcname, att + " not found in " + str_infile)
    result = "missing"
  end if

  leave_msg(scriptname, funcname)
  return(result)

end

; #############################################################################
undef("copy_CoordNames_n")
function copy_CoordNames_n(var_from,
                           var_to,
                           n)
;
; Arguments
;    var_from: an array of any dimensionality. Must have named dimensions and
;              coordinate variables.
;    var_to: an array of any dimensionality.
;    n: index of coordinate which is not to be copied.
;
; Return value
;    An array with the same dimensionality of var_from.
;
; Description
;:    Copies all named dimensions and coordinate variables except the n-th one.
;
; Caveats
;
; References
;
; Modification history
;    20130419-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, var_from, var_to, n, rank
begin

  funcname = "copy_CoordNames_n"
  scriptname = "interface_scripts/auxiliary.ncl"
  enter_msg(scriptname, funcname)

  rank = dimsizes(dimsizes(var_from))

  if (rank.le.n) then
    error_msg("f", scriptname, funcname, "index is not in the " + \
              "range of the variables")
  else
    names = getvardims(var_from)
    do ii = 0, rank - 1
      if (.not.ismissing(names(ii))) then
        var_to!ii = names(ii)
        if (ii.ne.n) then
          var_to&$names(ii)$ = var_from&$names(ii)$
        end if
      end if
    end do
  end if

  leave_msg(scriptname, funcname)
  return(var_to)

end

; #############################################################################
undef("extend_var_at")
function extend_var_at(var,
                       coord:integer,
                       index:integer)
;
; Arguments
;    var : array of variables to be extended .
;    coord : index of coordinate to be extended.
;    index : index where additional entry is created.
;
; Return value
;    An array with an extra dimension.
;
; Description
;    Adds to var in coord at index a new empty entry by retaining metadata.
;
; Caveats
;
; References
;
; Modification history
;    20130419-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, var, coord, index, var_new, var_dim, \
  var_dim_plus, rank
begin

  funcname = "extend_var_at"
  scriptname = "interface_scripts/auxiliary.ncl"
  enter_msg(scriptname, funcname)

  var_dim = dimsizes(var)
  var_dim_plus = var_dim
  var_dim_plus(coord) = var_dim(coord) + 1

  ; Define new variables with one more entry for diagnostics
  rank = dimsizes(var_dim)
  var_temp = new((/var_dim_plus/), typeof(var))

  ; Name cooridnates
  var_new = copy_CoordNames_n(var, var_temp, coord)

  ; Copy old entries
  if (rank.eq.1) then
    if (index - 1.ge.0) then
      var_new(:index - 1) = var(:index - 1)
    end if
    if (index.lt.var_dim(0)) then
      var_new(index + 1:) = var(index:)
    end if
  elseif (rank.eq.2) then
    if (coord.eq.0) then
      if (index - 1.ge.0) then
        var_new(:index - 1, :) = var(:index - 1, :)
      end if
      if (index.lt.var_dim(0)) then
        var_new(index + 1:, :) = var(index:, :)
      end if
    else
      if (index - 1.ge.0) then
        var_new(:, :index - 1) = var(:, :index - 1)
      end if
      if (index.lt.var_dim(1)) then
        var_new(:, index + 1:) = var(:, index:)
      end if
    end if
  elseif (rank.eq.3) then
    if (coord.eq.0) then
      if (index - 1.ge.0) then
        var_new(:index - 1, :, :) = var(:index - 1, :, :)
      end if
      if (index.lt.var_dim(0)) then
        var_new(index + 1:, :, :) = var(index:, :, :)
      end if
    elseif (coord.eq.1) then
      if (index - 1.ge.0) then
        var_new(:, :index - 1, :) = var(:, :index - 1, :)
      end if
      if (index.lt.var_dim(1)) then
        var_new(:, index + 1:, :) = var(:, index:, :)
      end if
    else
      if (index - 1.ge.0) then
        var_new(:, :, :index - 1) = var(:, :, :index - 1)
      end if
      if (index.lt.var_dim(2)) then
        var_new(:, :, index + 1:) = var(:, :, index:)
      end if
    end if
  else
    error_msg("f", scriptname, funcname, "implemented only for " + \
              "variables with dimension less or equal 3!")
  end if

  leave_msg(scriptname, funcname)
  return(var_new)

end

; #############################################################################
undef("remove_index")
function remove_index(array,
                      idx:integer)
;
; Arguments
;    array: an array of any type
;    index: the index(es) to be removed.
;
; Description
;    Remove the selected index(es) from the array's dimension 0.
;
; Caveats
;    Only ranking upto 5 is supported.
;
; References
;
; Modification history
;    20160704-righi_mattia: written.
;
local funcname, scriptname, dsize, rank, iarr
begin

  funcname = "remove_index"
  scriptname = "interface_scripts/auxiliary.ncl"
  enter_msg(scriptname, funcname)

  dsize = dimsizes(array)
  rank = dimsizes(dsize)
  dsize := dsize(0)
  iarr = ispan(0, dsize - 1, 1)
  iarr@_FillValue = default_fillvalue(typeof(iarr))

  iarr(idx) = iarr@_FillValue
  iarr := iarr(ind(.not.ismissing(iarr)))

  if (rank.eq.1) then
    out = array(iarr)
  end if
  if (rank.eq.2) then
    out = array(iarr, :)
  end if
  if (rank.eq.3) then
    out = array(iarr, :, :)
  end if
  if (rank.eq.4) then
    out = array(iarr, :, :, :)
  end if
  if (rank.eq.5) then
    out = array(iarr, :, :, :, :)
  end if
  if (rank.gt.5) then
    error_msg("f", funcname, scriptname, "rank greater than 5 not supported")
  end if

  leave_msg(scriptname, funcname)
  return(out)

end

; #############################################################################
undef("set_default_att")
procedure set_default_att(info:logical,
                          att:string,
                          default)
;
; Arguments
;    info: a logical containing all attributes (usually diag_script_info)
;    att: the attribute to be checked for
;    value: the default value to be assigned if the attribute is not present
;
; Description
;    Check for a given attribute and assign the given default value if not
;    present
;
; References
;
; Modification history
;    20171214-righi_mattia: written.
;
local funcname, scriptname
begin

  funcname = "set_default_att"
  scriptname = "interface_scripts/auxiliary.ncl"
  enter_msg(scriptname, funcname)

  if (.not.isatt(info, att)) then
    info@$att$ = default
  end if

  leave_msg(scriptname, funcname)

end

; #############################################################################
undef("empty_str")
function empty_str(nn[1]:integer)
;
; Arguments
;     nn: length of the returned string.
;
; Return value
;     An empty string of length nn.
;
; Description
;     Returns an empty string of length nn (useful for text formatting).
;
; Caveats
;
; References
;
; Modification history
;    20140228-righi_mattia: written.
;
local funcname, scriptname, ii
begin

  funcname = "empty_str"
  scriptname = "interface_scripts/auxiliary.ncl"
  enter_msg(scriptname, funcname)

  str = " "
  do ii = 1, nn - 1
    str = str + " "
  end do

  leave_msg(scriptname, funcname)
  return(str)

end
