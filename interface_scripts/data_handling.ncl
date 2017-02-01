;;#############################################################################
;; GENERAL ROUTINES FOR DATA HANDLING
;;#############################################################################
;; Please consider using of extend
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    procedure add_data_var
;;    function get_varName_from_file
;;    function extract_data_extended
;;    function extract_data
;;    procedure modify_data
;;    function read_data
;;
;;#############################################################################

load "./interface_data/ncl.interface"
load "./interface_scripts/messaging.ncl"

; a_laue_ax+
; global variables to save the names of all files read by "read_data"
gFilesRead0 = new(10000, string)
gNFiles0 = 0
; a_laue_ax-

;;#############################################################################
undef("add_data_var")
procedure add_data_var(index[1]:integer,
                       output_data[1]:logical,
                       new_data:numeric,
                       variable:string)
;;
;; Arguments
;;    output_data: The data structure into which 'new_data' will be packed.
;;                   @<filepath>_filepath: a string with the file path.
;;                   @<filepath>_dimNames: a list with the dimension names.
;;                   @<filepath>_inMemo: = 0 if data is on disk.
;;                                       = 1 if data is in memory.
;;                   @<filepath>: if inMemo set to 1, this attribute contains
;;                                the data.
;;                   @<filepath>_dimName(j)_values: values of dimension j.
;;                   @<filepath>_dimName(j)_attNames: attribute names of
;;                                                    dimension j.
;;                   @<filepath>_dimName(j)_attributes_attNames(k): attribute k
;;                                                                  of dim j.
;;                   @<filepath>_attributes_attNames(j): global attribute
;;                                                       number j.
;;    new_data: the data to add.
;;    index: index pointing to the current model.
;;    variable: the variable to add.
;;
;; Description
;;    Adds data as attributes to the internal data structure ('output_data')
;;    used by the tool. This structure is a dictionary where the keys are
;;    attribute names to the logical variable 'output_data' (see above),
;;    pointing to their respective attribute value.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;
local funcname, scriptname, verbosity, maxDataSize, field_type, ds, rank, \
    varType, datasize, i, out_file, cn, cnt, dimNames, fout, attNames, fAtt, \
    comment, file_atts
begin

    funcname = "add_data_var"
    scriptname = "interface_scripts/data_handling.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 2)

    maxDataSize = max_data_filesize(0)
    field_type = derived_field_type

    output_data = True
    ds = dimsizes(new_data)
    rank = dimsizes(ds)

    ;; Estimate the data size of the current variable
    varType = typeof(new_data)
    if (varType.eq."double") then
        datasize = 8.
    else
        datasize = 4.
    end if
    do j = 0, rank - 1, 1
        datasize = datasize * ds(j)
    end do
    delete(ds)
    datasize = datasize / 1000000  ; in Mbyte

    ;; Set output filename
    ;; 'get_fullpath' is loaded from  'interface_data/ncl.interface'
    out_file = interface_get_fullpath(variable, field_type, index)

    ;; 'cn' is the basic building block for the keys used in 'output_data'
    ;; see 'interface_data/ncl.interface'
    cn = interface_get_dictkeys_no_var(index)

    ;; 'cnt' is used (and reused) to specific designate specific keys
    ;; in the 'output' data structure. Here it points out the file path
    cnt = cn + "_filepath"
    output_data@$cnt$ = out_file
    dimNames = getvardims(new_data)
    cnt = cn + "_dimNames"
    output_data@$cnt$ = dimNames

    if (datasize .gt. maxDataSize) then
        ;; Only on hdd, not in current memory
        cnt = cn + "_inMemo"
        output_data@$cnt$ = 0
        output_data@$cn$ = 0
    else
        ;; Read into memory
        cnt = cn + "_inMemo"
        output_data@$cnt$ = 1
        output_data@$cn$ = (/new_data/)
    end if

    ;; Output data to file
    info_output("adding file " + out_file, verbosity, 2)
    fout = addfile(out_file, "c")
    filedimdef(fout, "time", -1, True)
    fout->$variable$ = new_data

    if (dimNames(0) .ne. "ncl_scalar") then
        do j = 0, rank - 1, 1
            cnt = cn + "_" + dimNames(j) + "_values"

            ;; Copy variable dimension data
            output_data@$cnt$ = new_data&$dimNames(j)$

            ;; Copy attribute for coord variables
            cnt = cn + "_" + dimNames(j) + "_attNames"
            attNames = getvaratts(new_data&$dimNames(j)$)
            output_data@$cnt$ = attNames
            do k = 0, dimsizes(attNames) - 1, 1
                cnt = cn + "_" + dimNames(j) + "_attributes_" + attNames(k)
                output_data@$cnt$ = new_data&$dimNames(j)$@$attNames(k)$
            end do
            delete(attNames)
        end do
    end if

    ;; Add global attributes
    attNames = getvaratts(new_data)
    cnt = cn + "_attNames"
    output_data@$cnt$ = attNames
    do j = 0, dimsizes(attNames) - 1, 1
        cnt = cn + "_attributes_" + attNames(j)
        output_data@$cnt$ = new_data@$attNames(j)$
    end do

    delete(dimNames)
    delete(attNames)

    fAtt = True  ; assign file attributes
    fAtt@fields = field_type
    comment = "created by add_data.ncl"

    file_atts = getvaratts(fAtt)
    do i = 0, dimsizes(file_atts) - 1
        if (fAtt@$file_atts(i)$ .eq. "") then
            delete(fAtt@$file_atts(i)$)
        end if
    end do

; a_laue_ax+
    do i = 0, gNFiles0 - 1
        fAtt@$sprinti("infile_%04d", i)$ = gFilesRead0(i)
    end do

    ; reset
    gNFiles0 = 0

;    flist = [/gFilesRead0/]
;    print("=========================================================================")
;    print("add_var_data")
;    print("---------")
;    print_table(flist, "%s")
;    print("=========================================================================")
;
;    gFilesRead0 = ""

    fAtt@version = getenv("0_ESMValTool_version")
; a_laue_ax-

    fileattdef(fout, fAtt)  ; copy file attributes

    delete(fout)
    delete(fAtt)
    delete(file_atts)
    info_output("added file " + out_file, verbosity, 2)

    leave_msg(scriptname, funcname, 2)

end

;;#############################################################################
undef("get_varName_from_file")
function get_varName_from_file(filename[1]:string)
;;
;; Arguments
;;    filename: the input filename.
;;
;; Return value
;;    The names of the variable in the file.
;;
;; Description
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;
local funcname, scriptname, verbosity, variable, delimiter, nfields, \
    shortfile, fil, varNames, i
begin

    funcname = "get_varName_from_file"
    scriptname = "interface_scripts/data_handling.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 2)

    variable = new(1, string)
    delimiter = "/"
    nfields = str_fields_count(filename, delimiter)
    shortfile = str_get_field(filename, nfields, delimiter)
    fil = addfile(filename, "r")
    varNames = getfilevarnames(fil)
    do i = 0, dimsizes(varNames) - 1
        if isStrSubset(shortfile, "_" + varNames(i) + "_") then
            variable = varNames(i)
        end if
    end do
    delete(varNames)

    if (ismissing(variable)) then
        error_msg("f", "data_handling.ncl", funcname, "non-standard " + \
                  "filename: it does not contain varName, check the " + \
                  "dignostic code: " + filename)
    end if

    leave_msg(scriptname, funcname, 2)
    return(variable)
end

;;#############################################################################
undef("extract_data_extended")
function extract_data_extended(index[1]:integer,
                               var[1]:string,
                               data_pointer[1]:logical,
                               dim_id,
                               lw_limit[1]:numeric,
                               up_limit[1]:numeric)
;;
;; Arguments
;;    index: index to current "dictionary keys" (see below).
;;    data_pointer: a pointer to the data which will be extracted.
;;                    @<filepath>_filepath: a string with the file path.
;;                    @<filepath>_dimNames:  a list with the dimension names.
;;                    @<filepath>_inMemo: = 0 if data is on disk
;;                                        = 1 if data is in memory
;;                    @<filepath>: if inMemo set to 1, this attribute contains
;;                                 the data.
;;                    @<filepath>_dimName(j)_values: values of dimension j.
;;                    @<filepath>_dimName(j)_attNames: attribute names of
;;                                                     dimension j.
;;                    @<filepath>_dimName(j)_attributes_attNames(k): attrib. k
;;                                                                   of dim j.
;;                    @<filepath>_attributes_attNames(j): global attribute #
;;                                                        number j.
;;    dim_id: which dimension to fetch, either an integer or a dimension
;;            name (-1=all dimensions).
;;    lw_limit: lower limit of values to fetch. If an integer, lw_limit is
;;              interpreted as an lower index in the array. If any other
;;              numeric type it is interpreted as a lower cut off value. For
;;              dimension 'time' the expected syntax is YYYY.fraction_of_year.
;;              Set to zero if not used.
;;    up_limit: similar to lw_limit.
;;
;; Return value
;;    The actual data array with all metadata.
;;
;; Description
;;    Uses the variable 'dictkeys' defined in the
;;    'interface_data/ncl.interface' file to access and extract data from
;;    variable 'data'. 'dim_id' and 'lw_limit'/'up_limit' are used to fetch a
;;    subset of the requested data. See description in arguments above for
;;    further details.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;
local funcname, scriptname, verbosity, cn, cnt, dimNames, rank, inMemo, \
    dim_id_num, lw_idx, up_idx, dim_val, attNames, dim_val_idx, dim_val_time, \
    modelfile, variable, temp_data, i, j
begin

    funcname = "extract_data_extended"
    scriptname = "interface_scripts/data_handling.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 6)

    ;; 'cn' is the basic building block for the keys used in
    ;; the 'data_pointer' dictionary
    ;; 'dictkeys' is loaded from 'interface_data/ncl.interface'
    cn = interface_get_dictkeys(var, index)

    ;; 'cnt' is used (and reused) to specify a key in the dictionary.
    ;; Here it is used to fetch the dimension names.
    cnt = cn + "_dimNames"
    dimNames = data_pointer@$cnt$
    rank = dimsizes(dimNames)

    cnt = cn + "_inMemo"
    inMemo = data_pointer@$cnt$

    if (dimNames(0) .ne. "ncl_scalar") then
        ;; Find out dimension/coordinate number.
        if isnumeric(dim_id) then
            dim_id_num = dim_id
        else
            dim_id_num = ind(dimNames.eq.dim_id)
            if (all(ismissing(dim_id_num))) then
                error_msg("f", "data_handling.ncl", funcname, \
                          "dimension not found: " + dim_id)
            end if
        end if

        ;; Find out upper/lower limit index of dimension/coordinate to select
        if (dim_id_num .ne. -1) then
            if (isinteger(lw_limit) .and. isinteger(up_limit)) then
                ;; lw_limit & up_limit are index
                lw_idx = lw_limit
                up_idx = up_limit
            else
                ;; lw_limit & up_limit are value limits
                cnt = cn + "_" + dimNames(dim_id_num) + "_values"
                dim_val = data_pointer@$cnt$

                if (dimNames(dim_id_num).eq."time") then
                    cnt = cn + "_" + dimNames(dim_id_num) + "_attNames"
                    attNames = data_pointer@$cnt$
                    do i = 0, dimsizes(attNames) - 1, 1
                        cnt = cn + "_" + dimNames(dim_id_num) + \
                            "_attributes_" + attNames(i)
                        dim_val@$attNames(i)$ = data_pointer@$cnt$
                    end do
                    delete(attNames)

                    ;; Option '4' gives time in format YYYY.fraction_of_year
                    dim_val_time = cd_calendar(dim_val, 4)
                    dim_val_idx = ind((dim_val_time.ge.lw_limit)\
                                      .and.(dim_val_time.le.up_limit))
                    delete(dim_val_time)
                else
                    dim_val_idx = ind((dim_val.ge.lw_limit)\
                                      .and.(dim_val.le.up_limit))
                end if
                if (all(ismissing(dim_val_idx))) then
                    error_msg("f", "data_handling.ncl", funcname, \
                              "outside of range (" + \
                              dimNames(dim_id_num) + "): " + \
                              flt2string(lw_limit) + ", " + \
                              flt2string(up_limit))
                end if
                lw_idx = dim_val_idx(0)
                up_idx = dim_val_idx(dimsizes(dim_val_idx) - 1)
                delete(dim_val)
                delete(dim_val_idx)
            end if
        end if

        if (inMemo .eq. 0) then
            ;; Read from hdd
            cnt = cn + "_filepath"
            modelfile = data_pointer@$cnt$

            ;; Get variable name from file
            variable = get_varName_from_file(modelfile)

            fm = addfile(modelfile, "r")
            if (dim_id_num .eq. -1) then  ; fetch all data
                data_new = fm->$variable$
            else  ; fetch a subset defined by lw_idx/up_dix
                if (rank.eq.4) then
                    if (dim_id_num.eq.0) then
                        data_new = fm->$variable$(lw_idx:up_idx, :, :, :)
                    end if
                    if (dim_id_num.eq.1) then
                        data_new = fm->$variable$(:, lw_idx:up_idx, :, :)
                    end if
                    if (dim_id_num.eq.2) then
                        data_new = fm->$variable$(:, :, lw_idx:up_idx, :)
                    end if
                    if (dim_id_num.eq.3) then
                        data_new = fm->$variable$(:, :, :, lw_idx:up_idx)
                    end if
                end if
                if (rank.eq.3) then
                    if (dim_id_num.eq.0) then
                        data_new = fm->$variable$(lw_idx:up_idx, :, :)
                    end if
                    if (dim_id_num.eq.1) then
                        data_new = fm->$variable$(:, lw_idx:up_idx, :)
                    end if
                    if (dim_id_num.eq.2) then
                        data_new = fm->$variable$(:, :, lw_idx:up_idx)
                    end if
                end if
                if (rank.eq.2) then
                    if (dim_id_num.eq.0) then
                        data_new = fm->$variable$(lw_idx:up_idx, :)
                    else
                        data_new = fm->$variable$(:, lw_idx:up_idx)
                    end if
                end if
                if (rank.eq.1) then
                    data_new = fm->$variable$(lw_idx:up_idx)
                end if
            end if
            delete(fm)
        else
            ;; Read from memory
            temp_data = data_pointer@$cn$

            ;; Loop over all dimensions and create dimensions in data_new
            ;; and copy dimension values to data_new
            do i = 0, dimsizes(dimNames) - 1, 1
                cnt = cn + "_" + dimNames(i) + "_values"
                temp_data!i = dimNames(i)
                temp_data&$dimNames(i)$ = data_pointer@$cnt$

                ;; Copy out coord attributes
                cnt = cn + "_" + dimNames(i) + "_attNames"
                attNames = data_pointer@$cnt$
                do j = 0, dimsizes(attNames) - 1, 1
                    cnt = cn + "_" + dimNames(i) + "_attributes_" + attNames(j)
                    temp_data&$dimNames(i)$@$attNames(j)$ = data_pointer@$cnt$
                end do
                delete(attNames)
            end do

            ;; Get attribute names
            cnt = cn + "_attNames"
            attNames = data_pointer@$cnt$
            ;; Create attributes of data_new
            do i = 0, dimsizes(attNames) - 1, 1
                cnt = cn + "_attributes_" + attNames(i)
                temp_data@$attNames(i)$ = data_pointer@$cnt$
            end do

            ;; Fetch data
            if (dim_id_num .eq. -1) then  ; fetch all data
                data_new = temp_data
            else  ; fetch a subset defined by lw_idx/up_dix
                if (rank.eq.4) then
                    if (dim_id_num.eq.0) then
                        data_new = temp_data(lw_idx:up_idx, :, :, :)
                    end if
                    if (dim_id_num.eq.1) then
                        data_new = temp_data(:, lw_idx:up_idx, :, :)
                    end if
                    if (dim_id_num.eq.2) then
                        data_new = temp_data(:, :, lw_idx:up_idx, :)
                    end if
                    if (dim_id_num.eq.3) then
                        data_new = temp_data(:, :, :, lw_idx:up_idx)
                    end if
                end if
                if (rank.eq.3) then
                    if (dim_id_num.eq.0) then
                        data_new = temp_data(lw_idx:up_idx, :, :)
                    end if
                    if (dim_id_num.eq.1) then
                        data_new = temp_data(:, lw_idx:up_idx, :)
                    end if
                    if (dim_id_num.eq.2) then
                        data_new = temp_data(:, :, lw_idx:up_idx)
                    end if
                end if
                if (rank.eq.2) then
                    if (dim_id_num.eq.0) then
                        data_new = temp_data(lw_idx:up_idx, :)
                    else
                        data_new = temp_data(:, lw_idx:up_idx)
                    end if
                end if
                if (rank.eq.1) then
                    data_new = temp_data(lw_idx:up_idx)
                end if
            end if
            delete(temp_data)
        end if

        ;; Reduce rank/dimension of return array if necessary
        if ((dim_id_num .ne. -1).and.(lw_idx.eq.up_idx)) then
            temp_data = dim_avg_n_Wrap(data_new, dim_id_num)
            delete(data_new)
            data_new = temp_data
            delete(temp_data)
        end if
    else
        ;; Copy data values into data_new
        data_new = data_pointer@$cn$

        ;; Get attribute names
        cnt = cn + "_attNames"
        attNames = data_pointer@$cnt$
        ;; Create attributes of data_new
        do i = 0, dimsizes(attNames) - 1, 1
            cnt = cn + "_attributes_" + attNames(i)
            data_new@$attNames(i)$ = data_pointer@$cnt$
        end do
    end if

    leave_msg(scriptname, funcname, 6)
    return(data_new)

end

;;#############################################################################
undef("extract_data")
function extract_data(index[1]:integer,
                      data_pointer[1]:logical,
                      dim_id,
                      lw_limit[1]:numeric,
                      up_limit[1]:numeric)
;;
;; Arguments
;;
;; Return value
;;    The actual data array with all metadata.
;;
;; Description
;;    Wrapper for the extended version below. Calling the extended
;;    version directly is necessary when different variables are
;;    linked to different input files.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;
local funcname, scriptname, verbosity, data_new
begin

    funcname = "extract_data"
    scriptname = "interface_scripts/data_handling.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 6)

    data_new = extract_data_extended(index, "default", data_pointer, dim_id, \
                                     lw_limit, up_limit)

    leave_msg(scriptname, funcname, 6)
    return(data_new)

end

;;#############################################################################
undef("modify_data")
procedure modify_data(index[1]:integer,
                      data_pointer[1]:logical,
                      data_new:numeric)
;;
;; Arguments
;;    index: index to current "dictionary keys" (see below)
;;    data_pointer: a pointer to the data which will be modified
;;    data_new: the data to insert into the data_pointer structure
;;
;; Description
;;    Uses the variable 'dictkeys' defined in 'interface_data/ncl.interface'
;;    to access and modify data describeb by the variable 'data_pointer'. The
;;    'data_pointer' structure is effectively a dictionary with the following
;;    key-value pairs.
;;        KEY: VALUE
;;        -------------------------------------------------------------
;;        data_pointer@<filepath>_filepath: a string with the file path
;;        data_pointer@<filepath>_dimNames: a list with the dimension names
;;        data_pointer@<filepath>_inMemo: =0 if data is on disk,
;;                                        =1 if data is in memory
;;        data_pointer@<filepath>: if inMemo = 1, this attribute contains the
;;                                 data
;;        data_pointer@<filepath>_dimName(j)_values: values of dimension j
;;        data_pointer@<filepath>_dimName(j)_attNames: attribute names of
;;                                                     dimension j
;;        data_pointer@<filepath>_dimName(j)_attributes_attNames(k): attribute
;;                                                                   k of
;;                                                                   dimension
;;                                                                   j
;;        data_pointer@<filepath>_attributes_attNames(j): global attribute
;;                                                        number j
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;
local funcname, scriptname, verbosity, cnt, inMemo, modelfile, fout, vNames, \
    file_atts, fin, variable, delimiter, file_parts, idxLast, tempfile
begin

    funcname = "modify_data"
    scriptname = "interface_scripts/data_handling.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 2)

    ;; 'cn' is the basic building block for the keys used in the 'data_pointer'
    ;; dictionary
    ;; (see 'interface_data/ncl.interface')
    cn = interface_get_dictkeys_no_var(index)

    cnt = cn + "_inMemo"
    inMemo = data_pointer@$cnt$

    if (inMemo .eq. 1) then
        ;; Copy data_new values into data
        data_pointer@$cn$ = data_new
    end if

    ;; Update attribute (@units)
    cnt = cn + "_attributes_units"
    data_pointer@$cnt$ = data_new@units

    ;; Write to hdd
    cnt = cn + "_filepath"
    modelfile = data_pointer@$cnt$

    ;; Modify filename
    delimiter = "/"
    file_parts = str_split(modelfile, delimiter)
    idxLast = dimsizes(file_parts) - 1
    file_parts(idxLast) = "junk.nc"
    tempfile = str_join(file_parts, delimiter)
    if (file_parts(0).ne.".") then
        tempfile = "/" + tempfile
    end if
    info_output("tempfile = " + str_sub_str(tempfile, pwd, "."), verbosity, 2)

    ;; Get variable name
    variable = get_varName_from_file(modelfile)

    ;; Output data to file
    fout = addfile(tempfile, "c")
    filedimdef(fout, "time", -1, True)
    fout->$variable$ = data_new

    fin = addfile(modelfile, "r")
    vNames = getfilevarnames(fin)
    if (sum(where(vNames .eq. "hyam", 1, 0)) .gt. 0) then
        fout->hyam = fin->hyam
    end if
    if (sum(where(vNames .eq. "hybm", 1, 0)) .gt. 0) then
        fout->hybm = fin->hybm
    end if
    if (sum(where(vNames .eq. "hyai", 1, 0)) .gt. 0) then
        fout->hyai = fin->hyai
    end if
    if (sum(where(vNames .eq. "hybi", 1, 0)) .gt. 0) then
        fout->hybi = fin->hybi
    end if
    if (sum(where(vNames .eq. "p0", 1, 0)) .gt. 0) then
        fout->p0 = fin->p0
    end if

    ;; Copy global attributes
    file_atts = getvaratts(fin)
    if (.not. all(ismissing(file_atts))) then
        do i = 0, dimsizes(file_atts) - 1
            fout@$file_atts(i)$ = fin@$file_atts(i)$
        end do
    end if
    delete(file_atts)

    delete(fin)
    delete(vNames)
    delete(fout)

    system("mv " + tempfile + " " + modelfile)

    leave_msg(scriptname, funcname, 2)

end

;;#############################################################################
undef("read_data")
function read_data(index[1]:integer,
                   variable[1]:string,
                   field_number:string)
;;
;; Arguments
;;    index: index to current infile defined in 'interface_data/ncl.interface'
;;    variable: the variable to read
;;    field_number: the type of data requested (see the 'doc/*.pdf')
;;
;; Return value
;;    The actual data array with all metadata.
;;
;; Description
;;    Uses the 'interface_get_fullpath' function defined in
;;    'interface_data/ncl.interface' to find the current indata file. Read the
;;    indata file into the 'data_new' structure. If the indata file is
;;    considered too big, only the meta data is read into the structure
;;    'data_new'.
;;    The data_new structure is described below. It is basically a
;;    dictionary where the keys are attribute names to the logical variable
;;    'data_new', pointing to the attribute values.
;;    The keys-values are:
;;        KEY: VALUE
;;        -------------------------------------------------------------
;;        data_new@<filepath>_filepath: a string with the file path
;;        data_new@<filepath>_dimNames: a list with the dimension names
;;        data_new@<filepath>_inMemo: =0 if data is on disk
;;                                    =1 if data is in memory
;;        data_new@<filepath>: if inMemo set to 1, this attribute contains the
;;                             data
;;        data_new@<filepath>_dimName(j)_values: values of dimension j
;;        data_new@<filepath>_dimName(j)_attNames: attribute names of dimension
;;                                                 j
;;        data_new@<filepath>_dimName(j)_attributes_attNames(k): attribute k
;;                                                               of dimension j
;;        data_new@<filepath>_attributes_attNames(j): global attribute number j
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;
;;   20151117-A-laue_ax: changed "in_file" in call to "inquire_and_save_fileinfo"
;;                       to "ls_in_file"
;;
local funcname, scriptname, verbosity, attNames, cn, cnt, data, data_new, \
    datasize, dimNames, dimt, fid, field_number, in_file, j, k, maxDataSize, \
    rank, temp, temp1, var0, varCoords, varType, verbosity, curr_var_idx
begin

    funcname = "read_data"
    scriptname = "interface_scripts/data_handling.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 4)

    maxDataSize = max_data_filesize(0)

    temp = stringtochar(field_number)
    temp1 = stringtochar(derived_field_type)
    if temp(0).eq."*" then
        temp(0) = temp1(0)
    end if
    if temp(2).eq."*" then
        temp(2) = temp1(2)
    end if
    field_number = chartostring(temp)
    delete(temp)
    delete(temp1)
    data_new = True

    ;; 'interface_get_fullpath' is defined in 'interface_data/ncl.interface'
    in_file = interface_get_fullpath(variable, field_number, index)

    ;; Use 'ls' to remove globbing

    ls_in_file = systemfunc("ls " + in_file)

    info_output("Reading in file = " + ls_in_file, verbosity, 1)
    fid = addfile(ls_in_file, "r")

; a_laue_ax+
    if (gNFiles0.lt.dimsizes(gFilesRead0)) then
        gFilesRead0(gNFiles0) = ls_in_file
        gNFiles0 = gNFiles0 + 1
    end if

    ; read attributes and write info to log-file
    inquire_and_save_fileinfo(ls_in_file, models@name(index), variable)
; a_laue_ax-

    ;; 'cn' is the basic building block for the keys used in the
    ;; 'data_new' dictionary
    ;; 'dictkeys' is loaded from 'interface_data/ncl.interface'
    cn = interface_get_dictkeys(variable, index)

    ;; 'cnt' is used (and reused) to specify a key in the dictionary,
    ;; here it is the file path
    cnt = cn + "_filepath"
    data_new@$cnt$ = in_file

    ;; ***********************************************
    ;; Estimate the data size of the current variable
    ;; ***********************************************
    varType = getfilevartypes(fid, variable)
    if (varType.eq."double") then
        datasize = 8.
    else
        datasize = 4.
    end if

    dimt = getfilevardimsizes(fid, variable)
    rank = dimsizes(dimt)
    do j = 0, rank - 1, 1
        datasize = datasize * dimt(j)
    end do
    delete(dimt)

    datasize = datasize / 1000000  ; In Mbyte

    ;; Add dimension names 'data_new'
    dimNames = getfilevardims(fid, variable)
    cnt = cn + "_dimNames"
    data_new@$cnt$ = dimNames

    ;; ***********************************
    ;; If file is too big, save only metadata in 'data_new'
    ;; ***********************************
    if (datasize .gt. maxDataSize) then
        ;; Only on hdd, not in current memory
        cnt = cn + "_inMemo"
        data_new@$cnt$ = 0
        data_new@$cn$ = 0
        if (dimNames(0) .ne. "ncl_scalar") then
            do j = 0, rank - 1, 1
                cnt = cn + "_" + dimNames(j) + "_values"
                varCoords = fid->$dimNames(j)$
                ;; Add coordinate variables to data_new
                data_new@$cnt$ = (/varCoords/)
                cnt = cn + "_" + dimNames(j) + "_attNames"
                ;; Copy attribute for coord variables
                attNames = getvaratts(varCoords)
                data_new@$cnt$ = attNames
                do k = 0, dimsizes(attNames) - 1, 1
                    cnt = cn + "_" + dimNames(j) + "_attributes_" + attNames(k)
                    data_new@$cnt$ = varCoords@$attNames(k)$
                end do
                delete(varCoords)
                delete(attNames)
            end do
        end if
        if (rank.eq.5) then
            var0 = fid->$variable$(0, 0, 0, 0, 0)
        end if
        if (rank.eq.4) then
            var0 = fid->$variable$(0, 0, 0, 0)
        end if
        if (rank.eq.3) then
            var0 = fid->$variable$(0, 0, 0)
        end if
        if (rank.eq.2) then
            var0 = fid->$variable$(0, 0)
        end if
        if (rank.eq.1) then
            var0 = fid->$variable$(0)
        end if
        ;; Copy global attributes to data_new
        attNames = getvaratts(var0)
        cnt = cn + "_attNames"
        data_new@$cnt$ = attNames(rank:dimsizes(attNames) - 1)
        do j = rank, dimsizes(attNames) - 1, 1
            cnt = cn + "_attributes_" + attNames(j)
            data_new@$cnt$ = var0@$attNames(j)$
        end do
        delete(var0)
    else
        ;; ***********************************
        ;; If file is not too big, read everything into 'data_new'
        ;; ***********************************
        cnt = cn + "_inMemo"
        data_new@$cnt$ = 1
        data = fid->$variable$
        data_new@$cn$ = (/data/)
        if (dimNames(0) .ne. "ncl_scalar") then
            do j = 0, rank - 1, 1
                cnt = cn + "_" + dimNames(j) + "_values"

                ;; Add coordinate variables names to data_new
                data_new@$cnt$ = data&$dimNames(j)$

                ;; Copy attribute for coord variables
                cnt = cn + "_" + dimNames(j) + "_attNames"
                attNames = getvaratts(data&$dimNames(j)$)
                data_new@$cnt$ = attNames
                do k = 0, dimsizes(attNames) - 1, 1
                    cnt = cn + "_" + dimNames(j) + "_attributes_" + attNames(k)
                    data_new@$cnt$ = data&$dimNames(j)$@$attNames(k)$
                end do
                delete(attNames)
            end do
        end if

        ;; Copy global attributes to data_new
        attNames = getvaratts(data)
        cnt = cn + "_attNames"
        data_new@$cnt$ = attNames
        do j = 0, dimsizes(attNames) - 1, 1
            cnt = cn + "_attributes_" + attNames(j)
            data_new@$cnt$ = data@$attNames(j)$
        end do
        delete(data)
    end if
    delete(attNames)
    delete(dimNames)
    delete(fid)

    leave_msg(scriptname, funcname, 4)
    return(data_new)

end
