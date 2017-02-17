;
;  This function takes in a data structure(data) with data for multiple models
;  that was created in read_data.ncl and extracts the data for one model.
;
;  2008-07-25 CAF
;  2008-11-06 CAF
;     added ensemble id
;  2009-01-12 CAF
;     fixed problem with same model but different time period
;  2009-01-28 CAF
;     added copy for coord. attributes
;  2010-11-17 HS
;     added parameters (dim_id,lw_limit,up_limit)
;     dim_id: all_data(-1), dim_index (e.g., 0), or dim_name (e.g., "time")
;     lw_limit/up_limit: lower/upper limits, either index or values
;     depends on the attribute (inMemo), data gets read from memory or file.
; 2012-08-28 ME
;     Added inline documentation
;

; Loading temporary file. This file is updated on the fly by Python and
; used for Python->NCL communication
load "./interface_data/ncl.interface"

load "./interface_scripts/info_output.ncl"
load "./interface_scripts/get_varName_from_file.ncl"

undef("extract_data_extended")
function extract_data_extended(index [1] : integer,
                               var [1] : string,
                               data_pointer [1] : logical,
                               dim_id,          ; Either integer or string
                               lw_limit [1] : numeric,
                               up_limit [1] : numeric)
;;                    return_val : numeric
;; Arguments:
;;    index    - Index to current "dictionary keys" (see below)
;;    data_pointer  - A pointer to the data which will be extracted
;;    dim_id   - Which dimenstion to fetch, either an integer or a dimension
;;               name (-1=all dimensions)
;;    lw_limit - Lower limit of values to fetch. If an interger, lw_limit is
;;               interpreted as an lower index in the array. If any other
;;               numeric type it is interpreted as a lower cut off value. For
;;               dimension 'time' the expected syntax is YYYY.fraction_of_year.
;;               Set to zero if not used.
;;    up_limit - Similar to lw_limit
;; Return value:
;;    data_new - The actual data array with all metadata
;; Description:
;;    Uses the variable 'dictkeys' defined in the 'interface_data/ncl.interface'
;;    file to access and extract data from variable 'data'. 'dim_id' and
;;    'lw_limit'/'up_limit' are used to to fetch a subset of the requested data.
;;    See description in arguments above for further details. The
;;    'data_pointer' is effectively a dictionary with the following
;;     key-value pairs.
;;        KEY:                                           VALUE
;;        -------------------------------------------------------------
;;        data_pointer@<filepath>_filepath               --  a string with the file path
;;        data_pointer@<filepath>_dimNames               --  a list with the dimension names
;;        data_pointer@<filepath>_inMemo                 -- =0 if data is on disk, =1 if data is in memory
;;        data_pointer@<filepath>                        -- if inMemo set to 1, this attribute contains the data
;;        data_pointer@<filepath>_dimName(j)_values      -- values of dimension 'j'
;;        data_pointer@<filepath>_dimName(j)_attNames    -- attribute names of dimension 'j'
;;        data_pointer@<filepath>_dimName(j)_attributes_attNames(k)
;;                                                   -- attribute 'k' of dimension 'j'
;;        data_pointer@<filepath>_attributes_attNames(j) -- global attribute number 'j'
local cn, cnt, dimNames, rank, inMemo, dim_id_num, lw_idx, up_idx, \
    dim_val, attNames, dim_val_idx, dim_val_time, modelfile, \
    variable, temp_data, i, j
begin

    funcname = "extract_data_extended"

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering extract_data.ncl", verbosity, 6)

    ;; 'cn' is the basic building block for the keys used in
    ;; the 'data_pointer' dictionary
    cn = interface_get_dictkeys(var, index)  ;; 'dictkeys' is loaded from
                                             ;; 'interface_data/ncl.interface'

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
                error_msg("f", "extract_data.ncl", funcname, \
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
                        cnt = cn + "_" + dimNames(dim_id_num) + "_attributes_" + attNames(i)
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
                    error_msg("f", "extract_data.ncl", funcname, \
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
            if (dim_id_num .eq. -1) then
                data_new = fm->$variable$   ; fetch all data
            else                       ; fetch a subset defined by lw_idx/up_dix
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
            if (dim_id_num .eq. -1) then
                data_new = temp_data        ; fetch all data
            else                       ; fetch a subset defined by lw_idx/up_dix
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

    info_output(">>>>>>>> Leaving extract_data.ncl", verbosity, 6)
    return(data_new)
end

undef("extract_data")
function extract_data(index [1] : integer,
                      data_pointer [1] : logical,
                      dim_id,          ; Either integer or string
                      lw_limit [1] : numeric,
                      up_limit [1] : numeric)
;;                    return_val : numeric
; Description
;    Wrapper for the extended version below. Calling the extended
;    version directly is necessary when different variables are
;    linked to different input files.
local data_new
begin
    data_new = extract_data_extended(index, "default", data_pointer, dim_id, lw_limit, up_limit)
    return data_new
end

