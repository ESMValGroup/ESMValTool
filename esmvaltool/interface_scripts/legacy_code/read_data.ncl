;
; 2008-06-27 CAF  function to read in date from netcdf files
;                 currently just a place holder
; 2008-07-18 CAF  modified to read in climo files from AMWG diag
; 2008-07-25 CAF  modified to create a structure that allows
;                 different dimensioned data to be read in.
; 2008-08-27 CAF  Improved file search some
; 2008-10-28 CAF  Added field_number
; 2008-11-06 CAF  Added ensemble id
; 2009-01-12 CAF  fixed problem with same model but different time period
; 2009-01-20 CAF  added coord variable attributes copying to structure
; 2009-01-27 CAF  change filename and climo directory structure
; 2009-03-10 CAF  if double there turn float value
; 2009-05-22 CAF  field_number that is passed in can have a format of
;                 *2*s *2*z *3*, then the C or T replace the first *
;                 I or M replaced the second *.  This is based on the
;                 original field number.  Leaving out the * forces
;                 that character to be used.
; 2009-05-22 CAF  check for filenames with and without date range in filename
; 2009-05-22 CAF  fixed deletion error
; 2009-05-27 CAF  error check for out of date range
; 2009-05-28 CAF  fixed error check for out of date range
; 2009-07-18 CAF  fixes from some models
; 2010-11-17 HS   added two attributes (filepath, inMemo) to hold information
;                 data exceeds maxDataSize, then don't get loaded at this stage
; 2010-11-26 HS   file formats problems are dealt in earlier stages now.
; 2012-08-28 ME   Added inline documentation
;
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"

; Loading temporary files. These are updated on the fly by Python and
; used for Python->NCL communication
load "./interface_data/ncl.interface"

load "./interface_scripts/info_output.ncl"

undef("read_data")
function read_data(index [1] : integer,
                   variable [1] : string,
                   field_number : string)
;;                 return_val [1] : logical
;; Arguments:
;;    index    - index to current infile defined in the
;;               'interface_data/ncl.interface'-file
;;    variable - the variable to read
;;    field_number - the type of data requested (see the 'doc/*.pdf'
;;                   for further details.
;; Description:
;;    Uses the 'interface_get_fullpath' function defined in the 'interface_data/ncl.interface'
;;    file to find the current indata file. Read the indata file into the
;;    'data_new' structure. If the indata file is considered too big, only
;;    the meta data is read into the structure 'data_new'.
;;
;;    The data_new structure is described below. It is basically a
;;    dictionary where the keys are attribute names to the logical variable
;;    'data_new', pointing to the attribute values.
;;    The keys-values are,
;;        KEY:                                           VALUE
;;        -------------------------------------------------------------
;;        data_new@<filepath>_filepath               --  a string with the file path
;;        data_new@<filepath>_dimNames               --  a list with the dimension names
;;        data_new@<filepath>_inMemo                 -- =0 if data is on disk, =1 if data is in memory
;;        data_new@<filepath>                        -- if inMemo set to 1, this attribute contains the data
;;        data_new@<filepath>_dimName(j)_values      -- values of dimension 'j'
;;        data_new@<filepath>_dimName(j)_attNames    -- attribute names of dimension 'j'
;;        data_new@<filepath>_dimName(j)_attributes_attNames(k)
;;                                                   -- attribute 'k' of dimension 'j'
;;        data_new@<filepath>_attributes_attNames(j) -- global attribute number 'j'
local attNames,         cn,        cnt,\
                                    data,\
    data_new,         datasize,  dimNames,\
    dimt,             fid,       field_number,\
    in_file,          j,         k,\
    maxDataSize,      rank,      temp,\
    temp1,            var0,      varCoords,\
    varType,          verbosity, curr_var_idx
begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering read_data.ncl", verbosity, 4)

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
    info_output("Reading in file = " + in_file, verbosity, 1)
    fid = addfile(in_file, "r")

    ;; 'cn' is the basic building block for the keys used in the
    ;; 'data_new' dictionary
    cn = interface_get_dictkeys(variable, index)  ;; 'dictkeys' is loaded from
                                                  ;; 'interface_data/ncl.interface'
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
                data_new@$cnt$ = (/varCoords/)  ; Add coordinate variables to data_new
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
        do j = rank, dimsizes(attNames)-1, 1
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
    info_output(">>>>>>>> Leaving read_data.ncl", verbosity, 4)
    return(data_new)
end
