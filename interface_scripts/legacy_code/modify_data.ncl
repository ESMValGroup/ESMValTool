;
;
;  This procedure modifies data in the data structure.
;  Only modifies data and not dimensions or attributes
;
;  2008-07-25 CAF
;  2008-11-06 CAF
;     added ensemble id
;  2009-01-12 CAF
;     fixed problem with same model but different time period
;  2010-11-29 HS
;     depends on the attribute (inMemo), data gets write to from
;     memory and/or file.
;     also variable attribute (@units) gets updated to match data content.
;  2012-05-16 HS
;     pass global attributes from initial file (fin)
;
;

; Loading temporary file. This file is updated on the fly by Python and
; used for Python->NCL communication
load "./interface_data/ncl.interface"

load "./interface_scripts/get_varName_from_file.ncl"

undef("modify_data")
procedure modify_data(index [1] : integer,
                      data_pointer [1] : logical,
                      data_new : numeric)
;; Arguments:
;;    index    - Index to current "dictionary keys" (see below)
;;    data_pointer  - A pointer to the data which will be modified
;;    data_new - The data to insert into the data_pointer structure
;; Description:
;;    Uses the variable 'dictkeys' defined in the 'interface_data/ncl.interface'
;;    file to access and modify data describeb by the variable
;;    'data_pointer'. The 'data_pointer' structure is effectively
;;    a dictionary with the following key-value pairs.
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
local verbosity, cnt, inMemo, modelfile, fout, vNames, file_atts, fin, variable, 
    delimiter, file_parts, idxLast, tempfile
begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering modify_data.ncl", verbosity, 2)

    ;; 'cn' is the basic building block for the keys used
    ;;  in the 'data_pointer' dictionary
    cn = interface_get_dictkeys_no_var(index)  ; see 'interface_data/ncl.interface'

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
    info_output("tempfile = " + str_sub_str(tempfile, pwd, "."), verbosity, 3)

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
    info_output(">>>>>>>> Leaving modify_data.ncl", verbosity, 2)
end
