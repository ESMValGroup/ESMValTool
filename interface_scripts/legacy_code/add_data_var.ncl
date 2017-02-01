;
; Function to create the data structure.  This is needed when deriving a field
; producing a result of a different rank, ie tropopause location
;
;  output_data output data
;  new_data input data
;  models  list of models
;  index   current model to have data added
;
; 2009-01-12 CAF created, based off of read_data.ncl
; 2010-11-30 HS  added two attributes (filepath, inMemo) to hold information
;                depends on the attribute (inMemo), data gets write to from
;                memory and/or file.
; 2012-08-28 ME  Added inline documentation
;
;

; Loading temporary files. These are updated on the fly by Python and
; used for Python->NCL communication
load "./interface_data/ncl.interface"

load "./interface_scripts/info_output.ncl"

undef("add_data_var")
procedure add_data_var(index [1] : integer,
                       output_data [1] : logical,
                       new_data : numeric,
                       variable : string)
;; Arguments:
;;    output_data - The data structure into which 'new_data' will be packed
;;    new_data  - the data to add
;;    index  - index pointing to the current model
;;    variable - the variable to add
;; Description:
;;    Adds data to the internal data structure ('output_data') used by the tool.
;;    This structure is a dictionary where the keys are attribute names to the
;;    logical variable 'output_data', pointing to their respecive attribute
;;    value. The keys-vaules are,
;;        KEY:                                           VALUE
;;        -------------------------------------------------------------
;;        output_data@<filepath>_filepath         --  a string with the
;;                                                    file path
;;        output_data@<filepath>_dimNames         --  a list with the
;;                                                    dimension names
;;        output_data@<filepath>_inMemo           -- =0 if data is on disk
;;                                                   =1 if data is in memory
;;        output_data@<filepath>                  -- if inMemo set to 1, this
;;                                                   attribute contains the data
;;        output_data@<filepath>_dimName(j)_values    -- values of dimension 'j'
;;        output_data@<filepath>_dimName(j)_attNames  -- attribute names of dimension 'j'
;;        output_data@<filepath>_dimName(j)_attributes_attNames(k)
;;                                                   -- attribute 'k' of dimension 'j'
;;        output_data@<filepath>_attributes_attNames(j) -- global attribute number 'j'
local maxDataSize, field_type, ds, rank, varType, datasize, i, \
    out_file, cn, cnt, dimNames, fout, attNames, fAtt, comment, file_atts

begin
    ;; Get environment variables
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering add_data_var.ncl", verbosity, 2)

    maxDataSize = max_data_filesize(0)
    field_type = derived_field_type

    output_data = True
    ds = dimsizes(new_data)
    rank = dimsizes(ds)

    ;; ***********************************************
    ;; Estimate the data size of the current variable
    ;; ***********************************************
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
    cn = interface_get_dictkeys_no_var(index)  ; see 'interface_data/ncl.interface'

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

    fAtt = True                         ; assign file attributes
    fAtt@fields = field_type
    comment = "created by add_data.ncl"

    file_atts = getvaratts(fAtt)
    do i = 0, dimsizes(file_atts) - 1
        if (fAtt@$file_atts(i)$ .eq. "") then
            delete(fAtt@$file_atts(i)$)
        end if
    end do
    fileattdef(fout, fAtt)            ; copy file attributes

    delete(fout)
    delete(fAtt)
    delete(file_atts)
    info_output("added file " + out_file, verbosity, 3)
    info_output(">>>>>>>> Leaving add_data_var.ncl", verbosity, 2)
end
