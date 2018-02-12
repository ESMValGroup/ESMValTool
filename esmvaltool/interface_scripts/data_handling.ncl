;;#############################################################################
;; GENERAL ROUTINES FOR DATA HANDLING
;;#############################################################################
;; Please consider using of extend
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    function read_data
;;
;;#############################################################################

; Global variables to save the names of all files read by "read_data"
gFilesRead0 = new(10000, string)
gNFiles0 = 0

;;#############################################################################
undef("read_data")
function read_data(index[1]:integer,
                   variable[1]:string,
                   field_type:string)
;;
;; Arguments
;;    index: index to current input file as defined in the <variable>_info.ncl
;;           file in the preproc directory of the current variable
;;    variable: name of the variable to read
;;    field_type: the variable's field type
;;
;; Return value
;;    The actual data array with all metadata.
;;
;; Description
;;    Retrieve the variable data and metadata based on the input_file_info
;;    logical defined in the <variable>_info.ncl file automatically generated
;;    by the preprocessor and stored in the preproc/<preproc_id>/ subdirectory.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;   20180212_A-righ_ma: written.
;;
local funcname, scriptname, infile, fin
begin

    funcname = "read_data"
    scriptname = "interface_scripts/data_handling.ncl"
    enter_msg(scriptname, funcname)

    ;; Check array boundaries first
    if (index.lt.0 .or. index.gt.(dimsizes(input_file_info@filename) - 1)) then
        error_msg("f", scriptname, funcname, "model index out of range")
    end if

    infile = input_file_info@filename(index)

    ;; Check existance and supported format
    if (.not.isfilepresent(infile)) then
        error_msg("f", scriptname, funcname, "input file " + infile + \
                  " not found or not a supported format")
    end if

    ;; Check field type correspondance
    if (input_file_info@field(index) .ne. field_type) then
        error_msg("f", scriptname, funcname, "requested field type " + \
                  field_type + " does not match")
    end if

    ;; Open file
    log_info("Reading in file = " + infile)
    fin = addfile(infile, "r")

    ;; Check that variable is in the file
    if (.not. isfilevar(fin, variable)) then
        error_msg("f", scriptname, funcname, "variable " + variable + \
                  "is not in the file")
    end if

    ;; Add to the list of read in files for traceability/provenance log
    if (gNFiles0.lt.dimsizes(gFilesRead0)) then
        gFilesRead0(gNFiles0) = infile
        gNFiles0 = gNFiles0 + 1
    end if

    ; Read attributes and write info to log-file
    inquire_and_save_fileinfo(infile, input_file_info@model(index), variable)

    ;; Read data
    data = fin->$variable$

    leave_msg(scriptname, funcname)
    return(data)

end
