;;#############################################################################
;; GENERAL ROUTINES FOR WRITING META DATA
;;#############################################################################
;; Please consider using of extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    procedure ESMValMD
;;
;;#############################################################################

load "./interface_scripts/messaging.ncl"
load "./interface_data/ncl.interface"

;;#############################################################################
undef("ESMValMD")
procedure ESMValMD(filename:string,
                  tags[*]:string,
                  caption:string,
                  id:string)
;;
;; Arguments
;;    filename: file name of the figure (including path)
;;    tags: list of tags
;;    caption: figure caption
;;    id: id string
;;
;; Return value
;;    0: success
;;    1: error
;;
;; Description
;;    Creates meta data string and calls Pzthon function to write meta data to
;;    figure files.
;;
;; Caveats
;;
;; References
;;
;; Modification history:
;;    20170509-A_laue_ax: written.
;;
local funcname, scriptname, verbosity
begin

    funcname = "ESMValMD"
    scriptname = "diag_scripts/lib/ncl/meta_data.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 6)

    n = dimsizes(tags)
    str = new(5, string)
    str(0) = filename
    str(1) = "both"  ; type
    str(2) = ""
    do i = 0, n - 2
        str(2) = str(2) + tags(i) + ","
    end do
    str(2) = str(2) + tags(n - 1)
    str(3) = caption
    str(4) = id

    ascii_file = filename + "_list.txt"
    asciiwrite(ascii_file, str)
    delete(str)

    system("python diag_scripts/lib/python/running_MD_for_ncl_with_file.py " + ascii_file)
    system("rm " + ascii_file)

    leave_msg(scriptname, funcname, 6)

end
