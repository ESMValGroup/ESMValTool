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
                  id:string,
                  varnames[*],
                  modelnames[*],
                  infiles[*],
                  diag_name:string,
                  contrib_authors[*]:string)
;;
;; Arguments
;;    filename: file name of the figure (including path)
;;    tags: list of tags
;;    caption: figure caption
;;    id: id string
;;    varnames: list of variable names
;;    modelnames: list of model names
;;    infiles: list of input files (= climo files)
;;    diag_name: name of diagnostic
;;    contib_authors: list of contributing authors
;;
;; Description
;;    Creates meta data string and calls Python function to write meta data to
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

    m_tags = "M_" + modelnames
    v_tags = "V_" + varnames
    tags_plus  = array_append_record(tags, m_tags, 0)
    tags_plus := array_append_record(tags_plus, v_tags, 0)
    delete(m_tags)
    delete(v_tags)

    n = dimsizes(tags_plus)
    str = new(8, string)
    str(0) = filename
    str(1) = "both"  ; type
    str(2) = ""
    do i = 0, n - 2
        str(2) = str(2) + tags_plus(i) + ","
    end do
    str(2) = str(2) + tags_plus(n - 1)
    str(3) = caption
    str(4) = id

    n = dimsizes(infiles)
    str(5) = ""
    do i = 0, n - 2
        str(5) = str(5) + infiles(i) + ","
    end do
    str(5) = str(5) + infiles(n - 1)
    str(6) = diag_name

    n = dimsizes(contrib_authors)
    str(7) = ""
    do i = 0, n - 2
        str(7) = str(7) + contrib_authors(i) + ","
    end do
    str(7) = str(7) + contrib_authors(n - 1)

    ascii_file = filename + "_list.txt"
    asciiwrite(ascii_file, str)
    delete(str)

    system("python diag_scripts/lib/python/running_MD_for_ncl_with_file.py " + ascii_file)
;    system("rm " + ascii_file)

    leave_msg(scriptname, funcname, 6)

end
