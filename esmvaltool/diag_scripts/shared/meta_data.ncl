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

;load "./interface_scripts/messaging.ncl"
load "./interface_data/ncl.interface"

;;#############################################################################
;undef("ESMValMD")
procedure ESMValMD(outfile:string,
                   caption:string,
                   statistics[*]:string,
                   domain:string,
                   plottype:string,
                   diag_description:string,
                   diag_authors[*]:string,
                   diag_ref[*]:string)
;;
;; Arguments
;;    outfile: file name of the figure (including path)
;;    caption: figure caption
;;    statistics: list of metrics/statistics calculated
;;    domain: of the data displayed
;;    plottype: plot type of the figure
;;    diag_description: description of the diagnostic
;;    diag_authors: authors related to the diagnostic
;;    diag_ref: references related to the diagnotic
;;
;; Description
;;    Creates a yaml file with meta data from the diagnostic
;;
;; Caveats
;;
;; References
;;
local funcname, scriptname, verbosity
begin

    funcname = "ESMValMD"
    scriptname = "diag_scripts/shared/lib/meta_data.ncl"
    enter_msg(scriptname, funcname)

    m_tags = "M_" + modelnames
    v_tags = "V_" + varnames
    tags_plus  = array_append_record(tags, m_tags, 0)
    tags_plus := array_append_record(tags_plus, v_tags, 0)
    delete(m_tags)
    delete(v_tags)


    str = new(8, string)
    str0 = "outfile: " + outfile 
    str1 = "caption: " + caption
    n = dimsizes(statistics)
    str2 = "statistics: "
    do i = 0, n - 2
        str2 = str2 + "- " + statistics(i) + ","
    end do
    str2 = str2 + "- " + statistics(n - 1)
    str3 = "domain: " + domain 
    str4 = "plottype: " + plottype
    str5 = "diag_description: " + diag_des
    n = dimsizes(diag_authors)
    str6 = "diag_authors: "
    do i = 0, n - 2
        str6 = str6 + "- " + diag_authors(i) + ","
    end do
    str6 = str6 + "- " + diag_authors(n - 1)
    n = dimsizes(diag_refs)
    str7 = "diag_refs: "
    do i = 0, n - 2
        str7 = str7 + "- " + diag_refs(i) + ","
    end do
    str7 = str7 + "- " + diag_refs(n - 1)

    str := array_append_record(str6,str7,0) 
    str := array_append_record(str5,str,0) 
    str := array_append_record(str4,str,0) 
    str := array_append_record(str3,str,0) 
    str := array_append_record(str2,str,0) 
    str := array_append_record(str1,str,0) 

    yaml_file = filename + "_diag_prov.yml"
    asciiwrite(yaml_file, str)
    delete(str)

    leave_msg(scriptname, funcname)

end
