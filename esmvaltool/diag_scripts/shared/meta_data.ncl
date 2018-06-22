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

load "interface_scripts/logging.ncl"
load "interface_scripts/interface.ncl"

;;#############################################################################
undef("ESMValMD")
procedure ESMValMD(outfile:string,
                   caption:string,
                   statistics[*]:string,
                   domain:string,
                   plottype:string,
                   diag_des:string,
                   diag_authors[*]:string,
                   diag_refs[*]:string)
;;
;; Arguments
;;    outfile: file name of the figure (including path)
;;    caption: figure caption
;;    statistics: list of metrics/statistics calculated
;;    domain: of the data displayed
;;    plottype: plot type of the figure
;;    diag_des: description of the diagnostic
;;    diag_authors: authors related to the diagnostic
;;    diag_refs: references related to the diagnotic
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

    file_des = "#provenance information coming from the diagnostic"
    str0 = "outfile: " + outfile 
    str1 = "caption: " + caption
    n = dimsizes(statistics)
    str2 = new(n+1, string)
    str2(0) = "statistics: "
    do i = 0, n-1
        str2(i+1) = "  - " + statistics(i) 
    end do
    str3 = "domain: " + domain 
    str4 = "plottype: " + plottype
    str5 = "diag_description: " + diag_des
    n = dimsizes(diag_authors)
    str6 = new(n+1, string)
    str6(0) = "diag_authors: "
    do i = 0, n-1
        str6(i+1) = "  - " + diag_authors(i) 
    end do
    n = dimsizes(diag_refs)
    str7 = new(n+1, string)
    str7(0) = "diag_references: "
    do i = 0, n-1
        str7(i+1) = "  - " + diag_refs(i) 
    end do

    str := array_append_record(str6,str7,0) 
    str := array_append_record(str5,str,0) 
    str := array_append_record(str4,str,0) 
    str := array_append_record(str3,str,0) 
    str := array_append_record(str2,str,0) 
    str := array_append_record(str1,str,0) 
    str := array_append_record(str0,str,0) 
    str := array_append_record(file_des,str,0) 

    str_list = [/str/]

    ;yaml_file = config_user_info@plot_dir + outfile + "_diag_prov.yml"
    yaml_file = outfile + "_diag_prov.yml"
    write_table(yaml_file, "w", str_list, "%s %s %s %s %s %s %s %s")
    delete(str)

    log_info(" write meta data to " + yaml_file)

    leave_msg(scriptname, funcname)

end
