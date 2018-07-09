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
load "interface_scripts/auxiliary.ncl"

;;#############################################################################
undef("ESMValMD")
procedure ESMValMD(nc_file:string,
                   outfile:string,
                   caption:string,
                   statistics[*]:string,
                   domain:string,
                   plottype:string,
                   diag_des:string,
                   diag_authors[*]:string,
                   diag_refs[*]:string,
                   infiles[*]:string)
;;
;; Arguments
;;    nc_file: netcdf file with data related to the plot
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
;;    script and save it in the plot directory
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

    str0 = systemfunc("basename " + nc_file) + ":"
    str1 = "    plotfile: " + systemfunc("basename " + outfile) 
    str2 = "    caption: '" + caption + "'"
    n = dimsizes(statistics)
    str3 = new(n+1, string)
    str3(0) = "    statistics: "
    do i = 0, n-1
        str3(i+1) = "      - " + statistics(i) 
    end do
    str4 = "    domain: " + domain 
    str5 = "    plottype: " + plottype
    str6 = "    diag_description: '" + diag_des + "'"
    n = dimsizes(diag_authors)
    str7 = new(n+1, string)
    str7(0) = "    diag_authors: "
    do i = 0, n-1
        str7(i+1) = "      - " + diag_authors(i) 
    end do
    n = dimsizes(diag_refs)
    str8 = new(n+1, string)
    str8(0) = "    diag_references: "
    do i = 0, n-1
        str8(i+1) = "      - " + diag_refs(i) 
    end do
    n = dimsizes(infiles)
    str9 = new(n+1, string)
    str9(0) = "    input_files: "
    do i = 0, n-1
        str9(i+1) = "      - ???" 
        ;str8(i+1) = "      - " + infiles(i) 
    end do

    str := array_append_record(str8,str9,0) 
    str := array_append_record(str7,str,0) 
    str := array_append_record(str6,str,0) 
    str := array_append_record(str5,str,0) 
    str := array_append_record(str4,str,0) 
    str := array_append_record(str3,str,0) 
    str := array_append_record(str2,str,0) 
    str := array_append_record(str1,str,0) 
    str := array_append_record(str0,str,0) 

    str_list = [/str/]

    yaml_file = outfile + "_diag_prov.yml"
    write_table(yaml_file, "a", str_list, "%s %s %s %s %s %s %s %s %s")
    delete(str)

    log_info(" write meta data to " + yaml_file)

    leave_msg(scriptname, funcname)

end
