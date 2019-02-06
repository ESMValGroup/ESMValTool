; #############################################################################
;  GENERAL ROUTINES FOR WRITING META DATA
; #############################################################################
;  Please consider using of extending existing routines before adding new ones.
;  Check the header of each routine for documentation.
;
;  Contents:
;     procedure ESMValMD
;
; #############################################################################

load "interface_scripts/logging.ncl"
load "interface_scripts/interface.ncl"
load "interface_scripts/auxiliary.ncl"

; #############################################################################
undef("ESMValMD")
procedure ESMValMD(nc_file:string,
                   outfile:string,
                   caption:string,
                   statistics[*]:string,
                   domain:string,
                   plottype:string,
                   diag_authors[*]:string,
                   diag_refs[*]:string,
                   infiles[*]:string)
;
;  Arguments
;     nc_file: netcdf file with data related to the plot
;     outfile: file name of the figure (including path)
;     caption: figure caption
;     statistics: list of metrics/statistics calculated
;     domain: of the data displayed
;     plottype: plot type of the figure
;     diag_authors: authors related to the diagnostic
;     diag_refs: references related to the diagnotic
;
;  Description
;     Creates a yaml file with meta data from the diagnostic
;     script and save it in the plot directory
;
;  Caveats
;
;  References
;
local funcname, scriptname, verbosity
begin

    funcname = "ESMValMD"
    scriptname = "diag_scripts/shared/lib/meta_data.ncl"
    enter_msg(scriptname, funcname)

    str0 = "? " + nc_file
    n = dimsizes(infiles)
    str1 = ": ancestors: ["
    str1 = str1 + infiles(0)
    if(n.ge.1) then
      do i = 1, n-1
          str1 = str1 + ", " + infiles(i)
      end do
    end if
    str1 = str1 + "]"
    n = dimsizes(diag_authors)
    str2 = "  authors: ["
    str2 = str2 + diag_authors(0)
    if(n.ge.1) then
      do i = 1, n-1
          str2 = str2 + ", " + diag_authors(i)
      end do
    end if
    str2 = str2 + "]"
    str3 = "  caption: '" + caption + "'"
    str4 = "  domain: " + domain
    str5 = "  plotfile: " + outfile
    str6 = "  plot_type: " + plottype
    n = dimsizes(diag_refs)
    str7 = "  references: ["
    str7 = str7 + diag_refs(0)
    if(n.ge.1) then
      do i = 1, n-1
          str7 = str7 + ", " + diag_refs(i)
      end do
    end if
    str7 = str7 + "]"
    n = dimsizes(statistics)
    str8 = "  statistics: ["
    str8 = str8 + statistics(0)
    if(n.ge.1) then
      do i = 1, n-1
          str8 = str8 + ", " + statistics(i)
      end do
    end if
    str8 = str8 + "]"

    str := array_append_record(str7, str8, 0)
    str := array_append_record(str6, str, 0)
    str := array_append_record(str5, str, 0)
    str := array_append_record(str4, str, 0)
    str := array_append_record(str3, str, 0)
    str := array_append_record(str2, str, 0)
    str := array_append_record(str1, str, 0)
    str := array_append_record(str0, str, 0)

    str_list = [/str/]

    yaml_file = config_user_info@run_dir + "diagnostic_provenance.yml"
    write_table(yaml_file, "a", str_list, "%s %s %s %s %s %s %s %s %s")
    delete(str)

    log_info(" write meta data to " + yaml_file)

    leave_msg(scriptname, funcname)

end
