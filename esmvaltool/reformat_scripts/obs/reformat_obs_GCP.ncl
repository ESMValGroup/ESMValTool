;;#############################################################################
;; REFORMAT SCRIPT FOR ESRL OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    http://www.globalcarbonproject.org/carbonbudget/
;;
;; Last access
;;    20140115
;;
;; Download and processing instructions
;;    Select the following filters:
;;        Category: Carbon Budget.
;;        Parameter name: fossil fuel and cement emissions
;;       	          land-use change emissions	
;;                        atmospheric growth	
;;                        ocean sink	
;;                        land sink	
;;                        sum land+ocean sinkCarbon Dioxide.
;;        Type: Reanalysis.
;;        Frequency: Annual.
;;    Convert the data to unix format using dos2unix.
;;    This list must be manually formatted in a comma-separated table, given as
;;    FILE_STAT below.
;;    Check the period covered by the data and give it below as YEAR1 and
;;    YEAR2.
;;
;; Modification history
;;    20151001-A_wenz_sa: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "GCP"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 1959
    YEAR2 = 2011

    ;; Selected variable (standard name)
    VAR = "co2flux"

    ;; Field
    FIELD = "T0M"

end

load "./interface_scripts/logging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    diag_script = "reformat_obs_ESACCI-GCP.ncl"

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Create time coordinate
    timec = create_timec(YEAR1, YEAR2)
    datec = ut_calendar(timec, 1)

    ;; Read station info
    log_info("Reading station information")

    infile = INDIR + OBSNAME + "_2012_v1.3.txt"

    if (fileexists(infile)) then

        head = readAsciiHead(infile, "#Year")
        data = readAsciiTable(infile, 7, "float", dimsizes(head))
        delete(head)

        log_info("Processing " + infile)
        ;; Read time and extract full years only
        year = data(:, 0)
        date = toint(year)
        var = data(:, 6)
        delete(data)

        ;; Check time range
        if (min(toint(year)).lt.YEAR1 .or. max(toint(year)).gt.YEAR2) then
            error_msg("f", diag_script, "", "given time range (" + \
                                       YEAR1 + "-" + YEAR2 + ")" + \
                                   " not consistent with data (" + \
                       min(toint(year)) + "-" + max(toint(year)) + ")")
        end if
        delete(year)
        
        ;; Create time coordinate
        timec = create_timec(YEAR1, YEAR2)
        datec = cd_calendar(timec, 1)
        
        ;; Create output variable
        output = new(dimsizes(timec), float)
        output!0 = "time"
        output&time = timec
        output@_FillValue = FILL

        do ii = 0, dimsizes(var) - 1
            output(6+ii*12) = (/var(ii)/)
        end do
        delete(var)
        delete(date)

        output = where(output.lt.0., output@_FillValue, output)

        ;; Set variable attributes (CMOR)
        tmp = var_attrib(output, VAR)
        delete(output)
        output = tmp
        delete(tmp)

        ;; Reformat time coordinate
        ctime = time_attrib(output&time)
        delete(output&time)
        output&time = ctime
        delete(ctime)

        ;; Set global attributes
        gAtt = True
        gAtt@history = "Created on " + systemfunc("date")
        gAtt@host = systemfunc("echo $HOST")
        gAtt@user = systemfunc("echo $USER")
        gAtt@period = YEAR1 + "-" + YEAR2
        gAtt@field = FIELD
        gAtt@tier = TIER
        gAtt@source = "http://www.globalcarbonproject.org/carbonbudget/"
        gAtt@title = OBSNAME + " data reformatted for the ESMValTool"
        gAtt@conventions = "CF/CMOR"
        
        ;; Outfile
        fout = OUTDIR + "OBS_" + OBSNAME + "_reanaly_1" \
                      + "_" + FIELD + "_" + VAR + ".nc"
        
        ;; Write variable
        write_nc(fout, VAR, output, gAtt)
        delete(gAtt)
        
    end if

    
end
