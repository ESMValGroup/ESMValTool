;;#############################################################################
;; REFORMAT SCRIPT FOR SRB SATELLITE DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    https://eosweb.larc.nasa.gov/project/srb/srb_table
;;
;; Last access
;;    20150331
;;
;; Download and processing instructions
;;    Select NetCDF Data products
;;        for Longwave: REL3.1 LW 3HRLY NC (GEWEX)
;;        for Shortwave: REL3.0 SW MONTHLY UTC NC (GEWEX)
;;    and Order data. This will redirect you to the NASA Reverb page.
;;    Click on "Search for granules" and add all the results to the cart.
;;    Go to the cart and order the data. Processing your order may take a few
;;    hours. A confirmation e-mail is sent when data is ready for download.
;;
;; Caveats
;;
;; Modification history
;;    20150330-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "SRB"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 1984
    YEAR2 = 2007

    ;; Selected variable (standard name)
    VARS = (/"rlutcs", "rlut", "rsut"/)

    ;; Name in the raw data
    NAME = (/"clr_lw_toa_up", "lw_toa_up", "sw_toa_up"/)

    ;; File name prefix
    FNAME = (/"srb_rel3.1_longwave_monthly_", \
              "srb_rel3.1_longwave_monthly_", \
              "srb_rel3.0_shortwave_monthly_utc_"/)

    ;; Initialize global variable
    FIELD = "T2Ms"

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_SRB.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Loop over variables
    do vID = 0, dimsizes(VARS) - 1

        info_output("", verbosity, 1)
        info_output("Processing " + VARS(vID), verbosity, 1)

        ;; Create timeseries
        do yy = YEAR1, YEAR2

            info_output("  " + yy, verbosity, 1)
            do mm = 1, 12

                fname = INDIR + FNAME(vID) + yy + sprinti("%0.2i", mm) + ".nc"
                f = addfile(fname, "r")

                ;; Read variable
                xx = f->$NAME(vID)$

                ;; Add to time series
                if (.not.isdefined("output")) then
                    output = xx
                else
                    tmp = array_append_record(output, xx, 0)
                    delete(output)
                    output = tmp
                    delete(tmp)
                end if
                delete(xx)
            end do
        end do

        ;; Set variable attributes
        tmp = var_attrib(output, VARS(vID))
        delete(output)
        output = tmp
        delete(tmp)

        ;; Format time coordinate
        ctime = time_attrib(output&time)
        delete(output&time)
        output&time = ctime
        delete(ctime)

        ;; Format latitude coordinate
        output!1 = "lat"
        if (isMonotonic(output&lat) .eq. 0) then
            error_msg("f", diag_script, "", \
                      "non-monotonic latitude coordinate")
        end if
        if (isMonotonic(output&lat) .eq. -1) then
            output = output(:, ::-1, :)
        end if
        clat = lat_attrib(output&lat)
        delete(output&lat)
        output&lat = clat
        delete(clat)

        ;; Format longitude coordinate
        output!2 = "lon"
        if (isMonotonic(output&lon) .eq. 0) then
            error_msg("f", diag_script, "", \
                      "non-monotonic longitude coordinate")
        end if
        if (any(output&lon.lt.0.)) then
            output = lonFlip(output)
        end if
        clon = lon_attrib(output&lon)
        delete(output&lon)
        output&lon = clon
        delete(clon)

        info_output("  Range: " + min(output) + "; " + \
                    max(output), verbosity, 1)

        ;; Set global attributes
        gAtt = True
        gAtt@history       = "Created on " + systemfunc("date")
        gAtt@host          = systemfunc("echo $HOST")
        gAtt@user          = systemfunc("echo $USER")
        gAtt@period        = YEAR1 + "-" + YEAR2
        gAtt@field         = FIELD
        gAtt@tier          = TIER
        gAtt@source        = \
            "https://eosweb.larc.nasa.gov/project/srb/srb_table"
        gAtt@reference = "GEWEX-news: February 2011, vol. 21, No. 1, 2011."
        gAtt@title = OBSNAME + \
            " satellite data reformatted for the ESMValTool"
        gAtt@conventions = "CF/CMOR"

        ;; Outfile
        fout = OUTDIR + "OBS_" + OBSNAME + "_sat_1_" + FIELD + \
            "_" + VARS(vID) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

        ;; Write variable
        write_nc(fout, VARS(vID), output, gAtt)
        delete(gAtt)
        delete(output)

    end do

    leave_msg(diag_script, "", 1)

end
