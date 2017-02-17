;;#############################################################################
;; REFORMAT SCRIPT FOR ERA-Interim REANALYSIS ACCUMULATED FLUXES AT THE SURFACE
;;#############################################################################
;;
;; Tier
;;    Tier 3: restricted dataset.
;;
;; Source
;;    http://apps.ecmwf.int/datasets/data/interim-full-mnth/?levtype=sfc
;;
;; Last access
;;    20150522
;;
;; Download and processing instructions
;;    Select the time period: ideally, data should be downloaded grouped by
;;    year (i.e., one file per year), to allow for an easy extension of the
;;    dataset.
;;    The script has been written to convert 12h accumulated fluxes, so in
;;    "Select time"click "00:00:00" and "12:00:00", and in "Select step" click
;;    "12". The downloaded files will have two values per month. These values
;     are then averaged and divided by the accumulation interval to get the
;;    amount per second. Note that we assume that the amount is 0 at the start
;;    of the accumulation interval.
;;    Select the variables of interest.
;;    Save the files as ERA-Interim_surfflux_YYYY.nc and specify the path as
;;    INDIR below.
;;    Script will also calculate precip+evap (evap negative) and save it as
;;    wfpe, and total heat flux hfds by summing ssr+str+slhf+sshf
;;    A registration is required for downloading the data.
;;
;; Caveats
;;
;; Modification history
;;    20150522-A_vanu_be: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "ERA-Interim"

    ;; Tier
    TIER = 3

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 2000
    YEAR2 = 2005

    ;; Accumulation period in hours, all input fields are divided by this
    ;; period to get the amount per second
    HACC = 12

    ;; Selected variable (standard name)
    VARS = (/"pr", "evspsbl", "hfls", "hfss", "rsns", "rlns"/)

    ;; Name in the raw data
    NAME = (/"tp", "e", "slhf", "sshf", "ssr", "str"/)

    ;; Units in the input  data and corresponding conversion factor
    ;; to CMOR standard units
    RAWUNITS = (/"m", "m of water equivalent", "J m**-2", "J m**-2", \
                "J m**-2", "J m**-2"/)
    CONVERSION = (/1000., 1000., 1., 1., 1., 1./)

    ;; Field
    FIELD = "T2Ms"

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_ERA-Interim-surffluxes.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Set global attributes that are the same for every variable
    gAtt = True
    gAtt@host   = systemfunc("echo $HOST")
    gAtt@user   = systemfunc("echo $USER")
    gAtt@period = YEAR1 + "-" + YEAR2
    gAtt@field  = FIELD
    gAtt@source = "http://apps.ecmwf.int/datasets/data/" + \
        "interim-full-mnth/?levtype=sfc"
    gAtt@title = OBSNAME + \
        " reanalysis data for the ESMVal evaluation tool"
    gAtt@conventions = "CF-1.1"
    gAtt@reference     = "Dee, D. P. et al., Q. J. Roy. Meteor. Soc., " + \
        "137, 553-597, doi:10.1002/qj.828, 2011."

    ;; Loop over variables
    do vID = 0, dimsizes(VARS) - 1

        info_output("", verbosity, 1)
        info_output("Processing " + VARS(vID) + " (" + FIELD + ") ", \
                    verbosity, 1)

        fill_start = 0

        ;; Create timeseries
        do yy = YEAR1, YEAR2
            info_output("  " + yy, verbosity, 1)

            ;; Read file
            fname = INDIR + "ERA-Interim_surfflux_" + yy + ".nc"
            f = addfile(fname, "r")

            ;; Check time coordinate
            date = cd_calendar(f->time, 0)
            if (any(date(:, 0) .ne. yy)) then
                error_msg("f", get_script_name(), "", \
                          "Time coordinate not corresponding to selected year")
            end if
            npermonth = 24 / HACC
            if (dimsizes(f->time) .ne. 12 * npermonth) then
                error_msg("f", get_script_name(), "", \
                          "Unexpected number of time steps for given HACC " + \
                          "in " + fname)
            end if

            ;; Read variable
            xx = f->$NAME(vID)$

            ;; Add to missing values count
            fill_start = fill_start + num(ismissing(xx))

            ;; Unpack variable according to metadata information
            tmp1 = todouble(xx)
            tmp1 = tmp1 * xx@scale_factor + xx@add_offset
            tmp2 = tofloat(tmp1)
            copy_VarMeta(xx, tmp2)
            delete(xx)
            delete(tmp1)

            ;; Divide values over accumulation period
            tmp2 = tmp2 / (HACC * 3600)

            ;; calculate average per month from given values
            dim_sizes = dimsizes(tmp2)
            monavg = new((/12, dim_sizes(1), dim_sizes(2)/), float)
            time = new(12, integer)
            latitude = f->latitude
            longitude = f->longitude
            monavg!0 = "time"
            monavg!1 = "latitude"
            monavg!2 = "longitude"
            copy_VarAtts(tmp2, monavg)
            copy_VarAtts(f->time, time)
            copy_VarAtts(f->latitude, latitude)
            copy_VarAtts(f->longitude, longitude)
            do mm = 1, 12
                istart = (mm - 1) * npermonth
                iend = mm * npermonth - 1
                monavg(mm - 1, :, :) = dim_avg_n(tmp2(istart:iend, :, :), 0)
                time(mm - 1) = f->time(istart)
            end do
            monavg&time = time
            monavg&latitude = latitude
            monavg&longitude = longitude
            delete(tmp2)

            ;; calculate average per month and add to time series
            if (.not.isdefined("output")) then
                output = monavg
                copy_VarMeta(output, monavg)
            else
                tmp3 = array_append_record(output, monavg, 0)
                delete(output)
                output = tmp3
                delete(tmp3)
            end if
            delete(monavg)
            delete(time)
            delete(latitude)
            delete(longitude)
        end do

        ;; Convert units
        if (output@units.ne.RAWUNITS(vID)) then
            error_msg("f", get_script_name(), "", \
                      "Unnexpected units for " + NAME(vID))
        else
            output = output * CONVERSION(vID)
        end if

        ;; Set variable attributes
        tmp = var_attrib(output, VARS(vID))
        delete(output)
        output = tmp
        delete(tmp)

        ;; Format time coordinate
        info_output("  Formatting time coordinate", verbosity, 1)
        ctime = time_attrib(output&time)
        delete(output&time)
        output&time = ctime
        delete(ctime)

        ;; Format latitude coordinate
        info_output("  Formatting latitude coordinate", verbosity, 1)
        output!1 = "lat"
        if (isMonotonic(output&lat) .eq. 0) then
            error_msg("f", get_script_name(), "", \
                      "Non-monotonic latitude coordinate")
        end if
        if (isMonotonic(output&lat) .eq. -1) then
            output = output(:, ::-1, :)
        end if
        clat = lat_attrib(output&lat)
        delete(output&lat)
        output&lat = clat
        delete(clat)

        ;; Format longitude coordinate
        info_output("  Formatting longitude coordinate", verbosity, 1)
        output!2 = "lon"
        if (isMonotonic(output&lon) .eq. 0) then
            error_msg("f", get_script_name(), "", \
                      "Non-monotonic longitude coordinate")
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
        gAtt@history  = "Created on " + systemfunc("date")

        ;; Check fill values
        fill_end = num(ismissing(output))
        if (fill_start.ne.fill_end) then
            error_msg("f", get_script_name(), "", \
                      "Missing values lost during conversion")
        end if
        delete(fill_start)
        delete(fill_end)

        ;; Outfile
        fout = OUTDIR + "OBS_" + OBSNAME + "_reanaly_1_" + FIELD + \
            "_" + VARS(vID) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

        ;; Write variable
        write_nc(fout, VARS(vID), output, gAtt)

        ;; Store variables that we need lateron
        if (VARS(vID) .eq. "pr") then
            pr = output
        end if
        if (VARS(vID) .eq. "evspsbl") then
            evspsbl = output
        end if
        if (VARS(vID) .eq. "hfls") then
            hfls = output
        end if
        if (VARS(vID) .eq. "hfss") then
            hfss = output
        end if
        if (VARS(vID) .eq. "rsns") then
            rsns = output
        end if
        if (VARS(vID) .eq. "rlns") then
            rlns = output
        end if

        ;; clean
        delete(output)

    end do

    ;; Prepare wfpe (precipitation + evaporation)
    varname = "wfpe"
    info_output("", verbosity, 1)
    info_output("Calculating " + varname, verbosity, 1)

    if (.not. isvar("pr") .or. .not. isvar("evspsbl")) then
        info_output(get_script_name() + ": Precipitation and/or " + \
                    "evaporation not available, not creating " + varname, \
                    verbosity, 1)
    else
        wfpe = pr + evspsbl
        copy_VarMeta(pr, wfpe)

        ;; Set variable attributes
        tmp = var_attrib(wfpe, varname)
        delete(wfpe)
        wfpe = tmp
        delete(tmp)

        ;; Set global attributes
        gAtt@history = "Created on " + systemfunc("date")

        ;; Outfile
        fout = OUTDIR + "OBS_" + OBSNAME + "_reanaly_1_" + FIELD + "_" + \
               varname + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

        ;; Write variable
        write_nc(fout, varname, wfpe, gAtt)

        ;; clean
        delete(wfpe)
    end if
    if (isvar("pr")) then
        delete(pr)
    end if
    if (isvar("evspsbl")) then
        delete(evspsbl)
    end if

    ;; Prepare total heat flux (SSR+STR+SLHF+SSHF)
    varname = "hfds"
    info_output("", verbosity, 1)
    info_output("Calculating " + varname, verbosity, 1)

    if (.not. isvar("rsns") .or. .not. isvar("rlns") .or. \
        .not. isvar("hfss") .or. .not. isvar("hfls")) then
        info_output(get_script_name() + ": Not all variables required for" + \
                    " the calculation of the total heat flux " + \
                    "(SSR+STR+SLHF+SSHF) are available, not creating " \
                    + varname, verbosity, 1)
    else
        hfds = hfls + hfss + rsns + rlns
        copy_VarMeta(hfls, hfds)

        ;; Set variable attributes
        tmp = var_attrib(hfds, varname)
        delete(hfds)
        hfds = tmp
        delete(tmp)

        ;; Set global attributes
        gAtt@history  = "Created on " + systemfunc("date")

        ;; Outfile
        fout = OUTDIR + "OBS_" + OBSNAME + "_reanaly_1_" + FIELD + "_" + \
            varname + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

        ;; Write variable
        write_nc(fout, varname, hfds, gAtt)

        ;; clean
        delete(hfds)
    end if
    if (isvar("rsns")) then
        delete(rsns)
    end if
    if (isvar("rlns")) then
        delete(rlns)
    end if
    if (isvar("hfss")) then
        delete(hfss)
    end if
    if (isvar("hfls")) then
        delete(hfls)
    end if

    ;; clean
    delete(gAtt)

    leave_msg(diag_script, "", 1)

end
