;;#############################################################################
;; REFORMAT SCRIPT FOR ERA-Interim REANALYISIS DATA
;;#############################################################################
;;
;; Tier
;;    Tier 3: restricted dataset.
;;
;; Source
;;    http://apps.ecmwf.int/datasets/data/interim-full-moda/
;;
;; Last access
;;    20150327
;;
;; Download and processing instructions
;;    Select Type of level ("Pressure levels" or "Surface") and 
;;    ERA Interim Fields ("Monthly means of daily means" or "Daily").
;;    For "Daily" select only the step 0 and all times (00, 06, 12, 18).
;;    Data shall be downloaded grouped by year and variable, to allow for an
;;    easy extension of the dataset. Downloaded file shall be saved in the
;;    INPUT dir given below and named as ERA-Interim_<var>_<mean>_<year>.nc,
;;    where <var> is the ERA-Interim variable name and <mean> is either monthly
;;    or daily.
;;    A registration is required for downloading the data.
;;
;; Caveats
;;
;; Modification history
;;    20160913-A_laue_ax: added variable tcc
;;    20150820-A_laue_ax: added variables lwp, iwp, sftlf
;;    20150327-A_righ_ma: merged with ERA-Interim-surf.ncl and 
;;                        ERA-Interim-surf_daily.ncl.
;;    20140128-A_righ_ma: written.
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
    YEAR1 = 1979
    YEAR2 = 2014

    ;; Selected variable (standard name)
    VARS = (/          "psl",         "ta",           "zg",\
                       "hus",         "ua",           "va",\
                       "tas",         "ps",          "tos",\
                      "tauu",       "tauv",        "clwvi",\
                     "clivi",      "sftlf",          "clt"/)
     
    ;; Name in the raw data
    NAME = (/          "msl",          "t",            "z",\
                         "q",          "u",            "v",\
                       "t2m",         "sp",          "sst",\
                      "iews",       "inss",      "p56.162",\
                   "p57.162",        "lsm",          "tcc"/)
          
    ;; Units in the input  data and corresponding conversion factor
    ;; to CMOR standard units
    RAWUNITS = (/       "Pa",          "K",   "m**2 s**-2",\
                 "kg kg**-1",    "m s**-1",      "m s**-1",\
                         "K",         "Pa",            "K",\
                   "N m**-2",    "N m**-2",     "kg m**-2",\
                  "kg m**-2",    "(0 - 1)",      "(0 - 1)"/)  
    
    CONVERSION = (/       1.,           1.,     1./9.80665,\
                          1.,           1.,             1.,\
                          1.,           1.,             1.,\
                          1.,           1.,             1.,\
                          1.,         100.,           100./)
     
    ;; Fields
    FLDARR = (/       "T2Ds",        "T3M",          "T3M",\
                       "T3M",        "T3M",          "T3M",\
                      "T2Ms",       "T2Ms",        "TO2Ms",\
                     "TO2Ms",      "TO2Ms",         "T2Ms",\
                      "T2Ms",       "T2Ms",         "T2Ms"/)

    ;; Initialize global variable
    FIELD = ""

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_ERA-Interim.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent_esmval(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Loop over variables
    do vID = 0, dimsizes(VARS) - 1

        ;; Set field
        FIELD = FLDARR(vID)

        info_output("", verbosity, 1)
        info_output("Processing " + VARS(vID) + " (" + FIELD + ") ", \
                    verbosity, 1)

        fill_start = 0

        ;; Create timeseries
        do yy = YEAR1, YEAR2

            info_output("  " + yy, verbosity, 1)

            ; land sea mask is invariant

            if (VARS(vID).eq."sftlf") then
                fname = INDIR + "ERA-Interim_" + NAME(vID) + ".nc"
                f = addfile(fname, "r")
            else
                ;; Read file
                if (isStrSubset(FIELD, "M")) then
                    fname = INDIR + "ERA-Interim_" + NAME(vID) + "_monthly_" + \
                        yy + ".nc"
                else if (isStrSubset(FIELD, "D")) then
                    fname = INDIR + "ERA-Interim_" + NAME(vID) + "_daily_" + \
                        yy + ".nc"
                end if
                end if

                f = addfile(fname, "r")

                ;; Check time coordinate
                date = cd_calendar(f->time, 0)
                if (any(date(:, 0).ne.yy)) then
                    error_msg("f", diag_script, "", "time coordinate not " + \
                              "corresponding to selected year")
                end if
                delete(date)
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

            ; -----------------------------------------------------------------

            ; land area fraction (sftlf) has only only record
            ; ---> add missing records and create corresponding time
            ;      variable

            if (VARS(vID).eq."sftlf") then  ; create 12 records in total per year
                out = tmp2
                do i = 1, 11
                    tmp3 = array_append_record(out, tmp2, 0)
                    delete(out)
                    out = tmp3
                    delete(tmp3)
                end do
                delete(tmp2)
                tmp2 = out
                delete(out)

                ; fill time records accordingly

                time = new(12, integer)

                dd = 15
                hh = 0
                mn = 0
                sc = 0
                units  = "days since 1950-01-01 00:00:00"

                do i = 1, 12
                    mm = i
                    time(i-1) = \
                        toint(cd_inv_calendar(yy,mm,dd,hh,mn,sc,units,0))
                end do

                tmp2&time = time
                tmp2&time@units = units
            end if

            ; -----------------------------------------------------------------

            ; =================================================================
            ; --- lwp is not a cmor variable ---> calculate clwvi = lwp+iwp ---
            ; =================================================================

            if (VARS(vID).eq."clwvi") then
                ;; Read 2nd file containing iwp (variable "p57.162")
                iwpname = "p57.162"

                if (isStrSubset(FIELD, "M")) then
                    fname2 = INDIR + "ERA-Interim_" + iwpname \
                             + "_monthly_" + yy + ".nc"
                else if (isStrSubset(FIELD, "D")) then
                    fname2 = INDIR + "ERA-Interim_" + iwpname \
                             + "_daily_" + yy + ".nc"
                end if
                end if

                f2 = addfile(fname2, "r")

                ;; Check time coordinate of 2nd file
                date = cd_calendar(f2->time, 0)
                if (any(date(:, 0).ne.yy)) then
                    error_msg("f2", diag_script, "", "time coordinate not " + \
                              "corresponding to selected year")
                end if
                delete(date)

                ;; Read iwp from 2nd file
                xx2 = f2->$iwpname$

                ;; Unpack variable according to metadata information
                tmp1a = todouble(xx2)
                tmp1a = tmp1a * xx2@scale_factor + xx2@add_offset
                tmp2a = tofloat(tmp1a)
                copy_VarMeta(xx2, tmp2a)
                delete(xx2)
                delete(tmp1a)

                ; add iwp to lwp to calculate cmor variable "clwvi"

                tmp2 = tmp2 + tmp2a
                delete(tmp2a)

            end if  ; if variable = clwvi

            ; =================================================================

            ;; Calculate daily mean
            if (isStrSubset(FIELD, "D")) then
                if (isStrSubset(FIELD, "T3")) then
                    dmean = tmp2(::4, :, :, :)
                    do ii = 0, dimsizes(dmean&time) - 1
                        dmean(ii, :, :, :) = \
                            dim_avg_n_Wrap(tmp2(4*ii:4*ii + 3, :, :, :), 0)
                    end do
                else
                    dmean = tmp2(::4, :, :)
                    do ii = 0, dimsizes(dmean&time) - 1
                        dmean(ii, :, :) = \
                            dim_avg_n_Wrap(tmp2(4*ii:4*ii + 3, :, :), 0)
                    end do
                end if
                delete(tmp2)
                tmp2 = dmean
                delete(dmean)    
            end if


            ;; Add to time series
            if (.not.isdefined("output")) then
                output = tmp2
            else
                tmp3 = array_append_record(output, tmp2, 0)
                delete(output)
                output = tmp3
                delete(tmp3)
            end if
            delete(tmp2)

        end do

        ;; Convert units
        if (output@units.ne.RAWUNITS(vID)) then
            error_msg("fatal",\
                      diag_script,\
                      "",\
                      "unexpected units for " + NAME(vID))
        else
            output = output * CONVERSION(vID)
        end if

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

        ;; Format vertical coordinate
        if (isStrSubset(FIELD, "T3")) then
            output!1 = "plev"
            if (isMonotonic(output&plev) .eq. 0) then
                error_msg("fatal",\
                          diag_script,\
                          "",\
                          "non-monotonic vertical coordinate")
            end if
            if (isMonotonic(output&plev).eq.1) then
                output = output(:, ::-1, :, :)
            end if
            if (output&plev@units.eq."millibars") then
                output&plev = output&plev * 100
            else
                error_msg("fatal",\
                          diag_script,\
                          "",\
                          "unexpected units for plev " + output&plev@units)
            end if
            cplev = plev_attrib(output&plev)
            delete(output&plev)
            output&plev = cplev
            delete(cplev)
        end if

        ;; Format latitude coordinate
        if (isStrSubset(FIELD, "T3")) then
            output!2 = "lat"
            if (isMonotonic(output&lat) .eq. 0) then
                error_msg("fatal",\
                          diag_script,\
                          "",\
                          "non-monotonic latitude coordinate")
            end if
            if (isMonotonic(output&lat) .eq. -1) then
                output = output(:, :, ::-1, :)
            end if
        else
            output!1 = "lat"
            if (isMonotonic(output&lat) .eq. 0) then
                error_msg("fatal",\
                          diag_script,\
                          "",\
                          "non-monotonic latitude coordinate")
            end if
            if (isMonotonic(output&lat) .eq. -1) then
                output = output(:, ::-1, :)
            end if
        end if
        clat = lat_attrib(output&lat)
        delete(output&lat)
        output&lat = clat
        delete(clat)

        ;; Format longitude coordinate
        if (isStrSubset(FIELD, "T3")) then
            output!3 = "lon"
        else
            output!2 = "lon"
        end if
        if (isMonotonic(output&lon) .eq. 0) then
            error_msg("fatal",\
                      diag_script,\
                      "",\
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
        gAtt@source        = "http://apps.ecmwf.int/datasets/data/" + \
            "interim_full_moda/"
        gAtt@reference     = "Dee, D. P. et al., Q. J. Roy. Meteor. Soc., " + \
            "137, 553-597, doi:10.1002/qj.828, 2011."
        gAtt@title = OBSNAME + \
            " reanalysis data reformatted for the ESMValTool"
        gAtt@conventions = "CF/CMOR"

        ;; Check fill values
        fill_end = num(ismissing(output))
        if (fill_start .ne. fill_end) then
            error_msg("fatal",\
                      diag_script,\
                      "",\
                      "missing values lost during conversion")
        end if
        delete(fill_start)
        delete(fill_end)

        ;; Outfile
        fout = OUTDIR + "OBS_" + OBSNAME + "_reanaly_1_" + FIELD + \
            "_" + VARS(vID) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

        ;; Add height coordinate to tas variable (required by the new backend)
        if (VARS(vID).eq."tas") then
            output@coordinates = "height"
        end if

        ;; Write variable
        write_nc(fout, VARS(vID), output, gAtt)
        delete(gAtt)
        delete(output)

        ;; Add height coordinate to tas variable (required by the new backend)
        if (VARS(vID).eq."tas") then
            height = 2.d
            height!0 = "ncl_scalar"
            height@units = "m"
            height@axis = "Z"
            height@positive = "up"
            height@long_name = "height"
            height@standard_name = "height"
            w = addfile(fout, "w")
            w->height = height
            delete(w)
        end if

    end do

    leave_msg(diag_script, "", 1)

end
