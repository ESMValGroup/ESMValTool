;;#############################################################################
;; REFORMAT SCRIPT FOR TRMM OBSERVATIONAL DATA
;;#############################################################################
;;
;; Full name of data set
;;    TRMM_3B42 Version 7: TRMM 3-Hourly 0.25 deg.
;;    TRMM and Others Rainfall Estimate Data 
;;
;; Tier
;;    Tier: 2 (other freely available data set)
;;
;; Source
;;    To cite the data in publications:
;;       Tropical Rainfall Measurement Mission Project (TRMM) (2011),
;;       TRMM 3-Hourly 0.25 deg. TRMM and Others Rainfall Estimate Data,
;;       version 7, Greenbelt, MD:Goddard Space Flight Center Distributed
;;       Active Archive Center (GSFC DAAC),
;;       Accessed <Enter User Data Access Date> at
;;       http://disc.sci.gsfc.nasa.gov/datacollection/TRMM_3B42_V7.html
;;    Key reference:
;;       Huffman, G.J., R.F. Adler, D.T. Bolvin, G. Gu, E.J. Nelkin,
;;       K.P. Bowman, Y. Hong, E.F. Stocker, D.B. Wolff: The TRMM
;;       Multi-satellite Precipitation Analysis: Quasi-Global, Multi-Year,
;;       Combined-Sensor Precipitation Estimates at Fine Scale,
;;       J. Hydrometeor., 8, 38-55, 2007.
;;
;; Last access
;;    September 2015
;;
;; Download and processing instructions
;;    Download: Go to http://disc.sci.gsfc.nasa.gov/datacollection/TRMM_3B42_V7.html
;;              and click on the tab "GetData". For easy automated download of
;;              multiple files (using wget) select "Simple Subset Wizard (SSW)"
;;              and follow the instructions.
;;    Processing: this script (reformat_obs_TRMM_3B42_3hourly.ncl)
;;
;; Caveats
;;    None
;;
;; Modification history
;;    20150904-A_laue_ax: written.
;;
;; ############################################################################

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity =  stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_TRMM_3B42_3hourly.ncl"
    enter_msg(diag_script, "", 1)

    inpath = getenv("ESMValTool_RAWOBSPATH") + "/Tier2/TRMM-3B42/"
    outpath = getenv("ESMValTool_OBSPATH") + "/Tier2/TRMM-3B42/"

    YEAR1 = 1998
    YEAR2 = 2014

    ;; Create output dir
    if (.not.isfilepresent(outpath)) then
        system("mkdir -p " + outpath)
    end if

    dayspermonth = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

    yr = YEAR1
    mo = 1
    dy = 1
    hr = 0
    nt = 0

    firstime = True

    ;; Loop over years

    do while (yr.le.YEAR2)

        leap = 0
        if ((yr % 4).eq.0) then
            leap = 1
        end if
        if ((yr % 100).eq.0) then
           leap = 0
        end if
        if ((yr % 400).eq.0) then
           leap = 1
        end if
        if (leap.eq.0) then
           dayspermonth(1) = 28
        else
           dayspermonth(1) = 29
        end if

        fname = inpath + "3B42." + sprinti("%04i", yr) + sprinti("%02i", mo) \
                + sprinti("%02i", dy) + "." + sprinti("%02i", hr) + ".7"
        if ((yr.ge.2000).and.(yr.le.2010)) then
            fname = fname + "A.HDF.Z"
        else
            fname = fname + ".HDF.Z"
        end if

        info_output("*** "+ fname + " ***", verbosity, 1)

        tmpfile = "./temp.hdf"
        system("gzip -dc <" + fname + " >" + tmpfile)
        f = addfile(tmpfile, "r")

        ; --------------------------------------------------------------------
        ; 1st time initialization
        ; --------------------------------------------------------------------

        if (firstime) then
            pr = f->precipitation
            dims = dimsizes(pr)
            nlon = dims(0)
            nlat = dims(1)
            ntime = 31 * 8

            lon = new((/nlon/), float)
            lat = new((/nlat/), float)
            time = new((/ntime/), double)

            ; define grid manually (values extracted from HDF header)

            ; Registration=CENTER
            LatitudeResolution = 0.25
            LongitudeResolution = 0.25
            NorthBoundingCoordinate = 50
            SouthBoundingCoordinate = -50
            EastBoundingCoordinate = 180
            WestBoundingCoordinate = -180

            lat = fspan(SouthBoundingCoordinate + 0.5*LatitudeResolution, \
                        NorthBoundingCoordinate - 0.5*LatitudeResolution, nlat)

            lon = fspan(WestBoundingCoordinate + 0.5*LatitudeResolution, \
                        EastBoundingCoordinate - 0.5*LatitudeResolution, nlon)

            lon@standard_name = "longitude"
            lon@units = "degrees_east"
            lon@axis = "X"
            lon@long_name = "longitude"
            lon!0 = "lon"

            lat@standard_name = "latitude"
            lat@units = "degrees_north"
            lat@axis = "Y"
            lat@long_name = "latitude"
            lat!0 = "lat"

            time@standard_name = "time"
            time@units = "days since 1950-01-01 00:00:00"
            time@axis = "T"
            time@long_name = "time"
            time@calendar = "standard"
            if (isatt(time, "_FillValue")) then
                delete(time@_FillValue)
            end if

            precip = new((/ntime, nlat, nlon/), float)
            msk = new(nlon, integer)

            precip@standard_name = "precipitation_flux"
            precip@units = "kg m-2 s-1"
            precip@long_name = "Precipitation"
            precip!0 = "time"
            precip!1 = "lat"
            precip!2 = "lon"
            precip&lon = lon
            precip&lat = lat

            precip = 0.0

            fillval = precip@_FillValue

            gAtt = True
            gAtt@host        = systemfunc("echo $HOST")
            gAtt@user        = systemfunc("echo $USER")
            gAtt@tier        = 2
            gAtt@field       = "T2Is"
            gAtt@source      = "GES DISC, http://disc.sci.gsfc.nasa.gov/datacollection/TRMM_3B42_V7.html"
            gAtt@reference   = "Huffman et al. (2007), J. Hydrometeor."
            gAtt@title       = "TRMM_3B42 Version 7" + \
                               " satellite data reformatted for the ESMValTool"
            gAtt@conventions = "CF/CMOR"

            firstime = False

        end if

        ; --------------------------------------------------------------------

        ; read and process current time step

        pr = f->precipitation

        ; swap longitudes and latidudes and convert from mm/hr to kg/m2/s

        do i = 0, nlat - 1
            msk = where(pr(:, i).ge.0.0, 1, 0)
            pr(:, i) = where(msk.eq.1, pr(:, i), 0.0)
            precip(nt, i, :) = pr(:, i) / 3600.0
        end do

        ; --------------------------------------------------------------------

        time(nt) = cd_inv_calendar(yr, mo, dy, hr, 0, 0, \
                   "days since 1950-01-01 00:00:00", 0)

        ; advance time step

        hr = hr + 3

        if (hr.ge.24) then
            hr = 0
            dy = dy + 1
            if (dy.gt.(dayspermonth(mo-1))) then
                ; save results to file (1 file per month)

                precip&time = time

                gAtt@history = "Created on " + systemfunc("date")
                gAtt@period  = yr

                datestr = sprinti("%04i", yr) + sprinti("%02i", mo) + "-" + \
                    sprinti("%04i", yr) + sprinti("%02i", mo)
                fout = outpath + "OBS_TRMM_sat_3B42-v7_T2Is_pr_" + datestr + ".nc"

                ; write variable
                write_nc(fout, "pr", precip(0:nt, :, :), gAtt)

                ; reset
                nt = -1

                dy = 1
                mo = mo + 1
                if (mo.gt.12) then
                    mo = 1
                    yr = yr + 1
                end if
            end if
        end if

        nt = nt + 1
    end do

    delete(lon)
    delete(lat)
    delete(time)
    delete(precip)
    delete(msk)

    system ("rm " + tmpfile)

    leave_msg(diag_script, "", 1)

end
