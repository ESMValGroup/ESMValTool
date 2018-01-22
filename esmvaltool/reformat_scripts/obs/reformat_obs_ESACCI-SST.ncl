;;#############################################################################
;; REFORMAT SCRIPT FOR ESACCI-SST SATELLITE DATA (analysis)
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    ftp://anon-ftp.ceda.ac.uk/neodc/esacci/sst/data/lt/Analysis/L4/v01.1/
;;
;; Last access
;;    20160818
;;
;; Download and processing instructions
;;    Download the data from the URL given above retaining the subdirectory
;;    structure, i.e. yyyy/mm/dd/*.nc
;;
;; Caveats
;;
;; Modification history
;;    20160818-A_laue_ax: written
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "ESACCI-SST"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period
    YR1 = 1992
    YR2 = 2010

    ;; Selected variable (standard name)
    VARS = (/"ts", "tsStderr"/)

    ;; Field
    FIELD = "T2Ms"

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin
    dayspermonth = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)


    verbosity = diag_script_info@verbosity
    diag_script = "reformat_obs_ESACCI-SST.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    yr = YR1
    mo = 1
    dy = 1

    firstime = True

    ;; Loop over years

    do while (yr.le.YR2)

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

        print(sprinti("%04i",yr)+"-"+sprinti("%02i",mo)+"-"+sprinti("%02i",dy))

        fname = INDIR + "lt/Analysis/L4/v01.1/" + sprinti("%04i", yr)   \
                + "/" + sprinti("%02i", mo) + "/" + sprinti("%02i", dy) + "/" \
                + sprinti("%04i", yr) + sprinti("%02i", mo) + sprinti("%02i", dy) \
                + "120000-ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_LT-v02.0-fv01.1.nc"

        ;; Read sea surface temperature data
        f=addfile(fname,"r")
        xx = f->analysed_sst
        OUTPUT = tofloat(xx * xx@scale_factor + xx@add_offset)
        if (firstime) then
            copy_VarAtts(xx, OUTPUT)
            copy_VarCoords(xx, OUTPUT)
        end if
        delete(xx)

        ;; Read sea surface temperature uncertainty data
        xx = f->analysis_error
        OUTPUT2 = tofloat(xx * xx@scale_factor + xx@add_offset)
        if (firstime) then
            copy_VarAtts(xx, OUTPUT2)
            copy_VarCoords(xx, OUTPUT2)
        end if
        delete(xx)

        ; mask anything that is not open sea water (i.e. land, ice, lakes)

        xxmask = f->mask
        OUTPUT = where(xxmask.eq.1, OUTPUT, OUTPUT@_FillValue)
        OUTPUT2 = where(xxmask.eq.1, OUTPUT2, OUTPUT2@_FillValue)

        ; debug output
        ;system ("rm debug.nc")
        ;debugfile = addfile("debug.nc","c")
        ;debugfile->mask = xxmask
        ;debugfile->output = OUTPUT
        ;debugfile->output2 = OUTPUT2
        ;exit

        delete(xxmask)

        ;fo=addfile("debug.nc","c")
        ;fo->sst=OUTPUT
        ;exit

        if (firstime) then
            dims = dimsizes(OUTPUT)
            nlat = dims(1)
            nlon = dims(2)
            ntime = 12

            time = new((/ntime/), double)
            sstdata = new((/ntime, nlat, nlon/), float)
            ssterrdata = new((/ntime, nlat, nlon/), float)

            ; use unsigned bytes to save some memory
            sstpts = new((/ntime, nlat, nlon/), ubyte)
            ssterrpts = new((/ntime, nlat, nlon/), ubyte)

            sstdata = 0.0
            ssterrdata = 0.0
            sstpts = 0B
            ssterrpts = 0B

            sstpts@_FillValue = 255B
            ssterrpts@_FillValue = 255B

            tmp = var_attrib(sstdata, VARS(0))
            delete(sstdata)
            sstdata = tmp
            delete(tmp)

            tmp = var_attrib(ssterrdata, VARS(1))
            delete(ssterrdata)
            ssterrdata = tmp
            delete(tmp)

            ;; Format latitude coordinate
            OUTPUT!1 = "lat"
            OUTPUT2!1 = "lat"
            if (isMonotonic(OUTPUT&lat) .eq. 0) then
                error_msg("f", diag_script, \
                    "", "non-monotonic latitude coordinate")
            end if
            if (isMonotonic(OUTPUT&lat) .eq. -1) then
                swap_lats = True
            else
                swap_lats = False
            end if
            clat = lat_attrib(OUTPUT&lat)
            delete(OUTPUT&lat)
            OUTPUT&lat = clat
            delete(OUTPUT2&lat)
            OUTPUT2&lat = clat
            delete(clat)

            ;; Format longitude coordinate
            OUTPUT!2 = "lon"
            if (isMonotonic(OUTPUT&lon) .eq. 0) then
                error_msg("f", diag_script, \
                "", "non-monotonic longitude coordinate")
            end if
            if (any(OUTPUT&lon.lt.0.)) then
                swap_lons = True
            else
                swap_lons = False
            end if
            clon = lon_attrib(OUTPUT&lon)
            delete(OUTPUT&lon)
            OUTPUT&lon = clon
            delete(OUTPUT2&lon)
            OUTPUT2&lon = clon
            delete(clon)

            sstdata!0 = "time"
            sstdata!1 = "lat"
            sstdata!2 = "lon"
            sstdata&lon = OUTPUT&lon
            sstdata&lat = OUTPUT&lat
            ssterrdata!0 = "time"
            ssterrdata!1 = "lat"
            ssterrdata!2 = "lon"
            ssterrdata&lon = OUTPUT2&lon
            ssterrdata&lat = OUTPUT2&lat

            time@standard_name = "time"
            time@units = "days since 1950-01-01 00:00:00"
            time@axis = "T"
            time@long_name = "time"
            time@calendar = "standard"
            if (isatt(time, "_FillValue")) then
                delete(time@_FillValue)
            end if

            sstdata&time = time
            ssterrdata&time = time

            firstime = False
        end if

        ; consider only grid cells with valid data
        sstdata(mo-1, :, :) = where(ismissing(OUTPUT(0, :, :)), sstdata(mo-1, :, :), \
                                   sstdata(mo-1, :, :) + OUTPUT(0, :, :))
        ssterrdata(mo-1, :, :) = where(ismissing(OUTPUT2(0, :, :)), ssterrdata(mo-1, :, :), \
                                      ssterrdata(mo-1, :, :) + OUTPUT2(0, :, :)^2)
        ; keep track of how many valid grid cells we have to calculate monthly means later
        sstpts(mo-1, :, :) = where(ismissing(OUTPUT(0, :, :)), sstpts(mo-1, :, :), \
                                  sstpts(mo-1, :, :) + 1B)
        ssterrpts(mo-1, :, :) = where(ismissing(OUTPUT2(0, :, :)), ssterrpts(mo-1, :, :), \
                                     ssterrpts(mo-1, :, :) + 1B)

        ; advance time step

        dy = dy + 1
        if (dy.gt.(dayspermonth(mo-1))) then
            dy = 1
            mo = mo + 1
            if (mo.gt.12) then

                ; -----------------------------------------------------------------------------------------
                ; calculate monthly means

                sstpts = where(sstpts.gt.0, sstpts, sstpts@_FillValue)
                ssterrpts = where(ssterrpts.gt.0, ssterrpts, ssterrpts@_FillValue)

                sstdata = where(ismissing(sstpts), sstdata@_FillValue, sstdata / sstpts)
                ssterrdata = sqrt(where(ismissing(ssterrpts), ssterrdata@_FillValue, ssterrdata / ssterrpts))

                ; -----------------------------------------------------------------------------------------

                if (swap_lats) then
                    sstdata = sstdata(:, ::-1, :)
                    ssterrdata = ssterrdata(:, ::-1, :)
                end if
                if (swap_lons) then
                    sstdata = lonFlip(sstdata)
                    ssterrdata = lonFlip(ssterrdata)
                end if

                ;; Set global attributes
                gAtt = True
                gAtt@comment = "ESACCI satellite data for sea surface temperature (analysis)"
                gAtt@history = "Created on " + systemfunc("date")
                gAtt@host    = systemfunc("echo $HOST")
                gAtt@user    = systemfunc("echo $USER")
                gAtt@period  = yr
                gAtt@field   = FIELD
                gAtt@tier    = TIER
                gAtt@source  = "ftp://anon-ftp.ceda.ac.uk/neodc/esacci/"
                gAtt@reference = "Merchant et al. (2014), Geoscience Data Journal, doi:10.1002/gdj3.20"
                gAtt@title = OBSNAME + " satellite data reformatted for the ESMValTool"
                gAtt@conventions = "CF/CMOR"

                sstdata&time@YEAR1 = yr
                sstdata&time@YEAR2 = yr
                ssterrdata&time@YEAR1 = yr
                ssterrdata&time@YEAR2 = yr

                do i = 0, 11
                    time(i) = cd_inv_calendar(yr, i+1, 15, 0, 0, 0, \
                        "days since 1950-01-01 00:00:00", 0)
                end do

                sstdata&time = time
                ssterrdata&time = time

                ;; Outfile
                fout = OUTDIR + "OBS_" + OBSNAME + "_sat_L4-GHRSST-SSTdepth-OSTIA-GLOB_" \
                    + FIELD + "_" + VARS(0) + "_" + sprinti("%04i", yr) + "01-" +        \
                    sprinti("%04i", yr) + "12.nc"

                ;; Write variable
                write_nc(fout, VARS(0), sstdata, gAtt)
                delete(OUTPUT)

                fout = OUTDIR + "OBS_" + OBSNAME + "_sat_L4-GHRSST-SSTdepth-OSTIA-GLOB_" \
                    + FIELD + "_" + VARS(1) + "_" + sprinti("%04i", yr) + "01-" +        \
                    sprinti("%04i", yr) + "12.nc"

                ;; Write variable
                write_nc(fout, VARS(1), ssterrdata, gAtt)
                delete(gAtt)
                delete(OUTPUT2)

                ; reset variables
                sstdata = 0.0
                ssterrdata = 0.0
                sstpts = 0B
                ssterrpts = 0B

                mo = 1
                yr = yr + 1
            end if
        end if

    end do

    leave_msg(diag_script, "", 1)

end
