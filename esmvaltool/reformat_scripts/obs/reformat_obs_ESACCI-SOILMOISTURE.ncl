;;#############################################################################
;; REFORMAT SCRIPT FOR ESACCI-SOILMOISTURE SATELLITE DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    ftp://anon-ftp.ceda.ac.uk/neodc/esacci/
;;
;; Last access
;;    20151124
;;
;; Download and processing instructions
;;    Download the data from:
;;       
;;       
;;
;; Caveats
;;
;; Modification history
;;    20160824-A_laue_ax: added processing of volumetric soil moisture
;;                        content (sm, smStderr)
;;    20160721-A_laue_ax: read original daily files
;;                        added processing of uncertainty
;;    20150523-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "ESACCI-SOILMOISTURE"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = "/export/pa_data01/ESMVal/obs/RAW/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = "/export/pa_data02/ESMVal/obs/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period
    YR1 = 1988
    YR2 = 2008

    ;; Selected variable (standard name)
    VARS = (/"sm", "smStderr", "dos", "dosStderr"/)

    ;; Field
    FIELD = "T2Ms"

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin
    dayspermonth = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

    sm_ind = ind(VARS.eq."sm")
    smStderr_ind = ind(VARS.eq."smStderr")
    dos_ind = ind(VARS.eq."dos")
    dosStderr_ind = ind(VARS.eq."dosStderr")

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_ESACCI-SOILMOISTURE.ncl"
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

        fname = INDIR + "ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-" + sprinti("%04i", yr) \
                + sprinti("%02i", mo) + sprinti("%02i", dy) + "000000-fv02.2.nc"

        ;; Read soil moisture data
        f=addfile(fname,"r")
        xx = f->sm
        OUTPUT_sm = tofloat(xx * xx@scale_factor)
        copy_VarAtts(xx, OUTPUT_sm)
        copy_VarCoords(xx, OUTPUT_sm)
        delete(xx)

        ;; Read soil moisture uncertainty data
        xx = f->sm_uncertainty
        OUTPUT_smerr = tofloat(xx * xx@scale_factor)
        copy_VarAtts(xx, OUTPUT_smerr)
        copy_VarCoords(xx, OUTPUT_smerr)
        delete(xx)

        if (firstime) then
            ;; Read porosity data
            f=addfile(INDIR + "ESACCI-SOILMOISTURE-POROSITY_V01.1.nc","r")
            tmp = f->porosity
            zz = conform(OUTPUT_sm, tmp, (/1, 2/))
            delete(tmp)

            dims = dimsizes(OUTPUT_sm)
            nlat = dims(1)
            nlon = dims(2)
            ntime = 12

            time = new((/ntime/), double)
            smdata = new((/ntime, nlat, nlon/), float)
            smerrdata = new((/ntime, nlat, nlon/), float)
            dosdata = new((/ntime, nlat, nlon/), float)
            doserrdata = new((/ntime, nlat, nlon/), float)

            smpts = new((/ntime, nlat, nlon/), integer)
            smerrpts = new((/ntime, nlat, nlon/), integer)

            smdata = 0.0
            smerrdata = 0.0
            smpts = 0
            smerrpts = 0
            dosdata = 0.0
            doserrdata = 0.0

            smpts@_FillValue = -9999
            smerrpts@_FillValue = -9999

            tmp = var_attrib(smdata, VARS(sm_ind))
            delete(smdata)
            smdata = tmp
            delete(tmp)

            tmp = var_attrib(dosdata, VARS(dos_ind))
            delete(dosdata)
            dosdata = tmp
            delete(tmp)

            tmp = var_attrib(smerrdata, VARS(smStderr_ind))
            delete(smerrdata)
            smerrdata = tmp
            delete(tmp)

            tmp = var_attrib(doserrdata, VARS(dosStderr_ind))
            delete(doserrdata)
            doserrdata = tmp
            delete(tmp)
        end if

        OUTPUT = OUTPUT_sm * 100. / zz         ; dos
        OUTPUT_err = OUTPUT_smerr * 100. / zz  ; dosStderr

        ;; Format latitude coordinate
        OUTPUT_sm!1 = "lat"
        OUTPUT_smerr!1 = "lat"
        if (isMonotonic(OUTPUT_sm&lat) .eq. 0) then
            error_msg("f", diag_script, \
                "", "non-monotonic latitude coordinate")
        end if
        if (isMonotonic(OUTPUT_sm&lat) .eq. -1) then
            OUTPUT_sm = OUTPUT_sm(:, ::-1, :)
            OUTPUT_smerr = OUTPUT_smerr(:, ::-1, :)
            OUTPUT = OUTPUT(:, ::-1, :)
            OUTPUT_err = OUTPUT_err(:, ::-1, :)
        end if
        clat = lat_attrib(OUTPUT_sm&lat)
        delete(OUTPUT_sm&lat)
        OUTPUT_sm&lat = clat
        delete(OUTPUT_smerr&lat)
        OUTPUT_smerr&lat = clat
        OUTPUT!1 = "lat"
        OUTPUT&lat = clat
        OUTPUT_err!1 = "lat"
        OUTPUT_err&lat = clat
        delete(clat)

        ;; Format longitude coordinate
        OUTPUT_sm!2 = "lon"
        OUTPUT_smerr!2 = "lon"
        if (isMonotonic(OUTPUT_sm&lon) .eq. 0) then
            error_msg("f", diag_script, \
            "", "non-monotonic longitude coordinate")
        end if
        if (any(OUTPUT_sm&lon.lt.0.)) then
            OUTPUT_sm = lonFlip(OUTPUT_sm)
            OUTPUT_smerr = lonFlip(OUTPUT_smerr)
            OUTPUT = lonFlip(OUTPUT)
            OUTPUT_err = lonFlip(OUTPUT_err)
        end if
        clon = lon_attrib(OUTPUT_sm&lon)
        delete(OUTPUT_sm&lon)
        OUTPUT_sm&lon = clon
        delete(OUTPUT_smerr&lon)
        OUTPUT_smerr&lon = clon
        OUTPUT!2 = "lon"
        OUTPUT&lon = clon
        OUTPUT_err!2 = "lon"
        OUTPUT_err&lon = clon
        delete(clon)

        if (firstime) then
            smdata!0 = "time"
            smdata!1 = "lat"
            smdata!2 = "lon"
            smdata&lon = OUTPUT_sm&lon
            smdata&lat = OUTPUT_sm&lat
            smerrdata!0 = "time"
            smerrdata!1 = "lat"
            smerrdata!2 = "lon"
            smerrdata&lon = OUTPUT_smerr&lon
            smerrdata&lat = OUTPUT_smerr&lat
            dosdata!0 = "time"
            dosdata!1 = "lat"
            dosdata!2 = "lon"
            dosdata&lon = OUTPUT&lon
            dosdata&lat = OUTPUT&lat
            doserrdata!0 = "time"
            doserrdata!1 = "lat"
            doserrdata!2 = "lon"
            doserrdata&lon = OUTPUT_err&lon
            doserrdata&lat = OUTPUT_err&lat

            time@standard_name = "time"
            time@units = "days since 1950-01-01 00:00:00"
            time@axis = "T"
            time@long_name = "time"
            time@calendar = "standard"
            if (isatt(time, "_FillValue")) then
                delete(time@_FillValue)
            end if

            smdata&time = time
            smerrdata&time = time
            dosdata&time = time
            doserrdata&time = time

            firstime = False
        end if

        ; consider only grid cells with valid data
        smdata(mo-1, :, :) = where(ismissing(OUTPUT_sm(0, :, :)), smdata(mo-1, :, :), \
                                   smdata(mo-1, :, :) + OUTPUT_sm(0, :, :))
        dosdata(mo-1, :, :) = where(ismissing(OUTPUT(0, :, :)), dosdata(mo-1, :, :), \
                                   dosdata(mo-1, :, :) + OUTPUT(0, :, :))
        smerrdata(mo-1, :, :) = where(ismissing(OUTPUT_smerr(0, :, :)), smerrdata(mo-1, :, :), \
                                      smerrdata(mo-1, :, :) + OUTPUT_smerr(0, :, :)^2)
        doserrdata(mo-1, :, :) = where(ismissing(OUTPUT_err(0, :, :)), doserrdata(mo-1, :, :), \
                                      doserrdata(mo-1, :, :) + OUTPUT_err(0, :, :)^2)

        ; keep track of how many valid grid cells we have to calculate monthly means later
        smpts(mo-1, :, :) = where(ismissing(OUTPUT_sm(0, :, :)), smpts(mo-1, :, :), \
                                  smpts(mo-1, :, :) + 1)
        smerrpts(mo-1, :, :) = where(ismissing(OUTPUT_smerr(0, :, :)), smerrpts(mo-1, :, :), \
                                     smerrpts(mo-1, :, :) + 1)

        ; advance time step

        dy = dy + 1
        if (dy.gt.(dayspermonth(mo-1))) then
            dy = 1
            mo = mo + 1
            if (mo.gt.12) then
                ; -----------------------------------------------------------------------------------------
                ; calculate monthly means

                smpts = where(smpts.gt.0, smpts, smpts@_FillValue)
                smerrpts = where(smerrpts.gt.0, smerrpts, smerrpts@_FillValue)

                smdata = smdata / smpts
                smerrdata = sqrt(smerrdata / smerrpts)
                dosdata = dosdata / smpts
                doserrdata = sqrt(doserrdata / smerrpts)
                ; -----------------------------------------------------------------------------------------

                ;; Set global attributes
                gAtt = True
                gAtt@comment = "ESACCI satellite data for soil moisture"
                gAtt@history = "Created on " + systemfunc("date")
                gAtt@host    = systemfunc("echo $HOST")
                gAtt@user    = systemfunc("echo $USER")
                gAtt@period  = yr
                gAtt@field   = FIELD
                gAtt@tier    = TIER
                gAtt@source  = "ftp://anon-ftp.ceda.ac.uk/neodc/esacci/"
                gAtt@reference = "Not available"
                gAtt@title = OBSNAME + " satellite data reformatted for the ESMValTool"
                gAtt@conventions = "CF/CMOR"

                smdata&time@YEAR1 = yr
                smdata&time@YEAR2 = yr
                smerrdata&time@YEAR1 = yr
                smerrdata&time@YEAR2 = yr
                dosdata&time@YEAR1 = yr
                dosdata&time@YEAR2 = yr
                doserrdata&time@YEAR1 = yr
                doserrdata&time@YEAR2 = yr

                do i = 0, 11
                    time(i) = cd_inv_calendar(yr, i+1, 15, 0, 0, 0, \
                        "days since 1950-01-01 00:00:00", 0)
                end do

                smdata&time = time
                smerrdata&time = time
                dosdata&time = time
                doserrdata&time = time

                ;; Outfile
                fout = OUTDIR + "OBS_" + OBSNAME + "_sat_L3S-SSMV-COMBINED_" + FIELD + \
                    "_" + VARS(dos_ind) + "_" + sprinti("%04i", yr) + "01-" + \
                    sprinti("%04i", yr) + "12.nc"

                ;; Write variable (dos)
                write_nc(fout, VARS(dos_ind), dosdata, gAtt)
                delete(OUTPUT)

                fout = OUTDIR + "OBS_" + OBSNAME + "_sat_L3S-SSMV-COMBINED_" + FIELD + \
                    "_" + VARS(dosStderr_ind) + "_" + sprinti("%04i", yr) + "01-" + \
                    sprinti("%04i", yr) + "12.nc"

                ;; Write variable (dosStderr)
                write_nc(fout, VARS(dosStderr_ind), doserrdata, gAtt)
                delete(OUTPUT_err)

                fout = OUTDIR + "OBS_" + OBSNAME + "_sat_L3S-SSMV-COMBINED_" + FIELD + \
                    "_" + VARS(sm_ind) + "_" + sprinti("%04i", yr) + "01-" + \
                    sprinti("%04i", yr) + "12.nc"

                ;; Write variable (sm)
                write_nc(fout, VARS(sm_ind), smdata, gAtt)
                delete(OUTPUT_sm)

                fout = OUTDIR + "OBS_" + OBSNAME + "_sat_L3S-SSMV-COMBINED_" + FIELD + \
                    "_" + VARS(smStderr_ind) + "_" + sprinti("%04i", yr) + "01-" + \
                    sprinti("%04i", yr) + "12.nc"

                ;; Write variable (smStderr)
                write_nc(fout, VARS(smStderr_ind), smerrdata, gAtt)
                delete(OUTPUT_smerr)

                delete(gAtt)

                ; reset variables
                smdata = 0.0
                smerrdata = 0.0
                dosdata = 0.0
                doserrdata = 0.0
                smpts = 0
                smerrpts = 0

                mo = 1
                yr = yr + 1
            end if
        end if

    end do

    leave_msg(diag_script, "", 1)

end
