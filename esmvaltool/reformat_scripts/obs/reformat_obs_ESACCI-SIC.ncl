;;#############################################################################
;; REFORMAT SCRIPT FOR ESACCI-SSMI AND -AMSR SIC DATA
;;#############################################################################
;;
;; Tier
;;   Tier 2: other freely-available dataset.
;;
;; Source
;;   ftp://anon-ftp.ceda.ac.uk/neodc/esacci/sea_ice/data/sea_ice_concentration/
;;
;; Last access
;;   20160216
;;
;; Download and processing instructions
;;   Download from above ftp server, login with "anonymous" and e-mail,
;;     download SSMI and AMSR *.nc.tar files and extract them.
;;
;; Caveats
;;   The data is only available for daily values. This script creates
;;     monthly means.
;;   Be aware that ESACCI sea ice data comes from two different instruments,
;;     SSMI and AMSRE. Both of them can be reformatted in this script.
;;
;; Modification history
;;   20160722-A_laue_ax: added processing of uncertainty
;;   20160310-A_senf_da: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = (/"ESACCI-SSMI", "ESACCI-AMSR"/)
    SATNAME = (/"SSMI", "AMSR"/)  ; corresponding to OBSNAME(iobs)

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/ESACCI-SIC/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/ESACCI-SIC/" + OBSNAME + "/"

    ;; Period
    YEAR1 = (/1992, 2003/)  ; corresponding to OBSNAME(iobs)
    YEAR2 = (/2008, 2010/)  ; corresponding to OBSNAME(iobs)

    ;; Selected variable (standard name)
    VAR = (/"sic", "sicStderr"/)

    ;; Name in the raw data
    NAME = (/"ice_conc", "total_standard_error"/)

    ;; Data is divided between the two hemispheres ("NH"/"SH")
    HEMISPHERE = (/"NH", "SH"/)

    ;; Fields
    FIELD = "T2Ms"

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_ESACCI-SIC.ncl"
    enter_msg(diag_script, "", 1)

    do iobs = 0, dimsizes(OBSNAME) - 1  ; loop over observation sets

        ;; Create output dir
        if (.not.isfilepresent_esmval(OUTDIR(iobs))) then
            system("mkdir -p " + OUTDIR(iobs))
        end if

        ;; Loop over hemispheres
        do ihem = 0, dimsizes(HEMISPHERE) - 1

            do ivar = 0, dimsizes(VAR) - 1

                info_output("", verbosity, 1)
                info_output("Processing " + OBSNAME(iobs) + " " + VAR(ivar) + " (" + \
                            FIELD + "), " + HEMISPHERE(ihem), verbosity, 1)

                time = create_timec(YEAR1(iobs), YEAR2(iobs))
                date = cd_calendar(time, 1)
                mm   = (/01, 02, 03, 04, 05, 06, 07, 08, 09, 10, 11, 12/)  ; months

                ;; Read files
                do yy = YEAR1(iobs), YEAR2(iobs)

                    info_output("Calculating year " + yy, verbosity, 1)

                    ldate = yy + sprinti("%0.2i", mm)
                    fyear = tostring(yy)
                    fname = INDIR(iobs) + "ESACCI-SEAICE-L4-SICONC-" + \
                            SATNAME(iobs) + "-" + HEMISPHERE(ihem) + \
                            "25kmEASE2-" + fyear + "0101-" + fyear + \
                            "1231-fv01.11.nc"
                    f = addfile(fname, "r")
                    xx = f->$NAME(ivar)$

                    ;; Set some standards
                    xx_tmp = tofloat(xx)
                    copy_VarMeta(xx, xx_tmp)
                    xx := xx_tmp
                    delete(xx_tmp)
                    xx = where(ismissing(xx), 1.e20, xx)
                    xx@_FillValue = 1.e20
                    xx@missing_value = 1.e20

                    ;; Correct factor
                    xx = xx / 10000.  ; data represented as fractions of 10,000

                    ;; Create monthly averages
                    if (VAR(ivar).eq."sicStderr") then
                        xx2 = xx^2
                        copy_VarMeta(xx, xx2)
                        xx_monthly = sqrt(calculate_monthly_values(xx2, "avg", 0, False))
                        delete(xx2)
                    else
                        xx_monthly = calculate_monthly_values(xx, "avg", 0, False)
                    end if

                    ;; Assign to global array
                    if (.not.isdefined("sic")) then
                        dims  = dimsizes(xx_monthly)
                        sic   = new((/dimsizes(time), dims(1), dims(2)/), \
                                    typeof(xx))
                        j     = ispan(0, dims(1) - 1, 1)
                        i     = ispan(0, dims(2) - 1, 1)
                        sic!0 = "time"
                        sic&time = time
                        sic!1 = "j"
                        sic!2 = "i"
                        sic&j = j
                        sic&i = i
                        sic&j@units     = "1"
                        sic&j@long_name = "cell index along first dimension"
                        sic&i@units     = "1"
                        sic&i@long_name = "cell index along second dimension"
                        lat   = f->lat
                        lat!0 = "j"
                        lat!1 = "i"
                        lat&j = j
                        lat&i = i
                        lon   = f->lon
                        lon!0 = "j"
                        lon!1 = "i"
                        lon&j = j
                        lon&i = i
                        sic@coordinates = "lon lat"
                        sic@units = "1"
                    end if
                    time0 = ind(toint(date).eq.min(toint(ldate)))
                    time1 = ind(toint(date).eq.max(toint(ldate)))
                    sic(time0:time1, :, :) = (/xx_monthly/)

                    ;; Clean up
                    delete([/ldate, fyear, fname, f,xx, xx_monthly, time0, time1/])

                end do  ; yy

                ;; Add global attributes
                gAtt         = True
                gAtt@history = "Created on " + systemfunc("date")
                gAtt@host    = systemfunc("echo $HOST")
                gAtt@user    = systemfunc("echo $USER")
                gAtt@period  = YEAR1(iobs) + "-" + YEAR2(iobs)
                gAtt@field   = FIELD
                gAtt@tier    = TIER
                gAtt@source  = "ftp://anon-ftp.ceda.ac.uk/neodc/esacci/" + \
                               "sea_ice/data/sea_ice_concentration/"
                gAtt@reference = "Not available"
                gAtt@title   = OBSNAME(iobs) + \
                               " satellite data reformatted for the ESMValTool"
                gAtt@conventions = "CF/CMOR"

                ;; Create and add areafile
                areacello = new((/dimsizes(sic&j), dimsizes(sic&i)/), typeof(sic))
                areacello = 6.25e+08  ; equal-area grid! (Thomas Lavergne)
                areacello!0 = "j"
                areacello!1 = "i"
                areacello&j = sic&j
                areacello&i = sic&i
                copy_VarMeta(lat, areacello)
                areacello@standard_name = "cell_area"
                areacello@long_name     = "Ocean Grid-Cell Area"
                areacello@units         = "m2"
                areacello@missing_value = sic@_FillValue
                areacello@_FillValue    = sic@_FillValue
                gAtt@areacello          = "added areacello to file, based on " + \
                                          "equal-area EASE grid [25km x 25km]"

                ;; Write file
                delete(sic&time@_FillValue)  ; check attributes
                outfile = OUTDIR(iobs) + "OBS_" + OBSNAME(iobs) + "_sat_" + \
                          HEMISPHERE(ihem) + "_" + FIELD + "_" + VAR(ivar) + \
                          "_" + YEAR1(iobs) + "01-" + YEAR2(iobs) + "12.nc"
                write_nc(outfile, VAR(ivar), sic, gAtt)
                fout = addfile(outfile, "w")
                fout->lat = lat
                fout->lon = lon
                fout->areacello = areacello
                info_output("Generated file " + outfile, verbosity, 1)

            end do ; ivar-loop

            ;; Clean up for ihem
            delete([/sic, dims, j, i, lon, lat, outfile, fout/])

        end do  ; ihem

        ;; Clean up for iobs
        delete([/time, date/])

    end do  ; iobs

    leave_msg(diag_script, "", 1)

end
