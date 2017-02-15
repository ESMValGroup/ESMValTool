;;#############################################################################
;; REFORMAT SCRIPT FOR NSIDC SEA-ICE OBSERVATIONS
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    ftp://sidads.colorado.edu/pub/DATASETS/NOAA/G02202_v2/
;;
;; Last access
;;    20150708
;;
;; Download and processing instructions
;;    You can bypass the optional registration.
;;    Download instructions can be found in the README file on above website.
;;    You need to download all files under north/monthly/ and south/monthly/
;;    Downloaded files need to be merged into two files issuing
;;      'cdo mergetime seaice_conc_monthly_nh*.nc NSIDC_all_monthly_NH.nc'
;;      'cdo mergetime seaice_conc_monthly_sh*.nc NSIDC_all_monthly_SH.nc'
;;    You can then delete the downloaded files.
;;
;; Additional information
;;    The files contain multiple NSIDC datasets. This script reformats two of
;;      them simultaneously: NSIDC NASA-Team Algortihm (NSIDC-NT) and
;;      NSIDC Bootstrap Algorithm (NSIDC-BT).
;;
;; Caveats
;;    Need additional files containing the area of each grid cell ("areacello")
;;
;; Modification history
;;    20150708-A_senf_da: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "NSIDC"

    ;; Tier
    TIER = 2

    ;; Input dir/files (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"
    INFILES = (/"NSIDC_all_monthly_NH.nc", "NSIDC_all_monthly_SH.nc"/)

    ;; Selected variable (standard name)
    VAR     = "sic"

    ;; Initialize global variable
    FIELD   = "T2Ms"

    ;; Output directory (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"
    outstem_end = (/"_sat_NH_" + FIELD + "_" + VAR, \
                    "_sat_SH_" + FIELD + "_" + VAR /)
                  ; must correspond to INFILES

    ;; Additional settings
    areapath = INDIR
    areafiles = (/"psn25area_v3.dat", "pss25area_v3.dat"/)

end

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_NSIDC.ncl"
    enter_msg(diag_script, "", 1)

    setnames = (/"NSIDC-NT", "NSIDC-BT"/)
    setabbr  = (/"nt", "bt"/)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Loop over files
    do ifile = 0, dimsizes(INFILES) - 1

        ;; Loop over Datasets
        do iset = 0, dimsizes(setnames) - 1

            ;; Read data
            f = addfile(INDIR + INFILES(ifile), "r")
            var = "goddard_" + setabbr(iset) + "_seaice_conc_monthly"
            sic = f->$var$

            ;; Convert from short to float & scale
            sic_tmp = tofloat(sic) / 100.
            copy_VarMeta(sic, sic_tmp)
            sic := sic_tmp
            delete(sic_tmp)

            ;; Reset missing values
            sic = where(sic.eq.sic@_FillValue, 1e+20, sic)
            sic@_FillValue    = 1e+20
            sic@missing_value = sic@_FillValue
            delete(sic@scale_factor)

            ;; Determine time range for output file name
            yyyymm       = tostring(cd_calendar(sic&time, -1))
            yyyymm_start = yyyymm(0)
            yyyymm_end   = yyyymm(dimsizes(yyyymm) - 1)
            start_year   = toint(yyyymm_start) / 100
            end_year     = toint(yyyymm_end) / 100

            ;; Add Coordinate variables
            s = dimsizes(sic)
            j = ispan(0, s(1) - 1, 1)
            i = ispan(0, s(2) - 1, 1)
            sic!1 = "j"
            sic!2 = "i"
            sic&j = j
            sic&i = i
            sic&j@units     = "1"
            sic&j@long_name = "cell index along first dimension"
            sic&i@units     = "1"
            sic&i@long_name = "cell index along second dimension"
            lat   = f->latitude
            lat!0 = "j"
            lat!1 = "i"
            lat&j = j
            lat&i = i
            lon   = f->longitude
            lon!0 = "j"
            lon!1 = "i"
            lon&j = j
            lon&i = i
            sic@coordinates = "lon lat"

            ;; Set global attributes
            gAtt = True
            gAtt@history = "Created on " + systemfunc("date")
            gAtt@host    = systemfunc("echo $HOST")
            gAtt@user    = systemfunc("echo $USER")
            gAtt@period  = start_year + "-" + end_year
            gAtt@field   = FIELD
            gAtt@tier    = TIER
            gAtt@source  = "ftp://sidads.colorado.edu/pub/" + \
                "DATASETS/NOAA/G02202_v2/"
            gAtt@reference = "Not available"
            gAtt@title   = setnames(iset) + \
                " satellite data reformatted for the ESMValTool"
            gAtt@conventions = "CF/CMOR"

            ;; Add areacello
            temp = fbindirread(areapath + areafiles(ifile), 0, (/s(1), s(2)/),\
                               "integer")
            areacello = 1000. * int2flt(temp)  ; m2
            copy_VarMeta(lat, areacello)
            areacello@standard_name = "cell_area"
            areacello@long_name     = "Ocean Grid-Cell Area"
            areacello@units         = "m2"
            areacello@missing_value = sic@_FillValue
            areacello@_FillValue    = sic@_FillValue
            gAtt@areacello          = "added areacello to file, based on " + \
                                      areafiles(ifile)
            delete([/temp/])

            ;; Write netCDF file
            outfile = OUTDIR + "OBS_" + setnames(iset) + outstem_end(ifile) \
                + "_" + yyyymm_start + "-" + yyyymm_end + ".nc"
            write_nc(outfile, VAR, sic, gAtt)
            fout = addfile(outfile, "w")
            fout->lat = lat
            fout->lon = lon
            fout->areacello = areacello

            ;; Clean up
            delete([/f, sic, yyyymm, s, j, i, lat, lon, gAtt, areacello, \
                    outfile, fout/])

        end do  ; iset

    end do  ; ifile

    leave_msg(diag_script, "", 1)

end
