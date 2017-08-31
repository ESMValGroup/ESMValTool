;;#############################################################################
;; REFORMAT SCRIPT FOR Tilmes OZONE SONDES DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    http://acd.ucar.edu/~tilmes/ozone.html
;;    See also:
;;    Tilmes, S., J.-F. Lamarque, L. K. Emmons, A. Conley, M. G. Schultz,
;;    M. Saunois, V. Thouret, A. M. Thompson, S. J. Oltmans, B. Johnson, and
;;    D. Tarasick (2011): Technical Note: Ozonesonde climatology between 1995
;;    and 2011: Description, Evaluation and Applications, Atmos. Chem. Phys.,
;;    12, 7475-7497, 2012, doi:10.5194/acp-12-7475-2012
;;
;; Last access
;;    20150720
;;
;; Download and processing instructions
;;    Select Monthly averaged ozone profiles averaged between 1995-2011, and
;;    Profiles in pressure altitudes for 42 sites (ncdf).
;;
;; Caveats
;;    Raw data contain 12 monthly values per file, averaged over 1995-2011.
;;    The time coordinate in the output files is arbitrarily set to the 15th
;;    of each month in 2000.
;;
;; Modification history
;;    20150720-A_righ_ma: updated to new data version (1995-2011).
;;    20140630-A_gott_kl: written
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "Tilmes"

    ;; Time range (data version)
    RANGE = "1995_2011"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/" + RANGE + "/stations/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Input file selectors:
    ;; Specify explicitely, since not all available files contribute to
    ;; regional averages. Compare to table 2 of Tilmes et al., 2012
    ;;  File, lat, lon(E)
    INFILES = (/"Alert" + RANGE + ".nc", "83", "-62", \
                "Ascension" + RANGE + ".nc", "-8", "-14", \
                "Churchill" + RANGE + ".nc", "59", "-94", \
                "Debilt" + RANGE + ".nc", "52", "5", \
                "Edmonton" + RANGE + ".nc", "53", "-114", \
                "Eureka" + RANGE + ".nc", "80", "-86", \
                "Fiji" + RANGE + ".nc", "-18", "178", \
                "Goose_bay" + RANGE + ".nc", "53", "-60", \
                "Hilo" + RANGE + ".nc", "19", "-155", \
                "Hohenpeissenberg" + RANGE + ".nc", "48", "11", \
                "Hongkong" + RANGE + ".nc", "22", "114", \
                "Huntsville" + RANGE + ".nc", "35", "87", \
                "Kagoshima" + RANGE + ".nc", "32", "131", \
                "Lauder" + RANGE + ".nc", "-45", "170", \
                "Legionowo" + RANGE + ".nc", "52", "21", \
                "Lerwick" + RANGE + ".nc", "60", "-1", \
                "Lindenberg" + RANGE + ".nc", "52", "14", \
                "Macquarie" + RANGE + ".nc", "-55", "159", \
                "Marambio" + RANGE + ".nc", "-64", "57", \
                "Naha" + RANGE + ".nc", "26", "123", \
                "Nairobi" + RANGE + ".nc", "-1", "37", \
                "Natal" + RANGE + ".nc", "-5", "-35", \
                "Neumayer" + RANGE + ".nc", "-71", "-8", \
                "Ny_Alesund" + RANGE + ".nc", "79", "12", \
                "Paramaribo" + RANGE + ".nc", "6", "-55", \
                "Payerne" + RANGE + ".nc", "47", "7", \
                "Praha" + RANGE + ".nc", "50", "15", \
                "Resolute" + RANGE + ".nc", "74", "-95", \
                "Samoa" + RANGE + ".nc", "-14", "-171", \
                "Sancristobal" + RANGE + ".nc", "-1", "-90", \
                "Scoresbysund" + RANGE + ".nc", "71", "-22", \
                "Syowa" + RANGE + ".nc", "-69", "40", \
                "Tateno" + RANGE + ".nc", "36", "140", \
                "Uccle" + RANGE + ".nc", "51", "4", \
                "Wallops_Island" + RANGE + ".nc", "40", "-76", \
                "Watukosek" + RANGE + ".nc", "-8", "113"/)

    ;; Regional averages
    LATMAX = (/20., 90., -20./)  ; (/90., 30., 0., -30./)
    LATMIN = (/-20., 20., -90./)  ; (/30., 0., -30., -90./)
    REGIONS = (/"Trop", "NHext", "SHext"/)  ; (/"NE", "NT", "ST", "SE"/)

    ;; Variable name in the raw files
    GAS = (/"o3_mean"/)

    YEAR1 = 2000  ; 1995
    YEAR2 = 2000  ; 2009

    ;; Selected variable (standard name): one per corresponding element of GAS
    VAR = (/"tro3"/)

    ;; Units conversion
    UNITS = (/"1e9"/)  ; target units
    CONVERSION = (/1./)  ; conversion factor to be multiplied to raw data

    FIELD = "T1M"

    soft_links = True

end

load "./interface_scripts/messaging.ncl"
load "./interface_scripts/auxiliary.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_Tilmes.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Input ##################################################################
    nfiles = (dimsizes(INFILES)) / 3
    in_files = new(nfiles, string)
    lat = new(nfiles, float)
    lon = new(nfiles, float)
    do i = 0, nfiles - 1
        in_files(i) = INFILES(i * 3)
        lat(i) = tofloat(INFILES((i * 3) + 1))
        lon(i) = tofloat(INFILES((i * 3) + 2))
    end do
    lon = where(lon.lt.0., 360. + lon, lon)

    ;; Output files ###########################################################
    do i = 0, dimsizes(REGIONS) - 1
        do j = 0, dimsizes(VAR) - 1

            spec = GAS(j)
            out_file = OUTDIR + "OBS_" + OBSNAME + "_insitu_" + \
                REGIONS(i) + "_" + FIELD + "_" + VAR(j) + "_" + \
                str_get_field(RANGE, 1, "_") + "01-" + \
                str_get_field(RANGE, 2, "_") + "12.nc"

            ;; Select region
            index = ind(lat.le.LATMAX(i) .and. lat.ge.LATMIN(i))
            if (all(ismissing(index))) then
                error_msg("f", diag_script, "", "no station found in " + \
                          "selected region (" + tostring(LATMAX(i)) + "-" + \
                          tostring(LATMIN(i)))
            end if

            ;; Loop over all files in that region
            rfiles = INDIR + in_files(index)
            nrfiles = dimsizes(index)
            rlat = lat(index)
            rlon = lon(index)

            ;; dimensions
            f = addfile(rfiles(0), "r")
            uu = f->$spec$(months|:, level|:)
            pp = f->press(level|:)
            np = dimsizes(pp)
            month = uu&months
            nt = dimsizes(month)
            data_full = new((/nt, np, nrfiles/), float)
            data_full(:, :, 0) = uu
            delete(f)
            delete(uu)

            ;; collect data
            if(nrfiles.gt.1) then
                do k = 1, nrfiles - 1
                    info_output("Reading " + rfiles(k), verbosity, 1)
                    f = addfile(rfiles(k), "r")
                    uu = f->$spec$(months|:, level|:)
                    data_full(:, :, k) = uu
                    delete(f)
                    delete(uu)
                end do
            end if

            ;; average & scale
            data = dim_avg(data_full)
            delete(data_full)
            data = data * CONVERSION(j)

            ;; time coordinates
            timec = create_timec(YEAR1, YEAR2)
            datec = ut_calendar(timec, 1)
            data!0 = "time"
            data&time = timec
            delete(timec)

            ;; pressure coordinates
            data!1 = "plev"
            if (isMonotonic(pp).eq.1) then
                plev = 100. * pp(::-1)
            else
                plev = 100. * pp
            end if
            plev@long_name     = "pressure"
            plev@positive      = "down"
            plev@axis          = "Z"
            plev@units         = "Pa"
            plev@standard_name = "air_pressure"
            data&plev = plev
            delete(plev)

            ;; Write data
            output = var_attrib(data, VAR(j))

            ;; Reformat time coordinate
            ctime = time_attrib(output&time)
            delete(output&time)
            output&time = ctime
            delete(ctime)

            ;; Format lat & lon
            slon = str_concat(str_squeeze(sprintf("%4.0f", rlon)) + ", ")
            slat = str_concat(str_squeeze(sprintf("%4.0f", rlat)) + ", ")
            slon = "(/" + str_get_cols(slon, 0, -3) + "/)"
            slat = "(/" + str_get_cols(slat, 0, -3) + "/)"

            ;; Set global attributes
            gAtt = True
            gAtt@history = "Created on " + systemfunc("date")
            gAtt@host = systemfunc("echo $HOST")
            gAtt@user = systemfunc("echo $USER")
            gAtt@longitudes = slon
            gAtt@latitudes = slat
            gAtt@raw_directory = INDIR
            gAtt@raw_files = str_concat(in_files(index) + ", ")
            gAtt@latmax = LATMAX(i)
            gAtt@latmin = LATMIN(i)
            gAtt@period = str_sub_str(RANGE, "_", "-")
            gAtt@field = FIELD
            gAtt@tier = TIER
            gAtt@source = "http://acd.ucar.edu/~tilmes/ozone.html"
            gAtt@reference = "Tilmes, S. et al., Atmos. Chem. Phys., " + \
                "12, 7475-7497, doi:10.5194/acp-12-7475-2012, 2012."
            gAtt@title = OBSNAME + " in-situ data " + \
                "reformatted for the ESMValTool"
            gAtt@conventions = "CF/CMOR"

            ;; Write variable
            write_nc(out_file, VAR(j), output, gAtt)
            info_output("  Wrote " + out_file, verbosity, 1)
            delete(gAtt)
            delete(output)
            delete(data)

            ;; Create soft links (workaround to deal with naming conventions)
            if (soft_links) then
                link_file = str_sub_str(out_file, VAR, VAR + "_" + REGIONS(i))
                system("ln -s " + out_file + " " + link_file)
            end if

        end do  ; j (variables)

        delete(index)
        delete(rlat)
        delete(rlon)
        delete(rfiles)
    end do  ; i (regions)

    leave_msg(diag_script, "", 1)

end
