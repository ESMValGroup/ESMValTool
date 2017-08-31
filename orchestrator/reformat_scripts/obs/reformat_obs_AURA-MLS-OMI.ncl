;;#############################################################################
;; REFORMAT SCRIPT FOR AURA-MLS-OMI DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    http://acdb-ext.gsfc.nasa.gov/Data_services/cloud_slice/new_data.html
;;
;; Last access
;;    20140415
;;
;; Download and processing instructions
;;
;; Caveats
;;
;; Modification history
;;    20140417-A_gott_kl wrote data reading algorithmn; adjusted cmor.
;;    20140415-A_fran_fr written based on reformat_obs_AURA-TES.ncl.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "AURA-MLS-OMI"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 2005
    YEAR2 = 2013

    ;; Selected variable (standard name)
    VARS = (/"tropoz"/)

    ;; Name in the raw data
    NAME = (/"O3"/)

    ;; Units in the input  data and corresponding conversion factor
    ;; to CMOR standard units
    RAWUNITS = (/"(Dobson Units) X 10"/)
    CONVERSION = (/0.1/)

    ;; Field
    FIELD = "T2Ms"

    ;; Gridtype
    GRIDTYPE = "Limb"

    ;; # of Header lines
    nhead = 3

    ;; missing values
    FillValue = -999  ; just something impossible
    threshold = (/99.9/)  ; set >= values to missing

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_AURA-MLS-OMI.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Creat time coordinate
    timec = create_timec(YEAR1, YEAR2)
    datec = ut_calendar(timec, 0)

    years = str_get_cols(tostring(ispan(YEAR1, YEAR2, 1)), 2, 3)
    months = (/"jan", "feb", "mar", "apr", "may", "jun", \
               "jul", "aug", "sep", "oct", "nov", "dec"/)

    ;; Grid
    nlon = 288
    nlat = 120
    lon = fspan(-179.375, 179.375, nlon)
    lat = fspan(-59.5, 59.5, nlat)

    ;; Loop over variables
    do vID = 0, dimsizes(VARS) - 1

        info_output("Processing " + VARS(vID), verbosity, 1)

        ;; Loop over time
        tt = 0
        do ii = 0, dimsizes(years) - 1
            do jj = 0, dimsizes(months) - 1

                ;; Find files
                fname = "L3_tropo_ozone_column_" + months(jj) + years(ii)
                files = systemfunc("ls " + INDIR + fname)

                if (.not.all(ismissing(files))) then

                    ;; collect file names -> global attribute in outfile
                    if(tt.eq.0) then
                        original_files = fname
                    else
                        original_files = original_files + ", " + fname
                    end if

                    ;; Read variable
                    data  = new((/nlon, nlat/), float)

                    info_output("Reading " + files, verbosity, 1)
                    head = readAsciiHead(files, nhead)
                    stuff = readAsciiTable(files, 1, "string", 3)
                    stuff = str_get_cols(stuff, 1, 75)

                    split13 = new(14, "integer")
                    split13(:) = 3
                    split13(13) = 14
                    split25 = new(25, "integer")
                    split25(:) = 3

                    i = -1
                    do ilat = 1, nlat
                        do row = 1, 11

                            i = i + 1
                            dat = str_split_by_length(stuff(i, 0), split25)

                            do ilon = 1, 25
                                data(25 * (row - 1) + ilon - 1, ilat - 1)\
                                    = tofloat(dat(ilon - 1))
                            end do

                        end do

                        i = i + 1
                        dat0 = str_split_by_length(stuff(i, 0), split13)

                        do ilon = 1, 13
                            data(275 + ilon - 1, ilat - 1) = \
                                tofloat(dat0(ilon - 1))
                        end do

                    end do
                    delete(stuff)

                    ;; Create array
                    if (.not.isdefined("outvar")) then
                        outvar = new(array_append_record(dimsizes(timec), \
                                                         dimsizes(data), 0),\
                                     float)
                    end if

                    outvar(tt, :, :) = data

                    tt = tt + 1
                    delete(data)

                else

                    info_output(" --> File missing, skipping this date", \
                                verbosity, 1)

                end if
            end do  ; months
        end do  ; years

        ;; Scale & filter missing values
        outvar = outvar * CONVERSION(vID)
        outvar@_FillValue = FillValue
        outvar = where(outvar.eq.threshold, FillValue, outvar)

        ;; Rename coordinate
        outvar!0 = "time"
        outvar!1 = "lon"
        outvar!2 = "lat"
        outvar&lon = lon
        outvar&lat = lat

        ;; Reorder
        output = outvar(time|:, lat|:, lon|:)
        output@_FillValue = FILL
        output&time = timec
        delete(outvar)

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
        clat = lat_attrib(output&lat)
        delete(output&lat)
        output&lat = clat
        delete(clat)

        ;; Format longitude coordinate
        clon = lon_attrib(output&lon)
        delete(output&lon)
        output&lon = clon
        delete(clon)

        ;; Set global attributes
        gAtt = True
        gAtt@comment       = "Values = 99.9 are flagged as missing"
        gAtt@history       = "Created on " + systemfunc("date")
        gAtt@host          = systemfunc("echo $HOST")
        gAtt@user          = systemfunc("echo $USER")
        gAtt@period        = YEAR1 + "-" + YEAR2
        gAtt@field         = FIELD
        gAtt@tier          = TIER
        gAtt@source        = "http://acdb-ext.gsfc.nasa.gov/Data_services/" + \
            "cloud_slice/new_data.html"
        gAtt@reference     = "Ziemke, J. R. et al., Atmos. Chem. Phys., " + \
            "11, 9237-9251, doi:10.5194/acp- 11-9237-2011, 2011."
        gAtt@title = OBSNAME + " satellite data reformatted for the ESMValTool"
        gAtt@conventions   = "CF/CMOR"

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
