;;#############################################################################
;; REFORMAT SCRIPT FOR DONG ET AL (2008) ARGO ASCII DATA
;;#############################################################################
;;
;; Source
;;  http://onlinelibrary.wiley.com/doi/10.1029/2006JC004051/full
;;  (then click on Supporting Information)
;;
;; Last access
;;    20150915
;;
;; Download and processing instructions
;;    Download following file:
;;      jgrc10600-sup-0002-ds01.txt
;;
;;    This script will convert the text file to a netcdf file that is accepted
;;    by the ESMValTool. Only MLDpd will be included.
;;
;; Modification history
;;    20150915-A_laue_ax: updated download information, added global
;;                        attributes TIER and REFERENCE.
;;    20150521-A_vanu_be: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "Dong08-ARGO-monthly"

    ;; Input dir (raw data) and file
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier2/" + OBSNAME + "/"
    fname = INDIR + "jgrc10600-sup-0002-ds01.txt"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier2/" + OBSNAME + "/"

    ;; Selected variable (standard name)
    VAR = "mlotst"

    ;; Name in the raw data
    NAME = "MLDpd"

    ;; Field type
    FIELD = "TO2Ms"

    ;; Tier
    TIER = 2

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_Dong08.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent_esmval(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; set properties of file
    ncol = 15  ; number of columns in file
    nlat = 35
    nlon = 360
    npoints = nlat * nlon
    nmonths = 12
    nhead = 5  ; skip the first lines

    ;; read the input file
    indata = readAsciiTable(fname, ncol, "float", nhead)

    ;; check we the expected number of rows
    nelems = dimsizes(indata)
    if (nelems(0) .ne. nmonths * npoints) then
        error_msg("f", get_script_name(), "", \
                  "Unexpected number of rows in " + INFILE)
    end if

    ;; take lat and lon from first month
    lat = indata(0:nlat - 1, 1)
    lat@units = "degrees_north"
    lat@long_name = "Latitude"
    lat@standard_name = "latitude"
    lon = indata(0:npoints - 1:nlat, 2)
    lon@units = "degrees_east"
    lon@long_name = "Longitude"
    lon@standard_name = "longitude"

    ;; create time variable, use middle of month
    time = indata(0:nmonths * npoints - 1:npoints, 0) - 0.5
    time@units = "months since 0000-1-1"

    ;; create the array from the required data column
    ;; note that lat and lon are swapped still
    mlotst = reshape(indata(:, 3), (/nmonths, nlon, nlat/))
    mlotst@units = "m"
    mlotst@standard_name = "ocean_mixed_layer_thickness"
    mlotst@long_name = "Ocean Mixed Layer Thickness"
    mlotst@_FillValue = 9999.0
    mlotst@comment = "MLDpd from Dong (2008)"

    ;; assign dimensions and coordinate variables
    mlotst!0 = "time"
    mlotst!1 = "lon"
    mlotst!2 = "lat"
    mlotst&time = time
    mlotst&lon = lon
    mlotst&lat = lat

    ;; now we have named dimensions we can re-arrange lat and lon
    ;; go from (time,lon,lat) to (time,lat,lon)
    mlotst_sw = mlotst(time|:, lat|:, lon|:)

    ;; global attributes
    gAtt = True
    gAtt@history       = "Created on " + systemfunc("date")
    gAtt@host          = systemfunc("echo $HOST")
    gAtt@user          = systemfunc("echo $USER")
    gAtt@period        = "climatology"
    gAtt@field         = FIELD
    gAtt@tier          = TIER
;    gAtt@source        = "http://onlinelibrary.wiley.com/doi/10.1029" + \
;                         "/2006JC004051/suppinfo"
    gAtt@source        = "http://onlinelibrary.wiley.com/doi/10.1029/2006JC004051/full"
    gAtt@reference     = "Dong, S. et al., J. Geophys. Res. (oceans), " + \
                         "113, C6, doi:10.1029/2006JC004051, 2008."

    gAtt@inputfile     = fname
    gAtt@title         = OBSNAME + \
        " climatology data for the ESMVal evaluation tool"

    ;; Write variable
    fout = OUTDIR + "OBS_" + OBSNAME + "_insitu_1_" + FIELD + "_" + VAR + \
        "_000001-000012.nc"
    write_nc(fout, VAR, mlotst_sw, gAtt)

    info_output("Done, output in " + fout, verbosity, 1)
    delete([/fout, indata, mlotst, mlotst_sw, time, lat, lon/])

    leave_msg(diag_script, "", 1)

end
