;;#############################################################################
;; REFORMAT SCRIPT FOR HadISST REANALYSIS DATA
;;#############################################################################
;;
;; Tier
;;    2 (freely available data set other than obs4MIPs and ana4MIPs)
;;
;; Source
;;    Reference:
;;       Rayner, N. A., D. E. Parker, E. B. Horton, C. K. Folland, L. V.
;;          Alexander, D. P. Rowell, E. C. Kent, A. Kaplan(2003), Global
;;          analyses of sea surface temperature, sea ice, and night marine
;;          air temperature since the late nineteenth century, J. Geophys.
;;          Res., 108(D14), 4407, doi: 10.1029/2002JD002670.
;;
;; Last access
;;    04/2015
;;
;; Download and processing instructions
;;    Download: http://www.metoffice.gov.uk/hadobs/hadisst/data/download.html
;;    Processing: this script (reformat_obs_HadISST.ncl)
;;
;; Caveats
;;    None
;;
;; Modification history
;;    20150422D-A_laue_ax: written.
;;
;; ############################################################################

begin

    ;; Source name
    OBSNAME = "HadISST"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Field
    FIELD = "T2Ms"

end

load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity =  stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_ESACCI-HadISST.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    initialize = 1

    fname_sst = INDIR + "HadISST_sst.nc"
    fname_ice = INDIR + "HadISST_ice.nc"

    info_output("*** " + fname_sst + " ***", verbosity, 1)
    info_output("*** " + fname_ice + " ***", verbosity, 1)

    fice = addfile(fname_ice, "r")
    fsst = addfile(fname_sst, "r")

    lat  = (/fsst->latitude/)
    lon  = (/fsst->longitude/)
    time = (/fsst->time/)

    lat@standard_name = "latitude"
    lat@units = "degrees_north"
    lat@axis = "Y"
    lat@long_name = "latitude"
    lat!0 = "lat"

    lon@standard_name = "longitude"
    lon@units = "degrees_east"
    lon@axis = "X"
    lon@long_name = "longitude"
    lon!0 = "lon"

    time@standard_name = "time"
    time@units = "days since 1870-01-01 00:00:00"
    time@axis = "T"
    time@long_name = "time"
    time@calendar = "gregorian"
    if (isatt(time, "_FillValue")) then
        delete(time@_FillValue)
    end if

    sst = (/fsst->sst/)

    sst@_FillValue = -1000
    sst@missing_value = -1000
    sst@units = "K"
    sst@long_name = "sea surface temperature"
    sst@standard_name = "sst"

    sst!0 = "time"
    sst!1 = "lat"
    sst!2 = "lon"
    sst&time = time
    sst&lat = lat
    sst&lon = lon

    ; convert from degrees C to K
    sst = sst + 273.15

    gAtt = True
    gAtt@history     = "Created on " + systemfunc("date")
    gAtt@host        = systemfunc("echo $HOST")
    gAtt@user        = systemfunc("echo $USER")
    gAtt@period      = "1870-2014"
    gAtt@tier        = TIER
    gAtt@field       = FIELD
    gAtt@source      = "http://www.metoffice.gov.uk/hadobs/hadisst/data/" + \
        "download.html"
    gAtt@reference   = "Rayner et al. (2013), J. Geophys. Res."
    gAtt@title       = "HadISST reanalysis data reformatted" + \
        " for the ESMValTool"
    gAtt@conventions = "CF/CMOR"

    ; outfile (file naming following Dominik Kunert, DLR)
    fout = OUTDIR + "OBS_HadISST_reanaly_1_" + FIELD + "_ts_187001-201412.nc"

    ; write sst
    write_nc(fout, "ts", sst, gAtt)

    delete(sst)

    ice = (/fice->sic/)

    ice@_FillValue = -1.0e30
    ice@missing_value = -1.0e30
    ice@units = "1"
    ice@long_name = "monthly 1 degree resolution sea ice concentration"
    ice@standard_name = "sea_ice_area_fraction"

    ice!0 = "time"
    ice!1 = "lat"
    ice!2 = "lon"
    ice&time = time
    ice&lat = lat
    ice&lon = lon

    gAtt = True
    gAtt@history     = "Created on " + systemfunc("date")
    gAtt@host        = systemfunc("echo $HOST")
    gAtt@user        = systemfunc("echo $USER")
    gAtt@period      = "1870-2014"
    gAtt@tier        = TIER
    gAtt@field       = FIELD
    gAtt@source      = "http://www.metoffice.gov.uk/hadobs/hadisst/data/" + \
        "download.html"
    gAtt@reference   = "Rayner et al. (2013), J. Geophys. Res."
    gAtt@title       = "HadISST reanalysis data reformatted" + \
        " for the ESMValTool"
    gAtt@conventions = "CF/CMOR"

    ; outfile (file naming following Dominik Kunert, DLR)
    fout = OUTDIR + "OBS_HadISST_reanaly_20130524_" + FIELD + \
        "_sic_187001-201412.nc"

    ; write sic
    write_nc(fout, "sic", ice, gAtt)

    delete(ice)

    leave_msg(diag_script, "", 1)

end
