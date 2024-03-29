; #############################################################################
; ESMValTool CMORizer for ESACCI-CLOUD data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    https://public.satproj.klima.dwd.de/data/ESA_Cloud_CCI/CLD_PRODUCTS/v3.0/
;
; Last access
;    20210428
;
; Download and processing instructions
;    Download the data from:
;      L3C/AVHRR-AM/ and L3C/AVHRR-PM, e.g.:
;      wget -r -nH -e robots=off --cut-dirs=9 --no-parent
;        --reject="index.html*"
;        https://public.satproj.klima.dwd.de/data/ESA_Cloud_CCI/
;          CLD_PRODUCTS/v3.0/L3C/AVHRR-AM/
;      wget -r -nH -e robots=off --cut-dirs=9 --no-parent
;        --reject="index.html*"
;        https://public.satproj.klima.dwd.de/data/ESA_Cloud_CCI/
;          CLD_PRODUCTS/v3.0/L3C/AVHRR-PM/
;
;    All files are expected in a single directory (no subdirectories
;    with years).
;
; Modification history
;    20230818-lauer_axel: added output of clwvi (in addition to iwp, lwp)
;    20210428-lauer_axel: AVHRR-AM and AVHRR-PM data are now averaged during
;                         the overlapping time; TOA radiative fluxes are now
;                         also processed
;    20190201-righi_mattia: adapted to v2.
;    20181116-lauer_axel: changed to use CLOUD-CCI v3.0 data (AVHRR-PM), gaps
;                         (1994/09 - 1995/01) are filled with AVHRR-AM data
;    20180522-righi_mattia: changed to use AVHRR-PM data.
;    20160519-laue_axel: written (reformat_obs_ESACCI-AEROSOL.ncl).
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/data/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "esacci_cloud.ncl"

  ; Source name
  OBSNAME = "ESACCI-CLOUD"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = get_year(start_year, 1982)
  YEAR2 = get_year(end_year, 2016)

  ; Selected variable (standard name)
  VAR = (/"clt", "cltStderr", "clivi", "lwp", "clwvi", "rlut", "rlutcs", \
          "rsut", "rsutcs", "rsdt", "rlus", "rsus", "rsuscs"/)

  ; Name in the raw data
  NAME = (/"cfc", "cfc_unc", "iwp_allsky", "lwp_allsky", "iwp_allsky", \
           "toa_lwup", "toa_lwup_clr", "toa_swup", "toa_swup_clr", \
           "toa_swdn", "boa_lwup", "boa_swup", "boa_swup_clr"/)

  ; Conversion factor
  CONV = (/100., 1., 0.001, 0.001, 0.001, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, \
          1.0, 1.0/)

  ; MIP
  MIP = (/"Amon", "Amon", "Amon", "Amon", "Amon", "Amon", "Amon", "Amon", \
          "Amon", "Amon", "Amon", "Amon", "Amon"/)

  ; Frequency
  FREQ = (/"mon", "mon", "mon", "mon", "mon", "mon", "mon", "mon", "mon", \
           "mon", "mon", "mon", "mon"/)

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + \
    (/"/cmip5/Tables/CMIP5_Amon", \
      "/custom/CMOR_cltStderr.dat", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/custom/CMOR_lwp.dat", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon"/)

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "AVHRR-AMPM-fv3.0"

  ; Global attributes
  SOURCE = "https://public.satproj.klima.dwd.de/data/ESA_Cloud_CCI/" + \
    "CLD_PRODUCTS/v3.0/"
  REF = \
    "Stengel et al., Earth Syst. Sci. Data, doi:10.5194/essd-9-881-2017, 2017"
  COMMENT = ""

end

begin

  firstime = True

  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    time = create_timec(YEAR1, YEAR2)
    date = cd_calendar(time, 1)

    ; Create timeseries
    do yy = YEAR1, YEAR2

      syear = sprinti("%i", yy)
      do mm = 1, 12

        smonth = sprinti("%0.2i", mm)

        ; Read file
        fname = systemfunc("ls " + input_dir_path + syear + smonth + \
                           "-ESACCI-L3C_CLOUD-CLD_PRODUCTS-AVHRR_*-" + \
                           str_sub_str(VERSION, "AVHRR-AMPM-", "") + ".nc")

        ; No files found
        if (all(ismissing(fname))) then
          continue
        end if

        ifiles = dimsizes(fname)
        log_info(syear + smonth + ": " + tostring(ifiles) + " input file(s)")

        do i = 0, ifiles - 1
          ; Extract data
          f = addfile(fname(i), "r")
          xx = f->$NAME(vv)$

          if (i .eq. 0) then
            xdims = array_append_record((/ifiles/), dimsizes(xx), 0)
            xx_all = new(xdims, float)
            delete(xdims)
          end if

          ; Convert units
          xx_all(i, :, :, :) = xx * CONV(vv)

          ; *** calculate clwvi (lwp + iwp) ***
          if (VAR(vv) .eq. "clwvi") then
            xx2 = f->lwp_allsky * 0.001
            xx_all(i, :, :, :) = xx_all(i, :, :, :) + xx2
            delete(xx2)
          end if

          if (firstime) then
            lat = f->lat
            lon = f->lon
            firstime = False
          end if
          delete(f)
        end do  ; loop over all files (am/pm) per date

        delete(fname)
        delete(xx)

        xx = dim_avg_n(xx_all, 0)  ; ignore missing values
        delete(xx_all)

        ; Assign to global array
        if (.not.isdefined("output")) then
          dims = dimsizes(xx)
          dims(0) = dimsizes(time)
          output = new(dims, float)
          output!0 = "time"
          output&time = time
          output!1 = "lat"
          output&lat = lat
          output!2 = "lon"
          output&lon = lon
        end if
        output(ind(toint(yy * 100 + mm).eq.date), :, :) = (/xx/)

      end do
    end do

    ; Set fill value
    output = where(output.eq.-999, output@_FillValue, output)

    ; Format coordinates
    output!0 = "time"
    output!1 = "lat"
    output!2 = "lon"
    format_coords(output, YEAR1 + "0101", YEAR2 + "1231", FREQ(vv))

    ; Set variable attributes
    tmp = format_variable(output, VAR(vv), CMOR_TABLE(vv))
    delete(output)
    output = tmp
    delete(tmp)

    ; Calculate coordinate bounds
    bounds = guess_coord_bounds(output, FREQ(vv))

    ; Set global attributes
    gAtt = set_global_atts(OBSNAME, TIER, SOURCE, REF, COMMENT)

    ; Output file
    DATESTR = YEAR1 + "01-" + YEAR2 + "12"
    fout = output_dir_path + \
      str_join((/"OBS", OBSNAME, TYPE, VERSION, \
                 MIP(vv), VAR(vv), DATESTR/), "_") + ".nc"

    ; Write variable
    write_nc(fout, VAR(vv), output, bounds, gAtt)
    delete(gAtt)
    delete(output)
    delete(bounds)

  end do

end
