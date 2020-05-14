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
;    20190201
;
; Download and processing instructions
;    Download the data from:
;      L3C/AVHRR-PM/
;    To fill the gap 199409-199501, also download:
;      L3C/AVHRR-AM/AVHRR_NOAA-12/1994/ \
;        199409-ESACCI-L3C_CLOUD-CLD_PRODUCTS-AVHRR_NOAA-12-fv3.0.nc
;      L3C/AVHRR-AM/AVHRR_NOAA-12/1994/ \
;        199410-ESACCI-L3C_CLOUD-CLD_PRODUCTS-AVHRR_NOAA-12-fv3.0.nc
;      L3C/AVHRR-AM/AVHRR_NOAA-12/1994/ \
;        199411-ESACCI-L3C_CLOUD-CLD_PRODUCTS-AVHRR_NOAA-12-fv3.0.nc
;      L3C/AVHRR-AM/AVHRR_NOAA-12/1994/ \
;        199412-ESACCI-L3C_CLOUD-CLD_PRODUCTS-AVHRR_NOAA-12-fv3.0.nc
;      L3C/AVHRR-AM/AVHRR_NOAA-12/1995/ \
;        199501-ESACCI-L3C_CLOUD-CLD_PRODUCTS-AVHRR_NOAA-12-fv3.0.nc
;    Put all files under a single directory (no subdirectories with years).
;
; Modification history
;    20190201-righi_mattia: adapted to v2.
;    20181116-lauer_axel: changed to use CLOUD-CCI v3.0 data (AVHRR-PM), gaps
;                         (1994/09 - 1995/01) are filled with AVHRR-AM data
;    20180522-righi_mattia: changed to use AVHRR-PM data.
;    20160519-laue_axel: written (reformat_obs_ESACCI-AEROSOL.ncl).
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "esacci_cloud.ncl"

  ; Source name
  OBSNAME = "ESACCI-CLOUD"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1982
  YEAR2 = 2016

  ; Selected variable (standard name)
  VAR = (/"clt", "cltStderr", "clivi", "clwvi"/)

  ; Name in the raw data
  NAME = (/"cfc", "cfc_unc", "iwp_allsky", "lwp_allsky"/)

  ; Conversion factor
  CONV = (/100., 1., 0.001, 0.001/)

  ; MIP
  MIP = (/"Amon", "Amon", "Amon", "Amon"/)

  ; Frequency
  FREQ = (/"mon", "mon", "mon", "mon"/)

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + \
    (/"/cmip5/Tables/CMIP5_Amon", \
      "/custom/CMOR_cltStderr.dat", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon"/)

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "AVHRR-fv3.0"

  ; Global attributes
  SOURCE = "https://public.satproj.klima.dwd.de/data/ESA_Cloud_CCI/" + \
    "CLD_PRODUCTS/v3.0/"
  REF = \
    "Stengel et al., Earth Syst. Sci. Data, doi:10.5194/essd-9-881-2017, 2017"
  COMMENT = ""

end

begin

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
                           "-ESACCI-L3C_CLOUD-CLD_PRODUCTS-AVHRR_NOAA-*-" + \
                           str_sub_str(VERSION, "AVHRR-", "") + ".nc")

        ; No files found
        if (all(ismissing(fname))) then
          continue
        end if

        ; Extract data
        f = addfile(fname, "r")
        xx = f->$NAME(vv)$

        ; Convert units
        xx = xx * CONV(vv)

        ; lwp is not a CMOR variable, derive as clwvi = lwp + iwp
        if (VAR(vv).eq."clwvi") then

          ; Read 2nd variable containing iwp (variable "iwp_allsky")
          iwpname = "iwp_allsky"
          xx2 = f->$iwpname$

          ; Convert units
          xx2 = xx2 * 0.001

          ; add iwp to lwp to calculate cmor variable "clwvi"
          xx = xx + xx2
          delete(xx2)
        end if

        ; Assign to global array
        if (.not.isdefined("output")) then
          dims = dimsizes(xx)
          dims(0) = dimsizes(time)
          output = new(dims, float)
          output!0 = "time"
          output&time = time
          output!1 = "lat"
          output&lat = f->lat
          output!2 = "lon"
          output&lon = f->lon
        end if
        output(ind(toint(yy * 100 + mm).eq.date), :, :) = (/xx/)

        delete(fname)
        delete(f)

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
