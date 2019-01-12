; #############################################################################
; ESMValTool CMORizer for ESACCI-LANDCOVER data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    ftp://anon-ftp.ceda.ac.uk/neodc/esacci/land_cover/data/land_cover_maps/
;
; Last access
;    20190110
;
; Download and processing instructions
;    Download the 3 NetCDF files for 2000, 2005 and 2010.
;    Download the CCI-LC Tools from:
;      http://maps.elie.ucl.ac.be/CCI/viewer/download/lc-user-tools-3.14.zip
;    Unpack and run the CCI-LC Tools on each of the NetCDF files as follows:
;      bash lc-user-tools/bin/aggregate-map.sh -PgridName=GEOGRAPHIC_LAT_LON \
;      -PnumMajorityClasses=1 -PoutputAccuracy=false -PoutputPFTClasses=true \
;      -PoutputLCCSClasses=false <inputfile.nc>
;    Put the resulting processed data in the input_dir_path.
;
; Caveat
;   The CCI-LC Tools must be applied before running this script.
;   The CCI-LC Tools require Java Version 7 or higher.
;   The input data are available for a single year and are copied over to
;   generate a time series over their time range of validity.
;
; Modification history
;    20190110-A_righ_ma: rewritten in NCL for v2.
;    20160714-A_muel_bn: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + "/cmor/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_ESACCI-LANDCOVER.ncl"

  ; Source name
  OBSNAME = "ESACCI-LANDCOVER"

  ; Tier
  TIER = 2

  ; Years
  YEARS = (/2000, 2005, 2010/)

  ; Version
  VERSION = "L4-LCCS-Map-300m-P5Y-aggregated-0.083333Deg"

  ; Field
  FIELD = "T2Ms"

  ; CMOR table
  CMOR_TABLE = \
    getenv("esmvaltool_root") + "/cmor/tables/cmip5/Tables/CMIP5_Lmon"

  ; Variable names
  VARIABLES = (/"baresoilFrac"/)

  ; Auxiliary coordinate
  TYPES = (/"bare_ground"/)

  ; Corresponding aggregation classes in the raw data
  CLASSES = [/"Bare_Soil"/]

end

begin

  do yy = 0, dimsizes(YEARS) - 1

    fname = \
      input_dir_path + "/ESACCI-LC-" + VERSION + "-" + YEARS(yy) + "-v1.6.1.nc"

    f = addfile(fname, "r")

    ; Create time coordinate
    YEAR1 = YEARS(yy) - 2
    YEAR2 = YEARS(yy) + 2
    time = create_timec(YEAR1, YEAR2)

    do vv = 0, dimsizes(VARIABLES) - 1

      log_info("Processing variable " + VARIABLES(vv))

      ; Sum up classes
      class = CLASSES[vv]
      do cc = 0, dimsizes(class) - 1
        log_info("  adding class " + class(cc))
        if (cc.eq.0) then
          xx = f->$class(cc)$
        else
          xx = xx + f->$class(cc)$
        end if
      end do
      delete(class)

      ; Replace NaNs with missing
      replace_ieeenan(xx, FILL, 0)
      xx@_FillValue = FILL

      ; Define output field
      output = \
        new((/dimsizes(time), dimsizes(xx&lat), dimsizes(xx&lon)/), float)
      output!0 = "time"
      output&time = time
      output!1 = "lat"
      output&lat = xx&lat
      output!2 = "lon"
      output&lon = xx&lon
      output = conform(output, xx, (/1, 2/))
      delete(xx)

      ; Set standard fill value
      output@_FillValue = FILL

      ; Convert units
      output = output * 100
      output@units = "%"

      ; Set variable attributes
      tmp = var_attrib(output, VARIABLES(vv), CMOR_TABLE)
      delete(output)
      output = tmp
      delete(tmp)

      ; Format time coordinate
      ctime = time_attrib(output&time, YEAR1, YEAR2, "M")
      delete(output&time)
      output&time = ctime
      delete(ctime)

      ; Format latitude coordinate
      output = output(:, ::-1, :)
      clat = lat_attrib(output&lat)
      delete(output&lat)
      output&lat = clat
      delete(clat)

      ; Format longitude coordinate
      output = lonFlip(output)
      clon = lon_attrib(output&lon)
      delete(output&lon)
      output&lon = clon
      delete(clon)

      ; Add auxiliary coordinate
      output@coordinates = "type"

      ; Set global attributes
      gAtt = True
      gAtt@history       = "Created on " + systemfunc("date")
      gAtt@host          = systemfunc("echo $HOST")
      gAtt@user          = systemfunc("echo $USER")
      gAtt@period        = YEAR1 + "-" + YEAR2
      gAtt@field         = FIELD
      gAtt@tier          = TIER
      gAtt@source        = "ftp://anon-ftp.ceda.ac.uk/neodc/esacci/" + \
        "land_cover/data/land_cover_maps/"
      gAtt@reference     = "Defourny, P., et al., 2015. ESA Land Cover " + \
        "Climate Change Initiative (ESA LC_cci) data"
      gAtt@title = OBSNAME + " satellite data CMORized for the ESMValTool v2.0"
      gAtt@conventions = "CF/CMOR"

      ; Outfile
      fout = output_dir_path + "/" + \
        "OBS_" + OBSNAME + "_sat_" + VERSION + "_" + FIELD + "_" + \
        VARIABLES(vv) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

      ; Write variable
      write_nc(fout, VARIABLES(vv), output, gAtt)
      delete(gAtt)
      delete(output)

      ; Add auxiliary coordinate
      type = tochar(TYPES(vv))
      type!0 = "strlen"
      type@long_name = "surface type"
      type@standard_name = "area_type"
      w = addfile(fout, "w")
      w->type = type
      delete(w)

    end do
  end do

end
