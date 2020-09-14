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
;      bash lc-user-tools-3.14/bin/aggregate-map.sh \
;        -PgridName=GEOGRAPHIC_LAT_LON -PnumMajorityClasses=1 \
;        -PoutputAccuracy=false -PoutputPFTClasses=true \
;        -PoutputLCCSClasses=false -PnumRows=360 <inputfile.nc>
;    Put the resulting processed data in input_dir_path.
;
; Caveat
;    The CCI-LC Tools must be applied before running this script.
;    The CCI-LC Tools require Java Version 7 or higher.
;    The input data are available for a single year and are copied over to
;    generate a time series over their time range of validity.
;
; Modification history
;    20190110-righi_mattia: rewritten in NCL for v2.
;    20160714-benjamin_mueller: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "esacci_landcover.ncl"

  ; Source name
  OBSNAME = "ESACCI-LANDCOVER"

  ; Tier
  TIER = 2

  ; Years
  YEARS = (/2000, 2005, 2010/)

  ; Variable names
  VAR = \
    (/"baresoilFrac", "cropFrac", "grassFrac", "shrubFrac", "treeFrac"/)

  ; Corresponding aggregation classes in the raw data
  CLASSES = [/"Bare_Soil", \
              "Managed_Grass", \
              "Natural_Grass", \
              (/"Shrub_Broadleaf_Deciduous",    \
                "Shrub_Broadleaf_Evergreen",    \
                "Shrub_Needleleaf_Evergreen"/), \
              (/"Tree_Broadleaf_Deciduous",  \
                "Tree_Broadleaf_Evergreen",  \
                "Tree_Needleleaf_Deciduous", \
                "Tree_Needleleaf_Evergreen"/)/]

  ; MIPs
  MIP = (/"Lmon", "Lmon", "Lmon", "Lmon", "Lmon"/)

  ; Frequency
  FREQ = (/"mon", "mon", "mon", "mon", "mon"/)

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + "/cmip5/Tables/CMIP5_Lmon"

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "L4-LCCS-Map-300m-P5Y-aggregated-0.500000Deg"

  ; Global attributes
  SOURCE = "ftp://anon-ftp.ceda.ac.uk/neodc/esacci/land_cover/data/"
  REF = "Defourny, P.: ESA Land Cover Climate Change Initiative " + \
    "(Land_Cover_cci): Global Land Cover Maps, Version 1.6.1. " + \
    "Centre for Environmental Data Analysis, " + \
    "http://catalogue.ceda.ac.uk/uuid/4761751d7c844e228ec2f5fe11b2e3b0, 2016."
  COMMENT = ""

end

begin

  do yy = 0, dimsizes(YEARS) - 1

    fname = \
      input_dir_path + "ESACCI-LC-" + VERSION + "-" + YEARS(yy) + "-v1.6.1.nc"

    f = addfile(fname, "r")

    ; Create time coordinate
    YEAR1 = YEARS(yy) - 2
    YEAR2 = YEARS(yy) + 2
    time = create_timec(YEAR1, YEAR2)

    do vv = 0, dimsizes(VAR) - 1

      log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

      ; Set classes to be added up
      class = CLASSES[vv]

      ; Save mask before adding up classes
      do cc = 0, dimsizes(class) - 1
        qq = f->$class(cc)$
        replace_ieeenan(qq, FILL, 0)
        qq@_FillValue = FILL
        tmp = ismissing(qq)
        delete(qq)
        if (cc.eq.0) then
          lmask = tmp
        else
          lmask := lmask .and. tmp
        end if
        delete(tmp)
      end do

      ; Add up classes
      do cc = 0, dimsizes(class) - 1
        log_info("  adding class " + class(cc))
        tmp = f->$class(cc)$
        replace_ieeenan(tmp, FILL, 0)
        tmp@_FillValue = FILL
        tmp = where(ismissing(tmp), 0., tmp)
        if (cc.eq.0) then
          xx = tmp
        else
          xx = xx + tmp
        end if
        delete(tmp)
      end do
      delete(class)

      ; Reapply mask of missing values
      xx = where(lmask, xx@_FillValue, xx)

      ; Define output array
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

      ; Format coordinates
      output!0 = "time"
      output!1 = "lat"
      output!2 = "lon"
      format_coords(output, YEAR1 + "0101", YEAR2 + "1231", FREQ(vv))

      ; Set variable attributes
      tmp = format_variable(output, VAR(vv), CMOR_TABLE)
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

      ; Special case for baresoilFrac: add auxiliary coordinate
      if (VAR(vv).eq."baresoilFrac") then
        output@coordinates = "type"
      end if

      ; Write variable
      write_nc(fout, VAR(vv), output, bounds, gAtt)
      delete(gAtt)
      delete(output)
      delete(bounds)

      ; Special case for baresoilFrac: add auxiliary coordinate
      if (VAR(vv).eq."baresoilFrac") then
        type = tochar("bare_ground")
        type!0 = "strlen"
        type@long_name = "surface type"
        type@standard_name = "area_type"
        w = addfile(fout, "w")
        w->type = type
        delete(w)
      end if

    end do
  end do

end
