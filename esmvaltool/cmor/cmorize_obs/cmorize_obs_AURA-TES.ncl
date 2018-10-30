; #############################################################################
; REFORMAT SCRIPT FOR AURA-TES OBSERVATIONAL DATA
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    http://reverb.echo.nasa.gov/reverb/
;
; Last access
;    20140107
;
; Download and processing instructions
;    Select the following search options:
;        Platforms & Instruments: AURA / TES / FTS
;        Processing levels: 3
;    Select the following datasets:
;        TES/Aura L3 O3 Monthly Gridded V001
;        TES/Aura L3 O3 Monthly Gridded V002
;        TES/Aura L3 O3 Monthly Gridded V003
;    Add all the files to the cart and download the script as a text file.
;    To download the data use: wget -i <textfile>.
;    From the downloaded data, pick only the *.he5 files and save them in the
;    path specified below as INDIR.
;
; Modification history
;    20180628-A_righ_ma: adapted for v2.
;    20140129-A_righ_ma: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + "/cmor/cmorize_obs/interface.ncl")

begin

  ; Source name
  OBSNAME = "AURA-TES"

  ; Tier
  TIER = 2

  ; Input dir (raw data)
  INDIR = input_dir_path + OBSNAME

  ; Output dir (CMOR-ized data)
  OUTDIR = output_dir_path

  ; Selected variable (standard name)
  VARS = (/"tro3"/)

  ; Name in the raw data
  NAME = (/"O3"/)

  ; Conversion factor
  CONVERSION = (/1.e9/)

  ; Gridtype
  GRIDTYPE = "Nadir"

  ; Period
  YEAR1 = 2005
  YEAR2 = 2009

  ; Field
  FIELD = "T3M"

end

begin

  print(input_dir_path)

  exit
;;   diag_script = "reformat_obs_AURA-TES.ncl"
;;   enter_msg(diag_script, "", 1)

;;   ; Create output dir
;;   if (.not.isfilepresent(OUTDIR)) then
;;     system("mkdir -p " + OUTDIR)
;;   end if

;;   ; Creat time coordinate
;;   timec = create_timec(YEAR1, YEAR2)
;;   datec = ut_calendar(timec, 0)

;;   ; Loop over variables
;;   do vID = 0, dimsizes(VARS) - 1

;;     info_output("Processing " + VARS(vID), verbosity, 1)

;;     ; Loop over time
;;     do tt = 0, dimsizes(timec) - 1

;;       yy = toint(datec(tt, 0))
;;       mm = toint(datec(tt, 1))

;;       ; Find files
;;       fnames = INDIR + "TES-Aura_L3-" + NAME(vID) + "-M" + yy + "m" + \
;;         sprinti("%0.2i", mm) + "_F01_*.he5"
;;       files = systemfunc("ls " + fnames)

;;       if (.not.all(ismissing(files))) then

;;         ; Select most recent file version
;;         vv = -1000
;;         do ii = 0, dimsizes(files) - 1
;;           vi = toint(str_get_field(str_sub_str( \
;;                      systemfunc("basename " + files(ii)), ".he5", ""), 4, "_"))
;;           if (vi.gt.vv) then
;;             fname = files(ii)
;;           end if
;;         end do

;;         ; Read variable
;;         f = addfile(fname, "r")
;;         str = NAME(vID) + "_" + GRIDTYPE + "Grid"
;;         var = f->$str$

;;         ; Reorder
;;         d1 = "XDim_" + GRIDTYPE + "Grid"
;;         d2 = "YDim_" + GRIDTYPE + "Grid"
;;         d3 = "nLevels_" + GRIDTYPE + "Grid"
;;         tmp = var($d1$|:, $d2$|:, $d3$|:)
;;         delete(var)
;;         var = tmp
;;         delete(tmp)

;;         ; Read vertical coordinate
;;         str = "Pressure_" + GRIDTYPE + "Grid"
;;         press = f->$str$
;;         var!2 = "plev"
;;         var&plev = todouble(press)

;;         ; Create array
;;         if (.not.isdefined("outvar")) then
;;           outdim = array_append_record(dimsizes(timec), dimsizes(var), 0)
;;           outvar = new(outdim, float)
;;           outvar@_FillValue = var@_FillValue
;;         end if

;;         outvar(tt, :, :, :) = var
;;         delete(var)

;;       else

;;         log_info(" --> File missing, skipping this date")
;;       end if

;;     end do  ; time

;;     ; Rename coordinates
;;     outvar!0 = "time"
;;     outvar!1 = "lon"
;;     outvar!2 = "lat"

;;     ; Convert units [mol/mol] --> [1e9]
;;     outvar = outvar * 1.e9
;;     outvar@units = "1e9"

;;     ; Reorder
;;     output = outvar(time|:, plev|:, lat|:, lon|:)
;;     output = output(:, ::-1, :, :)
;;     output@_FillValue = FILL
;;     output&time = timec
;;     delete(outvar)

;;     ; Set variable attributes
;;     tmp = var_attrib(output, VARS(vID))
;;     delete(output)
;;     output = tmp
;;     delete(tmp)

;;     ; Format time coordinate
;;     ctime = time_attrib(output&time)
;;     delete(output&time)
;;     output&time = ctime
;;     delete(ctime)

;;     ; Format vertical coordinate
;;     output!1 = "plev"
;;     output&plev = output&plev * 100.  ; [hPa] --> [Pa]
;;     output&plev@units = "Pa"
;;     cplev = plev_attrib(output&plev)
;;     delete(output&plev)
;;     output&plev = cplev
;;     delete(cplev)

;;     ; Format latitude coordinate
;;     clat = lat_attrib(output&lat)
;;     delete(output&lat)
;;     output&lat = clat
;;     delete(clat)

;;     ; Format longitude coordinate
;;     clon = lon_attrib(output&lon)
;;     delete(output&lon)
;;     output&lon = clon
;;     delete(clon)

;;     ; Set global attributes
;;     gAtt = True
;;     gAtt@history       = "Created on " + systemfunc("date")
;;     gAtt@host          = systemfunc("echo $HOST")
;;     gAtt@user          = systemfunc("echo $USER")
;;     gAtt@period        = YEAR1 + "-" + YEAR2
;;     gAtt@field         = FIELD
;;     gAtt@tier          = TIER
;;     gAtt@source        = "http://reverb.echo.nasa.gov/reverb/"
;;     gAtt@reference     = "Beer, R., IEEE Trans. Geosci. Rem. Sens., " + \
;;       "44, 1102-1105, doi:10.1109/TGRS.2005.863716, 2006."
;;     gAtt@title = OBSNAME + " satellite data reformatted for the ESMValTool"
;;     gAtt@conventions = "CF/CMOR"

;;     ; Outfile
;;     fout = OUTDIR + "OBS_" + OBSNAME + "_sat_1_" + FIELD + \
;;       "_" + VARS(vID) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

;;     ; Write variable
;;     write_nc(fout, VARS(vID), output, gAtt)
;;     delete(gAtt)
;;     delete(output)

;;   end do

;;   leave_msg(diag_script, "", 1)

end
