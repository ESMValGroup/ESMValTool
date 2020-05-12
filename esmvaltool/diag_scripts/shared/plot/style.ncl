; #############################################################################
; GENERAL ROUTINES FOR PLOT STYLE MANAGEMENT
; #############################################################################
; Please consider using of extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
; Contents:
;    function project_style
;    function place_debuginfo
;    procedure place_description
;    function gsnColorRange
;    function format_units
;    function set_log_ticks
;    function sort_alphabetically
;    function unique_labels_min
;    function unique_labels_all
;    function get_unique_items
;
; #############################################################################

load "$diag_scripts/../interface_scripts/auxiliary.ncl"
load "$diag_scripts/../interface_scripts/logging.ncl"

; #############################################################################
undef("project_style")
function project_style(items: list,
                       info[1]:logical,
                       flag[1]:string)
;
; Arguments
;    items: list of input_file_info items
;    info: diag_script_info logical
;    flag: string determining the type of array requested:
;            "annots": annotation strings.
;            "colors": colors (named colors, RGB or RGBA codes)
;            "dashes": line dash patterns.
;            "thicks": line thicknesses.
;            "markers": marker indexes.
;            "avgstd": average/standard deviation flags (0 = takes part in
;                      the calculation of mean and standard deviation, 1 =
;                      does not take part; usually 0 is for models and 1 for
;                      observations and reanalyses).
;
; Return value
;    An array of the same size of items list, with the style
;    information for the given flag. The type depends on the flag.
;
; Description
;    Retruns style informations (annotations, colors, line dash patterns, line
;    thicknesses, marker indexes and avgstd flat) based on a given styleset.
;    The styleset is determined based on the following priority list:
;      1st: style information for the given flag explicitely set as
;           diag_script_info@$flag$
;      2nd: styleset explicitely set as diag_script_info@styleset
;      3rd: styleset not defined, set to DEFAULT
;
; Caveats
;
; References
;
; Modification history
;    20150512-righi_mattia: modified to read style info from external style
;                           files, instead of using hard-coded values in the
;                           code. Functionalities of the project_styleset and
;                           project_style_<styleset> functions porteed here.
;    20130419-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, result, info, flag, styleset, stylefile, tmp, sty
begin

  funcname = "project_style"
  scriptname = "diag_scripts/shared/plot/style.ncl"
  enter_msg(scriptname, funcname)

  ; Check for the available dictonary keys to be used for the annotations
  dkeys = (/"dataset", "ensemble", "exp", "mip"/)

  ; Style information for this flag available as diag_script_info@$flag$
  if (isatt(info, flag)) then
    result = info@$flag$
    leave_msg(scriptname, funcname)
    return(result)
  end if

  ; Use diag_script_info@styleset or set to DEFAULT
  if (isatt(info, "styleset")) then
    styleset = diag_script_info@styleset
  else
    styleset = "DEFAULT"
  end if

  ; Styleset defined, search for corresponding style file
  styleset = str_lower(styleset)
  stylefile = "$diag_scripts/shared/plot/styles/" + styleset + ".style"
  if (.not.fileexists(stylefile)) then
    error_msg("f", scriptname, funcname, \
              "style file " + stylefile + " not found")
  end if
  tmp = readAsciiTable(stylefile, 1, "string", \
                       dimsizes(readAsciiHead(stylefile, "MODEL NAME")) + 1)
  sty = tmp(:, 0)
  delete(tmp)

  ; Model annotations, return unique model names based on dict keys
  if (flag.eq."annots") then

    ; Start with highest priority
    result = metadata_att_as_array(items, dkeys(0))

    unique = get_unique_values(result)
    iprio = 1
    do while (dimsizes(unique) .ne. dimsizes(result))
      if (iprio .eq. dimsizes(dkeys)) then
        error_msg("w", scriptname, funcname, "non-unique labels in dataset" + \
                  "annotations")
        unique := result
      else
        do i = 0, dimsizes(unique) - 1
          index = ind(result .eq. unique(i))
          if (dimsizes(index) .gt. 1) then  ; More than one occurence
            append = metadata_att_as_array(items, dkeys(iprio))
            result(index) = result(index) + "_" + append(index)
          end if
          delete(index)
        end do
        delete(unique)
        unique = get_unique_values(result)
        iprio = iprio + 1
      end if
    end do

    leave_msg(scriptname, funcname)
    return(result)
  end if

  ; Read model names
  modelsinfile = str_squeeze(str_get_field(sty, 1, "|"))

  ; Set flag names
  flagnames = (/"annots", "colors", "dashes", "thicks", "markers", "avgstd"/)
  colid = ispan(1, dimsizes(flagnames), 1)

  ; Read column
  col = ind(flag.eq.flagnames) + 1
  if (ismissing(col)) then
    error_msg("f", scriptname, funcname, "flag " + flag + " not allowed")
  end if
  datasets = metadata_att_as_array(items, "dataset")
  ndatasets = dimsizes(datasets)
  column = new(ndatasets, string)
  wildcard_cnt = 0
  do ii = 0, ndatasets - 1
    row = ind(datasets(ii).eq.modelsinfile)

    ; Multiple entries, raise error
    if (dimsizes(row).gt.1) then
      error_msg("f", scriptname, funcname, "Multiple entries for dataset " + \
                datasets(ii) + " in style file " + stylefile)
    end if

    ; Missing dataset
    if (ismissing(row)) then

      ; Search for "default"
      row  = ind("default".eq.modelsinfile)
      if (.not.ismissing(row)) then
        log_debug("Dataset " + datasets(ii) + \
                  " not found in " + stylefile + ", using default")
        column(ii) = str_squeeze(str_get_field(sty(row), col, "|"))
        continue
      end if

      ; Search for wildcards "*"
      rows = ind("*".eq.modelsinfile)
      if (all(ismissing(rows))) then
        error_msg("f", scriptname, funcname, "Dataset " + \
                  datasets(ii) + " not found in stylefile" + \
                  stylefile)
      end if
      if (wildcard_cnt.ge.dimsizes(rows)) then
        error_msg("f", scriptname, funcname, "style file " + stylefile + \
                  " contains too few wildcard (*) " + "entries (" + \
                  wildcard_cnt + "), please add more")
      end if
      log_debug("Dataset " + datasets(ii) + " not found in " + \
                stylefile + ", using wildcard")
      row = rows(wildcard_cnt)
      wildcard_cnt = wildcard_cnt + 1
      column(ii) = str_squeeze(str_get_field(sty(row), col, "|"))
      continue

    end if
    column(ii) = str_squeeze(str_get_field(sty(row), col, "|"))
  end do

  ; Colors: normalize to 1 and covert to string
  if (flag.eq."colors") then
    result = column
    do ii = 0, dimsizes(result) - 1
      nc = str_fields_count(column(ii), ",")
      ; Named color
      if (nc.eq.1) then
        result(ii) = column(ii)
      else if (nc.eq.3) then
        c1 = tofloat(str_get_field(column(ii), 1, ",")) / 255.
        c2 = tofloat(str_get_field(column(ii), 2, ",")) / 255.
        c3 = tofloat(str_get_field(column(ii), 3, ",")) / 255.
        result(ii) = "(/" + \
          sprintf("%7.5f", c1) + "," + \
          sprintf("%7.5f", c2) + "," + \
          sprintf("%7.5f", c3) + "/)"
      else if (nc.eq.4) then
        c1 = tofloat(str_get_field(column(ii), 1, ",")) / 255.
        c2 = tofloat(str_get_field(column(ii), 2, ",")) / 255.
        c3 = tofloat(str_get_field(column(ii), 3, ",")) / 255.
        c3 = tofloat(str_get_field(column(ii), 4, ","))
        result(ii) = "(/" + \
          sprintf("%7.5f", c1) + "," + \
          sprintf("%7.5f", c2) + "," + \
          sprintf("%7.5f", c3) + "," + \
          sprintf("%7.5f", c4) + "/)"
      else
        error_msg("f", scriptname, funcname, \
                  "Invalid definition of color for dataset " + \
                  datasets(ii) + ", check style file " + \
                  stylefile)
      end if
      end if
      end if
    end do
    leave_msg(scriptname, funcname)
    return(result)
  end if

  ; Line dash patterns, marker indexes, avgstd flag: convert to integer
  if (any(flag.eq.(/"dashes", "markers", "avgstd"/))) then
    result = toint(column)
    leave_msg(scriptname, funcname)
    return(result)
  end if

  ; Line thicknesses: convert to float
  if (flag.eq."thicks") then
    result = tofloat(column)
    leave_msg(scriptname, funcname)
    return(result)
  end if

end

; #############################################################################
undef("place_debuginfo")
function place_debuginfo(wks[1]:graphic,
                         debugstring[1]:string,
                         res[1]:logical,
                         plot[1]:graphic)
;
; Arguments
;     wks: current workstation.
;     debugstring: string to attach.
;     res: resource settings for display box.
;     plot: graphic object to draw text onto
;
; Return value
;
; Description
;    Places the text string debugstring onto wks.
;
; Caveats
;
; References
;
; Modification history
;
local i, txvpWidth, txvpHeight, tx, x, y, parent_width, parent_height, \
  annoid, mvxf, mvyf, debugstring
begin

  funcname = "place_debuginfo"
  scriptname = "diag_scripts/shared/plot/style.ncl"
  enter_msg(scriptname, funcname)

  res@gsnDraw               = False
  res@gsnFrame              = False
  res@txFontHeightF         = 0.012
  res@txBackgroundFillColor = 0
  res@txFuncCode            = "~"
  res@txPerimOn             = True

  ; First debugstring
  tx = gsn_create_text(wks, debugstring, res)
  getvalues tx
    "vpWidthF" : txvpWidth
    "vpHeightF" : txvpHeight
  end getvalues
  getvalues plot
    "vpWidthF" : parent_width
    "vpHeightF" : parent_height
  end getvalues
  mvxf = txvpWidth / parent_width
  mvyf = txvpHeight / parent_height

  amres = True
  amres@amZone           = 2
  amres@amSide           = "Right"
  amres@amParallelPosF   = mvyf / 2
  amres@amOrthogonalPosF = -mvxf
  annoid = gsn_add_annotation(plot, tx, amres)
  leave_msg(scriptname, funcname)
  return(annoid)

end

; #############################################################################
undef("place_description")
procedure place_description(wks[1]:graphic,
                            description[1]:string,
                            y_ndc_coord[1]:float)
;
; Arguments
;    wks: current workstation
;    description: string to attach
;    y_ndc_coord: vertical placement in ndc space (-1 for default)
;
; Return value
;
; Description
;     Places the text strings in array debugboxes onto wks
;
; Caveats
;
; References
;
; Modification history
;
local tx, txvpWidth, txvpHeight, parent_x, parent_y, amres, wks_dummy, \
  annoid, txres, y_ndc_coord
begin

  funcname = "place_description"
  scriptname = "diag_scripts/shared/plot/style.ncl"
  enter_msg(scriptname, funcname)

  wks_dummy = gsn_open_wks("png", "dummy.png")
  txres = True
  txres@gsnDraw               = True
  txres@gsnFrame              = False
  txres@txFontHeightF         = 0.012
  txres@txBackgroundFillColor = 0
  txres@txPerimOn             = True
  tx = gsn_create_text(wks_dummy, description, txres)
  getvalues tx
    "vpWidthF" : txvpWidth
    "vpHeightF" : txvpHeight
  end getvalues

  if (y_ndc_coord .lt. 0) then
    y_ndc_coord = 0.8
  end if
  gsn_text_ndc(wks, description, txvpWidth / 2, y_ndc_coord, txres)
  leave_msg(scriptname, funcname)

end

; #############################################################################
undef("gsnColorRange")
function gsnColorRange(lower:numeric,
                       upper:numeric,
                       step:numeric,
                       center:numeric,
                       color_end:integer,
                       center_color:integer)
;
; Arguments
;    lower: cnMinLevelValF.
;    upper: cnMaxLevelValF.
;    step: cnLevelSpacingF.
;    center: The numerical value the colormap is centered on. For
;            anomalies or trends, it's common to use 0.0, so blue means cold
;            or cooling and red means warm or warming.
;    color_end: The number of colors in colormap (ex. 97 for BlRe, 253 for
;               BlueRed).
;    center_color = Color value on the left of the "center" value (see above).
;
; Description
;    Sets the gsnSpreadXXX resources necessary to correctly span a two-color
;    colortable.
;
; Return value
;    A logical variable with the start and end colors as attributes.
;
; Caveats
;
; References
;    http://www.ncl.ucar.edu/Applications/Scripts/contoursym_4.ncl
;
; Modification history
;    20130422-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, nboxes_left, nboxes_right, ncolors_left, \
  ncolors_right, output, color_start
begin

  funcname = "gsnColorRange"
  scriptname = "diag_scripts/shared/plot/style.ncl"
  enter_msg(scriptname, funcname)

  color_start = 2  ; All of the color maps begin with 2.

  ; Calculate the number of color bar boxes to the left and right of
  ; the "center" value.
  nboxes_left  = (center - (lower - step)) / step
  nboxes_right = ((upper + step) - center) / step

  ; Calculate the number of colors in the map on
  ; the left and right hand sides.
  ncolors_left = (center_color - color_start + 1)
  ncolors_right = (color_end - center_color)

  output = True
  ; Either the lower or upper extent of the color map will be adjusted. If
  ; the magnitude of the lower limit is less than the magnitude of the
  ; upper limit, then the lower limit has to be "moved in" towards the
  ; center color. Oppositely, the upper limit will be moved. If both the
  ; lower and upper numerical values are the same, then pass back 2 as the
  ; lower extent and color_end (number of colors) as the upper extent (use
  ; the whole map in other words).
  if (abs(lower) .lt. abs(upper)) then
    output@ColorStart = \
      round(center_color - (ncolors_right / nboxes_right) * \
            nboxes_left, 3)
      output@ColorEnd   = color_end
  else if (abs(lower) .gt. abs(upper)) then
    output@ColorStart = 2
    output@ColorEnd = \
      round(center_color + (ncolors_left / nboxes_left) * \
            nboxes_right, 3)
  else
    output@ColorStart = 2
    output@ColorEnd   = color_end
  end if
  end if

  leave_msg(scriptname, funcname)
  return(output)

end

; #############################################################################
undef("format_units")
function format_units(str[1]: string)
;
; Arguments
;    str: a string.
;
; Return value
;    A string.
;
; Description
;    Reformats units string to properly display superscripts
;    (e.g. m^2 --> m~S1~2)
;
; Caveats
;    Currently convering only very few cases, to be extended.
;
; References
;
; Modification history
;    20140320-righi_mattia: written.
;
local funcname, scriptname, loc_str
begin

  funcname = "format_units"
  scriptname = "diag_scripts/shared/plot/style.ncl"
  enter_msg(scriptname, funcname)

  ; Create a local copy
  loc_str = str

  if (isStrSubset(loc_str, "-STP")) then
    loc_str = str_sub_str(loc_str, "-STP", "")
    loc_str = str_sub_str(loc_str, "-1", "~S2~-1~B3~STP    ")
    loc_str = str_sub_str(loc_str, "-2", "~S2~-2~B3~STP    ")
    loc_str = str_sub_str(loc_str, "-3", "~S2~-3~B3~STP    ")
  else
    loc_str = str_sub_str(loc_str, "-1", "~S2~-1  ")
    loc_str = str_sub_str(loc_str, "-2", "~S2~-2  ")
    loc_str = str_sub_str(loc_str, "-3", "~S2~-3  ")
  end if

  if (isStrSubset(loc_str, "umol").or.isStrSubset(loc_str, "ug")) then
    loc_str = str_sub_str(loc_str, "u", "~F33~m~F21~")
  end if

  if (isStrSubset(loc_str, "degC")) then
    loc_str = str_sub_str(loc_str, "deg", "~F35~J~F~")
  end if

  leave_msg(scriptname, funcname)
  return(loc_str)

end

; #############################################################################
undef("set_log_ticks")
function set_log_ticks(fmin[1]:numeric,
                       fmax[1]:numeric,
                       opt[1]:string)
;
; Arguments
;    fmin: the minimum axis value (as specified by trXMinF or trYMinF)
;    fmax: the maximum axis value (as specified by trXMaxF or trYMaxF)
;    opt: "major" or "minor", to return major or minor ticks.
;
; Return value
;    An array of tickmark values, to be used in trXBValues, trYLValues and
;    trXBLabels, trYLLabels
;
; Description
;    Since NCL only allows for exponential notation in plots with logarithmic
;    axes, this function creates an explicit sets of tickmarks with float
;    notation.
;
; Caveats
;
; References
;
; Modification history
;    20141003-righi_mattia: written.
;
local funcname, scriptname, omin, omax, rmin, rmax, om, tmp, arr
begin

  funcname = "set_log_ticks"
  scriptname = "diag_scripts/shared/plot/style.ncl"
  enter_msg(scriptname, funcname)

  ; Find order of magnitude of min
  omin = toint(log10(fmin))
  if (fmin.lt.1) then
    omin = omin - 1
  end if

  ; Find order of magnitude of max
  omax = toint(log10(fmax))
  if (fmax.lt.1) then
    omax = omax - 1
  end if

  if (opt.eq."major") then
    out = 10 ^ ispan(omin + 1, omax, 1)
    leave_msg(scriptname, funcname)
    return(out)
  end if

  if (opt.eq."minor") then
    rmin = (toint(fmin / 10 ^ omin) + 1) * 10 ^ omin
    rmax = (toint(fmax / 10 ^ omax)) * 10 ^ omax
    do om = omin, omax
      step = ispan(1, 9, 1) * 10 ^ om
      if (isdefined("arr")) then
        tmp = array_append_record(arr, step, 0)
        delete(arr)
        arr = tmp
        delete(tmp)
      else
        arr = step
      end if
    end do
    out = arr(ind(arr.ge.rmin .and. arr.le.rmax))
    leave_msg(scriptname, funcname)
    return(out)
  end if

  error_msg("f", scriptname, funcname, "Unrecognized option " + opt)

end

; #############################################################################
undef("sort_alphabetically")
function sort_alphabetically(orig_names[*], idx_exclude, dest_exclude)
;
; Arguments
;    orig_names: the array of model names prior to sorting
;    idx_exclude: the index(es) to be excluded from sorting, -1 to include
;                 everything
;    dest_exclude: the position where to put the excluded values after
;                  sorting ("begin" or "end")
;
; Return value
;    An integer array of the sime size of orig_names, with the permutation
;    index to be used to sort the array in alphabetical order.
;
; Description
;    Given an array of model names, this function returns the permutation
;    indexes which can be used to sort the array in alphabetical order.
;    Certain elements of the array can be excluded from the sorting and
;    placed either at the beginning or at the end of the sorted array (e.g.,
;    for sorting model alphabetically but leaving observations at the end,
;    or multi-model mean at the beginning).
;    The function itself does NOT perform any sorting, it just returns the
;    permutation indexes. These have to be applied to both the data AND
;    the model coordinate to get consistent results.
;    For example:
;
;        data(models|:, lat|:, lon|:)
;        pid = sort_alphabetically(data&models, -1, "")
;        sorted_data = data(pid, :, :)
;        sorted_data&models = data&models(pid)
;
; Caveats
;    Overwriting the original data can lead to incorrect results:
;        data = data(pid, :, :)          ; THIS IS WRONG!
;        data&models = data&models(pid)  ; THIS IS WRONG!
;
; References
;
; Modification history
;    20151028-righi_mattia: written.
;
begin

  funcname = "sort_alphabetically"
  scriptname = "diag_scripts/shared/plot/style.ncl"
  enter_msg(scriptname, funcname)

  ; Create a local copy of input names and set to upper case
  loc_names = str_upper(orig_names)
  loc_names@_FillValue = default_fillvalue(typeof(loc_names))
  nn = dimsizes(loc_names)

  ; Check if there are entries to be excluded
  if (all(ismissing(idx_exclude))) then
    l_excl = False
  else if (all(idx_exclude.lt.0)) then
    l_excl = False
  else
    l_excl = True
  end if
  end if

  if (l_excl) then  ; exclude given entries
    sub_names = loc_names
    sub_names(idx_exclude) = sub_names@_FillValue
    exc_names = loc_names(idx_exclude)
    inc_names = sub_names(ind(.not.ismissing(sub_names)))
    sqsort(inc_names)
    if (dest_exclude.eq."begin") then
      reord_names = array_append_record(exc_names, inc_names, 0)
    else if (dest_exclude.eq."end") then
      reord_names = array_append_record(inc_names, exc_names, 0)
    else
      error_msg("f", scriptname, funcname, "unrecognized option" + \
                dest_exclude)
    end if
    end if
  else  ; consider all
    reord_names = loc_names
    sqsort(reord_names)
  end if

  ; Calculate permutation indexes
  perm_idx = new(nn, integer)
  do ii = 0, nn - 1
    perm_idx(ii) = ind(reord_names(ii).eq.loc_names)
  end do

  enter_msg(scriptname, funcname)
    return(perm_idx)

end

; #############################################################################
undef("unique_labels_min")
function unique_labels_min(items: list, \
                           prio: string)
;
; Arguments
;    items: list of input_file_info metadata
;    prio: string vector with attribute names of elements of items, ordered by
;          priority for annotation (starting with highest).
;
; Return value
;    A vector (string) with one element for each item element, each label
;    contains the least possible attribute strings.
;
; Description
;
; Caveats
;
; References
;
; Modification history
;    20181217-schlund_manuel: adapated to new NCL interface of v2.0.
;    20130422-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, prio, result, unique, iprio, index
begin

  funcname = "unique_labels_min"
  scriptname = "diag_scripts/shared/plot/style.ncl"
  enter_msg(scriptname, funcname)

  ; Start with highest priority
  result = metadata_att_as_array(items, prio(0))
  unique = get_unique_values(result)

  ; Iterate over priorities
  iprio = 1
  do while (dimsizes(unique) .ne. dimsizes(result))

    ; Not enough prios given
    if (iprio .eq. dimsizes(prio)) then
      error_msg("w", scriptname, funcname, "Add more attributes to prio " + \
                "to make labels unique! Continuing with non-unique labels")
      break

    ; Add further prio
    else
      do i = 0, dimsizes(unique) - 1
        index := ind(result .eq. unique(i))
        new_prio = metadata_att_as_array(items, prio(iprio))

        ; More than one occurence
        if (dimsizes(index) .gt. 1) then
          result(index) = result(index) + "_" + new_prio(index)
        end if
      end do

      ; Reset
      unique := get_unique_values(result)
      iprio = iprio + 1
    end if
  end do

  leave_msg(scriptname, funcname)
  return(result)

end

; #############################################################################
undef("unique_labels_all")
function unique_labels_all(items: list, \
                           prio:string)
;
; Arguments
;    items: list of input_file_info metadata
;    prio: string vector with attribute names of elements of items, ordered by
;          priority for annotation (starting with highest).
;
; Return value
;    A vector (string) with one element for each item element, all labels
;    contain the same (least possible) number of attribute strings.

; Description
;
; Caveats
;
; References
;
; Modification history
;    20181217-schlund_manuel: adapated to new NCL interface of v2.0.
;    20130422-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, prio, result, unique, iprio, index
begin

  funcname = "unique_labels_all"
  scriptname = "diag_scripts/shared/plot/style.ncl"
  enter_msg(scriptname, funcname)

  ; Start with highest priority
  result = metadata_att_as_array(items, prio(0))
  unique = get_unique_values(result)

  ; Iterate over priorities
  iprio = 1
  do while (dimsizes(unique) .ne. dimsizes(result))

    ; Not enough prios given
    if (iprio .eq. dimsizes(prio)) then
      error_msg("w", scriptname, funcname, "Add more attributes to prio " + \
                "to make labels unique! Continuing with non-unique labels")
      break

    ; Add further prio
    else
      do i = 0, dimsizes(unique) - 1
        index := ind(result .eq. unique(i))
        new_prio = metadata_att_as_array(items, prio(iprio))
        result(index) = result(index) + "_" + new_prio(index)
      end do

      ; Reset
      unique := get_unique_values(result)
      iprio = iprio + 1
    end if
  end do

  leave_msg(scriptname, funcname)
  return(result)

end

; #############################################################################
undef("get_unique_items")
function get_unique_items(items: list,
                          prio: string)

;
; Arguments
;    items: list of input_file_info metadata
;    prio: string vector with attribute names of elements of items, ordered by
;          priority for annotation (starting with highest).
;
; Description
;    Return all unique list elements with the given prio.
;
; Caveats
;
; References
;
; Modification history
;    20181218-schlund_manuel: written
;
local items, prio, funcname, scriptname, unique_labels, ilist, unique_list, \
  atts, info
begin

  funcname = "get_unique_items"
  scriptname = "diag_scripts/shared/plot/style.ncl"
  enter_msg(scriptname, funcname)

  ; Get unique lables
  unique_labels = unique_labels_min(items, prio)
  do ilist = 0, ListCount(items) - 1
    items[ilist]@label = unique_labels(ilist)
  end do

  ; Select unique labels (always return 1st element)
  unique_list = NewList("lifo")
  unique_labels := get_unique_values(unique_labels)
  do ilab = 0, dimsizes(unique_labels) - 1
    atts := True
    atts@label = unique_labels(ilab)
    info = select_metadata_by_atts(items, atts)
    ListAppend(unique_list, info[0])
  end do

  leave_msg(scriptname, funcname)
  return(unique_list)

end
