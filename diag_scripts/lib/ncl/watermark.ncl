;;#############################################################################
;; GENERAL ROUTINES FOR ADDING WATERMARKS
;;#############################################################################
;; Please consider using of extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    procedure make_watermark
;;
;;#############################################################################

load "./interface_scripts/messaging.ncl"
load "./interface_data/ncl.interface"

;;#############################################################################

undef("make_watermark")
procedure make_watermark(filepath[1]: string,
                        folder[1]: string,
                         watermark[1]: string)
;;
;; Arguments
;;    filepath: filename of the picture to watermark
;;    folder: output folder, if "same" -> overwrites original image
;;
;; Return Value
;;    None
;;
;; Description
;;    Applies an ESMValTool watermark to given picture
;;    TODO: add optional outfile name? something to add to text stamp?
;;
;; Modification history
;;    20170821-A_gier_be: written

local alpha, wdth, hght, cmd_width, cmd_height, cmd_resize, cmd_tag, tag, \
      date, cmd_stmp_fgnd, cmd_stmp_mask, cmd_stamp_compose, cmd_stamp_trim, \
      cmd_apply_logo, cmd_apply_stmp, infile, outfile, fld_test, cmd_date, \
      wm,cmd_rotate

begin

    funcname = "make_watermark"
    scriptname = "diag_scripts/lib/ncl/watermark.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 6)


  ;;Check output folder and make it if necessary
  infile = str_split(filepath, "/")

  if folder.eq."same" then
    outfile = filepath
  else
    fld_test = systemfunc("test -d " + folder + "; echo $?")
    if fld_test.ne.0 then
      system("mkdir " + folder)
    end if
    outfile = folder + "/" + infile(dimsizes(infile)-1)
  end if

  ;;Retrieve width and height of image to be watermarked
  dq = str_get_dq()
  cmd_width = "identify -format " + dq + "%w\n" + dq + " " + filepath
  wdth_cmd = toint(systemfunc(cmd_width))
  wdth = wdth_cmd(0)

  cmd_height = "identify -format " + dq + "%h\n" + dq + " " + filepath
  hght_cmd = toint(systemfunc(cmd_height))
  hght = hght_cmd(0)

  ;;Compute angle for watermark to go across the diagonal (radian -> deg)
  alpha = atan(tofloat(hght)/tofloat(wdth)) / 3.14159 * 180

  ;;Check which watermark is desired, no angle on cmip6 one
  ;;CHANGE FILEPATHS IF NECESSARY
  wm_folder = "diag_scripts/lib/ncl/watermarks/"
  if watermark.eq."cmip6_prelim" then
    wm = wm_folder + "wm_esmvaltool_logo+text.png"
  else
    wm = wm_folder + "wm_esmvaltool_logo.png"
  end if

  ;;Rotate and Resize Watermark according to image
  cmd_rotate = "convert -background " + dq + "#00000000" + dq + \
               " -rotate -" + toint(alpha) + " " + wm + " interim1.png"
  system(cmd_rotate)
  system("mogrify -trim +repage interim1.png")
  
  cmd_resize = "convert interim1.png -resize " + wdth + "x" + hght + \
               "   interim1.png"
  system(cmd_resize)

  ;;Set Date + Version
  cmd_date = "date " + dq + "+%Y-%m-%d" + dq
  date = systemfunc(cmd_date)
  ;; version! If on git repository which has version tags, it's easy
  cmd_tag = "git describe --abbrev=0 --tags"
  tag = systemfunc(cmd_tag)
  ;;Set default value if no git tag available
  ;tag = "v1.1.0"
  if ismissing(tag) then
    tag = "v1.1.0"
  end if
  
  
  fontheight = hght/24
  ;;Make Version + Date stamp akin to example in
  ;;http://www.imagemagick.org/Usage/annotating/#wmark_text
  cmd_stmp_fgnd = "convert -size " + wdth + "x" + hght + \
                  " xc:grey30 -pointsize " + fontheight + " -gravity center" + \
                  " -draw " + dq + "fill grey70  text 0,0  'ESMValTool " + \
                  tag + "'" + dq + " -draw " + dq + "fill grey70  text 0," + \
                  fontheight + " '" + date + dq + " stamp_fgnd.png"
  system(cmd_stmp_fgnd)

  cmd_stmp_mask = "convert -size " + wdth + "x" + hght + \
                  " xc:black -pointsize " + fontheight + " -gravity center" + \
                  " -draw " + dq + " fill white  text  1,1  'ESMValTool " + \
                  tag + "'" + "text  1," + (fontheight + 1) + " '" + date + \
                  "' text  0,0  'ESMValTool " + tag + "'" +  \
                  " fill black text  0," + (fontheight - 1) + " '" + date + \
                  "' text -1,-1 'ESMValTool " + tag + "'" + dq + \
                  " +matte stamp_mask.png"
  system(cmd_stmp_mask)
  
  cmd_stamp_compose = "composite -compose CopyOpacity  stamp_mask.png" + \
                      " stamp_fgnd.png  stamp.png"
  system(cmd_stamp_compose)

  cmd_stamp_trim = "mogrify -trim +repage stamp.png"
  system(cmd_stamp_trim)

  ;;clean up interim stamp images
  system("rm stamp_fgnd.png")
  system("rm stamp_mask.png")


  cmd_apply_logo = "composite -dissolve 35% -gravity center interim1.png " + \
                   filepath + " " + outfile
  system(cmd_apply_logo)

  cmd_apply_stmp = "composite -gravity SouthEast -geometry +5+5 stamp.png " +\
                   outfile + " " + outfile
  system(cmd_apply_stmp)
  
  ;;Clean up all the interim files
  system("rm interim1.png")
  system("rm stamp.png")

  leave_msg(scriptname, funcname, 6)
end
