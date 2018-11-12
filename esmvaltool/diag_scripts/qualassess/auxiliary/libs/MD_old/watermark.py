import datetime
import sys
import os
import math


def make_watermark(filepath, folder="same", type="default"):
    """Function adding watermarks to image.

    Parameters:
    filepath: string containing filepath including file extension
    folder: folder for watermark picture. "same" overwrites old file (optional)
    type: Type of watermark. "default" is the ESMValTool Logo (optional)"""
    # Modification history
    #    2010824-A_gier_be: written

    # Check output folder and create if necessary
    infile = filepath.split("/")
    if folder == "same":
        outfile = filepath
    else:
        if not os.path.isdir(folder):
            os.mkdir(folder)
        outfile = folder + "/" + infile[-1]

    # Check which watermark to use
    wm_folder = "diag_scripts/lib/ncl/watermarks/"
    if type == "cmip6_prelim":
        wm = wm_folder + "wm_esmvaltool_logo+text.png"
    else:
        wm = wm_folder + "wm_esmvaltool_logo.png"

    # Retrieve width and height of image to be watermarked
    try:
        wdth = int(os.popen('identify -format "%w\n" ' + filepath).read())
    except (ValueError):
        return

    hght = int(os.popen('identify -format "%h\n" ' + filepath).read())

    # Computer angle for diagonal
    alpha = int(math.degrees(math.atan(float(wdth) / float(hght))))

    # Rotate and Resize watermark to fit image
    os.system('convert -background "#00000000" -rotate -' +
              str(int(alpha)) + ' ' + wm + ' interim1.png')
    os.system("mogrify -trim +repage interim1.png")
    os.system("convert interim1.png -resize " + str(wdth) +
              "x" + str(hght) + " interim1.png")

    # Set date + version
    date = datetime.date.today().strftime('%Y-%m-%d')
    version = os.environ['0_ESMValTool_version']
    #tag = os.popen('git describe --abbrev=0 --tags 2>/dev/null').read()
    # if len(tag)==0:
    #    #Set a default version if no tags/no git respository
    #    print("No version tag set/Not in git repository. Assuming v.1.1.0")
    #    tag = "v1.1.0"

    # Make Version + Date stamp akin to example in
    # http://www.imagemagick.org/Usage/annotating/#wmark_text
    fontheight = hght / 24
    os.system('convert -size ' + str(wdth) + "x" + str(hght) +
              " xc:grey30 -pointsize " + str(fontheight) +
              " -gravity center -draw " +
              '"fill grey70 text 0,0 ' + "'ESMValTool " +
              version + "'" + '" -draw "fill grey70  text 0,' +
              str(fontheight) + " '" + date + "'" + '" stamp_fgnd.png')
    os.system("convert -size " + str(wdth) + "x" + str(hght) +
              " xc:black -pointsize " + str(fontheight) +
              " -gravity center -draw " + '"' +
              " fill white  text  1,1  'ESMValTool " +
              version + "'" + " text  1," + str(fontheight + 1) +
              " '" + date + "' text  0,0  'ESMValTool " +
              version + "' fill black text  0," + str(fontheight - 1) +
              " '" + date + "' text -1,-1 'ESMValTool " +
              version + "'" + '"' + " +matte stamp_mask.png")
    os.system("composite -compose CopyOpacity  stamp_mask.png" +
              " stamp_fgnd.png  stamp.png")
    os.system("mogrify -trim +repage stamp.png")

    # Clean up interim stamp images
    os.system("rm stamp_fgnd.png")
    os.system("rm stamp_mask.png")

    # Apply logo and stamp to image
    os.system("composite -dissolve 35% -gravity center interim1.png " +
              filepath + " " + outfile)
    os.system("composite -gravity SouthEast -geometry +5+5 stamp.png " +
              outfile + " " + outfile)

    # Clean up interim image files
    os.system("rm interim1.png")
    os.system("rm stamp.png")
