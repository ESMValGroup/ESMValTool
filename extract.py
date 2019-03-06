"""
call like:
    python catch_bad_files.py FOLDER
to check if imagefiles in FOLDER has metadata tag and delete if not.

"""
from PIL import Image
from PIL import PngImagePlugin
from PIL.ExifTags import TAGS, GPSTAGS
import argparse
import os
from os.path import join, getsize

parser = argparse.ArgumentParser(description='Extract Metadata')
parser.add_argument('indir', metavar='INDIR', type=str, nargs=1,
                                    help='indir containing files with metadata')

args = parser.parse_args()
indir= args.indir[0]
bad_files = list()

for root, dirs, files in os.walk(indir):
    for f in files:
        absolut_path = os.path.join(root, f)
        image = Image.open(absolut_path)
        info = image.info
        if len(info) == 0:
            print(f + ": MISSING KEY")
            # os.system("rm " + os.path.join(root,f))
        for tag, value in info.items():
            try:
                key = TAGS.get(tag, tag)
            except:
                key = "BAD"
            print(f + " : " + key + " " + str(value))
            #if key == "Exif.Image.ImageDescription":
            #    #print(f + ":" + key + " " + str(value))
            #    #print(f + ":" + key)
            #    pass
            #else:
            #    #print(f + " : NO KEY!")
            #    print(f + ":" + key + " " + str(value))
            #    #if os.path.islink(absolut_path):
            #    #    os.remove(absolut_path)
            #    #    print("Removing link: " + absolut_path)

