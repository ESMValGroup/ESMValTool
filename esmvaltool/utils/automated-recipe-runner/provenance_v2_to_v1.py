#!/usr/bin/env python
# -*- coding: utf-8 -*-
import datetime
import argparse
import math
import os

from collections import defaultdict
from PIL import Image
from PIL.PngImagePlugin import PngInfo
from xml.etree import ElementTree as ET
import yaml
from PIL import Image
from PIL import ImageFile
from PIL import PngImagePlugin
from PIL.PngImagePlugin import PngInfo

PngImagePlugin.MAX_TEXT_CHUNK = PngImagePlugin.MAX_TEXT_CHUNK * 12
ImageFile.SAFEBLOCK = ImageFile.SAFEBLOCK * 10
PngImagePlugin.MAX_TEXT_MEMORY = PngImagePlugin.MAX_TEXT_MEMORY * 10

REF_FILE = 'tags.yml'

PREFIXES = [
    '{https://www.esmvaltool.org/attribute}',
    '{http://www.w3.org/ns/prov#}',
]
TAGS = {
    'domains': 'DM',
    'statistics': 'ST',
    'short_name': 'V',
    'themes': 'T',
    'dataset': 'M',
    'references': 'D',
    'realms': 'R',
    'plot_types': 'PT',
}

REF_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), REF_FILE)
with open(REF_PATH, 'r') as FILE:
    TAGS_REF = yaml.safe_load(FILE)


def get_args():
    """Define the command line."""
    parser = argparse.ArgumentParser(
        description="Prov new to old converter.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '-i',
        '--inpath',
        required=True,
        help=
        'Path to the directory that contains the ESMValTool output png files.')
    parser.add_argument('-o',
                        '--outpath',
                        required=True,
                        help='Path to the output directory.')
    parser.add_argument('-r',
                        '--recipe',
                        required=True,
                        help='Name of recipe.')
    args = parser.parse_args()

    if not os.path.exists(args.inpath):
        raise Exception(f"Inpath {args.inpath} does not exist.")
    if not os.path.exists(args.outpath):
        os.makedirs(args.outpath)
    if isinstance(args.recipe, str) and args.recipe.startswith('recipe_'):
        pass
    else:
        raise Exception(
            f"Recipe name {args.recipe} must be a string starting with `recipe_`."
        )
    return args


def replace_tag(value, category):
    """Find and replace tag key for a value in a reference file."""
    tags = TAGS_REF.get(category, {})
    if value in tags.values():
        value = list(tags.keys())[list(tags.values()).index(value)]
    return value


def dict_to_etree(dict_):
    """Convert python dictionary to XML Element."""
    if not isinstance(dict_, dict):
        raise TypeError("Expected dict, got {}", type(dict_))
    elems = []
    for (key, val) in dict_.items():
        elem = ET.Element(key)
        if isinstance(val, dict):
            elem.extend(dict_to_etree(val))
        else:
            elem.text = str(val)
        elems.append(elem)
    return elems


def etree_to_dict(tree):
    """Convert XML structure to python dictionary."""
    d = {tree.tag: {} if tree.attrib else None}
    children = list(tree)
    if children:
        dd = defaultdict(list)
        for dc in map(etree_to_dict, children):
            for k, v in dc.items():
                dd[k].append(v)
        d = {tree.tag: {k: v[0] if len(v) == 1 else v for k, v in dd.items()}}
    if tree.attrib:
        d[tree.tag].update(('@' + k, v) for k, v in tree.attrib.items())
    if tree.text:
        text = tree.text.strip()
        if children or tree.attrib:
            if text:
                d[tree.tag]['#text'] = text
        else:
            d[tree.tag] = text
    return d


def remove_prefixes(string):
    """Remove prefixes."""
    for prefix in PREFIXES:
        string = string.replace(prefix, '')
    return string


def remove_prefixes_from_dict(dict_):
    """Remove prefixes in dictionary keys."""
    new_dict = {}
    for (key, val) in dict_.items():
        new_dict[remove_prefixes(key)] = val
    return new_dict


def extract_authors(list_):
    """Extract authors from list of dictionaries."""
    authors = []
    for elem in list_:
        elem = remove_prefixes_from_dict(elem)
        for val in elem.values():
            if 'author:' in val:
                authors.append(val.replace('author:', ''))
    return authors


def extract_diags(list_):
    """Extract diags from list of dictionaries."""
    diags = []
    for elem in list_:
        elem = remove_prefixes_from_dict(elem)
        if 'script_file' in elem:
            diags.append(elem['script_file'])
    return diags


def extract_files(dict_):
    """Extract files from dictionary."""
    dict_ = remove_prefixes_from_dict(dict_)
    if '@id' in dict_:
        file_ = os.path.basename(dict_['@id'])
        if file_.startswith('recipe'):
            return None
        else:
            return file_
    else:
        return None


def extract_projects(list_):
    """Extract project from list of dictionaries."""
    projects = []
    for elem in list_:
        elem = remove_prefixes_from_dict(elem)
        for val in elem.values():
            if 'project:' in val:
                val = val.replace('project:', '')
                val = replace_tag(val, 'projects')
                projects.append('P_{}'.format(val))
    return projects


def extract_recipes(list_):
    """Extract recipe from list of dictionaries."""
    recipes = []
    for elem in list_:
        for dict_ in elem.values():
            dict_ = remove_prefixes_from_dict(dict_)
            for val in dict_.values():
                if 'recipe:' in val:
                    recipes.append('namelist_{}'.format(
                        val.replace('recipe:', '')))
    return recipes


def extract_tags(dict_):
    """Extract tags from dictionary."""
    tags = []
    dict_ = remove_prefixes_from_dict(dict_)
    for (key, val) in dict_.items():
        if key in TAGS:
            val = val.strip()
            if '(' in val and ')' in val:
                val = val.replace('(', '')
                val = val.replace(')', '')
                vals = val.split("',")
                for v in vals:
                    v = v.replace("'", '')
                    if v == '':
                        continue
                    v = v.strip()
                    v = replace_tag(v, key)
                    tags.append('{}_{}'.format(TAGS[key], v))
            else:
                val = replace_tag(val, key)
                tags.append('{}_{}'.format(TAGS[key], val))
    return tags


def extract_tasks(list_):
    """Extract task from list of dictionaries."""
    tasks = []
    for elem in list_:
        for dict_ in elem.values():
            dict_ = remove_prefixes_from_dict(dict_)
            for val in dict_.values():
                if 'task:' in val:
                    tasks.append(val.replace('task:', ''))
    return tasks


def convert_to_v1(info):
    """Convert provencance information of v2 to v1."""
    content = info['ImageHistory']
    xml = ET.XML(content)
    xml_dict = etree_to_dict(xml)
    sub_dict = xml_dict[list(xml_dict.keys())[0]]
    sub_dict = remove_prefixes_from_dict(sub_dict)

    # Tags, files and scripts
    entity = sub_dict.get('entity', {})
    diags = extract_diags(entity)
    tags = []
    files = []
    for elem in entity:
        tags.extend(extract_tags(elem))
        new_file = extract_files(elem)
        if new_file is not None:
            files.append(new_file)

    # Authors and project
    agent = sub_dict.get('agent', [])
    authors = extract_authors(agent)

    # Project
    tags.extend(extract_projects(agent))

    # Recipe
    was_started_by = sub_dict.get('wasStartedBy', [])
    if not isinstance(was_started_by,list): was_started_by = [was_started_by]
    tags.extend(extract_recipes(was_started_by))

    # Software and caption
    software = info.get('Software', '')
    caption = info.get('ImageDescription', '')

    # Tasks
    was_derived_from = sub_dict.get('wasDerivedFrom', [])
    tasks = extract_tasks(was_derived_from)

    # Create strings
    for (idx, tag) in enumerate(tags):
        tags[idx] = tag.replace(',', '')
    tags = list(set(tags))
    tags = '|'.join(tags)
    files = list(set(files))
    files = '|'.join(files)
    authors = list(set(authors))
    authors = ','.join(authors)
    tasks = list(set(tasks))
    tasks = ','.join(tasks)
    diags = list(set(diags))
    diags = '|'.join(diags)

    # Creat dictionary
    prov = {
        'DataIDs': files,
        'Provenance': {
            'Contrib_authors': authors,
            'Diag_name': diags,
            'Software_versions': {
                'ESMValTool': software,
            },
        },
        'block': tasks,
        'built': 'N.A.',
        'caption': caption,
        'tags': tags,
    }
    return {'ESMValTool': prov}


def make_watermark(infile,
                   outfile,
                   root,
                   version,
                   watermark_file='wm_esmvaltool_logo.png'):
    """Function adding watermarks to image."""
    dir_path = os.path.dirname(os.path.realpath(__file__))
    wm_file = os.path.join(dir_path, watermark_file)
    interim_file = os.path.join(root, "interim.png")
    stamp_fgnd = os.path.join(root, "stamp_fgnd.png")
    stamp_mask = os.path.join(root, "stamp_mask.png")
    stamp = os.path.join(root, "stamp.png")
    # Retrieve width and height of image to be watermarked
    try:
        width = int(os.popen('identify -format "%w\n" ' + infile).read())
    except ValueError as exc:
        print(str(exc))
        return
    height = int(os.popen('identify -format "%h\n" ' + infile).read())

    # Computer angle for diagonal
    alpha = int(math.degrees(math.atan(float(width) / float(height))))

    # Rotate and Resize watermark to fit image
    os.system('convert -background "#00000000" -rotate -' + str(int(alpha)) +
              ' ' + wm_file + ' ' + interim_file)
    os.system("mogrify -trim +repage " + interim_file)
    os.system("convert " + interim_file + " -resize " + str(width) + "x" +
              str(height) + " " + interim_file)

    # Set date + version
    date = datetime.date.today().strftime('%Y-%m-%d')

    # Make Version + Date stamp akin to example in
    # http://www.imagemagick.org/Usage/annotating/#wmark_text
    fontheight = height / 24.0
    os.system('convert -size ' + str(width) + "x" + str(height) +
              " xc:grey30 -pointsize " + str(fontheight) +
              " -gravity center -draw " + '"fill grey70 text 0,0 ' +
              "'ESMValTool " + version + "'" +
              '" -draw "fill grey70  text 0,' + str(fontheight) + " '" + date +
              "'" + '" ' + stamp_fgnd)
    os.system("convert -size " + str(width) + "x" + str(height) +
              " xc:black -pointsize " + str(fontheight) +
              " -gravity center -draw " + '"' + " fill white  text  1,1  '" +
              version + "'" + " text  1," + str(fontheight + 1) + " '" + date +
              "' text  0,0  '" + version + "' fill black text  0," +
              str(fontheight - 1) + " '" + date + "' text -1,-1 '" + version +
              "'" + '"' + " +matte " + stamp_mask)
    os.system("composite -compose CopyOpacity  " + stamp_mask + " " +
              stamp_fgnd + "  " + stamp)
    os.system("mogrify -trim +repage " + stamp)

    # Clean up interim stamp images
    os.system("rm " + stamp_fgnd)
    os.system("rm " + stamp_mask)

    # Apply logo and stamp to image
    os.system("composite -dissolve 35% -gravity center  " + interim_file +
              " " + infile + " " + outfile)
    os.system("composite -gravity SouthEast -geometry +5+5 " + stamp + " " +
              outfile + " " + outfile)

    # Clean up interim image files
    os.system("rm " + interim_file)
    os.system("rm " + stamp)
    print("Added watermark '{}'".format(wm_file))


def main(args):
    """Main function."""
    for (root, dirs, files) in os.walk(args.inpath):
        new_root = os.path.join(args.outpath, args.recipe)
        if not os.path.exists(new_root):
            os.makedirs(new_root)
            print("Created '{}'".format(new_root))

        # Save header of files
        for file_ in files:
            old_file = os.path.join(root, file_)
            print("current file : " + old_file)
            if os.path.splitext(file_)[1] != '.png':
                print("Skipping '{}'".format(old_file))
                print("")
                continue
            new_file = os.path.join(new_root, file_)
            yml_file = os.path.join(root, os.path.splitext(file_)[0] + '.yml')
            with Image.open(old_file) as image:
                info = image.info

            # v1
            if 'Exif.Image.ImageDescription' in info:
                print("File '{}' already contains v1 provencance".format(
                    old_file))

            # v2
            elif 'ImageHistory' in info:
                print("Saving header of '{}'".format(old_file))
                xml = dict_to_etree(convert_to_v1(info))[0]
                xml_str = ET.tostring(xml).decode('utf-8')
                info['Exif.Image.ImageDescription'] = xml_str

                # Write yml
                with open(yml_file, 'w') as outfile:
                    yaml.dump(info, outfile, default_flow_style=False)
                    print("Wrote '{}'".format(yml_file))

                # Watermark
                make_watermark(old_file, new_file, root,
                               info.get('Software', 'ESMValTool v2'))

                # Write header
                with open(yml_file, 'r') as infile:
                    header = yaml.safe_load(infile)
                os.remove(yml_file)
                pnginfo = PngInfo()
                for (key, val) in header.items():
                    pnginfo.add_text(key, val)
                with Image.open(new_file) as new_image:
                    new_image.save(new_file, pnginfo=pnginfo)
                    print("Wrote '{}'".format(new_file))

            # Invalid file
            else:
                print("File '{}' does not contain v2 provenance information, "
                      "skipping".format(old_file))
            print("")


if __name__ == '__main__':
    args = get_args()
    main(args)
