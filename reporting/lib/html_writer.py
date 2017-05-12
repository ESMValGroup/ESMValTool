# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 09:46:13 2016

@author: bmueller
"""

import numpy as np
import fnmatch
import os
import subprocess
from itertools import compress
from METAdata import METAdata
from difflib import SequenceMatcher
import index_blocks


class HTML_writer(object):

    def __init__(self, **kwargs):
        super(HTML_writer, self).__init__(**kwargs)
        """
        Default values
        """
        self.lines = []
        self._default_Dict = {'ESMValTool': {
            'built': '',
            'tags': [],
            'caption': 'No caption available',
            'block': '#IDUnknown'
        }}

    def _get_files_in_directory(self, directory, pattern, asstring=True):
        """ returns list and number of files with pattern in directory """

#        if directory[-1] != os.sep:
#            directory += os.sep
#
#        L = glob.glob(directory + pattern)
#        N = len(L)

        L = []
        for root, dirs, files in os.walk(directory):
            for l_directory in dirs:
                subL, subN = self._get_files_in_directory(os.path.join(root, l_directory),pattern)
                L.append(subL)
            for filename in fnmatch.filter(files, pattern):
                L.append(os.path.join(root, filename))

        N = len(L)

        if asstring:
            L = ' '.join(L)
        return L, N

    def _get_img_files(self, directory):

        L = []

        for pat in ['png', 'jpg', 'jpeg', 'tiff']:  # any other formats?

            L_loc, N = self._get_files_in_directory(
                directory, '*.' + pat, False)

            L = L + L_loc

        return L

    def _get_pdf_files(self, directory):

        L = []

        for pat in ['pdf', 'eps', 'ps']:  # any other formats?

            L_loc, N = self._get_files_in_directory(
                directory, '*.' + pat, False)

            L = L + L_loc

        return L

    def tocblock(self, text_lines):
        "setup for block content1"

        lines = ['{% block toc %}']
        lines.extend(text_lines)
        lines.extend(['{% endblock toc %}'])

        return lines

    def block1(self, text_lines):
        "setup for block content1"

        lines = ['{% block content1 %}']
        lines.extend(text_lines)
        lines.extend(['{% endblock content1 %}'])

        return lines

    def block2(self, text_lines):
        "setup for block content2"

        lines = ['{% block content2 %}']
        lines.extend(text_lines)
        lines.extend(['{% endblock content2 %}'])

        return lines

    def block3(self, text_lines):
        "setup for block content3"

        lines = ['{% block content3 %}']
        lines.extend(text_lines)
        lines.extend(['{% endblock content3 %}'])

        return lines

    def write_html(self, name):
        "setup for write html"
        filename = "./templates/" + name + ".html"
        np.savetxt(filename, self.lines, fmt="%s")

    def append_html(self, name, txt):
        "setup for write html"
        filename = "./templates/" + name + ".html"
        with open(filename, 'a') as f:
            np.savetxt(f, txt, fmt="%s")

    def make_PDF_list(self, folder, name, restrictor, host):
        "produce list of pdf files for html"
        # This is just a rough idea!
        # TODO: implement host, horizontal lines, blocks, captions etc.
        file_list = self._get_pdf_files(folder)
        print(file_list)

        file_list = self.correct_file_list(file_list, name, restrictor)

        file_list = ['<object data=".' + f +
                     '" type="application/pdf"></object>' for f in file_list]

        return file_list

    def make_IMG_list(self, folder, name, restrictor, host, time=None):
        "produce list of image files for html"

        print(folder)
        file_list = self._get_img_files(folder)
        print(file_list)

        file_list = self.correct_file_list(file_list, name, restrictor)
        print(file_list)

        if time is not None:
            file_list = self.correct_file_list(file_list, name, time)

        additional_files = []
        for fl in file_list:
            MD = METAdata("xml", fl)
            MD.__adjust_xml_file__()
            additional_files.append(MD.get_file())
        file_list = file_list + additional_files

        ofolder = "./static/images/" + name

        subprocess.call(['mkdir', '-p', ofolder])

        # gather images and xml into static folder
        [subprocess.call(['ln', '-sfn', f, ofolder]) for f in file_list]

        file_list = self._get_img_files(ofolder)

        MD = METAdata()

        blocks = [MD.read(f).get_dict()['ESMValTool']['block'][0]
                  for f in file_list]

        index = sorted(range(len(blocks)), key=blocks.__getitem__)

        file_list = [file_list[i] for i in index]
        blocks = [blocks[i] for i in index]

        breaks = self._similar_follower(blocks)

        breaks = [b < 0.9 for b in breaks]

        breaks = ['<hr />' if e else '' for e in breaks]

        file_list = ['<figure><a href=' + "http://" + host + "/" + f +
                     ' title=' +
                     MD.read(f).get_dict()['ESMValTool']['block'][0] +
                     '><img src=".' + f + '"></a><figcaption>' +
                     MD.read(f).get_dict()['ESMValTool']['caption'][0] +
                     ' (V. ' +
                     MD.read(f).get_dict()['ESMValTool']['built'][0] + ')' +
                     '</figcaption></figure>' for f in file_list]

        file_list = [item for sublist in map(
            None, file_list, breaks) for item in sublist][:-1]

        return file_list

    def correct_file_list(self, file_list, name, restrictor_or_time):
        "correct list of available files concerning given files"

        MD = METAdata()

        file_list.sort(key=lambda x: MD.read(
            x).get_dict()['ESMValTool']['block'][0])

        if isinstance(restrictor_or_time, (str, unicode)):
            actual_file_list = [(restrictor_or_time in MD.read(x).get_dict()[
                                 'ESMValTool']['tags']) for x in file_list]
        elif isinstance(restrictor_or_time, list):
            actual_file_list = [all(r in
                                    MD.read(x).get_dict()['ESMValTool']['tags']
                                    for r in restrictor_or_time)
                                for x in file_list]
        elif isinstance(restrictor_or_time, float):
            actual_file_list = [
                    MD.read(x).get_dict()['ESMValTool']['built'][0] >
                    restrictor_or_time for x in file_list]
        else:
            assert False, ("This restrictor object is not of expected type! " +
                           "(Neither str, list or date!)")

        file_list = list(compress(file_list, actual_file_list))

        return file_list

    def diag_html(self, name="tab", folder="./", add_cfg="config",
                  restrictor=[""], host="0.0.0.0:5000", time=None):
        "setup for the html page for a specific diagnostic"
        self.lines = ['{% extends "index.html" %}']
        # TODO input folder here

        IMG_file_list = self.make_IMG_list(
            folder, name, restrictor, host, time)
        PDF_file_list = []
        # PDF_file_list = self.make_PDF_list(folder,name,mindate,host) #TODO
        # MetaData with PDF

        self.lines.extend(self.block1(PDF_file_list + IMG_file_list))

        self.lines.extend(self.block2(['<p>This is the summary of the ' +
                                       'preliminary reporting service!</p>',
                                       '<p>Original namelist version. </p>']))

        add_info = ["<h2>Config:</h2>"]
        add_info.extend(add_cfg)

        self.lines.extend(self.block3(add_info))

        self.write_html(name)
        return

    def index_html(self, list_of_vars=["var1", "var2"],
                   list_of_keys=["key1", "key2"], version="0.0alpha",
                   list_of_authors=["person1", "person2"]):
        "setup for the html page for the index"

        self.lines = []

        # write head (static)
        Head = index_blocks.index_head(list_of_vars)
        self.lines.extend(Head)

        # write index (dynamic)
        for key in list_of_keys:
            if list_of_keys[key] is None:
                toc = index_blocks.index_tabs(list_of_keys, key)
                act_toc = self.tocblock(toc)
                self.append_html(key, act_toc)
            else:
                for diag_key in list_of_keys[key]:
                    toc = index_blocks.index_tabs(list_of_keys,
                                                  diag_key['key'])
                    act_toc = self.tocblock(toc)
                    self.append_html(diag_key['key'], act_toc)

        Tabs = self.tocblock("")
        self.lines.extend(Tabs)

        # write foot (static)
        Foot = index_blocks.index_foot(version, list_of_authors)
        self.lines.extend(Foot)

        self.write_html("index")
        return

    def home_html(self, left_box, middle_box, right_box):
        "setup for the html page for the home page"
        self.lines = ['{% extends "index.html" %}']

        terminal = ['<div class="console">']
        terminal.extend(['<p>'])
        terminal.extend(middle_box.split("\n"))
        terminal.extend(['</p>'])
        terminal.append('</div>')

        self.lines.extend(self.block1(terminal))
        self.lines.extend(self.block2([left_box]))

        add_info = ["<h2>Namelist(s):</h2>"]
        add_info.extend(["<p>" + l + "</p>" for l in right_box])
#        add_info2=["<h2>Config:</h2>"]
#        add_info2.extend(add_cfg)

#        add_info=add_info1+add_info2

        self.lines.extend(self.block3(add_info))

        self.write_html("home")
        return

    def tag_html(self, name, folderlist, tags, host):
        "setup for the html page for a specific tag"
        self.lines = ['{% extends "index.html" %}']
        # TODO input folder here

        tags = [tags] if isinstance(tags, (str, unicode)) else tags

        IMG_file_list = []
        for folder in folderlist:
            # for i in range(len(tags)):
            IMG_file_list.extend(self.make_IMG_list(folder, name, tags, host))
            PDF_file_list = []
            # PDF_file_list = self.make_PDF_list(folder,name,mindate,host)
            # #TODO MetaData with PDF

        all_files = PDF_file_list + IMG_file_list

        self.lines.extend(self.block1(all_files if len(all_files) > 0
                                      else ["<p>Nothing to show!</p>",
                                            "<p>Could not find any files " +
                                            "for tag " + ", ".join(tags) +
                                            "!</p>"]))

        self.lines.extend(self.block2(['<p>This is the summary of the ' +
                                       'preliminary reporting service!</p>',
                                       '<p>Report namelist version. </p>']))

        add_info = ["<h2>Placeholder!</h2>"]

        self.lines.extend(self.block3(add_info))

        self.write_html(name)

        self.lines = ['{% extends "index.html" %}']

        return

    def _similar_follower(self, v):

        if len(v) == 0:
            return v
        else:
            v1 = list(v)
            v1.pop(0)
            v2 = list(v)
            v2.pop(-1)
            return list(map(self._similar, v1, v2))

    def _similar(self, a, b):
        return SequenceMatcher(None, a, b).ratio()
