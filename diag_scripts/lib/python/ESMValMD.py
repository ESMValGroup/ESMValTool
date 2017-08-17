#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 17:12:54 2017

@author: bmueller
"""
from METAdata import METAdata
import datetime
import string
import random
import sys
import os

class ESMValMD(METAdata):

    def __init__(self, dtype="xml", modfile=None, tags=['tag'],
                 caption='caption', blockID="#ID",
                 DataIDs="NO IDs found.", diag_name=None, contrib_authors=None,
                 directwrite=True, **kwargs):
        super(ESMValMD, self).__init__(**kwargs)

        self.provenance = None
        self.diag_name = diag_name
        self.contrib_authors = contrib_authors
        self.set_provenance()

        self._make_dict(tags, caption, blockID, DataIDs)

        self.set_type(dtype)

        self.set_file(modfile)

        if directwrite:
            self.write()

    def set_provenance(self):
        d = {}
        d['Software_versions'] = {}
        if '0_ESMValTool_version' in os.environ.keys():
            d['Software_versions']['ESMValTool'] = os.environ['0_ESMValTool_version']
        d['Software_versions']['Python'] = sys.version
        #de = {}
        #for t in dict(os.environ).iteritems():
        #    de[str(t[0])] = str(t[1])
        #d['Environment'] = de
        d['Diag_name'] = self.diag_name
        d['Contrib_authors'] = self.contrib_authors
        self.provenance = d

    def _make_dict(self, tags, caption, blockID, DataIDs):
        DICT = {'ESMValTool': {
            'built': str(datetime.datetime.utcnow()),
            'tags': tags,
            'caption': caption,
            'block': blockID,
            'DataIDs': DataIDs , 'Provenance': self.provenance
            }}
        self.set_dict(DICT)


class nclFileMD(ESMValMD):

    def __init__(self, MDfile=None, **kwargs):

        with open(MDfile, 'rU') as f:
            file_lines = f.readlines()
            file_lines = [l.strip() for l in file_lines]
            try:
                file_name = file_lines[0]
            except:
                file_name = None
            try:
                dtype = file_lines[1]
            except:
                dtype = "both"
            try:
                tags = file_lines[2].split(",")
            except:
                tags = None
            try:
                caption = file_lines[3]
            except:
                caption = "No caption found for file: " + file_name + "!"
            try:
                blockID = file_lines[4]
            except:
                blockID = "#IDrand_" + ''.join(random.SystemRandom().
                                               choice(string.ascii_uppercase +
                                                      string.digits) for
                                               _ in range(8))
            try:
                DataIDs = file_lines[5]
            except:
                DataIDs = "No IDs found!"
            try:
                diag_name = file_lines[6]
            except:
                diag_name = "No diagnostic name found!"
            try:
                contrib_authors = file_lines[7]
            except:
                contrib_authors = "No contributing authors found!"

        super(nclFileMD, self).__init__(
            dtype=dtype,
            modfile=file_name,
            tags=tags,
            caption=caption,
            blockID=blockID,
            DataIDs=DataIDs,
            diag_name=diag_name,
            contrib_authors=contrib_authors
        )
