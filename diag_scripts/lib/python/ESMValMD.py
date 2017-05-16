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


class ESMValMD(METAdata):

    def __init__(self, dtype="xml", modfile=None, tags=['tag'],
                 caption='caption', blockID="#ID",
                 DataIDs="NO IDs found.", directwrite=True, **kwargs):
        super(ESMValMD, self).__init__(**kwargs)

        self._make_dict(tags, caption, blockID, DataIDs)

        self.set_type(dtype)

        self.set_file(modfile)

        if directwrite:
            self.write()

    def _make_dict(self, tags, caption, blockID, DataIDs):
        DICT = {'ESMValTool': {
            'built': str(datetime.datetime.utcnow()),
            'tags': tags,
            'caption': caption,
            'block': blockID,
            'DataIDs': DataIDs
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

        super(nclFileMD, self).__init__(
            dtype=dtype,
            modfile=file_name,
            tags=tags,
            caption=caption,
            blockID=blockID,
            DataIDs=DataIDs
        )
