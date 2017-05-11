#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 17:12:54 2017

@author: bmueller
"""
from METAdata import METAdata
import datetime


class ESMValMD(METAdata):

    def __init__(self, dtype="xml", modfile=None, tags=['tag'],
                 caption='caption', blockID="#ID",
                 DataIDs="NO ID found.", directwrite=True, **kwargs):
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
            file_lines = [l.split()[0] for l in file_lines]
            file_name = file_lines[0]
            dtype = file_lines[1]
            tags = file_lines[2].split(",")
            caption = file_lines[3]
            blockID = file_lines[4]
            DataIDs = file_lines[5]

        super(nclFileMD, self).__init__(
                dtype=dtype,
                modfile=file_name,
                tags=tags,
                caption=caption,
                blockID=blockID,
                DataIDs=DataIDs
                )
