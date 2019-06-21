#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 13:05:25 2019

@author: bmueller
"""
import iris
import logging
import os

logger = logging.getLogger(os.path.basename(__file__))

def glob_temp_mean(data):
    return