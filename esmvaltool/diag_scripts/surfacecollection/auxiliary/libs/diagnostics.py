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
    """
    produces global_temporal_mean
    -----------------------------
    returns a list of mean cubes
    """
    cubes = []
    
    for c in data.get_all():
        cubes.append(c.collapsed("time",iris.analysis.MEAN))
    
    return cubes