#!/usr/bin/env python
"""Python example diagnostic."""
import os

import yaml


# Get metrics dict, get dataset color(s) or get datasetnames
def get_mpqb_cfg(cfgtype, cfgkey):
    """cfgtype: colormap, datasetcolor, datasetname"""
    # First read cfg file
    cfg_filename = os.path.join(os.path.split(__file__)[0],
                                'mpqb_cfg_xch4.yml')
    with open(cfg_filename, 'r', encoding="utf-8") as handle:
        mpqb_cfg = yaml.safe_load(handle)
    # Defaults to specified defaults in yml file
    if cfgtype == 'colormap':
        if cfgkey in mpqb_cfg['colormaps']:
            return mpqb_cfg['colormaps'][cfgkey]
        return mpqb_cfg['colormaps']['default']
    # Defaults to alias (provided as cfgkey)
    if cfgtype == 'datasetname':
        if cfgkey in mpqb_cfg['datasetnames']:
            return mpqb_cfg['datasetnames'][cfgkey]
        return cfgkey
    # Defaults to black.
    if cfgtype == 'datasetcolor':
        if cfgkey in mpqb_cfg['datasetcolors']:
            return mpqb_cfg['datasetcolors'][cfgkey]
        return mpqb_cfg['datasetcolors']['default']
    # Defaults to solid (provided as cfgkey)
    if cfgtype == 'linestyle':
        if cfgkey in mpqb_cfg['linestyles']:
            return mpqb_cfg['linestyles'][cfgkey]
        mpqb_cfg['linestyles']['default']
    # Defaults to 0.5 (provided as cfgkey)
    if cfgtype == 'linewidth':
        if cfgkey in mpqb_cfg['linewidths']:
            return mpqb_cfg['linewidths'][cfgkey]
        return mpqb_cfg['linewidths']['default']

    return None
