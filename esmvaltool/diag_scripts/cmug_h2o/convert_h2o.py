#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Functions to convert water vapour variables.

###############################################################################
convert_h2o.py
Author: Katja Weigel
EVal4CMIP + ESA-CMUG project
###############################################################################

Description
-----------
    Functions for:
        Convertion between vmr (H2O) [ppmV] and specific humidity (hus) [g/kg].
        Based on IDL routines from Dominik Brunner.

Configuration options
---------------------
    None

###############################################################################

"""


def h2o_to_hus(h2o):
    """Calculate H2O in specific humidity instead of vmr."""
    mda = 28.966  # molecular mass of dry air
    mwv = 18.016  # molecular mass of water vapour
    helpval = h2o * (mwv / mda)
    hus = helpval / (1.0 + helpval)

    return hus


def hus_to_h2o(hus):
    """Calculate H2O in vmr instead of specific humidity."""
    mda = 28.966  # molecular mass of dry air
    mwv = 18.016  # molecular mass of water vapour
    h2o = mda / mwv * hus / (1.0 - hus)

    return h2o
