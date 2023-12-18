"""
ESMValTool CMORizer for CFSR data.

Tier
    Tier 2: other freely-available dataset.

Source
    Research Data Archive (RDA):
    https://rda.ucar.edu/datasets/ds093.2/

Last access
    20230411

Download and processing instructions
    see download script cmorizers/data/downloaders/datasets/cfsr.py
"""

from .cfsv2 import cmorization

# The following line makes it clear that the above import is not an error
cmorization
