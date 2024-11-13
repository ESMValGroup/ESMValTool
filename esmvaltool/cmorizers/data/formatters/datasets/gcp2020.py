"""ESMValTool CMORizer for GCP2020 data.

Tier
    Tier 2: other freely-available dataset.

Source
    https://www.icos-cp.eu/science-and-impact/global-carbon-budget/2020

Last access
    20210908

Download and processing instructions
    Download the following file: '2020 Global Budget v1.0'
"""

from .gcp2018 import cmorization

# The following line makes it clear that the above import is not an error
cmorization
