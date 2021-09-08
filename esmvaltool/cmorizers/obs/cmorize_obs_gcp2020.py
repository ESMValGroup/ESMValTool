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

from .cmorize_obs_gcp2018 import cmorization


# This script is not supposed to be called directly; the main purpose of the
# following lines is to avoid FLAKE8 errors
if __name__ == '__main__':
    print(cmorization)
