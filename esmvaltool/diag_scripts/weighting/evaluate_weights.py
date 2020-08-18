"""
Implementation of the climwip weighting scheme

Lukas Brunner et al. section 2.4
https://iopscience.iop.org/article/10.1088/1748-9326/ab492f
"""
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import xarray as xr
from scipy.spatial.distance import pdist, squareform
import pandas as pd
from pathlib import Path

from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            get_plot_filename, run_diagnostic,
                                            select_metadata)


def read_csv(filename: str):
    """Read a `.csv` file into a Pandas Dataframe."""
    return pd.read_csv(filename)


def main(cfg):
    # TODO: Finish this script

    input_files = cfg['input_files']
    filename = cfg['weights']

    csv_path = Path(input_files[0]) / filename

    df = read_csv(csv_path)



if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)

