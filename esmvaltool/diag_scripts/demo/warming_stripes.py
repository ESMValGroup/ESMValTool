"""Python script that visualizes a timeseries of global temperature anomalies
in the form of the popular warming stripes figure by Ed Hawkins.
"""

from pathlib import Path

import iris
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared._base import get_plot_filename


def plot_warming_stripes(cube, name, cmap):
    """Plot timeseries of global temperature as warming stripes.

    cube: iris cube of temperature anomalies with only time dimension.
    """
    temperature = cube.data.data
    nx = len(temperature)
    x = np.arange(len(temperature))
    y = np.array([0, 1])
    temperature = np.vstack([temperature, temperature])

    figure, axes = plt.subplots()
    axes.pcolormesh(x, y, temperature, cmap=cmap, shading='auto')
    axes.set_axis_off()
    axes.set_title(f'Warming stripes for {name}')
    return figure


def main(cfg):
    """Compute the time average for each input dataset."""
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    # Example of how to loop over variables/datasets in alphabetical order
    for dataset in input_data:
        # Load the data
        input_file = dataset['filename']
        name = dataset['variable_group']
        cube = iris.load_cube(input_file)

        # Plot as warming stripes
        cmap = cfg['colormap']
        fig = plot_warming_stripes(cube, name, cmap)

        # Save output
        output_file = Path(input_file).stem.replace('tas', name)
        output_path = get_plot_filename(output_file, cfg)
        fig.savefig(output_path)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)

