# Diagnostics script for ozone variability over tropics

import logging
from pathlib import Path
import matplotlib.pyplot as plt
import iris
from esmvalcore.preprocessor import (
    anomalies,
    extract_region,
    extract_levels,
    area_statistics,
)
from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostics,
    save_data,
    save_figure,
    select_metadata,
)
from esmvaltool.diag_scripts.shared._base import (
    get_plot_filename, )

logger = logging.getLogger(Path(__file__).stem)


def main(cfg):
    for name, dataset in cfg['input_data'].items():
        logger.info(f"Processing dataset: {name}")

        #Loading data
        cube = iris.load_cube(dataset['filename'])

        #Calculation of ozone anomalies
        ozone_anomalies = calculate_ozone(cube)

        #Plotting
        fig, ax = plt.subplots()
        time_points = cube.coord('time').points
        plt.plot(time_points, ozone_anomalies.data)
        ax.set_xlabel('Time')
        ax.set_ylabel('Ozone Anomalies')

        save_figure(fig, "o3_figure")

if __name__ == '__main__':
    with run_diagnostics() as config:
        main(config)
