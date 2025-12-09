import matplotlib.pyplot as plt
import numpy as np
import iris
import logging
from pathlib import Path

from esmvaltool.diag_scripts.shared import (
    run_diagnostic,
    save_figure,
)

logger = logging.getLogger(Path(__file__).stem)


def main(cfg):
    """Create a polar spiral plot of Arctic sea ice area over time."""
    # Read the time period from the config
    start_year = cfg["period"]["extract_time"]["start_year"]
    end_year = cfg["period"]["extract_time"]["end_year"]
    plot_title = f"Arctic Sea Ice Area (million km²)\n{start_year} to {end_year}"

    # Set the months correct for the number of years
    n_years = end_year - start_year + 1
    months = np.linspace(0, 2 * np.pi, 12, endpoint=False)  # Needed as linspace includes start AND stop by default
    all_months = []
    for i in range(n_years):
        all_months.extend(months + 2 * np.pi * i)

    # Read the data from the config
    input_data = cfg["input_data"].values()

    # Set up the figure
    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={'projection': 'polar'})
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)

    # Plot the data
    for dataset in input_data:
        cube = iris.load_cube(dataset["filename"])
        # print(cube)
        ax.plot(all_months, cube.data, alpha=0.5, label=dataset['alias'])

    # Tidy up the plot
    ax.set_xticks(months)
    ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
    ax.legend(loc='lower left', bbox_to_anchor=(0, 0), bbox_transform=fig.transFigure)
    ax.set_title(plot_title)

    # Save figure
    save_figure(
        'Arctic_sea_ice_area',
        {},
        cfg,
        figure=fig,
        close=True,
    )


if __name__ == "__main__":
    with run_diagnostic() as config:
        main(config)
