"""Plot maps."""
import logging
from copy import deepcopy
from pathlib import Path
from pprint import pformat

import cartopy.crs as ccrs
import iris
import matplotlib.colors as colors
import matplotlib.pyplot as plt

from esmvaltool.diag_scripts.shared import (
    get_plot_filename,
    group_metadata,
    plot,
    run_diagnostic,
    select_metadata,
)

logger = logging.getLogger(Path(__file__).stem)


def get_reference_dataset(cfg, input_data):
    """Extract reference dataset."""
    if 'reference' in cfg:
        ref_kwargs = cfg['reference']
        logger.info("Using kwargs %s to extract reference dataset", ref_kwargs)
        ref_datasets = select_metadata(input_data, **ref_kwargs)
        if len(ref_datasets) != 1:
            raise ValueError(
                f"Expected exactly one reference dataset with kwargs "
                f"{ref_kwargs}, found {len(ref_datasets):d}:\n"
                f"{pformat(ref_datasets)}")
        ref_cube = load_and_preprocess(ref_datasets[0])
        return (ref_datasets[0], ref_cube)
    return (None, None)


def load_and_preprocess(dataset):
    """Load and preprocess data."""
    filename = dataset['filename']
    logger.info("Loading %s", filename)
    cube = iris.load_cube(filename)

    # Tmporal mean
    if cube.coords('time', dim_coords=True):
        cube = cube.collapsed('time', iris.analysis.MEAN)

    # Convert units of some variables
    if cube.var_name == 'tas':
        cube.convert_units('celsius')
    if cube.var_name == 'pr':
        cube.units = 'mm s-1'
        cube.convert_units('mm day-1')
    return cube


def plot_dataset_with_ref(cfg, dataset, ref_dataset, ref_cube):
    """Plot single dataset with a reference dataset."""
    logger.info("Plotting %s with reference %s", dataset[cfg['title_key']],
                ref_dataset[cfg['title_key']])
    cube = load_and_preprocess(dataset)

    # Create single figure with multiple axes
    (fig, axes_list) = plt.subplots(
        2, 2, subplot_kw={'projection': ccrs.Robinson(central_longitude=10)})

    # Plot absolute values of dataset and reference dataset
    plot_kwargs = cfg.get('plot_kwargs_abs', {})
    map_plot = iris.plot.contourf(cube, axes=axes_list[0][0], **plot_kwargs)
    iris.plot.contourf(ref_cube, axes=axes_list[0][1], **plot_kwargs)
    colorbar = plt.colorbar(map_plot, ax=list(axes_list[0]),
                            orientation='horizontal', aspect=50)
    fontsize = 8
    colorbar.set_label(f"{cube.var_name} [{cube.units}]", labelpad=0.0,
                       fontsize=fontsize)
    colorbar.ax.tick_params(labelsize=fontsize)

    # Plot bias
    bias_cube = cube - ref_cube
    plot_kwargs = dict(norm=colors.CenteredNorm(), cmap='bwr')
    plot_kwargs.update(cfg.get('plot_kwargs_bias', {}))
    map_plot = iris.plot.contourf(bias_cube, axes=axes_list[1][0],
                                  **plot_kwargs)
    colorbar = plt.colorbar(map_plot, ax=axes_list[1][0],
                            orientation='horizontal', aspect=30)
    colorbar.set_label(rf"$\Delta${cube.var_name} [{cube.units}]",
                       labelpad=0.0, fontsize=fontsize)
    colorbar.ax.tick_params(labelsize=fontsize)

    # Plot appearance
    title_kwargs = dict(pad=3.0, fontdict=dict(fontsize=10))
    title_key = cfg['title_key']
    axes_list[0][0].set_title(dataset[title_key], **title_kwargs)
    axes_list[0][1].set_title(ref_dataset[title_key], **title_kwargs)
    axes_list[1][0].set_title(
        f"{dataset[cfg['title_key']]} - {ref_dataset[cfg['title_key']]}",
        **title_kwargs)
    for axes in list(axes_list.ravel()):
        axes.gridlines(color='lightgrey', alpha=0.5)
        axes.coastlines()
        axes.set_global()
    axes_list[1][1].set_visible(False)

    # Return figure and basename for plot
    tag = dataset[cfg['title_key']].replace(' ', '_')
    tag_ref = ref_dataset[cfg['title_key']].replace(' ', '_')
    return (fig, f"{tag}_vs_{tag_ref}")


def plot_dataset_without_ref(cfg, dataset):
    """Plot single dataset without using a reference dataset."""
    logger.info("Plotting %s", dataset[cfg['title_key']])
    cube = load_and_preprocess(dataset)
    plot_kwargs = dict(levels=10, cbar_label=f"{cube.var_name} [{cube.units}]")
    plot.global_contourf(cube, **plot_kwargs)
    plt.title(dataset[cfg['title_key']])
    basename = dataset[cfg['title_key']].replace(' ', '_')
    return (plt.gcf(), basename)


def main(cfg):
    """Run diagnostic."""
    cfg = deepcopy(cfg)
    cfg.setdefault('title_key', 'dataset')
    logger.info("Using key '%s' to create titles for datasets",
                cfg['title_key'])

    # This diagnostic supports only a single variable, raise error otherwise
    input_data = list(cfg['input_data'].values())
    all_vars = list(group_metadata(input_data, 'short_name'))
    if len(all_vars) != 1:
        raise ValueError(
            f"Expected exactly 1 variable, got {len(all_vars):d}: {all_vars}")

    # Load reference dataset if possible
    (ref_dataset, ref_cube) = get_reference_dataset(cfg, input_data)

    # Create individual plot for each dataset
    for dataset in input_data:

        # Do not create plot for ref_dataset
        if dataset == ref_dataset:
            continue

        # Plot
        if ref_dataset is None:
            (fig, basename) = plot_dataset_without_ref(cfg, dataset)
        else:
            (fig, basename) = plot_dataset_with_ref(cfg, dataset, ref_dataset,
                                                    ref_cube)

        # General plot appearance
        if 'suptitle' in cfg:
            fig.suptitle(cfg['suptitle'].format(**dataset))
        else:
            fig.suptitle(f"{dataset['long_name']} ({dataset['start_year']}-"
                         f"{dataset['end_year']})")

        # Save plot
        plot_path = get_plot_filename(basename, cfg)
        plt.savefig(plot_path, bbox_inches='tight', orientation='landscape',
                    dpi=200)
        logger.info("Wrote %s", plot_path)
        plt.close()


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
