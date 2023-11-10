"""Plot 2d profile."""
import logging
import cf_units
from copy import deepcopy
from pathlib import Path

import iris
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np

import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    get_plot_filename,
    group_metadata,
    run_diagnostic,
)

logger = logging.getLogger(Path(__file__).stem)

def get_plot_kwargs(cfg, dataset, **kwargs):
    """Get plot kwargs for different datasets."""
    plot_kwargs = dict(kwargs)
    key = dataset.get(cfg['title_key'])
    label= dataset.get(cfg['facet_used_for_labels'])
    if key is None:
        return plot_kwargs
    return plot_kwargs

def load_and_preprocess(dataset):
    """Load and preprocess data."""
    filename = dataset['filename']
    logger.info("Loading %s", filename)
    cube = iris.load_cube(filename)
    if cube.coords('air_pressure', dim_coords=True):
        z_coord = cube.coord('air_pressure', dim_coords=True)
        z_coord.attributes['positive'] = 'down'
        z_coord.convert_units('hPa')
    elif cube.coords('altitude', dim_coords=True):
        z_coord = cube.coord('altitude')
        z_coord.attributes['positive'] = 'up'
    if cube.var_name == 'pr' and cube.units == 'kg m-2 s-1':
        cube.units = 'mm s-1'
        cube.convert_units('mm day-1')
        dataset['units'] = 'mm day-1'
    elif cube.var_name == 'hus':
        cube.convert_units('g kg-1')
        
    iris.coord_categorisation.add_month_number(cube, 'time')
    cube_mean=cube.aggregated_by(['month_number'],iris.analysis.MEAN)
    cube_std=cube.aggregated_by(['month_number'],iris.analysis.STD_DEV)

    cube_std_neg= cube_mean-cube_std
    cube_std_pos= cube_mean+cube_std

    return cube_mean, cube_std_neg, cube_std_pos

def main(cfg):
    """Run diagnostic."""
    cfg = deepcopy(cfg)
    cfg.setdefault('title_key', 'dataset')
    logger.info("Using key '%s' to create titles for datasets",
                cfg['title_key'])

    # This diagnostic supports only a single variable, raise error otherwise
    input_data = list(cfg['input_data'].values())
#    all_vars = list(group_metadata(input_data, 'short_name'))
#    if len(all_vars) != 1:
#        raise ValueError(
#            f"Expected exactly 1 variable, got {len(all_vars):d}: {all_vars}")

    # Set fixed size for figure
    plt.figure(figsize=(8, 4), dpi=150)

    # Plot all datasets in single figure
    cubes_mean=iris.cube.CubeList()
    ancestors = []
    colors=['C0','C1','C2','C3','C4']
    i=0
    for dataset in input_data:
        ancestors.append(dataset['filename'])
        cube_mean, cube_std_neg, cube_std_pos = load_and_preprocess(dataset)

            # Plot hourly
        plot_kwargs = get_plot_kwargs(cfg, dataset)
        std_neg=np.array(cube_std_neg.data)
        std_pos=np.array(cube_std_pos.data)
        cubes_mean.append(cube_mean)
        plt.plot(cube_mean.coord("month_number").points,cube_mean.data, label=dataset[cfg['facet_used_for_labels']],color=colors[i],                 **plot_kwargs)
        plt.fill_between(cube_mean.coord("month_number").points,std_pos,std_neg, color=colors[i], alpha=0.2)
        i+=1

    # Plot appearance
    long_name = input_data[0]['long_name']
    short_name = input_data[0]['short_name']
    units = cube_mean.units
    plt.title(f" {long_name}")
    plt.xlabel("month number")
    plt.ylabel(f"{short_name} [{units}]")
    plt.legend()
    
    cube_ICON160=cubes_mean[0]
    cube_ICON80=cubes_mean[1]
    cube_ICON40=cubes_mean[2]
    
    cube_ref=cubes_mean[3]
 
    cube_ICON160_bias=cube_ICON160-cube_ref
    cube_ICON80_bias=cube_ICON80-cube_ref
    cube_ICON40_bias=cube_ICON40-cube_ref
    
    cube_ICON160_coll=cube_ICON160_bias.collapsed(["month_number"],iris.analysis.MEAN)
    cube_ICON160_rmse=cube_ICON160_bias.collapsed(["month_number"],iris.analysis.RMS)
    
    cube_ICON80_coll=cube_ICON80_bias.collapsed(["month_number"],iris.analysis.MEAN)
    cube_ICON80_rmse=cube_ICON80_bias.collapsed(["month_number"],iris.analysis.RMS)
    
    cube_ICON40_coll=cube_ICON40_bias.collapsed(["month_number"],iris.analysis.MEAN)
    cube_ICON40_rmse=cube_ICON40_bias.collapsed(["month_number"],iris.analysis.RMS)
    
    print(cube_ICON160_coll.data, cube_ICON160_rmse.data)
    print(cube_ICON80_coll.data, cube_ICON80_rmse.data)
    print(cube_ICON40_coll.data, cube_ICON40_rmse.data)     

    # Save plot
    plot_path = get_plot_filename(short_name, cfg)
    plt.savefig(plot_path, bbox_inches='tight', orientation='landscape')
    logger.info("Wrote %s", plot_path)
    plt.close()

    # Provenance tracking
    caption = (f"time series of {input_data[0]['long_name']} for various "
               f"datasets.")
    provenance_record = {
        'ancestors': ancestors,
        'authors': ['schlund_manuel'],
        'caption': caption,
        'plot_types': ['line'],
    }
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(plot_path, provenance_record)


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)