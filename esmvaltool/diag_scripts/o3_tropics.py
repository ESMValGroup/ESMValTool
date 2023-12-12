"""Diagnostics example for ozone variability."""

import logging
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import iris
from esmvalcore.preprocessor import (
    anomalies,
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

from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(Path(__file__).stem)

def get_provenance_record(attributes, ancestor_files):
    caption = attributes['caption'].format(**attributes)
    record = {
        'authors': ['gulletutan_gulcin'],
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_types': ['zonal'],
        'ancestors': ancestor_files,
    }
    return record

#List of variables per dataset will be processed
#proc_vars = {"BCC-CSM2-MR": ['o3'] }

def calculate_ozone(dict_item):
    """Calculation of ozone."""
     
    var = iris.load_cube(dict_item['filename'])
    var = var.collapsed('o3', iris.analysis.MEAN)
    return var

def plotting_support(cube, key, **kwargs):
    if cube.coords('time', dim_coords=True):
        ih.unify_time_coord(cube)
    iris.quickplot.plot(cube, label=key, **kwargs)
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.ylabel("Ozone anomalies over tropics")
    plt.xlabel("Time (years)")
    plt.title(f"Time series of monthly mean ozone  anomalies")

def plot_timeseries(dictionary, var, cfg):
    """Timeseries plot."""
    
    fig = plt.figure(figsize=(10, 4))
    sns.set_style('whitegrid')
    colors = plt.cm.viridis(np.linspace(0, 1, len(dictionary.keys())))
    baseplotname = f"Timeseries_ozone_anomalies"
    filename = get_plot_filename(baseplotname, cfg)

def main(cfg):

    """Compute the time average for each input dataset."""
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()
    # Demonstrate use of metadata access convenience functions.
    selection = select_metadata(input_data, short_name='tas', project='CMIP5')
    logger.info("Example of how to select only CMIP5 temperature data:\n%s"
                pformat(selection))
    selection = sorted_metadata(selection, sort='dataset')
    logger.info("Example of how to sort this selection by dataset:\n%s",
                pformat(selection))
    grouped_input_data = group_metadata(input_data,
                                        'variable_group',
                                        sort='dataset')
    logger.info(
        "Example of how to group and sort input data by variable groups from "
        "the recipe:\n%s", pformat(grouped_input_data))
    # Example of how to loop over variables/datasets in alphabetical order
    groups = group_metadata(input_data, 'variable_group', sort='dataset')
    for group_name in groups:
        logger.info("Processing variable %s", group_name)
        for attributes in groups[group_name]:
            logger.info("Processing dataset %s", attributes['dataset'])
            input_file = attributes['filename']
            cube = compute_diagnostic(input_file)
            
            output_basename = Path(input_file).stem
            if group_name != attributes['short_name']:
                output_basename = group_name + '_' + output_basename
            if "caption" not in attributes:
                attributes['caption'] = input_file
             provenance_record = get_provenance_record(
                 attributes, ancestor_files=[input_file])
             plot_diagnostic(cube, output_basename, provenance_record, cfg)

if __name__ == '__main__':
    with run_diagnostics() as config:
        main(config)

