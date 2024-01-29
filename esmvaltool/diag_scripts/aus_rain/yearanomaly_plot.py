"""based on warming_stripes.py

"""

from pathlib import Path

import iris
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import dask
import dask.array as da
import os
import logging
from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared._base import get_plot_filename


# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))

def plot_time_anomaly(ds,colour,longname,outname):
    """
    cube: iris cube of precipitaion anomalies with only time dimension.
    change to xarray
    """
    ## resample time to yearly? -calc means?

    values = ds['pr'].values
    # x = np.arange(len(values))
    # y = np.array([0, 1])
    # values = np.vstack([values, values])
    time = ds['time'].values

    figure, axes = plt.subplots()
    # axes.pcolormesh(x, y, values, cmap=cmap, shading='auto')
    # axes.set_axis_off()
    # val_name = cube.long_name
    axes.set_title(f'{longname} anomaly. {outname}')
    
    axes.bar(time, values, width=2, edgecolor=colour, linewidth=2)
    # axes.spines['bottom'].set_position('center')
    axes.axhline(0,color='black',lw=1)

    # axes.set_xticklabels([])
    ary = rollingwindow(ds)
    axes.plot(ary['time'].values, ary.values,color='grey')

    return figure

def compute_anom(data_dict): ##full and ref files
    ## from _compute_anomalies() in esmvalcore/preprocessor/_time.py ## doesn't work..
    ds1 =  xr.open_dataset(data_dict['full']) #iris.load_cube(data_dict['full'])
    ref_mean = xr.open_dataset(data_dict['ref'])['pr'].values.mean() #iris.load_cube(data_dict['ref']) ##calc mean of array values

    # cube calc difference of value from mean #or other way around?
    ds1['pr'] = ref_mean - ds1['pr']

    # cube = cube.copy(data)
    # cube.remove_coord(cube_coord)
    return ds1

def rollingwindow(dataset):

    # dataarray from dataset? (time/year and summed pr)?
    r = dataset['pr'].rolling(time=11)
    ary = r.mean()

    return ary ##array for adding to plot?

def main(cfg):
    """Compute the time average for each input dataset."""
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    datadict = {}
    # logger.info(f"!!! input_data: {input_data}..")
    # Example of how to loop over variables/datasets in alphabetical order
    for dataset in input_data:
        # Load the data
        input_file = dataset['filename']
        name = dataset['variable_group']
        output_file = dataset['diagnostic']        
        longname = dataset['long_name']
        # cube = iris.load_cube(input_file)
        datadict[name.split('_')[-1]] = input_file


    ## run anomaly with 2 cubes #fnc
    cube = compute_anom(datadict)
    # ds1 =  xr.open_dataset(datadict['full'])
    # ary = rollingwindow(ds1)
    
    colour = cfg['colour']
    fig = plot_time_anomaly(cube, colour,longname,output_file)

    # Save output
    # output_file = input_data['diagnostic'] #Path(input_file).stem#.replace('pr', name)
    output_path = get_plot_filename(output_file, cfg)
    fig.savefig(output_path) 

if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
