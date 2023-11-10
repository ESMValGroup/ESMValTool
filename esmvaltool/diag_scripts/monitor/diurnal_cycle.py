"""Plot 2d profile."""
import logging
import cf_units
from copy import deepcopy
from pathlib import Path

import iris
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import datetime

import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvalcore.iris_helpers import add_leading_dim_to_cube, date2num
from cf_units import Unit

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
        
#    iris.coord_categorisation.add_hour(cube, 'time')
#    cube_mean=cube.aggregated_by(['hour'],iris.analysis.MEAN)
#    cube_std=cube.aggregated_by(['hour'],iris.analysis.STD_DEV)

#    cube_std_neg= cube_mean-cube_std
#    cube_std_pos= cube_mean+cube_std

    return cube

def localTimeApprox(cube):
    
   """Returns local hour approximation"""
   cubes_mean=iris.cube.CubeList()
   cubes_std_pos=iris.cube.CubeList()
   cubes_std_neg=iris.cube.CubeList()

   for lon in cube.coord("longitude").points:
        
        cube_lon=cube.extract(iris.Constraint(longitude=lon))
        time_coord=cube_lon.coord("time")

        local_datetime = time_coord.units.num2date(time_coord.points) + datetime.timedelta(hours=((lon + 180) % 360 - 180 ) / 180 * 12)
        t_units = Unit(
            'days since 1850-01-01', calendar='proleptic_gregorian'
         )
        local_dt_points = date2num(np.array(local_datetime), t_units)
  
        # Modify time coordinate in place
        time_coord.points = local_dt_points
        time_coord.units = t_units

        cube_lon.replace_coord(time_coord)
        if cube_lon.coords("hour"):
            cube_lon.remove_coord("hour")
        iris.coord_categorisation.add_hour(cube_lon, 'time')
        cube_lon_mean=cube_lon.aggregated_by(['hour'],iris.analysis.MEAN)
        cube_lon_std=cube_lon.aggregated_by(['hour'],iris.analysis.STD_DEV)
        
# Extract the coordinate values and corresponding data
        coord_values = cube_lon_mean.coord('hour').points
        data_mean = cube_lon_mean.data
        data_std = cube_lon_std.data
# Sort the coordinate values and the data
        sort_indices = np.argsort(coord_values)
        sorted_coord_values = coord_values[sort_indices]
        sorted_data_mean = data_mean[sort_indices]
        sorted_data_std = data_std[sort_indices]
        
        if sorted_coord_values[0]==1:
            sorted_coord_values=sorted_coord_values-1
        if sorted_coord_values[0]==2:
            sorted_coord_values=sorted_coord_values+1
        for i in range(len(sorted_coord_values)):
            if sorted_coord_values[i] == 24:
                sorted_coord_values[i] = 0
                
                
        sort_indices = np.argsort(sorted_coord_values)
        sorted_coord_values = sorted_coord_values[sort_indices]
        sorted_data_mean = data_mean[sort_indices]
        sorted_data_std = data_std[sort_indices]

# Create a new cube with the sorted values
        sorted_cube_lon_mean = cube_lon_mean.copy(data=sorted_data_mean)
        sorted_cube_lon_std = cube_lon_std.copy(data=sorted_data_std)
        sorted_cube_lon_mean.coord("hour").points = sorted_coord_values
        sorted_cube_lon_std.coord("hour").points = sorted_coord_values
        
        cubes_mean.append(sorted_cube_lon_mean)
        cubes_std_pos.append(sorted_cube_lon_mean+sorted_cube_lon_std)
        cubes_std_neg.append(sorted_cube_lon_mean-sorted_cube_lon_std)
 
   cube_mean=cubes_mean.merge_cube()
   cube_std_pos=cubes_std_pos.merge_cube()
   cube_std_neg=cubes_std_neg.merge_cube()
 
   grid_areas_mean = iris.analysis.cartography.area_weights(cube_mean)
   grid_areas_std_pos = iris.analysis.cartography.area_weights(cube_std_pos)
   grid_areas_std_neg = iris.analysis.cartography.area_weights(cube_std_neg)

   cube_mean = cube_mean.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_mean)
   cube_std_pos = cube_std_pos.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_std_pos)
   cube_std_neg = cube_std_neg.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_std_neg)
   for i in range(len(cube_std_pos.data)):
        if cube_std_pos.data[i] < 0:
            cube_std_pos.data[i] = 0
        if cube_std_neg.data[i] < 0:
            cube_std_neg.data[i] = 0     
   return cube_mean, cube_std_pos, cube_std_neg

def main(cfg):
    """Run diagnostic."""
    cfg = deepcopy(cfg)
    cfg.setdefault('title_key', 'dataset')
    logger.info("Using key '%s' to create titles for datasets",
                cfg['title_key'])

    # This diagnostic supports only a single variable, raise error otherwise
    input_data = list(cfg['input_data'].values())
    all_vars = list(group_metadata(input_data, 'short_name'))
#    if len(all_vars) != 1:
#        raise ValueError(
#            f"Expected exactly 1 variable, got {len(all_vars):d}: {all_vars}")

    # Set fixed size for figure
    plt.figure(figsize=(8, 4), dpi=150)

    # Plot all datasets in single figure
    ancestors = []
    colors=['C0','C1','C2','C3','C4']
    i=0
    cubes_mean=iris.cube.CubeList()
    for dataset in input_data:
        ancestors.append(dataset['filename'])
        cube = load_and_preprocess(dataset)
        cube_mean, cube_std_pos, cube_std_neg = localTimeApprox(cube)
        cubes_mean.append(cube_mean)
        # Plot hourly
        plot_kwargs = get_plot_kwargs(cfg, dataset)
        plt.plot(cube_mean.coord("hour").points,cube_mean.data, label=dataset[cfg['facet_used_for_labels']],color=colors[i], **plot_kwargs)
        plt.fill_between(cube_mean.coord("hour").points,cube_std_pos.data,cube_std_neg.data, color=colors[i], alpha=0.1)
        i+=1
    
    
    
    # Plot appearance
    long_name = input_data[0]['long_name']
    short_name = input_data[0]['short_name']
    units = cube.units
    
    cube_ICON160=cubes_mean[0]
    cube_ICON80=cubes_mean[1]
    cube_ICON40=cubes_mean[2]
    
    cube_ref=cubes_mean[3]
 
    cube_ICON160_bias=cube_ICON160-cube_ref
    cube_ICON80_bias=cube_ICON80-cube_ref
    cube_ICON40_bias=cube_ICON40-cube_ref
    
    cube_ICON160_coll=cube_ICON160_bias.collapsed(["hour"],iris.analysis.MEAN)
    cube_ICON160_rmse=cube_ICON160_bias.collapsed(["hour"],iris.analysis.RMS)
    
    cube_ICON80_coll=cube_ICON80_bias.collapsed(["hour"],iris.analysis.MEAN)
    cube_ICON80_rmse=cube_ICON80_bias.collapsed(["hour"],iris.analysis.RMS)
    
    cube_ICON40_coll=cube_ICON40_bias.collapsed(["hour"],iris.analysis.MEAN)
    cube_ICON40_rmse=cube_ICON40_bias.collapsed(["hour"],iris.analysis.RMS)
    
    print(cube_ICON160_coll.data, cube_ICON160_rmse.data)
    print(cube_ICON80_coll.data, cube_ICON80_rmse.data)
    print(cube_ICON40_coll.data, cube_ICON40_rmse.data)
    
#    for i in points:
#        cube_ref_slice=cube_ref.extract("hour"=i)
#        cube_ICON160_slice=cube_ICON160.extract("hour"=i)
#        cube_ICON80_slice=cube_ICON80.extract("hour"=i)
#        cube_ICON40_slice=cube_ICON40.extract("hour"=i)
        
#        cube_ICON160_bias=cube_ICON160_slice-cube_ref_slice
#        cube_ICON80_bias=cube_ICON80_slice-cube_ref_slice
#        cube_ICON40_bias=cube_ICON40_slice-cube_ref_slice

    plt.title(f" {long_name}")
    plt.xlabel("hour")
    plt.ylabel(f"{short_name} [{units}]")
    plt.legend()

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