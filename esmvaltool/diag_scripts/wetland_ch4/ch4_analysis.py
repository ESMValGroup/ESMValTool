#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 11:46:08 2022

@author: bergmant
"""

import logging
import os
import cartopy.crs as ccrs
import iris
import iris.coord_categorisation
import iris.plot as iplt
import matplotlib.pyplot as plt
import numpy as np
from esmvalcore.preprocessor import extract_region, extract_season

from esmvaltool.diag_scripts.shared import (
    get_diagnostic_filename,
    apply_supermeans,
    save_data,
    get_plot_filename,
    save_figure,
    get_control_exper_obs,
    group_metadata,
    run_diagnostic,
)
from esmvaltool.diag_scripts.shared._base import ProvenanceLogger
# from esmvaltool.diag_scripts.shared.plot import quickplot
# from esmvaltool.diag_scripts.shared.plot import global_contourf

logger = logging.getLogger(os.path.basename(__file__))


def scaling():
    """
    scale kg s-1 to g yr-1

    Returns
    -------
    scale_kg_per_m2_per_s_to_g_per_m2_per_year : float
        mass per m2 per s converted to g per m2 per year

    """
    scale_kg_per_s_to_g_per_year = 1E3 * (365 * 24 * 60 * 60)
    return scale_kg_per_s_to_g_per_year


def plot_map(axes, cube):
    """
    plot map of time mean of cube

    Parameters
    ----------
    axes : figure axes
        axes to plot on.
    cube : iris.cube
        cube with data to plot.

    Returns
    -------
    None.


    """
    #fig=plt.figure(figsize=(8,8))
    #projection=ccrs.NorthPolarStereo())
    #for index, yx_slice in enumerate(cube.slices(['latitude', 'longitude'])):
    #logger.info('Plotting year: %s',cube.coord('time'))
    #ax=fig.add_subplot(111,projection=ccrs.NorthPolarStereo())

    # Northern Hemisphere from 23 degrees north:
    axes.set_extent([-180, 180, 45, 90], crs=ccrs.PlateCarree())
    #newcube, extent = iris.analysis.cartography.project(yx_slice,ccrs.NorthPolarStereo(), nx=400, ny=200)
    #ax.set_global()
    #plt.subplot(1,cube.coord('time').shape[0],index+1)
    #print(repr(yx_slice))
    #logger.info('config quickplot %s',**cfg['quickplot'])
    #quickplot(newcube, **cfg['quickplot'])
    #global_contourf(yx_slice, cbar_range=[-1e-9,10e-9,7])
    scale = scaling()

    #con0=iplt.contourf(yx_slice*scale,cmap='hot_r',vmin=-1e-9, vmax=10e-9,levels=np.linspace(-1e-9,10e-9,11),extend='both')
    con0 = iplt.contourf(cube.collapsed(
        'time',iris.analysis.MEAN)*scale,
        cmap='hot_r',
        vmin=-1, vmax=20,
        levels = np.linspace(-1, 20, 22),
        extend='both'
        )
    axes.coastlines()
    gl = axes.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=2,
        color='gray',
        alpha=0.5,
        linestyle='--'
        )
    #cax0 = fig.add_axes([0.15, 0.1, 0.3, 0.05])
    cbar = plt.colorbar(
        con0,
        label=cube.long_name + ' (g m-2 yr-1)',
        orientation='horizontal'
        )
    return


def plot_timeseries(axes, cubes):
    """
    Plot area mean time series of cube.

    Parameters
    ----------
    axes : figure axes
        axes to plot on.
    cubes : iris.cube
        cube with data to plot.

    Returns
    -------
    lines : array
        handles for plotted lines.
    labels : list
        titles for plotted lines.


    """
    #plt.subplots(111)
    lines=[]
    labels=[]
    for model in cubes:
        logger.info(model)
        logger.info(type(model))
        lines.append(iplt.plot(model,label=model.attributes['source_id'],axes=axes)[0])
        logger.info(lines[0])
        labels.append(model.attributes['source_id'])
        #u#nit=model.attributes['units']
        #time=model.dimension['time']
        name = model.long_name
        #unit=model.units
    unit = 'kg'
    logger.info(lines)
    logger.info(labels)
    axes.xaxis.set_label_text('time')
    axes.yaxis.set_label_text(name + '\n[' + unit + ']')
    #axes.xaxis.set_label_text(time)
    return (lines,labels)
def plot_seasonal(axes,cubes):
    """
    Plot latitudinal seasonal mean across longitudes.

    Parameters
    ----------
    axes : TYPE
        DESCRIPTION.
    cubes : TYPE
        DESCRIPTION.

    Returns
    -------
    lines : array
        handles for plotted lines.
    labels : list
        titles for plotted lines.

    """
    seasonlen={}
    seasonlen['DJF'] = (31 + 31 + 28) * 3600 * 24
    seasonlen['MAM'] = (31 + 30 + 31) * 3600 * 24
    seasonlen['JJA'] = (30 + 31 + 31) * 3600 * 24
    seasonlen['SON'] = (30 + 31 + 30) * 3600 * 24
    seasonit = ['DJF', 'MAM', 'JJA', 'SON']
    lines = []
    labels = []
    meanseason_s = 3 * 30 * 24 * 3600
    unit = 'kg m-2'  #season.attributes['unit']
    for model in cubes:
        logger.info(model.collapsed('season', iris.analysis.MEAN))
        lines.append(iplt.plot(
            model.collapsed(
                'season',
                iris.analysis.MEAN),
            label=model.attributes['source_id'],
            axes=axes[0]))
        axes[0].xaxis.set_label_text('longitude')
        axes[0].set_xlim(0, 360)
        axes[0].yaxis.set_label_text(unit)
        name=model.long_name
        axes[0].set_title(name + ' ')
        for i,season in enumerate(model.slices('longitude')):
            #logger.info(model.slices('season'),i)
            #logger.info(season.season)
            logger.info(type(season))
            #lines.append(iplt.plot(season,label=season.attributes['source_id'],axes=axes[i+1]))
            iplt.plot(
                season * meanseason_s,
                label=season.attributes['source_id'],
                axes=axes[i + 1]
                )
            logger.info(lines[0])
            labels.append(model.attributes['source_id'])
            #u#nit=model.attributes['units']
            #time=model.dimension['time']
            name = season.long_name
            unit = 'kg m-2'#season.attributes['unit']
            axes[i + 1].xaxis.set_label_text('longitude')
            axes[i + 1].set_xlim(0,360)
            axes[i + 1].yaxis.set_label_text(unit)
            #if i!=0:
            #    axes[i+1].set_title(name+' in: yearly'
            #else:
            axes[i + 1].set_title(name + ' in: ' + seasonit[i])

    logger.info(lines)
    logger.info(labels)
    return (lines, labels)
def plot_scatter(axes, xcubes, ycubes):
    """
    Plot scatter plot with tas from xcube and ycube data.

    Parameters
    ----------
    axes : figure axes
        DESCRIPTION.
    xcubes : iris.cube
        cube with tas.
    ycubes : iris.cube
        cube data for scatter plot.

    Returns
    -------
    lines : list
        handles for plotted lines.
    labels : list
        titles for lines.

    """

    #plt.subplots(111)

    lines = []
    labels = []
    for model_x, model_y in zip(xcubes, ycubes):
        logger.info(model_x)
        logger.info(type(model_x))
        # convert kg s-1 to Tg
        emi_2_TG= 3600 * 24 * 30 * 1e-9
        logger.debug('change scling to actual months')
        lines.append(iplt.scatter(
            model_x,
            model_y*emi_2_TG,
            label=model_x.attributes['source_id'],
            axes=axes)
            )
        labels.append(model_x.attributes['source_id'])
        #u#nit=model.attributes['units']
        #time=model.dimension['time']
        name_x = model_x.long_name
        name_y = model_y.long_name
        unit_x = 'K'#model_x.attributes['unit']
        unit_y = 'Tg'#model_y.attributes['unit']
    if unit_y == 'kg m-2 s-1':
        unit_y = 'kg'
    axes.xaxis.set_label_text(name_x + ' (' + unit_x + ')')
    axes.yaxis.set_label_text(name_y + ' (' + unit_y + ')')
    #axes.xaxis.set_label_text(time)

    return (lines,labels)


def plot_diagnostics(data, cfg):
    """


    Parameters
    ----------
    data : dict of cubes
        dictionary of data cubes for plots.
    cfg : dict
        configuration.

    Returns
    -------
    path : STR
        path for plot.

    """
    #linestyle = {'linestyle':'-','linewidth':2}
    path=""
    if 'wetlandCH4_map' in data.keys():
        mapdata = data['wetlandCH4_map']
        logger.info(mapdata)
        nplots=len(mapdata)
        fig = plt.figure(figsize = (8,nplots * 8))
        #axes = fig.add_subplot(111,projection=ccrs.NorthPolarStereo())
        for i,model in enumerate(mapdata):
            axes = fig.add_subplot(
                nplots, 1, i + 1, projection=ccrs.NorthPolarStereo()
                )
            plot_map(axes,model)
            axes.set_title(model.attributes['source_id'])
        path = get_plot_filename('ch4_map', cfg)
    elif'wetlandCH4_timeseries' in data.keys():
        fig,axes = plt.subplots()
        plotdata = data['wetlandCH4_timeseries']
        logger.info(plotdata)
        logger.info(type(plotdata))
        lines, labels = plot_timeseries(axes, plotdata)
        legend = draw_legend(fig, lines, labels)
        #fig.legend(labels, handles=lines,
        #                  loc='upper left',
        #                  fontsize=8.)
        path = get_plot_filename('ch4_timeseries', cfg)

    elif 'wetlandCH4_scatter' in data.keys():
        logger.info(data.keys())
        fig,axes=plt.subplots()
        ydata=data['wetlandCH4_scatter']
        xdata=data['tas']
        logger.info(xdata)
        markers,labels = plot_scatter(axes, xdata, ydata)
        legend = draw_legend(fig, markers, labels)
        path = get_plot_filename('ch4_scatter', cfg)
    elif 'wetlandCH4_seasonal' in data.keys():
        logger.info(data.keys())
        fig,axes = plt.subplots(5, figsize=(12,24))
        logger.info(data)
        season = {
            12: 'DJF',
            1: 'DJF',
            2: 'DJF',
            3: 'MAM',
            4: 'MAM',
            5: 'MAM',
            6: 'JJA',
            7: 'JJA',
            8: 'JJA',
            9: 'SON',
            10: 'SON',
            11: 'SON'
        }
        plotdata=[]
        for cube in data['wetlandCH4_seasonal']:
            if cube.coords('month_number'):
                points = [
                    season[point] for point in cube.coord('month_number').points
                    ]
                cube.add_aux_coord(iris.coords.AuxCoord(
                    points, var_name='season'), cube.coord_dims('month_number')
                    )
                cube = cube.aggregated_by('season', iris.analysis.MEAN)
            logger.info(cube)
            plotdata.append(cube)
        logger.info(plotdata)
        #plotdata = data['wetlandCH4_seasonal']
        # logger.info(plotdata)
        # logger.info(type(plotdata))
        lines,labels = plot_seasonal(axes, plotdata)
        # legend = draw_legend(fig,lines,labels)
        # #fig.legend(labels, handles=lines,
        # #                  loc='upper left',
        # #                  fontsize=8.)
        path = get_plot_filename('ch4_seasonal', cfg)

    else:
        logger.info('No diagnostic yet for '+ data.keys())
    fig.savefig(path)
    return path


def draw_legend(fig, lines, labels):
    """Draw the legend."""
    return fig.legend(lines,
                      labels,
                      loc='upper left',
                      fontsize=6.,
                      bbox_to_anchor=(.81, .92))

def plot_diagnostic(cube, basename, provenance_record, cfg,variable_group):
    """Create diagnostic data and plot it."""

    # Save the data used for the plot
    save_data(basename, provenance_record, cfg, cube)
    time = cube.coord('time')
    logger.info('times %s',time)
    logger.info('times %s',time[0])
    logger.info(cfg.get('quickplot'))
    #logger.info(cfg['variable_group'])
    logger.info(cfg)
    #logger.info(cfg['variables'])
    if cfg.get('quickplot').get('plot_type')=='polar' and variable_group=='wetlandCH4_map':
        logger.info(cfg.get('quickplot'))
        logger.info(cfg)
        logger.info(str(cube.coord('time')))
        # Create the plot
        logger.info(cube.coord('time'))
        logger.info(cube.coord('time').shape)
        #plt.subplots(nrows=1,ncols=cube.coord('time').shape[0],projection=ccrs.NorthPolarStereo())
        plot_map(cube)
# =============================================================================
#         fig=plt.figure(figsize=(8,8))
#         #projection=ccrs.NorthPolarStereo())
#         for index, yx_slice in enumerate(cube.slices(['latitude', 'longitude'])):
#             logger.info('Plotting year: %s',cube.coord('time'))
#             ax=fig.add_subplot(111,projection=ccrs.NorthPolarStereo())
#             # Northern Hemisphere from 23 degrees north:
#             ax.set_extent([-180, 180, 45, 90], crs=ccrs.PlateCarree())
#             #newcube, extent = iris.analysis.cartography.project(yx_slice,ccrs.NorthPolarStereo(), nx=400, ny=200)
#             #ax.set_global()
#             #plt.subplot(1,cube.coord('time').shape[0],index+1)
#             print(repr(yx_slice))
#             #logger.info('config quickplot %s',**cfg['quickplot'])
#             #quickplot(newcube, **cfg['quickplot'])
#             #global_contourf(yx_slice, cbar_range=[-1e-9,10e-9,7])
#             scale=scaling()
#
#             #con0=iplt.contourf(yx_slice*scale,cmap='hot_r',vmin=-1e-9, vmax=10e-9,levels=np.linspace(-1e-9,10e-9,11),extend='both')
#             con0=iplt.contourf(yx_slice*scale,cmap='hot_r',vmin=-1, vmax=20,levels=np.linspace(-1,20,22),extend='both')
#             ax.coastlines()
#             #cax0 = fig.add_axes([0.15, 0.1, 0.3, 0.05])
#             cbar = plt.colorbar(con0,  label='g m-2 year-1',
#                                 orientation='horizontal')
# =============================================================================
    elif cfg.get('quickplot').get('plot_type')=='times' and variable_group=='wetlandCH4_timeseries':
        logger.info('timeseries plot')
        plot_timeseries(cube)
        #iplt.plot(cube)
        # And save the plot
    else:
        logger.info('no plot')
    save_figure(basename, provenance_record, cfg)
def load_data(config):
    """Load cubes into config dict."""
    for key in config['input_data'].keys():
        filename = config['input_data'][key]['filename']
        config['input_data'][key]['cube'] = iris.load_cube(filename)


def load_cube(filename):
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)

    #logger.debug("Running example computation")
    #cube = iris.util.squeeze(cube)
    return cube


def get_provenance_record( ancestor_files, caption):
    """Create a provenance record describing the diagnostic data and plot."""
    #caption = ("Multiyear mean of emission of wetland CH4 between 2001 and 2005 "
    #           ". (b) monthly mean emissions of wetland CH4.")

    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['nhpolar'],
        'plot_types': ['polar','times'],
        'authors': [
            'tommi_bergman',
        ],
        'references': [
            'CRESCENDO',
        ],
        'ancestors': ancestor_files,
    }
    return record

def setup_figure():
    pass

def prepare_data(cfg):
    """

    Parameters
    ----------
    cfg : TYPE
        DESCRIPTION.

    Returns
    -------
    data : list of cubes
        get data cubes from input.

    """
    groups = group_metadata(config['input_data'].values(), 'variable_group')
    logger.info(groups.keys())
    #wetlanch4_map=groups['wetlandCH4_map']
    #wetlanch4_ts=groups['wetlandCH4_timeseries']
    data = {}
    for i in groups:
        logger.info(groups[i])
        logger.info(type(groups[i]))
        data[i] = [ds['cube'] for ds in groups[i]]
    return data


def write_data(conf, data):
    """Write all the calculated data to output file."""
    logger.info(data)
    if 'wetlandCH4_map' in data.keys():
        varid='wetlandCH4_map'
    elif 'wetlandCH4_timeseries' in data.keys():
        varid='wetlandCH4_timeseries'
    elif 'wetlandCH4_scatter' in data.keys():
        varid='wetlandCH4_scatter'
    elif 'wetlandCH4_seasonal' in data.keys():
        varid='wetlandCH4_seasonal'
    else:
        logger.info('Unknown variable in the cube (write_data).')
    cubes = iris.cube.CubeList(data[varid])
    #logger.info(cubes[0])
    path = get_diagnostic_filename(varid, conf)
    iris.save(cubes, path)
    return path


def create_caption(cases):
    for case in cases:
        logger.info(case)
        if case == 'methane_diagnostic_map':
            caption = 'Multiyear mean of polar [NN-NN] wetland emissions of CH4.'
        elif 'methane_diagnostic_timeseries' in case:
            caption = 'Monthly mean  polar [NN-NN] wetland emissions of CH4.'
        elif 'methane_diagnostic_temperature' in case:
            caption = 'Monthly mean polar [NN-NN] wetland'
            ' emissions of CH4 dependency on temperature.'
        elif 'methane_diagnostic_seasonal' in case:
            caption = 'Seasonal mean wetland emissions of CH4.'
        else:
            caption = 'Caption not available for this case: '+ case
    return caption


def main(cfg):
    """Execute  analysis and plotting.



    """
    logger.setLevel(cfg['log_level'].upper())
    #input_data, grouped_input_data = do_preamble(cfg)
    load_data(cfg)
    data = prepare_data(cfg)
    plot_path = plot_diagnostics(data, cfg)
    #input_data=cfg['input_data'].values()
    #my_files_dict=group_metadata(input_data,sor'dataset')
    my_files_dict = group_metadata(cfg['input_data'].values(),'diagnostic')
    #how to get diagnostic name (variable_group)
    #logger.info(cfg['variable_group'])
    logger.info(my_files_dict)
    caption = create_caption(my_files_dict.keys())
    ancestor_files = list(config['input_data'].keys())
    provenance_record = get_provenance_record(ancestor_files,caption)
    netcdf_path = write_data(cfg, data)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)
        provenance_logger.log(plot_path, provenance_record)

    # groups=group_metadata(input_data,'variable_group',sort='dataset')
    # for group_name in groups:
    #     logger.info('Processing variable %s',group_name)
    #     for attributes in groups[group_name]:
    #         logger.info('processing dataset %s',attributes['dataset'])
    #         input_file=attributes['filename']
    #         data_cube=load_cube(input_file)

    #         output_basename = Path(input_file).stem
    #         if group_name != attributes['short_name']:
    #             output_basename = group_name + '_' + output_basename
    #         provenance_record = get_provenance_record(
    #             attributes, ancestor_files=[input_file])
    #         plot_diagnostic(data_cube, output_basename, provenance_record, cfg,group_name)
if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)