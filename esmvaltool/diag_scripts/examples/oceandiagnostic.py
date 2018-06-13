"""Python example diagnostic."""
import inspect
import logging
import os
import sys

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import yaml

from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def folder(name):
    """
        This snippet takes a string, makes the folder and the string.
        It also accepts lists of strings.
        """
    if isinstance(name, list):
        name = '/'.join(name)
    if name[-1] != '/':
        name = name + '/'
    if os.path.exists(name) is False:
        os.makedirs(name)
        logger.info('Making new directory:', name)
    return name


def get_cfg():
    """Read diagnostic script configuration from settings.yml."""
    settings_file = sys.argv[1]
    with open(settings_file) as file:
        cfg = yaml.safe_load(file)
    return cfg, settings_file


def get_input_files(cfg, index=0):
    """Get a dictionary with input files from metadata.yml files."""
    metadata_file = cfg['input_files'][index]
    with open(metadata_file) as file:
        metadata = yaml.safe_load(file)
    return metadata


def sensibleUnits(cube,name):
    """
    Convert the cubes into some friendlier units for the range of 
    values typically seen in BGC.
    """
    
    new_units = ''
    
    if name in ['tos','thetao']:
        new_units = 'celsius'
        
    if name in ['no3',]:
        new_units = 'mmol m-3'

    if name in ['chl',]:
        new_units = 'mg m-3'
    
    if new_units != '':
        logger.info(' '.join(["Changing units from ",str(cube.units),'to', new_units]))
        cube.convert_units(new_units)
        
    return cube
        

def get_image_path(cfg,
                   md,
                   prefix='',
                   suffix='',
                   image_extention='png',
                   basenamelist=[
                       'project', 'model', 'mip', 'exp', 'ensemble', 'field',
                       'short_name', 'start_year', 'end_year'
                   ]):
    """
        This produces a path to the final location of the image.
        The cfg is the opened global config,
        md is the metadata dictionairy (for the individual model file)
        """

    path = folder(cfg['plot_dir'])
    if prefix:
        path += prefix + '_'
    path += '_'.join([str(md[b]) for b in basenamelist])
    if suffix:
        path += '_' + suffix
    path += '.' + image_extention
    logger.info("Image path will be: %s", path)
    return path

def make_cube_layer_dict(cube):
    """
        This method takes a cube and return a dictionairy with a cube for each layer 
        as it's own item. ie:
          cubes[depth] = cube from specific layer
        Also, cubes with no depth component are returns as:
          cubes[''] = cube with no depth component.
        """
        
    # Check layering:
    depth = cube.coords('depth')
    cubes = {}
   
    if depth == []:
        cubes[''] = cube
    else: 
        # iris stores coords as a list with one entry:
        depth = depth[0]
        if len(depth.points) in [1,]:
            cubes[''] = cube
        else:
            for l,layer in enumerate(depth.points):
                cubes[layer] = cube[:,l]
    return cubes
    
def TimeSeriesPlots(cfg, md, fn,):
    """
        This function makes a simply plot for an indivudual model.
        The cfg is the opened global config,
        md is the metadata dictionairy (for the individual model file).
        Mean: the temporal mean of the entire time range
        """
    # Load cube and set up units
    cube = iris.load_cube(fn)  
    cube = sensibleUnits(cube, md['short_name'])

    # Is this data is a multi-model dataset?
    multiModel = md['model'].find('MultiModel') > -1

    # Collapse cube into a time series.
    weights = iris.analysis.cartography.area_weights(cube, normalize=False)
    coords = ['longitude', 'latitude']
    cube = cube.collapsed(
        coords,
        iris.analysis.MEAN,
        weights=weights,
    )
    
    # Make a dict of cubes for each layer.
    cubes = make_cube_layer_dict(cube)
                
    # Making plots for each layer
    for l,(layer,c) in enumerate(cubes.items()):
        layer = str(layer)
        
        if multiModel:
            qplt.plot(c, label=md['model'], ls=':')
        else:
            qplt.plot(c, label=md['model'])

        # Add title, legend to plots
        title = ' '.join([md['model'], md['long_name']])
        if layer:
            title = ' '.join([title,'(',layer,str(c.coords('depth')[0].units),')'])
        plt.title(title)
        plt.legend(loc='best')
    
        # Determine filename:
        if multiModel:
            path = folder(cfg['plot_dir']) + os.path.basename(fn).replace(
                '.nc', '_timeseries_'+str(l)+'.png')
        else:
            path = get_image_path(
                cfg,
                md,
                suffix='timeseries_'+str(l),
                image_extention='png',
            )
            
        # Saving files:
        if cfg['write_plots']:

            logger.info('Saving plots to %s', path)
            plt.savefig(path)

        plt.close()


def MapOfMeanOverTime(cfg, md, fn):
    """
        This function makes a simply plot for an indivudual model.
        The cfg is the opened global config,
        md is the metadata dictionairy (for the individual model file).
        """
    # Load cube and set up units
    cube = iris.load_cube(fn)  
    cube = sensibleUnits(cube, md['short_name'])
    
    # take mean over entire time range
    cube = cube.collapsed('time', iris.analysis.MEAN)

    # Is this data is a multi-model dataset? 
    multiModel = md['model'].find('MultiModel') > -1


    # Make a dict of cubes for each layer.
    cubes = make_cube_layer_dict(cube)
                
    # Making plots for each layer
    for l,(layer,c) in enumerate(cubes.items()):
        layer = str(layer)

        # Making plots
        qplt.contourf(c, 25)
        plt.gca().coastlines()
        
        # Add title to plots
        title = ' '.join([md['model'], md['long_name']])
        if layer:
            title = ' '.join([title,'(',layer,str(depth.units),')'])
        plt.title(title)
    
        # Determine filename:
        if multiModel:
             path = folder(cfg['plot_dir']) + os.path.basename(fn).replace(
                '.nc', '_map_'+str(l)+'.png')
        else:
            path = get_image_path(
                cfg,
                md,
                suffix='map_'+str(l),
                image_extention='png',
            )
            
        # Saving files:
        if cfg['write_plots']:

            logger.info('Saving plots to %s', path)
            plt.savefig(path)

        plt.close()
    

def multiModelTimeSeries(cfg, metadata, plotType='WeightedMean'):
    """
        This method makes a time series plot showing several models.
        """
    mdfiles = []
    for fn in sorted(metadata.keys()):

        print('-----------------')
        print(
            'model filenames:',
            fn,
        )
        # print('Varbiable name:', key)
        print('metadata[fn]:', metadata[fn])
        cube = iris.load_cube(fn)
        cube = sensibleUnits(cube, metadata[fn]['short_name'])

        multiModel = metadata[fn]['model'].find('MultiModel') > -1

        if plotType == 'WeightedMean':
            weights = iris.analysis.cartography.area_weights(
                cube, normalize=False)
            coords = ['longitude', 'latitude']
            cube = cube.collapsed(
                coords,
                iris.analysis.MEAN,
                weights=weights,
            )
            print('multiModelTimeSeries:', metadata[fn]['model'],
                  cube.coords('time'))

            if multiModel:
                qplt.plot(cube, label=metadata[fn]['model'], ls=':')
            else:
                qplt.plot(cube, label=metadata[fn]['model'])
            plt.legend(loc='best')

    if cfg['write_plots']:
        path = get_image_path(
            cfg,
            metadata[fn],
            prefix='MultiModel',
            suffix=plotType,
            image_extention='png',
            basenamelist=['field', 'short_name', 'start_year', 'end_year'],
        )
        logger.info('Saving plots to', path)
        plt.savefig(path)

    plt.close()


def main(cfg):
    #####

    input_files = get_input_files(cfg)

    print('cfg:\tContents:')
    for k in cfg.keys():
        print('CFG:\t', k, '\t', cfg[k])

    for i, metadatafilename in enumerate(cfg['input_files']):
        print(
            '\nmetadata filename:',
            metadatafilename,
        )
        
        metadata = get_input_files(cfg, index=i)

        #######
        # Multi model time series
        try:multiModelTimeSeries(cfg, metadata, plotType='WeightedMean')
        except:pass

        for fn in sorted(metadata.keys()):

            print('-----------------')
            print(
                'model filenames:\t',
                fn,
            )

            ######
            # Time series of individual model
            TimeSeriesPlots(cfg, metadata[fn], fn)

            ######
            # Map of temporal mean
            MapOfMeanOverTime(cfg, metadata[fn], fn)

    logger.debug("\n\nThis works\n\n")
    print('Success')


# cfg, settings = get_cfg()
# metadata = get_input_files(cfg)

if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
