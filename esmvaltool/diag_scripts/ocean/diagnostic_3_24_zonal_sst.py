"""
Diagnostics for fig 3.24 in Chapter 3 of IPCC AR6 WGI
========================


Author: Lee de Mora (PML)
        ledm@pml.ac.uk
Revised and corrected (15.01.2021): Elizaveta Malinina (CCCma)
                        elizaveta.malinina-rieger@canada.ca

"""
import logging
import os
import sys

import iris
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic
import esmvaltool.diag_scripts.shared.plot as eplot

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def titlify(title):
    """
    Check whether a title is too long then add it to current figure.

    Parameters
    ----------
    title: str
        The title for the figure.
    """
    cutoff = 40
    if len(title) > cutoff:
        # Find good mid point
        titles = title.split(' ')
        length = 0
        for itr, word in enumerate(titles):
            length += len(word)
            if length > cutoff:
                titles[itr] += '\n'
                length = 0.
        title = ' '.join(titles)
    plt.title(title)



def make_single_zonal_mean_plots(
        cfg,
        metadata,
        filename,
        obs_metadata={},
        obs_filename='',
):
    """
    Make a zonal mean error plot for an individual model.

    The optional observational dataset must be added.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    metadata: dict
        The metadata dictionairy for a specific model.
    filename: str
        The preprocessed model file.
    obs_metadata: dict
        The metadata dictionairy for the observational dataset.
    obs_filename: str
        The preprocessed observational dataset file.

    """
  # Load cube and set up units
    cube = iris.load_cube(filename)
    cube = diagtools.bgc_units(cube, metadata['short_name'])
    short_name = metadata['short_name']

    # Is this data is a multi-model dataset?
    multi_model = metadata['dataset'].find('MultiModel') > -1

   # Add observational data.
    if obs_filename:
        obs_cube = iris.load_cube(obs_filename)
        obs_cube = diagtools.bgc_units(obs_cube, metadata['short_name'])

        obs_key = obs_metadata['dataset']

    new_cube = cube - obs_cube

    # Zonal_mean_error
    if new_cube.data.shape == new_cube.coord('latitude').points.shape:
            plt.plot(new_cube.coord('latitude').points, new_cube.data )
            plt.ylabel('Latitude ('+r'$^\circ$'+'N)')
            key_word = 'Zonal mean SST error'

    # Equatorial_mean_error
    if new_cube.data.shape == new_cube.coord('longitude').points.shape:
            plt.plot(new_cube.coord('longitude').points, new_cube.data )
            plt.ylabel('Longitude ('+r'$^\circ$'+'E)')
            key_word = 'Equatorial SST error'

    plt.axhline(0., linestyle=':', linewidth=0.2, color='k')
    plt.ylabel('SST error ('+r'$^\circ$'+'C)')

    # Add title to plot
    if multi_model:
        title = ' '.join([key_word, 'multimodel mean', ])
    else:
        title = ' '.join([key_word, metadata['dataset'], ])
    plt.title(title)

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    if multi_model:
        path = diagtools.folder(
            cfg['plot_dir']) + os.path.basename(filename).replace(
                '.nc', key_word.replace(' ','') + image_extention)
    else:
        path = diagtools.get_image_path(
            cfg,
            metadata,
            suffix= key_word.replace(' ','')+short_name + image_extention,
        )

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)

    plt.close()


def match_model_to_key(
        model_type,
        cfg_dict,
        input_files_dict,
        variable_groups
):
    """
    Match up model or observations dataset dictionairies from config file.

    This function checks that the control_model, exper_model and
    observational_dataset dictionairies from the recipe are matched with the
    input file dictionairy in the cfg metadata.

    Arguments
    ---------
    model_type: str
        The string model_type to match (only used in debugging).
    cfg_dict: dict
        the config dictionairy item for this model type, parsed directly from
        the diagnostics/ scripts, part of the recipe.
    input_files_dict: dict
        The input file dictionairy, loaded directly from the get_input_files()
         function, in diagnostics_tools.py.

    Returns
    ---------
    dict
        A dictionairy of the input files and their linked details.
    """
    for input_file, intput_dict in input_files_dict.items():
        if intput_dict['variable_group'] not in variable_groups: continue
        intersect_keys = intput_dict.keys() & cfg_dict.keys()
        match = True
        for key in intersect_keys:
            if intput_dict[key] == cfg_dict[key]:
                continue
            match = False
        if match:
            return input_file
    logger.warning("Unable to match model: %s", model_type)
    return ''


def plot_zonal_cube(cube, plot_details):
    # Zonal_mean_error

    if cube.data.shape == cube.coord('latitude').points.shape:
        plt.plot(cube.coord('latitude').points, cube.data,
             c = plot_details['c'],
             ls = plot_details['ls'],
             )
        xlabel = 'Latitude ('+r'$^\circ$'+'N)'
        key_word = 'Zonal mean SST bias'

    # Equatorial_mean_error
    if cube.data.shape == cube.coord('longitude').points.shape:
        plt.plot(cube.coord('longitude').points, cube.data,
             c = plot_details['c'],
             ls = plot_details['ls'],
             )
        xlabel = 'Longitude ('+r'$^\circ$'+'E)'
        key_word = 'Equatorial SST bias'

    return key_word, xlabel


def fill_between_two_cubes(cube1, cube2, color):
    # Zonal_mean_error
    if cube1.data.shape == cube1.coord('latitude').points.shape:
        plt.fill_between(cube1.coord('latitude').points,
                        cube1.data,
                        cube2.data,
                        color=color,
                        alpha = 0.2,
                        linewidth = 0
             )

    if cube1.data.shape == cube1.coord('longitude').points.shape:
        plt.fill_between(cube1.coord('longitude').points,
                        cube1.data,
                        cube2.data,
                        color=color,
                        alpha = 0.2,
                        linewidth = 0
             )



def make_mean_of_cube_list(cube_list):
    """
    Takes the mean of a list of cubes (not an iris.cube.CubeList).

    Assumes all the cubes are the same shape.

    Careful, it will try to do these operators in place, need to make a copy.
    """
    cube_mean = cube_list[0].copy()
    for cube in cube_list[1:]:
        cube_mean+=cube
    cube_mean = cube_mean/ float(len(cube_list))
    cube_mean.units = 'celsius'
    return cube_mean


def make_std_of_cube_list(cube_list, axis):
    """
    Makes the standard deviation of a list of cubes (not an iris.cube.CubeList).

    Assumes all the cubes are the same shape and 1D.
    """
    cube_std = cube_list[0].copy()
    out_data = np.zeros_like(cube_std.data)

    if axis.lower() in ['lats', 'lat', 'latitude']:
        coords = cube_std.coord('latitude').points
    if axis.lower() in ['lons', 'lon', 'longitude']:
        coords = cube_std.coord('longitude').points

    out_dict = {coord:[] for coord in coords}

    for cube in cube_list:
        try:
            model = cube.attributes['model_id']
        except:
            print('No model id:', cube.attributes)
            print(cube.data)
            model = 'no model'

        print(model, cube.data[0])
        for l, coord in np.ndenumerate(coords):
            dat = cube.data.squeeze()[l]
            if np.isnan(dat) or np.ma.is_masked(dat):
                dat = 1.e20
            out_dict[coord].append(dat)

    for l, coord in np.ndenumerate(coords):
        dat = np.ma.array(out_dict[coord])
        dat = np.ma.masked_where((dat > 50.)+(dat < -20.), dat).compressed()
        print('itr:', l, 'coord:', coord, 'length:', len(dat), 'std:', np.std(dat), '\tdata:', dat)
        out_data[l] = np.std(dat)

    out_data = np.ma.masked_where(cube_std.data.mask, out_data)
    cube_std.data = out_data
    return cube_std

def make_perc_cube_list(cube_list, perc, axis):
    """
    Makes the percentiles of a list of cubes (not an iris.cube.CubeList).

    Assumes all the cubes are the same shape and 1D.
    """
    cube_perc = cube_list[0].copy()
    out_data = np.zeros_like(cube_perc.data)

    if axis.lower() in ['lats', 'lat', 'latitude']:
        coords = cube_perc.coord('latitude').points
    if axis.lower() in ['lons', 'lon', 'longitude']:
        coords = cube_perc.coord('longitude').points

    out_dict = {coord:[] for coord in coords}

    for cube in cube_list:
        try:
            model = cube.attributes['model_id']
        except:
            print('No model id:', cube.attributes)
            print(cube.data)
            model = 'no model'

        print(model, cube.data[0])
        for l, coord in np.ndenumerate(coords):
            dat = cube.data.squeeze()[l]
        #     print(model, l, coord, dat, cube.data.shape)
            if np.isnan(dat) or np.ma.is_masked(dat):
                dat = 1.e20
            out_dict[coord].append(dat)

    for l, coord in np.ndenumerate(coords):
        dat = np.ma.array(out_dict[coord])
        dat = np.ma.masked_where((dat > 50.)+(dat < -20.), dat).compressed()
        if dat.size == 0:
            print('itr:', l, 'coord:', coord, 'length:', len(dat),
                  str(perc) + ' perc: none', '\tdata:',
                  dat)
            if not(np.ma.is_masked(out_data[l])):
                out_data[l].mask = np.array(1, dtype=bool)
        else:
            print('itr:', l, 'coord:', coord, 'length:', len(dat),
                  str(perc)+' perc:', np.percentile(dat, perc), '\tdata:', dat)
            out_data[l] = np.percentile(dat, perc)

    out_data = np.ma.masked_where(cube_perc.data.mask, out_data)
    cube_perc.data = out_data
    return cube_perc

def load_obs(cfg, groups):
    """
    Load the observations.
    """
    obs_key = 'observational_dataset'
    obs_filename = ''
    obs_metadata = {}
    metadatas = diagtools.get_input_files(cfg)
    if obs_key in cfg:
        obs_filename = match_model_to_key(obs_key, #using local copy
                                                    cfg[obs_key],
                                                    metadatas,
                                                    groups)
        obs_metadata = metadatas[obs_filename]
        obs_cube = iris.load_cube(obs_filename)
        obs_cube = diagtools.bgc_units(obs_cube, obs_metadata['short_name'])
        obs_key = obs_metadata['dataset']
    return obs_cube, obs_key, obs_filename


def make_multimodle_zonal_mean_plots(
        cfg,
        pane = 'a',
        save = True,
        shortname = 'thetao',
):
    """
    Make a zonal mean error plot for an individual model.

    The optional observational dataset must be added.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.

    """
    metadatas = diagtools.get_input_files(cfg)
    plot_details = {}
    if pane in ['a']:
        groups = ['thetao_zonal', 'tos_zonal', ]
        obs_cube, obs_key, obs_filename = load_obs(cfg, groups)
        mask_arr = np.zeros(obs_cube.shape, dtype=bool)
    if pane in ['b', 'c']:
        groups = ['thetao_equator', 'tos_equator', ]
        obs_cube, obs_key, obs_filename = load_obs(cfg, groups)
        mask_cube, mask_key, mask_filename = load_obs(cfg, ['tos_eq_mask', ])
        mask_arr = [mask_cube.data[:,i].count()<=0.5*len(mask_cube.data[:,i])
                     for i in range(mask_cube.data.shape[1])]

    obs_cube.data.mask = obs_cube.data.mask | mask_arr

    #####
    # Load obs data and details
    # obs_cube, obs_key, obs_filename = load_obs(cfg, groups)

    plot_details = {}
    cmap = plt.cm.get_cmap('jet')

    #####
    # calculate the number of models
    number_models = {}
    projects = {}
    for i, filename in enumerate(sorted(metadatas)):
        metadata = metadatas[filename]
        if filename == obs_filename: continue
        if metadata['variable_group'] not in groups: continue
        number_models[metadata['dataset']] = True
        # highlight only the HiResMIP models
        if 'activity' in metadata:
            if metadata['activity'] == 'HighResMIP':
                projects['HighRes'] = True
            else:
                projects[metadata['project']] = True
        else:
            projects[metadata['project']] = True

    model_numbers = {model:i for i, model in enumerate(sorted(number_models))}
    print (number_models, model_numbers)
    number_models = len(number_models)

    #####
    # List of cubes to make means/stds.
    project_cubes = {project:{} for project in projects}

    for i, filename in enumerate(sorted(metadatas)):
        if filename == obs_filename: continue

        metadata = metadatas[filename]
        short_name = metadata['short_name']
        dataset = metadata['dataset']
        # correct selection of HiResMIP models
        if 'activity' in metadata:
            if metadata['activity'] == 'HighResMIP':
                project = 'HighRes'
            else:
                project = metadata['project']
        else:
            project = metadata['project']
        if metadata['variable_group'] not in groups:
            continue

        cube = iris.load_cube(filename)
        cube.data.mask = cube.data.mask | mask_arr

        coord_names = [coord.var_name for coord in cube.coords()]
        if len(cube.coord('longitude').points)>1:
            cube = cube.intersection(longitude=(20., 380.))
            # to make the pacific and atlantic continuous
            obs_cube = obs_cube.intersection(longitude=(20., 380.))
        cube = diagtools.bgc_units(cube, short_name)

        if number_models == 1:
            color = 'black'
        else:
            value = float(model_numbers[dataset]) / (number_models - 1.)
            color = cmap(value)

        # Is this data is a multi-model dataset?
        if dataset.find('MultiModel') > -1: continue

        new_cube = cube - obs_cube
        if new_cube.data.mean() < -200. :
                print ("Warning: this model is broken:", dataset, 'mean:', new_cube.data.mean())
                continue

        ####
        # Calculate the project lines
        project_cubes[project][dataset] = cube

    # Plot the project data range.
    for project in projects:
        for ds, cube in project_cubes[project].items():
            print(ds, '\t', project, '\t', cube.data.min(), '->', cube.data.min())

    # Plot the project means.
    for project in projects:
        if project in ['OBS', 'obs4mip']: continue

        ####
        # Calculate error
        errorcubeslist = [cube - obs_cube for cube in project_cubes[project].values()]
        project_mean_error = make_mean_of_cube_list(errorcubeslist)

        if project == 'CMIP5':
                mip_color  = (37/255, 81/255, 204/255)
        elif project == 'CMIP3':
                mip_color  = 'dodgerblue'
        elif project == 'CMIP6':
                mip_color  = (204/255, 35/255, 35/255)
        elif project == 'HighRes':
                mip_color = (30/256, 148/255, 130/255)
        else:  assert 0
        plot_details[project] = {'c': mip_color, 'ls': '-', 'label': project}
        if pane in 'ab':
                key_word, xlabel = plot_zonal_cube(project_mean_error, plot_details[project])
                ylabel = 'SST bias ('+r'$^\circ$'+'C)'

        if pane in 'a':
                cube_perc_5 = make_perc_cube_list(errorcubeslist, 5, 'lat')
                cube_perc_95 = make_perc_cube_list(errorcubeslist, 95, 'lat')
                plt.text(-75, -2.0 - list(projects.keys()).index(project)*0.45, project, c= mip_color)
                plt.ylim(-3.1, 2.75)
                fill_between_two_cubes(cube_perc_5, cube_perc_95, mip_color)

        if pane in 'ab':
                cube_perc_5 = make_perc_cube_list(errorcubeslist, 5, 'lon')
                cube_perc_95 = make_perc_cube_list(errorcubeslist, 95, 'lon')
                fill_between_two_cubes(cube_perc_5, cube_perc_95, mip_color)

        if pane in 'c':
                cubeslist = [cube  for cube in project_cubes[project].values()]
                project_mean = make_mean_of_cube_list(cubeslist)
                key_word, xlabel = plot_zonal_cube(project_mean, plot_details[project])

                cube_perc_5 = make_perc_cube_list(cubeslist, 5, 'lon')
                cube_perc_95 = make_perc_cube_list(cubeslist, 95, 'lon')
                print('project_mean', project_mean.data.min(), project_mean.data.max())
                fill_between_two_cubes(cube_perc_5, cube_perc_95, mip_color)

                plot_details[obs_key] = {'c': 'black', 'ls': '-', 'label': obs_key}
                key_word, xlabel = plot_zonal_cube(obs_cube, plot_details[obs_key])
                key_word = 'Equatorial SST'
                ylabel = 'SST ('+r'$^\circ$'+'C)'

    if pane in 'bc':
        plt.xticks(np.arange(30,380,30))

    #####
    # title and axis lables
    if pane in 'ab':
        plt.axhline(0., linestyle='dashed', linewidth=1, color='grey',
                    zorder=1, alpha=0.7)
    else:
        plt.text(240, 29, 'OBS: ' + obs_key+' v1')
        plt.ylim(23.9, 30.2)
        plt.vlines(130, 23, 27.5, colors='silver', linestyles='dashed',
                   linewidth=1)
        plt.text(60, 24.2, 'Indian\nocean', color='grey')
        plt.vlines(285, 23, 25, colors='silver', linestyles='dashed',
                   linewidth=1)
        plt.text(180, 24.2, 'Pacific\nocean', color='grey')
        plt.text(320, 24.2, 'Atlantic\nocean', color='grey')
    if pane in 'b':
        plt.ylim(-2.5, 3.1)
        plt.text(60, -2.3, 'Indian\nocean', color='grey')
        plt.vlines(130, -3.3, -0.8, colors='silver', linestyles='dashed',
                   linewidth=1)
        plt.text(180, -2.3, 'Pacific\nocean', color='grey')
        plt.vlines(285, -3.3, -0.7, colors='silver', linestyles='dashed',
                   linewidth=1)
        plt.text(320, -2.15, 'Atlantic\nocean', color='grey')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    title = ' '.join(['(', pane, ')', key_word ])
    plt.title(title)

    if not save:
        return plt.gca(), plot_details

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    path = diagtools.get_image_path(
        cfg,
        metadata,
        suffix= key_word.replace(' ','') + pane + short_name+ image_extention,
        )


def main(cfg):
    """
    Load the config file and some metadata, then pass them the plot making
    tools.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.

    """
    #####
    # individual panes
    for pane in ['a', 'b', 'c', ]:
        make_multimodle_zonal_mean_plots(cfg, pane=pane, save = True)

    st_file = eplot.get_path_to_mpl_style(cfg.get('mpl_style'))

    plt.style.use(st_file)

    #####
    # Altogether
    fig = plt.figure()
    fig.set_size_inches(6, 8)
    plot_details = {}
    axes = []
    for pane,sbpt in zip([ 'a', 'b', 'c', ], [311,312,313]):
        axes.append(plt.subplot(sbpt))
        ax, pt_dets = make_multimodle_zonal_mean_plots(cfg, pane=pane, save=False)
        plot_details.update(pt_dets)

    fig.suptitle('Sea Surface Temperature (SST)')

    fig.subplots_adjust(top=0.93, left=0.15, right=0.9, bottom=0.06, hspace=0.4)

    # Load image format extention and path
    image_extention = diagtools.get_image_format(cfg)
    path = cfg['plot_dir'] + '/fig_3.24'+image_extention

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        fig.savefig(path)
        fig.savefig(path[:-4] + '.png', dpi=250)

    plt.close()

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)