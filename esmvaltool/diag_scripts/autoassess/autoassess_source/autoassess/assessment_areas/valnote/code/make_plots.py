'''
This module contains code to make all the validation note plots.
'''

# Import modules
import csv
import os
import pdb
import string
import sys

import matplotlib
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import iris
import iris.analysis
import iris.analysis.maths
import iris.plot as iplt
import iris.quickplot as qplt

import ImageMetaTag as imt

import globalvar                   # This is the global variable holder
import rms
import valmod as vm


def load_one_dict(file_dict, temp_dir):
    """
    Loads all the data in either source_dict or obs_dict
    file_dict = either source_dict or obs_dict
    """

    # Set the global variable that will hold all the data
    global all_data_dict

    # Loop over items in this dictionary
    for model_key in file_dict.iterkeys():
        for info_key in file_dict[model_key].iterkeys():
            if vm.key_valid(info_key):

                # What is the merged key and pickle file name
                merged_key = model_key+'_'+info_key
                pickle_file = os.path.join(temp_dir, merged_key+'.pkl')
                pickle_contents_file = \
                    os.path.join(temp_dir, merged_key+'.contents')
                input_file = file_dict[model_key][info_key][0].strip()

                # Should we load the pickle file
                load_pickle = False
                if globalvar.pkl:
                    if os.path.isfile(pickle_file) and \
                       os.path.isfile(pickle_contents_file):
                        pickle_contents = \
                            vm.read_info_file(pickle_contents_file)['filename']
                        if pickle_contents == input_file:
                            load_pickle = True

                # Load the pickle file if you need to
                if load_pickle:
                    print 'Loading data '+merged_key+' from file '+pickle_file
                    cube_list = vm.load_pickle(pickle_file)

                # Otherwise load the data from the source file
                else:
                    print 'Loading data '+merged_key+' from file '+input_file
                    cube_list = iris.load(input_file)

                    # Save to pickle file
                    if globalvar.pkl:
                        vm.save_pickle(cube_list, pickle_file)
                        vm.write_info_file(pickle_contents_file,
                                           {'filename': input_file})

                # Save the data in all_data_dict
                all_data_dict[merged_key] = cube_list


def load_all_data(do_main, temp_dir):
    """Loads in all the data into one big dictionary"""

    # Initialise variable
    global all_data_dict
    all_data_dict = {}

    # Get the data from the source file
    load_one_dict(source_dict, temp_dir)

    # Get the data from the obs file
    if do_main:
        load_one_dict(obs_dict, temp_dir)


def add_constraint(all_constraints, axis_str, input_info):
    """
    Adds a constraint to all_constraints based on axis_str which is one of:
    pressure, model_level_number, soil_model_level_number, latitude, longitude,
    grid_latitude, grid_longitude
    """

    axis_value_str = globalvar.item_dict[input_info[2]][axis_str][0]

    # You can set a range of a coordinate to load between. E.g. a range of
    #  latitudes be used for meridional means.
    if len(globalvar.item_dict[input_info[2]][axis_str]) > 1:
        axis_max_str = globalvar.item_dict[input_info[2]][axis_str][1]
        lamfn = lambda cell: int(axis_value_str) <= cell <= int(axis_max_str)
        this_constraint = iris.Constraint(coord_values={axis_str: lamfn})

    # Or just select a single point for this axis.
    #  E.g. just the first model level
    else:
        ival = int(axis_value_str)
        this_constraint = iris.Constraint(coord_values={axis_str: ival})

    # Combine the constraints. Seeing as the stash constraints is a list
    #  we will need to keep the combination in list format.
    n_constraints = len(all_constraints)
    for c_num in range(n_constraints):
        all_constraints[c_num] = all_constraints[c_num] & this_constraint


def get_data(components):
    """Gets all the data needed for all the components of an equation

    components = a list of components where each item is a string of the
                 form source_season_field
    """

    # Initialise data dictionary
    global all_data_dict
    data_dict = {}

    # Loop over components
    for comp in components:

        # Split the component into info about what to load
        input_info = string.split(comp, '_')

        # Make sure input_info exists in either the obs_dict or source_dict
        if not (input_info[0] in obs_dict.keys() or
                input_info[0] in source_dict.keys()):
            print "ERROR: input info[0] not found in either"
            print "obs_file.dat or source_file.dat"
            print "input_info[0]="+input_info[0]
            print "obs_dict.keys()=", obs_dict.keys()
            print "source_dict.keys()=", source_dict.keys()
            if globalvar.debug:
                pdb.set_trace()
            else:
                globalvar.error = True
                return 1

        # Set up extra variables to control the constraints used to extract
        # fields, constraint_limiter is a list of strings and can contain
        # one or more of:
        #  'use_all'      = use all the constraints as listed in item_file.dat
        #                   for your field (default)
        #  'ignore_stash' = do not constrain by stash code or standard name.
        #  'single_level' = do not constrain by model_level_number
        constraint_type = 'stash'
        constraint_limiter = ['use_all']

        if input_info[0] not in source_dict.keys():

            # Check to see if this is a file that contains a single field
            if 'constraint_limiter' in obs_dict[input_info[0]].keys():
                # Limit the constraints in some way.
                constraint_limiter = \
                    obs_dict[input_info[0]]['constraint_limiter']

            # Check to see if the input data is in netcdf format
            filenames = []
            for season in ['djfm', 'mamm', 'jjam', 'sonm', 'annm']:
                if season in obs_dict[input_info[0]].keys():
                    filenames.append(obs_dict[input_info[0]][season])

            # What is the constraint_type
            if filenames[0][0][-3:] == '.nc':
                constraint_type = 'standard_name'

        # Get the information about the constraint_str from the item dictionary
        # constraint_str can either contain the stash code or the standard_name
        # (depending on whether it is pp or netcdf format)
        if 'ignore_stash' not in constraint_limiter:

            # Initialise a list of field constraints
            all_constraints = []

            try:
                n_constraints = \
                    len(globalvar.item_dict[input_info[2]][constraint_type])
            except:
                print "ERROR encountered extracting key " + input_info[2] + \
                      " from globalvar.item_dict."
                print "ERROR: ", sys.exc_info()[0]
                if globalvar.debug:
                    pdb.set_trace()
                else:
                    globalvar.error = True
                    return 1

            # Loop over constraints
            for c_num in range(n_constraints):

                # Get the constraint
                constraint_str = \
                    globalvar.item_dict[input_info[2]][constraint_type][c_num]

                if constraint_type == 'standard_name':
                    # Setup the name constraints. The second constraint is not
                    #  needed in this case but needs to be there to be
                    #  compatible with the stash version of all_constraints
                    all_constraints.append(iris.Constraint(name=constraint_str))

                else:
                    # Setup the stash loading constraints
                    (section, item) = divmod(int(constraint_str), 1000)
                    all_constraints.extend([
                        iris.AttributeConstraint(STASH=(01, section, item)),
                        iris.AttributeConstraint(STASH=(None, section, item))])

        # Set up blank constraints for the single field case.
        else:
            all_constraints = [iris.Constraint()]

        # Setup the extra loading constraints based on pressure,
        #  latitude and longitude
        if 'pressure' in globalvar.item_dict[input_info[2]]:
            add_constraint(all_constraints, 'pressure', input_info)
        if 'model_level_number' in globalvar.item_dict[input_info[2]] and \
           'single_level' not in constraint_limiter:
            add_constraint(all_constraints, 'model_level_number', input_info)
        if 'soil_model_level_number' in globalvar.item_dict[input_info[2]] and \
           'single_level' not in constraint_limiter:
            add_constraint(all_constraints, 'soil_model_level_number', input_info)
        if 'use_grid' in constraint_limiter:
            if 'grid_longitude' in globalvar.item_dict[input_info[2]]:
                add_constraint(all_constraints, 'grid_longitude', input_info)
            if 'grid_latitude' in globalvar.item_dict[input_info[2]]:
                add_constraint(all_constraints, 'grid_latitude', input_info)
        else:
            if 'longitude' in globalvar.item_dict[input_info[2]]:
                add_constraint(all_constraints, 'longitude', input_info)
            if 'latitude' in globalvar.item_dict[input_info[2]]:
                add_constraint(all_constraints, 'latitude', input_info)

        # Print to screen what we are loading
        merged_key = input_info[0]+'_'+input_info[1]
        print 'Extracting field from '+merged_key

        # Check that the key exists
        if merged_key not in all_data_dict:
            print ('ERROR: Cannot find key {} in all_data_dict. Look in the ' +
                   'file source_file.dat and check that you have means for ' +
                   'all the seasons for both experiment and ' +
                   'control.').format(merged_key)
            if globalvar.debug:
                pdb.set_trace()
            else:
                raise Exception

        # Load the data
        data_dict[comp] = vm.extract_no_diurnal(all_data_dict[merged_key],
                                                all_constraints=all_constraints)

        # Skip if an error has occured
        if globalvar.error:
            return 1

        # What are the stash code you finally retrieved
        if 'STASH' in data_dict[comp].attributes:
            section = data_dict[comp].attributes['STASH'].section
            item = data_dict[comp].attributes['STASH'].item

        # Apply the heavyside function
        if (input_info[0] in source_dict) and \
           section == 30 and (200 < item < 400):

            # Load the heavyside function using the same constraints as
            #  the model data
            all_constraints = [iris.AttributeConstraint(STASH=(01, 30, 301))]
            if 'pressure' in globalvar.item_dict[input_info[2]]:
                add_constraint(all_constraints, 'pressure', input_info)
            if 'longitude' in globalvar.item_dict[input_info[2]]:
                add_constraint(all_constraints, 'longitude', input_info)
            if 'latitude' in globalvar.item_dict[input_info[2]]:
                add_constraint(all_constraints, 'latitude', input_info)
            heavyside_cube = vm.extract_no_diurnal(all_data_dict[merged_key],
                                                   all_constraints=all_constraints)

            # Skip if an error has occured
            if globalvar.error:
                return 1

            # Sometimes the shape of the loaded cube is not the same as
            #  the shape of the heavyside. When this happens try to match the
            #  shapes.
            if data_dict[comp].shape != heavyside_cube.shape:
                data_dict_reshape = data_dict[comp]
                if data_dict_reshape.shape == heavyside_cube.shape:
                    data_dict[comp] = data_dict_reshape
                else:
                    print "ERROR: Can't match cube shape with heavyside " + \
                        "shape for stash "+constraint_str+" in "+input_info[0]
                    if globalvar.debug:
                        pdb.set_trace()
                    else:
                        globalvar.error = True
                        return 1

            print 'Applying heavyside'
            data_dict[comp] = vm.apply_heavyside(data_dict[comp],
                                                 heavyside_cube)

        # Scale up the model data if requested
        if (input_info[0] in source_dict) and \
           'scale_model' in globalvar.item_dict[input_info[2]]:
            scale_model = \
                eval(globalvar.item_dict[input_info[2]]['scale_model'][0])
            data_dict[comp] = data_dict[comp] * scale_model

    return data_dict


def define_page(page_title='', num_plots=4):
    """
    Set up the page for plotting

    page_title = (string) the title to go at the top of the page
    num_plots = (int) the number of plots on the page
    """

    global ax1, location_dict

    # Set up defaults
    figsize = (12, 9)

    # Set up margins
    if globalvar.plot_type != 'lat_lon':
        left = 0.08
        right = 0.95
        bottom = 0.05
        top = 0.94
        wspace = 0.25
        hspace = 0.25
    else:
        left = 0.05
        right = 0.95
        bottom = 0.05
        top = 0.95
        wspace = 0.17
        hspace = 0.1

    # Make an array of locations
    if num_plots == 1:
        location_dict = {'a': 111}
    elif num_plots == 2:
        location_dict = {'a': 121, 'b': 122}
    elif num_plots == 3:
        location_dict = {'a': 221, 'b': 222, 'c': 223}
    elif num_plots == 4:
        location_dict = {'a': 221, 'b': 222, 'c': 223, 'd': 224}
    elif num_plots >= 5 and num_plots <= 9:
        location_dict = {'a': 331, 'b': 332, 'c': 333, 'd': 334,
                         'e': 335, 'f': 336, 'g': 337, 'h': 338, 'i': 339}

    # Set up the font size
    if num_plots <= 4:
        font_size = 12
    else:
        font_size = 8

    if 'layout' in title_dict[title_type]:
        layout = title_dict[title_type]['layout'][0]
        if layout == '1by3':
            location_dict = {'a': 311, 'b': 312, 'c': 313}
            figsize = (5, 9)
            if globalvar.plot_type != 'lat_lon':
                left += 0.1
            font_size = 12
        elif layout == '2by3':
            location_dict = {'a': 321, 'b': 322, 'c': 323,
                             'd': 324, 'e': 325, 'f': 326}
            figsize = (10, 9)
            font_size = 12
        elif layout == '2by4':
            location_dict = {'a': 421, 'b': 422, 'c': 423, 'd': 424,
                             'e': 425, 'f': 426, 'g': 427, 'h': 428}
            figsize = (10, 12)
            font_size = 12
        else:
            print 'WARNINBG: layout ' + layout + \
                ' in title_file.dat not supported. Ignoring and using defaults.'

    # Define a linestyle for negative values (this is a global setting and
    #  will require a change to matplotlib to make it local)
    matplotlib.rcParams['contour.negative_linestyle'] = 'dashed'
    format_string = '%d'

    # Set up the page
    if 'ax1' in globals():
        ax1.set_size_inches(figsize)
    else:
        ax1 = plt.figure(1, figsize=figsize)
    plt.subplots_adjust(left=left, right=right, bottom=bottom,
                        top=top, wspace=wspace, hspace=hspace)

    matplotlib.rcParams.update({'font.size': font_size})


def contour_plot(toplot_cube_list, location, page_title, key, equation_list,
                 filename, rms_float=None):
    """
    Make a contour plot

    toplot_cube_list = (cube) a 2D cube (or list of cubes) to plot
    location = (int) the location of the plot on the page
    page_title = the title of the page from valorder.dat
    key = the key (a, b, c or d) for the plot
    equation = the equation that made the plot
    filename = the filename of the file that this plot is going into.
               This is needed for the rms table.
    rms_float = root mean square value to add to the bottom of the plot
    """

    global ax1

    # Set up the levels
    level_type = globalvar.plots_dict[page_title][key][-1]
    if 'levels' in globalvar.levels_dict[level_type]:
        levels_list = eval(globalvar.levels_dict[level_type]['levels'][0])

    else:
        # Levels not supplied. Set levels to go between min and max input data.
        max_level = numpy.amax(toplot_cube_list[0].data)
        min_level = numpy.amin(toplot_cube_list[0].data)
        print 'Dynamically assigning levels {} {}'.format(min_level, max_level)
        inc = (max_level-min_level)/15
        levels_list = numpy.arange(min_level, max_level+inc, inc).tolist()

    # Do we extend the colour bar. Valid values are:
    #    both = extend in both directions (defualt)
    #    neither = do not extend colour bar
    #    min = extend colour bar to minimum value
    #    max = extend colour bar to maximum value
    if 'extend' in globalvar.levels_dict[level_type]:
        extend = globalvar.levels_dict[level_type]['extend'][0]
    else:
        extend = 'both'

    # Set up the colours by first generating a colour map ...
    red_list = eval(globalvar.levels_dict[level_type]['red'][0])
    green_list = eval(globalvar.levels_dict[level_type]['green'][0])
    blue_list = eval(globalvar.levels_dict[level_type]['blue'][0])
    cmap = vm.make_cmap(red_list, green_list, blue_list)

    # ... and then mapping those colours into a colours list
    #     (using a list comprehension)
    colors_list = [cmap(float(contour_num)/len(levels_list))
                   for contour_num in range(len(levels_list)+1)]

    # Get the contour type
    contour_type = 'default'
    if 'contour_type' in globalvar.levels_dict[level_type]:
        contour_type = globalvar.levels_dict[level_type]['contour_type'][0]

    # Load the basic title
    title_append = ' '
    if (globalvar.plot_type in ('y_pressure', 'y_height')) and not globalvar.pub:
        title_append = ' Zonal mean '
    title_top = title_append+title_dict[title_type][key][0]
    title_bottom = title_dict[title_type][key][1]
    observations = 'obs'
    title_model = 'model'

    # Alter the title to add information about this plot
    comp = vm.split_equation(equation_list[0])[0]
    input_info = string.split(comp, '_')
    model_title_dict = {}
    if globalvar.pub:
        for model in source_dict:
            model_title_dict[model] = source_dict[model]['name'][0]
    else:
        for model in source_dict:
            model_title_dict[model] = source_dict[model]['jobid'][0].upper() \
                + ': '+source_dict[model]['name'][0]

    # Get the field title from either the plots file or the item file
    if 'title' in globalvar.plots_dict[page_title]:
        field_title = globalvar.plots_dict[page_title]['title'][0]
    else:
        field_title = globalvar.item_dict[input_info[2]]['title'][0]

    # Scale the plot if required
    n_plots = len(toplot_cube_list)
    if 'scale_plot' in globalvar.levels_dict[level_type]:
        scale_plot = eval(globalvar.levels_dict[level_type]['scale_plot'][0])
        for plot_num in range(n_plots):
            toplot_cube_list[plot_num] *= scale_plot

    # Shift the grid so that it runs from 180W to 180E
    if globalvar.plot_type == 'lat_lon':
        for plot_num in range(n_plots):
            toplot_cube_list[plot_num] = \
                vm.roll_cube(toplot_cube_list[plot_num])

    # If height is a vertical axis convert to km
    if globalvar.plot_type == 'y_height':
        for plot_num in range(n_plots):
            aux_coord_names = [aux_coord.long_name for aux_coord in
                               toplot_cube_list[plot_num].aux_coords]
            height_index = aux_coord_names.index('level_height')
            toplot_cube_list[plot_num].aux_coords[height_index].convert_units('km')

    # Set up the projection
    projection=None
    if 'projection' in globalvar.plots_dict[page_title]:
        projection = eval('ccrs.'+globalvar.plots_dict[page_title]['projection'][0])
        globalvar.plot_type = 'lat_lon_regional'
    if globalvar.plot_type == 'lat_lon':
        projection = ccrs.PlateCarree()

    # Set up the location of the sub plot
    thisplot = ax1.add_subplot(location, projection=projection)

    # For some reason the use of bespoke contour levels requires the following
    #  line with Orthographic projections. Go figure!!!!
    if isinstance(projection, ccrs.Orthographic):
        thisplot.set_global()
    else:
        if 'region' in globalvar.plots_dict[page_title]:
            region = eval(globalvar.plots_dict[page_title]['region'][0])
            thisplot.set_extent(region, crs=ccrs.PlateCarree())

    # Make the plot
    try:
        for plot_num in range(n_plots):

            toplot_cube = toplot_cube_list[plot_num]
            equation = equation_list[plot_num]

            # Determine what the label is
            component_set = vm.get_component_set(equation)
            for comp in component_set:
                input_info = string.split(comp, '_')
                if input_info[0] not in source_dict.keys():
                    observations = obs_dict[input_info[0]]['name'][0]
                else:
                    title_model = source_dict[input_info[0]]['jobid'][0]
                    if not globalvar.pub:
                        title_model = title_model.upper() + ': ' + \
                            source_dict[input_info[0]]['name'][0]

            # The first plot is a filled contour plot
            if plot_num == 0:
                if contour_type == 'block':
                    plot_object = iplt.pcolormesh(toplot_cube, cmap=cmap,
                                                  vmin=levels_list[0],
                                                  vmax=levels_list[-1])
                else:
                    plot_object = iplt.contourf(toplot_cube,
                                                levels=levels_list,
                                                colors=colors_list,
                                                extend=extend)

            # The second plot is either a contour plot or saved for
            #  a future vector plot
            if plot_num == 1:
                if n_plots == 2:
                    plot_object = iplt.contour(toplot_cube)
                else:
                    u_cube = toplot_cube.copy()
                    print 'u_cube =', u_cube.data[100, 100]

            # The third plot is a vector plot
            if plot_num == 2:
                v_cube = toplot_cube
                print 'v_cube =', v_cube.data[100, 100]
                vector_scale = 30
                if 'vector_scale' in globalvar.levels_dict[level_type]:
                    vector_scale = eval(globalvar.levels_dict[level_type]['vector_scale'][0])
                vm.quiver(u_cube, v_cube, scale=vector_scale)

            # Any more plots cause an error
            if plot_num >= 3:
                print "ERROR. Only 3 overlying plots allowed with 2d fields"
                raise exception

    except:
        print "ERROR encountered plotting "+page_title+", plot number "+key
        print "ERROR: ", sys.exc_info()
        if globalvar.debug:
            pdb.set_trace()
        else:
            globalvar.error = True

            # Add a blank plot with an error message
            plt.cla()
            plt.plot([0, 0, 0, 0])
            plt.title("Error making plot")

            return 1

    # Replace parts of the title top with new values
    title_top = title_top.replace('field', field_title)
    title_top = title_top.replace('season', input_info[1][0:3])
    for model_name in model_title_dict.keys():
        title_top = title_top.replace(model_name, model_title_dict[model_name])
    title_top = title_top.replace('observations', observations)
    title_top = title_top.replace('model', title_model)
    title_bottom = title_bottom.replace('field', field_title)
    title_bottom = title_bottom.replace('season', input_info[1][0:3])
    for model_name in model_title_dict.keys():
        title_bottom = \
            title_bottom.replace(model_name, model_title_dict[model_name])
    title_bottom = title_bottom.replace('observations', observations)
    title_bottom = title_bottom.replace('model', title_model)

    # Make the plot go from high to low. This should ideally come
    #  from the Cube's metadata. i.e. pressure coord direction
    if globalvar.plot_type in ('y_pressure', 'y_height'):
        if globalvar.plot_type == 'y_pressure':
            plt.gca().invert_yaxis()
        plt.gca().invert_xaxis()
        if 'logy' in globalvar.plots_dict[page_title]:
            plt.gca().semilogy()
            plt.ylim([1000, 10])

    # Plot contour lines (as long as no_contours is not set in contour
    #  type in the levels_file.dat file
    if contour_type != 'no_contours' and contour_type != 'block':
        cs = iplt.contour(toplot_cube, levels=plot_object.levels,
                          colors='k', extend=extend)
        format_string = '%d'
        # Note that although the MatPlotLib documentation states you can use
        #  relative fontsizes (small, large, etc.) this does not appear to
        #  work hence the use of an explicit fontsize.
        try:
            plt.clabel(cs, colors='black', inline_spacing=0, fmt=format_string,
                       fontsize=8)
        except ValueError:
            print 'Error labelling contours for plot ' + key + ' of page ' + \
                page_title
            print "ERROR: ", sys.exc_info()
            print 'Not labelling contours and continuing onwards.'

    # Set the plot title and type
    if globalvar.pub:
        plt.title(title_top+'\n'+key+') '+title_bottom, fontsize=10)
    else:
        plt.title(key+') '+title_top+'\n'+title_bottom, size='medium')
    if globalvar.plot_type == 'y_pressure':
        plt.ylabel('Pressure / hPa', size='small')
    if globalvar.plot_type == 'y_height':
        plt.ylabel('Height / km', size='small')

    # Add the RMS values to the bottom of the plot
    if rms_float is not None:
        if rms_float != globalvar.missing_data:
            rms_format = '{0:.2e}'
            if rms_float >= 0.1:
                rms_format = '{0:.3f}'
            if rms_float >= 1.0:
                rms_format = '{0:.2f}'
            if rms_float >= 10.0:
                rms_format = '{0:.1f}'
            if rms_float >= 100.0:
                rms_format = '{0:.0f}'
            if rms_float >= 1000.0:
                rms_format = '{0:.2e}'
            rms_string = 'Area-weighted rms diff = '+rms_format
            plt.text(0.5, -0.17, rms_string.format(rms_float),
                     transform=plt.gca().transAxes,
                     horizontalalignment='center', size='small')

    # Add the coastlines
    if globalvar.plot_type in ['lat_lon', 'lat_lon_regional']:
        plt.gca().coastlines()

    # Add the colour bar
    units = ''
    if 'units' in globalvar.levels_dict[level_type]:
        units = eval(globalvar.levels_dict[level_type]['units'][0])
    orientation='horizontal'
    if 'colour_bar' in globalvar.plots_dict[page_title]:
        orientation = globalvar.plots_dict[page_title]['colour_bar'][0].strip()
    vm.add_colorbar(plot_object, units, orientation=orientation)

    # Correct the longitude and latitude ticks
    if globalvar.plot_type in ('y_pressure', 'y_height'):
        thisplot.xaxis.set_major_formatter(LATITUDE_FORMATTER)
        thisplot.xaxis.set_ticks(np.linspace(-80, 80, num=9))
        thisplot.xaxis.set_tick_params(labelsize='small')
        thisplot.yaxis.set_tick_params(labelsize='small')
    elif globalvar.plot_type in ('lat_lon', ):
        glines = thisplot.gridlines(draw_labels=True)
        glines.xlines = None
        glines.ylines = None
        glines.xlabels_top = False
        glines.ylabels_right = False
        #glines.xlocator = mticker.FixedLocator(np.linspace(-180, 180, num=5))
        #glines.ylocator = mticker.FixedLocator(np.linspace(-80, 80, num=5))
        glines.xformatter = LONGITUDE_FORMATTER
        glines.yformatter = LATITUDE_FORMATTER
        glines.xlabel_style = {'size': 'small'}
        glines.ylabel_style = {'size': 'small'}

    # Alter the size of the plots to better fit everything on
    if globalvar.plot_type in ('y_pressure', 'y_height'):
        fig = plt.figure(1, figsize=(7, 4))


def line_plot(toplot_cube_list, location, page_title, key, equation_list,
              filename, rms_float=None):
    """
    Make a line plot

    toplot_cube_list = (cube) a 1D cube (or list of cubes) to plot
    location = (int) the location of the plot on the page
    page_title = the title of the page from valorder.dat
    key = the key (a, b, c or d) for the plot
    equation = the equation that made the plot
    filename = the filename of the file that this plot is going into.
               This is needed for the rms table.
    rms_float = root mean square value to add to the bottom of the plot
    """
    global ax1

    # Set up the levels
    level_type = globalvar.plots_dict[page_title][key][-1]
    if 'levels' in globalvar.levels_dict[level_type]:
        levels_list = eval(globalvar.levels_dict[level_type]['levels'][0])

    # Load the basic title
    title_top = title_dict[title_type][key][0]
    title_bottom = title_dict[title_type][key][1]
    observations = 'obs'

    # Alter the title to add information about this plot
    comp = vm.split_equation(equation_list[0])[0]
    input_info = string.split(comp, '_')
    model_title_dict = {}
    if globalvar.pub:
        for model in source_dict:
            model_title_dict[model] = source_dict[model]['name'][0]
    else:
        for model in source_dict:
            model_title_dict[model] = source_dict[model]['jobid'][0].upper() \
                + ': '+source_dict[model]['name'][0]

    # Get the field title from either the plots file or the item file
    if 'title' in globalvar.plots_dict[page_title]:
        field_title = globalvar.plots_dict[page_title]['title'][0]
    else:
        field_title = globalvar.item_dict[input_info[2]]['title'][0]

    # Scale the plot if required
    n_plots = len(toplot_cube_list)
    if 'scale_plot' in globalvar.levels_dict[level_type]:
        scale_plot = eval(globalvar.levels_dict[level_type]['scale_plot'][0])
        for plot_num in range(n_plots):
            toplot_cube_list[plot_num] *= scale_plot

    # Make a list to hold line style information.
    # See fmt in http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot
    line_style_list = ['-', '--', ':', '-.', 'o', '^']

    # Set up the location of the sub plot
    thisplot = ax1.add_subplot(location)

    # Make the plots
    try:
        for plot_num in range(n_plots):

            toplot_cube = toplot_cube_list[plot_num]
            equation = equation_list[plot_num]

            # Determine what the label is
            component_set = vm.get_component_set(equation)
            input_info = string.split(component_set.pop(), '_')
            if input_info[0] not in source_dict.keys():
                observations = obs_dict[input_info[0]]['name'][0]
                label = observations
            else:
                label = model_title_dict[input_info[0]]

            plot_object = iplt.plot(toplot_cube, line_style_list[plot_num],
                                    label=label)
    except:
        print "ERROR encountered plotting "+page_title+", plot number "+key
        print "ERROR: ", sys.exc_info()
        if globalvar.debug:
            pdb.set_trace()
        else:
            globalvar.error = True

            # Add a blank plot with an error message
            plt.cla()
            plt.plot([0, 0, 0, 0])
            plt.title("Error making plot")
            return 1

    # Replace parts of the title top with new values
    title_top = title_top.replace('field', field_title)
    title_top = title_top.replace('season', input_info[1][0:3])
    for model_name in model_title_dict:
        title_top = title_top.replace(model_name, model_title_dict[model_name])
    title_top = title_top.replace('observations', observations)
    title_bottom = title_bottom.replace('field', field_title)
    title_bottom = title_bottom.replace('season', input_info[1][0:3])
    for model_name in model_title_dict:
        title_bottom = title_bottom.replace(model_name,
                                            model_title_dict[model_name])
    title_bottom = title_bottom.replace('observations', observations)

    # Add the legend
    plt.legend()

    # Limit the y axis
    if 'levels' in globalvar.levels_dict[level_type]:
        plt.ylim(levels_list)

    # Set the plot title and type
    if globalvar.pub:
        plt.title(title_top+'\n'+key+') '+title_bottom, fontsize=10)
    else:
        plt.title(key+') '+title_top+'\n'+title_bottom)
    if globalvar.plot_type == 'x_longitude':
        plt.xlabel('Longitude')
    if globalvar.plot_type == 'x_latitude':
        plt.xlabel('Latitude')
    units = ''
    if 'units' in globalvar.levels_dict[level_type]:
        units = eval(globalvar.levels_dict[level_type]['units'][0])
    plt.ylabel(units)

    # Add the RMS values to the bottom of the plot
    if rms_float is not None:
        if rms_float != globalvar.missing_data:
            rms_format = '{0:.2e}'
            if rms_float >= 0.1:
                rms_format = '{0:.3f}'
            if rms_float >= 1.0:
                rms_format = '{0:.2f}'
            if rms_float >= 10.0:
                rms_format = '{0:.1f}'
            if rms_float >= 100.0:
                rms_format = '{0:.0f}'
            if rms_float >= 1000.0:
                rms_format = '{0:.2e}'
            rms_string = 'Area-weighted rms diff = '+rms_format
            plt.text(0.5, -0.17, rms_string.format(rms_float),
                     transform=plt.gca().transAxes,
                     horizontalalignment='center', size='small')

    # Correct the longitude and latitude ticks
    if globalvar.plot_type == 'x_latitude':
        thisplot.xaxis.set_ticks(range(-90, 90, 30))
        thisplot.xaxis.set_major_formatter(LATITUDE_FORMATTER)
    if globalvar.plot_type == 'x_longitude':
        # thisplot.xaxis.set_ticks(range(-180, 181, 90))
        thisplot.xaxis.set_major_formatter(LONGITUDE_FORMATTER)

    # Alter the size of the plots
    plt.subplots_adjust(left=0.09, right=0.91, top=0.91,
                        bottom=0.09, hspace=0.3)


def create_tags_from_filename(filename):
    '''
    Stripping code from vnc.py to get image tags
    '''
    # Split the string into two halves. The first contains the number,
    # variable name and season. The second contains the observations.
    if filename.count('_v_') == 1:
        sections_list = filename.split('_v_')
        first_part = sections_list[0]
        second_part = sections_list[1]
    else:
        first_part = filename.split('.png')[0]
        second_part = "none.png"

    # Get the variable name after the first underscore and before the last
    # underscore in this section.
    where_first_under = first_part.find('_')
    where_last_under = first_part.rfind('_')
    variable_name = first_part[where_first_under+1:where_last_under]
    variable_name = variable_name.replace('_', ' ')

    # Get the season which is after the last underscore in this section.
    season_name = first_part[where_last_under+1:]

    # Get the observations which is in the second section
    where_dot = second_part.find('.')
    obs_name = second_part[:where_dot]
    obs_name = obs_name.replace('_', ' ')
    return dict(field=variable_name, season=season_name,
                obs=obs_name, proc='mean')


def main(images_dir, control_dir, thumbs_dir, csv_dir, temp_dir,
         num_to_save=None):
    """
    This is the main part of the validation note code

    images_dir = (str) a directory where the final images will go
    control_dir = (str) a directory that contains the control files
    csv_dir = (str) the directory that contains the csv files of rms errors
    temp_dir = (str) the directory used for temporary files
    num_to_save = (int) plot number to save as a pickle file for use in
                        test_colours.py
    """

    # Set up global variables
    global valorder_dict, source_dict, obs_dict, title_dict, ax1, title_type
    globalvar.valnote = True

    # Read the control files
    print 'Reading control files'
    valorder_dict = vm.read_info_file(globalvar.valorder_file)
    globalvar.plots_dict = vm.read_control_file('plots_file.dat',
                                                ctldir=control_dir)
    globalvar.item_dict = vm.read_control_file('item_file.dat',
                                               ctldir=control_dir)
    title_dict = vm.read_control_file('title_file.dat', ctldir=control_dir)
    globalvar.levels_dict = vm.read_control_file('levels_file.dat',
                                                 ctldir=control_dir)

    # Get the names of the experiment and control keys for the RMS table.
    # This is a bit of a fudge and needs sorting out.
    exper_key = source_dict.keys()[0]
    control_key = source_dict.keys()[1]
    for key in source_dict.keys():
        if key[0] == 'e':
            exper_key = key
        if key[0] == 'c':
            control_key = key

    # Initialise the rms list
    rms_list = rms.start(exper=source_dict[exper_key]['jobid'][0],
                         control=source_dict[control_key]['jobid'][0])

    # Loop over the pages
    for page_num, page_title in iter(sorted(valorder_dict.iteritems())):

        print 'Plotting page '+page_num+' '+page_title
        globalvar.error = False

        # Save the page title (for use in test_colours.py)
        if 'num_to_save' in locals():
            if int(page_num) == num_to_save:
                with open(os.path.join(temp_dir, 'page_title.txt'), 'w') as f:
                    f.write(page_title)

        # Collect together all the components on the page
        component_list = []
        equation_dict = {}
        for key in sorted(globalvar.plots_dict[page_title].keys()):
            if len(key) == 1:

                # How many equations are there and put this list in a
                #  dictionary
                n_equations = len(globalvar.plots_dict[page_title][key])-1
                equation_dict[key] = \
                    globalvar.plots_dict[page_title][key][0:n_equations]

                for equation in globalvar.plots_dict[page_title][key][0:n_equations]:

                    # Get the components of the equation and put them in a list
                    for comp in vm.split_equation(equation):
                        component_list.append(comp)

        # Extract the data from all_data_dict to data_dict
        component_set = set(component_list)
        data_dict = get_data(component_set)

        # Check to make sure there have not been any errors so far
        #  before continuing.
        if not globalvar.error:

            title_type = '4up'
            if 'type' in globalvar.plots_dict[page_title]:
                title_type = globalvar.plots_dict[page_title]['type'][0]

            # Set up the filename
            filename, thumb_file = vm.file_name(page_num, page_title)

            # What letter do we start at for calculating rms values.
            if 'rms' in globalvar.plots_dict[page_title]:
                globalvar.rms_first_letter = \
                    globalvar.plots_dict[page_title]['rms'][0][0]

            # Loop through all the plots on the page
            for key in sorted(equation_dict.keys()):

                # See if we want to calculate rms values
                calc_rms = False
                if 'rms' in globalvar.plots_dict[page_title]:
                    if key in globalvar.plots_dict[page_title]['rms'][0]:
                        calc_rms = True

                # Setup some variables for use in the equations
                n_equations = len(equation_dict[key])
                toplot_cube = range(n_equations)
                rms_float = range(n_equations)
                same_levels = True
                level_type = globalvar.plots_dict[page_title][key][-1]
                if 'same_levels' in globalvar.levels_dict[level_type]:
                    if globalvar.levels_dict[level_type]['same_levels'][0].lower() == 'false':
                        same_levels = False

                # Perform the equations
                for eq_num in range(n_equations):
                    toplot_cube[eq_num], rms_float[eq_num] = \
                        vm.perform_equation(equation_dict[key][eq_num],
                                            data_dict,
                                            page_title,
                                            key,
                                            rms_list=rms_list,
                                            filename=filename,
                                            calc_rms=calc_rms,
                                            same_levels=same_levels)
                if globalvar.error:
                    continue

                # Save to pickle file (for use in test_colours.py)
                if 'num_to_save' in locals():
                    if int(page_num) == num_to_save:

                        # Define your pickle file
                        pickle_file = os.path.join(temp_dir, key+'.pkl')

                        # Save to pickle file
                        print 'Saving plots for page ' + page_title + \
                            ' into file '+pickle_file
                        vm.save_pickle(toplot_cube, pickle_file)

                # What is the plot type
                globalvar.plot_type = vm.plot_type_test(toplot_cube[0])

                # If this is the first plot then define the page
                if key == 'a':
                    define_page(page_title, num_plots=len(equation_dict))

                # Make the plots
                if globalvar.plot_type in ('y_pressure', 'y_height',
                                           'lat_lon', 'lat_lon_regional'):
                    contour_plot(toplot_cube, location_dict[key], page_title,
                                 key, equation_dict[key], filename,
                                 rms_float=rms_float[0])
                elif globalvar.plot_type in ('x_longitude', 'x_latitude'):
                    line_plot(toplot_cube, location_dict[key], page_title,
                              key, equation_dict[key], filename,
                              rms_float=rms_float[0])
                else:
                    break

            # Save the plot
            out_file = os.path.join(images_dir, filename)
            print 'Saving plot to file '+out_file
            # Using ImageMetaTag to add metadata to images
            #  img_converter=1 : Reduces file size
            #  img_tags=dbtags : tags to add to file
            #  db_file=XX : Database file to write to
            #  do_thumb=YY : Create thumbnails of size YY
            #  keep_open=True : Do NOT close figure
            dbtags = create_tags_from_filename(filename)
            imt.savefig(out_file, img_converter=1,
                        db_file=os.path.join(images_dir, 'img.db'),
                        img_tags=dbtags, do_thumb=461, keep_open=True)
            ax1.clear()

            # Convert to thumbnail
            thumbs_file = os.path.join(thumbs_dir, thumb_file)
            cmd = 'convert -resize 476x346 %s %s' % (out_file, thumbs_file)
            os.system(cmd)

    print 'Finished producing plots. Outputting csv files...'

    # Output the csv files
    rms.end(rms_list, csv_dir)
    print 'Finished outputting csv files.'

    # Print an error if there are any
    if globalvar.error:
        print 'WARNING: The make_plots.py program encountered errors.'
        print '         Scroll up to see them.'
