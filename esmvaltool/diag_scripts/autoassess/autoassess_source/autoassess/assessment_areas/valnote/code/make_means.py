'''
This module contains code to make the means table
(previously called the radiation table)
'''

import copy
import csv
import os
import pdb
import sys

import numpy

import iris

from markup import oneliner as line
import globalvar
import valmod as vm

import make_plots
import vnc


def area_mean(cube, analysis_method=iris.analysis.MEAN):
    """Make an area mean of a cube"""

    # Get the cube ready
    vm.get_cube_ready(cube)

    # If there is a model_level_number dimension then just use the first level
    # This is hard coded but a future improvement would be to use a level
    # as specified by the user in the means_file.dat
    dims = [cor.standard_name for cor in cube.dim_coords]
    if 'model_level_number' in dims:
        cube = cube.extract(iris.Constraint(model_level_number=1))

    # Generate weights for area averaging
    if not cube.coord('longitude').has_bounds():
        cube.coord('longitude').guess_bounds()
    if not cube.coord('latitude').has_bounds():
        cube.coord('latitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(cube)

    # Collapse the longitude and latitude dimensions
    area_average_cube = cube.collapsed(['longitude', 'latitude'],
                                       analysis_method,
                                       weights=grid_areas)

    if area_average_cube.ndim > 0:
        msg = ('ERROR: area_mean routine returned {0} dimensions when ' +
               ' there should be none.').format(area_average_cube.ndim)
        print msg
        if globalvar.debug:
            raise Exception(msg)
        globalvar.error = True
        return globalvar.missing_data
    else:
        return area_average_cube.data

def add_data_to_table(page, data):
    """ Adds data to a html table, replacing missing data with a question mark """
    
    if float(data) == globalvar.missing_data:
        page.td('?')
    else:
        page.td(data)


def csv2html(page_top, csv_dir=None):
    """Convert the csv file to a html file"""

    region_dict = {'global': 'Global means',
                   'land': 'Land only means',
                   'ocean': 'Ocean only means'}

    # Make a fileinfo list which contains a list of tables
    # and links to put near the top of the page
    fileinfo = list()
    for region in region_dict.keys():
        csvfile = os.path.join(csv_dir, 'means_'+region+'.csv')
        csv_file_end = os.path.split(csvfile)[1]
        filedict = {'region': region,
                    'linkname': 'Region: '+region,
                    'csv': csvfile,
                    'html': csv_file_end.replace('.csv', '.html')}
        fileinfo.append(filedict)

    # Loop over csvfiles, converting them all.
    for filedict in fileinfo:
        title = region_dict[filedict['region']] + ' table'
        print "Creating " + title + " html page"
        if not os.path.isfile(filedict['csv']):
            continue

        # Read the csv file
        columns = ['description', 'obs', 'control', 'exper']
        with open(filedict['csv'], 'r') as file:
            dictreader = csv.DictReader(file, columns)

            # Start the web page
            page = copy.deepcopy(page_top)
            page.init(title=title)

            # Make a row of links to the other means tables
            page.table()
            page.tr()
            for dfile in fileinfo:
                page.td(line.a(dfile['linkname'], href=dfile['html']))
            page.tr.close()
            page.table.close()

            # Make a title
            page.p()
            page.h2(title)

            # loop over rows in the csv file
            table_open = False
            for d in dictreader:
                if d['obs'].strip() == "obs":

                    # Close any old tables
                    if table_open:
                        page.table.close()

                    # Make a new table
                    page.p()
                    page.h3(d['description'])
                    page.table(border=1)
                    table_open = True

                    # Make a title row
                    page.tr()
                    page.td('Description')
                    page.td('Observations')
                    page.td('Control: '+d['control'])
                    page.td('Experiment: '+d['exper'])
                    page.tr.close()

                else:

                    # Make an information row
                    page.tr()
                    page.td(d['description'])
                    page.td(d['obs'])
                    add_data_to_table(page, d['control'])
                    add_data_to_table(page, d['exper'])
                    page.tr.close()

            # Close any old tables
            if table_open:
                page.table.close()

        # Output to html file
        with open(filedict['html'], 'w') as ofile:
            print >> ofile, page


def get_land_sea_frac(merged_key):
    """
    Get the land and sea fractions. Will return 4 values:
      land_frac, sea_frac, land_frac_mean, sea_frac_mean

    where:
    land_frac = land fraction cube on model grid
    sea_frac = sea fraction cube on model grid
    land_frac_mean = area weighted mean of the land fration
    sea_frac_mean = area weighted mean of the sea fration

    inputs:
    merged_key = (string) one of:
        'control_annm'
        'exper_annm'
    """

    # Try to get the land fraction from the supermeans
    field_constraints = [iris.AttributeConstraint(STASH=(01, 0, 505)),
                         iris.AttributeConstraint(STASH=(01, 3, 395))]
    land_frac = vm.extract_no_diurnal(make_plots.all_data_dict['control_annm'],
                                      all_constraints=field_constraints)

    if isinstance(land_frac, int):
        # Could not get the land fraction from the supermeans so
        #  therefore getting it from extra data
        land_frac = vm.get_land_frac(control_cube)

    # Calculate the sea fraction:
    # This should be 1.0 - land_frac but it has had to be re-written
    #  as iris equations have to start with a cube.
    sea_frac = land_frac*(-1.0) + 1.0

    # Calculate the fraction of total area
    land_frac_mean = float(area_mean(land_frac))
    sea_frac_mean = float(area_mean(sea_frac))

    return land_frac, sea_frac, land_frac_mean, sea_frac_mean


def calculate_mean_stash(item_list, all_data_cube, field_constraints,
                         region, data_frac_mask, data_frac_mean):

    """Calculate the mean values from input stash data"""

    # Extract the data you want
    data_cube = vm.extract_no_diurnal(all_data_cube,
                                      all_constraints=field_constraints)

    # If something has gone wrong then set to missing data
    if isinstance(data_cube, int):
        return globalvar.missing_data

    data_cube = vm.get_cube_ready(data_cube)

    # Set missing data to zero if required
    if item_list[6] == 'replace_missing' and \
       numpy.ma.is_masked(data_cube.data):
        newdata = numpy.where(data_cube.data.mask, 0,
                              data_cube.data.data)
        data_cube.data = newdata

    if region == 'land' or region == 'ocean':

        # Regrid the masks if needed
        n_lats = data_cube.coords(axis='y')[0].shape[0]
        if n_lats != data_frac_mask.coords(axis='y')[0].shape[0]:
            print 'Regridding land fraction for ' + item_list[0]

            data_frac_mask_temp = vm.regrid_cube(data_frac_mask, data_cube)
        else:
            data_frac_mask_temp = data_frac_mask

        # Mask out any land or ocean
        data_cube = data_cube * data_frac_mask_temp

    # Set up the analysis method
    if item_list[5].strip() == 'mean':
        analysis_method = iris.analysis.MEAN
        data_frac_mean_temp=data_frac_mean
    elif item_list[5].strip() == 'total':
        analysis_method = iris.analysis.SUM
        data_frac_mean_temp=1.0
    else:
        msg = 'WARNING: analysis method not supported for item ' + item_list[0]
        print msg
        if globalvar.debug:
            raise Exception(msg)
        print 'Using mean method'
        analysis_method = iris.analysis.MEAN
        data_frac_mean_temp=data_frac_mean
        globalvar.error = True

    # Do an area mean
    data_mean = float(area_mean(data_cube, analysis_method=analysis_method)) \
        / data_frac_mean_temp

    # Apply any multiplications
    data_mean = data_mean * float(item_list[4])

    return data_mean


def calculate_mean_equation(item_list, component_list, abbrev_data_dict):
    """Calculate means based on equations from previously calculated values"""

    equation_data = vm.expand_equation(item_list[2], component_list,
                                       dict_name='abbrev_data_dict')

    # Perform the equation
    try:
        data_mean = float(eval(equation_data))
    except:
        print 'ERROR executing equation ', equation_data
        print "ERROR: ", sys.exc_info()
        if globalvar.debug:
            pdb.set_trace()
        else:
            globalvar.error = True
            data_mean = globalvar.missing_data
    else:
        # Apply any multiplications
        data_mean = data_mean * float(item_list[4])

    return data_mean


def main(control_dir, csv_dir):

    print 'Making global means table'

    # Read the control files
    print 'Reading control files'
    means_dict = vm.read_control_file('means_file.dat', ctldir=control_dir)

    # Get the control and experiment jobids
    control_jobid = make_plots.source_dict['control']['jobid'][0]
    exper_jobid = make_plots.source_dict['exper']['jobid'][0]

    # Get the land and sea fractions
    (control_land_frac, control_sea_frac,
     control_land_frac_mean, control_sea_frac_mean) = \
        get_land_sea_frac('control_annm')
    (exper_land_frac, exper_sea_frac,
     exper_land_frac_mean, exper_sea_frac_mean) = \
        get_land_sea_frac('exper_annm')

    # Loop over the regions
    for region in ['global', 'land', 'ocean']:

        # Create a new file
        csv_file = csv_dir+'/means_'+region+'.csv'
        with open(csv_file, 'w') as f:

            # Initialise dictionaries to hold the abbreviations
            abbrev_control_dict = {}
            abbrev_control_dict['lc'] = 2.501e6
            abbrev_control_dict['lf'] = 0.334e6
            abbrev_exper_dict = {}
            abbrev_exper_dict['lc'] = 2.501e6
            abbrev_exper_dict['lf'] = 0.334e6

            # Create a fractional mask
            if region == 'global':
                control_frac_mask = None
                exper_frac_mask = None
                control_frac_mean = 1.0
                exper_frac_mean = 1.0
            if region == 'land':
                control_frac_mask = control_land_frac
                exper_frac_mask = exper_land_frac
                control_frac_mean = control_land_frac_mean
                exper_frac_mean = exper_land_frac_mean
            if region == 'ocean':
                control_frac_mask = control_sea_frac
                exper_frac_mask = exper_sea_frac
                control_frac_mean = control_sea_frac_mean
                exper_frac_mean = exper_sea_frac_mean

            # Loop over sections (in order)
            for section_key in sorted(means_dict.keys()):

                # Get the dictionary for this section
                section_dict = means_dict[section_key]

                # Add a line to the csv file for this section
                if 'section_title' in section_dict:
                    # Make a section title
                    section_title=section_dict['section_title'][0]
                    if region == 'land':
                        section_title = section_title.replace('Global', 'Land only')
                    if region == 'land':
                        section_title = section_title.replace('global', 'land only')
                    if region == 'ocean':
                        section_title = section_title.replace('Global', 'Ocean only')
                    if region == 'ocean':
                        section_title = section_title.replace('global', 'ocean only')
            
                    f.write('{0:45s},{1:10s},{2:10s},{3:10s}\n'.format(section_title, 'obs', control_jobid, exper_jobid))
                else:
                    msg = ('ERROR: no section_title detected for section {0}' +
                           ' of the means_file.dat').format(section_key)
                    print msg
                    if globalvar.debug:
                        raise Exception(msg)
                    globalvar.error = True
                    return globalvar.missing_data

                # Loop over items
                for item_key in sorted(section_dict.keys()):

                    # Get the list for this key
                    item_list = section_dict[item_key]

                    # Go to the next loop if we don't need to process this
                    if item_key in ['section_title']:
                        continue
                    if len(item_list) < 9:
                        print ('ERROR: means.file.dat needs to have 9 items ' +
                               'after the first key number for key number ' +
                               '{0} of section {1}'
                               ).format(item_key, section_dict['section_title'][0])
                    if item_list[8].strip() != region and \
                       item_list[8].strip() != 'all':
                        continue

                    print 'Processing {0}: {1} {2}'.format(section_dict['section_title'][0],
                                                           item_key,
                                                           item_list[0])
                    globalvar.error = False

                    # Start the next row of the csv file with the title for
                    #  this item
                    f.write('{0:45s},'.format(item_list[0]))

                    # Add the obs value for this
                    obs_value = ' '
                    if item_list[1].count(',') == 2:
                        # Found two commas therefore this string must contain
                        #  all three obs values
                        if region == 'global':
                            obs_value = item_list[1].split(',')[0]
                        if region == 'land':
                            obs_value = item_list[1].split(',')[1]
                        if region == 'ocean':
                            obs_value = item_list[1].split(',')[2]

                    else:
                        # If not two commas then obs value must be for global
                        #  mean
                        if region == 'global':
                            obs_value = item_list[1]
                    obs_value = obs_value.strip()
                    f.write('{0:10s},'.format(obs_value))

                    # Try to get the stash code you will load
                    from_stash = True
                    try:
                        constraint_int = int(item_list[2])
                    except:
                        from_stash = False

                    # Load the stash code
                    if from_stash:
                        section = constraint_int/1000
                        item = constraint_int-section*1000

                        # Setup the stash loading constraints
                        field_constraints = [iris.AttributeConstraint(STASH=(01, section, item)),
                                             iris.AttributeConstraint(STASH=(None, section, item))]

                        # Calculate the means
                        control_mean = calculate_mean_stash(item_list, make_plots.all_data_dict['control_annm'], field_constraints, region, control_frac_mask, control_frac_mean)
                        exper_mean = calculate_mean_stash(item_list, make_plots.all_data_dict['exper_annm'], field_constraints, region, exper_frac_mask, exper_frac_mean)

                    else:
                        print 'Executing equation ', item_list[2]

                        # Rearrange the equation so that it includes the
                        # abbrev_control_dict and abbrev_exper_dict
                        # dictionaries
                        component_list = vm.split_equation(item_list[2])

                        # Calculate the means
                        control_mean = calculate_mean_equation(item_list, component_list, abbrev_control_dict)
                        exper_mean = calculate_mean_equation(item_list, component_list, abbrev_exper_dict)

                    # Add this to the csv file
                    unformatted_string = '{{0:10.{0}f}},{{1:10.{0}f}}\n'.format(item_list[7])
                    f.write(unformatted_string.format(control_mean, exper_mean))

                    # If there is an abbreviation add this to the abbreviation
                    #  dictionaries for future use
                    if len(item_list[3].strip()) > 1:
                        abbrev_control_dict[item_list[3].strip()] = control_mean
                        abbrev_exper_dict[item_list[3].strip()] = exper_mean

            # Close the csv file
            print 'Finished making global means table. Outputting to CSV ' + \
                'file ' + csv_file
