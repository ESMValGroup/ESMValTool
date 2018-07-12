'''
This file contains an rms class which contains all the information needed
to make a rms table of a particular region.
'''

import csv
import math
import pdb

import numpy
import numpy.ma as ma
import iris


class rms_list_class(list):
    '''
    This is the class for a list of rms_class (i.e. for lots of regions)
    '''

    def __init__(self, *args, **kwargs):
        if len(args) == 0:
            super(rms_list_class, self).__init__()
        else:
            super(rms_list_class, self).__init__(args[0])

    def __repr__(self):
        # This defines how this class is shown on the screen if you print it
        # See http://stackoverflow.com/questions/1535327/python-how-to-print-a-class-or-objects-of-class-using-print
        rms_out = "["
        for rms_item in self:
            rms_out += "rms.rms_class for "+rms_item.region+", \n"
        rms_out = rms_out[0:-3]
        rms_out += "]"

        return rms_out

    def __call__(self, region=False):
        # This is what happens when this is called.
        # Return an rms item that matches the region requested.
        # http://www.rafekettler.com/magicmethods.html
        rms_found = False
        region_list = []
        for rms_item in self:
            region_list.append(rms_item.region)
            if region:
                if rms_item.region == region:
                    rms_returned = rms_item
                    rms_found = True
        if not region:
            print("Please supply a region using the region='xxxx' input. " + \
                "Available regions are:")
        elif not rms_found:
            print("ERROR: Requested region {0} not found. Available " + \
                "regions are:".format(region))
        if not rms_found:
            print(region_list)
            raise Exception
        return rms_returned


# This is the class for one set of rms values (i.e. for one region)
class rms_class:

    def __init__(self, region, exper='experiment', control='control'):
        """
        This not only creates instances of this class but also starts making
        html files that will contain all the rms data.

        region = the region name
        exper = experiment jobid
        control = control jobid
        """

        # Store the region name, experiment and control
        self.region = region
        self.exper = exper
        self.control = control

        # TODO: Store following information in dictionary

        # Store the region boundaries
        if region == 'north':
            self.region_bounds = [-180, 30, 180, 90]
        if region == 'south':
            self.region_bounds = [-180, -90, 180, -30]
        if region == 'tropical_land':
            self.region_bounds = [-180, -30, 180, 30]
        if region == 'tropical_ocean':
            self.region_bounds = [-180, -30, 180, 30]
        if region == 'east_asia':
            self.region_bounds = [100, 20, 145, 50]
        if region == 'natl_europe':
            self.region_bounds = [-45, 25, 60, 75]
        if region == 'australia_land':
            self.region_bounds = [100, -45, 155, -10]

        # Store the end of the mask key
        if region == 'tropical_land':
            self.mask_end = 'land_gt_50pc'
        if region == 'tropical_ocean':
            self.mask_end = 'ocean_gt_50pc'
        if region == 'australia_land':
            self.mask_end = 'land_gt_50pc'
        if region == 'sahara_n_africa':
            self.mask_end = 'SaharaNC'
        if region == 'tropical_n_africa':
            self.mask_end = 'TNAfrica'
        if region == 'east_africa':
            self.mask_end = 'EAfrica'
        if region == 'central_africa':
            self.mask_end = 'CAfrica'
        if region == 'southern_africa':
            self.mask_end = 'SAfrica'
        if region == 'africa_land':
            self.mask_end = 'Africa'

        # Make a blank dictionary to store the values
        self.data_dict = {}

        # What is the filename
        # filename = 'summary_'+region+'.csv'

        # Open a file
        # self.file_unit = open('/data/local2/hadco/temp/'+filename, 'w')
        # self.csv_writer = csv.writer(self.file_unit, delimiter=',')

        # Write out the header to the csv file
        # header=['Description',
        #         'RMS1 (expt:'+exper+' vs control:'+control+')',
        #         'RMS2 (control:'+control+' vs obs)',
        #         'RMS3 (expt:'+exper+' vs obs)']
        # self.csv_writer.writerow(header)

    # **************** Other class defined functions ****************

    # Allow iterations over this
    def __iter__(self):
        return(self)

    # This defines how this class is shown on the screen if you print it
    def __repr__(self):
        rms_out = "rms.rms_class for {0}".format(self.region)
        return rms_out

    # This defines what it returns if you call this class
    def __call__(self, description=None, letter=None):
        if letter:
            index_value = ord(letter)-ord(globalvar.rms_first_letter)+1
        if description and letter:
            return self.data_dict[description][index_value]
        elif description:
            rms_dict = {}
            for index in range(len(self.data_dict[description])-1):
                indexp1 = index+1
                this_letter = chr(index+ord(globalvar.rms_first_letter))
                rms_dict[this_letter] = self.data_dict[description][indexp1]
            return rms_dict
        elif letter:
            rms_dict = {}
            for description, rms_values in self.data_dict.iteritems():
                rms_dict[description] = rms_values[index_value]
            return rms_dict
        else:
            print("Please call rms.rms_class with either description or " + \
                  "letter inputs. e.g. my_rms_class(description='TOA " + \
                  "longwave vs EBAF', letter='a') or my_rms_class" + \
                  "(description='TOA longwave vs EBAF') or " + \
                  "my_rms_class(letter='a')")
            raise Exception

    def calc(self, toplot_cube):
        """Calculate the rms value of a cube for this region.

        toplot_cube = (cube) cube that is to be plotted
        """

        # Make a copy of the input cube
        working_cube = toplot_cube.copy()

        # What type of plot is this
        plot_type = 'lat_lon'
        if len(toplot_cube.coords(axis='x')) == 0:
            plot_type = 'zonal_mean'
        else:
            if len(toplot_cube.coords(axis='x')[0].points) == 1:
                plot_type = 'zonal_mean'
        if len(toplot_cube.coords(axis='y')) == 0:
            plot_type = 'meridional_mean'
        else:
            if len(toplot_cube.coords(axis='y')[0].points) == 1:
                plot_type = 'meridional_mean'

        # Apply the mask but only for lat_lon plots
        if hasattr(self, 'mask_end'):

            if plot_type == 'lat_lon':

                # What resultion mask should we be using
                #if working_cube.shape[1]/2 > 364:
                #    mask_key = 'n512_'+self.mask_end
                #elif working_cube.shape[1]/2 > 156:
                #    mask_key = 'n216_'+self.mask_end
                #elif working_cube.shape[1]/2 > 72:
                #    mask_key = 'n96_'+self.mask_end
                #else:
                #    mask_key = 'n48_'+self.mask_end

                # Extract the mask from the extra data dictionary
                #mask_cube = globalvar.extra_data_dict[mask_key]
                mask_cube = globalvar.extra_data_dict['1deg_1deg_landsea_frac']

                # Roll the mask
                mask_cube = vm.roll_cube(mask_cube)

                # Is the mask in need of regridding. If so regrid the mask.
                if mask_cube.shape != toplot_cube.shape:
                    if not mask_cube.coord(axis='x').has_bounds():
                        mask_cube.coord(axis='x').guess_bounds()
                    if not mask_cube.coord(axis='y').has_bounds():
                        mask_cube.coord(axis='y').guess_bounds()
                    mask_cube = vm.regrid_cube(mask_cube, working_cube,
                                               scheme_str='Linear')

                # Apply the mask
                try:
                    working_cube.data = \
                        ma.masked_array(working_cube.data,
                                        mask=(mask_cube.data > 0.5))
                except:
                    print('ERROR: Failed to assign mask')
                    if globalvar.debug:
                        pdb.set_trace()
                    else:
                        return globalvar.missing_data

            else:
                # If there is a mask but we are using zonal mean or meridional
                # mean data then just return missing data
                return globalvar.missing_data

        # Extract a region
        if hasattr(self, 'region_bounds'):

            # Extract just the latitudes you want
            lonc = iris.Constraint()
            latc = iris.Constraint()
            if plot_type == 'lat_lon' or plot_type == 'meridional_mean':
                lamfn = lambda lon: self.region_bounds[0] <= lon \
                    <= self.region_bounds[2]
                lonc = iris.Constraint(longitude=lamfn)
            if plot_type == 'lat_lon' or plot_type == 'zonal_mean':
                lamfn = lambda lat: self.region_bounds[1] <= lat \
                    <= self.region_bounds[3]
                latc = iris.Constraint(latitude=lamfn)
            working_cube = working_cube.extract(lonc & latc)

        # Check to see if we have any data left.
        # If not then apply a missing data number.
        amount_of_data = len(working_cube.data)
        if hasattr(working_cube.data, 'compressed'):
            amount_of_data = len(working_cube.data.compressed())
        if amount_of_data == 0:
            rms_float = globalvar.missing_data
        else:
            print('Calculating RMS for '+self.region)

            # Square the values
            squared_cube = working_cube ** 2

            # Mean the values
            area_average = vm.area_avg(squared_cube)

            # Square root the answer
            rms_float = math.sqrt(area_average)

        return rms_float

    def calc_wrapper(self, toplot_cube, letter, page_title,
                     filename='no plot available'):
        """
        Gets the RMS value and adds it to its own data array.

        toplot_cube = (cube) cube that is to be plotted
        letter = (str) letter number of this calculation. For validation
                       notes this is normally the plot letter. i.e.
             a = experiment
             b = experiment minus control
             c = control minus obs
             d = experiment minus obs
        page_title = (str) the page title for this plot
        filename = (str) the filename of the corresponding png plot.
        """

        rms_float = self.calc(toplot_cube)

        # If this is the first time you have calculated rms values then add
        #  a key to the dictionary. The first entry should be the filename.
        if letter == globalvar.rms_first_letter:
            self.data_dict[page_title] = []

        # Test to make sure that an rms entry exists for your page title.
        # This will fail if you have started making an RMS entry but not
        # using the first letter first.
        # remaining_letters='abcdefghijklmnop'.replace(globalvar.rms_first_letter, '')
        # if letter in remaining_letters:
        try:
            rms_len = len(self.data_dict[page_title])
        except:
            print('ERROR: Could not assign RMS value as incorrect first ' + \
                'letter used (in rms.py) in '+page_title)
            if globalvar.debug:
                pdb.set_trace()
            else:
                return globalvar.missing_data

        # Subsequent entries are the rms values
        self.data_dict[page_title].append(rms_float)
        return rms_float

    def add(self, letter, page_title, rms_float):
        filename = 'blank'
        if letter == globalvar.rms_first_letter:
            self.data_dict[page_title] = [filename]
        self.data_dict[page_title].append(rms_float)

    def tofile(self, csv_dir):
        """
        Output all the RMS statistics to csv files

        csv_dir = (str) the directory that contains the csv files of rms errors
        """

        # Create a new file
        #with open(csv_dir+'/summary_'+self.region+'.csv', 'w') as f:
        with open(csv_dir+'/summary_'+self.region+'_'+self.exper+'.csv', 'w') as f:

            '''
            # Add the heading row
            if globalvar.valnote:
                f.write('Description, Filename, RMS1 (expt:{1} vs control:{0}), RMS2 (control:{0} vs obs), RMS3 (expt:{1} vs obs)\n'.format(self.control, self.exper))
            else:
                f.write('Description, Filename, control:{0} vs obs, expt:{1} vs obs\n'.format(self.control, self.exper))
            '''

            # Loop over all the fields, outputting their values to the csv file
            for page_title, rms_list in self.data_dict.iteritems():
                f.write('{0},'.format(page_title))
                for index in range(len(rms_list)):
                    f.write('{0}'.format(rms_list[index]))
                    if index < len(rms_list)-1:
                        f.write(',')     # Next comma
                    else:
                        f.write('\n')    # Next line


def start(exper='experiment', control='control'):
    """
    Make some instances of the rms class

    exper = experiment jobid (optional)
    control = control jobid (optional)
    """

    # Loop over all regions. Regions are:
    # 0 = globe
    # 1 = north of 30N
    # 2 = south of 30S
    # 3 = tropical land
    # 4 = tropical ocean
    # 5 = east asia
    # 6 = north atlantic and europe
    # 7 = australian land
    # 8 = sahara and north african coast
    # 9 = tropical northern africa
    # 10 = east africa
    # 11 = central africa
    # 12 = southern africa
    # 13 = african land

    # Make a list of the regions
    region_list = ['global', 'north', 'south', 'tropical_land',
                   'tropical_ocean', 'east_asia', 'natl_europe',
                   'australia_land', 'sahara_n_africa', 'tropical_n_africa',
                   'east_africa', 'central_africa', 'southern_africa',
                   'africa_land']

    # Make a blank list that will hold the rms classes
    rms_list = rms_list_class()

    # Make the rms classes. This will also start making the summary web pages.
    for region in region_list:
        rms_list.append(rms_class(region, exper=exper, control=control))

    return rms_list


def calc_all(rms_list, toplot_cube, letter, page_title,
             filename='no plot available'):
    """
    Loop through all the regions, calculating rms values and storing
    them in the class.

    rms_list = list of rms classes that stores all the information to do
               with the rms regions and the resulting answers.
    toplot_cube = (cube) cube that is to be plotted
    letter = (str) letter number of this plot on the page.
             b = experiment minus control
             c = control minus obs
             d = experiment minus obs
    page_title = (str) the page title for this plot
    filename = (str) the filename of the png plot
    """

    # Run through the loop, calculating rms values for each region
    rms_float_list = []
    n_rms = len(rms_list)
    for i in range(n_rms):
        rms_float = rms_list[i].calc_wrapper(toplot_cube, letter, page_title,
                                             filename=filename)
        rms_float_list.append(rms_float)

    # Return the global rms value
    return rms_float_list[0]


def add_example(rms_list, letter, page_title):
    """Example code to add numbers to your rms_list"""
    n_rms = len(rms_list)
    rms_float = 1.2
    for i in range(n_rms):
        rms_list[i].add(letter, page_title, rms_float)
        rms_float = rms_float+1.0


def end(rms_list, csv_dir):
    """
    Finish using the rms class

    rms_list = list of rms classes that stores all the information to do with
               the rms regions and the resulting answers.
    """
    for rms_instance in rms_list:
        rms_instance.tofile(csv_dir)
