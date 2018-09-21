"""
Port for ESMValTool v2 from v1.

Uses: ESMValTool v2, Python 3.x
Valeriu Predoi, UREAD, July 2018

Functionality: computes root mean squares for a bunch of geographical
regions;

Original docstring:
This file contains an rms class which contains all the information needed
to make a rms table of a particular region.
"""

import os
import math
import logging
import numpy.ma as ma
import iris
from esmvaltool.diag_scripts.autoassess._valmod_radiation import area_avg

logger = logging.getLogger(os.path.basename(__file__))


class RMSLISTCLASS(list):
    """
    Construct the regions class.

    This is the class for a list of RMSCLASS (i.e. for lots of regions).
    """

    def __init__(self, *args):
        """Init."""
        if not args:
            super(RMSLISTCLASS, self).__init__()
        else:
            super(RMSLISTCLASS, self).__init__(args[0])

    def __repr__(self):
        """Repr."""
        rms_out = "["
        for rms_item in self:
            rms_out += "rms.RMSCLASS for " + rms_item.region + ", \n"
        rms_out = rms_out[0:-3]
        rms_out += "]"

        return rms_out

    def __call__(self, region=False):
        """Call."""
        rms_found = False
        region_list = []
        for rms_item in self:
            region_list.append(rms_item.region)
            if region:
                if rms_item.region == region:
                    rms_returned = rms_item
                    rms_found = True
        if not region:
            logger.warning(
                "Please supply a region using the region='xxx' input. " +
                "Available regions are:")
        elif not rms_found:
            logger.warning("ERROR: Requested region not found.")
        if not rms_found:
            logger.error(region_list)
            raise Exception
        return rms_returned


# This is the class for one set of rms values (i.e. for one region)
class RMSCLASS:
    """Class per region."""

    def __init__(self, region, exper='experiment', control='control'):
        """
        Create instances of this class but also start making.

        html files that will contain all the rms data. (old)

        region = the region name
        exper = experiment jobid
        control = control jobid.
        """
        # Store the region name, experiment and control
        self.region = region
        self.exper = exper
        self.control = control

        # This could be a dictionary in the future; not now tho

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

    # Allow iterations over this
    def __iter__(self):
        """Iter."""
        return self

    # This defines how this class is shown on the screen if you print it
    def __repr__(self):
        """Repr."""
        rms_out = "rms.RMSCLASS for {0}".format(self.region)
        return rms_out

    def calc(self, toplot_cube, mask_cube):
        """Calculate the rms value of a cube for this region.

        toplot_cube = (cube) cube that is to be plotted
        mask_cube = (cube) the mask to be applied (land/sea)
        """
        # Make a copy of the input cube
        working_cube = toplot_cube.copy()

        # What type of plot is this
        plot_type = 'lat_lon'
        if not toplot_cube.coords(axis='x'):
            plot_type = 'zonal_mean'
        else:
            if len(toplot_cube.coords(axis='x')[0].points) == 1:
                plot_type = 'zonal_mean'
        if not toplot_cube.coords(axis='y'):
            plot_type = 'meridional_mean'
        else:
            if len(toplot_cube.coords(axis='y')[0].points) == 1:
                plot_type = 'meridional_mean'

        # Apply the mask but only for lat_lon plots
        if hasattr(self, 'mask_end'):
            if plot_type == 'lat_lon':
                # Apply the mask
                working_cube.data = \
                    ma.masked_array(working_cube.data,
                                    mask=(mask_cube.data > 0.5))

            else:
                # If there is a mask but we are using zonal
                # mean or meridional mean, return missing
                return 1e+20

        # Extract a region
        if hasattr(self, 'region_bounds'):

            # Extract just the latitudes you want
            lonc = iris.Constraint()
            latc = iris.Constraint()
            if plot_type == 'lat_lon' or plot_type == 'meridional_mean':
                lonc = iris.Constraint(
                    longitude=lambda lon:
                    self.region_bounds[0] <= lon <= self.region_bounds[2]
                )
            if plot_type == 'lat_lon' or plot_type == 'zonal_mean':
                latc = iris.Constraint(
                    latitude=lambda lat:
                    self.region_bounds[1] <= lat <= self.region_bounds[3]
                )
            working_cube = working_cube.extract(lonc & latc)

        # Check to see if we have any data left.
        # If not then apply a missing data number.
        amount_of_data = len(working_cube.data)
        if hasattr(working_cube.data, 'compressed'):
            amount_of_data = len(working_cube.data.compressed())
        if amount_of_data == 0:
            rms_float = 1e+20
        else:
            logger.info('Calculating RMS for %s', self.region)

            # Square the values
            squared_cube = working_cube**2

            # Mean the values
            area_average = area_avg(
                squared_cube, coord1='latitude', coord2='longitude')

            # Square root the answer
            rms_float = math.sqrt(area_average.data)

        return rms_float

    def calc_wrapper(self, toplot_cube, mask_cube, page_title):
        """
        Get the RMS value and adds it to its own data array.

        toplot_cube = (cube) cube that is to be plotted
        mask_cube = (cube) mask land/sea
        page_title = (str) the page title for this plot
        """
        rms_float = self.calc(toplot_cube, mask_cube)
        self.data_dict[page_title] = []
        if rms_float:
            self.data_dict[page_title].append(rms_float)
        return rms_float

    def tofile(self, csv_dir):
        """Output all the RMS statistics to csv files."""
        csv_file = 'summary_' + self.region + '_RMS_' + self.exper + '.csv'
        csv_path = os.path.join(csv_dir, csv_file)
        with open(csv_path, 'a') as out_file:
            for page_title, rms_list in self.data_dict.items():
                out_file.write('{0}: '.format(page_title))
                for rms_val in rms_list:
                    out_file.write('{0}'.format(str(rms_val)))
                    out_file.write('\n')


def start(exper='experiment', control='control'):
    """
    Make some instances of the rms class.

    exper = experiment jobid (optional)
    control = control jobid (optional).
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
    region_list = [
        'global', 'north', 'south', 'tropical_land', 'tropical_ocean',
        'east_asia', 'natl_europe', 'australia_land', 'sahara_n_africa',
        'tropical_n_africa', 'east_africa', 'central_africa',
        'southern_africa', 'africa_land'
    ]

    # Make a blank list that will hold the rms classes
    rms_list = RMSLISTCLASS()

    # Make the rms classes. This will also start making the summary web pages.
    for region in region_list:
        rms_list.append(RMSCLASS(region, exper=exper, control=control))

    return rms_list


def calc_all(rms_list, toplot_cube, mask_cube, page_title):
    """
    Loop through all the regions.

    Calculate rms values and store them in the class.
    rms_list = list of rms classes that stores all the information to do
               with the rms regions and the resulting answers.
    toplot_cube = (cube) cube that is to be plotted
    page_title = (str) the page title for this plot.
    """
    # Run through the loop, calculating rms values for each region
    rms_float_list = []
    n_rms = len(rms_list)
    for i in range(n_rms):
        rms_float = rms_list[i].calc_wrapper(toplot_cube, mask_cube,
                                             page_title)
        rms_float_list.append(rms_float)

    # Return the global rms value
    return rms_float_list[0]


def end(rms_list, csv_dir):
    """
    Finish using the rms class.

    rms_list = list of rms classes that stores all the information to do with
               the rms regions and the resulting answers.
    """
    for rms_instance in rms_list:
        rms_instance.tofile(csv_dir)
