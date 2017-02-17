"""

"""

import ConfigParser
import os
import pdb
import sys
import projects
import numpy as np

from netCDF4 import Dataset


class ESMValProject(object):
    """
    class to easily retrieve informations for ESMVal project_info
    in alphabetical order (mixed general and diagnostics specific)
    """
    def __init__(self, project_info):
        """
        Parameters
        ----------
        project_info : dict
            project information like provided by the launcher
        """
        self.project_info = project_info
        self.firstime = True
        self.oldvar = ""
#        self.version = os.environ['0_ESMValTool_version']

    def _get_path_with_sep(self, p):
        """ ensure that a pathname has the path separator at the end """
        if p[-1] == os.sep:
            return p
        else:
            return p + os.sep

    def average_data(self, data, dim_index):
        """Returns the mean values over certain dimensions. """
        # Usually the input array is three dimensional (time, lats, lons)
        if   (type(dim_index) == int):
            means = np.zeros(data.shape[dim_index])
            for mean in xrange(len(means)):
                if (dim_index == 0):
                    means[mean] = np.mean(data[mean, :, :])
                if (dim_index == 1):
                    means[mean] = np.mean(data[:, mean, :])
                if (dim_index == 2):
                    means[mean] = np.mean(data[:, :, mean])
        elif (dim_index == 'monthly'):
            means = np.zeros(12)
            for month in xrange(12):
                means[month] = np.mean(data[month::12, :, :])
        elif (dim_index == 'annual'):
            new_shape = data.shape
            new_shape = list(new_shape)
            new_shape[0] = 12
            means = np.zeros((new_shape))
            for month in xrange(12):
                for lat in xrange(means.shape[1]):
                    for lon in xrange(means.shape[2]):
                        means[month, lat, lon] = np.mean(data[month::12, lat, lon])
        elif (dim_index == 'annual'):
            new_shape = data.shape
            new_shape = list(new_shape)
            new_shape[0] = 12
            means = np.zeros((new_shape))
            for month in xrange(12):
                for lat in xrange(means.shape[1]):
                    for lon in xrange(means.shape[2]):
                        means[month, lat, lon] = np.mean(data[month::12, lat, lon])
        return means

    def check_model_instances(self, first_set, second_set):
        """Checks that two model sets have same elemenents. """
        # Not pretty but effective: we do a loop for both sets

        for model in first_set:
            try:
                second_set[model]
            except KeyError:
                print("PY  ERROR: I am getting inconsistent model sets for " +
                      "precipitation and temperature")
                print("PY  ERROR: All models must have both variables present.")
                print("PY  ERROR: This error is caused by " + model)
                print("PY  ERROR: Stopping the script and exiting")
                sys.exit()

        for model in second_set:
            try:
                first_set[model]
            except KeyError:
                print("PY  ERROR: I am getting inconsistent model sets for " +
                      "precipitation and temperature")
                print("PY  ERROR: All models must have both variables present.")
                print("PY  ERROR: This error is caused by " + model)
                print("PY  ERROR: Stopping the script and exiting")
                sys.exit()

    def ensure_directory(self, path):
        """ Checks if a given directory exists and creates it if necessary. """
        if not (os.path.exists(path)):
            os.makedirs(path)

    def ensure_looping(self, array):
        """ Checks the first and last element of array,
        currently only done for lons. FIXME if you need lats as well. """
        array.flags.writeable = True
        if (array[0] > array[-1]):
            for i in xrange(len(array) - 1):
                if (array[i + 1] < array[i]):
                    array[i + 1] += 360
        return array

    def extract_seasonal_mean_values(self,
                                     modelconfig,
                                     data,
                                     experiment,
                                     season,
                                     monthly=False):
        """Returns the season specific mean values for each lat, lon from data.
        We assume the usual indexing of time, lat, lon"""
        season_key = experiment + '_season_' + season
        data_shape = data.shape

        if (season == 'annual'):
            # For annual season we merely copy the data
            masked_values = data
        else:
            # For a specific season we mask the undesired values
            season_months = modelconfig.get(season_key, 'season_months').split()
            mask = np.ones((data_shape))
            for month in season_months:
                month_loc = int(month) - 1
                mask[month_loc::12, :, :] = 0
            masked_values = np.ma.masked_array(data, mask)

        if (monthly):
            return masked_values

        mean_values = np.zeros((data_shape[1], data_shape[2]))
        # Seasonal mean values
        for lat in xrange(data_shape[1]):
            for lon in xrange(data_shape[2]):
                mean_values[lat, lon] = masked_values[:, lat, lon].mean()
        return mean_values

    def find_nearest_value(self, array, value):
        """ Finds the nearest value in an array. """
        return np.abs(array - value).argmin()

    def get_area_coordinates(self, modelconfig, experiment, area):
        """Returns the coordinates (lat/lon) of the area of interest. """
        config_file = self.get_configfile()
        area_key = experiment + '_' + area

        # Checking that all needed values are in place in the config file
        if not modelconfig.has_section(area_key):
            print("PY  ERROR: Terminating in function 'get_area_coordinates'")
            print("PY  ERROR: undefined area specification: '" + area + "'")
            print("PY  ERROR: missing section: '[" + area_key + "]'")
            print("PY  ERROR: check your configuration file: " + config_file)
            sys.exit()
        for key in ['lat_min', 'lat_max', 'lon_min', 'lon_max']:
            if not modelconfig.has_option(area_key, key):
                print("PY  ERROR: Terminating in function 'get_area_coordinates'")
                print("PY  ERROR: undefined option " + key
                      + " in section '[" + area_key + "]'")
                print("PY  ERROR: check your configuration file: "
                      + config_file)
                sys.exit()
        lat_min = modelconfig.getint(area_key,  'lat_min')
        lat_max = modelconfig.getint(area_key,  'lat_max')
        lon_min = modelconfig.getint(area_key,  'lon_min')
        lon_max = modelconfig.getint(area_key,  'lon_max')
        return lat_min, lat_max, lon_min, lon_max

    def get_array_indices(self, lats, lons, coords):
        """Returns array indices of nearest values for given coordinates. """
        lat_min = self.find_nearest_value(lats, coords[0])
        lat_max = self.find_nearest_value(lats, coords[1])
        lon_min = self.find_nearest_value(lons, coords[2])
        lon_max = self.find_nearest_value(lons, coords[3])
        return lat_min, lat_max, lon_min, lon_max

    def get_clim_dir(self):
        return self._get_path_with_sep(self.project_info['GLOBAL']['climo_dir'])

    def get_clim_model_filenames(self, variable=None, monthly=True):
        """
        get names of models

        # TODO this is currently hardcoded for CMIP5 and should be more flexible
        perhaps one can even replace this code totally with code in the
        interface dir

        Parameters
        ----------
        climatology : bool
            return climatological mean values
        monthly : bool
            returned values are monthly (e.g. monthly=True)
            returns 12 fields of data
        """
        import projects
        if variable is None:
            raise ValueError('You need to specify a variable!')

        res = {}
        for currDiag in self.project_info['DIAGNOSTICS']:
            variables = currDiag.get_variables()
            field_types = currDiag.get_field_types()
            mip = currDiag.get_var_attr_mip()
            exp = currDiag.get_var_attr_exp()
            for idx in range(len(variables)):
                for model in self.project_info['MODELS']:
                    currProject = getattr(vars()['projects'], model.split_entries()[0])()
                    fullpath = currProject.get_cf_fullpath(self.project_info,
                                                           model,
                                                           field_types[idx],
                                                           variables[idx],
                                                           mip[idx],
                                                           exp[idx])
                    if variable == variables[idx] and os.path.isfile(fullpath):
                        if monthly:
                            if 'M' in field_types[idx]:
                                name = currProject.get_model_name(model)
                                res.update({name: fullpath})
        return res

    def get_clim_model_and_obs_filenames(self, variable=None):
        """Returns variable specific model and observation filenames with full
        paths. Also checks existance. """
        if variable is None:
            raise ValueError('You need to specify a variable!')
            sys.exit()

        obs = ''
        # The input variable is the actual variable but the obs id may differ
        if   (variable == 'pr' or variable == 'pr-mmday'):
            obs_id = 'pr_obs'
        elif (variable == 'ts'):
            obs_id = 'ts_obs'
        elif (variable == 'ua' or variable == 'ua-1000'):
            obs_id = 'ua_obs'
        else:
            obs_id = 'obs'

        # Check which model is the observations and extract that from the rest
        model_filenames = self.get_clim_model_filenames(variable=variable)
        for model in model_filenames:
            if (obs_id == self.get_model_id(model)
                    or variable + '_obs' == self.get_model_id(model)):
                obs = model

        if (len(obs) == 0):
            print("PY  ERROR: I didn't find observations for: '" + variable + "'")
            print("PY  ERROR: You should use id tags such as 'pr_obs', "
                  + "'ts_obs' or 'obs' for observations in the namelist file.")
            print("PY  ERROR: Stopping the script and exiting")
            sys.exit()

        # obs is the name of the obs model, obs_file is it's location
        # model_filenames contain the actual models (and paths)
        obs_file = model_filenames[obs]
        del model_filenames[obs]
        return obs, obs_file, model_filenames

    def get_configfile(self):
        """ returns the cfg file full location """
        currDiag = self.project_info['RUNTIME']['currDiag']
        cfg_file = currDiag.get_diag_script_cfg()
        return cfg_file

    def get_configfile_name(self):
        """ returns the cfg file filename in basic form """
        configfile_fullname = os.path.basename(self.get_configfile())
        configfile_name = os.path.splitext(configfile_fullname)[0]
        return configfile_name

    def get_currVars(self):
        """ returns the diagnostic variables """
        currDiag = self.project_info['RUNTIME']['currDiag']
        currVars = currDiag.get_variables()
        return currVars

    def get_diag_script_name(self):
        """ returns the current diagnostic script name """
        currDiag = self.project_info['RUNTIME']['currDiag']
        diag_script_name = os.path.splitext(currDiag.diag_script)[0]
        return diag_script_name

    def get_field_type(self):
        """ returns the first (and often only) field type """
        currDiag = self.project_info['RUNTIME']['currDiag']
        field_type = currDiag.get_field_types()[0]
        return field_type

    def get_graphic_format(self):
        """
        returns plotting directory path of project
        """
        return self.project_info['GLOBAL']['output_file_type']

    def get_model_data(self, modelconfig, experiment, area,
                       datakey, datafile, extend=''):
        """Extracts desired data for a specific area for a model from all
        model data. Also extends lat/lon coordinates (and therefore required
        values) if specified (these area "ghostlayers" for interpolation)."""
        lats = datafile.variables['lat'][:]
        lons = datafile.variables['lon'][:]

        # We take coordinates and transfer them into array indices
        coords = self.get_area_coordinates(modelconfig, experiment, area)
        indices = self.get_array_indices(lats, lons, coords)

        # Next define ghostlayers if needed
        # Latitudes extension assumes that you're not at poles
        # Longitutes extension works with all kinds of values
        # This is a nasty way to change tuple values but the easiest at hand
        lon_global = False
        if (extend == 'lats' or extend == 'both'):
            ilist = list(indices)
            ilist[0] += -1
            ilist[1] += 1
            indices = tuple(ilist)
        if (extend == 'lons' or extend == 'both'):
            # If we have global longtitude values then assign ghost layers later
            if   (indices[3] - indices[2] + 1 == len(lons)):
                lon_global = True
            else:
                ilist = list(indices)
                ilist[2] += -1
                ilist[3] += 1
                indices = tuple(ilist)

        # Next extract required lats and lons, the breaks are needed for
        # special cases when we have to loop over lats and/or lons without
        # taking all the values
        lat_break = len(lats) - indices[0]
        lon_break = len(lons) - indices[2]

        if   (indices[0] > indices[1] and indices[2] > indices[3]):
            # Both are looped
            lats_req = np.zeros(len(lats) - indices[0] + indices[1] + 1)
            lats_req[:lat_break] = lats[indices[0]:]
            lats_req[lat_break:] = lats[:indices[1] + 1]

            lons_req = np.zeros(len(lons) - indices[2] + indices[3] + 1)
            lons_req[:lon_break] = lons[indices[2]:]
            lons_req[lon_break:] = lons[:indices[3] + 1]

        elif (indices[0] > indices[1] and indices[2] <= indices[3]):
            # Only lats are looped
            lats_req = np.zeros(len(lats) - indices[0] + indices[1] + 1)
            lats_req[:lat_break] = lats[indices[0]:]
            lats_req[lat_break:] = lats[:indices[1] + 1]

            lons_req = lons[indices[2]:indices[3] + 1]

        elif (indices[0] <= indices[1] and indices[2] > indices[3]):
            # Only lons are looped
            lats_req = lats[indices[0]:indices[1] + 1]

            lons_req = np.zeros(len(lons) - indices[2] + indices[3] + 1)
            lons_req[:lon_break] = lons[indices[2]:]
            lons_req[lon_break:] = lons[:indices[3] + 1]

        else:
            # Neither are looped
            lats_req = lats[indices[0]:indices[1] + 1]
            lons_req = lons[indices[2]:indices[3] + 1]

        # Read in the required data for a specific model.
        # We have to define some special cases when there's looping over
        # latitudes or longtitudes.
        data = np.zeros((datafile.variables[datakey].shape[0],
                         lats_req.shape[0], lons_req.shape[0]))

        if (indices[0] > indices[1] and indices[2] > indices[3]):
            # Extract the four corners from the nc file
            data[:, :lat_break, :lon_break] = \
                datafile.variables[datakey][:, indices[0]:, indices[2]:]

            data[:, :lat_break, lon_break:] = \
                datafile.variables[datakey][:, indices[0]:, :indices[3] + 1]

            data[:, lat_break:, :lon_break] = \
                datafile.variables[datakey][:, :indices[1] + 1, indices[2]:]

            data[:, lat_break:, lon_break:] = \
                datafile.variables[datakey][:, :indices[1] + 1, :indices[3] + 1]

        elif (indices[0] > indices[1] and indices[2] <= indices[3]):
            # Looping only over lats
            data[:, :lat_break, :] = \
                datafile.variables[datakey][:, indices[0]:,
                                            indices[2]:indices[3] + 1]

            data[:, lat_break:, :] = \
                datafile.variables[datakey][:, :indices[1] + 1,
                                            indices[2]:indices[3] + 1]

        elif (indices[0] <= indices[1] and indices[2] > indices[3]):
            # Looping only over lons
            data[:, :, :lon_break] = \
                datafile.variables[datakey][:, indices[0]:indices[1] + 1,
                                            indices[2]:]

            data[:, :, lon_break:] = \
                datafile.variables[datakey][:, indices[0]:indices[1] + 1,
                                            :indices[3] + 1]
        else:
            data = datafile.variables[datakey][:, indices[0]:indices[1] + 1,
                                               indices[2]:indices[3] + 1]

        # Now we mask unwanted values if specified
        mask = modelconfig.getboolean('general', 'mask_unwanted_values')
        if (mask is True):
            # enquire what to do
            llow, lhigh = False, False
            if modelconfig.has_option('general', 'mask_limit_low'):
                llow = True
                low = modelconfig.getfloat('general', 'mask_limit_low')
            if modelconfig.has_option('general', 'mask_limit_high'):
                lhigh = True
                high = modelconfig.getfloat('general', 'mask_limit_high')
            # call the masking process
            if (llow and lhigh):
                data = self.mask_unwanted_values(data, low=low, high=high)
            elif (lhigh):
                data = self.mask_unwanted_values(data, high=high)
            elif (llow):
                data = self.mask_unwanted_values(data, low=low)

        # Ghost layers for global longtitude values
        if (lon_global is True):
            newlons = np.zeros(len(lons_req) + 2)
            newdata = np.zeros((data.shape[0], data.shape[1], data.shape[2] + 2))
            newlons[0] = lons_req[-1] - 360
            newlons[1:-1] = lons_req
            newlons[-1] = lons_req[0] + 360
            newdata[:, :, 0] = data[:, :, -1]
            newdata[:, :, 1:-1] = data[:, :, :]
            newdata[:, :, -1] = data[:, :, 0]
            lons_req = newlons
            data = newdata

        # Specify what to return based on the experiment
        # Zonal means want latitudes as well
        if   (experiment == 'zonal_means'):
            return lats_req, data

        # Equatorial and Southern Hemisphere need all
        elif (experiment == 'equatorial'
                or experiment == 'SouthernHemisphere'):
            return lats_req, lons_req, data

        # The scatterplots are fine with just the data
        else:
            return data

    def get_model_id(self, inmodel):
        """ Returns the id tag of the model if defined and empty string if not.
        The input is just the model name."""
        model_id = ''
        models = self.project_info['MODELS']
        for model in models:
            # loop over all entries and take the one we're after
            if (inmodel == model.split_entries()[1]):
                if ("id" in model.attributes):
                    model_id = model.attributes["id"]
                else:
                    pass
            else:
                pass
        return model_id

    def get_model_plot_style(self, model):
        """ Returns the style in which to plot the model:
        color (rgb), dashes and linewidth
        Check TropicalVariability.py for an example on the use of dashes."""
        # First we check that the style file for python is in place

        style_file = './diag_scripts/lib/python/style.cfg'
        if os.path.isfile(style_file):
            styleconfig = ConfigParser.ConfigParser()
            styleconfig.read(style_file)
        else:
            print("PY  ERROR: I didn't find the file for plotting styles.")
            print("PY  ERROR: This should be in: " + style_file)
            print("PY  ERROR: Stopping the script and exiting")
            sys.exit()

        if styleconfig.has_section(model):
            r, g, b = styleconfig.get(model, 'RGB').split()
            dash = styleconfig.getint(model, 'dash')
            width = 0.7 * styleconfig.getint(model, 'thickness')
        else:
            print("PY  ERROR: Couldn't find style for model " + model)
            print("PY  ERROR: Using the type for unknown style")
            r, g, b = styleconfig.get('unknown', 'RGB').split()
            dash = styleconfig.getint('unknown', 'dash')
            width = 0.7 * styleconfig.getint('unknown', 'thickness')

        color = [np.float_(r), np.float_(g), np.float_(b)]

        # Specifying the dashes
        if (dash == 0):
            dashes = []
        elif (dash == 1):
            dashes = [1, 2]
        elif (dash == 2):
            dashes = [3, 2]
        elif (dash == 3):
            dashes = [5, 2]
        elif (dash == 4):
            dashes = [1, 2, 3, 2]
        elif (dash == 5):
            dashes = [1, 2, 5, 2]
        elif (dash == 6):
            dashes = [3, 2, 5, 2]
        elif (dash == 7):
            dashes = [1, 2, 3, 2, 3, 2]
        elif (dash == 8):
            dashes = [1, 2, 3, 2, 5, 2]
        elif (dash == 9):
            dashes = [8, 4, 2, 4, 2, 4]
        elif (dash == 10):
            dashes = [8, 4, 8, 4, 8, 4]
        elif (dash == 11):
            dashes = [1, 6, 3, 6, 1, 3]
        elif (dash == 12):
            dashes = [7, 4, 1, 1, 5, 1]
        else:
            dashes = []
            print("PY  warning: You should specify more dashes in function")
            print("PY  warning: get_model_plot_style, located in file")
            print("PY  warning: " + os.getcwd())

        return color, dashes, width

    def get_plot_dir(self):
        """
        returns plotting directory path of project
        """
        return self._get_path_with_sep(self.project_info['GLOBAL']['plot_dir'])

    def get_plot_output_filename(self,
                                 diag_name='',
                                 variable='',
                                 model='',
                                 specifier='',
                                 end_specifier=''):
        """ Returns the output filename with optional keywords. Variable should
        be current variable but with multivariable plots define your own.
        Use specifier to further separate plots from oneanother. If model is
        specified, this script searches for additional keywords [2]
        [1]: <diag_script_name>_<specifier>_<variable>_<field_type>
        [2]: <model name>_<MIP>_<scenario>_<ensemble>_<start_year>_<end_year>
        """

        if (len(diag_name) == 0):
            diag_name = self.get_configfile_name()
        field_type = self.get_field_type()
        separator = "_"

        if (len(specifier) > 0):
            base = separator.join([diag_name,
                                   variable,
                                   field_type,
                                   specifier])
        else:
            base = separator.join([diag_name,
                                   variable,
                                   field_type])
        if (len(end_specifier) > 0):
            base = separator.join([base,
                                   end_specifier])

        if (len(model) > 0):
            """ This part is harcoded for CMIP5 runs and IT IS NOT ERROR FREE.
            It's slightly hazardous so use at your own risk i.e. you shouldn't
            use two instances of same model (i.e. two separate time periods -
            this loops over the models and takes the last one that hits) """
            for key in self.project_info['MODELS']:
                if (key.split_entries()[1] == model):
                    project, name, MIP, scenario, ensemble, start, end = key.split_entries()[0:7]
            years = "-".join([start, end])
            base = separator.join([base,
                                   project,
                                   model,
                                   MIP,
                                   scenario,
                                   ensemble,
                                   years])

        base = ".".join([base,
                         self.get_graphic_format()])
        return base

    def get_ticks(self, maxticks, array1, array2=''):
        """Returns suitable ticks for plots - we only deal with integers here.
        Takes in maxticks and max two arrays. This is because we want it was
        coded for scatterplotting observations and values simultaneously. """
        # getting the int limits
        if (len(array2) > 0):
            start = int(np.min([array1[0], array2[0]]))
            end_float = np.max([array1[-1], array2[-1]])
        else:
            start = int(array1[0])
            end_float = array1[-1]
        if (int(end_float) >= end_float):
            end = int(end_float)
        else:
            end = int(end_float + 1)

        # Next we define the ticks - always at least 2 (start and end)
        # The first if is here for correct behavior when looping over lons

        if (start <= end):
            length = end - start
        else:
            length = 360 - start + end

        for i in xrange(maxticks - 1):
            num_ticks = maxticks - i
            if ((length) % (num_ticks - 1) == 0):
                ticks = (np.linspace(start, start + length, num_ticks)).astype(int)
                break

        return ticks

    def get_ticks_labels(self, ticks, direction):
        """Returns tick labels with S/N, W/E definition with EQ if present. """
        labels = ticks.tolist()
        directions = []
        if (direction == 'lats'):
            directions = ['S', 'EQ', 'N']
            for item in xrange(len(ticks)):
                tick = ticks[item]
                if   (-90 <= tick < 0):
                    labels[item] = str(-tick) + directions[0]
                elif (tick == 0):
                    labels[item] = directions[1]
                elif (0 < tick <= 90):
                    labels[item] = str(tick) + directions[2]

        elif (direction == 'lons'):
            directions = ['W', 'EQ', 'E']
            for item in xrange(len(ticks)):
                tick = ticks[item]
                if (tick == 0 or tick == 360):
                    labels[item] = str(0)
                elif (tick == 180 or tick == 540):
                    labels[item] = str(180)
                elif (0 < tick < 180):
                    labels[item] = str(tick) + directions[2]
                elif (180 < tick < 360):
                    labels[item] = str(360 - tick) + directions[0]
                elif (360 < tick < 540):
                    labels[item] = str(tick - 360) + directions[2]
                elif (540 < tick < 720):
                    labels[item] = str(360 - tick) + directions[0]

        elif (direction == 'monthly'):
            labels = ['J', 'F', 'M', 'A', 'M', 'J',
                      'J', 'A', 'S', 'O', 'N', 'D']
        else:
            print("PY  ERROR: Unsupported direction in: 'get_ticks_labels'")
            print("PY  ERROR: You're only getting regular labels.")
        return labels

    def get_title_basename(self, datakey):
        """Return the title basename based on datakey."""
        if   (datakey == 'clt'):
            title = "Total cloud cover"
        elif (datakey == 'clivi'):
            title = "Cloud ice path"
        elif (datakey == 'clwvi'):
            title = "Cloud liquid water path"
        elif (datakey == 'hfls'):
            title = "Surface latent heat flux"
        elif (datakey == 'hfss'):
            title = "Surface sensible heat flux"
        elif (datakey == 'rlut'):
            title = "TOA outgoing longwave radiation"
        elif (datakey == 'rlutcs'):
            title = "TOA clear-sky outgoing longwave radiation"
        elif (datakey == 'rsut'):
            title = "TOA outgoing shortwave radiation"
        elif (datakey == 'rsutcs'):
            title = "TOA clear-sky outgoing shortwave radiation"
        elif (datakey == 'rlds'):
            title = "Surface incoming longwave radiation"
        elif (datakey == 'rldscs'):
            title = "Surface clear-sky incoming longwave radiation"
        elif (datakey == 'rsds'):
            title = "Surface incoming shortwave radiation"
        elif (datakey == 'rsdscs'):
            title = "Surface clear-sky incoming shortwave radiation"
        else:
            print("")
            print("PY  INFO: Unexpected datakey, no title selected.")
            print("PY  INFO: Add your own to function: 'get_title_basename'")
            print("PY  INFO: in diag_scripts/lib/python/esmval_lib.py")
            title = ""
        return title

    def get_verbosity(self):
        """ returns verbosity from namelist """
        return self.project_info['GLOBAL']['verbosity']

    def get_work_dir(self):
        """ returns the work directory """
        return self._get_path_with_sep(self.project_info['GLOBAL']['wrk_dir'])

    def mask_unwanted_values(self, array, low=None, high=None):
        """ Returns the given array masked with values outside the limits.
        Useful for throwing out missing values (i.e. sea surface temperature
        over land). Thus gives you a warning when converting missing values to
        nan (happens with masked arrays) - don't worry about it. """
        if (low is not None):
            out_array = np.ma.masked_less(array, low)
            if (high is not None):
                out_array = np.ma.masked_greater(array, high)
        elif (high is not None):
            out_array = np.ma.masked_greater(array, high)
        elif (low is None and high is None):
            print("PY  ERROR: You should specify at least one of the limits")
            print("PY  ERROR: low/high when using function:")
            print("PY  ERROR: 'mask_unwanted_values'")
            out_array = array
        return out_array

    def get_raw_inputfile(self):
        """
        get name of raw input file
        This is usefull if entire processing (including preprocessing)
        shall be done in the diag script itself
        """
        import projects


        res={}
        for currDiag in self.project_info['DIAGNOSTICS']:
            variables = currDiag.get_variables()
            field_types = currDiag.get_field_types()
            mip = currDiag.get_var_attr_mip()
            exp = currDiag.get_var_attr_exp()

            requested_vars = currDiag.get_variables_list()

            for idx in range(len(variables)):
                for model in self.project_info['MODELS']:
                        #~ print model.split_entries()[0] # gives 'JSBACH'
                        #~ print vars()
                        currProject = getattr(vars()['projects'], model.split_entries()[0])()
                        variable_defs_base_vars = currDiag.add_base_vars_fields(requested_vars, model)

                        base_vars = currDiag.select_base_vars(variable_defs_base_vars, model,
                                      currProject, self.project_info)
                        name = currProject.get_model_name(model)
                        res.update({name: {}})
                        for base_var in base_vars:

                            tmp_dir, tmp_files = currProject.get_cf_infile(self.project_info, model, base_var.fld, base_var.var, base_var.mip, base_var.exp)  # variable actually not used for e.g. JSBACH
                            #currProject.get_cf_infile(self.project_info, model, base_var, variable, mip, exp)  # variable actually not used for e.g. JSBACH

                            # store results
                            #return a dictionary by model and variable
                            res[name][base_var.var] =  {'directory' : tmp_dir,
                                                        'files' : tmp_files}
                            del tmp_dir, tmp_files

        return res

    # ##################################################################
    # write info for call to "write_references" (NCL) to temporary file,
    # then run NCL script "interface_scripts/write_references.ncl" which
    # will create the ESMValTool log-file
    # ##################################################################

    # Requires that "ncl" is in the search path.

    def write_references(self,
                         ref_script,
                         ref_auth,
                         ref_contr,
                         ref_diag,
                         ref_obs,
                         ref_proj,
                         project_info,
                         verbosity,
                         exit_on_warning):

        self.project_info['TEMPORARY'] = {}
        self.project_info['TEMPORARY']['ref_auth'] = ref_auth
        self.project_info['TEMPORARY']['ref_contr'] = ref_contr
        self.project_info['TEMPORARY']['ref_obs'] = ref_obs
        self.project_info['TEMPORARY']['ref_proj'] = ref_proj
        self.project_info['TEMPORARY']['ref_diag'] = ref_proj
        self.project_info['TEMPORARY']['ref_script'] = ref_script

        projects.run_executable("interface_scripts/write_references.ncl",
                                project_info,
                                verbosity,
                                exit_on_warning)

    # ###################################
    # write info on file read to log file
    # ###################################

    def add_to_filelist(self, filename):

        logfile = self.project_info['RUNTIME']['out_refs']
        log = open(logfile, "a+")

        f = Dataset(filename, 'r')

        path, fname = os.path.split(filename)

        try:
            var = f.variable
        except:
            var = ""

        try:
            mod = f.model
        except:
            mod = ""

        try:
            ver = f.version
        except:
            ver = "unknown"

        try:
            fix = f.fixfile
        except:
            fix = ""

        ref = []

        try:
            reftmp = f.reference
            for r in reftmp.split("\n"):
                ref.append(r)
        except:
            pass

        if (self.firstime is True):
            log.write("PREPROCESSING/REFORMATTING (ESMValTool v" + ver + ")\n\n")
            self.firstime = False

        if len(var) > 0:
            if (self.oldvar is not var):
                log.write("  Variable: " + var + "\n\n")
                self.oldvar = var

        if len(mod) > 0:
            log.write("    Model: " + mod + "\n")

        log.write("    Input path: " + path + "\n")
        log.write("    Input file(s):\n")
        log.write("      (1) " + fname + "\n")

        i = 0
        while i < 10000:
            strg = "infile_" + str(i).zfill(4)
            try:
                sfile = getattr(f, strg)
            except:
                sfile = ""
                break

            if len(sfile) > 0:
                try:
                    fs = Dataset(sfile, 'r')
                    try:
                        tid = fs.tracking_id
                    except:
                        tid = ""
                    try:
                        ref.append(fs.reference)
                    except:
                        pass
                    fs.close()
                except:
                    tid = ""
                    print("***** info: could not open original source file: " + sfile + " *****")
                    #pass

                if i is 0:
                    log.write("      Original source file(s) of all input file(s):\n")

                if len(tid) > 0:
                    log.write("        -S- (" + str(i+1) + ") " + sfile + " (tracking_id: " + tid + ")\n")
                else:
                    log.write("        -S- (" + str(i+1) + ") " + sfile + "\n")

            i = i + 1

        f.close()

        if len(fix) > 0:
            log.write("      Fixes applied to original source file(s): " + fix + "\n")

        i = 1
        for r in ref:
            if i == 1:
                log.write("    Reference(s) of original source file(s):\n")
            log.write("       (" + str(i) + ") " + r + "\n")
            i = i + 1

        log.write("\n")
        log.close()

