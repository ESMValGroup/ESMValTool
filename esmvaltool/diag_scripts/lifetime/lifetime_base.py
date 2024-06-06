"""Base class for lifetime diagnostics."""

import logging
import os
import re

from copy import deepcopy
# import cartopy
import iris
import numpy as np
import dask.array as da
import matplotlib.pyplot as plt
import yaml
from iris.analysis import MEAN
from iris.util import broadcast_to_shape
from mapgenerator.plotting.timeseries import PlotSeries
from scipy.constants import g, N_A, R

from esmvaltool.diag_scripts.shared import ProvenanceLogger, names

# set standard molar masses
M_AIR = 28.970  # [g_air/mol_air]
M_H2O = 18.02   # [g_air/mol_air]

logger = logging.getLogger(__name__)

def create_press(var):
    """Create a pressure variable."""
    resolver = iris.common.resolve.Resolve(var, var)
    if var.coord('air_pressure').points.shape == var.shape:
        press = resolver.cube(var.coord('air_pressure').points)
    else:
        press = resolver.cube(broadcast_to_shape(
            var.coord('air_pressure').points,
            var.shape,
            var.coord_dims('air_pressure')
        ))
    press.long_name = 'air_pressure'
    press.var_name = 'air_pressure'
    press.units = 'Pa'

    return press


def calculate_gridmassdry(press, hus, z_coord):
    """Calculate the dry gridmass according to pressure and humidity.

    Used constants:
    g  standard acceleration of gravity from scipy
    """
    # calculate minimum in cubes
    if isinstance(press, iris.cube.Cube):
        pmin = da.nanmin(press.core_data())
    else:
        pmin = da.nanmin(press)

    # surface air pressure as cube
    pmax = press[:, 0, :, :].copy()
    pmax.data = press.coord('surface_air_pressure').points

    # delta pressure > kg m-1 s-2
    delta_p = dpres_plevel_4d(press,
                              pmin,
                              pmax,
                              z_coord=z_coord)
    # grid mass per square meter > kg m-2
    delta_m = delta_p / g
    # grid area > m2
    area = calculate_area(hus.coord('latitude'),
                          hus.coord('longitude'))
    # grid mass per grid > kg
    delta_gm = area * delta_m
    # to dry air
    gridmassdry = delta_gm * (1. - hus)
    # attach metadata
    gridmassdry.long_name = 'gridmassdry'
    gridmassdry.var_name = 'gridmassdry'
    gridmassdry.units = 'kg'

    return gridmassdry


def calculate_rho(variables):
    """Calculate number density (rho).

    Calculate number density with respect to the given input
    variables.
    """
    logger.info("Calculate number density (rho)")

    # model levels
    if 'grmassdry' in variables and 'grvol' in variables:
        rho = _number_density_dryair_by_grid(
            variables['grmassdry'],
            variables['grvol'])
    # pressure levels
    elif ('ta' in variables and
          'hus' in variables):
        rho = _number_density_dryair_by_press(
            variables['ta'],
            variables['hus'])
    else:
        raise NotImplementedError("The necessary variables"
                                  " to calculate number"
                                  " density of dry air"
                                  " are not provided.\n"
                                  "Provide either:\n"
                                  " - grmassdry and grvol\n"
                                  " or\n"
                                  " - ta and hus")

    return rho


def _number_density_dryair_by_press(temp, hus, press=None):
    """
    Calculate number density of dry air.

    Used to convert from mol / mol_dry into molec / cm3
    by using present temperature and humidity.

    ##qqq
    Should there be an option to provide the simulated pressure field
    rather than the derived pressure from interpolation?
    ###

    Used constants:
    - N_A    Avogrado constant from scipy
    - R      gas constant from scipy
    - M_AIR   Molarmass of Air
    - M_H2O   Molarmass of watervapor
    """
    logger.info('Calculate number density of dry air by pressure')

    if not press:
        logger.info('Pressure not given')
        press = create_press(temp)

    rho = N_A / 10.**6
    rho = rho * iris.analysis.maths.divide(press, R * temp)
    rho = rho * iris.analysis.maths.divide(1. - hus,
                                           1. + hus *
                                           (M_AIR
                                            / M_H2O - 1.))

    # correct metadata
    rho.var_name = 'rho'
    rho.units = 'cm-3'
    # [ 1 / cm^3 ]

    return rho


def _number_density_dryair_by_grid(grmassdry, grvol):
    """
    Calculate number density of dry air.

    Used to convert from mol / mol_dry into molec / cm3
    by using present gridmass of dry air and gridvolume.

    Since gridvolume might not be appropriate after interpolation
    to pressure coordinates, this version should only be used
    on model levels.

    Used constants:
    - N_A    Avogrado constant from scipy
    - M_AIR   Molarmass of Air
    """
    logger.info('Calculate number density of dry air by grid information')
    rho = ((grmassdry / grvol)
           * (N_A / M_AIR) * 10**(-3))  # [ 1 / cm^3 ]
    # correct metadata
    rho.var_name = 'rho'
    rho.units = 'cm-3'

    return rho


def guess_interfaces(coordinate):
    """
    Calculate the interfaces of a coordinate given its midpoints.

    Parameters
    ----------
    coordinate : iris.cube
        The coordinate array.

    Returns
    -------
    numpy.array
        An array with the lenght of the coordinate plus one.
    """
    interfaces = 0.5 * (coordinate[0:-1] + coordinate[1:])
    first = 0.5 * (3 * coordinate[0] - coordinate[1])
    last = 0.5 * (3 * coordinate[-1] - coordinate[-2])

    # Check limits
    # if coordinate.name.lower() in ['lat', 'latitude']:
    #     first = np.sign(first) * min([abs(first), 90.])
    #     last = np.sign(last) * min([abs(last), 90.])

    interfaces = da.insert(interfaces, 0, first)
    interfaces = da.append(interfaces, last)

    return interfaces


def calculate_area(latitude, longitude):
    """
    Calculate the area of each grid cell on a rectangular grid.

    Parameters
    ----------
    latitude : iris.cube.Coord
        The latitude coordinate of the grid.

    longitude : iris.cube.Coord
        The longitude coordinate of the grid.

    Returns
    -------
    iris.cube
        An array with the area (in m2) and the input latitude and longitude as
        coordinates.
    """
    r_earth = 6378100.0  # from astropy.constants

    lat_i = da.deg2rad(guess_interfaces(latitude.points))
    lon_i = da.deg2rad(guess_interfaces(longitude.points))

    delta_x = abs(lon_i[1:] - lon_i[:-1])
    delta_y = abs(da.sin(lat_i[1:]) - da.sin(lat_i[:-1]))

    output = da.outer(delta_y, delta_x) * r_earth**2
    output = output.astype('float32')

    result = iris.cube.Cube(output,
                            standard_name='cell_area',
                            long_name='cell_area',
                            var_name='cell_area',
                            units='m2',
                            dim_coords_and_dims=[(latitude, 0),
                                                 (longitude, 1)])

    return result


def dpres_plevel_4d(plev, pmin, pmax, z_coord='air_pressure'):
    """Calculate delta pressure levels.

    The delta pressure levels are based
    on the given pressure level as a
    four dimensional cube.

    """
    cubelist_dplev = [plev_slice.copy()
                      for plev_slice in plev.slices(['time',
                                                     'latitude',
                                                     'longitude'],
                                                    ordered=True)]
    cubelist_plev = [plev_slice.copy()
                     for plev_slice in plev.slices(['time',
                                                    'latitude',
                                                    'longitude'],
                                                   ordered=True)]

    increasing = (plev.coord(z_coord,
                             dim_coords=True).attributes['positive'] == 'down')
    last = plev.coords(z_coord)[0].shape[0] - 1

    for i, lev in enumerate(cubelist_plev):
        if increasing:
            increment = [i + 1, i - 1]
        else:
            increment = [i - 1, i + 1]

        if i == 0:
            cube = (cubelist_plev[increment[0]] - lev) / 2. + (lev - pmin)
            cubelist_dplev[i].data = cube.core_data()
        elif i == last:
            cube = (pmax - lev) + (lev - cubelist_plev[increment[1]]) / 2.
            cubelist_dplev[i].data = cube.core_data()
        else:
            cube = ((lev - cubelist_plev[increment[0]]) / 2.
                    + (cubelist_plev[increment[1]] - lev) / 2.)
            cubelist_dplev[i].data = cube.core_data()

        cubelist_dplev[i].add_aux_coord(iris.coords.AuxCoord(
            lev.coords(z_coord)[0].points[0]))

    dplev = iris.cube.CubeList(cubelist_dplev).merge_cube()
    dplev.transpose(new_order=[1, 0, 2, 3])
    dplev = iris.util.reverse(dplev, 1)
    return dplev


def dpres_plevel_1d(plev, pmin, pmax):
    """Calculate delta pressure levels.

    The delta pressure levels are based
    on the given pressure level coordinate
    as numpy array (one-dimensional).

    pmax can be set as float or as a multidimensional
    cube. The output of dpres_plevel will have the
    same dimensionality.
    """
    increasing = (da.diff(plev) >= 0).all()
    decreasing = (da.diff(plev) <= 0).all()

    if isinstance(pmax, float):

        dplev = plev.copy()

        for i, lev in enumerate(plev):
            if increasing:
                if lev == min(plev):
                    dplev[i] = (plev[i + 1] - lev) / 2. + (lev - pmin)
                elif lev == max(plev):
                    dplev[i] = (pmax - lev) + (lev - plev[i - 1]) / 2.
                else:
                    dplev[i] = ((plev[i + 1] - lev) / 2.
                                + (lev - plev[i - 1]) / 2.)
            elif decreasing:
                if lev == min(plev):
                    dplev[i] = (lev - pmin) + (plev[i - 1] - lev) / 2.
                elif lev == max(plev):
                    dplev[i] = (lev - plev[i + 1]) / 2. + (pmax - lev)
                else:
                    dplev[i] = ((lev - plev[i + 1]) / 2.
                                + (plev[i - 1] - lev) / 2.)

    elif isinstance(pmax, iris.cube.Cube):

        cubelist = [pmax.copy() for lev in plev]

        for i, lev in enumerate(plev):
            if increasing:
                if lev == min(plev):
                    cubelist[i].data = (lev - pmin) + (plev[i + 1] - lev) / 2.
                    cubelist[i] = (((lev - pmin)
                                   + (plev[i + 1] - lev) / 2.)
                                   * cubelist[i] / cubelist[i])
                elif lev == max(plev):
                    cubelist[i] = (pmax - lev) + (lev - plev[i - 1]) / 2.
                else:
                    cubelist[i] = (((plev[i + 1] - lev) / 2.
                                   + (lev - plev[i - 1]) / 2.)
                                   * cubelist[i] / cubelist[i])
            elif decreasing:
                if lev == min(plev):
                    cubelist[i] = (((lev - pmin)
                                   + (plev[i - 1] - lev) / 2.)
                                   * cubelist[i] / cubelist[i])
                elif lev == max(plev):
                    cubelist[i] = (lev - plev[i + 1]) / 2. + (pmax - lev)
                else:
                    cubelist[i] = (((lev - plev[i + 1]) / 2.
                                   + (plev[i - 1] - lev) / 2.)
                                   * cubelist[i] / cubelist[i])
            cubelist[i].units = 'Pa'
            cubelist[i].var_name = 'air_pressure'
            cubelist[i].add_aux_coord(iris.coords.AuxCoord(lev))
        # merge to single cube
        dplev = iris.cube.CubeList(cubelist).merge_cube()

    else:
        raise NotImplementedError("Function not implemented"
                                  f" for type {type(pmax)}")

    return dplev


def calculate_lifetime(dataset, plot_type, region):
    """Calculate the lifetime for the given plot_type and region."""

    # extract region from weights and reaction
    reaction = extract_region(dataset, region, case='reaction')
    weight = extract_region(dataset, region, case='weight')

    # calculate nominator and denominator
    # and sum of nominator and denominator via plot_type dimensions
    nominator = sum_up_to_plot_dimensions(weight, plot_type)
    denominator = sum_up_to_plot_dimensions(weight * reaction, plot_type)

    # division
    division = nominator / denominator

    return division


def extract_region(dataset, region, case='reaction'):
    """Return cube with everything outside region set to zero.


    Current aware regions:
    - TROP: troposphere (excl. tropopause)
    - STRA: stratosphere (incl. tropopause)
    """
    var = dataset[case]
    use_z_coord = dataset['use_z_coord']

    # mask regions outside
    if region in ['TROP', 'STRA']:

        z_4d = broadcast_to_shape(
            var.coord(use_z_coord).points,
            var.shape,
            var.coord_dims(use_z_coord)
        )

        tp_4d = broadcast_to_shape(
            dataset['tropopause'].data,
            var.shape,
            tuple(( var.coord_dims(item)[0]
                    for item in var.dim_coords
                    if not item == dataset['z_coord'] )),
        )

        if region == 'TROP':
            var.data = da.ma.masked_array(
                var.core_data(),
                mask=(z_4d <= tp_4d),
            )
        elif region == 'STRA':
            var.data = da.ma.masked_array(
                var.core_data(),
                mask=(z_4d > tp_4d),
            )
    else:
        raise NotImplementedError(f"region '{region}' is not supported")

    return var


def climatological_tropopause(cube):
    """Return cube with climatological tropopause pressure."""
    if not cube.coords('latitude', dim_coords=True):
        raise NotImplementedError("The provided cube must"
                                  " have a latitude cooridnate")

    tpp = (300. - 215. * (
        da.cos(da.deg2rad(cube.coord('latitude').points)) ** 2)) * 100.

    tp_clim = cube.copy()
    tp_clim.data = broadcast_to_shape(
        tpp,
        cube.shape,
        cube.coord_dims('latitude')
    )
    tp_clim.var_name = 'tp_clim'
    tp_clim.long_name = 'climatological tropopause pressure'
    tp_clim.units = 'Pa'

    return tp_clim


def sum_up_to_plot_dimensions(var, plot_type):
    """Return the cube summed over the appropriate dimensions."""
    if plot_type in ['timeseries', 'annual_cycle']:
        if var.coords('air_pressure', dim_coords=True):
            z_coord = var.coords('air_pressure', dim_coords=True)[0]
        elif var.coords('lev', dim_coords=True):
            z_coord = var.coords('lev', dim_coords=True)[0]
        elif var.coords('atmosphere_hybrid_sigma_pressure_coordinate',
                        dim_coords=True):
            z_coord = var.coords('atmosphere_hybrid_sigma_pressure_coordinate',
                                 dim_coords=True)[0]

    if plot_type == 'timeseries':
        cube = var.collapsed(['longitude', 'latitude', z_coord],
                             iris.analysis.SUM)
    elif plot_type == 'zonalmean':
        cube = var.collapsed(['longitude'], iris.analysis.SUM)
    elif plot_type == '1d_profile':
        cube = var.collapsed(['longitude', 'latitude'], iris.analysis.SUM)
    elif plot_type == 'annual_cycle':
        # TODO!
        # not use iris.analysis.SUM but some kind of mean
        # cube = var.collapsed(['longitude', 'latitude', z_coord],
        #                      iris.analysis.SUM)
        raise NotImplementedError("The sum to plot dimensions for plot_type"
                                  f" {plot_type} is currently not implemented")

    return cube


def calculate_reaction_rate(temp, reaction_type,
                            coeff_a, coeff_er, coeff_b=None):
    """Calculate the reaction rate.

    Calculated in Arrhenius form or in a given special form
    depending on the oxidation partner.
    """
    reaction_rate = deepcopy(temp)
    reaction_rate.units = 'unknown'

    # special reaction rate
    if coeff_b is not None:
        reaction_rate = coeff_a * iris.analysis.maths.exp(
            coeff_b * iris.analysis.maths.log(reaction_rate)
            - coeff_er / reaction_rate )
    else:
        # standard reaction rate (arrhenius)
        reaction_rate = coeff_a * iris.analysis.maths.exp(
            coeff_er / reaction_rate)

    # set units
    reaction_rate.units = 'cm3 s-1'
    reaction_rate.var_name = 'reaction_rate'
    reaction_rate.long_name = f'Reaction rate of {reaction_type}'
    return reaction_rate


def _replace_tags(paths, variable):
    """Replace tags in the config-developer's file with actual values."""
    if isinstance(paths, str):
        paths = set((paths.strip('/'), ))
    else:
        paths = set(path.strip('/') for path in paths)
    tlist = set()
    for path in paths:
        tlist = tlist.union(re.findall(r'{([^}]*)}', path))
    if 'sub_experiment' in variable:
        new_paths = []
        for path in paths:
            new_paths.extend(
                (re.sub(r'(\b{ensemble}\b)', r'{sub_experiment}-\1', path),
                 re.sub(r'({ensemble})', r'{sub_experiment}-\1', path)))
            tlist.add('sub_experiment')
        paths = new_paths

    for tag in tlist:
        original_tag = tag
        tag, _, _ = _get_caps_options(tag)

        if tag == 'latestversion':  # handled separately later
            continue
        if tag in variable:
            replacewith = variable[tag]
        else:
            raise ValueError(f"Dataset key '{tag}' must be specified for "
                             f"{variable}, check your recipe entry")
        paths = _replace_tag(paths, original_tag, replacewith)
    return paths


def _replace_tag(paths, tag, replacewith):
    """Replace tag by replacewith in paths."""
    _, lower, upper = _get_caps_options(tag)
    result = []
    if isinstance(replacewith, (list, tuple)):
        for item in replacewith:
            result.extend(_replace_tag(paths, tag, item))
    else:
        text = _apply_caps(str(replacewith), lower, upper)
        result.extend(p.replace('{' + tag + '}', text) for p in paths)
    return list(set(result))


def _get_caps_options(tag):
    lower = False
    upper = False
    if tag.endswith('.lower'):
        lower = True
        tag = tag[0:-6]
    elif tag.endswith('.upper'):
        upper = True
        tag = tag[0:-6]
    return tag, lower, upper


def _apply_caps(original, lower, upper):
    if lower:
        return original.lower()
    if upper:
        return original.upper()
    return original


class LifetimeBase():
    """Base class for lifetime diagnostic.

    It contains the common methods for path creation, provenance
    recording, option parsing and to create some common plots.

    """

    def __init__(self, config):
        self.cfg = config
        plot_folder = config.get(
            'plot_folder',
            '{plot_dir}/../../{dataset}/{exp}/{modeling_realm}/{real_name}',
        )
        plot_folder = plot_folder.replace('{plot_dir}',
                                          self.cfg[names.PLOT_DIR])
        self.plot_folder = os.path.abspath(plot_folder)
        self.plot_filename = config.get(
            'plot_filename',
            '{plot_type}_{real_name}_{dataset}_{mip}_{exp}_{ensemble}')
        self.plots = config.get('plots', {})
        default_config = os.path.join(os.path.dirname(__file__),
                                      "lifetime_config.yml")
        # cartopy_data_dir = config.get('cartopy_data_dir', None)
        # if cartopy_data_dir:
        #     cartopy.config['data_dir'] = cartopy_data_dir
        with open(config.get('config_file', default_config),
                  encoding='utf-8') as config_file:
            self.config = yaml.safe_load(config_file)

    def _add_file_extension(self, filename):
        """Add extension to plot filename."""
        return f"{filename}.{self.cfg['output_file_type']}"

    # def _get_proj_options(self, map_name):
    #     return self.config['maps'][map_name]

    def _get_variable_options(self, variable_group, map_name):
        options = self.config['variables'].get(
            variable_group, self.config['variables']['default'])
        if 'default' not in options:
            variable_options = options
        else:
            variable_options = options['default']
            if map_name in options:
                variable_options = {**variable_options, **options[map_name]}

        if 'bounds' in variable_options:
            if not isinstance(variable_options['bounds'], str):
                variable_options['bounds'] = [
                    float(n) for n in variable_options['bounds']
                ]
        logger.debug(variable_options)
        return variable_options

    def plot_timeseries(self, cube, var_info, period='', **kwargs):
        """Plot timeseries from a cube.

        It also automatically smoothes it for long timeseries of monthly data:
            - Between 10 and 70 years long, it also plots the 12-month rolling
              average along the raw series
            - For more than ten years, it plots the 12-month and 10-years
              rolling averages and not the raw series

        """
        if 'xlimits' not in kwargs:
            kwargs['xlimits'] = 'auto'
        length = cube.coord("year").points.max() - cube.coord(
            "year").points.min()
        filename = self.get_plot_path(f'timeseries{period}', var_info,
                                      add_ext=False)
        caption = ("{} of "
                   f"{var_info[names.LONG_NAME]} of dataset "
                   f"{var_info[names.DATASET]} (project "
                   f"{var_info[names.PROJECT]}) from "
                   f"{var_info[names.START_YEAR]} to "
                   f"{var_info[names.END_YEAR]}.")
        if length < 10 or length * 11 > cube.coord("year").shape[0]:
            self.plot_cube(cube, filename, **kwargs)
            self.record_plot_provenance(
                self._add_file_extension(filename),
                var_info,
                'timeseries',
                period=period,
                caption=caption.format("Time series"),
            )
        elif length < 70:
            self.plot_cube(cube, filename, **kwargs)
            self.record_plot_provenance(
                self._add_file_extension(filename),
                var_info,
                'timeseries',
                period=period,
                caption=caption.format("Time series"),
            )

            # Smoothed time series (12-month running mean)
            plt.gca().set_prop_cycle(None)
            self.plot_cube(cube.rolling_window('time', MEAN, 12),
                           f"{filename}_smoothed_12_months",
                           **kwargs)
            self.record_plot_provenance(
                self._add_file_extension(f"{filename}_smoothed_12_months"),
                var_info,
                'timeseries',
                period=period,
                caption=caption.format(
                    "Smoothed (12-months running mean) time series"),
            )
        else:
            # Smoothed time series (12-month running mean)
            self.plot_cube(cube.rolling_window('time', MEAN, 12),
                           f"{filename}_smoothed_12_months",
                           **kwargs)
            self.record_plot_provenance(
                self._add_file_extension(f"{filename}_smoothed_12_months"),
                var_info,
                'timeseries',
                period=period,
                caption=caption.format(
                    "Smoothed (12-months running mean) time series"),
            )

            # Smoothed time series (10-year running mean)
            self.plot_cube(cube.rolling_window('time', MEAN, 120),
                           f"{filename}_smoothed_10_years",
                           **kwargs)
            self.record_plot_provenance(
                self._add_file_extension(f"{filename}_smoothed_10_years"),
                var_info,
                'timeseries',
                period=period,
                caption=caption.format(
                    "Smoothed (10-years running mean) time series"),
            )

    def record_plot_provenance(self, filename, var_info, plot_type, **kwargs):
        """Write provenance info for a given file."""
        with ProvenanceLogger(self.cfg) as provenance_logger:
            prov = self.get_provenance_record(
                ancestor_files=[var_info['filename']],
                plot_type=plot_type,
                long_names=[var_info[names.LONG_NAME]],
                **kwargs,
            )
            provenance_logger.log(filename, prov)

    def plot_cube(self, cube, filename, linestyle='-', **kwargs):
        """Plot a timeseries from a cube.

        Supports multiplot layouts for cubes with extra dimensions
        `shape_id` or `region`.

        """
        plotter = PlotSeries()
        plotter.filefmt = self.cfg['output_file_type']
        plotter.img_template = filename
        region_coords = ('shape_id', 'region')

        for region_coord in region_coords:
            if cube.coords(region_coord):
                if cube.coord(region_coord).shape[0] > 1:
                    plotter.multiplot_cube(cube, 'time', region_coord,
                                           **kwargs)
                    return
        plotter.plot_cube(cube, 'time', linestyle=linestyle, **kwargs)

    @staticmethod
    def get_provenance_record(ancestor_files, **kwargs):
        """Create provenance record for the diagnostic data and plots."""
        record = {
            'authors': [
                'vegas-regidor_javier',
            ],
            'references': [
                'acknow_project',
            ],
            'ancestors': ancestor_files,
            **kwargs
        }
        return record

    def get_plot_path(self, plot_type, var_info, add_ext=True):
        """Get plot full path from variable info.

        Parameters
        ----------
        plot_type: str
            Name of the plot
        var_info: dict
            Variable information from ESMValTool
        add_ext: bool, optional (default: True)
            Add filename extension from configuration file.

        """
        return os.path.join(
            self.get_plot_folder(var_info),
            self.get_plot_name(plot_type, var_info, add_ext=add_ext),
        )

    def get_plot_folder(self, var_info):
        """Get plot storage folder from variable info.

        Parameters
        ----------
        var_info: dict
            Variable information from ESMValTool

        """
        info = {
            'real_name': self._real_name(var_info['variable_group']),
            **var_info
        }
        folder = os.path.expandvars(
            os.path.expanduser(
                list(_replace_tags(self.plot_folder, info))[0]
            )
        )
        if self.plot_folder.startswith('/'):
            folder = '/' + folder
        if not os.path.isdir(folder):
            os.makedirs(folder, exist_ok=True)
        return folder

    def get_plot_name(self, plot_type, var_info, add_ext=True):
        """Get plot filename from variable info.

        Parameters
        ----------
        plot_type: str
            Name of the plot
        var_info: dict
            Variable information from ESMValTool
        add_ext: bool, optional (default: True)
            Add filename extension from configuration file.

        """
        info = {
            "plot_type": plot_type,
            'real_name': self._real_name(var_info['variable_group']),
            **var_info
        }
        file_name = list(_replace_tags(self.plot_filename, info))[0]
        if add_ext:
            file_name = self._add_file_extension(file_name)
        return file_name

    @staticmethod
    def _set_rasterized(axes=None):
        """Rasterize all artists and collection of axes if desired."""
        if axes is None:
            axes = plt.gca()
        if not isinstance(axes, list):
            axes = [axes]
        for single_axes in axes:
            for artist in single_axes.artists:
                artist.set_rasterized(True)
            for collection in single_axes.collections:
                collection.set_rasterized(True)

    @staticmethod
    def _real_name(variable_group):
        for subfix in ('Ymean', 'Ysum', 'mean', 'sum'):
            if variable_group.endswith(subfix):
                variable_group = variable_group.replace(subfix, '')
        return variable_group
