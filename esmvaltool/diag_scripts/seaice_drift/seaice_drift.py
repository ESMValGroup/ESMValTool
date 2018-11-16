"""Sea ice drift diagnostic"""
import os
import logging
import math
import csv
import calendar

import numpy as np
from scipy import stats
from matplotlib import pyplot as plt

import iris
import iris.cube
import iris.analysis
import iris.coords

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n

logger = logging.getLogger(os.path.basename(__file__))

MONTHS_PER_YEAR = 12

class SeaIceDrift(object):

    def __init__(self, config):
        self.cfg = config
        self.datasets = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.variables = esmvaltool.diag_scripts.shared.Variables(self.cfg)

        self.siconc = {}
        self.sivol = {}
        self.sispeed = {}
        self.lat_mask = {}

        self.slope_drift_siconc = {}
        self.intercept_drift_siconc = {}
        self.slope_ratio_drift_siconc = {}
        self.error_drift_siconc = {}

        self.slope_drift_sivol = {}
        self.intercept_drift_sivol = {}
        self.slope_ratio_drift_sivol = {}
        self.error_drift_sivol = {}

    def compute(self):
        logger.info('Loading sea ice concentration')
        siconc_original= {}
        siconc_files = self.datasets.get_path_list(
            standard_name='sea_ice_area_fraction')
        for filename in siconc_files:
            reference_dataset = self._get_reference_dataset(
                self.datasets.get_info('reference_dataset', filename)
            )
            alias = self._get_alias(filename, reference_dataset)
            siconc = self._load_cube(filename, 'sea_ice_area_fraction')
            siconc_original[alias] = siconc

            self.siconc[alias] = self._compute_mean(
                siconc, self._get_mask(siconc, filename)
            )


        logger.info('Loading sea ice thickness')
        sithick_files = self.datasets.get_path_list(
            standard_name='sea_ice_thickness')
        for filename in sithick_files:
            reference_dataset = self._get_reference_dataset(
                self.datasets.get_info('reference_dataset', filename)
            )
            alias = self._get_alias(filename, reference_dataset)
            sithick = self._load_cube(filename, 'sea_ice_thickness')
            self.sivol[alias] = self._compute_mean(
                sithick * siconc_original[alias],
                self._get_mask(sithick, filename)
            )
            del sithick

        logger.info('Load sea ice velocities')
        sispeed_files = self.datasets.get_path_list(
            standard_name='sea_ice_speed')
        for filename in sispeed_files:
            reference_dataset = self._get_reference_dataset(
                self.datasets.get_info('reference_dataset', filename)
            )
            alias = self._get_alias(filename, reference_dataset)
            sispeed = self._load_cube(filename, 'sea_ice_speed')
            self.sispeed[alias] = self._compute_mean(
                sispeed, self._get_mask(sispeed, filename)
            )

        self.compute_metrics()
        self.results()
        self.save()
        self.plot_results()

    def _get_reference_dataset(self, reference_dataset):
        for filename in self.datasets:
            dataset = self.datasets.get_info(n.DATASET, filename)
            if dataset == reference_dataset:
                return filename
        raise ValueError(
            'Reference dataset "{}" not found'.format(reference_dataset)
        )

    def _get_mask(self, data,filename):
        lat_threshold = self.cfg['latitude_treshold']
        mask = data.coord('latitude').points > lat_threshold
        mask = mask.astype(np.int8)
        dataset_info = self.datasets.get_dataset_info(filename)
        area_cello = iris.load_cube(dataset_info[n.FX_FILES]['areacello'])
        return area_cello.data * mask

    def _load_cube(self, filepath, standard_name):
        cube = iris.load_cube(filepath, standard_name)
        cube.remove_coord('day_of_month')
        cube.remove_coord('day_of_year')
        cube.remove_coord('year')
        return cube

    def compute_metrics(self):
        for dataset in self.siconc.keys():
            logger.info('Compute diagnostics for %s', dataset)
            logger.info('Metrics drift-concentration')
            logger.info('Slope ratio (no unit)')
            slope, intercept, sd, sig = self._get_slope_ratio(
                self.siconc[dataset], self.sispeed[dataset])
            self.slope_drift_siconc[dataset] = slope
            self.intercept_drift_siconc[dataset] = intercept

            logger.info('Metrics drift-thickness')
            logger.info('Slope ratio (no unit)')
            slope, intercept, sd, sig = self._get_slope_ratio(
                self.sivol[dataset], self.sispeed[dataset]
            )
            self.slope_drift_sivol[dataset] = slope
            self.intercept_drift_sivol[dataset] = intercept

        for dataset in self.siconc.keys():
            if dataset == 'reference':
                continue
            logger.info('Compute metrics for %s', dataset)
            logger.info('Compute mean errors (%)')
            self.error_drift_siconc[dataset] = self._compute_error(
                self.siconc[dataset], self.siconc['reference'],
                self.sispeed[dataset], self.siconc['reference']
            )
            self.error_drift_sivol[dataset] = self._compute_error(
                self.sivol[dataset], self.sivol['reference'],
                self.sispeed[dataset], self.siconc['reference']
            )

            logger.info('Compute relative slope ratios ')
            self.slope_ratio_drift_siconc[dataset] = \
                self.slope_drift_siconc[dataset] / self.slope_drift_siconc['reference']
            self.slope_ratio_drift_sivol[dataset] = \
                self.slope_drift_sivol[dataset] / self.slope_drift_sivol['reference']


    def _get_alias(self, filename, reference_dataset):
        filename = self._get_alias_name(filename)
        reference_dataset = self._get_alias_name(reference_dataset)
        if filename == reference_dataset:
            return 'reference'
        else:
            return filename

    def _get_alias_name(self, filename):
        info = self.datasets.get_dataset_info(filename)
        template = '{project}_{dataset}_{experiment}_{ensemble}_{start}_{end}'
        return template.format(
            project=info[n.PROJECT],
            dataset=info[n.DATASET],
            experiment=info[n.EXP],
            ensemble=info[n.ENSEMBLE],
            start=info[n.START_YEAR],
            end=info[n.END_YEAR],
        )

    def _compute_mean(self, data, mask):
        domain_mean = iris.cube.CubeList()

        for map_slice in data.slices_over('time'):
            domain_mean.append(map_slice.collapsed(['latitude', 'longitude'],
                                                   iris.analysis.MEAN,
                                                   weights=mask))
        domain_mean = domain_mean.merge_cube()
        return domain_mean.aggregated_by('month_number', iris.analysis.MEAN)

# def main():
#     cfg = get_cfg()
#     logger.setLevel(cfg['log_level'].upper())
#
#     input_files = get_input_files(cfg)
#
#     for variable_name, filenames in input_files.items():
#         logger.info("Processing variable %s", variable_name)
#         for filename, attributes in filenames.items():
#             plot_filename = os.path.join(
#                 cfg['plot_dir'],
#                 os.path.splitext(os.path.basename(filename))[0] + '.png',
#             )
#             cube = iris.load_cube(filename)
#             cube = cube.collapsed('time', iris.analysis.MEAN)
#
#
# SCICEX = 'scicex'
# LAT = 'lat'
# RAW = 'raw'
#
# DOMAINS = (SCICEX, LAT)
#
#
#
# class SeaIceDrift(object):
#
#     def compute(self):
#         self.load_models()
#         self.load_observations()
#         self.compute_monthly_mean_over_domain()
#         self.compute_multi_year_monthly_mean()
#         self.compute_metrics()
#         self.results()
#         self.save()
#         self.plot_results()
#
#
#     def load_observations(self):
#         logger.info('Load observations')
#         self._load_piomas()
#         self._load_IABP_bouys()
#         self._load_osisaf_sic()
#
#     def _load_IABP_bouys(self):
#         logger.info('Sea ice drift (IABP buoys)')
#         if not self.recalculate and os.path.isfile(self.drift_IABP_file):
#             self.drift_IABP[SCICEX] = iris.load_cube(self.drift_IABP_file,
#                                                      'Drift mean over '
#                                                      'SCICEX domain')
#             self.drift_IABP[LAT] = iris.load_cube(self.drift_IABP_file,
#                                                   'Drift mean over latitude '
#                                                   'threshold')
#             return
#
#         # List sea ice drift speed data files -
#         # retrieved from NSIDC (https://nsidc.org/data/NSIDC-0116)
#         filelist = glob.glob(os.path.join(self.dir_iabp, '*v3.txt'))
#         filelist.sort()
#         n_files = len(filelist)  # number of files
#
#         # Check for empty files
#         def checkIfEmpty(fname, header_cutoff):
#             return os.path.getsize(fname) < header_cutoff
#
        # # Read EASE-Grid file - retrieved from NSIDC
        # # (https://nsidc.org/data/NSIDC-0116)
        # grid_file_path = os.path.join(self.dir_iabp, 'north_x_y_lat_lon.txt')
        # x_grid, y_grid, lat_grid, lon_grid = np.loadtxt(grid_file_path,
        #                                                 usecols=(0, 1, 2, 3),
        #                                                 unpack=True)

        # # Longitude and latitude of SCICEX vertices
        # lon_scicex = np.array((-15., -60., -130., -141, -141, -155, 175,
        #                        172, 163, 126, 110, 80, 57, 33, 8),
        #                       dtype=float)
        # lat_scicex = np.array((87, 86.58, 80, 80, 70, 72, 75.5, 78.5,
        #                        80.5, 78.5, 84.33, 84.42, 85.17, 83.8,
        #                        84.08),
        #                       dtype=float)

#         # Find closest EASE-Grid x and y coordinates of SCICEX vertices
#         x_scicex = np.zeros(15)
#         y_scicex = np.zeros(15)
#         # parameter to take a range around latitude and longitude
#         variation = 0.5
#         for j_scicex in np.arange(len(lon_scicex)):
#             lon0 = lon_scicex[j_scicex] + 2 * variation
#             lat0 = lat_scicex[j_scicex] + 2 * variation
#             for j in np.arange(len(x_grid)):
#                 if (lon_scicex[j_scicex] - variation) <= lon_grid[j] <= \
#                     (lon_scicex[j_scicex] + variation) and \
#                     (lat_scicex[j_scicex] - variation) <= lat_grid[j] <= \
#                     (lat_scicex[j_scicex] + variation):
#                     lon1 = lon_grid[j]
#                     lat1 = lat_grid[j]
#                     if np.abs(lon1 - lon_scicex[j_scicex]) < \
#                         np.abs(lon0 - lon_scicex[j_scicex]) and \
#                         np.abs(lat1 - lat_scicex[j_scicex]) < \
#                         np.abs(lat0 - lat_scicex[j_scicex]):
#                         x_scicex[j_scicex] = x_grid[j]
#                         y_scicex[j_scicex] = y_grid[j]
#                         lon0 = lon_grid[j]
#                         lat0 = lat_grid[j]
#
#         # Put SCICEX vertices (EASE-Grid x and y) into a path
#         scicex_box = matplotlib.path.Path(
#             [(x_scicex[0], y_scicex[0]), (x_scicex[1], y_scicex[1]),
#              (x_scicex[2], y_scicex[2]), (x_scicex[3], y_scicex[3]),
#              (x_scicex[4], y_scicex[4]), (x_scicex[5], y_scicex[5]),
#              (x_scicex[6], y_scicex[6]), (x_scicex[7], y_scicex[7]),
#              (x_scicex[8], y_scicex[8]), (x_scicex[9], y_scicex[9]),
#              (x_scicex[10], y_scicex[10]), (x_scicex[11], y_scicex[11]),
#              (x_scicex[12], y_scicex[12]), (x_scicex[13], y_scicex[13]),
#              (x_scicex[14], y_scicex[14])])
#
#         # Find corners (xmin, xmax, ymin, ymax) of domain north
#         # of latitude threshold
#         indexes = np.where(lat_grid >= self.lat_threshold)
#         xmin = np.min(x_grid[indexes])
#         xmax = np.max(x_grid[indexes])
#
#         ymin = np.min(y_grid[indexes])
#         ymax = np.max(y_grid[indexes])
#
#         # Compute multi-year mean sea ice drift averaged over whole domain
#         # daily mean over SCICEX domain
#         drift_s = np.zeros(n_files)
#         # daily mean over domain north of latitude threshold
#         drift_l = np.zeros(n_files)
#         k = 0
#         # multi-year monthly mean over SCICEX domain
#         drift_scicex = np.zeros(MONTHS_PER_YEAR)
#         # multi-year monthly mean over domain north of latitude threshold
#         drift_lat = np.zeros(MONTHS_PER_YEAR)
#         # number of observations to compute drift_scicex
#         n_scicex = np.zeros(MONTHS_PER_YEAR)
#         # number of observations to compute drift_lat
#         n_lat = np.zeros(MONTHS_PER_YEAR)
#         scicex_index = []
#         lat_index = []
#
#         for i in filelist:
#             if not checkIfEmpty(i, 15):
#                 number = os.path.basename(i)
#                 number = re.sub('icemotion[.]vect[.]buoy[.]([0-9]+)'
#                                 '[.]n[.]v3[.]txt$', r'\1', number)
#                 date = datetime.datetime.strptime(number, '%Y%j')
#                 month_index = date.month - 1
#                 x, y, u, v, t = np.loadtxt(i, skiprows=1, usecols=range(5),
#                                            unpack=True)
#                 if len(np.atleast_1d(x)) > 1:
#                     drift = np.sqrt(u ** 2. + v ** 2.) * 86400. / 100000
#                     n_s = 0  # number of observations to compute drift_s
#                     n_l = 0  # number of observations to compute drift_l
#                     for j in np.arange(len(x)):
#                         # SCICEX box
#                         if scicex_box.contains_points([(x[j] / 5, y[j] /
#                                                                   5)]):
#                             drift_s[k] = drift_s[k] + drift[j]
#                             n_s += + 1
#                             scicex_index.append((k, j))
#                         # Domain north of latitude threshold
#                         if xmin <= round(x[j] / 5) <= xmax and \
#                             ymin <= round(y[j] / 5) <= ymax:
#                             drift_l[k] = drift_l[k] + drift[j]
#                             n_l += 1
#                             lat_index.append((k, j))
#                     if n_s > 0:
#                         drift_s[k] = drift_s[k] / n_s
#                     else:
#                         drift_s[k] = np.nan
#                     if n_l > 0:
#                         drift_l[k] = drift_l[k] / n_l
#                     else:
#                         drift_l[k] = np.nan
#                     if not np.isnan(drift_s[k]):
#                         drift_scicex[month_index] = \
#                             drift_scicex[month_index] + drift_s[k]
#                         n_scicex[month_index] += 1
#                     if not np.isnan(drift_l[k]):
#                         drift_lat[month_index] = drift_lat[month_index] + \
#                                                  drift_l[k]
#                         n_lat[month_index] += 1
#             k += 1
#
#         drift_scicex = np.divide(drift_scicex, n_scicex)
#         drift_lat = np.divide(drift_lat, n_lat)
#
#         drift_scicex = iris.cube.Cube(drift_scicex,
#                                       long_name="Drift mean over SCICEX "
#                                                 "domain",
#                                       var_name="scicex_drift",
#                                       units='km day-1')
#         drift_scicex.add_dim_coord(
#             iris.coords.DimCoord(range(1, 13), var_name='month_number',
#                                  long_name='month_number', units='1'), 0)
#
#         drift_lat = iris.cube.Cube(drift_lat,
#                                    long_name="Drift mean over latitude "
#                                              "threshold",
#                                    var_name="lat_drift",
#                                    units='km day-1')
#         drift_lat.attributes['lat_threshold'] = self.lat_threshold
#         drift_lat.add_dim_coord(
#             iris.coords.DimCoord(range(1, 13), var_name='month_number',
#                                  long_name='month_number', units='1'), 0)
#         iris.save((drift_scicex, drift_lat), self.drift_IABP_file, zlib=True)
#         self.drift_IABP[SCICEX] = drift_scicex
#         self.drift_IABP[LAT] = drift_lat
#
#     def _load_osisaf_sic(self):
#         logger.info('Sea ice concentration (OSI SAF)')
#
#         if self.recalculate or not os.path.isfile(self.siconc_OSISAF_file):
#             fname = 'siconc_OImon_OSISAF_10km*.nc'
#
#             with iris.FUTURE.context(cell_datetime_objects=True):
#                 cubes = iris.load(os.path.join(self.dir_osisaf, fname),
#                                   self.years_constraint &
#                                   iris.Constraint(name=
#                                                   'sea_ice_area_fraction'))
#             equalise_attributes(cubes)
#             iris.util.unify_time_units(cubes)
#             cube = cubes.concatenate_cube()
#             iris.coord_categorisation.add_day_of_year(cube, 'time')
#             iris.coord_categorisation.add_year(cube, 'time')
#             iris.coord_categorisation.add_month_number(cube, 'time')
#             tmp_file = os.path.join(self.dir, 'tmp.nc')
#             iris.save(cube, tmp_file, zlib=True)
#             Cdo().remapbil(self.models[0].siconc_file, input=tmp_file,
#                            output=self.siconc_OSISAF_file)
#             os.remove(tmp_file)
#
#         self.siconc_OSISAF[RAW] = iris.load_cube(self.siconc_OSISAF_file,
#                                                  'sea_ice_area_fraction')
#         return
#
#     def _load_piomas(self):
#         logger.info('Sea ice thickness (PIOMAS)')
#         if not self.recalculate and os.path.isfile(self.sivol_PIOMAS_file):
#             self.sivol_PIOMAS[SCICEX] = iris.load_cube(
#                 self.sivol_PIOMAS_file,
#                 'Sea ice volume mean over SCICEX domain')
#             self.sivol_PIOMAS[LAT] = iris.load_cube(
#                 self.sivol_PIOMAS_file,
#                 'Sea ice volume mean over latitude threshold')
#             return
#         filelist = glob.glob(os.path.join(self.dir_piomas, 'heff.txt*'))
#         # multi-year monthly mean thickness over SCICEX domain
#         sithic = np.zeros(MONTHS_PER_YEAR)
#         # multi-year monthly mean thickness over domain north of threshold
#         h_lat = np.zeros(MONTHS_PER_YEAR)
#         h = [None] * MONTHS_PER_YEAR
#         nyears = 0
#         for i in sorted(filelist):
#             year, lat, lon, h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7], \
#             h[8], h[9], h[10], h[11] = np.loadtxt(i, unpack=True)
#             if self.start_year <= year[0] <= self.end_year:
#                 mask_scicex = self.scicex_domain_1d(lon, lat)
#                 mask_lat = self.lat_domain(lat)
#                 n_s = 0  # number of observations for SCICEX domain
#                 # monthly mean thickness over SCICEX domain
#                 h_s = np.zeros(MONTHS_PER_YEAR)
#                 # monthly mean thickness over domain north of
#                 # latitude threshold
#                 h_l = np.zeros(MONTHS_PER_YEAR)
#                 # number of observations for domain north of
#                 # latitude threshold
#                 n_l = 0
#                 for j in np.arange(len(year)):
#                     if mask_scicex[j]:
#                         for x in range(MONTHS_PER_YEAR):
#                             h_s[x] += h[x][j]
#                         n_s = n_s + 1
#                     if mask_lat[j]:
#                         for x in range(MONTHS_PER_YEAR):
#                             h_l[x] += h[x][j]
#                         n_l = n_l + 1
#                 sithic += (h_s/n_s)
#                 h_lat += (h_l/n_l)
#                 nyears += 1
#         sithic /= nyears
#         h_lat /= nyears
#
#         self.sivol_PIOMAS[SCICEX] = iris.cube.Cube(sithic,
#                                                    long_name="Sea ice "
#                                                              "volume mean
#                                                              "over "
#                                                              "SCICEX domain",
#                                                    var_name="scicex_sivol",
#                                                    units='m')
#         self.sivol_PIOMAS[SCICEX].add_dim_coord(
#             iris.coords.DimCoord(range(1, 13), var_name='month_number',
#                                  long_name='month_number', units='1'), 0)
#
#         self.sivol_PIOMAS[LAT] = iris.cube.Cube(h_lat,
#                                                 long_name="Sea ice volume "
#                                                           "mean over "
#                                                           "latitude "
#                                                           "threshold",
#                                                 var_name="lat_sivol",
#                                                 units='m')
#         self.sivol_PIOMAS[LAT].attributes['lat_threshold'] = \
#             self.lat_threshold
#         self.sivol_PIOMAS[LAT].add_dim_coord(
#             iris.coords.DimCoord(range(1, 13), var_name='month_number',
#                                  long_name='month_number', units='1'), 0)
#         iris.save((self.sivol_PIOMAS[SCICEX], self.sivol_PIOMAS[LAT]),
#                   self.sivol_PIOMAS_file, zlib=True)
#
#     def lat_domain(self, lat):
#         return lat >= self.lat_threshold
#
#     def scicex_domain(self, lon, lat):
#         self._prepare_scicex_box()
#
#         x, y = self.map(lon, lat)
#         # Check which model grid points fall into domain and create mask
#         # (1: inside domain; 0: outside domain)
#         mask_scicex = np.zeros(shape=lon.shape)
#         for jx in np.arange(lon.shape[0]):
#             for jy in np.arange(lon.shape[1]):
#                 if self.scicex_box.contains_points([(x[jx, jy], y[jx, jy])]):
#                     mask_scicex[jx, jy] = 1
#         return mask_scicex
#
#     def scicex_domain_1d(self, lon, lat):
#         # lon: longitude
#         # lat: latitude
#         self._prepare_scicex_box()
#         n = len(lon)
#         x, y = self.map(lon, lat)
#
#         # Check which model grid points fall into domain and create mask
#         # (1: inside domain; 0: outside domain)
#         mask_scicex = np.zeros(n)
#         for j in range(n):
#             if self.scicex_box.contains_points([(x[j], y[j])]):
#                 mask_scicex[j] = 1
#
#         # Return value of mask
#         return mask_scicex
#
#     def _prepare_scicex_box(self):
#         if self.scicex_box is not None:
#             return
#         # Store SCICEX vertices (lon,lat)
#         scicex_vertices = (
#             (-15, 87), (-60, 86.58), (-130, 80), (-141, 80), (-141, 70),
#             (-155, 72), (175, 75.5), (172, 78.5),(163, 80.5),(126, 78.5),
#             (110, 84.33), (80, 84.42), (57, 85.17), (33, 83.8), (8, 84.08))
#         # Map projection
#         boundlat = 50.
#         l0 = 0.
#         self.map = Basemap(projection='nplaea', boundinglat=boundlat,
#                            lon_0=l0, resolution='c')
#         path_list = [self.map(i[0], i[1]) for i in scicex_vertices]
#         self.scicex_box = matplotlib.path.Path(path_list)
#
#     def compute_monthly_mean_over_domain(self):
#         logger.info('Compute monthly mean averaged over whole domain')
#         logger.info('Models')
#         for model in self.models:
#             self._average_model_variables(model)
#         self._average_osisaf()
#
#     def _average_osisaf(self):
#         logger.info('OSISAF averages')
#         siconc_raw = self.siconc_OSISAF[RAW]
#         lon = siconc_raw.coord('longitude').points
#         lat = siconc_raw.coord('latitude').points
#         grid_area = self._get_grid_area(self.models[0])
#         mask_scicex = (siconc_raw.data >=
#                        self.siconc_threshold[SCICEX]) * \
#                       self.scicex_domain(lon, lat) * grid_area
#         mask_lat = (siconc_raw.data >=
#                     self.siconc_threshold[LAT]) * \
#                    self.lat_domain(lat) * grid_area
#         self.siconc_OSISAF[SCICEX] = \
#             siconc_raw.collapsed(['latitude','longitude'],
#                                  iris.analysis.MEAN,
#                                  weights=mask_scicex)
#         self.siconc_OSISAF[LAT] = \
#             siconc_raw.collapsed(['latitude','longitude'],
#                                  iris.analysis.MEAN,
#                                  weights=mask_lat)
#         del siconc_raw
#
#     def _average_model_variables(self, model):
#         logger.info(model.name)
#         lon = model.siconc[RAW].coord('longitude').points
#         lat = model.siconc[RAW].coord('latitude').points
#
#         for dom in DOMAINS:
#             if dom == SCICEX:
#                 mask = self.scicex_domain(lon, lat)
#             else:
#                 mask = self.lat_domain(lat)
#             grid_area = self._get_grid_area(model)
#             mask = mask * grid_area * (model.siconc[RAW].data >=
#                                        self.siconc_threshold[dom])
#
#             model.drift[dom] = self._compute_domain_mean(model.drift[RAW],
#                                                          mask)
#             model.siconc[dom] = self._compute_domain_mean(model.siconc[RAW],
#                                                           mask)
#             model.sivol[dom] = self._compute_domain_mean(model.sivol[RAW],
#                                                          mask)
#         del model.siconc[RAW]
#         del model.drift[RAW]
#         del model.sivol[RAW]
#
#     def _get_grid_area(self, model):
#         # Compute area of grid cells from grid file
#         mesh_file = os.path.join(model.path, model.mesh_file)
#         e1t = iris.load_cube(mesh_file, 'e1t')
#         e2t = iris.load_cube(mesh_file, 'e1t')
#         grid_area = e1t * e2t
#         grid_area.add_dim_coord(iris.coords.DimCoord([0], 'time'), 0)
#         grid_area = grid_area.collapsed('time', iris.analysis.MEAN)
#         return grid_area.data
#
#     def _compute_domain_mean(self, data, mask):
#         return data.collapsed(['latitude', 'longitude'], iris.analysis.MEAN,
#                               weights=mask)
#
#     def compute_multi_year_monthly_mean(self):
#         logger.info('Compute multi-year monthly mean averaged over '
#                          'whole domain')
#
#         for model in self.models:
#             for domain in DOMAINS:
#                 model.drift[domain] = self._multiyear_mean(
#                     model.drift[domain])
#                 model.siconc[domain] = self._multiyear_mean(
#                     model.siconc[domain])
#                 model.sivol[domain] = self._multiyear_mean(
#                     model.sivol[domain])
#
#         for domain in DOMAINS:
#             self.drift_IABP[domain] = self._multiyear_mean(
#                 self.drift_IABP[domain])
#             self.siconc_OSISAF[domain] = self._multiyear_mean(
#                 self.siconc_OSISAF[domain])
#             self.sivol_PIOMAS[domain] = self._multiyear_mean(
#                 self.sivol_PIOMAS[domain])
#
#     def _multiyear_mean(self, data):
#         try:
#             data.coord('month_number')
#         except iris.exceptions.CoordinateNotFoundError:
#             iris.coord_categorisation.add_month_number(data, 'time')
#         return data.aggregated_by('month_number', iris.analysis.MEAN)
#
#     def compute_metrics(self):
#         for model in self.models:
#             logger.info('Compute metrics for {0}'.format(model.name))
#             for domain in DOMAINS:
#                 logger.info('Domain {0}'.format(domain))
#                 logger.info('Metrics drift-concentration')
#                 logger.info('Slope ratio (no unit)')
#                 slope, intercept, sd, sig = self._get_slope_ratio(
#                     model.siconc[domain], model.drift[domain])
#                 slope_obs, intercept_obs, sd_obs, sig_obs = \
#                     self._get_slope_ratio(self.siconc_OSISAF[domain],
#                                           self.drift_IABP[domain])
#                 model.slope_drift_siconc[domain] = slope
#                 model.intercept_drift_siconc[domain] = intercept
#                 self.slope_drift_siconc[domain] = slope_obs
#                 self.intercept_drift_siconc[domain] = intercept_obs
#                 model.slope_ratio_drift_siconc[domain] = slope / slope_obs
#
#                 logger.info('Mean error (%)')
#                 model.error_drift_siconc[domain] = \
#                     self._compute_error(model.siconc[domain],
#                                         self.siconc_OSISAF[domain],
#                                         model.drift[domain],
#                                         self.drift_IABP[domain])
#
#                 logger.info('Metrics drift-thickness')
#                 logger.info('Slope ratio (no unit)')
#
#                 slope, intercept, sd, sig = \
#                     self._get_slope_ratio(model.sivol[domain],
#                                           model.drift[domain])
#                 slope_obs, intercept_obs, sd_obs, sig_obs = \
#                     self._get_slope_ratio(self.sivol_PIOMAS[domain],
#                                           self.drift_IABP[domain])
#
#                 model.slope_drift_sivol[domain] = slope
#                 model.intercept_drift_sivol[domain] = intercept
#                 self.slope_drift_sivol[domain] = slope_obs
#                 self.intercept_drift_sivol[domain] = intercept_obs
#                 model.slope_ratio_drift_sivol[domain] = slope / slope_obs
#
#                 logger.info('Mean error (%)')
#                 model.error_drift_sivol[domain] = \
#                     self._compute_error(model.sivol[domain],
#                                         self.sivol_PIOMAS[domain],
#                                         model.drift[domain],
#                                         self.drift_IABP[domain])
#
    def _compute_error(self, var, var_obs, drift, drift_obs):
        var = var.data
        var_obs = var_obs.data
        drift = drift.data
        drift_obs = drift_obs.data
        return 100. * np.nanmean((np.absolute(var - var_obs) /
                                  np.nanmean(var_obs)) ** 2. +
                                 (np.absolute(drift - drift_obs) /
                                  np.nanmean(drift_obs)) ** 2.)

    def _get_slope_ratio(self, siconc, drift):
        slope, intercept = np.polyfit(siconc.data, drift.data, 1)
        sd, sig = self._sd_slope(slope, intercept, siconc.data, drift.data)
        return slope, intercept, sd, sig

    def _sd_slope(self, slope, intercept, sivar, drift):
        # Parameters
        alpha = 0.05  # significance level
        nfreedom = MONTHS_PER_YEAR - 2  # number of degrees of freedom
        t_crit = stats.t.ppf(1 - alpha / 2, nfreedom)  # critical Student's t

        # Compute standard deviation of slope
        lreg = slope * sivar + intercept  # linear regression
        s_yx = np.sum((drift - lreg) ** 2) / (MONTHS_PER_YEAR - 2)
        SS_xx = np.sum((sivar - np.mean(sivar)) ** 2)
        sd_slope = np.sqrt(s_yx / SS_xx)  # Standard deviation of slope

        # Significance
        ta = slope / sd_slope  # Student's t
        sig_slope = 0
        if np.abs(ta) > t_crit:
            sig_slope = 1

        # Return value of mask
        return sd_slope, sig_slope
#
    def results(self):
        logger.info('Results')
        for model in self.siconc.keys():

            self._print_results(model)

    def _print_results(self, model):
        if model == 'reference':
            return
        logger.info('Dataset {0}'.format(model))
        logger.info(
            'Metrics computed over domain north of %s',
            self.cfg['latitude_treshold']
        )
        logger.info('Slope ratio Drift-Concentration = {0:.3}'
                    ''.format(self.slope_ratio_drift_siconc[model]))
        logger.info('Mean error Drift-Concentration (%) = {0:.4}'
                    ''.format(self.error_drift_siconc[model]))
        logger.info('Slope ratio Drift-Thickness = {0:.3}'
                    ''.format(self.slope_ratio_drift_sivol.get(model, math.nan)))
        logger.info('Mean error Drift-Thickness (%) = {0:.4}'
                    ''.format(self.error_drift_sivol.get(model, math.nan)))

    def save(self):
        if not True:
            return

        logger.info('Save variables')
        for dataset in self.siconc.keys():
            self.save_climatologies(dataset)
            self.save_slope(dataset)

    def save_climatologies(self, dataset):
        base_path = os.path.join(self.cfg[n.WORK_DIR], dataset)
        if not os.path.isdir(base_path):
            os.makedirs(base_path)
        iris.save(
            self.sispeed[dataset],
            os.path.join(
                base_path,
                'sispeed_clim.nc'
            ),
            zlib=True
        )
        iris.save(
            self.siconc[dataset],
            os.path.join(
                base_path,
                'siconc_clim.nc'
            ),
            zlib=True
        )
        iris.save(
            self.sivol[dataset],
            os.path.join(
                base_path,
                'sivol_clim.nc'
            ),
            zlib=True
        )
    def save_slope(self, dataset):
        base_path = os.path.join(self.cfg[n.WORK_DIR], dataset)
        if not os.path.isdir(base_path):
            os.makedirs(base_path)
        siconc_path = os.path.join(base_path, 'metric_drift_siconc.csv')
        sivol_path = os.path.join(base_path, 'metric_drift_sivol.csv')
        with open(siconc_path, 'w') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(('slope', 'intercept',
                                 'slope_ratio', 'error'))
            csv_writer.writerow((
                self.slope_drift_siconc[dataset],
                self.intercept_drift_siconc[dataset],
                self.slope_ratio_drift_siconc.get(dataset, None),
                self.error_drift_siconc.get(dataset, None)
            ))

        with open(sivol_path, 'w') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(('slope', 'intercept',
                                 'slope_ratio', 'error'))
            csv_writer.writerow((
                self.slope_drift_sivol[dataset],
                self.intercept_drift_sivol[dataset],
                self.slope_ratio_drift_sivol.get(dataset, None),
                self.error_drift_sivol.get(dataset, None)
            ))

    def plot_results(self):
        logger.info('Plot results')
        for model in self.siconc.keys():
            if model == 'reference':
                continue
            logger.info('Results for {0}'.format(model))
            self._plot_domain(model)

    def _plot_domain(self, dataset):
        fig, ax = plt.subplots(1, 2, figsize=(18, 6))
        self._plot_drift_siconc(ax[0], dataset)
        self._plot_drift_sivol(ax[1], dataset)

        # Save figure
        if True:
            base_path = os.path.join(self.cfg[n.PLOT_DIR], dataset)
            if not os.path.isdir(base_path):
                os.makedirs(base_path)
            fig.savefig(os.path.join(
                base_path,
                'Drift-Strength.{0}'.format(self.cfg[n.OUTPUT_FILE_TYPE])
                )
            )

    def _plot_drift_sivol(self, ax, dataset):

        drift = self.sispeed[dataset].data
        sivol = self.sivol[dataset].data

        slope_sivol = self.slope_drift_sivol[dataset]
        slope_ratio_sivol = self.slope_ratio_drift_sivol[dataset]
        intercept_sivol = self.intercept_drift_sivol[dataset]
        error_sivol = self.error_drift_sivol[dataset]

        slope_sivol_obs = self.slope_drift_sivol[dataset]
        intercept_sivol_obs = self.intercept_drift_sivol[dataset]

        drift_obs = self.sispeed['reference'].data
        sivol_obs = self.sivol['reference'].data

        ax.plot([sivol[-1], sivol[0]], [drift[-1], drift[0]], 'r-',
                linewidth=2)
        ax.plot(sivol, drift, 'ro-', label='NEMO-LIM3.6', linewidth=2)
        ax.plot(sivol, slope_sivol * sivol + intercept_sivol, 'r:',
                linewidth=2)

        ax.plot([sivol_obs[-1], sivol_obs[0]], [drift_obs[-1], drift_obs[0]],
                'b-', linewidth=2)
        ax.plot(sivol_obs, drift_obs, 'bo-',
                label='reference volume / speed ($s_h$=' +
                      str(np.round(slope_ratio_sivol, 1)) +
                      '; $\epsilon_h$=' + str(np.round(error_sivol, 1)) +
                      '$\%$)', linewidth=2)
        ax.plot(sivol_obs, slope_sivol_obs * sivol_obs + intercept_sivol_obs,
                'b:', linewidth=2)
        ax.legend(loc='lower right', shadow=True, frameon=False, fontsize=12)
        ax.set_xlabel('Sea ice thickness (m)', fontsize=18)
        ax.set_ylabel('Sea ice drift speed (km d$^{-1}$)', fontsize=18)
        ax.tick_params(axis='both', labelsize=14)
        high_sivol, low_sivol = self._get_plot_limits(sivol, sivol_obs)
        high_drift, low_drift = self._get_plot_limits(drift, drift_obs)
        ax.axis([low_sivol, high_sivol, low_drift, high_drift])
        self._annotate_points(ax, sivol, drift)
        self._annotate_points(ax, sivol_obs, drift_obs)
        ax.grid()
        ax.set_title('Seasonal cycle 1979-2013 {0}'.format(dataset),
                     fontsize=18)

    def _plot_drift_siconc(self, ax, dataset):
        drift = self.sispeed[dataset].data
        siconc = self.siconc[dataset].data

        slope_siconc = self.slope_drift_siconc[dataset]
        slope_ratio_siconc = self.slope_ratio_drift_sivol[dataset]
        intercept_siconc = self.intercept_drift_siconc[dataset]
        error_siconc = self.error_drift_siconc[dataset]

        slope_siconc_obs = self.slope_drift_siconc[dataset]
        intercept_siconc_obs = self.intercept_drift_siconc[dataset]

        drift_obs = self.sispeed['reference'].data
        siconc_obs = self.siconc['reference'].data

        ax.plot(siconc, drift, 'ro-', label=dataset)
        ax.plot(siconc, slope_siconc * siconc + intercept_siconc, 'r:',
                linewidth=2)
        ax.plot(siconc_obs, drift_obs, 'bo-',
                label='reference ($s_A$=' +
                      str(np.round(slope_ratio_siconc, 1)) +
                      '; $\epsilon_A$=' + str(np.round(error_siconc, 1)) +
                      '$\%$)')
        ax.plot(siconc_obs, slope_siconc_obs * siconc_obs +
                intercept_siconc_obs,
                'b:', linewidth=2)
        ax.legend(loc='lower left', shadow=True, frameon=False, fontsize=12)
        ax.set_xlabel('Sea ice concentration', fontsize=18)
        ax.set_ylabel('Sea ice drift speed (km d$^{-1}$)', fontsize=18)
        ax.tick_params(axis='both', labelsize=14)
        high_drift, low_drift = self._get_plot_limits(drift, drift_obs)
        ax.axis([0.5, 1, low_drift, high_drift])
        self._annotate_points(ax, siconc, drift)
        self._annotate_points(ax, siconc_obs, drift_obs)
        ax.grid(linewidth=0.01)
        ax.set_title('Seasonal cycle 1979-2013 {0}'.format(dataset),
                     fontsize=18)

    def _annotate_points(self, ax, xvalues, yvalues):
        for x, y, z in zip(xvalues, yvalues, range(1, 12 + 1)):
            ax.annotate(calendar.month_abbr[z][0], xy=(x, y), xytext=(10, 5),
                        ha='right', textcoords='offset points')

    def _get_plot_limits(self, sivol, sivol_obs):
        low = min(min(sivol), min(sivol_obs)) - 0.4
        low = 0.5 * math.floor(2.0 * low)
        high = max(max(sivol), max(sivol_obs)) + 0.4
        high = 0.5 * math.ceil(2.0 * high)
        return high, low
#
#
# class Model(object):
#     SICONC = 'siconc'
#     SIVOL = 'sivol'
#     DRIFT = 'drift'
#
#     def __init__(self, name, path_to_data, mesh_file):
#         self._slope_drift_sivol = None
#         self._slope_drift_siconc = None
#         self.name = name
#         self.path = path_to_data
#         self.mesh_file = mesh_file
#         self._scratch_dir = None
#         self.start_year = None
#         self.end_year = None
#
#         self._data = {Model.SICONC: {}, Model.SIVOL: {}, Model.DRIFT: {}}
#
#         self._slope_drift_sivol = {}
#         self._slope_drift_siconc = {}
#         self._intercept_drift_siconc = {}
#         self._intercept_drift_sivol = {}
#
#         self._slope_ratio_drift_sivol = {}
#         self._slope_ratio_drift_siconc = {}
#         self._error_drift_sivol = {}
#         self._error_drift_siconc = {}
#     @property
#     def scratch_dir(self):
#         return self._scratch_dir
#

#
#     @property
#     def scratch_dir(self):
#         return self._scratch_dir
#
#     @scratch_dir.setter
#     def scratch_dir(self, value):
#         self._scratch_dir = value
#         if self._scratch_dir and \
#             not os.path.isdir(os.path.join(self.scratch_dir, self.name)):
#             os.mkdir(os.path.join(self.scratch_dir, self.name))
#
#     @property
#     def siconc(self):
#         return self._data[Model.SICONC]
#
#     @property
#     def sivol(self):
#         return self._data[Model.SIVOL]
#
#     @property
#     def drift(self):
#         return self._data[Model.DRIFT]
#
#     @property
#     def slope_drift_siconc(self):
#         return self._slope_drift_siconc
#
#     @property
#     def slope_drift_sivol(self):
#         return self._slope_drift_sivol
#
#     @property
#     def intercept_drift_siconc(self):
#         return self._intercept_drift_siconc
#
#     @property
#     def intercept_drift_sivol(self):
#         return self._intercept_drift_sivol
#
#     @property
#     def slope_ratio_drift_siconc(self):
#         return self._slope_ratio_drift_siconc
#
#     @property
#     def slope_ratio_drift_sivol(self):
#         return self._slope_ratio_drift_sivol
#
#     @property
#     def error_drift_siconc(self):
#         return self._error_drift_siconc
#
#     @property
#     def error_drift_sivol(self):
#         return self._error_drift_sivol
#
#     @property
#     def drift_file(self):
#         return os.path.join(self.scratch_dir, self.name,
#                             'drift_{0.start_year}-{0.end_year}'
#                             '.nc'.format(self))
#
#     @property
#     def siconc_file(self):
#         return os.path.join(self.scratch_dir, self.name,
#                             'siconc_{0.start_year}-{0.end_year}'
#                             '.nc'.format(self))
#
#     @property
#     def sivol_file(self):
#         return os.path.join(self.scratch_dir, self.name,
#                             'sivol_{0.start_year}-{0.end_year}'
#                             '.nc'.format(self))
#
#
#     def drift_clim_file(self, domain):
#         return os.path.join(self.scratch_dir, self.name,
#                             'drift_clim_{1}_{0.start_year}-{0.end_year}'
#                             '.nc'.format(self, domain))
#
#     def siconc_clim_file(self, domain):
#         return os.path.join(self.scratch_dir, self.name,
#                             'siconc_clim_{1}_{0.start_year}-{0.end_year}'
#                             '.nc'.format(self, domain))
#
#     def sivol_clim_file(self, domain):
#         return os.path.join(self.scratch_dir, self.name,
#                             'sivol_clim_{1}_{0.start_year}-{0.end_year}'
#                             '.nc'.format(self, domain))
#
#
# if __name__ == '__main__':
#     nemo = Model('NEMO',
#                  '/group_workspaces/jasmin2/primavera1/WP2/OCE/'
#                  'NEMO-LIM3-6/HIST_ORCA1/',
#                  'ORCA1_mesh_mask.nc')
#     start_time = datetime.datetime.now()
#     SeaIceDrift(1979, 2013, (nemo,)).compute()
#     print ('Ellapsed time: {0}'.format(datetime.datetime.now() - start_time))


def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        SeaIceDrift(config).compute()


if __name__ == '__main__':
    main()
