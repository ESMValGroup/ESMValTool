"""Sea ice drift diagnostic."""
import os
import logging
import math
import csv
import warnings

import numpy as np
from scipy import stats
from matplotlib import pyplot as plt

import iris
import iris.cube
import iris.analysis
import iris.analysis.cartography
import iris.coords
from iris.util import broadcast_to_shape
from iris.aux_factory import AuxCoordFactory
from pyproj import Transformer
from shapely.geometry import Polygon, Point


import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared._base import ProvenanceLogger

logger = logging.getLogger(os.path.basename(__file__))

MONTHS_PER_YEAR = 12

warnings.filterwarnings("ignore", category=DeprecationWarning)


class SeaIceDrift():
    """Class to compute SeaIce Drift metric."""

    def __init__(self, cfg):
        self.cfg = cfg
        self.datasets = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.variables = esmvaltool.diag_scripts.shared.Variables(self.cfg)

        self.references = {}
        self.siconc = {}
        self.sivol = {}
        self.sispeed = {}
        self.region_mask = {}

        self.slope_drift_sic = {}
        self.intercept_drift_siconc = {}
        self.slope_ratio_drift_siconc = {}
        self.error_drift_siconc = {}

        self.slope_drift_sivol = {}
        self.intercept_drift_sivol = {}
        self.slope_ratio_drift_sivol = {}
        self.error_drift_sivol = {}

    def compute(self):
        """Compute metric"""
        logger.info('Loading sea ice concentration')
        siconc_original = {}
        siconc_files = self.datasets.get_path_list(
            standard_name='sea_ice_area_fraction')
        for filename in siconc_files:
            reference_dataset = self._get_reference_dataset(
                'sic',
                self.datasets.get_info('reference_dataset', filename)
            )
            alias = self._get_alias(filename, reference_dataset)
            siconc = iris.load_cube(filename, 'sea_ice_area_fraction')
            siconc.convert_units('1.0')
            siconc_original[alias] = siconc

            self.siconc[alias] = self._compute_mean(
                siconc, self._get_mask(siconc, filename)
            )

        logger.info('Loading sea ice thickness')
        sithick_files = self.datasets.get_path_list(
            standard_name='sea_ice_thickness')
        for filename in sithick_files:
            reference_dataset = self._get_reference_dataset(
                'sithick',
                self.datasets.get_info('reference_dataset', filename)
            )
            alias = self._get_alias(filename, reference_dataset)
            sithick = iris.load_cube(filename, 'sea_ice_thickness')
            self.sivol[alias] = self._compute_mean(
                sithick,
                self._get_mask(sithick, filename)
            )

        logger.info('Load sea ice velocities')
        sispeed_files = self.datasets.get_path_list(
            standard_name='sea_ice_speed')
        obs_file = self.cfg.get('sispeed_obs', '')
        for filename in sispeed_files:
            reference_dataset = self._get_reference_dataset(
                'sispeed',
                self.datasets.get_info('reference_dataset', filename)
            )
            alias = self._get_alias(filename, reference_dataset)
            if obs_file and alias == 'reference':
                obs_data = np.load(obs_file)
                obs_data = obs_data.reshape((12, 35), order='F')
                logger.debug(obs_data)
                sispeed = iris.cube.Cube(
                    obs_data,
                    'sea_ice_speed',
                    units='km day-1'
                )
                sispeed.add_dim_coord(
                    iris.coords.DimCoord(
                        range(1, 13), var_name='month_number'
                    ),
                    0
                )
                sispeed.add_dim_coord(
                    iris.coords.DimCoord(
                        range(1979, 1979 + 35), var_name='year'
                    ),
                    1
                )
                sispeed.extract(
                    iris.Constraint(year=lambda c: 1979 <= c <= 2005)
                )
                sispeed = sispeed.collapsed('year', iris.analysis.MEAN)
                logger.debug(sispeed)
                self.sispeed[alias] = sispeed
            else:
                sispeed = iris.load_cube(filename, 'sea_ice_speed')
                sispeed.convert_units('km day-1')
                self.sispeed[alias] = self._compute_mean(
                    sispeed, self._get_mask(sispeed, filename)
                )

        self._compute_metrics()
        self._results()
        self._save()
        self._plot_results()

    def _get_reference_dataset(self, var, reference_dataset):
        for filename in self.datasets:
            dataset = self.datasets.get_info(n.DATASET, filename)
            if dataset == reference_dataset:
                self.references[var] = self.datasets.get_info(
                    n.ALIAS, filename)
                return filename
        raise ValueError(f'Reference dataset {reference_dataset} not found')

    def _get_mask(self, data, filename):
        if 'latitude_treshold' in self.cfg:
            lat_threshold = self.cfg['latitude_treshold']
            mask = data.coord('latitude').points > lat_threshold
            mask = mask.astype(np.int8)
        else:
            polygon = self.cfg['polygon']
            factory = InsidePolygonFactory(
                polygon,
                data.coord('latitude'),
                data.coord('longitude'),
            )
            data.add_aux_factory(factory)
            mask = data.coord('Inside polygon').points
            mask = mask.astype(np.int8)
            coord = data.coord('Inside polygon')
            dim_coords = data.coord_dims(coord)
            data.remove_aux_factory(factory)
            data.add_aux_coord(coord, dim_coords)
            data.remove_coord('Inside polygon')

        dataset_info = self.datasets.get_dataset_info(filename)
        var_info = esmvaltool.diag_scripts.shared.group_metadata(
            self.cfg['input_data'].values(), 'alias'
        )[dataset_info[n.ALIAS]]
        var_info = esmvaltool.diag_scripts.shared.group_metadata(
            var_info, 'short_name')
        if 'areacello' in var_info:
            area_file = var_info['areacello'][0]['filename']
            area_cello = iris.load_cube(area_file)
        else:
            area_cello = iris.analysis.cartography.area_weights(data)

        return area_cello.data * mask

    def _compute_metrics(self):
        for dataset in self.siconc:
            logger.info('Compute diagnostics for %s', dataset)
            logger.info('Metrics drift-concentration')
            logger.debug('Siconc: %s', self.siconc[dataset].data)
            logger.debug('Sispeed: %s', self.sispeed[dataset].data)
            logger.info('Slope ratio (no unit)')
            slope, intercept, _, _ = self._get_slope_ratio(
                self.siconc[dataset], self.sispeed[dataset])
            self.slope_drift_sic[dataset] = slope
            self.intercept_drift_siconc[dataset] = intercept

            logger.info('Metrics drift-thickness')
            logger.debug('sivol: %s', self.sivol[dataset].data)
            logger.debug('Sispeed: %s', self.sispeed[dataset].data)
            logger.info('Slope ratio (no unit)')
            slope, intercept, _, _ = self._get_slope_ratio(
                self.sivol[dataset], self.sispeed[dataset]
            )
            self.slope_drift_sivol[dataset] = slope
            self.intercept_drift_sivol[dataset] = intercept

        for dataset in self.siconc:
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
                self.slope_drift_sic[dataset] / \
                self.slope_drift_sic['reference']
            self.slope_ratio_drift_sivol[dataset] = \
                self.slope_drift_sivol[dataset] / \
                self.slope_drift_sivol['reference']

    def _get_alias(self, filename, reference_dataset):
        filename = self._get_alias_name(filename)
        reference_dataset = self._get_alias_name(reference_dataset)
        if filename == reference_dataset:
            return 'reference'
        return filename

    def _get_alias_name(self, filename):
        info = self.datasets.get_dataset_info(filename)
        return info[n.ALIAS]

    @staticmethod
    def _compute_mean(data, weights):
        mapping = set(
            data.coord_dims('latitude') + data.coord_dims('longitude')
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return data.collapsed(
                ('latitude', 'longitude'),
                iris.analysis.MEAN,
                weights=broadcast_to_shape(weights, data.shape, mapping)
            )

    @staticmethod
    def _compute_error(var, var_obs, drift, drift_obs):
        var = var.data
        var_obs = var_obs.data
        drift = drift.data
        drift_obs = drift_obs.data

        return 100. * np.nanmean(np.sqrt(
            SeaIceDrift._var_error(var, var_obs) +
            SeaIceDrift._var_error(drift, drift_obs)
        ))

    @staticmethod
    def _var_error(var, obs):
        var_error = np.absolute(var - obs)
        var_mean = np.nanmean(obs)
        var_error_normal = var_error / var_mean
        return var_error_normal ** 2

    @staticmethod
    def _get_slope_ratio(siconc, drift):
        slope, intercept = np.polyfit(siconc.data, drift.data, 1)
        std_dev, sig = SeaIceDrift._sd_slope(
            slope, intercept, siconc.data, drift.data)
        return slope, intercept, std_dev, sig

    @staticmethod
    def _sd_slope(slope, intercept, sivar, drift):
        # Parameters
        alpha = 0.05  # significance level
        nfreedom = MONTHS_PER_YEAR - 2  # number of degrees of freedom
        t_crit = stats.t.ppf(1 - alpha / 2, nfreedom)  # critical Student's t

        # Compute standard deviation of slope
        lreg = slope * sivar + intercept  # linear regression
        s_yx = np.sum((drift - lreg) ** 2) / (MONTHS_PER_YEAR - 2)
        ss_xx = np.sum((sivar - np.mean(sivar)) ** 2)
        sd_slope = np.sqrt(s_yx / ss_xx)  # Standard deviation of slope

        # Significance
        t_student = slope / sd_slope
        sig_slope = 0
        if np.abs(t_student) > t_crit:
            sig_slope = 1

        return sd_slope, sig_slope

    def _results(self):
        logger.info('Results')
        for model in self.siconc:
            self._print_results(model)

    def _print_results(self, model):
        if model == 'reference':
            return
        logger.info('Dataset %s', model)
        if 'latitude_treshold' in self.cfg:
            logger.info(
                'Metrics computed over domain north of %s',
                self.cfg['latitude_treshold']
            )
        else:
            logger.info(
                'Metrics computed inside %s region',
                self.cfg.get('polygon_name', 'SCICEX')
            )
        logger.info('Slope ratio Drift-Concentration = {0:.3}'
                    ''.format(self.slope_ratio_drift_siconc[model]))
        logger.info('Mean error Drift-Concentration (%) = {0:.4}'
                    ''.format(self.error_drift_siconc[model]))
        logger.info('Slope ratio Drift-Thickness = {0:.3}'.format(
            self.slope_ratio_drift_sivol.get(model, math.nan)
        ))
        logger.info('Mean error Drift-Thickness (%) = {0:.4}'
                    ''.format(self.error_drift_sivol.get(model, math.nan)))

    def _save(self):
        if not True:
            return
        logger.info('Save variables')
        for dataset in self.siconc:
            self._save_slope(dataset)

    def _save_slope(self, dataset):
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
                self.slope_drift_sic[dataset],
                self.intercept_drift_siconc[dataset],
                self.slope_ratio_drift_siconc.get(dataset, None),
                self.error_drift_siconc.get(dataset, None)
            ))

        sic_data, ancestors_sic = self._get_data_and_ancestors(dataset, 'sic')
        _, ancestors_sispeed = self._get_data_and_ancestors(dataset, 'sispeed')
        sithick_data, ancestors_sithick = self._get_data_and_ancestors(
            dataset, 'sithick')
        caption = (
            f"Drift - siconc metric between {sic_data[n.START_YEAR]} and "
            f"{sic_data[n.END_YEAR]} according to {dataset}"
        )
        self._create_prov_record(
            siconc_path, caption, ancestors_sic + ancestors_sispeed)

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
        caption = (
            f"Drift - sithick metric between {sithick_data[n.START_YEAR]} and "
            f"{sithick_data[n.END_YEAR]} according to {dataset}"
        )
        self._create_prov_record(
            sivol_path, caption, ancestors_sithick + ancestors_sispeed)

    def _plot_results(self):
        if not self.cfg[n.WRITE_PLOTS]:
            return
        logger.info('Plotting results')
        for model in self.siconc:
            if model == 'reference':
                continue
            logger.info('Results for %s', model)
            self._plot_domain(model)

    def _plot_domain(self, dataset):
        fig, axes = plt.subplots(1, 2, figsize=(18, 6))
        plt.suptitle('Seasonal cycle {0}'.format(dataset), fontsize=18)
        self._plot_drift_siconc(axes[0], dataset)
        self._plot_drift_sivol(axes[1], dataset)
        base_path = os.path.join(self.cfg[n.PLOT_DIR], dataset)
        if not os.path.isdir(base_path):
            os.makedirs(base_path)
        plot_path = os.path.join(
            base_path,
            'drift-strength.{0}'.format(self.cfg[n.OUTPUT_FILE_TYPE])
        )
        fig.savefig(plot_path)

        sic_data, ancestors_sic = self._get_data_and_ancestors(dataset, 'sic')
        _, ancestors_sispeed = self._get_data_and_ancestors(dataset, 'sispeed')
        _, ancestors_sithick = self._get_data_and_ancestors(dataset, 'sithick')

        _, ancestors_sic_ref = self._get_data_and_ancestors(
            'reference', 'sic')
        _, ancestors_sispeed_ref = self._get_data_and_ancestors(
            'reference', 'sispeed')
        _, ancestors_sithick_ref = self._get_data_and_ancestors(
            'reference', 'sithick')

        caption = (
            "Drift - sithick and drift - siconc plot between "
            f"{sic_data[n.START_YEAR]} and "
            f"{sic_data[n.END_YEAR]} for {dataset} and reference"
        )
        self._create_prov_record(
            plot_path,
            caption,
            ancestors_sic + ancestors_sispeed + ancestors_sithick +
            ancestors_sic_ref + ancestors_sispeed_ref + ancestors_sithick_ref
        )

    def _plot_drift_sivol(self, axes, dataset):
        drift = self.sispeed[dataset].data
        sivol = self.sivol[dataset].data

        slope_sivol = self.slope_drift_sivol[dataset]
        intercept_sivol = self.intercept_drift_sivol[dataset]

        slope_sivol_obs = self.slope_drift_sivol['reference']
        intercept_sivol_obs = self.intercept_drift_sivol['reference']

        # slope_ratio_sivol = self.slope_ratio_drift_sivol[dataset]
        # error_sivol = self.error_drift_sivol[dataset]

        drift_obs = self.sispeed['reference'].data
        sivol_obs = self.sivol['reference'].data

        axes.plot([sivol[-1], sivol[0]], [drift[-1], drift[0]], 'r-',
                  linewidth=2)
        axes.plot(sivol, drift, 'ro-', label='model', linewidth=2)
        axes.plot(sivol, slope_sivol * sivol + intercept_sivol, 'r:',
                  linewidth=2)

        axes.plot([sivol_obs[-1], sivol_obs[0]], [drift_obs[-1], drift_obs[0]],
                  'b-', linewidth=2)
        axes.plot(
            sivol_obs,
            drift_obs,
            'bo-',
            label=r'reference',
            #   str(np.round(slope_ratio_sivol, 1)) +
            # #   r'; $\epsilon_h$=' +
            # #   str(np.round(error_sivol, 1)) +
            #   r')',
            linewidth=2
        )
        axes.plot(sivol_obs, slope_sivol_obs * sivol_obs + intercept_sivol_obs,
                  'b:', linewidth=2)

        axes.set_xlabel('Sea ice thickness (m)', fontsize=18)
        axes.set_ylabel('Sea ice drift speed (km d$^{-1}$)', fontsize=18)
        axes.tick_params(axis='both', labelsize=14)
        high_sivol, low_sivol = self._get_plot_limits(sivol, sivol_obs, 0.2)
        high_drift, low_drift = self._get_plot_limits(drift, drift_obs)
        axes.axis([low_sivol, high_sivol, low_drift, high_drift])
        axes.legend(loc='lower left', shadow=True, frameon=False, fontsize=12)
        self._annotate_points(axes, sivol, drift)
        self._annotate_points(axes, sivol_obs, drift_obs)
        axes.grid()

    def _plot_drift_siconc(self, axes, dataset):
        drift = self.sispeed[dataset].data
        siconc = self.siconc[dataset].data

        slope_siconc = self.slope_drift_sic[dataset]
        # slope_ratio_siconc = self.slope_ratio_drift_sivol[dataset]
        intercept_siconc = self.intercept_drift_siconc[dataset]
        # error_siconc = self.error_drift_siconc[dataset]

        slope_siconc_obs = self.slope_drift_sic['reference']
        intercept_siconc_obs = self.intercept_drift_siconc['reference']

        drift_obs = self.sispeed['reference'].data
        siconc_obs = self.siconc['reference'].data

        axes.plot(siconc, drift, 'ro', label='model')
        axes.plot(siconc, slope_siconc * siconc + intercept_siconc, 'r:',
                  linewidth=2)

        axes.plot(
            siconc_obs,
            drift_obs,
            'bo',
            label='reference'
        )
        axes.plot(siconc_obs, slope_siconc_obs * siconc_obs +
                  intercept_siconc_obs,
                  'b:', linewidth=2)

        axes.set_xlabel('Sea ice concentration', fontsize=18)
        axes.set_ylabel('Sea ice drift speed (km d$^{-1}$)', fontsize=18)
        axes.tick_params(axis='both', labelsize=14)
        high_drift, low_drift = self._get_plot_limits(drift, drift_obs)
        _, low_siconc = SeaIceDrift._get_plot_limits(siconc, siconc_obs, 0.1)
        axes.axis([low_siconc, 1.01, low_drift, high_drift])
        axes.legend(loc='lower left', shadow=True, frameon=False, fontsize=12)
        SeaIceDrift._annotate_points(axes, siconc, drift)
        SeaIceDrift._annotate_points(axes, siconc_obs, drift_obs)
        axes.grid()

    @staticmethod
    def _annotate_points(axes, xvalues, yvalues):
        for i, j, k in zip(xvalues, yvalues, range(1, 12 + 1)):
            axes.annotate(k, xy=(i, j), xytext=(10, 5),
                          ha='right', textcoords='offset points')

    @staticmethod
    def _get_plot_limits(sivol, sivol_obs, step=0.55):
        low = min(min(sivol), min(sivol_obs)) - 0.5 * step
        low = step * math.floor(low / step)
        low = max(low, 0)
        high = max(max(sivol), max(sivol_obs)) + 0.5 * step
        high = step * math.ceil(high / step)
        return high, low

    def _create_prov_record(self, filepath, caption, ancestors):
        record = {
            'caption': caption,
            'domains': ['nhpolar'],
            'autors': ['docquier_david'],
            'references': ['docquier2017cryo'],
            'ancestors': ancestors
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(filepath, record)

    def _get_data_and_ancestors(self, dataset, var):
        if dataset == 'reference':
            dataset = self.references[var]

        data = esmvaltool.diag_scripts.shared.group_metadata(
            self.cfg['input_data'].values(), 'alias'
        )[dataset]
        data = esmvaltool.diag_scripts.shared.group_metadata(
            data, 'short_name'
        )
        return (
            data[var][0],
            [data[var][0]['filename'], data['areacello'][0]['filename']]
        )


class InsidePolygonFactory(AuxCoordFactory):
    """Defines a coordinate."""

    def __init__(self, polygon=None, lat=None, lon=None):
        """
        Args:
        * polygon: List
            List of (lon, lat) tuples defining the polygon
        * lat: Coord
            The coordinate providing the latitudes.
        * lon: Coord
            The coordinate providing the longitudes.
        """
        super(InsidePolygonFactory, self).__init__()
        self.lat = lat
        self.lon = lon
        self.standard_name = None
        self.long_name = 'Inside polygon'
        self.var_name = 'inpoly'
        self.units = '1.0'
        self.attributes = {}

        polygon.append(polygon[0])
        self.transformer = Transformer.from_crs(
            "WGS84",
            "North_Pole_Stereographic"
        )

        transformed = []
        for lon_val, lat_val in polygon:
            transformed.append(self.transformer.transform(lon_val, lat_val))
        self.polygon = Polygon(transformed)

    @property
    def dependencies(self):
        """
        Return a dict mapping from constructor names to coordinates.
        """
        return {'lat': self.lat, 'lon': self.lon}

    def _derive(self, lat, lon):
        def in_polygon(lat, lon):
            """Check if point is inside polygon"""
            if lon > 180:
                lon -= 360
            point = self.transformer.transform(lon, lat)
            if self.polygon.contains(Point(point[0], point[1])):
                return 1.
            return np.nan
        vectorized = np.vectorize(in_polygon)
        return vectorized(lat, lon)

    def make_coord(self, coord_dims_func):
        """
        Returns a new :class:`iris.coords.AuxCoord`

        Args:
        * coord_dims_func:
            A callable which can return the list of dimensions relevant
            to a given coordinate.
            See :meth:`iris.cube.Cube.coord_dims()`.
        """
        # Which dimensions are relevant?
        derived_dims = self.derived_dims(coord_dims_func)
        dependency_dims = self._dependency_dims(coord_dims_func)

        # Build the points array.
        nd_points_by_key = self._remap(dependency_dims, derived_dims)
        points = self._derive(nd_points_by_key['lat'],
                              nd_points_by_key['lon'],)

        in_polygon = iris.coords.AuxCoord(
            points,
            standard_name=self.standard_name,
            long_name=self.long_name,
            var_name=self.var_name,
            units=self.units,
            bounds=None,
            attributes=self.attributes,
            coord_system=self.coord_system
        )
        return in_polygon

    def update(self, old_coord, new_coord=None):
        """
        Notify factory about the removal/replacement of a coordinate


        Args:
        * old_coord:
            The coordinate to be removed/replaced.
        * new_coord:
            If None, any dependency using old_coord is removed, otherwise
            any dependency using old_coord is updated to use new_coord.
        """
        if self.lat is old_coord:
            self.lat = new_coord
        elif self.lon is old_coord:
            self.lon = new_coord


if __name__ == '__main__':
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        SeaIceDrift(config).compute()
