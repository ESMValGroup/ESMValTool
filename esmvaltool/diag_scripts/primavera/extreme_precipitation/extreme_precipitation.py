"""
Generalised extreme value analysis.

Apply GEV analysis grid point by grid point.
"""

import os
import logging

import numpy as np
import iris
import iris.cube
import iris.analysis
import iris.util

from rpy2.robjects import r, numpy2ri, IntVector
from rpy2.robjects.packages import importr

import esmvaltool.diag_scripts.shared
from esmvaltool.diag_scripts.shared.plot import quickplot
import esmvaltool.diag_scripts.shared.names as n

logger = logging.getLogger(os.path.basename(__file__))


class ExtremeValuesAnalysis():
    """Extreme values analysis."""

    def __init__(self, config):
        self.cfg = config
        self.filenames = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.return_period = self.cfg['return_period']
        self.confidence_interval = self.cfg['confidence_interval']

        self.gev_param_symbols = ['mu', 'sigma', 'xi']
        self.gev_param_name = dict(mu='location', sigma='scale', xi='shape')
        self.return_period = IntVector(self.return_period)
        self.extremes = importr('extRemes')

    @staticmethod
    def _get_return_period_name(period):
        return '{}-year level'.format(period)

    def compute(self):
        """Compute extreme values metric."""
        numpy2ri.activate()
        r['options'](warn=-1)
        for filename in self.filenames:
            logger.info('Processing %s', filename)
            cube = iris.load_cube(filename)
            logger.debug(cube)
            model_cube = cube[0, ...]

            for season in set(cube.coord('clim_season').points):
                # Clean R objects
                r('rm(list = ls())')
                logger.info('Processing season %s', season)
                season_cube = cube.extract(iris.Constraint(clim_season=season))
                fevd, ret_values = self._compute_season(
                    season_cube, model_cube)
                self._plot_results(filename, season, fevd, ret_values)
                self._save_results(filename, season, fevd, ret_values)

    def _compute_season(self, season_cube, model_cube):
        logger.info('GEV analysis...')
        shape = model_cube.shape
        units = model_cube.units
        fevd = dict()
        for par in self.gev_param_symbols:
            fevd[par] = np.full(shape, np.nan)

        ret_values = dict()
        for period in self.return_period:
            ret_values[period] = np.full(shape, np.nan)

        for i in range(shape[0]):
            for j in range(shape[1]):
                data = season_cube.data[..., i, j]
                if np.any(data):
                    self._compute_metric(
                        data, units, fevd, ret_values, i, j
                    )

        for par, data in fevd.items():
            fevd[par] = self._create_cube(data, par, model_cube)

        for period, data in ret_values.items():
            ret_values[period] = self._create_cube(
                data,
                self._get_return_period_name(period),
                model_cube
            )
        return fevd, ret_values

    def _compute_metric(self, data, units, fevd, ret_values, i, j):
        evdf = self.extremes.fevd(data, units=units.origin)
        results = evdf.rx2('results').rx2('par')
        # -ve mu/sigma invalid
        if results.rx2('location')[0] > 0. and results.rx2('scale')[0] > 0.:
            for par in self.gev_param_symbols:
                fevd[par][i, j] = results.rx2(self.gev_param_name[par])[0]
            r_level = self.extremes.return_level(
                evdf,
                return_period=self.return_period,
                qcov=self.extremes.make_qcov(evdf)
            )
            for values, period in zip(r_level, self.return_period):
                ret_values[period][i, j] = values

    @staticmethod
    def _create_cube(data, name, model_cube):
        cube = model_cube.copy(data)
        cube.standard_name = None
        cube.long_name = name
        cube.var_name = None
        return cube

    def _plot_results(self, filename, season, fevd, ret_values):
        if not self.cfg[n.WRITE_PLOTS]:
            return

        default_plot_config = {
            'plot_type': 'pcolormesh',
            'cmap': 'Blues',
            'coastlines': True,
        }
        logger.info('Plot results')
        results_subdir = os.path.join(
            self.cfg[n.PLOT_DIR],
            self.filenames.get_info(n.ALIAS, filename),
            season,
        )
        os.makedirs(results_subdir, exist_ok=True)
        for par in self.gev_param_symbols:
            par_ffp = os.path.join(
                results_subdir,
                '{}.{}'.format(
                    self.gev_param_name[par],
                    self.cfg[n.OUTPUT_FILE_TYPE]
                )
            )
            quickplot(
                fevd[par],
                filename=par_ffp,
                **(self.cfg.get('gev_quickplot', default_plot_config))
            )
        for return_period in ret_values:
            return_periods_path = os.path.join(
                results_subdir,
                'return_period_{}years.{}'.format(
                    return_period,
                    self.cfg[n.OUTPUT_FILE_TYPE]
                )
            )
            quickplot(
                ret_values[return_period],
                filename=return_periods_path,
                **(self.cfg.get(
                    'return_period_quickplot', default_plot_config))
            )

    def _save_results(self, filename, season, fevd, return_periods):
        if not self.cfg[n.WRITE_NETCDF]:
            return
        logger.info('Saving data...')
        results_subdir = os.path.join(
            self.cfg[n.WORK_DIR],
            self.filenames.get_info(n.PROJECT, filename),
            self.filenames.get_info(n.DATASET, filename),
            season,
        )
        os.makedirs(results_subdir, exist_ok=True)
        for par in self.gev_param_symbols:
            par_ffp = os.path.join(results_subdir, '{}.nc'.format(par))
            iris.save(fevd[par], par_ffp)

        return_periods_path = os.path.join(
            results_subdir, 'return_periods.nc'
        )
        iris.save(list(return_periods.values()), return_periods_path)


def main():
    """Execute analysis."""
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        ExtremeValuesAnalysis(config).compute()


if __name__ == "__main__":
    main()
