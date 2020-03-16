"""
Generalised extreme value analysis

Apply GEV analysis  grid point-by-grid point
"""

import os
import logging

import numpy as np
import iris
import iris.cube
import iris.analysis
import iris.util
import rpy2.robjects

import esmvaltool.diag_scripts.shared
from esmvaltool.diag_scripts.shared.plot import quickplot
import esmvaltool.diag_scripts.shared.names as n

logger = logging.getLogger(os.path.basename(__file__))


class ExtremePrecipitation(object):

    def __init__(self, config):
        self.cfg = config
        self.filenames = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.return_period = self.cfg['return_period']
        self.confidence_interval = self.cfg['confidence_interval']

        self.gev_param_symbols = ['mu', 'sigma', 'xi']
        self.gev_param_name = dict(mu='location', sigma='scale', xi='shape')
        self.return_period = rpy2.robjects.IntVector(self.return_period)

    def _get_return_period_name(self, period):
        return '{}-year level'.format(period)

    def compute(self):
        from rpy2.robjects import r as R, numpy2ri
        from rpy2.robjects.packages import importr

        extRemes = importr('extRemes')
        numpy2ri.activate()
        R['options'](warn=-1)
        for filename in self.filenames:
            logger.info('Processing %s', filename)
            cube = iris.load_cube(filename)
            logger.info(cube)

            model_cube = cube[0, ...]
            shape = model_cube.shape
            units = model_cube.units

            for season in set(cube.coord('clim_season').points):
                # Clean R objects
                R('rm(list = ls())')
                logger.info('Processing season %s', season)
                season_cube = cube.extract(iris.Constraint(clim_season=season))

                logger.info('GEV analysis...')
                fevd = dict()
                for par in self.gev_param_symbols:
                    fevd[par] = np.full(shape, np.nan)

                rl = dict()
                for period in self.return_period:
                    rl[period] = np.full(shape, np.nan)

                for x in range(shape[0]):
                    for y in range(shape[1]):
                        data = cube.data[..., x, y]
                        if np.any(data):
                            self._compute_metric(
                                data, units, fevd, rl, extRemes, x, y
                            )

                for par, data in fevd.items():
                    fevd[par] = self._create_cube(data, par, model_cube)

                for period, data in rl.items():
                    rl[period] = self._create_cube(
                        data,
                        self._get_return_period_name(period),
                        model_cube
                    )
                self._plot_results(filename, season, fevd, rl)
                self._save_results(filename, season, fevd, rl)

    def _compute_metric(self, data, units, fevd, rl, extRemes, x, y):
        evdf = extRemes.fevd(data, units=units.origin)
        results = evdf.rx2('results').rx2('par')
        # -ve mu/sigma invalid
        if results.rx2('location')[0] > 0. and results.rx2('scale')[0] > 0.:
            for par in self.gev_param_symbols:
                fevd[par][x, y] = results.rx2(self.gev_param_name[par])[0]
            r_level = extRemes.return_level(
                evdf,
                return_period=self.return_period,
                qcov=extRemes.make_qcov(evdf)
            )
            for data, period in zip(r_level, self.return_period):
                rl[period][x, y] = data

    def _create_cube(self, data, name, model_cube):
        cube = model_cube.copy(data)
        cube.standard_name = None
        cube.long_name = name
        cube.var_name = None
        return cube

    def _plot_results(self, filename, season, fevd, rl):
        if not self.cfg[n.WRITE_PLOTS]:
            return

        DEFAULT_PLOT_CONFIG = {
            'plot_type': 'pcolormesh',
            'cmap': 'Blues',
            'coastlines': True,
        }
        logger.info('Plot results')
        results_subdir = os.path.join(
            self.cfg[n.PLOT_DIR],
            self.filenames.get_info(n.PROJECT, filename),
            self.filenames.get_info(n.DATASET, filename),
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
                **(self.cfg.get('gev_quickplot', DEFAULT_PLOT_CONFIG))
            )
        for return_period in rl:
            return_periods_path = os.path.join(
                results_subdir,
                'return_period_{}years.{}'.format(
                    return_period,
                    self.cfg[n.OUTPUT_FILE_TYPE]
                )
            )
            quickplot(
                rl[return_period],
                filename=return_periods_path,
                **(self.cfg.get(
                    'return_period_quickplot', DEFAULT_PLOT_CONFIG))
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
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        ExtremePrecipitation(config).compute()


if __name__ == "__main__":
    main()
