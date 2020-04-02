"""
Diagnostic to evaluate the negative ice growth-ice thickness feedback

The codes presented here is derived from
TECLIM's GitHub code developed by F. Massonnet:
http://www.climate.be:3000/TECLIM/ClimateData.git
branch develop-fmasson
"""

import os
import logging
import math
import warnings
import numpy as np

import scipy.stats
import iris
import matplotlib.pyplot as plt

import esmvaltool.diag_scripts.shared
from esmvaltool.diag_scripts.shared import group_metadata
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared._base import ProvenanceLogger

logger = logging.getLogger(os.path.basename(__file__))


class NegativeSeaIceFeedback(object):
    """
    Diagnostic to evaluate the negative ice growth-ice thickness feedback

    Parameters
    ----------
    conf : dict
        Diagnostic execution info

    """

    def __init__(self, conf):
        self.cfg = conf
        self.datasets = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.variables = esmvaltool.diag_scripts.shared.Variables(self.cfg)

    def compute(self):
        """
        Compute diagnostic

        """
        negative_feedback = list()
        p_value = list()
        datasets = list()
        grouped_input_data = group_metadata(
            self.cfg['input_data'].values(), 'alias', sort='alias')
        for alias, dataset in grouped_input_data.items():
            try:
                feedback, p_val = self._compute_dataset(alias, dataset)
            except Exception as ex:
                logger.error('Failed to compute for %s', alias)
                logger.exception(ex)
            else:
                negative_feedback.append(feedback)
                p_value.append(p_val)
                datasets.append(alias)

        if self.cfg[n.WRITE_PLOTS]:
            self._plot_comparison(negative_feedback, datasets)
            self._plot_comparison(p_value, datasets, p_values=True)

    def _compute_dataset(self, alias, dataset):
        var_info = group_metadata(dataset, 'short_name')
        logger.info('Computing %s', alias)
        area_cello = iris.load_cube(
            var_info['areacello'][0]['filename']
        )
        cellarea = area_cello.data
        sit = iris.load_cube(var_info['sit'][0]['filename'])
        mask = np.asarray(
            sit.coord('latitude').points > 80.0,
            dtype=np.int8
        )
        try:
            mask = np.broadcast_to(mask, cellarea.shape)
        except ValueError:
            try:
                mask = np.broadcast_to(np.expand_dims(mask, -1),
                                       cellarea.shape)
            except ValueError:
                mask = np.broadcast_to(np.expand_dims(mask, 0),
                                       cellarea.shape)
        volume = self.compute_volume(sit, cellarea, mask=mask)
        del cellarea, sit

        neg_feedback, stats, _ = self.negative_seaice_feedback(
            var_info['sit'][0], volume, period=12, order=2
        )
        del volume
        logger.info("Negative feedback: %10.4f", neg_feedback)
        logger.info("P-Value:           %10.4f", stats[1])
        return (neg_feedback, stats[1])

    @staticmethod
    def compute_volume(avg_thick, cellarea, mask=1):
        """
        Compute sea ice volume

        Parameters
        ----------
        avg_thick : iris.Cube
            Sea ice or snow volume per unit cell area, in meters
        cellarea : [type]
            Grid cell area (sq. meters)
        mask : int, optional
            mask (1 on ocean, 0 on continent) (the default is 1)

        Raises
        ------
        ValueError
            avg_thick has not 2 nor 3 dimensions or mask not between 0 and 1

        Returns
        -------
        numpy.array
            Sea ice or snow volume in the region defined by the mask
        """
        if np.max(mask) != 1.0 or np.min(mask) < 0.0:
            raise ValueError("Mask not between 0 and 1")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            max_thick = avg_thick.collapsed(
                avg_thick.coords(), iris.analysis.MAX
            )
        if float(max_thick.data) > 20.0:
            logger.warning("Large sea ice thickness:"
                           "Max = %f",
                           max_thick.data)

        if avg_thick.coords('time'):
            vol = []
            for thick_slice in avg_thick.slices_over('time'):
                vol.append(
                    np.sum(
                        thick_slice.data * cellarea.data * mask.data
                    ) / 1e12
                )
            vol = np.asarray(vol)
        elif len(avg_thick.shape) == 2:
            vol = np.sum(avg_thick * cellarea * mask) / 1e12
        else:
            raise ValueError("avgthickness has not 2 nor 3 dimensions")
        return vol

    @staticmethod
    def detrend(data, order=1, period=None):
        """
        Detrend signal

        Parameters
        ----------
        data : numpy.array
            Data to detrend. Assumed to be sampled at evenly spaced times
        order : int, optional
            rder of the polynomial for detrending (the default is 1)
        period : int, optional
            possible existing periodicity of the signal, coming e.g.
            from external forcing,  expressed in units time steps.
            That is, data[i] and data[i + period] correspond
            to two realizations of the process at times where the
            forcing might be similar. Common examples include
            the seasonal cycle forcing, or the diurnal forcing.

            If "period" is not None, the detrending is performed
            separately for each time step (e.g., all 1st of January,
            all 2nd of January, ..., all 31st of December in case
            of annual cycle.

            If "period" is None, the detrending is performed on
            the given time series
            The default is None

        Raises
        ------
        ValueError
            [description]

        Returns
        -------
        numpy.array
            the signal detrended using a least-square polynomial
            regression of order "order"
        """

        if len(data.shape) != 1:
            raise ValueError("Non-conform input data")

        # Remove possible nans from the data. All the regression (ie polyfit)
        # parameters will be estimated based on the no-nan data but the
        # residuals will be computed from the original data in order
        # to keep the same size and to restitute NaNs where they appeared

        data_nonan = data[~np.isnan(data)]

        # If the signal has no periodicity, we just make a linear regression
        if period is None:
            time_nonan = np.arange(len(data_nonan))
            time = np.arange(len(data))
            polynom = np.polyfit(time_nonan, data_nonan, order)
            residuals = data - np.sum([polynom[i] * time ** (order - i)
                                       for i in range(order + 1)], axis=0)

        # If the signal contains a periodical component, we do the regression
        # time step per time step
        else:
            residuals = np.empty([n])

            # For each time step of the period, detrend
            # Note that another common option is to first remove a seasonal
            # cycle and then detrend the anomalies. However this assumes that
            # a cycle can be estimated, which in presence of a trend is tricky
            # because the trend component interferes with the mean. I have
            # tried that and it gives ugly step-wise anomalies. Detrending day
            # per day seems the most natural way to do, at least as long as we
            # assume that the raw signal at some time is the result of a
            # seasonal cycle depending on the position of the time step in the
            # period, plus a common trend, plus some noise.
            for i in np.arange(period):
                raw = data[np.arange(i, n, period)]
                raw_nonan = raw[~np.isnan(raw)]
                time = np.arange(len(raw))
                time_nonan = np.arange(len(raw_nonan))
                polynom = np.polyfit(time_nonan, raw_nonan, order)
                residuals[np.arange(i, n, period)] = \
                    raw - np.sum([polynom[i] * time ** (order - i)
                                  for i in range(order + 1)], axis=0)
        return residuals

    def negative_seaice_feedback(self, dataset_info, volume, period, order=1):
        """
        Function to estimate the negative ice-thickness ice growth feedback
        and its significance.

        Parameters
        ----------
        volume : iris.Cube
            Time series of sea ice volume
        period : int
            period of the signal (period expressed in time steps:
            12 for monthly, 365 for daily...)
        order : int, optional
            order of the polynomial detrending (>=0) (the default is 1)

        Raises
        ------
        ValueError
            If volume is not 1D or timeseries length is not a multiple of
            the period

        Returns
        -------
        nv: float
            Feedback parameter expressed as the regression
            between dV on V_min

        corr: tuple(float, float, float)
            Correlation between those two, the p-value under the null
            hypothesis of no correlation between dV and V_min and standard
            deviation

        volume: tuple(float, float(
            [V_min, dV]: detrended time series of annual minimum of sea ice
            volume, detrended series of wintertime volume production
        """

        if len(volume.shape) != 1:
            raise ValueError("Volume is not 1-D")

        if volume.size % period != 0:
            raise ValueError(
                "Length of volume series is not multiple of period"
            )

        # 1. Locate the minima for each year
        imin = [t + np.nanargmin(volume.data[t:t + period])
                for t in range(0, volume.size, period)]

        # 2. Locate the maxima for each year
        imax = [t + np.nanargmax(volume.data[t:t + period])
                for t in np.arange(0, volume.size, period)]

        # 3. Detrend series. A one-year shift is introduced to make sure we
        #    compute volume production *after* the summer minimum
        vol_min = self.detrend(volume[imin[:-1]], order=order)
        dvol = self.detrend(volume[imax[1:]] - volume[imin[:-1]], order=order)

        # 4. Compute diagnostics
        # If all Vmins are zero or all dVs are zero, return Nan
        # (pathological case)
        if np.max(vol_min) == 0.0 or np.max(dvol == 0.0):
            fit = np.nan
            fit_complete = [np.nan, np.nan]
            corr = np.nan
            pval = np.nan
            std = np.nan
        else:
            corr = np.corrcoef(vol_min, dvol)[0, 1]
            # The t-statistic.
            tstat = corr / np.sqrt((1 - corr ** 2) / (len(vol_min) - 2))
            # Under the null hypothesis of no correlation,
            # tstat follows a student's law with  N - 2 dof.
            pval = 1.0 - scipy.stats.t.cdf(np.abs(tstat), len(vol_min) - 2)

            if pval > 0.05:
                logger.warning(
                    "Check the scatterplot of dV versus V_min, it is most "
                    "likely suspicious, and the feedback factor likely "
                    "meaningless: p-value: %f", pval
                )

            try:
                fit_complete, cov = np.polyfit(vol_min, dvol, 1, cov=True)
                logger.info(fit_complete)
                fit = fit_complete[0]  # Fit parameter
                std = np.sqrt(cov[0, 0])  # Standard deviation on it
            except ValueError:
                logger.error("(negative_seaice_feedback) PROBLEM,"
                             "series badly conditioned: "
                             "Input volume: %f Vmin: %f dv: %f",
                             volume, vol_min, dvol)
                raise
        if self.cfg[n.WRITE_PLOTS]:
            self._plot_ife(dataset_info, vol_min, dvol, fit_complete)
        return [fit, [corr, pval, std], [vol_min, dvol]]

    def _plot_ife(self, dataset_info, vol_min, dvol, fit_complete):
        path = os.path.join(
            self.cfg[n.PLOT_DIR],
            f'ife_{dataset_info[n.ALIAS]}.{self.cfg[n.OUTPUT_FILE_TYPE]}'
        )
        plot_options = self.cfg.get('plot', {})
        fig = plt.figure()
        plt.scatter(
            vol_min,
            dvol,
            plot_options.get('point_size', 8),
            color=plot_options.get('point_color', 'black'),
        )
        minx, maxx = plt.xlim()
        xvals = np.linspace(minx, maxx)
        plt.plot(xvals, xvals * fit_complete[0] + fit_complete[1])
        axes = plt.gca()
        axes.set_title(
            f'Evaluation of the IFE \n{dataset_info[n.ALIAS]} '
            f'({dataset_info[n.START_YEAR]}-{dataset_info[n.END_YEAR]})'
        )
        axes.set_ylabel('Wintertime volume range \n(anomalies) [10³ km³]')
        axes.set_xlabel('Volume at minimum\n(anomalies) [10³ km³]')
        plt.grid(True, 'both', 'both')
        plt.tight_layout()
        fig.savefig(path)
        plt.close(fig)
        self._create_prov_record(
            path, f'Evaluation of IFE for {dataset_info[n.ALIAS]}',
            [info['filename'] for info in group_metadata(
                self.cfg['input_data'].values(), n.ALIAS
            )[dataset_info[n.ALIAS]]
            ]
        )

    def _plot_comparison(self, data, datasets, p_values=False):
        if p_values:
            filename = 'feedback_p_values'
        else:
            filename = 'feedback'

        path = os.path.join(
            self.cfg[n.PLOT_DIR],
            f'{filename}.{self.cfg[n.OUTPUT_FILE_TYPE]}'
        )

        plot_options = self.cfg.get('plot', {})
        fig = plt.figure()
        index = np.arange(len(data))
        plt.scatter(
            index,
            data,
            plot_options.get('point_size', 8),
            color=plot_options.get('point_color', 'black'),
        )
        if p_values:
            plt.hlines(0.05, -1, index[-1] + 1, colors='red')
        axes = plt.gca()
        logger.debug(data)
        max_limit = math.ceil(max(data))
        if max_limit < 0:
            max_limit = 0
        min_limit = math.floor(min(data))
        separation = max_limit - min_limit

        if plot_options.get('show_values', False):
            def _get_y_position(value):
                if value > min_limit + separation * 0.75:
                    return value - separation * 0.05
                return value + separation * 0.10

            for i, value in enumerate(data):
                axes.annotate(
                    f'{value:.2f}',
                    xy=(index[i], value),
                    xycoords='data',
                    textcoords='data',
                    xytext=(index[i], _get_y_position(value)),
                    rotation=90,
                )

        # axes and labels
        axes.set_ylim(min_limit, max_limit)
        if p_values:
            axes.set_ylabel('P-value [log]')
            plt.ylim(0, max(0.25, max(data)))
        else:
            axes.set_ylabel('IFE')
        axes.set_title('IFE comparison')
        _, xtick_names = plt.xticks(index, datasets)
        plt.xlim(index[0] - 0.5, index[-1] + 0.5)
        plt.setp(xtick_names, rotation=90, fontsize=10)
        plt.grid(True, 'both', 'y')
        plt.tight_layout()
        fig.savefig(path)
        plt.close(fig)
        self._create_prov_record(
            path, f'IFE {filename} comparison for all datasets',
            group_metadata(self.cfg['input_data'].values(), n.ALIAS)
        )

    def _create_prov_record(self, filepath, caption, ancestors):
        record = {
            'caption': caption,
            'domains': ['nhpolar'],
            'ancestors': ancestors
        }
        with ProvenanceLogger(self.cfg) as provenance_logger:
            provenance_logger.log(filepath, record)


def main():
    """Run diagnostic"""
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        NegativeSeaIceFeedback(config).compute()


if __name__ == '__main__':
    main()
