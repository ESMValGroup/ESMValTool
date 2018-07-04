import os
import csv
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from geoval.core.mapping import SingleMap
from ESMValMD import ESMValMD
from diagnostic import BasicDiagnostics


class SoilMoistureDiagnostic(BasicDiagnostics):
    """
    class to implement soil moisture diagnostics,
    like e.g. global means, global differences, RMSD etc.

    TODO implement testing for this diagnostic

    """

    def __init__(self, **kwargs):
        super(SoilMoistureDiagnostic, self).__init__(**kwargs)

        self._project_info = {}
        self._modtype = None
        self._reftype = None
        self._plot_dir = '.' + os.sep
        self._work_dir = '.' + os.sep

        self._vartype = 'soil moisture'  # default value as there must be one
        self.output_type = 'png'  # default ouput file type
        self._changed = False

#        self._diag_name = 'Soil Moisture'
#        self.server = 'login.zmaw.de'
#        self.remote_dir = '/data/icdc/land/esa_cci_soilmoisture/DATA/'
#        self.user = 'm300028'
#        self.raw_reference_directory = './cci_sm/raw'
#        self.reference_preproc_directory = './cci_sm/preprocessed'
#        self.targetgrid = kwargs.pop('targetgrid', 't63grid')
#        self._projected_files = []

    def run_diagnostic(self):
        """
        running the diagnostics
        """
        super(SoilMoistureDiagnostic, self).run_diagnostic()

        self._specific_diag()

    def _specific_diag(self):
        """
        Diagnostic management
        """

        if "percentile" in self.cfg.__dict__.keys() and self.cfg.percentile:
            self._percentile_comparison(
                plist=np.arange(self.cfg.percentile_pars[0],
                                self.cfg.percentile_pars[1] +
                                0.1 * self.cfg.percentile_pars[2],
                                self.cfg.percentile_pars[2]))

        if "anomaly" in self.cfg.__dict__.keys() and self.cfg.anomaly:
            self._anomaly_correlation()

    def write_data(self, plot=True):
        """
        write data
        """
        super(SoilMoistureDiagnostic, self).write_data()

        if self.cfg.regionalization and '_regions' not in self.__dict__.keys():
            self._regions = self._mod_data.get_regions(
                self._reg_shape, self.cfg.shapeNames - 1)
            self._write_regionalization_header()
            print "this should not be calculated here"

        if '_percentile_list' in self.__dict__.keys():
            # write statistic file
            plist = np.array(np.arange(self.cfg.percentile_pars[0],
                                       self.cfg.percentile_pars[1] +
                                       0.1 * self.cfg.percentile_pars[2],
                                       self.cfg.percentile_pars[2]))

            self._write_percentile_correlations(plist)
            self._plot_percentile_correlation(plist, np.array(self._r_list))

            for p in np.arange(len(self._percentile_list)):
                # generate map plots for percentiles
                self._plot_percentile_maps(plist[p],
                                           self._percentile_list[p][0],
                                           self._percentile_list[p][1],
                                           np.array(self._r_list)[p])
                if self.cfg.regionalization:
                    self._percentile_list[p][0].get_shape_statistics(
                        self._regions)
                    self._percentile_list[p][1].get_shape_statistics(
                        self._regions)
                    self._write_shape_statistics(
                        self._percentile_list[p][0].regionalized,
                        'percentile_' + str(int(plist[p] * 100)).zfill(3),
                        self.refname)
                    self._write_shape_statistics(
                        self._percentile_list[p][1].regionalized,
                        'percentile_' + str(int(plist[p] * 100)).zfill(3),
                        self.modname)
        else:
            print 'No percentile data to plot!'

        if '_ref_anomaly' and '_mod_anomaly' in self.__dict__.keys():
            self._plot_anomaly_correlation(
                self._ref_anomaly, self._RA_pval,
                self._mod_anomaly, self._RM_pval)
            if self.cfg.regionalization:
                self._ref_anomaly.get_shape_statistics(self._regions)
                self._mod_anomaly.get_shape_statistics(self._regions)
                self._write_shape_statistics(
                    self._ref_anomaly.regionalized,
                    'anomaly_correlation', self.refname)
                self._write_shape_statistics(
                    self._mod_anomaly.regionalized,
                    'anomaly_correlation', self.modname)
        else:
            print 'No anomaly correlation to plot!'

    def _write_percentile_correlations(self, plist):
        """ writing percentile correlations as csv """
        oname = self._get_output_rootname() + '_percentile_correlation.csv'
        if os.path.exists(oname):
            os.remove(oname)
        f = open(oname, 'w')
        try:
            writer = csv.writer(f, quoting=csv.QUOTE_NONNUMERIC)
            writer.writerow(('percentile', 'correlation'))
            writer.writerows(zip(np.char.array(plist), np.array(self._r_list)))
        finally:
            f.close()

        ESMValMD("xml",
                 oname,
                 self._basetags + ['DM_global', 'M_' + self.refname]
                                + ['M_{0}'.format(str(item)) for item in self.modelnames]
                                + ['ST_corr', 'ST_perc'],
                 'Development of global pattern correlation values over ' +
                 'different percentile levels for ' + self.refname + ' and ' +
                 self.modname + ' ' + self._vartype + ' data.',
                 '#ID' + 'devcorrperctab' + self.var,
                 ','.join(self._infiles))

    def _percentile_comparison(self, plist=np.arange(0.0, 1.01, 0.05)):
        """
        calculate percentiles for model and observational dataset
        and compare these

        Parameters
        ----------
        plist : list
            list of percentile values to analyze. A value of e.g. 0.05
            corresponds to a 5% percentile

        References
        ----------
        * Loew et al. (2013): doi: 10.5194/hess-17-3523-2013

        """

        print('   percentile analysis ...')

        # calculate for each percentile the spatial distribution for both
        # the model and reference data
        self._r_list = []
        self._percentile_list = []
        for p in plist:
            pmod = self._mod_data.get_percentile(p)  # model percentile map
            pref = self._ref_data.get_percentile(p)  # ref data percentile map

            # TODO I cannot find the error... data behaves weird
            # THIS IS HARDCODED!
            perc_mask = ((pmod.data.data > 1.5) & (pref.data.data > 1.5))
            pref.data.mask = pref.data.mask | perc_mask
            pmod.data.mask = pmod.data.mask | perc_mask

            self._percentile_list.append([pmod, pref])

            # calculate spatial correlation
            r_value, p_value = self._calc_spatial_correlation(pmod, pref)
            self._r_list.append(r_value)

    def _plot_percentile_maps(self, p, mod, ref, r):
        """
        plot percentile maps
        """
        f = plt.figure(figsize=(20, 6))
        ax1 = f.add_subplot(121)
        ax2 = f.add_subplot(122)

        def submap(data, ax, title, vmin, vmax, cmap,
                   ctick={'ticks': None, 'labels': None}):
            Map = SingleMap(data,
                            backend=self.plot_backend,
                            show_statistic=True,
                            savefile=None,
                            ax=ax,
                            show_unit=True)
            Map.plot(title=title,
                     show_zonal=False,
                     show_histogram=False,
                     show_timeseries=False,
                     nclasses=self.plot_nclasses,
                     colorbar_orientation=self.plot_cborientation,
                     show_colorbar=self.plot_cbshow,
                     cmap=cmap,
                     vmin=min(self.cfg.mima_percentile),
                     vmax=max(self.cfg.mima_percentile),
                     proj_prop=self.cfg.projection,
                     ctick_prop=ctick,
                     drawparallels=True,
                     titlefontsize=self.plot_tfont)

        submap(ref, ax=ax1, title=self.refname, vmin=0, vmax=1, cmap='Blues')
        submap(mod, ax=ax2, title=self.modname, vmin=0, vmax=1, cmap='Blues')
        f.suptitle('Percentile: ' + str(int(p*100)) + "% (r=" +
                   str(round(r, 2)) + ")")

        oname = self._get_output_rootname() + '_percentile_' + \
            str(int(p * 100)).zfill(3) + '.' + self.output_type
        if os.path.exists(oname):
            os.remove(oname)
        f.savefig(oname)

        plt.close(f.number)  # close figure for memory reasons!
        del f

        ESMValMD("both",
                 oname,
                 self._basetags + ['DM_global', 'PT_geo', 'M_' + self.refname]
                                + ['M_{0}'.format(str(item)) for item in self.modelnames]
                                + ['ST_perc'],
                 'Comparison of global patterns of ' + self._vartype +
                 ' for ' + str(int(p * 100)) + 'th-percentile of ' +
                 self.refname + ' and ' + self.modname + ' data. ' +
                 'The spatial correlation (r) is noted in the title.',
                 '#ID' + 'perc' + str(int(p * 100)).zfill(3) + self.var,
                 ','.join(self._infiles))

    def _plot_percentile_correlation(self, p, r):
        """
        plot correlation as function of percentile

        Parameters
        ----------
        p : ndarray
            percentiles
        r : ndarray
            pearson correlation coefficient
        """
        oname = self._get_output_rootname() + \
            '_percentile_spatial_correlation' + \
            '.' + self.output_type
        if os.path.exists(oname):
            os.remove(oname)

        f = plt.figure()
        f.suptitle(self.refname + "-" + self.modname +
                   ' percentile spatial correlation', fontsize=14)
        ax = f.add_subplot(111)
        ax.plot(p*100, r, linestyle='-', marker='d')
        ax.set_xlabel('percentile [%]')
        ax.set_ylabel('correlation')
        ax.grid()
        ax.set_ylim(-1., 1.)
        f.savefig(oname)
        plt.close(f.number)

        ESMValMD("both",
                 oname,
                 self._basetags + ['DM_global', 'PT_pro', 'M_' + self.refname]
                                + ['M_{0}'.format(str(item)) for item in self.modelnames]
                                + ['ST_corr', 'ST_perc'],
                 'Development of global pattern correlation values over ' +
                 'different percentile levels for ' + self.refname + ' and ' +
                 self.modname + ' ' + self._vartype + ' data.',
                 '#ID' + 'devcorrperc' + self.var,
                 ','.join(self._infiles))

    def _anomaly_correlation(self):
        """
        Anomaly correlation analysis with precipitation
        """

        print('   anomaly analysis ...')

        # read additional precipitation data
        _aux_pr_data = self._mod_pr_data

        _aux_pr_data.apply_temporal_subsetting(
            self._start_time, self._stop_time)

        # calculate anomalies from precipitation data
        anomaly_pr = _aux_pr_data.get_deseasonalized_anomaly(base='current')

        # calculate anomalies from model soil moisture data
        # There are nans in the respective clim object during calculation of
        # anomalies
        anomaly_mod = self._mod_data.get_deseasonalized_anomaly(base='current')

        # calculate anomalies from reference soil moisture data
        anomaly_ref = self._ref_data.get_deseasonalized_anomaly(base='current')

        # correlate anomalies from precipitation with soil moistures
        self._ref_anomaly, self._RA_pval = self._mapping_tau(
            anomaly_ref, anomaly_pr)
        self._mod_anomaly, self._RM_pval = self._mapping_tau(
            anomaly_mod, anomaly_pr)

    def _plot_anomaly_correlation(self, anomalyref, pref, anomalymod, pmod):
        """
        plot anomaly maps
        """
        f = plt.figure(figsize=(20, 10))
        ax1 = f.add_subplot(221)
        ax2 = f.add_subplot(222)
        ax3 = f.add_subplot(223)
        ax4 = f.add_subplot(224)

        def submap(data, ax, title, vmin, vmax, cmap,
                   ctick={'ticks': None, 'labels': None}):
            Map = SingleMap(data,
                            backend=self.plot_backend,
                            stat_type='mean',
                            savefile=None,
                            ax=ax,
                            show_unit=False)
            Map.plot(title=title,
                     show_zonal=False,
                     show_histogram=False,
                     show_timeseries=False,
                     nclasses=self.plot_nclasses,
                     colorbar_orientation=self.plot_cborientation,
                     show_colorbar=self.plot_cbshow,
                     cmap=cmap,
                     vmin=vmin,
                     vmax=vmax,
                     proj_prop=self.cfg.projection,
                     ctick_prop=ctick,
                     drawparallels=True,
                     titlefontsize=self.plot_tfont)

        submap(anomalyref, ax=ax1, title=self.refname,
               vmin=-1, vmax=1, cmap='RdBu')
        submap(anomalymod, ax=ax2, title=self.modname,
               vmin=-1, vmax=1, cmap='RdBu')

        tick_def = np.arange(0, 1.01, 0.1)
        submap(pref, ax=ax3, title="p-value " + self.refname, vmin=0, vmax=1,
               cmap='summer',
               ctick={'ticks': tick_def,
                      'labels': np.append(tick_def[:-1].astype('string'),
                                          '> 1.0')})
        submap(pmod, ax=ax4, title="p-value " + self.modname, vmin=0, vmax=1,
               cmap='summer',
               ctick={'ticks': tick_def,
                      'labels': np.append(tick_def[:-1].astype('string'),
                                          '> 1.0')})
        f.suptitle(self.refname + " and " + self.modname +
                   " soil moisture anomaly correlations with precipitation")

        oname = self._get_output_rootname() + '_anomaly_with_percipitation' + \
            '.' + self.output_type
        if os.path.exists(oname):
            os.remove(oname)
        f.savefig(oname)

        plt.close(f.number)  # close figure for memory reasons!
        del f

        ESMValMD("both",
                 oname,
                 self._basetags + ['pr', 'ST_anomaly', 'ST_correlation',
                                   'DM_global', 'PT_geo', 'M_' + self.refname]
                                + ['M_{0}'.format(str(item)) for item in self.modelnames],
                 'Pixelwise anomaly correlation with precipitation for ' +
                 self.modname + ' and ' + self.refname + ' ' + self._vartype +
                 ' data. The left column shows reference data and the right ' +
                 'column shows model data, while top row shows correlation ' +
                 'values and bottom row shows p-values, accordingly. The ' +
                 'p-values higher than 1.0 are not shown separately.',
                 '#ID' + 'ACorr' + self.var + 'pr',
                 ','.join([f for f in self._allfiles if "Shapefiles" not in f])
                 )

    def _load_model_data(self):
        """ load soil moisture model data """
        orig_mod_file = self._mod_file

        mod_info = Dataset(self._mod_file)
        try:
            lat = mod_info.dimensions['lat'].size
            lon = mod_info.dimensions['lon'].size
        except:  # regridding required in any case
            lat = -1
            lon = -1
        mod_info.close()

        if not (lat == self.reso["lat"] and lon == self.reso["lon"]):

            grid = self.resoLUT[str(self.reso["lat"]) + "-" +
                                str(self.reso["lon"])]

            newfile = self._mod_file + "." + grid + "built.nc"
            newfile = newfile.split("/")
            newdir = (self._work_dir if self._work_dir[-1] ==
                      os.sep else self._work_dir + os.sep) +\
                "AUX_Files_sm_ESACCI"
            newfile = newdir + os.sep + newfile[-1]

            if not os.path.exists(newfile):
                tempfile = self._aggregate_resolution(
                    self._mod_file, grid, remove=False)
                if not os.path.exists(newdir):
                    os.makedirs(newdir)
                shutil.copy2(tempfile, newfile)
                os.remove(tempfile)

            self._mod_file = newfile

        # load data
        self._mod_data = self._load_cmip_generic(
                self._mod_file, self._project_info['RUNTIME']['currDiag'].
                get_variables()[0])

        if "anomaly" in self.cfg.__dict__.keys() and self.cfg.anomaly:
        # A_laue_ax   if self.cfg.anomaly:

            self._mod_pr_file = orig_mod_file
            self._mod_pr_file = self._mod_pr_file.split("/")
            self._mod_pr_file[-1] = self._mod_pr_file[-1].replace(
                "_sm_", "_pr_").replace("_Lmon_", "_Amon_")
            self._mod_pr_file = "/".join(self._mod_pr_file)

            mod_info = Dataset(self._mod_pr_file)
            try:
                lat = mod_info.dimensions['lat'].size
                lon = mod_info.dimensions['lon'].size
            except:  # regridding required in any case
                lat = -1
                lon = -1
            mod_info.close()

            if not (lat == self.reso["lat"] and lon == self.reso["lon"]):

                grid = self.resoLUT[str(self.reso["lat"]) + "-" +
                                    str(self.reso["lon"])]

                newfile = self._mod_pr_file + "." + grid + "built.nc"
                newfile = newfile.split("/")
                newdir = (self._work_dir if self._work_dir[-1] ==
                          os.sep else self._work_dir + os.sep) +\
                    "AUX_Files_sm_ESACCI"
                newfile = newdir + os.sep + newfile[-1]

                if not os.path.exists(newfile):
                    tempfile = self._aggregate_resolution(
                        self._mod_pr_file, grid, remove=False)
                    if not os.path.exists(newdir):
                        os.makedirs(newdir)
                    shutil.copy2(tempfile, newfile)
                    os.remove(tempfile)

                self._mod_pr_file = newfile

            # load data
            self._mod_pr_data = self._load_cmip_generic(
                    self._mod_pr_file, "pr")

    def _load_observation_data(self):
        """ load obs data """

        mod_info = Dataset(self._ref_file)
        try:
            lat = mod_info.dimensions['lat'].size
            lon = mod_info.dimensions['lon'].size
        except:  # regridding required in any case
            lat = -1
            lon = -1
        mod_info.close()

        reso = {"lat": lat, "lon": lon}

        if str(lat) + "-" + str(lon) not in self.resoLUT.keys():

            grid = self.resoLUT[str(self.reso["lat"]) + "-" +
                                str(self.reso["lon"])]

            newfile = self._ref_file + "." + grid + "built.nc"
            newfile = newfile.split("/")
            newdir = (self._work_dir if self._work_dir[-1] ==
                      os.sep else self._work_dir + os.sep) +\
                "AUX_Files_sm_ESACCI"
            newfile = newdir + os.sep + newfile[-1]

            if not os.path.exists(newfile):
                tempfile = self._aggregate_resolution(
                    self._ref_file, grid, remove=False)
                if not os.path.exists(newdir):
                    os.makedirs(newdir)
                shutil.copy2(tempfile, newfile)
                os.remove(tempfile)

            self._ref_file = newfile

            reso = self.reso

        self.reso = reso

        proj_var = self._project_info['RUNTIME']['currDiag'].get_variables()[0]
        if proj_var == "sm":
            self._ref_data = self._load_cci_generic(self._ref_file, proj_var)
        else:
            assert False, 'Not supported yet'
