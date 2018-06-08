'''
Module to deal with plotting Normalised Assessment Criteria Plots
'''

import os

import matplotlib.pyplot as mplt
import numpy as np
import numpy.ma as ma

# Define some colours
BLACK = '#000000'
RED = '#FF0000'
AMBER = '#FF8C00'
GREEN = '#7CFC00'
OBS_GREY = '#000000'
ACC_GREY = '#00FFFF'
STD_GREY = '#EEEEEE'
NOOBS_GREY = '#A9A9A9'

# Set available markers
# Those towards end of list are not really suitable but do extend list
# - How about custom symbols?
MARKERS = 'ops*dh^v<>+xDH.,'

# Create fakelines for legend using MARKERS above
FAKELINES = [mplt.Line2D([0, 0], [0, 1], marker=marker, color=BLACK,
                         linestyle='') for marker in MARKERS]

class MetricData:
    '''Class to hold data from a metrics file'''

    def __init__(self, obj, desc=None, mdi=-10000.0, normalised=False):
        '''
        Initialise object with filename and optional MDI value
        Read file to get initialise data
        '''
        self.desc = desc
        self.mdi = mdi
        self.normalised = normalised

        # What inputs are acceptable?
        #  1) String filename
        #  2) Numpy record array
        #  3) Metrics dictionary
        #  4) List of metric names

        if isinstance(obj, list):
            # If object is list then convert to records array just
            # to get label names
            recs = [(label,) for label in obj]
            self._data = np.core.records.fromrecords(recs)
        elif isinstance(obj, dict):
            # If object is dictionary then convert to records array by writing
            # data into records format
            recs = [tuple([label] + list(value))
                    for (label, value) in obj.items()]
            self._data = np.core.records.fromrecords(recs)
        elif isinstance(obj, np.core.records.recarray):
            # If object is records array then just copy into data
            self._data = obj
        elif isinstance(obj, basestring):
            # If object is a string then assume it is a filename
            self.filename = obj
            try:
                print 'Creating metrics from file {}'.format(self.filename)
                with open(self.filename) as inf:
                    self._data = np.recfromtxt(inf,
                                               delimiter=',',
                                               autostrip=True)
            except IOError:
                print 'File {} not found'.format(self.filename)

    def __eq__(self, other):
        '''
        Test "equality" of instances. Actually tests to make sure that
        objects have same metrics in same order.
        If either object is blank, i.e. _data is None, then get False.
        '''
        assert isinstance(other, MetricData)
        samesize = False
        sameorder = False
        if self.haslabel and other.haslabel:
            samesize = (self.label.size == other.label.size)
            if samesize:
                sameorder = np.all(self.label == other.label)
        return (samesize and sameorder)

    def __ne__(self, other):
        '''Inverse of "equality"'''
        return not (self == other)

    def __str__(self):
        '''String representation of metric data'''
        printlines = []
        printlines.append('description: {}'.format(self.desc))
        printlines.append('normalised: {}'.format(self.normalised))
        printlines.append('mdi: {}'.format(self.mdi))
        if hasattr(self, 'filename'):
            printlines.append('filename: {}'.format(self.filename))
        printlines.append('label: {}'.format(self.label))
        printlines.append('data: {}'.format(self.data))
        return '\n'.join(printlines)

    def __repr__(self):
        '''Write out representation as string'''
        return str(self)

    @property
    def label(self):
        '''Returns labels as ndarray from records array'''
        if hasattr(self, '_data'):
            if self._data.dtype.names:
                dtype_names = list(self._data.dtype.names)
                names = self._data[dtype_names[0]]
            else:
                names = self._data
            return np.array(names)
        else:
            return None

    @property
    def data(self):
        '''Returns data as ndarray from records array'''
        if hasattr(self, '_data'):
            if self._data.dtype.names:
                dtype_names = list(self._data.dtype.names)
                tmpdata = tuple([self._data[name] for name in dtype_names[1:]])
                return ma.masked_values(np.vstack(tmpdata), self.mdi)
            else:
                return None
        else:
            return None

    @property
    def hasdata(self):
        '''Determine if object has data'''
        return self.data is not None

    @property
    def haslabel(self):
        '''Determine if object has list of metrics'''
        return self.label is not None

    @property
    def plotready(self):
        '''Returns if data object is ready to be plotted'''
        data_not_exist = not self.hasdata
        data_exist_and_norm = self.hasdata and self.normalised
        return data_not_exist or data_exist_and_norm

    def filter(self, labels, nonmdi=False, nonzero=False):
        '''
        Filter data to only those fields given in labels (in that order)
        Can also filter out model data that is zero or missing
        '''
        # Filter out by label and sort in given order
        # If metric in label but not in self then add missing data
        if isinstance(labels, MetricData):
            labels = labels.label
        if (labels is not None) and self.haslabel:
            newrecs = []
            for lab in labels:
                if lab in self.label:
                    index = np.argwhere(self.label == lab)[0]
                    newrecs.extend(self._data[index].tolist())
                else:
                    rec = list(self._data[0].tolist())
                    newrec = [lab] + [self.mdi for val in rec[1:]]
                    newrecs.append(tuple(newrec))
            self._data = np.core.records.fromrecords(newrecs)

        # Only do data filters with shape[0]==1
        if (self.hasdata) and (self.data.shape[0] == 1):
            # Filter out missing data
            if nonmdi and np.any(self.data.mask):
                indices = np.argwhere(np.logical_not(self.data.mask[0]))[:, 0]
                self._data = self._data[indices]
            # Filter out zero fields
            if nonzero:
                indices = np.argwhere(self.data.data[0])[:, 0]
                self._data = self._data[indices]

    def normalise(self, norm, scale=1.0):
        '''Normalise data using given norm'''
        try:
            assert self.hasdata, "no data in data object"
            assert not self.normalised, "data already normalised"
            assert self == norm, "data objects do not match"
        except AssertionError, msg:
            print 'Cannot normalise, {}.'.format(msg)
            return
        # Normalise data in numpy ndarray form
        data = self.data / norm.data
        data *= scale
        # Write data into records format and convert to records array
        recs = [tuple([label] + list(data[:, j].filled()))
                for (j, label) in enumerate(self.label)]
        self._data = np.core.records.fromrecords(recs)
        self.normalised = True
        return

    def addstd(self, std):
        '''Add model uncertainty to data'''
        try:
            assert self.hasdata, "no data in data object"
            assert self.normalised, "data not normalised"
            assert self == std, "data objects do not match"
            assert std.hasdata, "no data in std data object"
            assert std.normalised, "std data not normalised"
        except AssertionError, msg:
            print 'Cannot normalise, {}.'.format(msg)
            return
        # Add std data in numpy ndarray form
        data = self.data + std.data
        # Write data into records format and convert to records array
        recs = [tuple([label] + list(data[:, j].filled()))
                for (j, label) in enumerate(self.label)]
        self._data = np.core.records.fromrecords(recs)
        return

    def split_obs(self):
        '''
        Split observations data into observational uncertainty and acceptable
        range. Acceptable range held in new object.
        '''

        # Split data into obs and acc
        if self.hasdata:
            # If acceptable range present
            if self.data.shape[0] == 4:
                newrecs = [(rec['f0'], rec['f3'], rec['f4'])
                           for rec in self._data]
                acc_data = np.core.records.fromrecords(newrecs)
            else:
                acc_data = None
            # If observational uncertainty range present
            if self.data.shape[0] >= 2:
                newrecs = [(rec['f0'], rec['f1'], rec['f2'])
                           for rec in self._data]
                obs_data = np.core.records.fromrecords(newrecs)
            # If only single observational uncertainty value present
            else:
                newrecs = [(rec['f0'], 0.0, rec['f1'])
                           for rec in self._data]
                obs_data = np.core.records.fromrecords(newrecs)
        else:
            obs_data = None
            acc_data = None

        obs = MetricData(obs_data, desc='Observational Uncertainty')
        acc = MetricData(acc_data, desc='Acceptable Range')

        return (obs, acc)

    # Need to make this routine obsolete, should be able to configure
    # obs ranges so that this doesn't matter.
    def absmod(self):
        '''Method to modify obs ranges if no range specified'''
        # If absolute plotting then copy top range into bottom range
        # Assumes bottom of range is zero
        nrecs = [(rec['f0'], rec['f2'], rec['f2']) for rec in self._data]
        self._data = np.core.records.fromrecords(nrecs)

    def min(self):
        '''Shortcut to data minimum'''
        return self.data.min()

    def max(self):
        '''Shortcut to data maximum'''
        return self.data.max()

    def absmax(self):
        '''Shortcut to abs(data) maximum'''
        return np.fabs(self.data).max()

    def abssum(self, dim):
        '''Shortcut to abs(data) sum'''
        return np.fabs(self.data[dim, :]).sum()


class NACPlot:
    '''
    Class to construct data for Normalised Assessment Criteria plot
    '''

    def __init__(self, norm, expts,
                 file_nrm, file_mod, file_obs, file_std, file_sub,
                 absolute=False, extend_y=False, no_min_exty=False,
                 alt=False):
        '''Initialise object with options provided'''

        # Check to make sure not too many experiments requested in plot
        assert len(expts) <= len(MARKERS), \
            "Too many experiments for number of usable symbols"

        # Plot options
        self.plot_options = dict(absolute=absolute, extend_y=extend_y,
                                 no_min_exty=no_min_exty, alt=alt)

        # Create data objects from each input file
        self.data_sub = MetricData(file_sub, desc='Subset')
        self.data_nrm = MetricData(file_nrm, desc=norm)
        self.data_obs = MetricData(file_obs, desc='Observations')
        self.data_std = MetricData(file_std, desc='Model Uncertainty')
        self.data_mod = [MetricData(mfile, desc=expt) for
                         (mfile, expt) in zip(file_mod, expts)]

        # Split obs into obs range and acceptable limits range
        (self.data_obs, self.data_acc) = self.data_obs.split_obs()

        # Filter control from subset file
        # Note that if there is no subset file this does nothing!
        # Also remove zero or missing data entries
        self.data_nrm.filter(self.data_sub, nonmdi=True, nonzero=True)

        # Filter data files based on list of fields left in control
        self.data_obs.filter(self.data_nrm)
        self.data_acc.filter(self.data_nrm)
        self.data_std.filter(self.data_nrm)
        for data_mod in self.data_mod:
            data_mod.filter(self.data_nrm)

        # Check to see if obs data is absolute or relative
        # Modify obs and acc data based on whether absolute plotting or not
        # Should be able to make this obsolete.
        self.check_obs_data()

        # Normalise obs/acc/stdev/model
        self.data_obs.normalise(self.data_nrm)
        self.data_acc.normalise(self.data_nrm)
        self.data_std.normalise(self.data_nrm, scale=1.5)
        for data_mod in self.data_mod:
            data_mod.normalise(self.data_nrm)

        # Normalise control data (=1!) and add model uncertainty for plotting
        self.data_nrm.normalise(self.data_nrm)
        self.data_nrm.addstd(self.data_std)

    def plot(self, ofile=None):
        '''Method to interface with plotting routine'''
        plot_nac(self.data_nrm, self.data_mod, self.data_obs,
                 self.data_acc, self.data_std, self.plot_options,
                 ofile=ofile)

    # Need to make this routine obsolete, should be able to configure
    # obs ranges so that this doesn't matter.
    def check_obs_data(self):
        '''
        Check to see if obs data is absolute or relative
        If absolute then switch plot to absolute
        Modify data if necessary
        '''
        # Check if observational uncertainty is absolute
        if self.data_obs.hasdata:
            if self.data_obs.abssum(0) > 0:
                self.plot_options['absolute'] = True

        # Check if acceptable range is absolute
        if self.data_acc.hasdata:
            if self.data_acc.abssum(0) == 0:
                self.plot_options['absolute'] = True

        # Modify observational uncertainty data if absolute range plotting
        if self.data_obs.hasdata:
            if self.data_obs.abssum(0) == 0:
                if self.plot_options['absolute']:
                    self.data_obs.absmod()

        # Modify acceptable range data based if absolute range plotting
        if self.data_acc.hasdata:
            if self.data_acc.abssum(0) == 0:
                if self.plot_options['absolute']:
                    self.data_acc.absmod()


def get_data_lim(data_exps, data_obs, absolute=False,
                 extend_y=False, no_min_exty=False):
    '''Determine data axis limits'''
    # Set extension factor
    extend = 1.05

    # Calculate max/min for experiments
    minval = min([data_exp.min() for data_exp in data_exps])
    maxval = max([data_exp.max() for data_exp in data_exps])

    # If using absolute plot ranges
    if absolute:
        # Calculate max abs value for experiments
        maxabs = max([data_exp.absmax() for data_exp in data_exps])
        # If want to extend beyond range of observations
        if extend_y:
            maxabs = max(maxabs, data_obs.absmax())
            minval = min(minval, data_obs.min())
            maxval = max(maxval, data_obs.max())
        # Add/Subtract a little bit extra to fully contain plot
        extra = maxabs * (extend - 1.0)
        minval -= extra
        maxval += extra
        # Make range no less than (0, 2) unless otherwise requested
        if not no_min_exty:
            maxval = max(maxval, 2.0)
            minval = min(minval, 0.0)
    # If using relative plot ranges
    else:
        minval = 0.0
        maxval = max(2, maxval * extend)
        # If want to extend beyond range of observations
        if extend_y:
            maxval = max(maxval, data_obs.max()) * extend

    return (minval, maxval)


def plot_std(data, color=STD_GREY, alt=False):
    '''Plot model uncertainty as filled bars about nac=1 line'''
    if data.hasdata:
        coord = np.arange(data.data.size) + 1
        std = data.data[0, :]
        if alt:
            mplt.barh(coord, 2.0*std, left=1.0-std, height=1.0,
                      align='center', color=color, linewidth=0, zorder=0)
        else:
            mplt.bar(coord, 2.0*std, bottom=1.0-std, width=1.0,
                     align='center', color=color, linewidth=0, zorder=0)


def plot_obs(data, color=OBS_GREY, alt=False):
    '''Plot obs range as error bars'''
    if data.hasdata:
        coord = np.arange(data.data.shape[1]) + 1
        orig = 0.5 * (data.data[1, :] + data.data[0, :])
        err = 0.5 * (data.data[1, :] - data.data[0, :])
        # Get user warnings here if any data is missing
        if alt:
            mplt.errorbar(orig, coord, xerr=err,
                          fmt=None, ecolor=color, capsize=5, zorder=1)
        else:
            mplt.errorbar(coord, orig, yerr=err,
                          fmt=None, ecolor=color, capsize=5, zorder=1)


def calc_obs_range(data_obs, data_acc):
    '''Calculate obs range to use in judging model performance'''
    # Set up missing range by default
    obs_range = ma.masked_values(np.ones((2, data_obs.data.shape[1])), 1.0)
    # Use observational uncertainty if available
    if data_obs.hasdata:
        obs_range = data_obs.data
    # Use acceptable range if available
    if data_acc.hasdata:
        obs_range = data_acc.data
    return obs_range


def data_location_check(model, obs_range, upper=None):
    '''Obtain position of model data relative to obs'''
    if upper and upper.hasdata:
        model += upper.data
    diff = model - obs_range
    inside = diff.prod(axis=0) <= 0
    err = np.fabs(diff).min(axis=0)
    return (inside, err)


def metric_colours(nrm, mod, obs_msk):
    '''Routine to determine colour of metric in plot'''
    (nrm_ins, nrm_err) = nrm
    (mod_ins, mod_err) = mod

    # RED is default color
    mcolour = np.array([RED]*nrm_ins.size)

    # GREEN is where model inside obs range color
    mcolour[np.where(mod_ins)] = GREEN

    # AMBER is where neither model or norm are in obs range,
    #  but model is closer than norm
    amber_logic = (~mod_ins) & (~nrm_ins) & (mod_err <= nrm_err)
    mcolour[np.where(amber_logic)] = AMBER

    # GREY is where there are no obs or acceptable ranges
    mcolour[np.where(obs_msk)] = NOOBS_GREY
    return mcolour


def create_legend(data_ctl, data_exps):
    '''
    Create legend by from a collection of fake plotlines with
    appropriate plot styles and labels
    '''
    labels = [data_exp.desc for data_exp in data_exps]
    legend = mplt.legend(FAKELINES[0:len(labels)], labels,
                         bbox_to_anchor=(1, 1), loc=2, numpoints=1,
                         fancybox=True, fontsize='small')
    legend.set_title('Vs %s' % data_ctl.desc, prop={'size': 'small'})
    return legend


def set_plot_niceties(data_ctl, data_exps, data_obs, plot_options):
    '''Routine to add plot niceties'''

    # Determine data axis limits
    limits = get_data_lim(data_exps, data_obs,
                          absolute=plot_options['absolute'],
                          extend_y=plot_options['extend_y'],
                          no_min_exty=plot_options['no_min_exty'])

    # Get current axes instance in order to add plot niceties
    axes = mplt.gca()

    # Set limits, label axes and add norm=1 line
    if plot_options['alt']:
        axes.axvline(1.0, color=BLACK)
        axes.set_yticks(np.arange(data_ctl.data.size)+1)
        axes.set_yticklabels(data_ctl.label, ha='right', fontsize='x-small')
        axes.set_ylim(data_ctl.data.size+0.5, 0.5)
        axes.set_xlim(limits)
        axes.set_xlabel('Normalised Assessment Criteria', fontsize='small')
        axes.tick_params(axis='x', labelsize='small')
    else:
        axes.axhline(1.0, color=BLACK)
        axes.set_xticks(np.arange(data_ctl.data.size)+1)
        axes.set_xticklabels(data_ctl.label, ha='left', fontsize='x-small',
                             rotation=315, stretch='ultra-condensed')
        axes.set_xlim(0.5, data_ctl.data.size+0.5)
        axes.set_ylim(limits)
        axes.set_ylabel('Normalised Assessment Criteria', fontsize='small')
        axes.tick_params(axis='y', labelsize='small')

    # Create legend
    return create_legend(data_ctl, data_exps)


def plot_nac(data_ctl, data_exps, data_obs, data_acc, data_std,
             plot_options, ofile=None):
    '''
    Generate Normalised Assessment Criteria plot
    Write to specified file or display on screen

    Input data must be MetricData and have been made plot ready

    '''
    try:
        assert isinstance(data_ctl, MetricData) and data_ctl.plotready
        assert isinstance(data_obs, MetricData) and data_obs.plotready
        assert isinstance(data_acc, MetricData) and data_acc.plotready
        assert isinstance(data_std, MetricData) and data_std.plotready
        for data_exp in data_exps:
            assert isinstance(data_exp, MetricData) and data_exp.plotready
    except AssertionError, msg:
        raise ValueError("One of the inputs is wrong \n"+str(msg))

    # Calculate obs range to use in judging model performance
    obs_range = calc_obs_range(data_obs, data_acc)

    # Obtain position of normalised control data relative to obs
    # This is normally 1, but if we have model uncertainty data then
    # this becomes 1+data_std
    norm_check = data_location_check(data_ctl.data, obs_range)

    # Create plot figure
    fig = mplt.figure()

    # Plot model uncertainty, obs ranges and acceptable ranges
    plot_std(data_std, color=STD_GREY, alt=plot_options['alt'])
    plot_obs(data_acc, color=ACC_GREY, alt=plot_options['alt'])
    plot_obs(data_obs, color=OBS_GREY, alt=plot_options['alt'])

    # Loop over experiments and MARKERS
    for (data_exp, marker) in zip(data_exps, MARKERS):
        # Obtain position of expt relative to obs
        modl_check = data_location_check(data_exp.data, obs_range)

        # Get metric colours
        mcolor = metric_colours(norm_check, modl_check, obs_range[1, :].mask)

        # Plot data using scatter plot as that is only way of having
        # different colour symbols
        # What does having fixed zorder do for multiple experiments?
        if plot_options['alt']:
            mplt.scatter(data_exp.data,
                         np.arange(data_exp.data.size)+1,
                         s=50, edgecolors=BLACK, c=mcolor,
                         marker=marker, zorder=3)
        else:
            mplt.scatter(np.arange(data_exp.data.size)+1,
                         data_exp.data,
                         s=50, edgecolors=BLACK, c=mcolor,
                         marker=marker, zorder=3)

    # Add plot niceties, including legend
    leg = set_plot_niceties(data_ctl, data_exps, data_obs, plot_options)

    fig.tight_layout()

    if ofile:
        # Note that bbox_inches only works for png plots
        mplt.savefig(ofile, bbox_extra_artists=(leg,), bbox_inches='tight')
    else:
        mplt.show()
    mplt.close()


def central_metrics(run, name, period=False):
    '''Construct name of centrally held metrics file from run dictionary'''
    store = run['mstore']
    area_id = run['area']
    if run['subarea']:
        area_id += '_{}'.format(run['subarea'])
    fname = os.path.join(store, '{0}_{1}'.format(area_id, run[name]))
    if period and run.period:
        fname += '_{}'.format(run.period)
    fname += '.csv'
    return fname


def nacfromarea(run_ctl, run_exps, plotfile, metrics=None,
                alt=False, absolute=False, extend_y=False, no_min_exty=False):
    '''Plot the metrics in standard format from given list of metrics'''
    if metrics:
        assert isinstance(metrics, list)

    # Get run ids
    ctl_id = '{1} ({0})'.format(run_ctl.id, run_ctl.title)
    exps_id = ['{1} ({0})'.format(run.id, run.title) for run in run_exps]

    # Get list of metrics
    file_sub = metrics

    # Use newly created metrics files or use centrally held versions
    file_ctl = '{}.csv'.format(run_ctl['metrics_file'])
    if not os.path.exists(file_ctl):
        file_ctl = central_metrics(run_ctl, 'metrics_file', period=True)

    # Use newly created metrics files or use centrally held versions
    file_exps = []
    for run_exp in run_exps:
        file_exp = '{}.csv'.format(run_exp['metrics_file'])
        if not os.path.exists(file_exp):
            file_exp = central_metrics(run_exp, 'metrics_file', period=True)
        file_exps.append(file_exp)

    # Use centrally held observations metrics files
    file_obs = central_metrics(run_ctl, 'metrics_obs_name', period=True)
    if not os.path.exists(file_obs):
        file_obs = central_metrics(run_ctl, 'metrics_obs_name', period=False)
        if not os.path.exists(file_obs):
            file_obs = None

    # Use centrally held model uncertainty metrics files
    if 'metrics_model_uncertainty' in run_ctl:
        file_std = central_metrics(run_ctl, 'metrics_model_uncertainty',
                                   period=True)
        if not os.path.exists(file_std):
            file_std = central_metrics(run_ctl, 'metrics_model_uncertainty',
                                       period=False)
            if not os.path.exists(file_std):
                file_std = None
    else:
        file_std = None

    # Fill argument tuple with required fields
    nac_args = (ctl_id, exps_id,
                file_ctl, file_exps, file_obs, file_std, file_sub)

    # Fill keywords dictionary with options
    nac_kwargs = dict(alt=alt,
                      absolute=absolute,
                      extend_y=extend_y,
                      no_min_exty=no_min_exty)

    # Produce plot
    try:
        nac = NACPlot(*nac_args, **nac_kwargs)
        nac.plot(ofile=plotfile)
    except:
        raise


def nacfromcl():
    '''Plot the metrics in standard format'''
    import argparse

    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Create Normalised Assessment Criteria Plots'
    )

    # Create arguments
    parser.add_argument('--plot', dest='plot', default=None,
                        help='Plot file to be created')
    parser.add_argument('--exp', dest='exp', default=None,
                        help='Experiments for which validation notes and \
assessment are to be done (commma separated)')
    parser.add_argument('--ctl', dest='ctl', default=None,
                        help='Experiment which --exp options are compared \
against')
    parser.add_argument('--file_obs', dest='file_obs', default=None,
                        help='Observations metrics file')
    parser.add_argument('--file_sub', dest='file_sub', default=None,
                        help='Metrics subset file')
    parser.add_argument('--file_std', dest='file_std', default=None,
                        help='Model uncertainty metrics file')
    parser.add_argument('--file_ctl', dest='file_ctl', default=None,
                        help='Control metrics file')
    parser.add_argument('--file_exp', dest='file_exp', default=None,
                        help='Experiment metrics files (commma separated)')
    parser.add_argument('--abs', dest='abs', default=False,
                        action='store_true',
                        help='Plot absolute metric ranges, not relative')
    parser.add_argument('--exty', dest='exty', default=False,
                        action='store_true',
                        help='Extend y axis in plots for obs uncertainties \
if required')
    parser.add_argument('--nominy', dest='nominy', default=False,
                        action='store_true',
                        help='Dont limit y axis to (0,2) if range is less \
than that')
    parser.add_argument('--alt', dest='alt', default=False,
                        action='store_true',
                        help='Alternative format plot')

    # Parse arguments
    options = parser.parse_args()

    # Check size of experiment inputs
    file_model = options.file_exp.split(',')
    expt_model = options.exp.split(',')
    if len(file_model) != len(expt_model):
        raise ValueError('Number of experiments and number of experiment \
files are different, this will cause plotting to fail')

    # Fill argument tuple with required fields
    nac_args = (options.ctl, expt_model, options.file_ctl, file_model,
                options.file_obs, options.file_std, options.file_sub)

    # Fill keywords dictionary with options
    nac_kwargs = dict(alt=options.alt,
                      absolute=options.abs,
                      extend_y=options.exty,
                      no_min_exty=options.nominy)

    # Set plot file name
    plotfile = options.plot

    # Produce plot
    try:
        nac = NACPlot(*nac_args, **nac_kwargs)
        nac.plot(ofile=plotfile)
    except:
        raise


# Capability to run plotting script from command line, use:
# > python2.7 plot_norm_ac.py <args>
if __name__ == '__main__':
    nacfromcl()
