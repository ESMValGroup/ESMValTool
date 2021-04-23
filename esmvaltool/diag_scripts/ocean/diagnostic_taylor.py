"""
Taylor diagram code adapted from https://gist.github.com/ycopin/3342888
Acknowledgment to Yannick Copin
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.projections import PolarAxes
import mpl_toolkits.axisartist.floating_axes as fa
import mpl_toolkits.axisartist.grid_finder as gf


class TaylorDiagram():
    """
    Taylor diagram.

    Plot model standard deviation and correlation to reference (data)
    sample in a single-quadrant polar plot, with r=stddev and
    theta=arccos(correlation).
    """

    def __init__(self, refstd,
                 fig=None, rect=111, label='_', srange=(0, 1.5), extend=False):
        """
        Set up Taylor diagram axes, i.e. single quadrant polar
        plot, using `mpl_toolkits.axisartist.floating_axes`.

        Parameters

        refstd: float
            reference standard deviation to be compared to
        fig: matplotlib figure object
            input Figure or None
        rect: integer
            subplot definition
        label: text
            reference label style
        srange: tuple or list
            stddev axis extension, in units of refstd
        extend: boolean
            extend diagram to negative correlations
        """

        self.refstd = refstd            # Reference standard deviation

        trpolar = PolarAxes.PolarTransform()

        # Correlation labels
        rlocs = np.array([0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1])
        if extend:
            # Diagram extended to negative correlations
            self.tmax = np.pi
            rlocs = np.concatenate((-rlocs[:0:-1], rlocs))
        else:
            # Diagram limited to positive correlations
            self.tmax = np.pi/2
        tlocs = np.arccos(rlocs)        # Conversion to polar angles
        gl1 = gf.FixedLocator(tlocs)    # Positions
        tf1 = gf.DictFormatter(dict(zip(tlocs, map(str, rlocs))))

        # Standard deviation axis extent (in units of reference stddev)
        self.smin = srange[0] * self.refstd
        self.smax = srange[1] * self.refstd

        ghelper = fa.GridHelperCurveLinear(
            trpolar,
            extremes=(0, self.tmax, self.smin, self.smax),
            grid_locator1=gl1, tick_formatter1=tf1)

        if fig is None:
            fig = plt.figure()

        axs = fa.FloatingSubplot(fig, rect, grid_helper=ghelper)
        fig.add_subplot(axs)

        # Adjust axes
        axs.axis["top"].set_axis_direction("bottom")   # "Angle axis"
        axs.axis["top"].toggle(ticklabels=True, label=True)
        axs.axis["top"].major_ticklabels.set_axis_direction("top")
        axs.axis["top"].label.set_axis_direction("top")
        axs.axis["top"].label.set_text("Correlation")

        axs.axis["left"].set_axis_direction("bottom")  # "X axis"
        axs.axis["left"].label.set_text("Standard deviation")

        axs.axis["right"].set_axis_direction("top")    # "Y-axis"
        axs.axis["right"].toggle(ticklabels=True)
        axs.axis["right"].major_ticklabels.set_axis_direction(
            "bottom" if extend else "left")

        if self.smin:
            axs.axis["bottom"].toggle(ticklabels=False, label=False)
        else:
            axs.axis["bottom"].set_visible(False)          # Unused

        self._ax = axs                         # Graphical axes
        self.axs = axs.get_aux_axes(trpolar)   # Polar coordinates

        # Add reference point and stddev contour
        val, = self.axs.plot([0], self.refstd, 'ko',
                             ls='', ms=8, label=label)
        tlocs = np.linspace(0, self.tmax)
        rlocs = np.zeros_like(tlocs) + self.refstd
        self.axs.plot(tlocs, rlocs, 'k--', label='_')

        # Collect sample points for latter use (e.g. legend)
        self.samplepoints = [val]


    def add_sample(self, stddev, corrcoef, *args, **kwargs):
        """
        Add sample (*stddev*, *corrcoeff*) to the Taylor
        diagram. *args* and *kwargs* are directly propagated to the
        `Figure.plot` command.
        """
        val, = self.axs.plot(np.arccos(corrcoef), stddev,
                             *args, **kwargs)  # (theta, radius)
        self.samplepoints.append(val)

        return val


    def add_grid(self, *args, **kwargs):
        """Add a grid."""
        self._ax.grid(*args, **kwargs)


    def add_contours(self, levels=5, **kwargs):
        """
        Add constant centered RMS difference contours, defined by *levels*.
        """
        rps, tps = np.meshgrid(np.linspace(self.smin, self.smax),
                               np.linspace(0, self.tmax))
        # Compute centered RMS difference
        rms = np.sqrt(self.refstd**2 + rps**2 - 2*self.refstd*rps*np.cos(tps))
        contours = self.axs.contour(tps, rps, rms, levels, **kwargs)

        return contours
