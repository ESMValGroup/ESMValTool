"""
Implementation for DUMMY diagnostics into ESMValTool
"""

import iris
from c3s_511_basic import Basic_Diagnostic
from ESMValMD import ESMValMD
import sys
import os
import matplotlib.pyplot as plt
[sys.path.insert(0, os.path.join(os.path.dirname(
    os.path.abspath(__file__)), dir)) for dir in ["lib", "plots"]]
from plot import Plot2D


class DUMMY_Diagnostic(Basic_Diagnostic):
    """
    class to implement dummy diagnostics, like e.g. global means,
    global differences, RMSD etc.
    """

    def run_diagnostic(self):
        """
        run parent diagnostic and the DUMMY specific diagnostic
        """
        super(DUMMY_Diagnostic, self).run_diagnostic()
        self.__my_diag__()

    def __my_diag__(self):
        """
        method
        """
        # the name of the diagnostic
        this_function = "my diag"

        # takes the first slice of a 3D cube
        my_data = self.sp_data[1, :, :]

        # make list of plots
        list_of_plots = []

        # define filename
        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + \
            "_mydiag" + "." + self.__output_type__

        # add current plot to list
        list_of_plots.append(filename)

        # produce the plot object
        x = Plot2D(my_data)

        # initialize the figure
        fig = plt.figure()

        # intialize the axis where we plot
        ax = [plt.subplot(1, 1, 1)]

        # call the plot
        x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [
               0, 2, 1, 3]]) + " (" + self.__time_period__ + ")")

        # save figure
        fig.savefig(filename)

        # close figure (memory reasons)
        plt.close(fig.number)

        # add metadata
        ESMValMD("meta",
                 filename,
                 self.__basetags__ + ['DM_global', 'C3S_mean_var'],
                 # caption below EDIT!
                 str("mydiag" +
                     ' values of ' +
                     self.__varname__ +
                     ' for the data set "' +
                     "_".join(self.__dataset_id__) +
                     '" (' +
                     self.__time_period__ +
                     ')'),
                 # identifier EDIT!
                 '#C3S' + "mydiag" + self.__varname__,
                 self.__infile__,
                 self.diagname,
                 self.authors)

        # produce report
        self.__do_report__(
            content={
                "plots": list_of_plots},
            filename=this_function.upper())
