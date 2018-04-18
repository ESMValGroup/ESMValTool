'''
Module to analyse the wallclock time files
'''

import copy
import csv
import os
import pdb

import globalvar

import vnc


class csvError(Exception):
    """ Error class to handle problems loading the csv file. """
    def __init__(self, message):
        print message


class SpeedInfo:
    """
    A class that contains all the key information about
    a job required for speed calculations.

    Contains:
    AverageWallclock = The average wallclock time (minutes)
    TotalCPU = The total cpu time (minutes)
    CPUprocessor = The cpu time per processor (minutes)
    MaxWallclock = The maximum wallclock time (minutes)
    timesteps = The number of timesteps in a resubmission period
    ResubDays = The number of days in a resubmission period
    SpeedYpDay = The speed of the model in years per day
    pes = The number of processors
    wWind = The maximum w component of wind
    """

    def __init__(self, csvfile):
        """ What gets initialised with this class. """

        # Initialise lists
        pes_l = []
        timesteps_l = []
        resub_mon = []
        resub_day = []
        wallclock = []
        cpu = []
        w = []

        # Load the csv file
        csvreader = csv.reader(open(csvfile, 'rb'))

        # Make a list of timesteps and processors
        self.timesteps = 0
        for row in csvreader:
            if len(row) >= 5:
                timesteps_l.append(int(row[3]))
                pes_l.append(int(row[0]))

        # Find the most used number of timesteps and the most often
        #  used number of processors
        self.timesteps = max_count(timesteps_l)
        self.pes = max_count(pes_l)

        # Populate the lists
        csvreader = csv.reader(open(csvfile, 'rb'))
        for row in csvreader:
            if len(row) >= 5:
                if int(row[3]) == self.timesteps and int(row[0]) == self.pes:
                    resub_mon.append(int(row[1]))
                    resub_day.append(int(row[2]))
                    try:
                        wallclock.append(float(row[4]))
                    except ValueError:
                        print 'Error converting wallclock time to float'
                    try:
                        cpu.append(float(row[5]))
                    except:
                        pass
                    try:
                        w.append(float(row[6]))
                    except:
                        pass

        # Check that we have some values
        if len(wallclock) == 0:
            raise csvError('Unable to gather any wallclock times from file '+csvfile)

        # Check that the resubmission period stayed the same
        if max(resub_mon) != min(resub_mon) and \
           max(resub_day) != min(resub_day):
            raise csvError('Resubmission periods vary throughout the run in file '+csvfile)

        # Generate output
        self.AverageWallclock = mean(wallclock)/60.0
        self.MaxWallclock = max(wallclock)/60.0
        self.ResubDays = resub_mon[0]*30 + resub_day[0]
        self.SpeedYpDay = (60.0*24.0)/(self.AverageWallclock*360.0/self.ResubDays)
        self.TotalCPU = mean(cpu)/60.0
        self.DayCPU = self.TotalCPU/self.ResubDays
        self.maxw = globalvar.missing_data
        if len(w) > 0:
            self.maxw = max(w)
        self.meanw = mean(w)


def mean(inlist):
    """ Calculate the mean of a list """
    thismean = globalvar.missing_data
    if len(inlist) > 0:
        thismean = sum(inlist)/len(inlist)
    return thismean


def max_count(inlist):
    """ Calculate the most frequent occurance of a number in a list """
    count = 0
    for item in set(inlist):
        thiscount = len([x for x in inlist if x == item])
        if thiscount >= count:
            count = thiscount
            returnnum = item
    return returnnum


def get_table_heading(jobid):
    """ Get the title of the jobid from source_file.dat """
    heading = ' '
    csvreader = csv.reader(open('source_file.dat', 'rb'), delimiter=':')
    for row in csvreader:
        try:
            if row[0] == jobid+'.name':
                heading = row[1]
        except IndexError:
            pass
    return heading


def add_table_heading(page, controlcol, expercol):
    """ Add a heading to the html table."""
    page.tr.open()
    page.th()
    page.th(controlcol)
    page.th(expercol)
    page.tr.close()


def add_table_line(page, title, controlcol, expercol, format):
    """ Add a line to the html table with information."""
    page.tr.open()
    page.td(title)
    page.td(format % controlcol)
    page.td(format % expercol)
    page.tr.close()


def create_page(page_top, exper, control):
    """
    Create a web page documenting the speed of your
    model compared to the control.

    Inputs:
    exper = string containing the experiment jobid
    control = string containing the control jobid
    valnote_dir = where the validation note code is
    """

    # Set up defaults
    user = os.environ['USER']
    tmpdir = os.environ['TMPDIR']

    # Set up lists (_l), dictionaries (_d)
    speedinfo_d = {}

    # Loop over jobs
    for job in [exper, control]:

        # Get the wallclock files
        os.system('scp -B '+user+'@hpc2e:/data/cr/cmde/hadco/wallclock/'+job+'_wallclock.list '+tmpdir)

        # Get the iteration count files
        os.system('scp -B '+user+'@hpc2e:/data/cr/cmde/hadco/iterations/'+job+'_iterations.list '+tmpdir)

        # Extract the data
        try:
            speedinfo_d[job] = SpeedInfo(tmpdir+'/'+job+'_wallclock.list')
        except csvError:
            print 'CSV reading error. Not creating speed page.'
            return
        except IOError:
            print 'Missing speed file. Not creating speed page.'
            return
        except ZeroDivisionError:
            print 'Divide by zero error. This may be due to errors in the wallclock csv file. Not creating speed page.'
            return

    # Make a plot of the iteration counts
    # os.system('tidl '+valnote_dir+'/idl/iterations.pro -args '+tmpdir+'/'+control+'_iterations.list '+tmpdir+'/'+exper+'_iterations.list')

    page = copy.deepcopy(page_top)
    page.br()
    page.div(style='text-align:center')
    page.table.open(width=800)

    # Get the table headings
    control_head = get_table_heading(control)
    exper_head = get_table_heading(exper)

    # Populate the table headings
    add_table_heading(page, control_head, exper_head)
    add_table_heading(page, control, exper)

    # Populate the table
    add_table_line(page, 'Number of processors', speedinfo_d[control].pes, speedinfo_d[exper].pes, '%i')
    add_table_line(page, 'Length of resubmission period (days)', speedinfo_d[control].ResubDays, speedinfo_d[exper].ResubDays, '%i')
    add_table_line(page, 'Number of timesteps in resubmission period', speedinfo_d[control].timesteps, speedinfo_d[exper].timesteps, '%i')
    if speedinfo_d[control].maxw > 0.0:
        add_table_line(page, 'Maximum w wind (over whole run)', speedinfo_d[control].maxw, speedinfo_d[exper].maxw, '%.3f')
        add_table_line(page, 'Average maximum w wind (average of w_max in each resubmission)', speedinfo_d[control].meanw, speedinfo_d[exper].meanw, '%.3f')
    add_table_line(page, 'Maximum wallclock time for resubmission period (minutes)', speedinfo_d[control].MaxWallclock, speedinfo_d[exper].MaxWallclock, '%.3f')
    add_table_line(page, 'Average wallclock time for resubmission period (minutes)', speedinfo_d[control].AverageWallclock, speedinfo_d[exper].AverageWallclock, '%.3f')
    if speedinfo_d[control].TotalCPU > 0.0:
        add_table_line(page, 'CPU time for resubmission period (minutes)', speedinfo_d[control].TotalCPU, speedinfo_d[exper].TotalCPU, '%.3f')
        add_table_line(page, 'CPU time for a day (minutes)', speedinfo_d[control].DayCPU, speedinfo_d[exper].DayCPU, '%.3f')
    add_table_line(page, 'Speed of the model if continuous (years/day)', speedinfo_d[control].SpeedYpDay, speedinfo_d[exper].SpeedYpDay, '%.3f')

    # Close the table
    page.table.close()

    # What is the percentage change in speed of a processor
    page.br()
    page.h1.open()
    if speedinfo_d[control].DayCPU > 0.0:
        speedup = speedinfo_d[control].DayCPU / speedinfo_d[exper].DayCPU * 100.0 - 100.0
        explanation = 'based on cpu time'
        if speedup > 0.1:
            page.add('%s is %.1f%% faster than %s (%s)' %
                     (exper, speedup, control, explanation))
        elif speedup < -0.1:
            page.add('%s is %.1f%% slower than %s (%s)' %
                     (exper, -1.0*speedup, control, explanation))
        else:
            page.add('%s is the same speed as %s (%s)' %
                     (exper, control, explanation))

    # What is the percentage change in speed of the model
    speedup = speedinfo_d[exper].SpeedYpDay / speedinfo_d[control].SpeedYpDay * 100.0 - 100.0
    explanation = 'based on wallclock time'
    page.br()
    if speedup > 0.1:
        page.add('%s is %.1f%% faster than %s (%s)' %
                 (exper, speedup, control, explanation))
    elif speedup < -0.1:
        page.add('%s is %.1f%% slower than %s (%s)' %
                 (exper, -1.0*speedup, control, explanation))
    else:
        page.add('%s is the same speed as %s (%s)' %
                 (exper, control, explanation))

    page.h1.close()
    page.br()
    page.br()

    # Plot the iterations plot
    if os.path.isfile('iterations.ps'):
        os.system('convert -rotate 270 iterations.ps iterations.png')
        page.br()
        page.h1.open()
        page.add('Histogram plot of number of iterations of the solver - ' +
                 'plotted as a fraction of time')
        page.br()
        page.img(src='iterations.png')
        page.h1.close()
        page.br()
        page.br()

    # Output to file
    with open('speed.html', 'w') as file:
        print >>file, page
