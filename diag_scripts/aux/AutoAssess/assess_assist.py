import os, sys, stat, subprocess, glob, config
import restore_data
from assess_obj import Lock, Cmd
from mksuper_mod import mksuper_main

def parse_config(options):
    conf = config.Config()

    if not os.path.exists(options.config_file):
        sys.stderr.write('Config file %s does not exist\n' % options.config_file)
        sys.exit(2)

    conf.read(options.config_file)

    try:
        options.valorder = conf.get('valnote', 'valorder')

        if not os.path.exists(options.valorder):
            sys.stderr.write('%s: Valnote file %s does not exist\n' % (options.config_file,options.valorder))
            sys.exit(2)

    except config.ConfigParser.NoOptionError, e:
        options.valorder = None


    try:
        options.valcolour = conf.get('valnote', 'valcolour')

        if not os.path.exists(options.valcolour):
            sys.stderr.write('%s: Valcolour file %s does not exist\n' % (options.config_file,options.valcolour))
            sys.exit(2)

    except config.ConfigParser.NoOptionError, e:
        options.valcolour = None


    try:
        options.itemfile = conf.get('valnote', 'itemfile')

        if not os.path.exists(options.itemfile):
            sys.stderr.write('%s: Itemfile file %s does not exist\n' % (options.config_file,options.itemfile))
            sys.exit(2)

    except config.ConfigParser.NoOptionError, e:
        options.itemfile = None


    try:
        datadir = conf.get('directories', 'datadir')

        # Check directory exists
        if not os.path.exists(datadir):
            sys.stderr.write('datadir %s does not exist\n' % datadir)
            sys.exit(2)

        options.datadir = datadir

    except config.ConfigParser.NoOptionError, e:
        sys.stderr.write('Unable to obtain datadir from config file %s\n' % options.config_file)
        sys.exit(2)


    try:
        auto_assess_results_dir = conf.get('directories', 'auto_assess_results_dir')

        # Check directory exists
        if not os.path.exists(auto_assess_results_dir):
            sys.stderr.write('auto_assess_results_dir %s does not exist\n' % auto_assess_results_dir)
            sys.exit(2)

        auto_assess_results_dir = os.path.join(auto_assess_results_dir, '%s_v_%s' % (options.expt, options.cntl))
        options.auto_assess_results_dir = os.path.abspath(auto_assess_results_dir)

    except config.ConfigParser.NoOptionError, e:
        sys.stderr.write('Unable to obtain auto_assess_results_dir from config file %s\n' % options.config_file)
        sys.exit(2)


    try:
        valnote_results_dir = conf.get('directories', 'valnote_results_dir')

        # Check directory exists
        if not os.path.exists(valnote_results_dir):
            sys.stderr.write('valnote_results_dir %s does not exist\n' % valnote_results_dir)
            sys.exit(2)

        valnote_results_dir = os.path.join(valnote_results_dir, '%s_v_%s' % (options.expt, options.cntl))
        options.valnote_results_dir = os.path.abspath(valnote_results_dir)

    except config.ConfigParser.NoOptionError, e:
        sys.stderr.write('Unable to obtain valnote_results_dir from config file %s\n' % options.config_file)
        sys.exit(2)


    try:
        options.flag_dir = conf.get('directories', 'flags')

    except config.ConfigParser.NoOptionError, e:
        sys.stderr.write('Unable to obtain flags from config file %s\n' % options.config_file)
        sys.exit(2)


    try:
        options.metrics_dir = conf.get('directories', 'metrics')

    except config.ConfigParser.NoOptionError, e:
        sys.stderr.write('Unable to obtain metrics from config file %s\n' % options.config_file)
        sys.exit(2)


    try:
        options.plots_dir = conf.get('directories', 'plots')

    except config.ConfigParser.NoOptionError, e:
        sys.stderr.write('Unable to obtain plots from config file %s\n' % options.config_file)
        sys.exit(2)


    try:
        os.mkdir(options.flag_dir)

    except OSError, e:
        if e.errno != 17:  # Directory exists
            sys.stderr.write('Cannot create flag directory %s: %s\n' % (options.flag_dir, str(e)) )
            sys.exit(2)


def def_meanfile(jobid, startyear, nyears, name_system, filestr):
   """ Generate a mean file name """

   # startyear = the start year
   # nyears = the number of years
   # name_system = the naming system to use
   # filestr: m = mean, v = variance, s = standard deviation

   stream=['s', 's', 's', 's', 'y']

   # Set the strings that are needed to generate file names
   yearabs=startyear+nyears-1
   decade=yearabs/10
   cyclic_decade=decade % 36
   if cyclic_decade > 9:
      decadestr=chr(cyclic_decade-10+97)
   else:
      decadestr=chr(cyclic_decade+48)
   yearindec=yearabs-decade*10
   yearindecstr=chr(yearindec+48)

   # When archiving what stream do you want to archive to?
   if nyears <= 9: meanstream=chr(nyears+48)
   if nyears == 10: meanstream='a'
   if nyears >= 11 and nyears <= 19: meanstream='b'
   if nyears >= 20 and nyears <= 24: meanstream='k'
   if nyears >= 25 and nyears <= 29: meanstream='s'
   if nyears >= 30 and nyears <= 39: meanstream='t'
   if nyears >= 40 and nyears <= 49: meanstream='q'
   if nyears >= 50 and nyears <= 99: meanstream='l'
   if nyears >= 100 and nyears <= 249: meanstream='u'
   if nyears >= 250 and nyears <= 499: meanstream='w'
   if nyears >= 500: meanstream='d'

   # Define the filename of the resulting mean, var and stdev files.
   if name_system == 1:
      yearstr = decadestr+yearindecstr
   elif name_system == 2:
      yearstr = str(yearabs)
   else:
      print 'Naming system not defined'

   filename_root = '%sa.%s%s%s' % (jobid, filestr, meanstream, yearstr)
   return filename_root


def setup_runs(options, areas, ctldir):
    ssexpr_expt = {}
    ssexpr_cntl = {}
    for area in areas:
        for filter_file in glob.glob(os.path.join(ctldir, '%s.*.filter'%area)):
            filterFilename = restore_data.FilterFilename(filter_file)
            f_expt = '%s/%s/%s.splitlev.%s' % (options.datadir, options.expt, options.expt, restore_data.hash_dirs[filterFilename.period()])
            f_cntl = '%s/%s/%s.splitlev.%s' % (options.datadir, options.cntl, options.cntl, restore_data.hash_dirs[filterFilename.period()])
            key = 'ss_%s' % filterFilename.period()
            if not ssexpr_expt.has_key(key):
                ssexpr_expt[key] = f_expt
            if not ssexpr_cntl.has_key(key):
                ssexpr_cntl[key] = f_cntl

    if ssexpr_expt.has_key('ss_monthly'):
        ssexpr_expt['ss_spatiot'] = ssexpr_expt['ss_monthly']
    if ssexpr_cntl.has_key('ss_monthly'):
        ssexpr_cntl['ss_spatiot'] = ssexpr_cntl['ss_monthly']

    keys = ['ss_annual', 'ss_seasonal', 'ss_monthly', 'ss_spatiot', 'ss_daily']

    config_parser = config.ConfigParser.ConfigParser()
    
    config_parser.add_section(options.expt)
    config_parser.set(options.expt, 'start', options.start_year)
    config_parser.set(options.expt, 'nyear', options.num_years)
    config_parser.set(options.expt, 'run_type', options.run_type)
    config_parser.set(options.expt, 'data_root', os.path.join(options.datadir, options.expt))

    use_keys = keys[:]
    for key, value in ssexpr_expt.iteritems():
        config_parser.set(options.expt, key, value)
        use_keys.remove(key)

    for key in use_keys:
        config_parser.set(options.expt, key, '')

    config_parser.set(options.expt, 'validation_note', options.valnote_results_dir)

    meanfile = def_meanfile(options.expt, options.start_year, options.num_years, 2, 'm')
    config_parser.set(options.expt, 'supermean_root',  os.path.join(os.path.join(os.path.join(options.datadir, options.expt), 'supermeans'), meanfile))


    #-----------------------------------------------------
    config_parser.add_section(options.cntl)
    config_parser.set(options.cntl, 'start', options.start_year)
    config_parser.set(options.cntl, 'nyear', options.num_years)

    config_parser.set(options.cntl, 'data_root', os.path.join(options.datadir, options.cntl))
    config_parser.set(options.cntl, 'run_type', options.run_type)

    config_parser.set(options.cntl, 'validation_note', options.valnote_results_dir)
    meanfile = def_meanfile(options.cntl, options.start_year, options.num_years, 2, 'm')
    config_parser.set(options.cntl, 'supermean_root',  os.path.join(os.path.join(os.path.join(options.datadir, options.cntl), 'supermeans'), meanfile))

    for key, value in ssexpr_cntl.iteritems():
        config_parser.set(options.cntl, key, value)

    use_keys = keys[:]
    for key, value in ssexpr_cntl.iteritems():
        config_parser.set(options.cntl, key, value)
        use_keys.remove(key)

    for key in use_keys:
        config_parser.set(options.cntl, key, '')


    with open('%s/runs.cfg'%options.tmpdir, 'wb') as configfile:
        config_parser.write(configfile)

    return configfile.name


def make_super(prefix, expt, start, nyear, retrieval_dir, use_query):

    print 'Running mksuper for %s' % expt
    
    # Setup defaults for whether or not we use a query file
    if use_query:
        query_file=os.path.join(prefix, 'maverick/etc/restore_data_filters/valnote_all.query')
    else:
        query_file='None'

    try:
        oldmask=os.umask(0)
        os.makedirs(retrieval_dir, stat.S_IRWXU|stat.S_IRWXG|stat.S_IRWXO)

    except OSError, e:
        if e.errno != 17:  # directory exists
            raise e

    finally:
        os.umask(oldmask)
    cwd = os.getcwd()
    os.chdir(retrieval_dir)
    
    # Create lock for data retrievals
    datalock = Lock(retrieval_dir, 'restore_supermean')
    user = datalock.isset()
    if user is not None:
        sys.stderr.write('There is a lock file in the data directory, cannot retrieve at this time\n')
        sys.stderr.write('File: %s\n' % datalock.filename())
        sys.stderr.write('User: %s\n' % user)
        sys.stderr.write('If you need to remove the lock file you can, but please check with user before hand!\n')
        sys.exit(2)
        
    # Put lock in place to prevent any one else from creating supermeans
    datalock.set()
    
    map(os.unlink, glob.glob('ftn*'))
    mksuper_main(expt, nyear, start, query_file=query_file)
    
    os.chdir(cwd)  
        
    try:
        cmd = Cmd()
        cmd.execute('chmod', '-R', 'a+w', retrieval_dir)
    except Exception, e:
        sys.stderr.write('Failed to give %s a+w permissions\n' % retrieval_dir)
        sys.stderr.write(str(e))
        
    # Remove lock
    datalock.clear()

