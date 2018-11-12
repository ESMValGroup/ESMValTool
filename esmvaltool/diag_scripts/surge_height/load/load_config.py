def load_config():
	import ConfigParser
	from datetime import datetime
	global config
	config = ConfigParser.ConfigParser()
	config.read('/usr/people/ridder/Documents/0_models/ESMValTool/esmvaltool/diag_scripts/surge_height/config/cfg_srg_estim.conf')

	global coastal_map, plt_tseries, SOIname, t0, tstart, tend, lat_min, lat_max, lon_min, lon_max 
	sect = 'GENERAL'
	coastal_map = config.getboolean(sect,'coastal_map')
	t0          = datetime.strptime(config.get(sect,'t0'),'%Y-%m-%d')
	plt_tseries = config.getboolean(sect,'plt_tseries')
	SOIname     = config.get(sect,'SOIname')
	tstart      = datetime.strptime(config.get(sect,'tstart'),'%Y-%m-%d')
	tend        = datetime.strptime(config.get(sect,'tend'),'%Y-%m-%d')
	lat_min     = float(config.get(sect,'lat_min'))
	lat_max     = float(config.get(sect,'lat_max'))
	lon_min     = float(config.get(sect,'lon_min'))
	lon_max     = float(config.get(sect,'lon_max'))
	
	global project, res, mapcol, savepath, plot_name
	sect = 'PLOT_PARAMS'
	project     = config.get(sect,'projection')
	res         = config.get(sect,'resolution')
	mapcol      = config.get(sect,'colormap')
	savepath    = config.get(sect,'savepath')
	plot_name   = config.get(sect,'plot_name')
