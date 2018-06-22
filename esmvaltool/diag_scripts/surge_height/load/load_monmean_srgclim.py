def load_monmean_srgclim(stat):
	from netCDF4 import Dataset
	PATHin = 'data/srgclim/'
	filename = 'monanom_ERAintWAQUA_surge_1979-2016_speed_'
	if not type(stat) == list:
		stat = [stat]
	#
	monmean_srgclim = {}
	for ifname in stat:
		nc = Dataset(PATHin + filename + str(ifname) + '.nc','r')
		monmean_srgclim[str(ifname)] = nc.variables['WAQUA_surge'][:]
		nc.close()
	#
	return monmean_srgclim
