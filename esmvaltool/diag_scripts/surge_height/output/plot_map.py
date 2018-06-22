def plot_map(dates, srg):
	import matplotlib.pyplot as plt
	#from ..load import load_config as llc
	from load import load_config as llc
	from mpl_toolkits.basemap import Basemap
	import numpy as np
	from datetime import datetime
	# coordinates of stations
	coords = {'aberdeen': [81, 111], 'aukalpha': [114, 102], 'bg2': [126, 46], 'borkums': [150, 68], 
		   'bremerha': [165, 68], 'cadzand': [123, 42], 'cromer': [108, 60], 'cuxhaven': [167, 72], 
		   'delfzijl': [153, 65], 'denhelde': [138, 65], 'denoever': [137, 61], 'devonpor': [64, 29], 
		   'dover': [108, 39], 'duinkerk': [116, 38], 'ekofisk': [123, 104], 'esbjerg': [163, 91], 
		   'europlat': [123, 49], 'f3': [134, 83], 'felixsto': [108, 48], 'goeree': [126, 48],  
		   'harlinge': [140, 64], 'helgeroa': [175, 133], 'helgolan': [160, 75], 'hoekvanh': [130, 49], 
		   'holyhead': [60, 66], 'huibertg': [148, 68], 'husum': [168, 78], 'ijmuiden': [133, 54], 
		   'ilfracom': [64, 40], 'immingha': [98, 69], 'innerdow': [101, 65], 'k13a': [123, 64], 
		   'kornwerd': [138, 61], 'lauwerso': [147, 66], 'leith': [72, 97], 'lerwick': [89, 151], 
		   'lowestof': [111, 55], 'meetpost': [131, 53], 'newhaven': [96, 34], 'newlyn': [53, 26],  
		   'northcor': [106, 160], 'northshi': [86, 85], 'offharwi': [109, 47], 'oostende': [120, 40], 
		   'os11': [125, 45], 'os15': [125, 44], 'oscarsbo': [182, 139], 'portsmou': [88, 35], 
		   'roompotb': [126, 45], 'scarboro': [94, 77], 'scheveni': [131, 51], 'scillyis': [45, 23], 
		   'sheernes': [104, 42], 'southend': [103, 43], 'stavange': [141, 133], 'stmarys': [47, 24], 
		   'stornowa': [51, 122], 'terschel': [139, 66], 'texelnoo': [135, 63], 'torsmind': [161, 101], 
		   'tregde': [158, 121], 'vidaa': [166, 85], 'vlaktevd': [122, 43], 'vlissing': [125, 42], 
		   'westkape': [124, 43],'westters': [138, 65], 'weymouth': [77, 32], 'wick': [73, 126],  
		   'zeebrugg': [122, 41]}	
	# create map
	m = Basemap(projection=llc.project, llcrnrlon=-7,
		urcrnrlon=11.5, 
		llcrnrlat=49.4,
		urcrnrlat=61.75,
		resolution=llc.res)
	lons = np.arange(-12.00,-12.00+(0.125*201),0.125).tolist()
	lats = np.arange(48.00,48.00+(0.08333*173),0.08333).tolist()
	x, y = m(lons,lats)
	X, Y = np.meshgrid(lons,lats)

	fig, axes = plt.subplots(nrows=1, ncols=1)
	ax1   = plt.subplot(1,1,1)
	coast = m.drawcoastlines(linewidth=1.)
	bndr  = m.drawmapboundary()
	par   = m.drawparallels(np.arange(47.5,67.5,2.5),color='grey',labels=[1,0,0,0],fontsize=12)
	mer   = m.drawmeridians(np.arange(-16,17.,4.),color='grey',labels=[0,0,0,1],fontsize=12)
	#
	for stat in srg.keys():
		scat = plt.scatter(lons[coords[stat][0]], lats[coords[stat][1]], c = srg[stat], edgecolors='k', 
				cmap = plt.get_cmap(llc.mapcol,20), vmin=-3, vmax=3,zorder=50)
	#
	cbar_scat = m.colorbar(scat, location="bottom",size = "5%", pad="7.5%")
	cbar_scat.set_label('surge height (m)')
	plt.title('North Sea coastal surge ' + datetime.strftime(dates[0],'(%d-%m-%Y)'))
	#plt.show()
	fdates = datetime.strftime(dates[0],'%Y-%m-%d')
	plt.savefig(llc.savepath + llc.plot_name + '_map_' + fdates + '.pdf', dpi=100, format = 'pdf')

