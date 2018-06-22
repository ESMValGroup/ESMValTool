def plot_tseries(dates, srg, stat):
	# 
	import matplotlib.pyplot as plt
	import numpy as np
	from datetime import datetime
	#from ..load import load_config as llc
	from load import load_config as llc
	#
	fig, axes = plt.subplots(1,1,figsize=(10,5))
	ax1 = plt.subplot(1,1,1)
	#
	ax1.plot(dates, srg)
	ax1.set_ylabel('surge (m)', fontsize=12)
	ax1.set_xlabel('dates', fontsize=12)
	ax1.set_title(stat)
	ax1.grid()
	#
	date_start = datetime.strftime(dates[0],'%Y-%m-%d')
	date_end   = datetime.strftime(dates[-1],'%Y-%m-%d')
	plt.savefig(llc.savepath + llc.plot_name + '_at_' + stat + '_' + date_start + '_to_' + date_end + 
			'.pdf', dpi=100, format = 'pdf')
