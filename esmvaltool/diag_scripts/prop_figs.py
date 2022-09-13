"""Python example diagnostic."""
import logging
import os

import matplotlib.pyplot as plt
import numpy as np
import iris
from iris.analysis import MEAN
from iris.analysis.stats import pearsonr

#from diagnostic import plot_diagnostic
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic

logger = logging.getLogger(os.path.basename(__file__))

def compute_diagnostic(filename):
    """Compute an example diagnostic."""
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)

    logger.debug("Running example computation")
    cube = iris.util.squeeze(cube)
    return cube

def get_provenance_record(attributes, ancestor_files, plot_type):
    """Create a provenance record describing the diagnostic data and plot."""
    caption = ("Correlation of {long_name} between {dataset} and "
               "{reference_dataset}.".format(**attributes))

    record = {
        'caption': caption,
        'statistics': ['corr'],
        'domains': ['global'],
        'plot_type': plot_type,
        'authors': [
            'andela_bouwe',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': ancestor_files,
    }
    return record


def main(cfg):
    """Compute the time average for each input dataset."""
    input_data = group_metadata(
        cfg['input_data'].values(),
        'dataset',
        sort='short_name',
    )

#    print ('input_data',input_data)
    for dataset in input_data:
        logger.info("Processing model %s", dataset)
        grouped_input_data=group_metadata(input_data[dataset],'variable_group',sort='short_name')
        for short_name in grouped_input_data:
            for attributes in grouped_input_data[short_name]:
                input_file = attributes['filename']
                cube = compute_diagnostic(input_file)

#                print ('cube',cube)
#                help (cube)
                print ('short_name',short_name)
                if short_name == 'tas':
                    tas=cube.data
                    print ('tas',tas)
                if short_name == 'nbp':
                    land_cumflux=cube.data
                    print ('land_cumflux',land_cumflux)
                if short_name == 'fgco2':
                    ocean_cumflux=cube.data
                    print ('ocean_cumflux',ocean_cumflux)
                atmos_co2=0.
                print (cube.time)    
#        ocean_cumflux=np.cumsum(atmos_ocean_flux[:])/1000.
#        land_cumflux=np.cumsum(atmos_land_flux[:])/1000.
        emiss=atmos_co2+ocean_cumflux+land_cumflux

        plt.figure(figsize=[5,10])
        plt.subplot(311)
        plt.plot(year[:],q[:],color='black')
        plt.xlabel('Year')
        plt.ylabel('Global average heat flux (W/m^2)')
        plt.title('(a) Atmosphere-ocean heat flux')
        plt.subplot(312)
        plt.xlabel('Year')
        plt.ylabel('Global atmos-ocean carbon flux (PgC/yr)')
        plt.title('(b) Atmosphere-ocean carbon flux')
        plt.plot(year[:],atmos_ocean_flux[:],color='black')
        plt.subplot(313)
        plt.xlabel('Year')
        plt.ylabel('Ratio of heat to carbon flux (W/m^2/(PgC/yr))')
        plt.title('(c) Ratio of heat to carbon flux')
        plt.plot(year[:],q[:]/atmos_ocean_flux[:],color='black')
        plt.tight_layout()
        plt.savefig('/home/rng/plots/tcre/fluxes.png')
        plt.close()
        
        plt.figure(figsize=[5,10])
        plt.subplot(311)
        plt.plot(year[:],tas[:],color="black")
        plt.xlabel('Year')
        plt.ylabel('GSAT anomaly (K)')
        plt.title('(a) Warming in 4xCO2')
        plt.axis([0,250,0,10])
        plt.subplot(312)
        plt.fill_between(year[:],np.zeros(250),atmos_co2[:],color="red")
        plt.xlabel('Year')
        plt.ylabel('Atmosphere/ocean/land CO2 (EgC)')
        plt.title('(b) Carbon pool anomalies in 4xCO2')
        plt.fill_between(year[:],atmos_co2[:],atmos_co2[:]+ocean_cumflux,color="blue")
        plt.fill_between(year[:],atmos_co2[:]+ocean_cumflux,atmos_co2[:]+ocean_cumflux+land_cumflux,color="green")
        plt.plot(year[:],emiss[:],color="black",label='Cumulative CO2 emissions')
        plt.text (100,1,'Atmosphere',fontsize =11,fontweight='bold', va='center')
        plt.text (100,2.1,'Ocean',fontsize =11,fontweight='bold', va='center')
        plt.text (100,3.0,'Land',fontsize =11,fontweight='bold', va='center')
        plt.legend(loc="upper left")
        plt.axis([0,250,0,5])
        
        plt.subplot(313)
        plt.xlabel('Year')
        plt.ylabel('Ratio of warming to carbon pools')
        plt.title('(c) Ratio')
        plt.plot(year[:],tas[:]/ocean_cumflux[:],linestyle='dashed',color='black',label='\\u0394T/(ocean carbon)')
        plt.plot(year[:],tas[:]/(ocean_cumflux[:]+atmos_co2[:]),linestyle='dashdot',color='black',label='\\u0394T/(atmos+ocean carbon)')
        plt.plot(year[:],tas[:]/emiss,color="black",label='\\u0394T/(atmos+ocean+land carbon)')
        plt.legend(loc="upper right")
        plt.axis([0,250,0,18])

        plt.tight_layout()
        plt.savefig('/home/rng/plots/tcre/tcre.png')
        plt.close()

        plt.show()


                
if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
