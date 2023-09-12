"""Python code to plot supplementary figures in Gillett (2023)."""
import logging
import os
import csv
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
            'nathan_gillett',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': ancestor_files,
    }
    return record


def main(cfg):
    input_data = group_metadata(
        cfg['input_data'].values(),
        'dataset',
        sort='short_name',
    )
    plot_dir=cfg['plot_dir']
    plt.ioff() #Turn off interactive plotting.
# Define variables.
    models=['CanESM5']
    TCRE=[2.09] #Arora et al. 2020, Table 4. Assume same TCRE for all CESM2 models.
    nmodel=len(models)
    plt.figure(1,figsize=[13,2.5])
    plt.figure(2,figsize=[13,2.5])

# Read in data (note only plots output for CanESM5 1pctCO2 and zero emissions simulation).
    for dataset in input_data:
        logger.info("Processing model %s", dataset)
        grouped_input_data=group_metadata(input_data[dataset],'variable_group',sort='exp')
        for short_name in grouped_input_data:
            logger.info("Processing variable %s", short_name)
            grouped_input_data_var=group_metadata(grouped_input_data[short_name],'exp',sort='exp')
            
            for exp in grouped_input_data_var:
                logger.info("Processing experiment %s", exp)
                for attributes in grouped_input_data_var[exp]:
                    input_file = attributes['filename']
                    if exp == '1pctCO2':
                        cube1 = compute_diagnostic(input_file)
                    if exp == 'esm-1pct-brch-1000PgC':
                        cube3 = compute_diagnostic(input_file)
                    if exp == 'piControl':
                        cube2 = compute_diagnostic(input_file)
            if short_name != 'co2':
               merged=np.zeros((166))
               merged[0:61]=cube1.data[:]
               merged[61:166]=cube3.data[:]
            else:
               merged=np.zeros((166))
               merged[0:61]=cube1.data[:,0]
               merged[61:166]=cube3.data[:,0]
                
            if short_name == 'tas':
                tas=merged-cube2.data
            if short_name == 'co2':
                atmos_co2=(merged-cube2.data[:,0])*2.13*1e3 #Convert to Eg
            if short_name == 'nbp':
                atmos_land_flux=(merged-cube2.data)*31556952/1e12 #in PgC/yr
            if short_name == 'fgco2':
                atmos_ocean_flux=(merged-cube2.data)*31556952/1e12 #in PgC/yr
            if short_name == 'hfds':
                q=(merged-cube2.data)/5.101e14 #Convert to W/m^2.

        ocean_cumflux=np.cumsum(atmos_ocean_flux[:])/1000.
        land_cumflux=np.cumsum(atmos_land_flux[:])/1000.
        emiss=atmos_co2[:]+ocean_cumflux[:]+land_cumflux[:]
        year=list(np.arange(166))
        mod_index=models.index(dataset)        

#Plot Supplementary Fig 2.
        plt.figure(1)
        plt.subplot(nmodel,3,mod_index*3+1)
        plt.plot(year[:],q[:],color='black')
        plt.ylabel('Heat flux (W/$m^2$)')

        plt.axis([0,150,0,2])
        plt.text (5,0.9*2,'a',fontweight='bold', va='center')
        if mod_index==nmodel-1:
          plt.xlabel('Year')
        plt.text (20,7,dataset,fontsize =10,fontweight='bold', va='center')
        if mod_index ==0:
          plt.title('Heat flux')
        if mod_index < nmodel-1:
          plt.xticks([])
        plt.subplot(nmodel,3,mod_index*3+2)
        plt.ylabel('Carbon flux (PgC/yr)')
        plt.plot(year[:],atmos_ocean_flux[:],color='black')
        plt.axis([0,150,0,7])
        plt.text (5,0.9*7,'b',fontweight='bold', va='center')
        if mod_index == 0:
          plt.title('Carbon flux')
        if mod_index == nmodel-1:
          plt.legend(loc="upper right",labelspacing=0.1,borderpad=0.1,fontsize=10)
          plt.xlabel('Year')
        if mod_index < nmodel-1:
          plt.xticks([])
        plt.subplot(nmodel,3,mod_index*3+3)
        plt.ylabel('(W/$m^2$)/(PgC/yr)')
        plt.plot(year[:],q[:]/(atmos_ocean_flux[:]),color='black')
        plt.axis([0,150,0,2])
        plt.text (5,0.9*2,'c',fontweight='bold', va='center')
        if mod_index == 0:
          plt.title('Ratio')
        if mod_index == nmodel-1:
          plt.xlabel('Year')
        if mod_index < nmodel-1:
          plt.xticks([])

#Plot Supplementary Fig 1.
        plt.figure(2)
        plt.subplot(nmodel,3,mod_index*3+1)
        plt.plot(year[:],tas[:],color="black",label='ESM')#,label='\u0394T')
        plt.ylabel('\u0394T (K)')
        if mod_index==nmodel-1:
          plt.xlabel('Year')
#          plt.legend(loc="center left",labelspacing=0.1,borderpad=0.1,fontsize=10)
        plt.text (5,10,dataset,fontsize =10,fontweight='bold', va='center')
        if mod_index ==0:
          plt.title('Warming')
        plt.axis([0,150,0,3])
        plt.text (15,0.9*3,'a',fontweight='bold', va='center')
        if mod_index < nmodel-1:
            plt.xticks([])
        plt.subplot(nmodel,3,mod_index*3+2)
        plt.fill_between(year[:],np.zeros(166),atmos_co2[:],color="red")
        plt.fill_between(year[:],atmos_co2[:],atmos_co2[:]+ocean_cumflux,color="blue")
        plt.fill_between(year[:],atmos_co2[:]+ocean_cumflux,atmos_co2[:]+ocean_cumflux+land_cumflux,color="green")
        plt.plot(year[:],emiss[:],color="black",label='Cumulative CO2 emissions')
        plt.ylabel('Carbon (EgC)')
        if mod_index == 0:
          plt.title('Cumulative emissions')
        plt.axis([0,150,0,2])
        plt.text (15,0.9*2,'b',fontweight='bold', va='center')
        if mod_index < nmodel-1:
          plt.xticks([])
        plt.subplot(nmodel,3,mod_index*3+3)        
        plt.ylabel('Ratios (K/EgC)')
        plt.plot(year[:],0.1*tas[:]/ocean_cumflux[:],color='blue',label='0.1\u0394T/\u0394$C_O$')
        plt.plot(year[:],tas[:]/(ocean_cumflux[:]+atmos_co2[:]),color='red',label='\u0394T/(\u0394$C_O$+\u0394$C_A$)')
        plt.plot(year[:],tas[:]/emiss,color="black",label='\u0394T/E')
        plt.plot([0,1000],[TCRE[mod_index],TCRE[mod_index]],color="gray")
        plt.axis([0,150,0,4])
        plt.text (15,0.9*4,'c',fontweight='bold', va='center')
        if mod_index == 0:
          plt.title('Ratios')
        if mod_index < nmodel-1:
          plt.xticks([])
        if mod_index == nmodel-1:        
          plt.xlabel('Year')
          plt.legend(loc="lower left",labelspacing=0.1,borderpad=0.1,fontsize=10)
#Write data out to a file, one per model.
        with open(plot_dir+'/gillett23_supp_figs1and2_data_'+dataset+'.csv', mode='w') as file:
          data_writer=csv.writer(file,delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
          data_writer.writerow(['Year','Delta T','Delta C_A','Delta C_A + Delta C_O',\
                                'E','0.1*Delta T / Delta C_O','Delta T/(Delta C_A+DeltaC_O)','Delta T/E','TCRE',\
                                'q','f','q/f'])
          for yy in range(150):
            data_writer.writerow([yy,tas[yy],atmos_co2[yy],atmos_co2[yy]+ocean_cumflux[yy],\
                                  emiss[yy],0.1*tas[yy]/ocean_cumflux[yy],tas[yy]/(ocean_cumflux[yy]+atmos_co2[yy]),tas[yy]/emiss[yy],TCRE[mod_index],\
                                  q[yy],atmos_ocean_flux[yy],q[yy]/atmos_ocean_flux[yy]])

    plt.figure(1)
    plt.tight_layout(pad=0.1)
    plt.savefig(plot_dir+'/ocean_fluxes.png')
    plt.close()

    plt.figure(2)
    plt.tight_layout(pad=0.1)
    plt.savefig(plot_dir+'/tcre.png')
    plt.close()

                
if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
