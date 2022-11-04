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
    plot_dir=cfg['plot_dir']
    plt.ioff() #Turn off interactive plotting.

#    print ('input_data',input_data)
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
                    if exp == 'abrupt-4xCO2':
                        cube1 = compute_diagnostic(input_file)
                    if exp == 'piControl':
                        cube2 = compute_diagnostic(input_file)

#                print ('cube',cube)
#            help (cube1)
#            help (cube2)
            print ('short_name',short_name)
            if short_name == 'tas':
                tas=cube1.data-cube2.data
                print ('tas',tas)
#            if short_name == 'co2':
#                atmos_co2=(cube1.data[:,0]-cube2.data[:,0])*2.13*1e3 #Convert to Eg
#                print ('atmos_co2',atmos_co2)
            if short_name == 'nbp':
                atmos_land_flux=(cube1.data-cube2.data)*31556952/1e12 #in PgC/yr
                print ('atmos_land_flux',atmos_land_flux)
            if short_name == 'fgco2':
                print ('cube1',cube1)
                print ('cube2',cube2)
                atmos_ocean_flux=(cube1.data-cube2.data)*31556952/1e12 #in PgC/yr
                print ('atmos_ocean_flux',atmos_ocean_flux)
                print ('len(atmos_ocean_flux),len(cube1.data),len(cube2.data)',len(atmos_ocean_flux),len(cube1.data),len(cube2.data))
            if short_name == 'hfds':
                q=(cube1.data-cube2.data)/5.101e14 #Convert to W/m^2.
                print ('q',q)
#                print (cube.time)
        atmos_co2=np.full([150],1.8163314) #Increase in atmospheric CO2 in 4xCO2 simulation.
        ocean_cumflux=np.cumsum(atmos_ocean_flux[:])/1000.
        land_cumflux=np.cumsum(atmos_land_flux[:])/1000.
        emiss=atmos_co2[:]+ocean_cumflux[:]+land_cumflux[:]
        year=list(np.arange(150))
        

        plt.figure(figsize=[5,10])
        plt.subplot(311)
        plt.plot(year[:],q[:],color='black')
        plt.xlabel('Year')
        plt.ylabel('Global mean heat flux (W/m^2)')
        plt.title('(a) Atmosphere-ocean heat flux')
#Calculate values to compare.
        TCRE=1.9 #K/EgC Source: Swart et al.
        F_4co2=3.82/0.476 #Source: AR6 WGI Table 7.6
        lambda_=F_4co2/(2*5.7) #Uses ECS from Swart et al.
        delta_ca=atmos_co2[0] 
        q0=F_4co2-lambda_*TCRE*delta_ca
#        alpha=(atmos_ocean_flux[0]+atmos_land_flux[0])*lambda_*TCRE/(q0*1000)
#        f_est=(atmos_ocean_flux[0]+atmos_land_flux[0])*np.exp(-1.*np.arange(150)*alpha)
#        emiss_est=atmos_co2[:]-((atmos_ocean_flux[0]+atmos_land_flux[0])/(alpha*1000.))*(np.exp(-1.*np.arange(150)*alpha)-1)
        alpha=(atmos_ocean_flux[0])*lambda_*TCRE/(q0*1000)
        f_est=(atmos_ocean_flux[0])*np.exp(-1.*np.arange(150)*alpha)
        emiss_est=atmos_co2[:]-((atmos_ocean_flux[0])/(alpha*1000.))*(np.exp(-1.*np.arange(150)*alpha)-1)
#        print ('ocean_cumflux_est/ocean_cumflux',((atmos_ocean_flux[0])/(alpha*1000.))*(np.exp(-1.*np.arange(150)*alpha)-1)*-1./ocean_cumflux)
        q_est=q0*np.exp(-1.*np.arange(150)*alpha)
#        print ('TCRE,F_4co2,lambda_,delta_ca,q0,alpha,f0',TCRE,F_4co2,lambda_,delta_ca,q0,alpha,atmos_ocean_flux[0])
#        print ('q_est',q_est)
#        print ('emiss_est',emiss_est)
        plt.plot(year[:],q_est[:],color='gray')
        plt.axis([0,150,0,6])
        plt.subplot(312)
        plt.xlabel('Year')
        plt.ylabel('Global atmos-ocean carbon flux (PgC/yr)')
        plt.title('(b) Atmosphere-ocean carbon flux')
        plt.plot(year[:],atmos_ocean_flux[:],color='black')
        plt.plot(year[:],f_est[:],color='gray')
        plt.axis([0,150,0,100])
        plt.subplot(313)
        plt.xlabel('Year')
        plt.ylabel('Ratio of heat to ocean carbon flux (W/(PgC/yr))')
        plt.title('(c) Ratio of heat to ocean carbon flux')
        plt.plot(year[:],q[:]/(atmos_ocean_flux[:]),color='black')
#        print ('q[:]/atmos_ocean_flux[:]',q[:]/(atmos_ocean_flux[:]))
        plt.plot([0,1000],[q0/atmos_ocean_flux[0],q0/atmos_ocean_flux[0]],color='gray')
        plt.axis([0,150,0,2])
        plt.tight_layout()
        plt.savefig(plot_dir+'/'+dataset+'_ocean_fluxes.png')
        plt.close()

        plt.figure(figsize=[5,10])
        plt.subplot(311)
        plt.plot(year[:],tas[:],color="black")
        plt.plot(year[:],(q[:]*-1.+F_4co2)/lambda_,color="black",linestyle='dashed')
        plt.plot(year[:],(q_est[:]*-1.+F_4co2)/lambda_,color="gray")
#        print ('T_est/T',(q_est[:]*-1.+F_4co2)/(tas[:]*lambda_))
        plt.xlabel('Year')
        plt.ylabel('GSAT anomaly (K)')
        plt.title('(a) Warming in 4xCO2')
        plt.axis([0,150,0,12])
        plt.subplot(312)
        plt.fill_between(year[:],np.zeros(150),atmos_co2[:],color="red")
        plt.xlabel('Year')
        plt.ylabel('Atmosphere/ocean/land CO2 (EgC)')
        plt.title('(b) Carbon pool anomalies in 4xCO2')
        plt.fill_between(year[:],atmos_co2[:],atmos_co2[:]+ocean_cumflux,color="blue")
        plt.fill_between(year[:],atmos_co2[:]+ocean_cumflux,atmos_co2[:]+ocean_cumflux+land_cumflux,color="green")
        plt.plot(year[:],emiss[:],color="black",label='Cumulative CO2 emissions')
        plt.plot(year[:],emiss_est[:],color="gray",label='Cumulative CO2 emissions - analytical')
#        print ('emiss_est/emiss',emiss_est[:]/emiss[:])
        plt.text (100,1,'Atmosphere',fontsize =11,fontweight='bold', va='center')
        plt.text (100,2.1,'Ocean',fontsize =11,fontweight='bold', va='center')
        plt.text (100,3.0,'Land',fontsize =11,fontweight='bold', va='center')
        plt.legend(loc="upper left")
        plt.axis([0,150,0,6])
        
        plt.subplot(313)
        plt.xlabel('Year')
        plt.ylabel('Ratio of warming to carbon pools')
        plt.title('(c) Ratio')
#        plt.plot(year[:],tas[:]/ocean_cumflux[:],linestyle='dashed',color='black',label='\u0394T/(ocean carbon)')
        plt.plot(year[:],tas[:]/(ocean_cumflux[:]+atmos_co2[:]),linestyle='dashdot',color='black',label='\u0394T/(atmos+ocean carbon)')
        plt.plot(year[:],tas[:]/emiss,color="black",label='\u0394T/(atmos+ocean+land carbon)')
        plt.plot([0,1000],[TCRE,TCRE],color="gray",label='TCRE')
        plt.legend(loc="upper right")
        plt.axis([0,150,0,5])

        plt.tight_layout()
        plt.savefig(plot_dir+'/'+dataset+'_tcre.png')
        plt.close()
#        print ('atmos_co2',atmos_co2)


                
if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
