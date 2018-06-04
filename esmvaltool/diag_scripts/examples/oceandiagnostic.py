"""Python example diagnostic."""
import logging
import inspect, os
import sys

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import yaml


#from esmvaltool.diag_scripts.shared import run_diagnostic
#from esmvaltool.diag_scripts.shared.plot import example_map_plot

logger = logging.getLogger(os.path.basename(__file__))



def folder(name):
        """
        This snippet takes a string, makes the folder and the string.
        It also accepts lists of strings.
        """
        if type(name) == type(['a','b','c']):
                name='/'.join(name)
        if name[-1] != '/':
                name = name+'/'
        if os.path.exists(name) is False:
                os.makedirs(name)
                logger.info('Making new directory:', name)
        return name


          
            
def get_cfg():
    """Read diagnostic script configuration from settings.yml."""
    settings_file = sys.argv[1]
    with open(settings_file) as file:
        cfg = yaml.safe_load(file)
    return cfg, settings_file

def get_input_files(cfg, index=0):
    """Get a dictionary with input files from metadata.yml files."""
    metadata_file = cfg['input_files'][index]
    with open(metadata_file) as file:
        metadata = yaml.safe_load(file)
    return metadata

def get_image_path(cfg, md, prefix='',suffix='',image_extention = 'png',
         basenamelist = ['project', 'model','mip','exp','ensemble','field','short_name','start_year','end_year']
         ):
        """
        This produces a path to the final location of the image.
        The cfg is the opened global config,
        md is the metadata dictionairy (for the individual model file)
        """
       
        path  = folder(cfg['plot_dir'])
        if prefix: path+=prefix+'_'
        path += '_'.join([str(md[b]) for b in basenamelist])
        if suffix: path += '_'+suffix
        path += '.'+ image_extention
        logger.info("Image path will be ", path)
        return path
       

def loadCubeAndSavePng(cfg,md, fn, plotType = 'Mean'):
        """
        This function makes a simplee plot for an indivudual model.
        The cfg is the opened global config,
        md is the metadata dictionairy (for the individual model file).
        """
        
        cube = iris.load_cube(fn)
        cube.convert_units('celsius')            
        
        print('Attempting to take time mean in iris.')
        multiModel = md['model'].find('MultiModel')>-1        
        if plotType == 'Mean':
                cube = cube.collapsed('time', iris.analysis.MEAN)               
                qplt.contourf(cube, 25)
                try:    plt.gca().coastlines()
                except: pass
                
                                
        if plotType == 'WeightedMean':          
                weights =  iris.analysis.cartography.area_weights(cube, normalize=False) 
                coords = ['longitude', 'latitude']
                cube = cube.collapsed(coords, iris.analysis.MEAN, weights=weights,)

                if multiModel: 
                    qplt.plot(cube,label = md['model'],ls=':')                    
                else:
                    qplt.plot(cube,label = md['model'])
                                    
                plt.legend(loc='best')
                
        try:    path = get_image_path(cfg, md, suffix=plotType, image_extention = 'png',)
        except: 
                path = folder(cfg['plot_dir']) + os.path.basename(fn).replace('.nc','_'+plotType+'.png')
        
        plt.title(md['model'] +' '+ md['long_name'])
        
        if cfg['write_plots']:               
                #plt.show()
                
                logger.info('Saving plots to',path)
                plt.savefig(path)
                
        plt.close()



def multiModelTimeSeries(cfg,metadata, plotType = 'WeightedMean'):
        """
        This method makes a time series plot showing several models.
        """
        mdfiles = []
        for fn in sorted(metadata.keys()):
            
            print('-----------------')
            print('model filenames:',fn,)
            #print('Varbiable name:', key)
            print('metadata[fn]:',metadata[fn])        
            cube = iris.load_cube(fn)
            cube.convert_units('celsius')
        
            multiModel = metadata[fn]['model'].find('MultiModel')>-1
            
            if plotType == 'WeightedMean':
                weights =  iris.analysis.cartography.area_weights(cube, normalize=False)
                coords = ['longitude', 'latitude']
                cube = cube.collapsed(coords, iris.analysis.MEAN, weights=weights,)

                if multiModel: 
                    qplt.plot(cube,label = metadata[fn]['model'],ls=':')                    
                else:
                    qplt.plot(cube,label = metadata[fn]['model'])
                plt.legend(loc='best')
           
                
        if cfg['write_plots']:               
                path = get_image_path(cfg, metadata[fn], prefix = 'MultiModel' , suffix= plotType, image_extention = 'png',
	                basenamelist = ['field','short_name','start_year','end_year'], )
                logger.info('Saving plots to',path)
                plt.savefig(path)
                
        plt.close()
                

def main():
    #####
    # This part sends debug statements to stdout
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))
         cfg, settings_file = get_cfg()
    #logger.setLevel(cfg['log_level'].upper())          
    
    input_files = get_input_files(cfg)

    print('CFG:\tsettings_file\t',settings_file)
    print('CFG:\tContents:')
    for k in cfg.keys(): 
        print('CFG:\t',k,'\t',cfg[k])

    for i,metadatafilename in enumerate(cfg['input_files']):
        print('\n\metadata filename:',metadatafilename,)
        metadata = get_input_files(cfg,index=i)

        #######
        # Multi model time series
        multiModelTimeSeries(cfg,metadata, plotType = 'WeightedMean')
            
        fns = sorted(metadata.keys())
        for fn in fns:
                
                print('-----------------')
                print('model filenames:',fn,)
                print('metadata[fn]:',metadata[fn])

                ######
                # Time series of individual model
                loadCubeAndSavePng(cfg,metadata[fn], fn, plotType = 'WeightedMean')            
                
                ######
                # Map of temperal mean
                loadCubeAndSavePng(cfg,metadata[fn], fn, plotType = 'Mean')                
                            
    logger.debug("\n\nThis works\n\n")
    print('Success')    

        
        

            

#cfg, settings = get_cfg()
#metadata = get_input_files(cfg)

if __name__ == '__main__':
    main()
    
    
