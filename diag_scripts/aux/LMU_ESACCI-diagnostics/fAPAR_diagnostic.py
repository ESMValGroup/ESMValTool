from diagnostic import *


class fAPARDiagnostic(BasicDiagnostics):
    """
    class to implement albedo diagnostics, like e.g. global means, global differences, RMSD etc.

    TODO implement testing for this diagnostic

    """
    def __init__(self, **kwargs):
        super(fAPARDiagnostic, self).__init__(**kwargs)
        
        self._project_info={}
        self._modtype = None 
        self._reftype = None
        self._plot_dir='.' + os.sep
        self._work_dir='.' + os.sep

        self._vartype = 'albedo' #default value as there must be one
        self.output_type = 'png'  #default ouput file type 
        self._changed=False

    def run_diagnostic(self):
        """
        running the diagnostics
        """
        super(fAPARDiagnostic, self).run_diagnostic()

    def _specific_diag(self, percentile=True):
        """
        Diagnostic management
        """

    def write_data(self,plot=True):
        """
        write data
        """
        super(fAPARDiagnostic, self).write_data()

        
    def _load_model_data(self):
        """ load albedo model data """
        
        edited=False
        
        newfile=self._mod_file+".T63built.nc"
        newfile=newfile.split("/")
        newdir=(self._work_dir if self._work_dir[-1]==os.sep else self._work_dir +os.sep) + "AUX_Files_fAPAR_QA4ECV"
        newfile=newdir + os.sep + newfile[-1]
         
        mod_info=Dataset(self._mod_file)
        lat=mod_info.dimensions['lat'].size
        lon=mod_info.dimensions['lon'].size
        mod_info.close()
        
        if not ((lat == 96 and lon == 192) or (lat == 6 and lon == 12)): #TODO add diffs
        
            if not os.path.exists(newfile):
                tempfile=self._aggregate_resolution(self._mod_file,"T63",remove=False)
                subprocess.call(["mkdir",newdir])
                subprocess.call(['cp',tempfile,newfile])
                os.remove(tempfile)
                
            self._mod_file_E=newfile
            edited=True
                
        #load data
        self._mod_data = self._load_cmip_generic(self._mod_file_E if edited else self._mod_file,self._project_info['RUNTIME']['currDiag'].get_variables()[0])
        self._mod_data.unit="-"   
        
        
    def _load_observation_data(self):
        """ load obs data """
        newfile=self._ref_file+".T63built.nc"
        newfile=newfile.split("/")
        newdir=(self._work_dir if self._work_dir[-1]==os.sep else self._work_dir +os.sep) + "AUX_Files_fAPAR_QA4ECV"
        newfile=newdir + os.sep + newfile[-1]
        
        mod_info=Dataset(self._ref_file)
        lat=mod_info.dimensions['lat'].size
        lon=mod_info.dimensions['lon'].size
        mod_info.close()
        
        if not ((lat == 96 and lon == 192) or (lat == 6 and lon == 12)): #TODO add diffs
            if not os.path.exists(newfile):
                tempfile=self._aggregate_resolution(self._ref_file,"T63",remove=False)
                subprocess.call(["mkdir",newdir])
                subprocess.call(['cp',tempfile,newfile])
                os.remove(tempfile)
                
            self._ref_file=newfile
            
        if self._project_info['RUNTIME']['currDiag'].get_variables()[0] == "fAPAR":
            self._ref_data=self._load_cci_generic(self._ref_file,self._project_info['RUNTIME']['currDiag'].get_variables()[0])
        else:
            assert False, 'Not supported yet'
        self._ref_data.unit="-" 

    def _load_regionalization_shape(self):
        """ load shape data """
        self._reg_shape = self._load_shape_generic(self._reg_file)
        
    def _preprocess_observations(self, infile, mod, var,check_f = None,force=False):
        """
        preprocess observations to adapt to temporal and spatial resolution needed
        Parameters:
        -----------
        
        infile : string with path to file
        mod : model data; the prepocessing should mirror its specifications
        var : name of the variable within the file
        check_f : alternativeley check this folder and write infile.built.nc
        """
        
        #choose timestep and resolution 
        #TODO get this from mod_data or proj_info
        if self._field_type in ["T2Ms", "T3M"]:
            timestep="monthly"
        else:
            assert False, "unnkown field_type"
            
        resolution="T63"
        
        ofile=infile+'.built.nc'
        
        if (os.path.isfile(infile) or os.path.isfile(ofile)) and not force:
            data = self._load_cci_generic(ofile if os.path.isfile(ofile) else infile,var)
            return data
        elif (os.path.isfile(infile) or os.path.isfile(ofile)) and force:
            
            #adjust timestep
            thisfile = self._aggregate_timestep(ofile if os.path.isfile(ofile) else infile,timestep)
            
            #adjust resolution
            thisfile = self._aggregate_resolution(thisfile,resolution)
            
            subprocess.call(['cp',thisfile,ofile])
            
            os.remove(thisfile)
            
            data = self._load_cci_generic(ofile,var)
            return data
            
        elif not check_f is None:
            file_list, list_length = self._get_files_in_directory(check_f,'*/*/*/*.nc',False)
            
            #choose necessary files
            def loc_timestamp_split(str):
                return str.split('/')[-1].split('-')[0]
                
            file_timestamps = np.asarray(map(int,map(loc_timestamp_split,file_list)))
            ys=file_timestamps/10000000000
            low_date = min(ys)
            high_date = max(ys)
            
            curyear=low_date
                       
            file_list_agg=[]
            
            while not curyear == (high_date + 1):
                
                use = np.array(np.where(np.all([file_timestamps/10000000000 == curyear],0)))
                use=use[0][(np.argsort(file_timestamps[use[0]]))]
                file_list_cur=np.array(file_list)[use].tolist()
                
                if len(file_list_cur)>0:
                    #concatenate files
                    thisfile = self._aggregate_obs_from_files(file_list_cur)
                    
                    #adjust mask
                    thisfile = self._apply_sst_flags(thisfile)
                
                    #adjust timestep
                    thisfile = self._aggregate_timestep(thisfile,timestep)
                    
                    #adjust resolution
                    thisfile = self._aggregate_resolution(thisfile,resolution)
                    
                    #add to new file list
                    file_list_agg.append(thisfile)
                
                curyear = curyear + 1
                
                print(file_list_agg)

            thisfile = self._aggregate_obs_from_files(file_list_agg)
            
            [os.remove(tf) for tf in file_list_agg]

            subprocess.call(['cp',thisfile,ofile])
            
            os.remove(thisfile)
            
            data = self._load_cci_generic(ofile,var)
            return data
            
        else:
            print infile
            assert False, "cannot find any files!" 
            
    def _apply_sst_flags(self,infile,remove=True):
        """ apply flag to file """
        """ currenty only monthly """
        cdo=Cdo()
        oname=self._work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1]
        tmpname1=self._work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1]
        tmpname2=self._work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1]
        tmpname4=self._work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1] #only if only sm is wanted
        cdo.selname("mask",input=infile,output=tmpname1,options='-L -f nc4 -b F32')
        cdo.selname("analysed_sst",input=infile,output=tmpname4,options='-L -f nc4 -b F32') #only if only sm is wanted
        cdo.setvrange(1,1,input=tmpname1,output=tmpname2,options='-L -f nc4 -b F32')
        os.remove(tmpname1)
        cdo.div(input=tmpname4 + " " + tmpname2,output=oname,options='-L -f nc4 -b F32') #only if only sm is wanted
        #cdo.div(input=infile + " " + tmpname3,output=oname,options='-L -f nc4 -b F32')
        os.remove(tmpname2)
        os.remove(tmpname4) #only if only sm is wanted
        
        if remove:
            os.remove(infile)
                
            
        return oname
        
###TODO CURRENTLY UNUSED

    def download_observations(self, year, force=True):
        """
        download observations and update local mirror with recent data

        The data is currently downloaded from the ICDC (http://icdc.zmaw.de/1.html)
        This requires a special account. An alternative approach would
        be to download the data directly from the SM CCI project
        and/or the CCI data protal

        TODO implement download from CCI or data portal
        """
        if not os.path.exists(self.raw_reference_directory):
            os.makedirs(self.raw_reference_directory)

        odir = self.raw_reference_directory + os.sep + year + os.sep
        if os.path.exists(odir):
            if force == False:
                print('No data download as data already available: ' + year)
                return

        if not os.path.exists(odir):
            os.makedirs(odir)
        else:  # empty directory
            os.system('rm -rf ' + odir + '*.nc')

        # download data via scp
        # if on the server, create only links, otherwise download using scp
        if os.path.exists(self.remote_dir):
            print('Linking data only as on server!')
            cmd = "ln -s " + self.remote_dir + year + '/*.nc' + ' ' + odir
            os.system(cmd)
        else:
            assert False, "not working yet with subprocess"
            cmd = self.user + '@' + self.server + ':' + self.remote_dir + year + '/*.nc' + ' ' + odir
            subprocess.call(["scp", "-r", cmd], shell=True)

    def update_observation(self, year, force=False):
        """
        download and update observations for a particular year

        Parameters
        ----------
        year : str
            year to process
        force : bool
            specifies if processing should be done in any case (true) or if
            the data is available, no processing will be made (false)
        """
        ofile_root = self.reference_preproc_directory + os.sep + 'CCI_SM_' + year

        self.download_observations(year, force=force)
        self.preprocess_observations(year, ofile_root, force=force)


