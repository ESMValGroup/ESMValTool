from diagnostic import *
    

class LandCoverDiagnostic(BasicDiagnostics):
    """
    class to implement soil moisture diagnostics, like e.g. global means, global differences, RMSD etc.

    TODO implement testing for this diagnostic

    """
    def __init__(self, **kwargs):
        super(LandCoverDiagnostic, self).__init__(**kwargs)
        
        self._project_info={}
        self._modtype = None 
        self._reftype = None
        self._plot_dir='.' + os.sep
        self._work_dir='.' + os.sep

        self._vartype = 'land cover' #default value as there must be one
        self.output_type = 'png'  #default ouput file type 
        self._changed=False
        

    def run_diagnostic(self):
        """
        running the diagnostics
        """
        
        if '_ref_data' not in self.__dict__.keys():
            return
        
        super(LandCoverDiagnostic, self).run_diagnostic(globmeants=self.cfg.globmeants, portrait=self.cfg.portrait, globmeandiff=self.cfg.globmeandiff, trend=self.cfg.trend)

        self._specific_diag(single_years=self.cfg.single_years)
            

    def _specific_diag(self,single_years=True):
        """
        Diagnostic management
        """
        
        if single_years:
            self._year_uncertainty()


    def write_data(self,plot=True):
        """
        write data
        """
        if '_ref_data' not in self.__dict__.keys():
            return
            
        super(LandCoverDiagnostic, self).write_data() 
        
        if self.cfg.regionalization and not '_regions' in self.__dict__.keys():
            self._regions=self._mod_data.get_regions(self._reg_shape,self.cfg.shapeNames-1)
            self._write_regionalization_header()
            print "this should not be calculated here"
            
        if '_yu_data' in self.__dict__.keys():
            
            for yu in self._yu_data.keys():
            # generate map plots for years
                self._plot_yearly_maps(self._yu_data[yu], yu)               
                if self.cfg.regionalization:
                    [self._yu_data[yu][n].get_shape_statistics(self._regions) for n in np.arange(3)]
                    self._write_shape_statistics(self._yu_data[yu][0].regionalized,'mean_m'+str(self.cfg.std_factor)+'std_'+yu,self.modname+"_"+self.refname)
                    self._write_shape_statistics(self._yu_data[yu][1].regionalized,'mean_'+yu,self.modname+"_"+self.refname)
                    self._write_shape_statistics(self._yu_data[yu][2].regionalized,'mean_p'+str(self.cfg.std_factor)+'std_'+yu,self.modname+"_"+self.refname)
        else:
            print 'No yearly data to plot!'
        
    def _year_uncertainty(self):
        """ get data for uncertainty plot """
        
        #comparison for the single years
        ref_save=self._ref_data.copy()
        mod_save=self._mod_data.copy()
        
        loc_dict=dict()
        
        years=[x.year for x in self._ref_data.date]
        
        
        #only the middle year
        theseyears=[years[int(len(years)/2)]]
        
        print(theseyears)
        
        for y in theseyears:
                    
            loc_dict[str(y)]=[]

            mod_data_std=self._mod_data.copy()
            mod_data_std.data=mod_data_std.data.std(axis=0)
            
            
            low_d_r=mod_data_std.copy()
            mid_d_r=mod_data_std.copy()
            high_d_r=mod_data_std.copy()
            
            low_d=self._mod_data.data - self.cfg.std_factor*mod_data_std.data - self._ref_data.data.mean(axis=0)
            low_d.mask=(self._mod_data.data.mask + mod_data_std.data.mask + self._ref_data.data.mask).astype(bool)
            low_d_r.data=low_d
            loc_dict[str(y)].insert(0,low_d_r)
            
            mid_d=self._mod_data.data - self._ref_data.data.mean(axis=0)
            mid_d.mask=(self._mod_data.data.mask + mod_data_std.data.mask + self._ref_data.data.mask).astype(bool)
            mid_d_r.data=mid_d
            loc_dict[str(y)].insert(1,mid_d_r)
            
            high_d=self._mod_data.data + self.cfg.std_factor*mod_data_std.data - self._ref_data.data.mean(axis=0)
            high_d.mask=(self._mod_data.data.mask + mod_data_std.data.mask + self._ref_data.data.mask).astype(bool)
            high_d_r.data=high_d
            loc_dict[str(y)].insert(2,high_d_r)
                        
            self._ref_data=ref_save.copy()
            self._mod_data=mod_save.copy()
        
        self._yu_data=loc_dict
        
    def _plot_yearly_maps(self,low_mid_high,year):
        """
        plot kendall's tau trend correlations
        """
        f = plt.figure(figsize=(30,6))
        ax1 = f.add_subplot(131)
        ax2 = f.add_subplot(132)
        ax3 = f.add_subplot(133)
        
        def submap(data,ax,title,vmin,vmax,cmap,ctick={'ticks': None, 'labels': None}):
            Map = SingleMap(data,
                            backend=self.plot_backend,
                            show_statistic=True,
                            savefile=None,
                            ax=ax,
                            show_unit=False)
            Map.plot(title= title,
                     show_zonal=False,
                     show_histogram=False,
                     show_timeseries=False,
                     nclasses=self.plot_nclasses,
                     colorbar_orientation=self.plot_cborientation,
                     show_colorbar=self.plot_cbshow, 
                     cmap=cmap, 
                     vmin=min(self.cfg.mima_single_year),
                     vmax=max(self.cfg.mima_single_year), 
                     proj_prop=self.cfg.projection, 
                     ctick_prop=ctick,
                     drawparallels=True,
                     titlefontsize=self.plot_tfont)

        submap(low_mid_high[0], ax=ax1, title="mean - " + str(self.cfg.std_factor) + "*std",vmin=-1,vmax=1,cmap='RdBu')
        submap(low_mid_high[1], ax=ax2, title="mean",vmin=0,vmax=1,cmap='RdBu')
        submap(low_mid_high[2], ax=ax3, title="mean + " + str(self.cfg.std_factor) + "*std",vmin=-1,vmax=1,cmap='RdBu')
        f.suptitle("Difference between " + self.modname + " and " + self.refname + " within the range of the " + str(self.cfg.std_factor) + "th standard deviation of the model for years around " + year) #TODO get the year range right

        oname = self._get_output_rootname() + "_diff_std_" + year + '.' + self.output_type
        if os.path.exists(oname):
            os.remove(oname)
        f.savefig(oname,dpi=self.plot_dpi)

        plt.close(f.number)  # close figure for memory reasons!
        del f
        
        
    def _load_model_data(self):
        """ load all land cover model data """

        edited=False
        
        newfile=self._mod_file+".T63built.nc"
        newfile=newfile.split("/")
        newdir=(self._work_dir if self._work_dir[-1]==os.sep else self._work_dir +os.sep) + "AUX_Files_lc_ESACCI"
        newfile=newdir + os.sep + newfile[-1]
         
        mod_info=Dataset(self._mod_file)
        lat=mod_info.dimensions['lat'].size #TODO this is not right from preprocessing!!!
        lon=mod_info.dimensions['lon'].size
        mod_info.close()
        
        if not (lat == 96 and lon == 192): #TODO add diffs
        
            if not os.path.exists(newfile):
                tempfile=self._aggregate_resolution(self._mod_file,"T63",remove=False)
                subprocess.call(["mkdir",newdir])
                subprocess.call(['cp',tempfile,newfile])
                os.remove(tempfile)
                
            self._mod_file_E=newfile
            edited=True
                
        #load data
        self._mod_data = self._load_cmip_generic(self._mod_file_E if edited else self._mod_file,self.var)
#        self._mod_data.data=np.nan_to_num(self._mod_data.data)

            
    def _load_observation_data(self):
        """ load obs data """
        newfile=self._ref_file+".T63built.nc"
        newfile=newfile.split("/")
        newdir=(self._work_dir if self._work_dir[-1]==os.sep else self._work_dir +os.sep) + "AUX_Files_lc_ESACCI"
        newfile=newdir + os.sep + newfile[-1]
        
        mod_info=Dataset(self._ref_file)
        lat=mod_info.dimensions['lat'].size
        lon=mod_info.dimensions['lon'].size
        mod_info.close()
        
        if not (lat == 96 and lon == 192): #TODO add diffs
            if not os.path.exists(newfile):
                tempfile=self._aggregate_resolution(self._ref_file,"T63",remove=False)
                subprocess.call(["mkdir",newdir])
                subprocess.call(['cp',tempfile,newfile])
                os.remove(tempfile)
                
            self._ref_file=newfile

        if self.var in ["baresoilFrac","grassNcropFrac","shrubNtreeFrac"]:
            self._ref_data=self._load_cci_generic(self._ref_file,self.var)
#            self._ref_data.data=np.nan_to_num(self._ref_data.data)
        else:
            assert False, 'Not supported yet'

    def _load_regionalization_shape(self):
        """ load shape data """
        self._reg_shape = self._load_shape_generic(self._reg_file)
      
    def _preprocess_data(self, infiles, mod,check_f = None,force=False):
        """
        preprocess observations and model to adapt to each others specification and temporal and spatial resolution needed
        
        Parameters:
        -----------
        
        infiles : string with path to files
        mod : exemplary model data; the prepocessing should mirror its specifications
        check_f : alternativeley check this folder and write infiles.built.nc
        force : also preprocess if file is there
        """
        
        #preprocess observation data
        varlist=self.cfg.translatorlist.keys()
        varlist=['_'.join(var.split(" ")) for var in varlist] 
    
        ofiles=infiles
        
        ofiles_built=[ofile + ".built.nc" for ofile in ofiles]
        
        #choose timestep and resolution 
        #TODO get this from mod_data or proj_info           
        resolution="T63"
        
        original_checker=[os.path.isfile(ofile) for ofile in ofiles]
        built_checker=[os.path.isfile(ofile) for ofile in ofiles_built]
        
        original_checker=len([x for x in original_checker if x]) == len(original_checker)
        built_checker=len([x for x in built_checker if x]) == len(built_checker)
        
        if ( original_checker or built_checker ) and not force:
            print("No LANDCOVER preprocessing needed!")
            
        elif (not (original_checker or built_checker) or force):
            #extracting with tool
            file_list_a, list_length_a = self._get_files_in_directory(check_f,'*aggregate*.nc',False)
            
            if list_length_a==0: 
                # test if there is any aggregated file. Only if there is none, the files will be produced!
            
                file_list, list_length = self._get_files_in_directory(check_f,'*.nc',False)
                
                if list_length==0:
                    print infiles
                    assert False, "cannot find any files!" 
                
                print('   running ' + path2lctool + ' for '+ str(list_length) + ' files ...')
    
                for locfile in file_list:
                            
                    preprocessingCommand= "bash " + self._work_dir + os.sep + self.cfg.path2lctool + os.sep + "bin" + os.sep + "aggregate-map.sh -PgridName=GEOGRAPHIC_LAT_LON -PnumMajorityClasses=1 -PoutputAccuracy=false -PoutputPFTClasses=true -PoutputLCCSClasses=false " + locfile
                
                    process = subprocess.Popen(preprocessingCommand.split(), stdout=subprocess.PIPE)
                    for line in iter(process.stdout.readline, b''):
                        print line,
                        process.stdout.close()
                    process.wait()                    
                    
                file_list_a, list_length_a = self._get_files_in_directory(self._work_dir + os.sep + self.cfg.path2ref + os.sep +self.cfg.ref_folder,'*aggregate*.nc',False)
            
    
            #aggregating for classes
            for it in range(len(self.cfg.translatorlist)):
                
                t_list=[]            
                
                for fl in file_list_a:
                
                    t_list.append(self._extract_variables(fl,self.cfg.translatorlist.items()[it][1][0],self.cfg.translatorlist.items()[it][0]))    
    
                cdo=Cdo()
                
                varfilename=self.cfg.ref_data
                varfilename=varfilename.split("*")
                varfilename[1]='_'.join(self.cfg.translatorlist.items()[it][0].split(" "))
                varfilename=''.join(varfilename)
                
                oname=self._work_dir + os.sep + "temp" + os.sep + varfilename
    
                cdo.mergetime(input=t_list,output=oname,options='-f nc4 -b F32')
                
                for ifi in t_list:
                    os.remove(ifi)
                
            #adjusting resolution and representation      
            file_list, list_length = self._get_files_in_directory(self._work_dir + os.sep + "temp" + os.sep,'*.nc',False)
            
            for locfile in file_list:
                self._adjust_lc_details(locfile,resolution)
            
                ofile=self._work_dir + os.sep + self.cfg.path2ref + os.sep + locfile.split(os.sep)[-1] + ".built.nc"
                
                subprocess.call(['cp',locfile,ofile])
                
                os.remove(locfile) 
                
                print("   Data written!")

        #load one file for details and test
        test_data=self._load_cci_generic(ofiles[0] if os.path.isfile(ofiles[0]) else ofiles_built[0],varlist[0])
        test_data_years=[test_data.date[x].year for x in range(len(test_data.date)) ]
        
        
        #preprocess model data
        new_mod_files_full=self._mod_file.split(os.sep)
        new_mod_files_sub=new_mod_files_full[-1].split("_")
        new_mod_files_full.remove(new_mod_files_full[-1])
        new_mod_files_sub.remove("lc")
        subprocess.call(['mkdir',"/".join(new_mod_files_full) + os.sep + "sub"])
        self._mod_files=["/".join(new_mod_files_full) + os.sep + "sub" + os.sep + var + "_" +  "_".join(new_mod_files_sub) for var in varlist]
        
        mod_checker=[os.path.isfile(files) for files in self._mod_files]
        mod_checker=len([x for x in mod_checker if x]) == len(mod_checker)
        
        if not mod_checker:
            for it in self.cfg.translatorlist.keys():
                alt_mod_files=self.cfg.translatorlist[it][1]
                alt_mod_files=[self.single_files[var].filename for var in alt_mod_files]
                if len(alt_mod_files)>1:
                    thisfile = self._sum_files(alt_mod_files,remove=False)
                else:
                    thisfile = alt_mod_files[0]
            
                #adjust timestep
                thisfile = self._rename_variable(thisfile,"newFrac")
                thisfile = self._aggregate_specific_years(thisfile,test_data_years,remove=False)
    
                
                #set new name for file
                ofile=self._mod_files[self.cfg.translatorlist.keys().index(it)]
    
                subprocess.call(['cp',thisfile,ofile])
                
                os.remove(thisfile)

                print("   Data written!")

        #run LC-Classes seperately
        for it in range(len(ofiles)):
            nml="namelist_" + varlist[it] + "_ESACCI.xml"
            print(nml + " started!")
            preprocessingCommand="python main.py nml/" + nml
            process = subprocess.Popen(preprocessingCommand.split(), stdout=subprocess.PIPE)
            for line in iter(process.stdout.readline, b''):
                print line,
            process.stdout.close()
            process.wait()
            
    def _adjust_lc_details(self,infile,resolution,remove=True):
        """ aggregate infile to resolution and adjust bounding box """
        """ currenty only T63 """
        cdo=Cdo()
        oname=infile
        tmp1name=self._work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1]
        tmp2name=self._work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1]
        if resolution=="T63":
            thisfile=self._aggregate_resolution(infile,resolution)
            cdo.invertlat(input=thisfile,output=tmp1name)
            cdo.sellonlatbox("0,360,-90,90",input=tmp1name,output=tmp2name)
            cdo.mulc(100,input=tmp2name,output=oname)
            if remove:
                os.remove(thisfile)
                os.remove(tmp1name)
                os.remove(tmp2name)
        else:
            assert False, "This resolution cannot be handled yet."
            
        return oname

     
