"""
;;#############################################################################
;; Sea Surface Temperature Diagnostics
;; Author: Benjamin Mueller (LMU Munich, GER)
;; ESA-CMUG project
;;#############################################################################
;; Description
;;    Produces various general diagnostic plots and statistics for the 
;;    reference data sets of the ESACCI Project
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    none
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    none
;;
;; Required variable_info attributes (variable specific)
;;    none
;;
;; Optional variable_info attributes (variable specific)
;;    none
;;
;; Caveats
;;
;; Modification history
;;    20161128-A_laue_ax: added call to write_references
;;    20160818-A_muel_bn: Routines written.
;;
;;#############################################################################
"""

# Basic Python packages
import sys
from copy import copy

# Add subfolder of the diagnostics to the path
sys.path.append('./diag_scripts/aux/LMU_ESACCI-diagnostics/')

from esmval_lib import ESMValProject

# Import full diagnostic routines
from sst_diagnostic import SeaSurfaceTemperatureDiagnostic

def main(project_info):
    import projects
    print(">>>>>>>> sst_ESACCI.py is running! <<<<<<<<<<<<")

# A_laue_ax+
    E = ESMValProject(project_info)

    verbosity = E.get_verbosity()
    diag_script = E.get_diag_script_name()

    res = E.write_references(diag_script,              # diag script name
                             ["A_muel_bn"],            # authors
                             [""],                     # contributors
                             [""],                     # diag_references
                             ["E_esacci-sst"],         # obs_references
                             ["P_cmug"],               # proj_references
                             project_info,
                             verbosity,
                             False)
# A_laue_ax-

    Diag=None
    
    for v in range(len(project_info['RUNTIME']['currDiag'].get_variables())):
    
        # read variable
        variable=project_info['RUNTIME']['currDiag'].get_variables()[v]
        
        #check if variable fits to diagnostics
        if variable == 'ts' or variable == 'tos':
            
            model_filelist=ESMValProject(project_info).get_clim_model_filenames(variable=variable)
    
            # only models are read
            for inc in range(len(project_info['MODELS'])):
                model=project_info['MODELS'][inc]
                currProject = getattr(vars()['projects'],
                                      model.split_entries()[0])()
                model_name = currProject.get_model_name(model)

                # only for non-reference models
                refference_model = project_info['RUNTIME']\
                                               ['currDiag'].\
                                               variables[v].\
                                               ref_model
                if not model_name == refference_model:
                
                    model_filename = model_filelist[model_name]
                    reference_filename = model_filelist[refference_model]
                
                    # copy old data to provide data that is needed again                                                                                                                                                           # copy old data to provide data that is needed again
                    D_old=copy(Diag)
                    
                    # initialize diagnostic
                    Diag = SeaSurfaceTemperatureDiagnostic()
                    # provide project_info to diagnostic
                    Diag.set_info(project_info,model,variable,reference_filename,model_filename,project_info['RUNTIME']['currDiag'].diag_script_cfg)
                    # reuse region info
                    if not D_old==None:
                        if "_regions" in D_old.__dict__.keys():
                            Diag._regions=D_old._regions
                    del(D_old)
                    # load the data
                    Diag.load_data()
                    # run the diagnostics defined by the import
                    Diag.run_diagnostic()
                    # write the results to the specific folder
                    Diag.write_data(project_info['GLOBAL']['write_plots'])
                
        if len(model_filelist)>2:
            Diag.write_overview(project_info['GLOBAL']['write_plots'])
    
    print(">>>>>>>> ENDED SUCESSFULLY!! <<<<<<<<<<<<")
    print('')

if __name__ == "__main__":
    main()
