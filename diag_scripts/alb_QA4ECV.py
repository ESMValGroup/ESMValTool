"""
;;#############################################################################
;; Albedo Diagnostics
;; Author: Benjamin Mueller (LMU Munich, GER)
;; QA4ECV project
;;#############################################################################
;; Description
;;    Produces various general diagnostic plots and statistics for the 
;;    reference data sets of the QA4ECV Project
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
#import imp

# Add subfolder of the diagnostics to the path
sys.path.append('./diag_scripts/aux/LMU_ESACCI-diagnostics/')

from esmval_lib import ESMValProject

# Import full diagnostic routine
from alb_diagnostic import AlbedoDiagnostic

def main(project_info):
    print('>>>>>>>> alb_QA4ECV.py is running! <<<<<<<<<<<<')

# A_laue_ax+
    E = ESMValProject(project_info)

    verbosity = E.get_verbosity()
    diag_script = E.get_diag_script_name()

#    res = E.write_references(diag_script,              # diag script name
#                             ["A_muel_bn"],            # authors
#                             [""],                     # contributors
#                             [""],                     # diag_references
#                             ["E_qa4ecv_albedo"],   # obs_references
#                             ["P_qa4ecv"],               # proj_references
#                             project_info,
#                             verbosity,
#                             False)
# A_laue_ax-

#    f = open(project_info['RUNTIME']['currDiag'].diag_script_cfg)
#    cfg = imp.load_source('cfg', '', f)
#    f.close()

    Diag=None

    for v in range(len(project_info['RUNTIME']['currDiag'].get_variables())):
    
        # read variable
        variable=project_info['RUNTIME']['currDiag'].get_variables()[v]
        
        if variable in ["alb"]:
            
            model_filelist=ESMValProject(project_info).get_clim_model_filenames(variable=variable)
            
            # only models are read
            for inc in range(len(project_info['MODELS'])):
        
                model=project_info['MODELS'][inc]
                
                # only for non-reference models
                if not model.model_line.split()[1] == project_info['RUNTIME']['currDiag'].variables[v].ref_model:
                
                    model_filename=model_filelist[model.model_line.split()[1]]
                    reference_filename=model_filelist[project_info['RUNTIME']['currDiag'].variables[v].ref_model]
                
                    # copy old data to provide data that is needed again
                    D_old=copy(Diag)
                    
                    # initialize diagnostic
                    Diag = AlbedoDiagnostic()

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
    
    print('>>>>>>>> ENDED SUCESSFULLY!! <<<<<<<<<<<<')
    print('')

if __name__ == "__main__":
    main()
