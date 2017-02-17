;;#############################################################################
;; carbon_constraint.ncl
;; Author: Sabrina Wenzel (DLR, Germany)    
;; EMBRACE project 
;;#############################################################################
;;                                                                            
;; Description:
;;    This is a dummy diagnostic to process data with the ESMValTool and save  
;;    in the climo directory
;;    
;;
;; Modification history
;;    * 20151027-A_laue_ax: added call to 'write_references'
;;    * 201510-A_Sabrina Wenzel: written 
;; #############################################################################
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"
load "./interface_scripts/constants.ncl"

begin
  diag_script = "carbon_dummy"
  
  verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
  info_output("<<<<<<<< Entering " + diag_script + ".ncl", verbosity, 4)
; A-laue_ax+
    ; ###########################################
    ; # references                              #
    ; ###########################################

    write_references(diag_script,      \  ; script name
                     "A_wenz_sa",      \  ; authors
                     "",               \  ; contributors
                     "",               \  ; diag_references
                     "",               \  ; obs_references
                     (/""/))              ; proj_references
; A-laue_ax-
  info_output("==> processed data using: " + diag_script + ".ncl", verbosity, 1)
  info_output(">>>>>>>> Leaving " + diag_script + ".ncl", verbosity, 4)

end 
;#################################################################################
