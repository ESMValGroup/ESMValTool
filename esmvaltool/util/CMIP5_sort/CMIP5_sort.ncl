;;##########################################################
;; Description:
;; This Skript sorts the CMIP5 files into subdirectories by 
;; using the substrings of the filename:
;;
;;	sub-strings: 0 = variable name
;;		     1 = realm
;;		     2 = model name
;;		     3 = experiment
;;		     4 = ensemble
;;
;; Modification history
;;    * 201501-A_wenz_sa: written   
;;##########################################################
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl" 

; user parameter: set output path "outpath"
; input files are expected in the current directory "."

  outpath = getenv(ESMValTool_MODELPATH) + "/ETHZ_CMIP5"
  fils = systemfunc ("ls *.nc")				
  dim_FILE = dimsizes(fils)
  
  do ff = 0,dim_FILE-1
     f_string = str_split(fils(ff),"_") 
     	    
     ;;Folder for experiment
     if (.not.isfilepresent(f_string(3))) then 
      system("mkdir " + outpath + "/" + f_string(3))	 
     end if
     f_path = outpath + "/" + f_string(3)+"/"
     ;print(f_path)

     ;;Folder for realm
     if (.not.isfilepresent(f_path+f_string(1))) then 
      system("mkdir " + f_path+f_string(1))	 
     end if
     f_path = f_path+f_string(1)+"/"
     ;print(f_path)

     ;;Folder for variable
     if (.not.isfilepresent(f_path+f_string(0))) then 
      system("mkdir " + f_path+f_string(0))	 
     end if
     f_path = f_path+f_string(0)+"/"
     ;print(f_path)

     ;;Folder for variable
     if (.not.isfilepresent(f_path+f_string(2))) then 
      system("mkdir " + f_path+f_string(2))	 
     end if
     f_path = f_path+f_string(2)+"/"
     ;print(f_path)

     ;;Folder for ensemble
     if (.not.isfilepresent(f_path+f_string(4))) then 
      system("mkdir " + f_path+f_string(4))	 
     end if
     f_path = f_path+f_string(4)+"/"

     print(f_path+" "+fils(ff))
     system("mv " + fils(ff) + " " + f_path)

     delete(f_string)
     delete(f_path)
  end do	
  delete(fils)
  delete(dim_FILE)
