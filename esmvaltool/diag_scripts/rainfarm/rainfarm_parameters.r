#####################################################################
#
#  RainFARM advanced configuration file
#
#  About: configuration file with advanced settings for ESMValTool RainFARM namelist.
#
#####################################################################

force_processing = T # TRUE to rewrite existing output file

# RainFARM options 
# Note - RainFARM usage:  rfarm [-s SLOPE] [-e NENS] [-n NF] [-r REGION] [-w WEIGHTS]
#                               [-o OUTFILE] [-v VARNAME] [-g] [-c]
# set FALSE to void argument
weights_climo =     # climatology file to be used to calculate weights (calling rfweights) 
          "/work/datasets/obs/unsorted/WORLDCLIM/prec/prec_06_16_central_europe.nc" 
varname = F         # input variable name - this is assigned by ESMValTool pre-processing 
conserv_glob = F    # conserve precipitation over full domain
conserv_smooth = T  # conserve precipitation using convolution

# NOTE: for more information type 'rfarm -h' from terminal
