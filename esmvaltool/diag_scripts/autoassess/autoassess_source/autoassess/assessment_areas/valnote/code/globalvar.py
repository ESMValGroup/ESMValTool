'''
This is a python config file used to store global variables used throughout
any code that uses:

make_plots.py, valmod.py, rms.py

http://dev.fyicenter.com/Interview-Questions/Python/How_do_I_share_global_variables_across_modules_.html
'''

# The defaults are:

# ************* FLAGS **************

error = False           # Flag to say that an error has occured somewhere
debug = False           # Flag to say we are using debug mode
pub = False             # Flag to say this is for a publication. Don't output
                        #  UMUI jobids or Rose suiteids in the final plot
                        #  headers.
valnote = False         # Flag to say we are doing a validation note
pkl = False             # Flag to use pkl files to speed up reruns
store_regrids=True      # Flag to store regridded data for future use.
                        #  Speeds up the validation note but uses more memory.

# ************ GENERAL VARIABLES *************

plot_type = 'lat_lon'   # What type of plot are we plotting
rms_first_letter = 'a'  # What is the first plot letter (key) to calculate
                        #  RMS values for
extra_data_dict = {}    # A dictionary containing extra data of land sea
                        #  masks and grids at different resolutions
missing_data = -10000   # A missing data indicator number

# *********** VARIABLES USED BY VALIDATION NOTES *****************

item_dict = {}          # A dictionary referencing items to their stash codes
                        #  and pre-processing options
levels_dict = {}        # A dictionary referencing plots to their levels,
                        #  colours and post-processing options
plots_dict = {}         # A dictionary referencing what equations are used
                        #  to make each plot
valorder_file = 'none'  # Full path to the valorder file used.
