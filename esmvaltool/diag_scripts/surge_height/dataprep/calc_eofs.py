import os
import pickle

from eofs.standard import Eof
from eofs.xarray import Eof as xrEof


def save_object(filename, obj):
	with open(filename, 'wb') as output:
		pickle.dump(obj,output,pickle.HIGHEST_PROTOCOL)


def calc_eofs(psl, uas, vas, gradslp, eof_savepath):
    # determine eofs
    psl_solver = xrEof(psl)
    fout_psl = os.path.join(eof_savepath,
       'eof_psl_solver_075x075_19790101-20151231.pkl'
    )
    save_object(fout_psl, psl_solver)


    gradslp_solver = Eof(gradslp)
    fout_lon = os.path.join(eof_savepath,
       'eof_gradslp_solver_075x075_19790101-20151231.pkl'
    )
    save_object(fout_lon, gradlat_solver)


    uas_solver = xrEof(uas)
    fout_uas = os.path.join(eof_savepath,
       'eof_uas_solver_075x075_19790101-20151231.pkl'
    )
    save_object(fout_uas, uas_solver)


    vas_solver = xrEof(vas)
    fout_vas = os.path.join(eof_savepath,
       'eof_vas_solver_075x075_19790101-20151231.pkl'
    )
    save_object(fout_vas, vas_solver) 

    return psl_solver, gradlon_solver, gradlat_solver, uas_solver, vas_solver
