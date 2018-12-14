import os
import pickle

import numpy as np


def check_solver(eof_inpath):
    psl_solver = os.path.join(eof_inpath,
      'eof_psl_solver_075x075_19790101-20151231.pkl')
    gradpsl_solver = os.path.join(eof_inpath, 
      'eof_gradpsl_solver_075x075_19790101-20151231.pkl')
    uas_solver = os.path.join(eof_inpath, 
      'eof_uas_solver_075x075_19790101-20151231.pkl')
    vas_solver = os.path.join(eof_inpath, 
      'eof_vas_solver_075x075_19790101-20151231.pkl')
    if (os.path.exists(psl_solver) and os.path.exists(gradpsl_solver) and
       os.path.exists(uas_solver) and os.path.exists(vas_solver)):
        for file in [psl_solver, gradpsl_solver, 
               uas_solver, vas_solver 
               ]:
            try:
                with open(file,'rb') as fp:
                    test = pickle.load(fp)
                    solvers_exist = True
            except:
                solvers_exist = False
    else:
        solvers_exist = False

    return solvers_exist
