import os
import pickle


def load_eofs(eof_inpath):
    #global SLPsolver, gradlonsolver, gradlatsolver, usolver, vsolver
    #
    with open(os.path.join(eof_inpath,
            'eof_psl_solver_075x075_19790101-20151231.pkl'),
            'rb') as input:
        SLPsolver = pickle.load(input)
    #
    with open(os.path.join(eof_inpath,
            'eof_gradpsl_solver_075x075_19790101-20151231.pkl'),
            'rb') as input:
        gradslpsolver = pickle.load(input)
    #
    with open(os.path.join(eof_inpath,
            'eof_uas_solver_075x075_19790101-20151231.pkl'),
            'rb') as input:
        usolver = pickle.load(input)
    #
    with open(os.path.join(eof_inpath,
            'eof_vas_solver_075x075_19790101-20151231.pkl'),
            'rb') as input:
        vsolver = pickle.load(input)

    return  SLPsolver, gradslpsolver, usolver, vsolver
