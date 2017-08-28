"""
Checks different mandatory fields that need to be passed in for analysis
Author: Valeriu Predoi, University of Reading, valeriu.predoi@ncas.ac.uk
First version: August 2017
"""
import sys

def models_checks(models_dict):
    """
    checking for standard keys and values e.g. for CMIP5*:
    {name: MPI-ESM-LR, project: CMIP5,  mip: Amon,  exp: historical,  ensemble: r1i1p1,  
     start_year: 2000,  end_year: 2002,  path: ./test_data}
    Fastest check is try/except
    """
    # models loop
    for m in models_dict:
        try:
            name = m['name']
            if name is None:
                print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR model name is None', m
                sys.exit(1)
        except KeyError as e:
            print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in model', m
            sys.exit(1)
        try:
            project = m['project']
            if project is None:
                print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR model project is None', m
                sys.exit(1)
        except KeyError as e:
            print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in model', m
            sys.exit(1)
        ####### check specific project name cases ########################################
        ####### CMIP5 and any other CMIP5 derivatives ####################################
        if project.startswith('CMIP5') is True:
            try:
                mip = m['mip']
                if mip is None:
                    print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR model mip is None', m
                    sys.exit(1)
            except KeyError as e:
                print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in model', m
                sys.exit(1)
            try:
                exp = m['exp']
                if exp is None:
                    print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR model exp is None', m
                    sys.exit(1)
            except KeyError as e:
                print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in model', m
                sys.exit(1)
        ######## EMAC and any other EMAC drivative ############################################
        # same as default
        ######## GFDL and any other GFDL derrivative ##########################################
        if project.startswith('GFDL') is True:
            try:
                realm = m['realm']
                if realm is None:
                    print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR model realm is None', m
                    sys.exit(1)
            except KeyError as e:
                print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in model', m
                sys.exit(1)
            try:
                shifty = m['shift']
                if shifty is None:
                    print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR model shift is None', m
                    sys.exit(1)
            except KeyError as e:
                print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in model', m
                sys.exit(1)
        ######### CCMVal and any CCMVal derrivative ###################################################
        if project.startswith('CCMVal') is True:
            try:
                exp = m['exp']
                if exp is None:
                    print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR model exp is None', m
                    sys.exit(1)
            except KeyError as e:
                print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in model', m
                sys.exit(1)
        ########## end specific project name cases ####################################################
        try:
            ensemble = m['ensemble']
            if ensemble is None:
                print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR model ensemble is None', m
                sys.exit(1)
        except KeyError as e:
            print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in model', m
            sys.exit(1)
        try:
            start_year = m['start_year']
            if start_year is None:
                print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR model start_year is None', m
                sys.exit(1)
        except KeyError as e:
            print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in model', m
            sys.exit(1)
        try:
            end_year = m['end_year']
            if end_year is None:
                print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR model end_year is None', m
                sys.exit(1)
        except KeyError as e:
            print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in model', m
            sys.exit(1)
# FIX-ME not needed anymore, can be deleted
#        try:
#            path = m['path']
#            if path is None:
#                print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR model path is None', m
#                sys.exit(1)
#        except KeyError as e:
#            print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in model', m
#            sys.exit(1)

def diags_checks(diags_dict):
    """
    checking for standard keys and values
    Fastest check is try/except
    """
    for c in diags_dict:
        D = diags_dict[c]
        try:
            variables = D.variables
            for v in variables:
                try:
                    vname = v['name']
                    if vname is None:
                        print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR variable name is None', v, c
                        sys.exit(1)
                except KeyError as e:
                    print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in variable', v, c
                    sys.exit(1)
                try:
                    vfield = v['field']
                    if vfield is None:
                        print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR variable field is None', v, c
                        sys.exit(1)
                except KeyError as e:
                    print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in variable', v, c
                    sys.exit(1)
                try:
                    vpreproc = v['preproc_id']
                    if vpreproc is None:
                        print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR variable preproc_id is None', v, c
                        sys.exit(1)
                except KeyError as e:
                    print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in variable', v, c
                    sys.exit(1)
                try:
                    vref = v['ref_model']
                    if vref is None:
                        print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR variable ref_model is None', v, c
                        sys.exit(1)
                    else:
                        if isinstance(vref, list) is False:
                            print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR variable ref_model needs to be a list', v, c
                            sys.exit(1)
                except KeyError as e:
                    print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in variable', v, c
                    sys.exit(1)
        except AttributeError as e:
            print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, ' in diagnostic', c
            sys.exit(1)
        try:
            scripts = D.scripts
            for s in scripts:
                try:
                    sscript = s['script']
                    if sscript is None:
                        print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR script is None', s, c
                        sys.exit(1)
                except KeyError as e:
                    print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in script', s, c
                    sys.exit(1)
                try:
                    scfg = s['cfg_file']
                    if scfg is None:
                        print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR cfg_file is None', s, c
                        sys.exit(1)
                except KeyError as e:
                    print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in cfg_file', s, c
                    sys.exit(1)
        except AttributeError as e:
            print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, ' in diagnostic', c
            sys.exit(1)

def preprocess_checks(preprocess):
    """
    Takes each dictionary object from PREPROCESS
    and looks for standard keys eg
    {id: pp1, select_level: None, target_grid: ref_model, regrid_scheme: linear,
     mask_fillvalues: True, mask_landocean: None, multimodel_mean: True}
    very similar to models_checks 
    """
    for m in preprocess:
        try:
            idp = m['id']
            if idp is None:
                print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR preprocess id is None', m
                sys.exit(1)
        except KeyError as e:
            print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in preprocess', m
            sys.exit(1)
        try:
            sl = m['select_level']
            if sl is None:
                print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR preprocess select_level is None', m
                sys.exit(1)
        except KeyError as e:
            print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in preprocess', m
            sys.exit(1)
        try:
            tg = m['target_grid']
            if tg is None:
                print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR preprocess target_grid is None', m
                sys.exit(1)
            else:
                if isinstance(tg, basestring) is False:
                    print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR preprocess target_grid needs to be a string: None, ref_model or XxY', m
                    sys.exit(1)
        except KeyError as e:
            print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in preprocess', m
            sys.exit(1)
        try:
            rs = m['regrid_scheme']
            if rs is None:
                print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR preprocess regrid_scheme is None', m
                sys.exit(1)
        except KeyError as e:
            print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in preprocess', m
            sys.exit(1)
        try:
            mf = m['mask_fillvalues']
            if mf is None:
                print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR preprocess mask_fillvalues is None', m
                sys.exit(1)
        except KeyError as e:
            print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in preprocess', m
            sys.exit(1)
        try:
            mlo = m['mask_landocean']
            if mlo is None:
                print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR preprocess mask_landocean is None', m
                sys.exit(1)
        except KeyError as e:
            print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in preprocess', m
            sys.exit(1)
        try:
            mmm = m['multimodel_mean']
            if mmm is None:
                print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR preprocess multimodel_mean is None', m
                sys.exit(1)
        except KeyError as e:
            print >> sys.stderr, 'PY  info:  >>> namelistchecks.py >>> ERROR ', e, 'is missing in preprocess', m
            sys.exit(1)
