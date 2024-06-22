import iris
import os
import re
import glob
from netCDF4 import num2date
from shelve import open as shopen
from itertools import product
import numpy as np

"""
Script to generate the additional datasets for the brightpost analysis.
"""


skip_models  = [      # because
    'AWI-CM-1-1-MR',  # triangular grid
    'AWI-ESM-1-1-LR', # triangular grid
    'IITM-ESM',       # No longitude bounds.
    'FGOALS',         # No longitude bounds in irregular grid.
    'NESM3',          # Branch time in parent is incorrect.
    'MIROC-ES2L',     # SSP data contains loads of other variables.
    ]


def load_exp_files(exp,
        short_names = ['tas', 'o2', 'ph', 'areacello', ],
        mips = ['Omon','Ofx', 'fx', ]
    ):
    """
    
    """
    shelve_fn = 'shelves/glob_test_'+exp+'.shelve'
    exists = [os.path.exists(fn) for fn in glob.glob(shelve_fn+'*')]

    if len(exists) == 3:
        sh = shopen(shelve_fn)
        file_path_dicts = sh['file_path_dicts']
        sh.close()
        return file_path_dicts

    file_path_dicts = {}
    for mip in mips:
        for short_name in short_names:
            glob_path = "/badc/cmip6/data/CMIP6/*/*/*/"+exp+"/r*/"+mip+"/"+short_name+"/gn/latest/*.nc"
            # /badc/cmip6/data/CMIP6/CMIP/MOHC/UKESM1-0-LL/piControl/r1i1p1f2/Ofx/areacello/gn/latest/
            files = glob.glob(glob_path)
            print("Globbing:", glob_path)
            print("Found:", files)
            for fn in files:
            
                paths = fn.split('/')
                index = tuple(paths[1:-1]) # skip start and final path
                if index in file_path_dicts:
                    file_path_dicts[index].append(fn)
                else:
                    file_path_dicts[index] = [fn, ]
    if len(file_path_dicts.keys())>1:
        sh = shopen(shelve_fn)
        sh['file_path_dicts'] = file_path_dicts
        sh.close()
        
    return file_path_dicts


def find_area_index(index, a_files, areaname = 'areacello'):
    """
    Figurge out the area index.
    """
    (badc, cmip6, data, CMIP6, CMIP, institute, model, exp, ensemble, mip, short_name, gn, latest) = index
    # Create a list of possible candidates.
    area_indices = []
    enss = ['r1i1p1f1', 'r1i1p1f2', 'r1i1p2f1', 'r1i1p1f3']
    enss.append(ensemble)
    for pien,damip, ofx, nexp  in product(enss, ['DAMIP', 'CMIP', 'ScenarioMIP'], ['Ofx', 'fx','Omon', ], [exp, 'piControl', 'historical']):
        area_indices.append((badc, cmip6, data, CMIP6, damip, institute, model,  nexp, pien, ofx, areaname, gn, latest))

    area_index = []
    areafn = None
    # iterqatre through the file of candidates and choose the first existing file:
    for area_index in area_indices:
        print("find_area_index:", '/' + '/'.join(area_index))
        if a_files.get(area_index, False):
            areafn = a_files[area_index]
            break
    return area_index, areafn


def convert_index_to_yml(index):
    """
    Takes a dataset index and outputs a string for the yaml.

    """
    #ot_it = (index[6], index[7])
    #f got_pair.get(got_it, False):
    #   print("printintg recipe, already got this one:", got_it)
    #   if single_members: continue
    #lse:
    #   got_pair[got_it] = True
    #f not index in complete_sets_hist_pi: continue
    #f index in missing_time:
    #   exclusion_reasons[(model, exp, ensemble)] = ''.join(['Missing some years.'])
    #   continue

    #count+=1
    # print('index:', count, '/'.join(index))
    #rint('        # dataset: ', index[6], index[8])
    #0     1      2     3     4     5          6      7    8         9    10          11    12 
    (badc, cmip6, data, proj, CMIP, institute, model, exp, ensemble, mip, short_name, grid, latest) = index

    if exp  in ['historical', 'hist-nat']:
        start_year='1990'
        end_year='2000'
    else:
        start_year='2040'
        end_year='2050'

    if short_name in ['areacello',]:
        start_year = '0'
        end_year = '1'

    txt = '\n            - {'
    txt += 'project: '+ proj  
    txt += ', mip: '+ mip   
    txt += ', grid: '+ grid   
    txt += ', institute: '+ institute
    txt += ', dataset: '+ model 
    txt += ', exp: '+ exp   
    txt += ', ensemble: '+ ensemble
    txt += ', start_year: '+ start_year
    txt += ', end_year: '+ end_year
    txt += '}'
    return txt


def main():
    """
    Load the data.
    """
    # load a dict of all files.
    all_files = {}
    short_names = ['areacello', 'tos', 'o2', 'ph', ]
    for exp in ['historical', 'ssp119', 'ssp126', 'ssp245', 'ssp370', 'ssp585']:
        all_files.update(load_exp_files(exp, short_names= short_names))
    all_files.update(load_exp_files('piControl', short_names = ['areacello', ]))


    # To fit out method, each model emsemble pair needs:
    #   tas, o2, ph
    # and each model needs one copy of both:  areacello
    # Unlike before, we don't care about pi control.
    final_datasets = {}
    exclusion_reasons = {}

    for index, hist_file in  all_files.items():
        (badc, cmip6, data, CMIP6, CMIP, institute, model, exp, ensemble, mip, short_name, gn, latest) = index
        if short_name != 'tos': continue
        
        if model in skip_models: continue

        # Look for O2:
        os_index = (badc, cmip6, data, CMIP6, CMIP, institute, model, exp, ensemble, 'Omon', 'o2', gn, latest)
        if not all_files.get(os_index, False):
            exclusion_reasons[(model, exp, ensemble)] = ' '.join(['No', exp, 'o2 data'])
            print('exclusion:', (model, exp, ensemble), 'no o2 data')
            continue 

        # Look for ph:
        ph_index = (badc, cmip6, data, CMIP6, CMIP, institute, model, exp, ensemble, 'Omon', 'ph', gn, latest)
        if not all_files.get(ph_index, False):
            exclusion_reasons[(model, exp, ensemble)] = ' '.join(['No', exp, 'ph data'])
            print('exclusion:', (model, exp, ensemble), 'no ph data')
            continue

        # look for  areacello
        areao_index, areaofn = find_area_index(index, all_files, areaname = 'areacello')
        if None in [areaofn,] :
            exclusion_reasons[(model, exp, ensemble)] = ''.join(['No areacello'])
            print('exclusion:', (model, exp, ensemble), 'no areacello')
            continue

        if int(re.split('(\d+)', ensemble)[1]) > 8:  
            exclusion_reasons[(model, exp, ensemble)] = ''.join(['too many ensembles'])
            print('exclusion:', (model, exp, ensemble), 'too many members')
            continue

        print('Found all five datasets:', index)
        final_datasets[index] = True
        final_datasets[os_index] = True
        final_datasets[ph_index] = True
        final_datasets[areao_index] = True

    print('\n\n')
    for index in sorted(exclusion_reasons.keys()): print(index, exclusion_reasons[index])

    print('\n\n')
    final_yml_txt = ''
    for short_name in short_names:
        final_yml_txt +='\n\n        short_name: '+ short_name+'\n'
        for index in sorted(final_datasets.keys()): 
            if short_name not in index: continue

            final_yml_txt += convert_index_to_yml(index) 

    print(final_yml_txt)

    outfn = 'shelves/output.yml'
    wrriting = open(outfn, 'w')
    wrriting.write(final_yml_txt)
    wrriting.close()
    print('saved as', outfn)
    return


main()
"""
# Look for complete sets of historical data/
# Complete set includes temperature, salinity and volume.

complete_sets_hists = []

linked_datasets = {}

pi_ens = ['r1i1p1f1', 'r1i1p1f2', 'r1i1p2f1', 'r1i1p1f3']
for index, hist_file in hist_files.items():
    (badc, cmip6, data, CMIP6, CMIP, institute, model, exp, ensemble, mip, short_name, gn, latest) = index
    if model in ['BCC-CSM2-MR', 'CNRM-CM6-1', 'IPSL-CM6A-LR', 'HadGEM3-GC31-LL', 'GISS-E2-1-G', 'GFDL-ESM4']: 
        print('found missing model', model)
    if short_name != 'thetao': continue
    so_index = (badc, cmip6, data, CMIP6, CMIP, institute, model, exp, ensemble, mip, 'so', gn, latest)
    if not hist_files.get(so_index, False): 
        exclusion_reasons[(model, exp, ensemble)] = ''.join(['No ', exp, 'salinity data'])
        continue
    #pi_exp = guess_PI_ensemble(dicts, keys, ens_pos = 8)

    vol_indices = []
    vol_indices.append((badc, cmip6, data, CMIP6, CMIP, institute, model, exp, ensemble, mip, 'volcello', gn, latest))
    for pien in pi_ens:
        vol_indices.append((badc, cmip6, data, CMIP6, CMIP, institute, model, exp, pien, 'Ofx', 'volcello', gn, latest))
        vol_indices.append((badc, cmip6, data, CMIP6, 'CMIP', institute, model, 'piControl', pien, 'Ofx', 'volcello', gn, latest))
        vol_indices.append((badc, cmip6, data, CMIP6, 'CMIP', institute, model, 'piControl', pien, 'fx', 'volcello', gn, latest))
        vol_indices.append((badc, cmip6, data, CMIP6, CMIP, institute, model, exp, pien,  'fx', 'volcello', gn, latest))
    vol_indices.append((badc, cmip6, data, CMIP6, CMIP, institute, model, exp, ensemble, mip, 'thkcello', gn, latest))
    #print(vol_indices)
    #assert 0
    volfn = []
    out__vol_index = []
    for vol_index in vol_indices:
        if hist_files.get(vol_index, False):
            volfn = hist_files[vol_index]
            out__vol_index = vol_index
            break
        if pi_files.get(vol_index, False):
            volfn = pi_files[vol_index]
            out__vol_index = vol_index
            break            
    if not volfn: 
        exclusion_reasons[(model, exp, ensemble)] = ''.join(['No ', exp, ' or Ofx volume data'])
        continue
    if 'thkcello' in out__vol_index:
        print("WARNING: this one doesn't have volume but does have thkcello")

        area_index, areafn = find_area_index(index, all_files)
        if not areafn:
            print("Found thkcello but not area:", out__vol_index)
            exclusion_reasons[(model, exp, ensemble)] = ''.join(["Found thkcello but not area:"])
            continue
    print("Found one!", index)
    complete_sets_hists.append(index)
    linked_datasets[index] = {}
    linked_datasets[index]['thetao_hist'] = index
    linked_datasets[index]['so_hist'] = so_index
    linked_datasets[index]['vol_hist'] = out__vol_index


print('\ncomplete_sets_hists\n:',complete_sets_hists)

complete_sets_hist_pi = []     
no_volumes = []


        
for index in sorted(complete_sets_hists):
    (badc, cmip6, data, CMIP6, CMIP, institute, model, exp, ensemble, mip, short_name, gn, latest) = index
    if model in ['BCC-CSM2-MR', 'CNRM-CM6-1', 'IPSL-CM6A-LR', 'HadGEM3-GC31-LL', 'GISS-E2-1-G', 'GFDL-ESM4']:
        print('found missing model', model)
#        assert 0
    piindex_ts = []
    piindex_ss = []
    for pien in pi_ens:
        piindex_ts.append((badc, cmip6, data, CMIP6, 'CMIP', institute, model, 'piControl', pien, 'Omon', 'thetao', gn, latest))
        piindex_ss.append((badc, cmip6, data, CMIP6, 'CMIP', institute, model, 'piControl', pien, 'Omon', 'so', gn, latest))

    pi_tmp = False
    pi_so = False
    for piindex_t in piindex_ts:
        pi_tmp = pi_files.get(piindex_t, False)
        if pi_tmp is not False: break   


    for piindex_s in piindex_ss:
        pi_so = pi_files.get(piindex_s, False)
        if pi_so is not False: break

    if False in [pi_tmp,]: 
        exclusion_reasons[(model, exp, ensemble)] = ''.join(['No pi Control Temperature data'])
        continue
    if False in [pi_so, ]:
        exclusion_reasons[(model, exp, ensemble)] = ''.join(['No pi Control salinity data'])
        continue

    pivol_indices = []
    for pien in pi_ens:
        pivol_indices.append((badc, cmip6, data, CMIP6, 'CMIP', institute, model, 'piControl', pien, 'Omon', 'volcello', gn, latest))
        pivol_indices.append((badc, cmip6, data, CMIP6, 'CMIP', institute, model, 'piControl', pien, 'Ofx', 'volcello', gn, latest))
        pivol_indices.append((badc, cmip6, data, CMIP6, CMIP, institute, model, exp         , pien, 'Ofx', 'volcello', gn, latest))
        pivol_indices.append((badc, cmip6, data, CMIP6, 'CMIP', institute, model, 'piControl', pien, 'fx', 'volcello', gn, latest))    
        pivol_indices.append((badc, cmip6, data, CMIP6, CMIP, institute, model,  exp        , pien, 'fx', 'volcello', gn, latest))
        pivol_indices.append((badc, cmip6, data, CMIP6, 'CMIP', institute, model, 'piControl', pien, 'Omon', 'thkcello', gn, latest))

    volfn = False
    for vol_index in pivol_indices:
        if pi_files.get(vol_index, False):
            volfn = pi_files[vol_index]
            break    
        if hist_files.get(vol_index, False):
            volfn = hist_files[vol_index]
            break                
    if not volfn: 
        exclusion_reasons[(model, exp, ensemble)] = ''.join(['No pi Control volume data'])
        continue
    if 'thkcello' in vol_index:
        print("WARNING: this one doesn't have volume but does have thkcello")
        area_index, areafn = find_area_index(index, all_files)

        if not areafn:
            print("Found thkcello but not area:", vol_index)
            exclusion_reasons[(model, exp, ensemble)] = ''.join(["Found thkcello but not area:"])
            continue    
        no_volumes.append(index)
    print("Found all files to run this analysis!", index)
    complete_sets_hist_pi.append(index)

    linked_datasets[index]['thetao_pi'] = piindex_t
    linked_datasets[index]['so_pi'] = piindex_s
    linked_datasets[index]['vol_pi'] = vol_index

#####
# Searching for the areacello file:
for index in linked_datasets.keys():
    (badc, cmip6, data, CMIP6, CMIP, institute, model, exp, ensemble, mip, short_name, gn, latest) = index

    area_index, areafn = find_area_index(index, all_files)

    if not areafn:
        exclusion_reasons[(model, exp, ensemble)] = ''.join(["No areacello"])
        continue
    linked_datasets[index]['areacello'] = area_index 

print('\n\n')
for nv in no_volumes:print('No Volume (thickcello only)', nv)

print('\ncomplete_sets_hist_pi\n:',complete_sets_hist_pi)

####
# Convert nc file path into list of years.
def calc_available_years(a_files):
    available_years = {}
    for index, files in a_files.items():
        #if not index in complete_sets_hist_pi: continue
        #print('/'.join(index), files)
        for fn in files:
            basename = fn.split('/')[-1]
            basename = basename.split('_')[-1]
            basename = basename.split('.')[0]
            years =  basename.split('-')
            print(years, basename, fn)
            # working here!
            if years == ['gn']: # Ofx file, no times!
                available_years[index] = None 
                # assert 0
                continue
            start = int(years[0][:-2])
            end = int(years[1][:-2]) + 1 
            if index in available_years:
                available_years[index] += list(range(start, end, 1))
            else:
                available_years[index] = list(range(start, end, 1))
        print('/'.join(index), files, available_years[index])
    return available_years 

Available_years = calc_available_years(all_files)
required_years = set(range(1860,2014,1))
missing_time = []

for index in linked_datasets.keys():
    if not index in complete_sets_hist_pi: 
        continue
    for key, ind in linked_datasets[index].items():
        if 'Ofx' in ind: continue
        if 'fx' in ind: continue
        print(key, ind)
        if Available_years.get(ind, None) == None: continue
        if ind[7] in [ 'historical', 'hist-nat']:
            years = set(Available_years[ind])
            intersect = years.intersection(required_years)
            if len(intersect) == len(required_years):
                print("no gaps!")
            else:
                print("time gaps!", '/'.join(index))
                missing_time.append(index)
                missing_time.append(ind)

####
# Calculate pi control years.
def calc_time(hist_fn, historical_range = [1860, 2014]): #, 2100]):
    cube= iris.load_cube(hist_fn)

    times = cube.coord('time')
    units = times.units.name
    calendar = times.units.calendar
    # We can assume a shared calendar!

    parent_branch_yr = num2date(cube.attributes['branch_time_in_parent'],
                                units=cube.attributes['parent_time_units'],
                                calendar=calendar ).year

    child_branch_yr = num2date(cube.attributes['branch_time_in_child'],
                               units=units, calendar=calendar ).year # A

    diff = child_branch_yr - parent_branch_yr

    #historical_range = [1860, 2014, 2100]
    print('date in parent:\t', historical_range, 'is', [h-diff for h in historical_range])
    return [h-diff for h in historical_range]

# Calculate PI control years.
pi_years = {}
pi_years_single = {}

for index in linked_datasets.keys():
    if not index in complete_sets_hist_pi:
        continue
    if index in missing_time:
        (badc, cmip6, data, CMIP6, CMIP, institute, model, exp, ensemble, mip, short_name, gn, latest) = index
        exclusion_reasons[(model, exp, ensemble)] = ''.join(['Missing some years.'])
        continue
    # Should only be hist thetao file index.
    hist_thetao_fn = all_files[index][0]
 
    linked_datasets[index]['thetao_hist'] = index
    pi_years[index] = calc_time(hist_thetao_fn, historical_range = [1860, 2014])



#### 
# Convert to recipe
count = 0

# Flag to set a single ensemble member only.
single_members = False

# tool to determinie the pi control time range for all model ensemble members.
if not single_members:
    for index, times in pi_years.items():
        index2= list(index)
        index2[8] = 'all' # 8 in ensemble member (r1i1p1f1)
        index2 = tuple(index2)
        old_times = pi_years_single.get(index2, times)
        pi_years_single[index2] = [np.min([old_times[0], times[0]]), np.max([old_times[1], times[1]])]
        print(index, times,'->', pi_years_single[index2])

    for index, times in pi_years.items():
        index2= list(index)
        index2[8] = 'all' # 8 in ensemble member (r1i1p1f1)
        index2 = tuple(index2)
        pi_years_single[index] = pi_years_single[index2]

got_pair = {}
recipe_txt = {'thetao':{}, 'so':{}, 'vol_Omon':{}, 'vol_fx':{}, 'areacello':{}, 'thkcello':{}}
for index in sorted(linked_datasets.keys()):
    got_it = (index[6], index[7])
    if got_pair.get(got_it, False):
        print("printintg recipe, already got this one:", got_it)
        if single_members: continue
    else:
        got_pair[got_it] = True
    if not index in complete_sets_hist_pi: continue
    if index in missing_time:
        exclusion_reasons[(model, exp, ensemble)] = ''.join(['Missing some years.'])
        continue

    count+=1
    # print('index:', count, '/'.join(index))
    print('        # dataset: ', index[6], index[8])
    for key, ind in linked_datasets[index].items():
        ofx_bool = ind[9] in ['Ofx', 'fx']
        if ind[7] in ['historical', 'hist-nat']:
            start_year=1850
            end_year=2014
        if ind[7] == 'piControl' and not ofx_bool:
            if single_members:
                start_year = pi_years[index][0]
                end_year = pi_years[index][1]
            else:
                start_year = pi_years_single[index][0]
                end_year = pi_years_single[index][1]

        if ofx_bool:
            start_year = 0
            end_year=1

        txt = '\n            - {'
        txt += 'project: '+ ind[3]
        txt += ', mip: '+ ind[9]
        txt += ', grid: '+ ind[11]
        txt += ', dataset: '+ ind[6]
        txt += ', exp: '+ ind[7]
        txt += ', ensemble: '+ ind[8]
        txt += ', start_year: '+ str(start_year)
        txt += ', end_year: '+ str(end_year)
        txt += '}'

        if 'thkcello' in ind:
             recipe_txt['thkcello'][txt] = True
             continue
        if ind[10] in ['thetao', 'so']:
            recipe_txt[ind[10]][txt] = True
        if ind[10] == 'volcello' and ind[9] == 'Omon' :
             recipe_txt['vol_Omon'][txt] = True
        if ind[10] == 'volcello' and ind[9] in ['Ofx', 'fx']:
            recipe_txt['vol_fx'][txt] = True
        if ind[10] == 'areacello':
            recipe_txt['areacello'][txt] = True


out_txt = ''
for key, truedict in recipe_txt.items():
    string = ''.join([txt for txt in sorted(truedict.keys())])
    print(key, string)
    out_txt+=key
    out_txt+=string

wrriting = open('shelves/output.yml', 'w')
wrriting.write(out_txt)
wrriting.close()
          # - {project: CMIP6, mip: Omon, grid: gn, dataset: ACCESS-ESM1-5,     exp: historical, ensemble: r1i1p1f1, start_year: 1860, end_year: 2014}



for (model, exp, ensemble) in sorted(exclusion_reasons.keys()):
    print(model, exp, ensemble, ':', exclusion_reasons[(model, exp, ensemble)])




"""
