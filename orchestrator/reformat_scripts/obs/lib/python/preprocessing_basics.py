# -*- coding: utf-8 -*-
"""
These are (cdo) wrapper functions that can be used by the obs reformatting 
scripts in python.

Most functions will write temporary files in a subdirectory "temp" that has to 
be written beforehand and deleted afterwards 
"""

import os, glob, tempfile, math, subprocess
from cdo import Cdo

        
def _get_files_in_directory(directory, pattern, asstring=True):
    """ returns list and number of files with pattern in directory """
    
    if directory[-1] != os.sep:
        directory += os.sep
        
    L = glob.glob(directory + pattern)
    N = len(L)
    if asstring:
        L = ' '.join(L)
    return L, N
    
def _get_subdirectories_in_directory(directory, asstring=True):
    """ returns list and number of directories in directory """
    
    if directory[-1] != os.sep:
        directory += os.sep
        
    L = glob.glob(directory)
    N = len(L)
    if asstring:
        L = ' '.join(L)
    return L, N
    
def _aggregate_obs_from_files(work_dir,file_list):
    """ reads different obs files into one file """
    cdo=Cdo()
    oname=work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1]
    try:
        cdo.cat(input=file_list,output=oname,options = '-f nc4')
        return oname
    except OSError:
        fl=len(file_list)
        if fl>1:
            fhalf=file_list[0:int(math.floor(fl/2))]
            shalf=file_list[int(math.floor(fl/2)):fl]
            foname=_aggregate_obs_from_files(work_dir,fhalf)
            soname=_aggregate_obs_from_files(work_dir,shalf)
            oname=_aggregate_obs_from_files(work_dir,[foname,soname])
            os.remove(foname)
            os.remove(soname)
            return oname
        else:
            print "package too small"
    
def _aggregate_timestep(work_dir,infile,timestep,remove=True):
    """ aggregate infile to timestep """
    """ currenty only monthly """
    cdo=Cdo()
    oname=work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1]
    if timestep=="monthly":
        cdo.monmean(input=infile,output=oname,options='-f nc4 -b F32')
        if remove:
            os.remove(infile)
    else:
        assert False, "This timestep cannot be handled yet."
        
    return oname
    
def _aggregate_specific_years(work_dir,infile,times,remove=True):
    """ aggregate infile to times with mean and sd"""
    
    cdo=Cdo()
    onameM=work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1]
    onameS=work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1]
    oname=work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1]
    tmpname=work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1]
    cdo.selyear(",".join([str(t) for t in times]),input=infile,output=tmpname,options='-f nc4 -b F32')
    cdo.timselmean(12,input=tmpname,output=onameM,options='-f nc4 -b F32')
    #cdo.timselstd(12,input=tmpname,output=onameS,options='-f nc4 -b F32')
    name=cdo.showname(input=onameM)
    cdo.setname(name[0] + "_std -timselstd,12",input=tmpname,output=onameS,options='-L -f nc4 -b F32')
    cdo.merge(input=[onameM,onameS],output=oname)
    if remove:
        os.remove(infile)
    os.remove(tmpname)
    os.remove(onameM)
    os.remove(onameS)
        
    return oname
    
def _aggregate_resolution(work_dir,infile,resolution,remove=True):
    """ aggregate infile to resolution """
    """ currenty only T63, T85 """
    cdo=Cdo()
    oname=work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1]
    if resolution=="T63":
        cdo.remapcon('t63grid',input=infile,output=oname,options='-f nc4 -b F32')
    elif resolution=="T85":
        cdo.remapcon('t85grid',input=infile,output=oname,options='-f nc4 -b F32')
    else:
        assert False, "This resolution cannot be handled yet."
        
    if remove:
        os.remove(infile)   
        
    return oname
    
def _select_variable(work_dir,infile,variablename,remove=False):
    """ select variables from infile """
    cdo=Cdo()
    oname=work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1]
    cdo.selname(variablename,input=infile,output=oname,options='-f nc4 -b F32')
    if remove:
        os.remove(infile)
    return oname
    
def _sum_files(work_dir,infiles,remove=True):
    """ sum up all infiles """
    cdo=Cdo()
    oname=work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1]
    cdo.enssum(input=" ".join(infiles),output=oname,options='-f nc4 -b F32')
    if remove:
        for ifi in infiles:
            os.remove(ifi)
    return oname
    
def _extract_variables(work_dir,infile,variablenames,newvarname,remove=True):
    """ select, sum up and rename variable(s) from infile """
    
    tmpname=work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1]
    oname=work_dir + os.sep + "temp" + os.sep + tempfile.NamedTemporaryFile().name.split('/')[-1]
    selfilenames=[]
    
    stop=0
    
    varnames=list(variablenames)
    
    while len(varnames)>stop:
        
        v=varnames[0]
        thisfile=_select_variable(infile,v)
        selfilenames.insert(0,thisfile)
        varnames.remove(v)
        
    thisfile=_sum_files(selfilenames)
    
    # adjust variable name and time axis
    cdo=Cdo()
    cdo.setname('_'.join(newvarname.split(" ")),input=thisfile,output=tmpname,options='-f nc4 -b F32')
    cdo.settaxis(infile.split("-")[-2]+"-01-01","00:00",input=tmpname,output=oname,options='-f nc4 -b F32')
    subprocess.call(['ncatted', '-O', '-a', 'units,' + '_'.join(newvarname.split(" ")) + ',c,c,%', oname])
    
    if remove:
        os.remove(thisfile)
        os.remove(tmpname)
    return oname
    
