#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 14:55:31 2018

@author: Valerio2
"""

import numpy as np
import math
import sys

class srvfile():
    
    def __init__(self ,filename=None, mode = 'r', hbyte = '8', fbyte = '8'):
        
        # open file
        self.fbyte = fbyte
        self.hbyte = hbyte
        self.fileobject = open(filename, mode=mode+'b') 
        # fix service format parameters
        self.define_service_format()
        # store filename
        self.filename = filename
        # read dpstep 
        #self.get_dpstep()
        
        

        
    """
    Define dtype which corresponds to SERVICE FORMAT
    """
    def define_service_format(self):
        self.headersize = int(self.hbyte)*8
        self.read_resolution() # define dim1 and dim2
        print(self.dim1,self.dim2,self.fbyte)
        self.fieldsize = self.dim1*self.dim2*int(self.fbyte) + 8       
        self.headtype       = np.dtype([('headstart','i4') , ('code','i'+self.hbyte ), ('level','i'+self.hbyte ), ('date','i'+self.hbyte ), ('time','i'+self.hbyte ), ('X','i'+self.hbyte ), ('Y','i'+self.hbyte ), ('head7','i'+self.hbyte ), ('expnum','i'+self.hbyte ),('headend','i4')])
        self.fieldtype      = np.dtype("("+str(self.dim2)+","+str(self.dim1)+")f"+self.fbyte )
        self.serviceformat  = np.dtype([('headstart','i4'),
                                  ('code','i'+self.hbyte ), ('level','i'+self.hbyte ), ('date','i'+self.hbyte ), ('time','i'+self.hbyte ), ('X','i'+self.hbyte ), ('Y','i'+self.hbyte ), ('head7','i'+self.hbyte ), ('expnum','i'+self.hbyte ),
                                  ('headend','i4'),
                                  ('fieldstart','i4'),
                                  ('field', "("+str(self.dim2)+","+str(self.dim1)+")f"+self.fbyte ),
                                  ('fieldend','i4')])
        self.recordlength=self.serviceformat.itemsize        
        
    """
    read resolution from header
    """        
    def read_resolution(self):
        
        self.fileobject.seek(0,0)
        self.headtype       = np.dtype([('headstart','i4'),('code','i'+self.hbyte ), ('level','i'+self.hbyte ), ('date','i'+self.hbyte ), ('time','i'+self.hbyte ), ('X','i'+self.hbyte ), ('Y','i'+self.hbyte ), ('head7','i'+self.hbyte ), ('expnum','i'+self.hbyte ),('headend','i4')])
        self.headtype333    = np.dtype([('headstart','i4'),('code','i'+self.hbyte ), ('zero','i'+self.hbyte ), ('date','i'+self.hbyte ), ('time','i'+self.hbyte ), ('nlon','i'+self.hbyte ), ('nlat','i'+self.hbyte ), ('nlev','i'+self.hbyte ), ('ntru','i'+self.hbyte ),('headend','i4')])
        self.firstheader=np.fromfile(self.fileobject,dtype=self.headtype, count = 1)   
        if self.firstheader['code'][0] == 333:
            self.fileobject.seek(0,0)
            self.firstheader=np.fromfile(self.fileobject,dtype=self.headtype333, count = 1)
            self.nlat = self.firstheader['nlat'][0]
            self.nlon = self.firstheader['nlon'][0]
            self.ntru = self.firstheader['ntru'][0]
            self.dtype333 = np.dtype(self.headtype333.descr + [('fieldstart','i4'),
                                 ('zsig',"("+str(self.nlon)+","+str(self.nlat)+")f"+self.fbyte),
                                 ('fieldend','i4')])
            raw = np.fromfile(self.fileobject,dtype=self.dtype333, count = 1)
            self.zsig = raw['zsig'][0]
            self.fileobject.seek(self.dtype333.itemsize, 0)
            self.offset = self.dtype333.itemsize
        else:
            self.offset = 0
        
        self.fileobject.seek(self.offset,0)
        self.firstheader=np.fromfile(self.fileobject,dtype=self.headtype, count = 1)   
        self.fileobject.seek(0,0)
        self.dim1=self.firstheader['X'][0]
        self.dim2=self.firstheader['Y'][0]
        #print("dim1 : ",self.dim1)
        #print("dim2 : ",self.dim2)
        
    """
    get timestep (and nlev)
    """        
    def get_dpstep(self):
        
        self.dt = {0:0}
        record = 1000
        while True:
            record = record + 1
            try: 
                self.read(start=1,count = record)
                break
            except ValueError: 
                pass
            
        try: self.dpstep = self.dt[0][0]
        except: 
            try: self.dpstep = self.dt[0]
            except: print(self.dt)
        if len(self.dt) > 1:
            self.dpstep = np.unique(self.dt)[0]
        
    def read_all(self):
        self.read(start=1,count=-1)

    def read_years(self,startyear,endyear, filtlevels = None):
        
        self.get_dpstep()
        self.startrecord=int((startyear-1)*360/int(self.dpstep))*self.nlev + 1
        self.endrecord=int(endyear*360/int(self.dpstep))*self.nlev
        self.nb_records = (self.endrecord - self.startrecord +1)
        self.read(start=self.startrecord,count=self.nb_records, filtlevels = filtlevels)
        
    """
    read everything (start is first recorded to be read in)
    """
    def read(self,start = 1,count = -1, filtlevels = None):
        if start < 1: raise ValueError("Record numbering starts with 1#")
        # read according to service format definition
        
        
        self.fileobject.seek(self.offset + (start-1)*self.recordlength ,0)
#        rawdata=np.fromfile(self.fileobject,count=count,dtype=self.serviceformat)   
#        rawdata=np.memmap(self.fileobject.name,dtype=self.serviceformat,offset = self.offset + (start-1)*self.recordlength)   
        rawdata=np.fromfile(self.fileobject.name,dtype=self.serviceformat)   
        self.levels = np.unique(rawdata['level'])
        self.codes = np.unique(rawdata['code'])
        self.rawdata = rawdata

        # consistency checks
        
        #self.rawdata = rawdata
        if filtlevels != None:        
            boollevel = np.array([(row['level'] in filtlevels) for row in rawdata])
            rawdata = rawdata[boollevel]
        else:
            self.nlev = len(self.levels)
            filtlevels = self.levels
        
        self.dt = {}
        self.field = {}
        for code in self.codes:
            # check time for all levels
            
            dtdiff = {}
            actuallevels = len(filtlevels)
            #filter for levels
            print('code: ',code,self.codes)
#            for lev in filtlevels:
#                print('level: ',lev)
#                time = np.asarray(list(map(int,rawdata['time'][np.logical_and(rawdata['level']==lev, rawdata['code']==code)])))
#                date = np.asarray(list(map(int,rawdata['date'][np.logical_and(rawdata['level']==lev, rawdata['code']==code)])))
#                
#                print(code,lev,len(time))
#                if len(time) == 0: 
#                    actuallevels = actuallevels - 1
#                    pass
#                elif len(time)==1:
#                    dtdiff[lev] = np.asarray([-1])
#                else:
#                    dtdiff[lev] = np.abs(self.timediff(date,time))
#                if len(np.unique(dtdiff)) > 1:
#                    self.dtdiff = dtdiff
#                    raise ValueError("Time Vector is not consistent for level "+str(lev)+". Different time steps are "+str(np.unique(dtdiff[lev] )))
#            dtcode = np.unique(np.asarray([np.unique(dtdiff[key]) for key in dtdiff.keys() if len(np.unique(dtdiff[key]))>=1]))
#            #print(dtcode,dtdiff)
#            self.dt[code] = dtcode[0]
#            if dtcode.shape[0] != 1:
#                raise ValueError("Time Vector is not consistent for code "+str(code)+". Different time steps are "+str(np.unique(self.dt[code] )))
#                
#            
            # check control bytes of header and field
            
            if (( rawdata['headstart'][code == rawdata['code']]!= self.headersize).any() or
                ( rawdata['headend'][code == rawdata['code']] != self.headersize).any() or 
                ( rawdata['fieldstart'][code == rawdata['code']] != self.dim1*self.dim2*int(self.fbyte )).any() or 
                ( rawdata['fieldend'][code == rawdata['code']] != self.dim1*self.dim2*int(self.fbyte )).any()):
                raise ValueError("Record Control Bytes are not consistent.")
            
            # read in data
            
         
            self.field[str(code)]=(rawdata['field'][code == rawdata['code']]).reshape((actuallevels,int((rawdata[code == rawdata['code']]).shape[0]/actuallevels),self.dim2,self.dim1),order='F')
        
        return rawdata
    
    """
    Compute time difference of time vector
    """    
    def timediff(self,date,time):
                years = np.asarray(list(map(lambda x: np.floor(x/10000),date)))
                date = date - years*10000
                months = np.asarray(list(map(lambda x: np.floor(x/100),date)))
                
                date = date - months*100
                days = np.asarray(list(map(lambda x: np.floor(x),date)))
                
                hours = np.asarray(list(map(lambda x: np.floor(x/100),time)))
                time = time - hours*100
                mins = np.asarray(list(map(lambda x: np.floor(x),time)))
                
                months_in_year = 12
                days_in_month = 30
                
                diffmins = np.diff(mins) % 60
                uebertragmins = -1*(np.diff(mins) < 0)
                diffhours = (np.diff(hours) + uebertragmins ) % 24
                uebertraghours = -1*(np.diff(hours) + uebertragmins  < 0)
                diffdays = (np.diff(days) + uebertraghours) % days_in_month
                uebertragdays = -1*(np.diff(days) + uebertraghours < 0)
                diffmonths = ( np.diff(months) + uebertragdays ) % months_in_year
                uebertragmonths = -1*(np.diff(months) + uebertragdays < 0)
                diffyears = np.diff(years) + uebertragmonths
                deltadays = diffmins/60/24 + diffhours/24 + diffdays + diffmonths * days_in_month + diffyears * months_in_year * days_in_month 
                return deltadays