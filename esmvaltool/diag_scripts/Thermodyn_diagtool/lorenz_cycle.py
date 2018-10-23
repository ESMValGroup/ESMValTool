#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 14:57:47 2018

@author: Valerio2
"""
import numpy as np
import sys
import math
import os
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import warnings
warnings.filterwarnings("ignore")
import imp
#sys.path.append('./diag_scripts/aux/Thermodynamics/')
import fluxogram as fluxogram
import srvfile_read as srv
from diagram_module import *
#import diag_scripts.aux.Thermodynamics.fluxogram as fluxogram
#import diag_scripts.aux.Thermodynamics.srvfile_read as srv
reload(srv)

g     = 9.81
R     = 287.00
cp    = 1003.5
aa    = 6.371E6
ps    = 101100.0
nwout = -999
nwoutf= -999
ntime = 1000000
nskip = 0
nw1   = 3
nw2   = 9
nw3   = 21

#log = open('lec_new.txt','w')    

class Lorenz_cycle():
    
    from lorenz_cycle import * 
        
    def lorenz(self, outpath, model, year, filesrv, filenc, plotfile, logfile):    
    #    PROGRAM MAIN
    #C
    #C***  *MAIN* MAIN PROGRAM
    #C
    #C     F.LUNKEIT         UNIHH        01.08.94
    #C
    #C     PURPOSE.
    #C     --------
    #C     READ NAMELIST AND CONTROL CALCULATION
    #C
    #C
    #C     INPUT.
    #C     ------
    #C     NAMELIST : STDIN
    #C     U,V,W,T  : FORT.10
    #C
    #C     DATA STRUCTURE: SERVIVE FORMAT (IH(8),FELD(NY,NRES))
    #C
    #C
    #C     OUTPUT.
    #C     -------
    #C     GLOBAL DIAGNOSTIC                              : STDOUT
    #C     TIME AVERAGED 3-D ENERGY TERMS (NY,NTP,NLEVEL) : FORT.20
    #C     TIME SERIES (G,N,S) MEANS (NTP,3)              : FORT.21  (SEE NWOUT)
    #C     TIME SERIES 3-D ENERGY TERMS (NY,NTP,NLEVEL)   : FORT.22  (SEE NWOUTF)
    #C
    #C     DATA STRUCTURE: SERVIVE FORMAT (IH(8),FELD(NTP,3) OR (NY,NTP))
    #C
    #C     CODES: 4001 KT       (TRANSIENT KIN. ENERGY)
    #C            4002 AT       (TRANSIENT AVAILABLE POT. ENERGY)
    #C            4003 A2KT     (A->K TRANSIENT)
    #C            4004 AE2AZT   (EDDY A -> ZONAL A TRANSIENT)
    #C            4005 KE2KZT   (EDDY K -> ZONAL K TRANSIENT)
    #C            4006 AT2AS    (TRANSIENT A -> STATIONARY A)
    #C            4007 KT2KS    (TRANSIENT K -> STATIONARY K)
    #C            4008 KS       (STATIONARY KIN ENERGY)
    #C            4009 AS       (STATIONARY AVAIL. POT. ENERGY)
    #C            4010 A2KS     (A -> K STATIONARY)
    #C            4011 AE2AZS   (EDDY A -> ZONAL A STATIONARY)
    #C            4012 KE2KZS   (EDDY K -> ZONAL K STATIONARY)
    #C
    #C     METHOD.
    #C     -------
    #C     READ NAMELIST FROM STDIN
    #C     GET DIMENSIONS BY *GETDIM*
    #C     CALL *ENERGIE*
    #C
    #C     EXTERNALS.
    #C     ----------
    #C     *GETDIM*  EXTRACT DIMENSIONS AND LEVELS FROM INPUT DATA
    #C     *ENERGIE* CALCULATE ENERGY TERMS
    #C
    #C     REFERENCES.
    #C     -----------
    #C     ULBRICH, U. AND P. SPETH (1991) METEOROL. ATMOS. PHYS. 45 125-138
    #C
    #C
    #C*    VARIABLE     TYPE     PURPOSE.
    #C     --------     ----     --------
    #C
    #C     GRAV         REAL     GRAVITATIONAL ACCELERATION  (9.81)
    #C     R            REAL     GAS CONSTANT (287.05)
    #C     CP           REAL     SPECIFIC HEAT OF DRY AIR (1005.46)
    #C     AA           REAL     PLANET RADIUS (6.371E6)
    #C     PS           REAL     MEAN SURFACE PRESSURE (1.E5)
    #C     NWOUT       INTEGET   INTERVAL FOR OUTPUT OF 1 D FIELDS  (NO OUTPUT)
    #C     NWOUTF      INTEGER   INTERVAL FOR OUTPUT OF FULL 3 D FIELDS (NO OUTPUT)
    #C     NTIME       INTEGER   NUMPER OF DATA TO BE USED (ALL; 100000)
    #C     NSKIP       INTEGER   NUMPER OF DATA TO BE SKIPED (NO)
    #C     NW1         INTEGER   LONG WAVES (1-NW1);FOR PRINT OUT (3)
    #C     NW2         INTEGER   SYNOPTIC WAVES (NW1+1-NW2); FOR PRINT OUT (10)
    #C     NW3         INTEGER   SHORT WAVES (NW2+1-NW3); FOR PRINT OUT (21)
    #C     TITLE       CHARACTER TO BE SET IN THE TITLE OF POSTSCRIPT OUTPUT
    #C
    
        lorenz = Lorenz_cycle()
        diagram = Diagram()
        
#        diagscript='diagram_module.py'
#        diagmod=imp.load_source('diagram_module',diagscript)
        
        log = open(logfile,'w')    
        
        log.write('########################################################\n')
        log.write('#                                                      #\n')
        log.write('#      LORENZ     ENERGY    CYCLE                      #\n')
        log.write('#                                                      #\n')
        log.write('########################################################\n')
     
        filename = filesrv
        filep    = filenc
        
    #    f=srv.srvfile(filename,hbyte="4",fbyte="4")
        f=srv.srvfile(filename,hbyte="8",fbyte="8")
    
        
        f.read_all()
    
    #    ta=f.field['130']
    #    ua=f.field['131']
    #    va=f.field['132']
    #    wap=f.field['135']
        ta=f.field['-1']
        ua=f.field['-2']
        va=f.field['-3']
        wap=f.field['-4']
    
        lev   = f.levels
        if(max(lev)<1000):
            lev=lev*100
            wap=wap*100
        #print(lev)
        nlev  = len(lev)
        nfc   = f.dim1
        nlat  = f.dim2
        lat = np.linspace(89.99,-89.99,num=nlat)
        dims  = np.shape(ta)
        ntime = dims[1]
        
        ta_r=ta[:,:,:,0::2]
        ta_i=ta[:,:,:,1::2]
        ua_r=ua[:,:,:,0::2]
        ua_i=ua[:,:,:,1::2]
        va_r=va[:,:,:,0::2]
        va_i=va[:,:,:,1::2]
        wap_r=wap[:,:,:,0::2]
        wap_i=wap[:,:,:,1::2]
        
        ta_c=ta_r + 1j * ta_i
        ua_c=ua_r + 1j * ua_i
        va_c=va_r + 1j * va_i
        wap_c=wap_r + 1j * wap_i
      
        ta_cinv=ta_c[::-1,:,:,:]
        ua_cinv=ua_c[::-1,:,:,:]
        va_cinv=va_c[::-1,:,:,:]
        wap_cinv=wap_c[::-1,:,:,:]
        ta_c=ta_cinv
        ua_c=ua_cinv
        va_c=va_cinv
        wap_c=wap_cinv
        
        log.write(' \n')
        log.write(' \n')
        log.write('INPUT DATA:\n')
        log.write('-----------\n')
        log.write(' \n')
        log.write('SPECTRAL RESOLUTION : {}\n'.format(nfc))
        log.write('NUMBER OF LATITUDES : {}\n'.format(nlat))
        log.write('NUMBER OF LEVEL : {}'.format(nlev))
        log.write('LEVEL  : {} Pa\n'.format(lev))
        log.write(' \n')
        
        ntp=nfc/2+1
        
        log.write('WAVES:\n')
        log.write(' \n')
        log.write('(1) : 1 - {}\n'.format(nw1))
        log.write('(2) : {} - {}\n'.format(nw1,nw2))
        log.write('(3) : {} - {}\n'.format(nw2,nw3))
    
    # Compute sigma level and dsigma  
        sig=np.zeros(len(lev))
        for jl in range(nlev):
            sig[jl]=lev[jl]/ps
        
        ds=np.zeros(len(lev))
        for jl in range(1,nlev-1,1):
            ds[jl]=0.5*abs(sig[jl+1]-sig[jl-1])
        ds[0]=sig[0]+0.5*abs(sig[1]-sig[0])
        ds[nlev-1]=1-sig[nlev-1]+0.5*abs(sig[nlev-1]-sig[nlev-2])
    #Compute Gaussian weights
        [pa,pw]=lorenz.gauaw(nlat)    
        gw=np.zeros(nlat)
        y=np.zeros(nlat)
        for jy in range(nlat): 
            gw[jy]=pw[jy]
            y[jy]=np.deg2rad(lat[jy])
        log.write(' \n')
        log.write('GLOBAL DIAGNOSTIC: \n')
        log.write('  \n')
        log.write('                            I GLOBAL I NORTH I SOUTH I\n')
        log.write('------------------------------------------------------\n')
        
    #Compute time mean
        ta_tmn  = np.nanmean(ta_c, axis=1)
        ua_tmn  = np.nanmean(ua_c, axis=1)
        va_tmn  = np.nanmean(va_c, axis=1)
        wap_tmn = np.nanmean(wap_c,axis=1)
    #Compute zonal mean of time means
        ta_zmn   = np.squeeze(np.real(ta_tmn[:,:,0])) 
        ua_zmn   = np.squeeze(np.real(ua_tmn[:,:,0])) 
        va_zmn   = np.squeeze(np.real(va_tmn[:,:,0]))
        wap_zmn  = np.squeeze(np.real(wap_tmn[:,:,0])) 
    #Compute global mean of time means
        ta_gmn=np.zeros(len(lev))
        ua_gmn=np.zeros(len(lev))
        va_gmn=np.zeros(len(lev))
        wap_gmn=np.zeros(len(lev))
        for j in range(nlev):
            #wap_gmn[j]  = np.nansum(wap_zmn[j,:]*gw)/np.nansum(gw)
            ta_gmn[j]   = np.nansum(ta_zmn[j,:]*np.cos(np.deg2rad(lat)))/np.nansum(np.cos(np.deg2rad(lat)))
            ua_gmn[j]   = np.nansum(ua_zmn[j,:]*np.cos(np.deg2rad(lat)))/np.nansum(np.cos(np.deg2rad(lat)))
            va_gmn[j]   = np.nansum(va_zmn[j,:]*np.cos(np.deg2rad(lat)))/np.nansum(np.cos(np.deg2rad(lat)))
            wap_gmn[j]  = np.nansum(wap_zmn[j,:]*np.cos(np.deg2rad(lat)))/np.nansum(np.cos(np.deg2rad(lat)))
    #Compute stability parameter    
        gam_tmn = lorenz.stabil(ta_gmn,lev,nlev)
        
        ek=np.zeros([nlev,ntime,nlat,ntp-1])
        ape=np.zeros([nlev,ntime,nlat,ntp-1])
        a2k=np.zeros([nlev,ntime,nlat,ntp-1])
        ae2az=np.zeros([nlev,ntime,nlat,ntp-1])
        ke2kz=np.zeros([nlev,ntime,nlat,ntp-1])
        at2as=np.zeros([nlev,ntime,nlat,ntp-1])
        kt2ks=np.zeros([nlev,ntime,nlat,ntp-1])
    
        for t in range(ntime):
            #print('Day {}'.format(t))
            ta_t     = ta_c[:,t,:,:]
            ua_t     = ua_c[:,t,:,:]
            va_t     = va_c[:,t,:,:]
            wap_t    = wap_c[:,t,:,:]
            ta_tan   = ta_c[:,t,:,:]-ta_tmn
            ua_tan   = ua_c[:,t,:,:]-ua_tmn
            va_tan   = va_c[:,t,:,:]-va_tmn
            wap_tan  = wap_c[:,t,:,:]-wap_tmn
        #Compute zonal means
            ta_tzmn   = np.squeeze(np.real(ta_t[:,:,0])) 
            ua_tzmn   = np.squeeze(np.real(ua_t[:,:,0])) 
            va_tzmn   = np.squeeze(np.real(va_t[:,:,0])) 
            wap_tzmn  = np.squeeze(np.real(wap_t[:,:,0])) 
            ta_tzan   = np.squeeze(np.real(ta_tan[:,:,0])) 
            ua_tzan   = np.squeeze(np.real(ua_tan[:,:,0])) 
            va_tzan   = np.squeeze(np.real(va_tan[:,:,0])) 
            wap_tzan  = np.squeeze(np.real(wap_tan[:,:,0])) 
        #Compute global means as a function of levels
            ta_tgmn=np.zeros(len(lev))
            ua_tgmn=np.zeros(len(lev))
            va_tgmn=np.zeros(len(lev))
            wap_tgmn=np.zeros(len(lev))
            ta_tgan=np.zeros(len(lev))
            ua_tgan=np.zeros(len(lev))
            va_tgan=np.zeros(len(lev))
            wap_tgan=np.zeros(len(lev))
            for j in range(nlev):
                ta_tgmn[j]   = np.nansum(ta_tzmn[j,:]*gw)/np.nansum(gw)
                ua_tgmn[j]   = np.nansum(ua_tzmn[j,:]*gw)/np.nansum(gw)
                va_tgmn[j]   = np.nansum(va_tzmn[j,:]*gw)/np.nansum(gw)
                wap_tgmn[j]  = np.nansum(wap_tzmn[j,:]*gw)/np.nansum(gw)
                ta_tgan[j]   = np.nansum(ta_tzan[j,:]*gw)/np.nansum(gw)
                ua_tgan[j]   = np.nansum(ua_tzan[j,:]*gw)/np.nansum(gw)
                va_tgan[j]   = np.nansum(va_tzan[j,:]*gw)/np.nansum(gw)
                wap_tgan[j]  = np.nansum(wap_tzan[j,:]*gw)/np.nansum(gw)
         #Compute kinetic energy   
            ek[:,t,:,:]=lorenz.makek(ua_tan,va_tan,nlat,ntp,nlev)
         #Compute available potential energy
            ape[:,t,:,:]=lorenz.makea(ta_tan,ta_tgan,gam_tmn,nlat,ntp,nlev)
         #Compute conversion between kin.en. and pot.en.
            a2k[:,t,:,:]=lorenz.mka2k(wap_tan,ta_tan,wap_tgan,ta_tgan,lev,nlat,ntp,nlev)
         #Compute conversion between zonal and eddy APE
            ae2az[:,t,:,:]=lorenz.mkaeaz(va_tan,wap_tan,ta_tan,ta_tmn,ta_gmn,lev,y,gam_tmn,nlat,ntp,nlev)
         #Compute conversion between zonal and eddy KE
            ke2kz[:,t,:,:]=lorenz.mkkekz(ua_tan,va_tan,wap_tan,ua_tmn,va_tmn,lev,y,nlat,ntp,nlev)
        #Compute conversion between stationary and transient eddy APE
            at2as[:,t,:,:]=lorenz.mkatas(ua_tan,va_tan,wap_tan,ta_tan,ta_tmn,gam_tmn,lev,y,nlat,ntp,nlev)   
        #Compute conversion between stationary and transient eddy KE
            kt2ks[:,t,:,:]=lorenz.mkktks(ua_tan,va_tan,wap_tan,ua_tmn,va_tmn,wap_tmn,lev,y,nlat,ntp,nlev)
        ##Average transient terms over time to obtain variables from 4001 to 4007 (transient eddies
        ## and zonal terms)
        #Variable 4001
        ek_tmn=np.nanmean(ek,axis=1)
        ek_tgmn=lorenz.globall_cg(ek_tmn,gw,ds,nlat,ntp,nlev)
        lorenz.table(ek_tgmn,ntp,'TOT. KIN. EN.    ',log)
        #Variable 4002
        ape_tmn=np.nanmean(ape,axis=1)
        ape_tgmn=lorenz.globall_cg(ape_tmn,gw,ds,nlat,ntp,nlev)
        lorenz.table(ape_tgmn,ntp,'TOT. POT. EN.   ',log)
        #Variable 4003
        a2k_tmn=np.nanmean(a2k,axis=1)
        a2k_tgmn=lorenz.globall_cg(a2k_tmn,gw,ds,nlat,ntp,nlev)
        lorenz.table_conv(a2k_tgmn,ntp,'KE -> APE (trans) ',log)
        #Variable 4004
        ae2az_tmn=np.nanmean(ae2az,axis=1)
        ae2az_tgmn=lorenz.globall_cg(ae2az_tmn,gw,ds,nlat,ntp,nlev)
        lorenz.table_conv(ae2az_tgmn,ntp,'AZ <-> AE (trans) ',log)
        #Variable 4005
        ke2kz_tmn=np.nanmean(ke2kz,axis=1)
        ke2kz_tgmn=lorenz.globall_cg(ke2kz_tmn,gw,ds,nlat,ntp,nlev)
        lorenz.table_conv(ke2kz_tgmn,ntp,'KZ <-> KE (trans) ',log)
        #Variable 4006
        at2as_tmn=np.nanmean(at2as,axis=1)
        at2as_tgmn=lorenz.globall_cg(at2as_tmn,gw,ds,nlat,ntp,nlev)
        lorenz.table_conv(at2as_tgmn,ntp,'ASE  <->  ATE   ',log)
        #Variable 4007
        kt2ks_tmn=np.nanmean(kt2ks,axis=1)
        kt2ks_tgmn=lorenz.globall_cg(kt2ks_tmn,gw,ds,nlat,ntp,nlev)
        lorenz.table_conv(kt2ks_tgmn,ntp,'KSE  <->  KTE   ',log)
     
        ##Use time averaged quantities to obtain variables from 4008 to 4012 (stationary terms)
        #Variable 4008
        ek_st=lorenz.makek(ua_tmn,va_tmn,nlat,ntp,nlev)
        ek_stgmn=lorenz.globall_cg(ek_st,gw,ds,nlat,ntp,nlev)
        lorenz.table(ek_stgmn,ntp,'STAT. KIN. EN.    ',log)
        #Variable 4009
        ape_st=lorenz.makea(ta_tmn,ta_gmn,gam_tmn,nlat,ntp,nlev)
        ape_stgmn=lorenz.globall_cg(ape_st,gw,ds,nlat,ntp,nlev)
        lorenz.table(ape_stgmn,ntp,'STAT. POT. EN.    ',log)
        #Variable 4010
        a2k_st=lorenz.mka2k(wap_tmn,ta_tmn,wap_gmn,ta_gmn,lev,nlat,ntp,nlev)
        a2k_stgmn=lorenz.globall_cg(a2k_st,gw,ds,nlat,ntp,nlev)
        lorenz.table_conv(a2k_stgmn,ntp,'KE -> APE (stat)',log)
        #Variable 4011
        ae2az_st=lorenz.mkaeaz(va_tmn,wap_tmn,ta_tmn,ta_tmn,ta_gmn,lev,y,gam_tmn,nlat,ntp,nlev)
        ae2az_stgmn=lorenz.globall_cg(ae2az_st,gw,ds,nlat,ntp,nlev)
        lorenz.table_conv(ae2az_stgmn,ntp,'AZ <-> AE (stat)',log)
        #Variable 4012
        ke2kz_st=lorenz.mkkekz(ua_tmn,va_tmn,wap_tmn,ua_tmn,va_tmn,lev,y,nlat,ntp,nlev)
        ke2kz_stgmn=lorenz.globall_cg(ke2kz_st,gw,ds,nlat,ntp,nlev)
        lorenz.table_conv(ke2kz_stgmn,ntp,'KZ <-> KE (stat)',log)
        
        #This part is for the flux diagram as in Ulbrich and Speth 1991
        apz    = '{:.2f}'.format(ape_tgmn[0,0]+ape_stgmn[0,0])
        az2kz  = '{:.2f}'.format(-1e5*(a2k_tgmn[0,0]+a2k_stgmn[0,0]))
        az2at  = '{:.2f}'.format(-1e5*np.nansum(ae2az_tgmn[0,1:ntp-1]))
        aps    = '{:.2f}'.format(np.nansum(ape_stgmn[0,1:ntp-1]))
        as2ks  = '{:.2f}'.format(1e5*np.nansum(a2k_stgmn[0,1:ntp-1]))
        apt    = '{:.2f}'.format(np.nansum(ape_tgmn[0,1:ntp-1]))
        at2kt  = '{:.2f}'.format(1e5*np.nansum(a2k_tgmn[0,1:ntp-1]))
        az2as  = '{:.2f}'.format(-1e5*np.nansum(ae2az_stgmn[0,1:ntp-1]))
        as2at  = '{:.2f}'.format(1e5*np.nansum(at2as_tgmn[0,1:ntp-1]))
        azin   = '{:.2f}'.format((float(az2at)+float(az2as)-float(az2kz)))
        asein  = '{:.2f}'.format((float(as2ks)+float(as2at)-float(az2as)))
        atein  = '{:.2f}'.format(float(at2kt)-float(az2at)-float(as2at))    
        kz     = '{:.2f}'.format(ek_tgmn[0,0]+ek_stgmn[0,0])
        kte    = '{:.2f}'.format(np.nansum(ek_tgmn[0,1:ntp-1]))
        kse    = '{:.2f}'.format(np.nansum(ek_stgmn[0,1:ntp-1]))
        kt2kz  = '{:.2f}'.format(1e5*np.nansum(ke2kz_tgmn[0,1:ntp-1]))
        kt2ks  = '{:.2f}'.format(-1e5*np.nansum(kt2ks_tgmn[0,1:ntp-1]))
        ks2kz  = '{:.2f}'.format(1e5*np.nansum(ke2kz_stgmn[0,1:ntp-1]))
        kteout = '{:.2f}'.format(float(at2kt)-float(kt2ks)-float(kt2kz))
        kseout = '{:.2f}'.format(float(kt2ks)+float(as2ks)-float(ks2kz))
        kzout  = '{:.2f}'.format(float(kt2kz)+float(ks2kz)-float(az2kz))    
        diagram.diagram(plotfile,azin,apz,asein,aps,atein,apt,as2ks,at2kt,kteout,
                               kte,kseout,kse,kzout,kz,az2kz,az2at,az2as,
                               as2at,kt2kz,kt2ks,ks2kz)
        lec_strength=float(kteout)+float(kseout)+float(kzout)
        
        #Print out to NetCDF files
        ek_aux=np.zeros([nlev,nlat,ntp-1])
        ape_aux=np.zeros([nlev,nlat,ntp-1])
        a2k_aux=np.zeros([nlev,nlat,ntp-1])
        ae2az_aux=np.zeros([nlev,nlat,ntp-1])
        ke2kz_aux=np.zeros([nlev,nlat,ntp-1])
        for l in range(nlev):
            ek_aux[l,:,:]    = ek_tmn[l,:,:]*ds[l]
            ape_aux[l,:,:]   = ape_tmn[l,:,:]*ds[l]
            a2k_aux[l,:,:]   = a2k_tmn[l,:,:]*ds[l]
            ae2az_aux[l,:,:] = ae2az_tmn[l,:,:]*ds[l]
            ke2kz_aux[l,:,:] = ke2kz_tmn[l,:,:]*ds[l]
        ek_vmn=np.nansum(ek_aux,axis=0)/np.nansum(ds)
        ape_vmn=np.nansum(ape_aux,axis=0)/np.nansum(ds)
        a2k_vmn=np.nansum(a2k_aux,axis=0)/np.nansum(ds)
        ae2az_vmn=np.nansum(ae2az_aux,axis=0)/np.nansum(ds)
        ke2kz_vmn=np.nansum(ke2kz_aux,axis=0)/np.nansum(ds)
    
        nc_f=outpath+'/ek_tmap_{}_{}.nc'.format(model,year)
        lorenz.removeif(nc_f)
        lorenz.pr_output(ek_vmn, 'ek', filep, nc_f, 1, verb=True)
        nc_f=outpath+'/ape_tmap_{}_{}.nc'.format(model,year)
        lorenz.removeif(nc_f)
        lorenz.pr_output(ape_vmn, 'ape', filep, nc_f, 1, verb=True)
        nc_f=outpath+'/a2k_tmap_{}_{}.nc'.format(model,year)
        lorenz.removeif(nc_f)
        lorenz.pr_output(a2k_vmn, 'a2k', filep, nc_f, 1, verb=True)
        nc_f=outpath+'/ae2az_tmap_{}_{}.nc'.format(model,year)
        lorenz.removeif(nc_f)
        lorenz.pr_output(ae2az_vmn, 'ae2az', filep, nc_f, 1, verb=True)
        nc_f=outpath+'/ke2kz_tmap_{}_{}.nc'.format(model,year)
        lorenz.removeif(nc_f)
        lorenz.pr_output(ke2kz_vmn, 'ke2kz', filep, nc_f, 1, verb=True)
        
        log.close() 
        
        return lec_strength
    
    def bsslzr(self, kdim):
            
        NDIM=50
        
        PI = math.pi
      
        zbes = [ 2.4048255577, 5.5200781103, 8.6537279129,  11.7915344391,  
                14.9309177086,  18.0710639679, 21.2116366299,  24.3524715308,  
                27.4934791320,  30.6346064684, 33.7758202136,  36.9170983537,  
                40.0584257646,  43.1997917132, 46.3411883717,  49.4826098974,  
                52.6240518411,  55.7655107550, 58.9069839261,  62.0484691902,  
                65.1899648002,  68.3314693299, 71.4729816036,  74.6145006437,  
                77.7560256304,  80.8975558711, 84.0390907769,  87.1806298436,  
                90.3221726372,  93.4637187819, 96.6052679510,  99.7468198587, 
                102.8883742542, 106.0299309165,109.1714896498, 112.3130502805, 
                115.4546126537, 118.5961766309, 121.7377420880, 124.8793089132, 
                128.0208770059, 131.1624462752, 134.3040166383, 137.4455880203, 
                140.5871603528, 143.7287335737, 146.8703076258, 150.0118824570, 
                153.1534580192, 156.2950342685 ]
    
        pbes=np.zeros(kdim)
        idim = min([kdim,NDIM])
        pbes[0:idim] = zbes[0:idim]
        for j in range(idim,kdim-1,1):
            pbes[j] = pbes[j-1] + PI
          
        return(pbes)
    
    
                
    def gauaw(self, ny):
        
        lorenz = Lorenz_cycle()
        
        c = (1-(2/math.pi)**2)/4
        eps = 0.00000000000001
      
        kk = ny/2
        #print(kk)
        pa=np.zeros(ny)
        #print(len(bsslzr(kk)))
        pa[0:kk]=lorenz.bsslzr(kk)
        #print(len(pa[0:kk]))
        pw=np.zeros(ny)
        for i in range(kk):
            xz = np.cos(pa[i]/math.sqrt((ny+0.5)**2+c))
            iterr = 0.
            zsp  = 1.0
            while (abs(zsp) > eps and iterr <= 10):
                #print(zsp)
                pkm1 = xz
                pkm2 = 1.0
                for n in range(2,ny,1):
                    pk   = ((n*2-1.0)*xz*pkm1 - (n-1.0) * pkm2) / n
                    pkm2 = pkm1
                    pkm1 = pk
                pkm1  = pkm2
                pkmrk = (ny * (pkm1 - xz * pk)) / (1.0 - xz**2)
                zsp   = pk / pkmrk
                xz    = xz - zsp
                iterr  = iterr + 1
            if iterr > 15:
                sys.exit("*** no convergence in gauaw ***")
            pa[i] = xz
            pw[i] = (2.0 * (1.0 - xz**2))/((ny**2)*(pkm1**2))
            pa[ny-1-i] = -pa[i]
            pw[ny-1-i] =  pw[i]
            
        psi = pa
        pgw = pw
      
        return psi,pgw
    
    
    
    def globall_cg(self, d3v,gw,ds,nlat,ntp,nlev):
        #C***  *GLOBAL* CALCULATE GLOBAL AND HEMISPHERIC MEANS
        #C
        #C     F.LUNKEIT         UNIHH        01.08.94
        #C
        #C     PURPOSE.
        #C     --------
        #C     CALCULATE GLOBAL AND HEMISPHERIC MEANS
        #C
        #C
    
        gmn=np.zeros([3,ntp-1])
        aux1=np.zeros([nlev,nlat/2,ntp-1])
        aux2=np.zeros([nlev,nlat/2,ntp-1])
        aux1v=np.zeros([nlev,ntp-1])
        aux2v=np.zeros([nlev,ntp-1])
        
        nhem=nlat/2
        fac=1/g*ps/1e5
        for l in range(nlev):
            for i in range(nhem):
                aux1[l,i,:]=fac*np.real(d3v[l,i,:])*gw[i]
                aux2[l,i,:]=fac*np.real(d3v[l,i+nhem-1,:])*gw[i+nhem-1]
            aux1v[l,:]=np.nansum(aux1[l,:,:],axis=0)/np.nansum(gw[0:nhem])*ds[l]
            aux2v[l,:]=np.nansum(aux2[l,:,:],axis=0)/np.nansum(gw[0:nhem])*ds[l]
        gmn[1,:]=(np.nansum(aux1v,axis=0)/np.nansum(ds))
        gmn[2,:]=(np.nansum(aux2v,axis=0)/np.nansum(ds))
        
        gmn[0,:]=0.5*(gmn[1,:]+gmn[2,:])
    
        return(gmn)
    
    
    def makek(self, u,v,nlat,ntp,nlev):
    #C
    #C
    #C
    #C***  *MAKEK* CALCULATE KINETIC ENERGY
    #C
    #C     F.LUNKEIT         UNIHH        01.08.94
    #C
    #C     PURPOSE.
    #C     --------
    #C     CALCULATE SPECTRAL COMPONENTS OF KINETIC ENERGY FRON U AND V
    #C
    #C
    
        ek=np.zeros([nlev,nlat,ntp-1])
        ck1=u*np.conj(u)
        ck2=v*np.conj(v)
        ek[:,:,0]=0.5*np.real(u[:,:,0]*u[:,:,0]+v[:,:,0]*v[:,:,0])
        ek=np.real(ck1+ck2)
        ek[:,:,0]=0.5*np.real(u[:,:,0]*u[:,:,0]+v[:,:,0]*v[:,:,0])
        
        return(ek)
    
        
    def makea(self, t,tg,gam,nlat,ntp,nlev):
    #C***  *MAKEA* CALCULATE AVAIL. POT. ENERGY
    #C
    #C     F.LUNKEIT         UNIHH        01.08.94
    #C
    #C     PURPOSE.
    #C     --------
    #C     CALCULATE SPECTRAL COMPONENTS OF AVAIL. POT. ENERGY FROM T
    #C
    #C
          
        a=gam[:,np.newaxis,np.newaxis]*np.real(t*np.conj(t))
        a[:,:,0]=gam[:,np.newaxis]*0.5*np.real((t[:,:,0]-tg[:,np.newaxis])*(t[:,:,0]-tg[:,np.newaxis]))
        
        return(a)
    
    def mka2k(self, wap,t,wg,tg,p,nlat,ntp,nlev):
    
    #C***  *MKA2K* CALCULATE CONVERSION A->K
    #C
    #C     F.LUNKEIT         UNIHH        01.08.94
    #C
    #C     PURPOSE.
    #C     --------
    #C     CALCULATE SPECTRAL COMPONENTS OF CONVERSION A->K
    #C
    #C
        
        a2k=-R/p[:,np.newaxis,np.newaxis]*np.real(t*np.conj(wap)+np.conj(t)*wap)
        a2k[:,:,0]=-R/p[:,np.newaxis]*np.real((t[:,:,0]-tg[:,np.newaxis])*(wap[:,:,0]-wg[:,np.newaxis]))
        
        return(a2k)
        
    
    def mkaeaz(self, v,wap,t,tt,ttg,p,lat,gam,nlat,ntp,nlev):
    #C
    #C
    #C***  *MKAEAZ* CALCULATE CONVERSION EDDY A->ZONAL A
    #C
    #C     F.LUNKEIT         UNIHH        01.08.94
    #C
    #C     PURPOSE.
    #C     --------
    #C     CALCULATE SPECTRAL COMPONENTS OF CONVERSION EDDY A->ZONAL A
    #C
    #C
    
        ae2az = np.zeros([nlev,nlat,ntp-1])
        dtdp  = np.zeros([nlev,nlat])
        dtdy  = np.zeros([nlev,nlat])
        for l in np.arange(nlev):
            if l == 0:
                t1   = np.real(tt[l,:,0])-ttg[l]
                t2   = np.real(tt[l+1,:,0])-ttg[l+1]
                dtdp[l,:] = (t2-t1)/(p[l+1]-p[l])
            elif l == nlev-1:
                t1   = np.real(tt[l-1,:,0])-ttg[l-1]
                t2   = np.real(tt[l,:,0])-ttg[l]
                dtdp[l,:] = (t2-t1)/(p[l]-p[l-1])
            else:
                t1   = np.real(tt[l,:,0])-ttg[l]
                t2   = np.real(tt[l+1,:,0])-ttg[l+1]
                dtdp1= (t2-t1)/(p[l+1]-p[l])
                t2   = t1
                t1   = np.real(tt[l-1,:,0])-ttg[l-1]
                dtdp2= (t2-t1)/(p[l]-p[l-1])
                dtdp[l,:] = (dtdp1*(p[l]-p[l-1])+dtdp2*(p[l+1]-p[l]))/(p[l+1]-p[l-1])
            dtdp[l,:] = dtdp[l,:]-R/(cp*p[l])*(tt[l,:,0]-ttg[l])
        
        for i in np.arange(nlat):        
            if i == 0:
                t1        = np.real(tt[:,i,0])
                t2        = np.real(tt[:,i+1,0])
                dtdy[:,i] = (t2-t1)/(lat[i+1]-lat[i])
            elif i == nlat-1:
                t1        = np.real(tt[:,i-1,0])
                t2        = np.real(tt[:,i,0])
                dtdy[:,i] = (t2-t1)/(lat[i]-lat[i-1])
            else:
                t1        = np.real(tt[:,i-1,0])
                t2        = np.real(tt[:,i+1,0])
                dtdy[:,i] = (t2-t1)/(lat[i+1]-lat[i-1])
        dtdy  = dtdy/aa
        c1 =  np.real(v*np.conj(t)+t*np.conj(v))
        c2 =  np.real(wap*np.conj(t)+t*np.conj(wap))
        ae2az = gam[:,np.newaxis,np.newaxis]*(dtdy[:,:,np.newaxis]*c1+dtdp[:,:,np.newaxis]*c2)      
        ae2az[:,:,0]=0.
    
        return(ae2az)
    
    
    
    def mkkekz(self, u,v,wap,ut,vt,p,lat,nlat,ntp,nlev):
    #C
    #C***  *MKKEKZ* CALCULATE CONVERSION EDDY K->ZONAL K
    #C
    #C     F.LUNKEIT         UNIHH        01.08.94
    #C
    #C     PURPOSE.
    #C     --------
    #C     CALCULATE SPECTRAL COMPONENTS OF CONVERSION EDDY K->ZONAL K
    #C
    #C
        dudp  = np.zeros([nlev,nlat])
        dvdp  = np.zeros([nlev,nlat])
        dudy  = np.zeros([nlev,nlat])
        dvdy  = np.zeros([nlev,nlat])
        
        for l in np.arange(nlev):
            if l == 0:
                dudp[l,:]  = (np.real(ut[l+1,:,0]-ut[l,:,0]))/(p[l+1]-p[l])
                dvdp[l,:]  = (np.real(vt[l+1,:,0]-vt[l,:,0]))/(p[l+1]-p[l])
            elif l == nlev-1:
                dudp[l,:]  = (np.real(ut[l,:,0]-ut[l-1,:,0]))/(p[l]-p[l-1])
                dvdp[l,:]  = (np.real(vt[l,:,0]-vt[l-1,:,0]))/(p[l]-p[l-1])
            else:
                dudp1  = (np.real(ut[l+1,:,0]-ut[l,:,0]))/(p[l+1]-p[l])
                dvdp1  = (np.real(vt[l+1,:,0]-vt[l,:,0]))/(p[l+1]-p[l])
                dudp2  = (np.real(ut[l,:,0]-ut[l-1,:,0]))/(p[l]-p[l-1])
                dvdp2  = (np.real(vt[l,:,0]-vt[l-1,:,0]))/(p[l]-p[l-1])
                dudp[l,:] = (dudp1*(p[l]-p[l-1])+dudp2*(p[l+1]-p[l]))/(p[l+1]-p[l-1])
                dvdp[l,:] = (dvdp1*(p[l]-p[l-1])+dvdp2*(p[l+1]-p[l]))/(p[l+1]-p[l-1])
           
        for i in np.arange(nlat):    
            if i == 0:
                dudy[:,i]  = (np.real(ut[:,i+1,0]-ut[:,i,0]))/(lat[i+1]-lat[i])
                dvdy[:,i]  = (np.real(vt[:,i+1,0]-vt[:,i,0]))/(lat[i+1]-lat[i])
            elif i == nlat-1:
                dudy[:,i]  = (np.real(ut[:,i,0]-ut[:,i-1,0]))/(lat[i]-lat[i-1])
                dvdy[:,i]  = (np.real(vt[:,i,0]-vt[:,i-1,0]))/(lat[i]-lat[i-1])
            else:
                dudy[:,i]  = (np.real(ut[:,i+1,0]-ut[:,i-1,0]))/(lat[i+1]-lat[i-1])
                dvdy[:,i]  = (np.real(vt[:,i+1,0]-vt[:,i-1,0]))/(lat[i+1]-lat[i-1])
        dudy  = dudy/aa
        dvdy  = dvdy/aa
        
        c1   = np.zeros([nlev,nlat,ntp-1])
        c2   = np.zeros([nlev,nlat,ntp-1])
        c3   = np.zeros([nlev,nlat,ntp-1])
        c4   = np.zeros([nlev,nlat,ntp-1])
        c5   = np.zeros([nlev,nlat,ntp-1])
        c6   = np.zeros([nlev,nlat,ntp-1])
        ke2kz= np.zeros([nlev,nlat,ntp-1])
        
        uu = u*np.conj(u)+u*np.conj(u)
        uv = u*np.conj(v)+v*np.conj(u)
        vv = v*np.conj(v)+v*np.conj(v)
        uw = u*np.conj(wap)+wap*np.conj(u)
        vw = v*np.conj(wap)+wap*np.conj(v)
        
        for i in np.arange(nlat):
            c1[:,i,:]=dudy[:,i][:,np.newaxis]*uv[:,i,:]
            c2[:,i,:]=dvdy[:,i][:,np.newaxis]*vv[:,i,:]
    #        c5[:,i,:]= np.tan(np.deg2rad(lat[i]))/aa*np.real(ut[:,i,0])[:,np.newaxis]*(uv[:,i,:])
    #        c6[:,i,:]=-np.tan(np.deg2rad(lat[i]))/aa*np.real(vt[:,i,0])[:,np.newaxis]*(uu[:,i,:])
            c5[:,i,:]= np.tan(lat[i])/aa*np.real(ut[:,i,0])[:,np.newaxis]*(uv[:,i,:])
            c6[:,i,:]=-np.tan(lat[i])/aa*np.real(vt[:,i,0])[:,np.newaxis]*(uu[:,i,:])    
        for l in np.arange(nlev):
            c3[l,:,:]=dudp[l,:][:,np.newaxis]*uw[l,:,:]
            c4[l,:,:]=dvdp[l,:][:,np.newaxis]*vw[l,:,:]
        ke2kz=(c1+c2+c3+c4+c5+c6)
        ke2kz[:,:,0]=0.
    
        return(ke2kz)
    
    
    def mkatas(self, u,v,wap,t,tt,gw,p,lat,nlat,ntp,nlev):
    #C***  *MKATAS* CALCULATE CONVERSION TRANSIENT A -> STATIONARY A
    #C
    #C     F.LUNKEIT         UNIHH        01.08.94
    #C
    #C     PURPOSE.
    #C     --------
    #C     CALCULATE SPECTRAL COMPONENTS OF CONVERSION TRANSIENT A -> STATIONARY A
    #C
    #C
    #C     REFERENCES.
    #C     -----------
    #C     ULBRICH, U. AND P. SPETH (1991) METEOROL. ATMOS. PHYS. 45 125-138
    
        tr=np.fft.irfft(t,axis=1)
        ur=np.fft.irfft(u,axis=1)
        vr=np.fft.irfft(v,axis=1)
        wr=np.fft.irfft(wap,axis=1)
        tur=tr*ur
        tvr=tr*vr
        twr=tr*wr
        tu=np.fft.rfft(tur,axis=1)
        tv=np.fft.rfft(tvr,axis=1)
        tw=np.fft.rfft(twr,axis=1)
        
        c1=tu*np.conj(tt)-tt*np.conj(tu)
        c6=tw*np.conj(tt)-tt*np.conj(tw)
        c2=np.zeros([nlev,nlat,ntp-1])
        c3=np.zeros([nlev,nlat,ntp-1])
        c5=np.zeros([nlev,nlat,ntp-1])
        for i in range(nlat):
            if i == 0:
                c2[:,i,:]  = tv[:,i,:]*np.conj(tt[:,i+1,:]-tt[:,i,:])/(aa*(lat[i+1]-lat[i]))
                c3[:,i,:]  = np.conj(tv[:,i,:])*(tt[:,i+1,:]-tt[:,i,:])/(aa*(lat[i+1]-lat[i]))
            elif i == nlat-1:
                c2[:,i,:]  = tv[:,i,:]*np.conj(tt[:,i,:]-tt[:,i-1,:])/(aa*(lat[i]-lat[i-1]))
                c3[:,i,:]  = np.conj(tv[:,i,:])*(tt[:,i,:]-tt[:,i-1,:])/(aa*(lat[i]-lat[i-1]))
            else:
                c2[:,i,:]  = tv[:,i,:]*np.conj(tt[:,i+1,:]-tt[:,i-1,:])/(aa*(lat[i+1]-lat[i-1]))
                c3[:,i,:]  = np.conj(tv[:,i,:])*(tt[:,i+1,:]-tt[:,i-1,:])/(aa*(lat[i+1]-lat[i-1]))
        
        for l in range(nlev):
            if l == 0:
                c5[l,:,:]  = (tt[l+1,:,:]-t[l,:,:])/(p[l+1]-p[l])
            elif l == nlev-1:
                c5[l,:,:]  = (tt[l,:,:]-t[l-1,:,:])/(p[l]-p[l-1])
            else:
                c51  = (tt[l+1,:,:]-t[l,:,:])/(p[l+1]-p[l])
                c52  = (tt[l,:,:]-t[l-1,:,:])/(p[l]-p[l-1])
                c5[l,:,:]  = (c51*(p[l]-p[l-1])+c52*(p[l+1]-p[l]))/(p[l+1]-p[l-1])
                
        K = np.arange(0,ntp-1)            
        at2as=gw[:,np.newaxis,np.newaxis]*((K-1)[np.newaxis,np.newaxis,:]*
             np.imag(c1)/(aa*np.cos(lat[np.newaxis,:,np.newaxis])) \
            + np.real(c2+c3)+np.real(tw*np.conj(c5)+\
            np.conj(tw)*c5)+R/(cp*p[:,np.newaxis,np.newaxis])*np.real(c6))        
        at2as[:,:,0]=0.
    
        return(at2as)    
       
    
    
    def mkktks(self, u,v,wap,ut,vt,wt,p,lat,nlat,ntp,nlev):    
    #C***  *MKKTKS* CALCULATE CONVERSION TRANSIENT K -> STATIONARY K
    #C
    #C     F.LUNKEIT         UNIHH        01.08.94
    #C
    #C     PURPOSE.
    #C     --------
    #C     CALCULATE SPECTRAL COMPONENTS OF CONVERSION TRANSIENT K -> STATIONARY K
    #C
    
        kt2ks = np.zeros([nlev,nlat,ntp-1])
        c1  = np.zeros([nlev,nlat,ntp-1])
        c21 = np.zeros([nlev,nlat,ntp-1])
        c22 = np.zeros([nlev,nlat,ntp-1])
        c3  = np.zeros([nlev,nlat,ntp-1])
        c41 = np.zeros([nlev,nlat,ntp-1])
        c42 = np.zeros([nlev,nlat,ntp-1])
        c5  = np.zeros([nlev,nlat,ntp-1])
        c6  = np.zeros([nlev,nlat,ntp-1])
        dut = np.zeros([nlev,nlat,ntp-1])
        dvt = np.zeros([nlev,nlat,ntp-1])
        dlat=np.zeros([nlat])
        
        ur=np.fft.irfft(u,axis=1)
        vr=np.fft.irfft(v,axis=1)
        uur=ur*ur
        uvr=ur*vr
        vvr=vr*vr
        uu=np.fft.rfft(uur,axis=1)
        vv=np.fft.rfft(vvr,axis=1)
        uv=np.fft.rfft(uvr,axis=1)
        
        c1=uu*np.conj(ut)-ut*np.conj(uu)
        c3=uv*np.conj(ut)+ut*np.conj(uv)
        c5=uu*np.conj(vt)+vt*np.conj(uu)
        c6=uv*np.conj(vt)-vt*np.conj(uv)
        for i in range(nlat):
            if i == 0:
                dut[:,i,:]=(ut[:,i+1,:]-ut[:,i,:])
                dvt[:,i,:]=(vt[:,i+1,:]-vt[:,i,:])
                dlat[i]=(lat[i+1]-lat[i])
            elif i == nlat-1:
                dut[:,i,:]=(ut[:,i,:]-ut[:,i-1,:])
                dvt[:,i,:]=(vt[:,i,:]-vt[:,i-1,:])
                dlat[i]=(lat[i]-lat[i-1])
            else:
                dut[:,i,:]=(ut[:,i+1,:]-ut[:,i-1,:])
                dvt[:,i,:]=(vt[:,i+1,:]-vt[:,i-1,:])
                dlat[i]=(lat[i+1]-lat[i-1])           
        c21 = np.conj(uu)*dut/dlat[np.newaxis,:,np.newaxis]
        c22 = uu*np.conj(dut)/dlat[np.newaxis,:,np.newaxis]
        c41 = np.conj(vv)*dvt/dlat[np.newaxis,:,np.newaxis]
        c42 = vv*np.conj(dvt)/dlat[np.newaxis,:,np.newaxis]
        
        K=np.arange(0,ntp-1)      
        kt2ks=(K-1)[np.newaxis,np.newaxis,:]/(aa*np.cos(lat)[np.newaxis,:,np.newaxis]) \
            *np.imag(c1+c6)+np.real(c21+c22+c41+c42)/aa+ \
            np.tan(lat)[np.newaxis,:,np.newaxis]*np.real(c1-c5)/aa
        kt2ks[:,:,0]=0.
        
        return(kt2ks)   
    
    
    def ncdump(self, nc_fid,key,verb):
        """
        Prints the NetCDF file attributes for a given key
    
        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
    
        # NetCDF global attributes
        nc_attrs = nc_fid.ncattrs()
        nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
        nc_vars = [var for var in nc_fid.variables]  # list of nc variables
        
        return nc_attrs, nc_dims, nc_vars
    
    
    
    def pr_output(self, varo, varname, filep, nc_f, opt, verb=True):
        '''
        NAME
            NetCDF with Python
        PROGRAMMER(S)
            Chris Slocum
        REVISION HISTORY
            20140320 -- Initial version created and posted online
            20140722 -- Added basic error handling to ncdump
                        Thanks to K.-Michael Aye for highlighting the issue
        '''
    
        lorenz = Lorenz_cycle()
        
        #I need a proxy dataset to extract the information on coordinates and global attributes
        nc_fid = Dataset(filep, 'r')  # Dataset is the class behavior to open the file
                                     # and create an instance of the ncCDF4 class
        nc_attrs, nc_dims, nc_vars = lorenz.ncdump(nc_fid,'var130',verb)
        # Extract data from NetCDF file
        lats = nc_fid.variables['lat'][:]  # extract/copy the data
        time = nc_fid.variables['time'][:]
        wave = nc_fid.variables['wave'][:]
        ntp=len(wave)/2
        
        # Writing NetCDF files
        w_nc_fid = Dataset(nc_f, 'w', format='NETCDF4')
        w_nc_fid.description = "NCEP/NCAR Reanalysis %s from its value at %s. %s"
        
        # Using our previous dimension info, we can create the new time dimension
        # Even though we know the size, we are going to set the size to unknown
        
        if opt ==0:
            w_nc_fid.createDimension('time', None)
            w_nc_dim = w_nc_fid.createVariable('time', nc_fid.variables['time'].dtype,\
                                               ('time',))
            for ncattr in nc_fid.variables['time'].ncattrs():
                w_nc_dim.setncattr(ncattr, nc_fid.variables['time'].getncattr(ncattr))
            # Assign the dimension data to the new NetCDF file.
            w_nc_fid.variables['time'][:] = time
            w_nc_fid.createDimension('lat', len(lats))
            w_nc_dim = w_nc_fid.createVariable('lat', nc_fid.variables['lat'].dtype,\
                                               ('lat',))
            for ncattr in nc_fid.variables['lat'].ncattrs():
                w_nc_dim.setncattr(ncattr, nc_fid.variables['lat'].getncattr(ncattr))
            w_nc_fid.variables['lat'][:] = lats
            w_nc_fid.createDimension('wave', ntp)
            w_nc_dim = w_nc_fid.createVariable('wave', nc_fid.variables['wave'].dtype,\
                                               ('wave',))
            w_nc_fid.variables['wave'][:] = wave[0:ntp]
            w_nc_var = w_nc_fid.createVariable(varname, 'f8', ('time','lat','wave'))
            lorenz.varatts(w_nc_var,varname,0,0)
            w_nc_fid.variables[varname][:] = varo
            w_nc_fid.close()  # close the new file
        elif opt ==1:
            w_nc_fid.createDimension('lat', len(lats))
            w_nc_dim = w_nc_fid.createVariable('lat', nc_fid.variables['lat'].dtype,\
                                               ('lat',))
            for ncattr in nc_fid.variables['lat'].ncattrs():
                w_nc_dim.setncattr(ncattr, nc_fid.variables['lat'].getncattr(ncattr))
            w_nc_fid.variables['lat'][:] = lats
            w_nc_fid.createDimension('wave', ntp)
            w_nc_dim = w_nc_fid.createVariable('wave', nc_fid.variables['wave'].dtype,\
                                               ('wave',))
            for ncattr in nc_fid.variables['wave'].ncattrs():
                w_nc_dim.setncattr(ncattr, nc_fid.variables['wave'].getncattr(ncattr))
            w_nc_fid.variables['wave'][:] = wave[0:ntp]
            w_nc_var = w_nc_fid.createVariable(varname, 'f8', ('lat','wave'))
            lorenz.varatts(w_nc_var,varname,1,0)
            w_nc_fid.variables[varname][:] = varo
            w_nc_fid.close()  # close the new file
        elif opt == 2:
            w_nc_fid.createDimension('lat', len(lats))
            w_nc_dim = w_nc_fid.createVariable('lat', nc_fid.variables['lat'].dtype,\
                                               ('lat',))
            for ncattr in nc_fid.variables['lat'].ncattrs():
                w_nc_dim.setncattr(ncattr, nc_fid.variables['lat'].getncattr(ncattr))
            w_nc_fid.variables['lat'][:] = lats
            w_nc_var = w_nc_fid.createVariable(varname, 'f8', ('lat'))
            lorenz.varatts(w_nc_var,varname,1,0)
            w_nc_fid.variables[varname][:] = varo
            w_nc_fid.close()  # close the new file
        elif opt == 3:
            w_nc_fid.createDimension('hem', 3)
            w_nc_dim = w_nc_fid.createVariable('hem', nc_fid.variables['hem'].dtype,\
                                               ('hem',))
            w_nc_fid.variables['hem'][:] = [0,1,2]
            w_nc_fid.createDimension('time', None)
            w_nc_dim = w_nc_fid.createVariable('time', nc_fid.variables['time'].dtype,\
                                               ('time',))
            for ncattr in nc_fid.variables['time'].ncattrs():
                w_nc_dim.setncattr(ncattr, nc_fid.variables['time'].getncattr(ncattr))
            w_nc_fid.variables['time'][:] = time
            w_nc_var = w_nc_fid.createVariable(varname, 'f8', ('hem','time'))
            lorenz.varatts(w_nc_var,varname,0,0)
            w_nc_fid.variables[varname][:] = varo
            w_nc_fid.close()  # close the new file
        elif opt == 4:
            w_nc_fid.createDimension('hem', 3)
            w_nc_dim = w_nc_fid.createVariable('hem', nc_fid.variables['hem'].dtype,\
                                               ('hem',))
            w_nc_fid.variables['hem'][:] = [0,1,2]
            w_nc_fid.createDimension('time', None)
            w_nc_dim = w_nc_fid.createVariable('time', nc_fid.variables['time'].dtype,\
                                               ('time',))
            for ncattr in nc_fid.variables['time'].ncattrs():
                w_nc_dim.setncattr(ncattr, nc_fid.variables['time'].getncattr(ncattr))
            w_nc_fid.variables['time'][:] = time
            w_nc_fid.createDimension('wave', ntp)
            w_nc_dim = w_nc_fid.createVariable('wave', nc_fid.variables['wave'].dtype,\
                                               ('wave',))
            for ncattr in nc_fid.variables['wave'].ncattrs():
                w_nc_dim.setncattr(ncattr, nc_fid.variables['wave'].getncattr(ncattr))
            w_nc_fid.variables['wave'][:] = wave[0:ntp]
            w_nc_var = w_nc_fid.createVariable(varname, 'f8', ('hem','time','wave'))
            lorenz.varatts(w_nc_var,varname,0,0)
            w_nc_fid.variables[varname][:] = varo
            w_nc_fid.close()  # close the new file
                
        nc_fid.close()
    
    
    
        
    def removeif(self, filename):
        try:
            os.remove(filename)
        except OSError:
            pass
    
    
    
    def stabil(self, ta_gmn,p,nlev):
    #C
    #C
    #C***  *STABIL* CALCULATE STABILITY PARAMETER
    #C
    #C     F.LUNKEIT         UNIHH        01.08.94
    #C
    #C     PURPOSE.
    #C     --------
    #C     CALCULATE STABILITY PARAMETER NEEDED FOR POT ENERGY
    #C
    #C     METHOD.
    #C     -------
    #C     GAMMA=CP/(T-P*DTDP*CP/R)
    #C
    #C     REFERENCES.
    #C     -----------
    #C     ULBRICH, U. AND P. SPETH (1991) METEOROL. ATMOS. PHYS. 45 125-138
    #C
    #C
    
        cpdr = cp/R
        
        t=ta_gmn
        
        
        gs=np.zeros(nlev)
        for i in range(nlev):
            if i == 0:
                dtdp=(t[i+1]-t[i])/(p[i+1]-p[i])
            elif i == nlev-1:
                dtdp =(t[i]-t[i-1])/(p[i]-p[i-1])
            else:
                dtdp1=(t[i+1]-t[i])/(p[i+1]-p[i])
                dtdp2=(t[i]-t[i-1])/(p[i]-p[i-1])
                dtdp = (dtdp1*(p[i]-p[i-1])+dtdp2*(p[i+1]-p[i]))/(p[i+1]-p[i-1])
            gs[i]=cp/(t[i]-p[i]*dtdp*cpdr)
    
        return gs
    
    
    def table(self, varin,ntp,name,log):
        
        varzon = varin[:,0]
        vared  = np.nansum(varin[:,1:ntp-1],axis=1)
        vared1 = np.nansum(varin[:,1:nw1-1],axis=1)
        vared2 = np.nansum(varin[:,nw1:nw2-1],axis=1)
        vared3 = np.nansum(varin[:,nw2:nw3-1],axis=1)
        vartot = varzon+vared
        
        log.write(' {} TOTAL    {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,vartot[0],vartot[1],vartot[2]))
        log.write('--------------------------------------\n')
        log.write(' {} ZONAL    {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,varzon[0],varzon[1],varzon[2]))
        log.write('--------------------------------------\n')
        log.write(' {} EDDY     {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,vared[0],vared[1],vared[2]))
        log.write('--------------------------------------\n')
        log.write(' {} EDDY(LW) {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,vared1[0],vared1[1],vared1[2]))
        log.write('--------------------------------------\n')
        log.write(' {} EDDY(SW) {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,vared2[0],vared2[1],vared2[2]))
        log.write('--------------------------------------\n')
        log.write(' {} EDDY(KW) {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,vared3[0],vared3[1],vared3[2]))
        log.write('--------------------------------------\n')
        
    
    def table_conv(self, varin,ntp,name,log):
        
        fac=1e5
        varin=fac*varin
        varzon = varin[:,0]
        vared  = np.nansum(varin[:,1:ntp-1],axis=1)
        vared1 = np.nansum(varin[:,1:nw1-1],axis=1)
        vared2 = np.nansum(varin[:,nw1:nw2-1],axis=1)
        vared3 = np.nansum(varin[:,nw2:nw3-1],axis=1)
        vartot = varzon+vared
        
        log.write(' {} TOTAL    {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,vartot[0],vartot[1],vartot[2]))
        log.write('--------------------------------------\n')
        log.write(' {} ZONAL    {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,varzon[0],varzon[1],varzon[2]))
        log.write('--------------------------------------\n')
        log.write(' {} EDDY     {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,vared[0],vared[1],vared[2]))
        log.write('-------------------------------------\n')
        log.write(' {} EDDY(LW) {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,vared1[0],vared1[1],vared1[2]))
        log.write('--------------------------------------\n')
        log.write(' {} EDDY(SW) {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,vared2[0],vared2[1],vared2[2]))
        log.write('--------------------------------------\n')
        log.write(' {} EDDY(KW) {: 4.3f}  {: 4.3f}  {: 4.3f}\n'.format(name,vared3[0],vared3[1],vared3[2]))
        log.write('--------------------------------------\n')   
    
    
    
    def varatts(self, w_nc_var,varname,tres,vres):
        if tres == 0:
            tatt= u"Daily\nM"
        elif tres == 1:
            tatt= u"Annual mean\nM"
        
        if vres == 0:
            vatt= u"Pressure levels\n"
        elif vres == 1:
            vatt = u"Vertically integrated\n"
        
        if varname == 'a':
            w_nc_var.setncatts({'long_name': u"Available Potential Energy",'units': u"W m-2", 'level_desc': vatt,\
                                'var_desc': u"APE -> KE",\
                                'statistic': tatt})
        elif varname == 'ek':
            w_nc_var.setncatts({'long_name': u"Kinetic Energy",'units': u"W m-2", 'level_desc': vatt,\
                                'var_desc': u"APE -> KE",\
                                'statistic': tatt})
        elif varname == 'a2k':
            w_nc_var.setncatts({'long_name': u"Conversion between APE and KE",'units': u"W m-2", 'level_desc': vatt,\
                                'var_desc': u"APE <-> KE",\
                                'statistic': tatt})
        elif varname == 'k':
            w_nc_var.setncatts({'long_name': u"Kinetic Energy",'units': u"W m-2", 'level_desc': vatt,\
                                'var_desc': u"APE -> KE",\
                                'statistic': tatt})
