#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 10:58:00 2018

@author: Valerio Lembo
"""
from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

import sys
sys.path.append('./diag_scripts/aux/Thermodynamics/')
import fluxogram as fluxogram
#import diag_scripts.aux.Thermodynamics.fluxogram as fluxogram

def diagram(filen,azin,apz,asein,aps,atein,apt,as2ks,at2kt,kteout,kte,kseout,kse,kzout,kz,
            az2kz,az2at,az2as,as2at,kt2kz,kt2ks,ks2kz):
    
    fl = fluxogram.Fluxogram(1000, 1000, grid_size = 20)
    
    fl.add_storage("AZ", 600, 0, 0)
    fl.add_storage("ASE", 600, 0.75, 0.25)
    fl.add_storage("ATE", 600, 1.5, 0)
    fl.add_storage("KTE", 600, 1.5, 1.5)
    fl.add_storage("KSE", 600, 0.75,1.25)
    fl.add_storage("KZ", 600, 0, 1.5)
    fl.add_storage("AZ+", 0, 0, -1)
    fl.add_storage("ASE+", 0, 0.75,-1)
    fl.add_storage("ATE+", 0, 1.5, -1)
    fl.add_storage("KTE-", 0, 1.5, 2.5)
    fl.add_storage("KSE-", 0, 0.75, 2.5)
    fl.add_storage("KZ-", 0, 0, 2.5)
    
    fl.add_flux("A2KZ", fl.storages[5], fl.storages[0], 100)
    fl.add_flux("AE2AZ", fl.storages[0], fl.storages[2], 150)
    fl.add_flux("AE2AS", fl.storages[0], fl.storages[1], 60)
    fl.add_flux("AE2AT", fl.storages[1], fl.storages[2], 60)
    fl.add_flux("A2KS", fl.storages[1], fl.storages[4], 60)
    fl.add_flux("A2KT", fl.storages[2], fl.storages[3], 100)
    fl.add_flux("KE2KS", fl.storages[3], fl.storages[4], 60)
    fl.add_flux("KS2KZ", fl.storages[4], fl.storages[5], 60)
    fl.add_flux("KE2KZ", fl.storages[3], fl.storages[5], 150)
    fl.add_flux("AZ+", fl.storages[6], fl.storages[0], 60)
    fl.add_flux("ASE+", fl.storages[7], fl.storages[1], 60)
    fl.add_flux("ATE+", fl.storages[8], fl.storages[2], 60)
    fl.add_flux("KTE-", fl.storages[3], fl.storages[9], 60)
    fl.add_flux("KSE-", fl.storages[4], fl.storages[10], 60)
    fl.add_flux("KZ-", fl.storages[5], fl.storages[11], 60)
    
    
    fl.draw(filen,azin,apz,asein,aps,atein,apt,as2ks,at2kt,kteout,kte,kseout,
            kse,kzout,kz,az2kz,az2at,az2as,as2at,kt2kz,kt2ks,ks2kz)
    #fl.show()
    #plt.savefig('prova.pdf')
