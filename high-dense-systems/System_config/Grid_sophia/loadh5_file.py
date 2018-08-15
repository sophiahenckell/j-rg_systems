#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 13:08:46 2018

@author: sophia
"""
"""
funktioniert NUR bei Ovitofiles, die bei Collision abgespeichert wurden!!!
"""


import h5py
filename = '/data/scc/sophia/Figure4_2_StatStrucFact_MSDCalculations1-n1372-t2000-phi0.494-k0-tb0.01-ts5/Ellipse-Benchmark-n1372-k0.0-phi0.494-id12-equilibrate.h5'
f = h5py.File(filename, 'r')
print('Keys', f.keys())

#get the data
data = list(f[u'particles']['statics'])




