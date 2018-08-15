#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 07:14:17 2018

@author: sophia
"""
#TODO 
# ist incoherent wirklich bei 7 ausgewertet?
# wo ist bei incoherent der cutoff?



import glob
import numpy as np
import matplotlib.pyplot as plt

####################################################################
# Erste Sortierphase: nur peakfreie kurven des stat struktur faktors
####################################################################

path = '/data/scc/sophia/Figure4_2_StatStrucFact_MSDCalculations1-n1372-t2000-phi0.494-k0-tb0.01-ts5/*.dump-static-sk_new.dat'
filenames = glob.glob(path)
counter = 0
counter_val=0
datay = {}
datax = {}
promising_systems={}
plt.figure(1)


for i in filenames:

    data2 = {counter: np.loadtxt(i, usecols = 1)}
    if any (t>12 for t in data2[counter]):
        pass
    else:
        good_system = {counter_val : i}
        promising_systems.update({counter: i})
        counter_val+=1
        
        data1 = {counter: np.loadtxt(i, usecols = 0)} 
        
        if counter == 0:
            sum = data2[counter]
        else:
            sum = np.add(sum, data2[counter])
        
        datay.update(data2)
        datax.update(data1)
        
    counter += 1
print('amount of detected Files ', counter) 
print('static melted systems ', counter_val)
values = sum/float(len(datax))
plt.plot(datax[0]*2 , values)
plt.xlabel('$ \mid q \mid$', fontsize = 15)
plt.ylabel('static structure factor', fontsize = 15)
plt.title('Static structure factor for $\phi = 0.494$')
plt.show()

#########################################################
# Zweite Sortierphase: incoherent structur factor decays?
#########################################################

path_incoherent = '/data/scc/sophia/Figure4_2_StatStrucFact_MSDCalculations1-n1372-t2000-phi0.494-k0-tb0.01-ts5/*.dump-incoherent-sk_new.dat' 
filenames_incoherent = glob.glob(path_incoherent)
counter_val2 = 0
counter = 0
data_inc_x = {}
data_inc_y = {}
suitable_files = {}
plt.figure(2)
for i in filenames_incoherent:
    data_inc_2 = {counter: np.loadtxt(i, usecols = 1)}
    #print(data_inc_2)
    if any ( j< 0.01  for j in data_inc_2[counter]):
        data_inc_1 = {counter: np.loadtxt(i, usecols = 0)}
        plt.plot(data_inc_1[counter], data_inc_2[counter])
        data_inc_x.update(data_inc_1)
        data_inc_y.update(data_inc_2)
        suitable_files.update({counter: i})
        counter_val2 += 1
    
    
    
    counter += 1
print('amount of detected files ', counter)
print('incoherent, promising ' , counter_val2)
plt.show()


###########################################################
# Ultimate Systemcheck
###########################################################
files = 0
Systems ={}
for id_number in range(counter):
    u1 = '/data/scc/sophia/Figure4_2_StatStrucFact_MSDCalculations1-n1372-t2000-phi0.494-k0-tb0.01-ts5/Ellipse-Benchmark-n1372-k0.0-phi0.494-id{}-equilibrate.h5.dump-static-sk_new.dat'.format(id_number)
    u2 = '/data/scc/sophia/Figure4_2_StatStrucFact_MSDCalculations1-n1372-t2000-phi0.494-k0-tb0.01-ts5/Ellipse-Benchmark-n1372-k0.0-phi0.494-id{}-equilibrate.h5.dump-incoherent-sk_new.dat'.format(id_number)
    if u1 in promising_systems.values():
        if u2 in suitable_files.values():
            another_system = {files: '/data/scc/sophia/Figure4_2_StatStrucFact_MSDCalculations1-n1372-t2000-phi0.494-k0-tb0.01-ts5/Ellipse-Benchmark-n1372-k0.0-phi0.494-id{}-equilibrate.h5.dump'.format(id_number)}
            files += 1
            Systems.update(another_system)
print('--> REALLY melted systems are ', files)