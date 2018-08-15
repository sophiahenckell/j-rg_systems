#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 15:09:39 2018

@author: sophia
"""

##TODO: will swell up to a given pacion fraction, last if-loop doesnt work however
## flatten is not the best way hugh?

import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d 
#n = 32 # must be n^3 * 4
#box_size = 100. #wichtig! mit punkt!

def _generate_particles_FCC_spheres(n, Phi, box_size):
    
    
    #check if n is a suitable particle number
    check = round((n /4)**(1./3),3) #check if particle number is cubic
    
    if not check**3 * 4 == n:
       print ('no valid particle number! \n')
       sys.exit('no valid number')
       
    
    #fcc setup
    cells_per_row = check 
    shift = box_size / (cells_per_row * 2)
    
    base = np.linspace(start = 0, stop = box_size, num = cells_per_row, endpoint = False) 
    xx, yy, zz = np.meshgrid(base, base, base)
    
    fcc_1 = xx + shift
    fcc_2 = yy + shift
    fcc_3 = zz + shift
    
    fcc_x = np.concatenate([xx, fcc_1 , fcc_1, xx], axis = 0).flatten()
    fcc_y = np.concatenate([yy, fcc_2, yy , fcc_2], axis = 0).flatten()
    fcc_z = np.concatenate([zz, zz, fcc_3, fcc_3] , axis = 0).flatten()

    a= fcc_x.reshape(1,-1).T
    b= fcc_y.reshape(1,-1).T
    c= fcc_z.reshape(1,-1).T
    
    fcc_lattice = np.concatenate((a,b,c), axis=1)

    
    ################
    #check via Plot
    ################
    plt.figure()
    ax = plt.axes(projection= '3d')
    ax.scatter3D(fcc_1, fcc_2, zz, marker = 'o', color = 'c')
    ax.scatter3D(fcc_1, yy, fcc_3, marker = 'o', color = 'c')
    ax.scatter3D(xx, fcc_2, fcc_3, marker = 'o', color = 'c')
    ax.scatter3D(xx,yy,zz, marker = 'o', color = 'b')
    #ax.scatter3D(fcc_1, fcc_2, fcc_3, marker = 'o')
    plt.show()
    
    
    ############
    #WRITE FILE
    ############
    with open('AleenaFCCtestfile.txt', 'w') as outfile:
        outfile.write('ITEM: TIMESTEP \n')
        outfile.write(str(1)+ '\n')
        outfile.write('ITEM: NUMBER OF ATOMS pp pp pp \n')
        outfile.write(str(n)+ '\n')
        outfile.write('ITEM: BOX BOUNDS \n ')
        outfile.write("{a}\t{b}\n".format(a=0,b=box_size)*3)
        outfile.write('ITEM: ATOMS id type x y z \n')
        for i in range(len(fcc_x-1)):
            outfile.write('{w:} \t {w:}\t {x: .4f} \t {y: .4f} \t {z: .4f} \n'.format(w = i, x = fcc_x[i], y = fcc_y[i], z= fcc_z[i]))
    outfile.close()       
    #############
    # Inflate particles onto a given packing fraction
    ############
    #Phi = 0.75
    r = (Phi * 3 /(4 * np.pi * 4) )**(1./3) * (box_size/ cells_per_row)
    r_max = np.sqrt(2) * (box_size/cells_per_row) /4
    
    if r > r_max :
        print('ATTENTION: choose lower packing fraction')
    
    print('Radius is ' + str(r))    
_generate_particles_FCC_spheres(32,0.6
                                ,100)