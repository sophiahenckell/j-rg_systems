#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 14:53:10 2018

@author: sophia
"""

#from visual import shpere, colour
#
#count = 3
#r = 0.3
#
#for x in range(-count, count+1):
#    for y in range(-count, count+1):
#        for z in range(-count, coutn+1):
#            sphere(pos= [x,y,z], radius = r, color = color.yellow)

import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d 

def _generate_particles_FCC_spheres(self, n, k, box_size):
    '''
    diese Klasse soll kugeln gewünschter Dichte in einem Fcc Gitter erzeugen
    '''
    self.particle_number = n
    self.aspectratio = k
    self.box_size = box_size
    
    #first check if n is a suitable particle number
    
    
    #n = 32 # n must be n^3 * 4
    #box_size = 100. #wichtig! mit punkt!
    
    check = round((n /4)**(1./3),3) #check if particle number is cubic
    if not check**3 * 4 == n:
        print ('no valid particle number! \n')
        sys.exit('no valid number')

 


    cells_per_row = check 
    shift_for_fcc = box_size / (cells_per_row * 2)

    grid_space = np.linspace(start=0, stop= box_size, num = cells_per_row, endpoint = False )
    grid_points = [[float(x), 
                    float(y), 
                    float(z)] for x in grid_space for y in grid_space for z in grid_space]


    #grid_space_for_fcc = grid_space[:-1] ##letzte zeile verschiebt sich sonst außerhalb des Limits


    ## 1.) bcc lattice

    #grid_points_bcc = [[float(x)+shift_for_fcc, 
    #          float(y)+shift_for_fcc, 
    #          float(z)+shift_for_fcc
    #          ] for x in grid_space_for_fcc for y in grid_space_for_fcc for z in grid_space_for_fcc]
    #
    #
    ##plotting for correction
    #x,y,z= zip(*grid_points) #fürs plotting
    #x_1, x_2 , x_3 = zip(*grid_points_bcc)
    #plt.scatter(x,y, color= 'c')
    #plt.scatter(x_1, x_2, color = 'b')


    # 2.) FCC lattice
    #setup for plotting
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    #gridpoints
    grid_points_fcc1 = [[
                    float(x) + shift_for_fcc,
                    float(y) + shift_for_fcc,
                    float(z)] for x in grid_points for y in grid_points for z in grid_space]
    
    grid_points_fcc2 = [[
                    float(x) + shift_for_fcc,
                    float(y) ,
                    float(z) + shift_for_fcc] 
                    for x in grid_points for y in grid_space for z in grid_points]
    
    grid_points_fcc3 = [[
                    float(x) ,
                    float(y) + shift_for_fcc,
                    float(z) + shift_for_fcc] 
                    for x in grid_space for y in grid_points for z in grid_points]
    
    fcc_data = np.concatenate((grid_points_fcc1, grid_points_fcc2, grid_points_fcc3, grid_points    ), axis = 0)

    #plotting for reasons of oversight
    x,y,z = zip(*grid_points) #grundmuster
    x_2, y_2, z_2 = zip(*grid_points_fcc2)
    x_3, y_3, z_3 = zip(*grid_points_fcc3)
    x_1, y_1, z_1 = zip(*grid_points_fcc1)
    
    #plt.scatter(x,y, color= 'b')
    ax.scatter(x_1, y_1, z_1, color= 'c')
    ax.scatter(x_2, y_2, z_2, color= 'c')
    ax.scatter(x_3, y_3, z_3, color= 'c')
    ax.scatter(x,y,z, color = 'b')
    
    #x_json, y_json, z_json = zip(*k)
    
    with open('test.txt', 'w') as outfile:
        outfile.write('ITEM: TIMESTEP \n')
        outfile.write(str(1)+ '\n')
        outfile.write('ITEM: NUMBER OF ATOMS pp pp pp \n')
        outfile.write(str(len(fcc_data))+ '\n')
        outfile.write('ITEM: BOX BOUNDS \n ')
        outfile.write("{a}\t{b}\n".format(a=0,b=box_size)*3)
        outfile.write('ITEM: ATOMS id type x y z \n')
        for i in range(len(fcc_data-1)):
            outfile.write('{w:} \t {w:}\t {x: .4f} \t {y: .4f} \t {z: .4f} \n'.format(w = i, x =    fcc_data[i][0], y = fcc_data[i][1], z= fcc_data[i][2], t = 0))
            
            
        
        outfile.close




