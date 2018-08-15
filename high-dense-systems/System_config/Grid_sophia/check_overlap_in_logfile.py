#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 19:03:08 2018

@author: sophia
"""



import argparse
import numpy as np

""" Funktioniert NUR für overlap-basierte Ovito files, bei denen auch wirklich die Kollision erfasst wurde! """

class ReadDumpFile(object):

    def __init__(self,file_name, brownian_timestep, saving_timestep, system_lifetime, packing_fraction, particle_number,aspect_ratio,simulation_id):
        self.filename = '/data/scc/sophia/SYSTEMERROR1-n1372-t10-phi0.4-k0-tb0.01-ts0.0001/Ellipse-Benchmark-n1372-k0.0-phi0.4-id%s-equilibrate.h5.dump' % (simulation_id)
        print(self.filename)
        #self.filename = file_name+"/Ellipse-Benchmark-n%s-k%s-phi%s-id%s-equilibrate.h5.dump" % (particle_number, aspect_ratio, packing_fraction, simulation_id)
        self.natoms = particle_number
    
    def getdata(self):
        """
        Read trajectories from an ovito file
        self.filename, str that contains the filename to be read
        n_atoms, int that gives the number of particles to be read
        returns the following
        L: str, size of the simulation box
        timelist: list of times that appears in the ovito file
        rbig: list of particle position per time in timelist
        qbig: list of particle quaternion per time in timelist 
        
        """
        
        f = open(self.filename,'r')
        print(self.filename)
        print("this is for", self.natoms, "particles")
        
        counter = 0
        time = 0
        
        id = []
        type = []
        r = []
        q = []
        
        r_init = []
        q_init = []
        
        rbig = []
        qbig = []
        timelist = []
        
        """
        Read given ovito file with particle number = self.natoms
        """
        
        for line in f:
            counter += 1
            if counter == 2:
                time = float(line.strip())
                timelist.append(time)
#                 print(time)
            if counter == 4:
                self.natoms = int(line.strip())
            if counter == 6:
                L = line.strip()
                L = float(L.split()[1])
            if counter > 9 and counter: #number of header lines
                line = line.strip()
                columns = line.split()
                id.append(int(columns[0]))
                type.append(int(columns[1]))
                r.append(np.array([float(columns[2]),float(columns[3]),float(columns[4])]))
#                 q.append(np.array([float(columns[5]),float(columns[6]),float(columns[7]),float(columns[8])]))
#                 aaxis = float(columns[9])
#                 baxis = float(columns[10])
#                 caxis = float(columns[11])
            if counter > self.natoms + 8 : #number of lines in one time step
                # start computing stuff
             
                if rbig == []:
                    r_init = r
                    q_init = q
                    
                rbig.append(r)
                qbig.append(q)
                
                counter = 0
                id = []
                type = []
                r = []
                q = []  

        return(L,timelist,rbig,qbig)

def calculate_norm(rbig, particle_number):
    overlaps = []
    for frame in rbig:
        for particle_id in range(particle_number):
            for neighbor_particle in range(particle_number):
                distance = np.linalg.norm((frame[particle_id]- frame[neighbor_particle]))
                if distance < 2 and neighbor_particle != particle_id :
                    overlaps.append(distance)


    return overlaps


if __name__ == '__main__':    
#    parser = argparse.ArgumentParser(description='Perform a Brownian simulation of the Woetzel benchmark ellipse system.')
#    simgroup = parser.add_argument_group('Simulation settings')
#    simgroup.add_argument('--path', metavar='path', help='Path to save simulation output to.', required=True)
#    simgroup.add_argument('--brownian-timestep', metavar="t_B", type=float, help='Brownian time step. If not set use Newtonian dynamics.', required=True)
#    simgroup.add_argument('--saving-timestep', metavar="t_s", type=float, help='Timestep for saving the system summary.', required=True)
#    simgroup.add_argument('--system-lifetime', metavar="T", type=float, help='Simulation duration.', required=True)
#    simgroup.add_argument('--packing-fraction', metavar="phi", type=float, help='Packing fraction of the system between 0 and pi/4.', required=True)
#    simgroup.add_argument('--particle-number', metavar="n", type=int, help='Number of particles in the system. Must be a square number.', required=True)
#    simgroup.add_argument('--aspect-ratio', metavar="k", type=float, help='Aspect ratio of the ellipses', required=True)
#    simgroup.add_argument('--simulation_id', metavar="id", type=int, help='Simulation number, useful for array jobs. Default is 0.', default=0)
#    simgroup.add_argument('--q-vector',metavar="q_vec", type=float, help='Q-space vector where incoherent density corr is evaluated',required = True)
#    parser.add_argument('--debug', action='store_true', help='print debug messages to stderr')
#    parser.add_argument('--verbose', help='Print simulation status to stdout.', action='store_true')
#    args = parser.parse_args()
 
    path = '/data/scc/sophia/Figure4_2_StatStrucFact_MSDCalculations1-n1372-t2000-phi0.494-k0-tb0.01-ts5'
    particle_number = 1372
    for i in range(4):
        k =ReadDumpFile(path, 0.01, 5, 2000, 0.4, particle_number, 0.0, (i+1)) #args.simulation_id
        # read the normal ovito file
        [L,timelist,rbig,qbig] = k.getdata()
        OVERLAPS = calculate_norm(rbig, particle_number)
        if OVERLAPS != []:
            print('minimum distance', np.amin(OVERLAPS))
        print('and all detected overlaps give \n', OVERLAPS)
        
    
    