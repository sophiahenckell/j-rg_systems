import argparse
import numpy as np
from collections import deque
from itertools import tee
from math import acos, sqrt
from cmath import exp



class ReadDumpFile(object):

    def __init__(self,file_name,n_atoms):
        self.filename = file_name
        self.natoms = n_atoms
    
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
            if counter > 9: #number of header lines
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
    

class StructureFactor(object):
    
    def __init__(self,box_size,n_atoms,big_array_r,big_array_q):
        self.L = box_size
        self.natoms = n_atoms
        self.bigarrayr = big_array_r
        self.smallarrayr = self.bigarrayr[-1] # get the last part of the file for the calculation of static structure factor
        self.bigarrayq = big_array_q
        self.maxq = []
        self.max_sq_static = [] # we normalize the correlation function with the max value of the structure factor
        
    def calculate_incoherent_sk(self):
        """
        Iterate of positions of the particle r_j(t) - r_j(0) for the
        calculation of the calculation of the incoherent part of the 
        correlation function
        """
        
        displacement = []
        r_init = self.bigarrayr[0]

        print("computing incoherent part")
        #Calculate displacement between r_j(t) - r_j(0)
        for time_frame in self.bigarrayr:
            disp_per_particle = []
            # r0 = r(t=0), rn = r(t=tn)
            for (r0,rn) in zip(r_init,time_frame):
             
                z = rn - r0
                disp_per_particle.append(z)
                 
            displacement.append(disp_per_particle)
        
        number_of_q_vectors = len(self.maxq) # note: the average over a binsize is considered here
        fac = 2.*np.pi/self.L
        incoherent_part = []
        for timeframe in displacement:
            particle_counter = -1
            incoherent_part_t = 0
            # get the contribution of each atoms
            for particledelta in timeframe:
                particle_counter +=1
                # get the S(q) contribution for each q vector that resides in the chosen bin
                inc_part = (np.sum([ exp(-1j * np.dot(q,particledelta) * fac) for q in self.maxq ])) # np.sum(np.dot(q,particledelta))
                incoherent_part_t += np.real(inc_part*np.conj(inc_part))/self.natoms # normalize with number of atoms
            # get the average over all atoms
            incoherent_part.append(incoherent_part_t/number_of_q_vectors)
            
        norm_sk = max(incoherent_part)
            
        return incoherent_part/norm_sk
        
    def calculate_structure_factor_static(self):
        """
        Calculate static structure factor, averaged over different systems used
        wave vector: q = 2pi/L(qx,qy,qz) because of PBC
        """
        
        # define the q-space coodinates  in 3D
        # define qnx, qny, qnz = 1,2,3,...20
        q_radius = 70 # radius in q-space
        q_space = np.linspace(start=-q_radius, stop=q_radius, num=q_radius + 1, endpoint=True).tolist() # num+1 because endpoint counts as well...
        q_data = [ [float(qx),float(qy),float(qz)] for qx in q_space for qy in q_space for qz in q_space]
        # remove the origin [0.,0.,0] because this cant be resolved dued to PBC
        q_origin = q_data[int(len(q_data)/2.)]
        if not all(o == 0. for o in q_origin):
            sys.exit("deleting",q_origin)
        np.delete(q_data, int(len(q_data)/2.))
        
        # get the actual magnitude of q-space vectors
        fac = 2.*np.pi/self.L
        q_data_norm = [np.linalg.norm(qmag)*fac for qmag in q_data] 
        # define the bin resolution dq by increasing num_bins
        res = 3 # resolved for dq/2
        numbins = int(sqrt(3.*q_radius**2)) + 1
        #q_bin = 0,dq,2*dq,3*dq,...n*dq
        q_bin = np.linspace(start = 0, stop = numbins, num = numbins*res, endpoint=True).tolist()
        q_binreal = [ (q*fac - 0.001) for q in q_bin]

        # initialize the array where the static structure is to be stored
        rho_q_squared = np.zeros(len(q_bin) + 1)
        num_data_points = np.zeros(len(q_bin) + 1)
        
        bins = np.digitize(q_data_norm,q_binreal) # gives the bin index for each q-magnitude
        print('len',len(rho_q_squared),max(q_data_norm),max(q_binreal))

#         binlist = bins.tolist()
#         data_per_bins = [binlist.count(ind) for ind in binlist] # gives the number of q-vectors that resides in each bin, used for averaging
# 
#         
#         # create a dictionary with the "bin index number" as keys and "number of data to be averaged" as values
#         dic = {}
#         countdata = zip(binlist,data_per_bins)
# 
#         for (ind,data) in countdata:
#             if ind not in dic:
#                 dic[ind] = data
        
        print("computing sk_static")        
        # get structure factor for one system


        for (ind,q) in zip(bins,q_data):
            # calculate density for each vector q, summed for all positions r
            rho_q = (np.sum([ exp(-1j * np.dot(q,r) * 2*np.pi/self.L) for r in self.smallarrayr ]))
            # calculate static structure factor and bin it
            rho_q_squared[ind] += np.real(rho_q*np.conj(rho_q))/self.natoms
            num_data_points[ind] += 1
            
        
        skk = zip(q_binreal,rho_q_squared,num_data_points)
        count = 0
        #sk = []
        sk_find_max = []
        q_find_max = []
        for (q,rho,counter) in skk:
            if counter != 0. and q > fac:
                #sk.append((q,rho/counter))
                sk_find_max.append(rho/counter)
                q_find_max.append(q)
#             count += 1
#             if count in dic:
#                 sk.append((q,rho/dic[count])) # averaged over the number of data points added per bin
#                 sk_find_max.append(rho/dic[count])
#                 q_find_max.append(q)
#             else:
#                 print('----',q,rho)

        #calculate the corresponding bin index of the maximum value of the static structure factor

        self.max_sq_static.append(max(sk_find_max[0:int(len(sk_find_max)/2)])) # get value of maximum to be used in the correlation function calculation
        # we only take half of the actual sk because numerical discrepancies are present at long q values
        #qind = sk_find_max.index(self.max_sq_static)
        qind = 13
        print("index",qind,q_find_max[qind],max(sk_find_max))
        print("range",q_find_max[qind-2],q_find_max[qind],q_find_max[qind+2])


        self.max_sq_static.append(q_find_max[qind-2])
        self.max_sq_static.append(q_find_max[qind-1])
        self.max_sq_static.append(q_find_max[qind+1])
        self.max_sq_static.append(q_find_max[qind+2])


        
        
        # find the q vectors that reside in this bin range
        for q in q_data:
            s = np.linalg.norm(q)*fac
            if (s >= q_binreal[qind-2]) and (s <= q_binreal[qind+2]): # get two-four bins before and after to set the range
                self.maxq.append(q)
                
        print("qs",self.maxq)
        
        return zip(q_find_max,sk_find_max)

def calc_sk(path, brownian_timestep, saving_timestep, system_lifetime, packing_fraction, particle_number, aspect_ratio, simulation_id, verbose):
    """
    calculate structure factor
    """
    filename = "%s/Ellipse-Benchmark-n%s-k%s-phi%s-id%s-equilibrate.h5.dump" % (path, particle_number, aspect_ratio, packing_fraction, simulation_id)

    get_dump = ReadDumpFile(filename,particle_number)
    # read the normal ovito file
    [L,timelist,rbig,qbig] = get_dump.getdata()
            
    # calculate the ff:
    calculate_sk = StructureFactor(L,particle_number,rbig,qbig)
    # 1. static structure factor
    sk_static = calculate_sk.calculate_structure_factor_static()
    # 2. incoherent part of the structure factor
    sk_incoherent = calculate_sk.calculate_incoherent_sk()

    # write the results in the file
    sk_static_dat = open(filename+'-static-sk.dat','w')
    sk_incoherent_part_dat = open(filename+'-incoherent-sk.dat','w')
    
    for (q,sk) in sk_static:
        sk_static_dat.write(str(format(q, '.12f')) + '\t'+ str(format(sk, '.12f')) + '\n')
        
    for (time,sk_i) in zip(timelist,sk_incoherent):
        sk_incoherent_part_dat.write(str(format(time, '.12f')) + '\t'+ str(format(sk_i, '.12f')) + '\n')



if __name__ == '__main__':    
    # parser = argparse.ArgumentParser(description='Calculates structure factor and Density correlation function.')
    # simgroup = parser.add_argument_group('Simulation settings')
    # simgroup.add_argument('--path', metavar='path', help='source and saving directory. Path to save simulation output to.', required=True)
    # simgroup.add_argument('--brownian-timestep', metavar="t_B", type=float, help='Brownian time step. If not set use Newtonian dynamics.', required=True)
    # simgroup.add_argument('--saving-timestep', metavar="t_s", type=float, help='Timestep for saving the system summary.', required=True)
    # simgroup.add_argument('--system-lifetime', metavar="T", type=float, help='Simulation duration.', required=True)
    # simgroup.add_argument('--packing-fraction', metavar="phi", type=float, help='Packing fraction of the system between 0 and pi/4.', required=True)
    # simgroup.add_argument('--particle-number', metavar="n", type=int, help='Number of particles in the system. Must be a square number.', required=True)
    # simgroup.add_argument('--aspect-ratio', metavar="k", type=float, help='Aspect ratio of the ellipses', required=True)
    # simgroup.add_argument('--simulation-id', metavar="id", type=int, help='Simulation number, useful for array jobs. Default is 0.', default=0)
    # parser.add_argument('--debug', action='store_true', help='print debug messages to stderr')
    # parser.add_argument('--verbose', help='Print simulation status to stdout.', action='store_true')
    # args = parser.parse_args()
    #
    # calc_sk(args.path, args.brownian_timestep, args.saving_timestep, args.system_lifetime, args.packing_fraction, args.particle_number, args.aspect_ratio, args.simulation_id, args.verbose)

     #calc_sk('/data/scc/sophia/Second_MSDCalculations2-n500-t4000-phi0.494-k1.2-tb0.01-ts10', 0.01, 10, 4000, 0.494, 500, 1.2, 5 , 'store_true')
     for i in range(149):
        
        calc_sk('/data/scc/sophia/Figure4_2_StatStrucFact_MSDCalculations1-n1372-t2000-phi0.494-k0-tb0.01-ts5', 0.01, 5, 2000, 0.494, 1372, 0.0, (i+201), 'store_true')