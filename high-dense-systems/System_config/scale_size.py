import argparse
import numpy as np
from collections import deque
from itertools import tee
from math import acos


class ReadDumpFile(object):

    def __init__(self,file_name,n_atoms,scaling_parameter):
        self.filename = file_name
        self.natoms = n_atoms
        self.scaling = scaling_parameter
        
    
    def getdata(self):
        """
        Read trajectories from an ovito file
        self.filename, str that contains the filename to be read
        n_atoms, int that gives the number of particles to be read
        returns the new scaled coordinates and sizes
        """
        
        f = open(self.filename,'r')
        print("this is for", self.natoms, "particles")
        
        counter = 0
        time = 0
        
        id = []
        type = []
        r = []
        q = []
        axis = []
        baxis = []
        caxis = []
        timelist = []
        s = self.scaling
        
        """
        Read given ovito file with particle number = self.natoms
        """
        
        for line in f:
            counter += 1
            if counter == 2:
                time = float(line.strip())
                timelist.append(time)
                if time == 0.:
                    writetype = 'w'
                else:
                    writetype = 'a+'
                print(time)
                OvitoFile = open(self.filename + '-scaled.dump', writetype)
                OvitoFile.write('ITEM: TIMESTEP \n')
                OvitoFile.write(str(time) + '\n')
            if counter == 4:
                self.natoms = int(line.strip())
                OvitoFile.write('ITEM: NUMBER OF ATOMS \n')
                OvitoFile.write(str(self.natoms) + '\n')
            if counter == 6:
                L = line.strip()
                L= s*float(L.split()[1])
                OvitoFile.write('ITEM: BOX BOUNDS pp pp pp \n')
                OvitoFile.write(str(0) + '\t' + str(L) + '\n')
                OvitoFile.write(str(0) + '\t' + str(L) + '\n')
                OvitoFile.write(str(0) + '\t' + str(L) + '\n')
                OvitoFile.write('ITEM: ATOMS id type x y z c_orient[1] c_orient[2] c_orient[3] c_orient[4] c_shape[1] c_shape[2] c_shape[3] \n')
                """
                Writing format for quarternion px,py,pz,s
                """

            if counter > 9: #number of header lines
                line = line.strip()
                columns = line.split()
                id.append(int(columns[0]))
                type.append(int(columns[1]))
                r.append(np.array([s*float(columns[2]),s*float(columns[3]),s*float(columns[4])]))
                q.append(np.array([float(columns[5]),float(columns[6]),float(columns[7]),float(columns[8])]))
                axis.append(s*float(columns[9]))
                baxis.append(s*float(columns[10]))
                caxis.append(s*float(columns[11]))
            if counter > self.natoms + 8 : #number of lines in one time step
                # start computing stuff

                traj = zip(id,r,q,axis,baxis,caxis)
                
                for (id_frame, r_frame, q_frame, axis_frame, baxis_frame, caxis_frame) in traj:
                    
                    OvitoFile.write(str(int(id_frame)) + '\t' + str(int(id_frame)) + '\t' + str(format(r_frame[0], '.4f')) + '\t' + str(format(r_frame[1], '.4f')) + '\t' + str(format(r_frame[2], '.4f')) + '\t' + str(format(q_frame[0], '.4f')) + '\t' + str(format(q_frame[1], '.4f')) + '\t' + str(format(q_frame[2], '.4f')) + '\t' + str(format(q_frame[3], '.4f')) + '\t' + str(format(axis_frame, '.4f')) + '\t' + str(format(baxis_frame, '.4f')) + '\t' + str(format(caxis_frame, '.4f')) + '\n')
                
                counter = 0
                id = []
                type = []
                r = []
                q = []
                axis = []
                baxis = []
                caxis = []
    

def main(filename,natoms,scaling):
    """
    Scales the ovito dump file to "args.scaling" its original size
    """
    
    # read the ovito file
    get_dump = ReadDumpFile(filename,natoms,scaling)
    get_dump.getdata()

if __name__ == '__main__':    
#    parser = argparse.ArgumentParser(description='Calculate mean square displacement and mean square angular displacement')
#    simgroup = parser.add_argument_group('Simulation settings')
#    simgroup.add_argument('--filename', metavar='filename', help='Ovito file to be read.', required=True)
#    simgroup.add_argument('--natoms', metavar="natoms", type=float, help='number of particles in the system', required=True)
#    simgroup.add_argument('--scaling',metavar="scaling", type=float,help='scaling parameter needed', required = True)
#    args = parser.parse_args()
    
#    main(args.filename,args.natoms,args.scaling)
    filename = '/data/scc/sophia/FreeBROWNIANParticle1-n1372-t1000-phi0.05-k0-tb0.01-ts0.001/Tails/Tail5.dump'
    natoms = 1372
    scaling_factor = 0.5
    main(filename, natoms, scaling_factor)
    
    
    
    