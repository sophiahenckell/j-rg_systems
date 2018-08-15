import argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import deque
from itertools import tee
from math import acos
import glob


class ReadDumpFile(object):

    def __init__(self,file_name1,file_name2,n_atoms):
        self.filename1 = file_name1
        self.filename2 = file_name2
        self.natoms = n_atoms
        
    def getfreedata(self):
        """
        Read trajectories from an cumulative position file with extension "free"
        self.filename, str that contains the filename to be read
        n_atoms, int that gives the number of particles to be read
        returns the following
        L: str, size of the simulation box
        timelist: list of times that appears in the ovito file
        rbig: list of particle position per time in timelist
        qbig: list of particle quaternion per time in timelist 
        
        """
        
        f = open(self.filename1,'r')
        
        counter = 0
        time = 0
        
        id = []
        r = []
        r_init = []
        rbig = []
        timelist = []
        
        """
        Read given ovito file with particle number = self.natoms
        """
        
        for line in f:
            counter += 1
            if counter == 2:
                time = float(line.strip())
                timelist.append(time)
            if counter > 4: #number of header lines
                line = line.strip()
                columns = line.split()
                id.append(int(columns[0]))
                r.append(np.array([float(columns[1]),float(columns[2]),float(columns[3])]))
            if counter > self.natoms + 2 : #number of lines in one time step
                # start computing stuff
             
                if rbig == []:
                    r_init = r
                rbig.append(r)
                counter = 0
                id = []
                r = []
                
        return(timelist,rbig)  
    
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
        
        f = open(self.filename2,'r')
        print("reading files for", self.natoms, "particles")
        
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
                q.append(np.array([float(columns[5]),float(columns[6]),float(columns[7]),float(columns[8])]))
                aaxis = float(columns[9])
                baxis = float(columns[10])
                caxis = float(columns[11])
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
    

class ConvertTo(object):
    # Needed if you want to get the body coordinates of the ellipsoids
    
    def __init__(self,vector,quaternion):
        self.vec = vector
        self.quat = quaternion
        
    def space_vector_norm(self):
        # convert to spatial coordinates
        
        [q0, q1, q2, q3] = self.quat
        [v0, v1, v2] = self.vec
        
        # corresponding rotation matrix of the quaternion self.quat
        A = np.array([[q3*q3 + q0*q0 - q1*q1 - q2*q2, 2.*(q0*q1 + q3*q2),2.*(q0*q2 - q3*q1)], [2.*(q0*q1 - q3*q2), q3*q3 - q0*q0 + q1*q1 - q2*q2, 2.*(q1*q2 + q3*q0)],[2.*(q0*q2 + q3*q1), 2.*(q1*q2 - q3*q0), q3*q3 - q0*q0 - q1*q1 + q2*q2]])
        AT = np.transpose(A)
        
        space_vector = np.matmul(AT,self.vec)
        norm = np.linalg.norm(space_vector)
        
        if norm != 0:
            space_vector = space_vector/norm
        
        return space_vector
    
    def body_vector_norm(self):
        # convert to body coordinates
        
        [q0, q1, q2, q3] = self.quat
        [v0, v1, v2] = self.vec
        
        # corresponding rotation matrix of the quaternion self.quat
        A = np.array([[q3*q3 + q0*q0 - q1*q1 - q2*q2, 2.*(q0*q1 + q3*q2),2.*(q0*q2 - q3*q1)], [2.*(q0*q1 - q3*q2), q3*q3 - q0*q0 + q1*q1 - q2*q2, 2.*(q1*q2 + q3*q0)],[2.*(q0*q2 + q3*q1), 2.*(q1*q2 - q3*q0), q3*q3 - q0*q0 - q1*q1 + q2*q2]])
        
        body_vector = np.matmul(A,self.vec)
        norm = np.linalg.norm(body_vector)
        
        if norm != 0:
            body_vector = body_vector/norm
        
        return body_vector
    
    def space_vector(self):
        # convert to spatial coordinates
        
        [q0, q1, q2, q3] = self.quat
        [v0, v1, v2] = self.vec
        
        # corresponding rotation matrix of the quaternion self.quat
        A = np.array([[q3*q3 + q0*q0 - q1*q1 - q2*q2, 2.*(q0*q1 + q3*q2),2.*(q0*q2 - q3*q1)], [2.*(q0*q1 - q3*q2), q3*q3 - q0*q0 + q1*q1 - q2*q2, 2.*(q1*q2 + q3*q0)],[2.*(q0*q2 + q3*q1), 2.*(q1*q2 - q3*q0), q3*q3 - q0*q0 - q1*q1 + q2*q2]])
        AT = np.transpose(A)
        
        space_vector = np.matmul(AT,self.vec)
        norm = np.linalg.norm(space_vector)
        
        return space_vector
    
    def body_vector(self):
        # convert to body coordinates
        
        [q0, q1, q2, q3] = self.quat
        [v0, v1, v2] = self.vec
        
        # corresponding rotation matrix of the quaternion self.quat
        A = np.array([[q3*q3 + q0*q0 - q1*q1 - q2*q2, 2.*(q0*q1 + q3*q2),2.*(q0*q2 - q3*q1)], [2.*(q0*q1 - q3*q2), q3*q3 - q0*q0 + q1*q1 - q2*q2, 2.*(q1*q2 + q3*q0)],[2.*(q0*q2 + q3*q1), 2.*(q1*q2 - q3*q0), q3*q3 - q0*q0 - q1*q1 + q2*q2]])
        
        body_vector = np.matmul(A,self.vec)
        norm = np.linalg.norm(body_vector)

        
        return body_vector


class MeanSquare(object):
    
    def __init__(self,box_size,n_atoms,big_array_r,big_array_q):
        """
        self.L : float, size of the simulation box
        self.natoms: float, number of particles inside the simulation box
        self.bigarray: list, list of time frame with list positions of each particles
        """
        self.L = box_size
        self.natoms = n_atoms
        self.bigarrayr = big_array_r
        self.bigarrayq = big_array_q
        
    def calculate_body_displacement(self):
        """
        Calculate mean square displacement in body coordinates of file in ovito format
        Returns an array that corresponds to mean square disp along the a,b,c axes resp
        """
        
        displacement = []
        a,b = tee(self.bigarrayr) # get time sequence
        next(b,None)
        bigr = zip(a,b)
        bigq = deque(self.bigarrayq) # get orientation for body coordinate transformation
        bigq.popleft() # remove the first time frame because it would not be used

        #Calculate mean square displacement for different time frames
        for (timer0,timer1) in bigr:
            
            disp_per_particle = []
            for (particlet0,particlet1) in zip(timer0,timer1):
                dr = particlet1 - particlet0
                disp_per_particle.append(dr)
            displacement.append(disp_per_particle)
        
        displacement2 = [] # get body coordinates
        for (disp_time,q_time) in zip(displacement,bigq):
            disp_per_particle2 = []
            for (disp_particle,q_particle) in zip(disp_time,q_time):
                QQ = ConvertTo(disp_particle,q_particle)
                drbody = QQ.body_vector()
                disp_per_particle2.append(drbody)
            displacement2.append(disp_per_particle2)
             
        msd = []
        
        for timeframe in displacement2:

            msd_x = 0. # major axis
            msd_y = 0. # minor1 axis
            msd_z = 0. # minor2 axis
            for particledelta in timeframe:
                msd_x += particledelta[0]**2
                msd_y += particledelta[1]**2
                msd_z += particledelta[2]**2
                
            msd.append([msd_x/self.natoms,msd_y/self.natoms,msd_z/self.natoms])
            
        return msd
    
    def calculate_displacement(self):
        """
        Calculate mesouran square displacement of file in ovito format in space coordinates (lab frame)
        Returns 1d MSD in lab coordinates
        """
        
        displacement = []
        r_init = self.bigarrayr[0]


        #Calculate mean square displacement
        for time_frame in self.bigarrayr:
            disp_per_particle = []
             
            for (r0,rn) in zip(r_init,time_frame):
             
                z = rn - r0
                disp_per_particle.append(z)
                 
            displacement.append(disp_per_particle)
             
        msd = []
        for timeframe in displacement:
            particle_counter = -1
            msd_t = 0
            for particledelta in timeframe:
                particle_counter +=1
                msd_t += np.dot(particledelta,particledelta)
            msd.append([msd_t/self.natoms])
            
        return msd
    
    def calculate_angular_displacement(self):
        """
        Calculate mean square angular displacement
        returns 2d array for the a,b respectively
        """
        
        displacement_major = []
        displacement_minor1 = []
        displacement_minor2 = []
        # get initial quaternion
        q_init = self.bigarrayq[0]
        
        vec_minor1_init = []
        vec_minor2_init = []
        vec_major_init = []
        # from the quaternion, get the initial vector representation
        
        # minor axis body coordinate: [0,1,0] or [0,0,1]
        # major axis body coordinate: [1,0,0] 
        
        a = [1.,0.,0.]
        b = [0.,1.,0.]
        c = [0.,0.,1.]
        
        for quat in q_init: # convert initial orientations to vector form
            orient_minor1 = ConvertTo(b,quat)
            orient_minor2 = ConvertTo(c,quat)
            orient_major = ConvertTo(a,quat)
            vec_minor1_init.append(orient_minor1.body_vector_norm())
            vec_minor2_init.append(orient_minor2.body_vector_norm())
            vec_major_init.append(orient_major.body_vector_norm())
            
        

        #Calculate all the vectors needed
        for time_frame in self.bigarrayq:
            disp_per_particle_major = []
            disp_per_particle_minor1 = []
            disp_per_particle_minor2 = []
            
            disp_major_sum = 0.
            disp_minor1_sum = 0.
            disp_minor2_sum = 0.
            # body coordinates corresponding to the major and minor axes

            # for each particle in each time frame, get the vector representation and get the angular 
            # difference wrt to the initial vector orientation
            for (theta0,thetan) in zip(vec_major_init,time_frame):
                
                orient_major=ConvertTo(a,thetan)
                delta_major = acos(round(np.dot(theta0,orient_major.body_vector_norm()),12))
                disp_major_sum += delta_major
                
            disp_major_sum = disp_major_sum/self.natoms
#                 disp_per_particle_major.append(delta_major)
                
            for (theta0,thetan) in zip(vec_minor1_init,time_frame):
                
                
                orient_minor1=ConvertTo(b,thetan)
                delta_minor1 = acos(round(np.dot(theta0,orient_minor1.body_vector_norm()),12))
                disp_minor1_sum += delta_minor1
#                 disp_per_particle_minor1.append(delta_minor1)
            disp_minor1_sum = disp_minor1_sum/self.natoms
                
            for (theta0,thetan) in zip(vec_minor2_init,time_frame):
                
                
                orient_minor2 = ConvertTo(c,thetan)
                delta_minor2 = acos(round(np.dot(theta0,orient_minor2.body_vector_norm()),12))
                disp_minor2_sum += delta_minor2
#                 disp_per_particle_minor2.append(delta_minor2)

            disp_minor2_sum = disp_minor2_sum/self.natoms
                
                
#             displacement_major.append(disp_per_particle_major)
#             displacement_minor1.append(disp_per_particle_minor1)
#             displacement_minor2.append(disp_per_particle_minor2)
            displacement_major.append(disp_major_sum**2)
            displacement_minor1.append(disp_minor1_sum**2)
            displacement_minor2.append(disp_minor2_sum**2)
                
        msr = []
        
        for (tminor1,tminor2,tmajor) in zip(displacement_minor1,displacement_minor2,displacement_major):
             
            msr.append([tminor1,tminor2,tmajor])
        
        
        return msr
        

class AverageOverFiles(object):
    
    def __init__(self,ms_quantity,no_timeframes,no_components,no_files):
        self.mean_square_quantity = ms_quantity
        self.number_of_timeframes = no_timeframes
        self.number_of_components = no_components
        self.number_of_files = no_files
        self.mean_square_quantity_ave = np.zeros([self.number_of_timeframes,self.number_of_components])
        
    def averaging(self):
        
        for msdfile in  self.mean_square_quantity:    
            for i in range(self.number_of_timeframes):
                for j in range(self.number_of_components):
                    self.mean_square_quantity_ave[i][j] += msdfile[i][j]
                    
                
                
#             for i in range(self.number_of_timeframes): # for each time frame in the msd per file
#                 for j in range(self.number_of_components):
#                     self.mean_square_quantity_ave[i][j] += msdfile[i][j]
                        
        return self.mean_square_quantity_ave/self.number_of_files

#def calc_msd(file1, file2 ,natoms):
def calc_msd(source_directory, natoms):

    """
    Calculates the MSD and MSAD_a and MSAD_b given a directory with dump files
    source_director: string, path or location of the files to be read
    natoms: number of atoms in the file
    Assumptions: same time frames are saved for all files in this directory 
    """

    file_list_free = glob.glob(source_directory + '/*comulative.dump') # stores all cumulative position
    file_list_quat = glob.glob(source_directory + '/*h5.dump') # ovito file with quat information

    
    no_files = len(file_list_free)
    if no_files != len(file_list_quat):
        sys.exit("some free or dump files are missing")
    
    msd_body_list = []
    msd_list = []
    msr_list = []
    timesteps = []
    counter2 = 0
    
    for (filename1,filename2) in zip(file_list_free,file_list_quat):
        counter2 += 1
        get_dump = ReadDumpFile(filename1,filename2,natoms)
        # read free bc positions file
        [timelist0,rfree] = get_dump.getfreedata()
        # read the normal ovito file
        [L,timelist,rbig,qbig] = get_dump.getdata()
        

        # calculate mean square displacement and mean square angular displacement --> MSD_a, MSD_b, MSD_c
        calculate_msd = MeanSquare(L,natoms,rfree,qbig)
        # body frame translational displacement
        msd_body = calculate_msd.calculate_body_displacement()
        # lab frame translational displacement --> MSD
        msd = calculate_msd.calculate_displacement()
        # lab frame rotational displacement (is there any other kind ? ) --> MSR_a, MSR_b
        msr = calculate_msd.calculate_angular_displacement()
        
        timesteps.append(len(msd))
        timesteps.append(len(msd_body))
        timesteps.append(len(msr))
        msd_body_list.append(msd_body)
        msd_list.append(msd)
        msr_list.append(msr)
        
    print(timesteps)
    no_timeframes = min(timesteps)
    
    avefilesmsdbody = AverageOverFiles(msd_body_list,no_timeframes,3,no_files)
    msdbodyave = avefilesmsdbody.averaging()
    
    avefilesmsd = AverageOverFiles(msd_list,no_timeframes,1,no_files)
    msdave = avefilesmsd.averaging()
#     
    avefilesmsr = AverageOverFiles(msr_list,no_timeframes,3,no_files)
    msrave = avefilesmsr.averaging()

    
#     counter = 0
#     msd_ave = np.zeros(min(filecounter))
#     for msd in msd_list: # for each file msd in the big list
#         counter += 1 # count the number of files read
#          # initialize the array where the msd will be stored
#         for i in range(min(filecounter)): # for each time frame in the msd per file
#             msd_ave[i] += msd[i] # add them here
#               
#     msd_ave = msd_ave/counter # this is the average listed in time
#     
#     counter2 = 0
#     msr_ave = np.zeros([min(filecounter),2])
#     for msr in msr_list:
#         counter2 += 1
#         for i in range(min(filecounter)):
#             msr_ave[i][0] += msr[i][0]
#             msr_ave[i][1] += msr[i][1]
#             
#     msr_ave = msr_ave/counter2
    
            
            
    # write the results in the file
    msdbody_dat = open(source_directory +'/msdbody.dat','w')
    msd_dat = open(source_directory + '/msd.dat','w')
    msr_dat = open(source_directory +'/msr.dat','w')
     
    for (time,[msdx,msdy,msdz]) in zip(timelist,msdbodyave):
        msdbody_dat.write(str(format(time, '.12f')) + '\t'+ str(format(msdx, '.12f')) + '\t'+ str(format(msdy, '.12f'))+ '\t'+ str(format(msdz, '.12f')) + '\n')
     
    for (time,[msdt]) in zip(timelist,msdave):
        msd_dat.write(str(format(time, '.12f')) + '\t'+ str(format(msdt, '.12f')) + '\n')
#          
    for (time,msrt) in zip(timelist,msrave):
        msr_dat.write(str(format(time, '.12f')) + '\t'+ str(format(msrt[0], '.12f')) + '\t'+ str(format(msrt[1], '.12f')) + '\n')



def create_image(directory):
    data1 = np.loadtxt(directory+'/msd.dat', dtype=float, usecols= 0)
    data2 = np.loadtxt(directory + '/msd.dat', dtype=float, usecols=1)
    plt.plot(data1, data2, marker = 'o')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    return None


if __name__ == '__main__':    
    # parser = argparse.ArgumentParser(description='Calculate mean square displacement and mean square angular displacement')
    # simgroup = parser.add_argument_group('Simulation settings')
    # simgroup.add_argument('--source_directory', metavar='source_directory', help='Directories where files are located', required=True)
    # simgroup.add_argument('--natoms', metavar="natoms", type=float, help='number of particles in the system', required=True)
    # args = parser.parse_args()
    #
    #
    # calc_msd(args.source_directory,args.natoms)
    directory = '/home/newton/sophia/Desktop/curr_MASTER/Project/data/BrownianCircle_radiusReduced0.5'
    calc_msd(directory, 1372)
    create_image(directory)