# checked -A
"""
Generate a trivial 2x2 ellipse grid
"""
import sys
import math
import cmath
from collections import deque
from functools import lru_cache

import numpy as np
from scipy.spatial import ConvexHull

from bedsim.files import BedfileWriter
from math import sqrt

import itertools
from mpmath import norm

"""
> import sysgen.ellipse_grid
> a = sysgen.ellipse_grid.EllipseGrid()
> a.save_to_file("abc.h5")
"""


class SystemGenerator(object):
    def generate(self, n, **kwargs):
        raise NotImplementedError()

    def save_to_file(self, filename):  ### FIXME: USE files.py API HERE!!!
        f = BedfileWriter(filename, 'ShortFormat', None)
        f.system.system_properties = {'localtime': 0}
        f.system.boundary = self.boundary_data
        f.system.grid = self.grid_data
        f.particles.save_statics(self.particle_data)
        #print('Check this: \n self.particle_data = {}'.format(self.particle_data)) ## SOPHIA JUNI

    def __init__(self):
        self.particle_data = None
        self.grid_data = []
        self.boundary_data = None


class EllipseBenchmark(SystemGenerator):
    def generate(self, n, **kwargs):
        """
        Edited by Aiyin
        Generate the system with n particles and 'PeriodicBox' boundary conditions.
        @param n Generate n particles
        @param k aspect ratio of the ellipses
        @param phi packing fraction of the ellipses
        """
        k = kwargs['k']
        phi = kwargs['phi']
        if phi > math.pi / 6:
            raise ValueError("Packing fraction phi=%s is too big! Maximum is pi/4." % phi)
            # box_size = n**(1./3.)*2*math.sqrt(k)
        box_size = ((n * 4. / 3.) * np.pi * major ** 3 / phi) ** (1. / 3.)
        print("boxsize, generate", box_size)
        cells_per_row = math.ceil(n ** (1. / 3))
        cell_size = box_size / cells_per_row

        self.__generate_particles(n, k, box_size)

        # generate the cellspace grid
        grid_space = np.linspace(start=0, stop=box_size, num=cells_per_row + 1,
                                 endpoint=True).tolist()  # num+1 because endpoint counts as well...
        self.grid_data = [[float(x), float(y), float(y)] for x in grid_space for y in grid_space for z in grid_space]

        # boundary conditions
        self.boundary_data = "PeriodicBox"

    def __generate_particles(self, n, k, box_size):
        # particle_data => Ellipse# => {position: [x,y], velocity=[vx, vy], ...}
        major = 1.
        minor = 1 / major  # we let two axes be equal to each other (for now)
        minor2 = minor

        print("boxsize generate particles", box_size)

        cubeN = math.ceil(n ** (1. / 3.))
        self.particle_data = []
        for x in range(cubeN):
            for y in range(cubeN):
                for z in range(cubeN):
                    counter += 1
                    self.particle_data.append({"type": "Ellipse", "id": counter,
                                               "position": np.array([2. * major * x, 2. * major * y, 2. * major * z],
                                                                    dtype=np.float64),
                                               "velocity": np.array([0., 0., 0.]), "angle": np.float64(math.pi / 4),
                                               "angvel": np.array([0., 0., 0.]), "major": np.float64(major),
                                               "minor": np.float64(minor), "minor2": np.float64(minor2),
                                               "pinned": False})


class EllipseBenchmarkRandom(EllipseBenchmark):
    def generate(self, n, **kwargs):
        """
        Edited by Aiyin
        Generate the system with n particles and 'PeriodicBox' boundary conditions.
        @param n Generate n particles
        @param k aspect ratio of the ellipses
        @param phi packing fraction of the ellipses
        """
        k = kwargs['k']
        phi = kwargs['phi']
        radius = kwargs['radius']
        filename = kwargs.get('filename', None)
        if phi > 0.2:
            raise ValueError("Packing fraction phi=%s is too big! Maximum is pi/4." % phi)

        box_size = ((n * 4. / 3.) * np.pi * radius ** 3 / phi) ** (1. / 3.)
        print("boxsize generate2", box_size, "volume fraction", phi)
        cells_per_row = math.ceil(n ** (1. / 3.))
        print('cells_per_row', cells_per_row)
        cell_size = box_size / cells_per_row

        if filename is None:
            self.__generate_particles(n, k, box_size, radius)
        else:
            self.__generate_particles_from_file(n, radius,  filename)


        # generate the cellspace grid
        grid_space = np.linspace(start=0, stop=box_size, num=cells_per_row + 1,
                                 endpoint=True).tolist()  # num+1 because endpoint counts as well...
        self.grid_data = [[float(x), float(y), float(z)] for x in grid_space for y in grid_space for z in grid_space]

        # boundary conditions
        self.boundary_data = "PeriodicBox"

    def __generate_particles(self, n, k, box_size, radius):
        # particle_data => Ellipse# => {position: [x,y], velocity=[vx, vy], ...}


        print("ELLIPSE GENERATE: box size", box_size, ", volume fraction", ((n * 4. / 3.) * np.pi * radius ** 3) / (box_size ** 3))
        counter = 0
        self.particle_data = []

        r0 = np.random.uniform(0, box_size, 3)
        rlist = [r0]
        rijlist2 = []

        while len(rlist) < n:
            i = np.random.uniform(0, box_size, 3)
            rijlist = []

            for j in rlist:
                rij = np.absolute(i - j)
                np.copyto(rij, box_size - rij, where=rij > 0.5 * box_size)
                rijlist.append(np.linalg.norm(rij))

            rijlist = np.array(rijlist)
            if (rijlist > 2. * radius).all():
                rlist.append(i)
                rijlist2.append(rijlist)

        now = 0
        ids = []
        posits = []
        angles = []

        #         self.system.boundary.unwrap(np.random.uniform(0,cubeN*2,3))

        for i in range(n):
            counter += 1
            ids.append(counter)
            angle = np.random.uniform(0, 1, 4)
            norm = np.linalg.norm(angle)
            angle = angle / norm
            angles.append(angle)
            veloc = np.zeros(3)
            angveloc = np.zeros(3)
            posit = rlist.pop()
            posits.append(posit)

            self.particle_data.append(
                {"type": "SphereOrientable", "id": i, "position": posits[i], "velocity": [0, 0, 0],
                 "angle": [0, 0, 0, 0],
                 "angvel": [0.5, 0.5, 0], "major": np.float64(1), "minor": np.float64(1),
                 "minor2": np.float64(1), "pinned": False, "radius": radius})

    def __generate_particles_from_file(self, n, radius, filename):
        """
        Reads the last frame of an ovito file format - A
        """
        print('reading file')
        f = open(filename, 'r')

        counter = 0
        ids = []
        type = []
        posits = []
        angles = []
        self.particle_data = []

        for line in f:
            counter += 1
            if counter == 2:
                time = float(line.strip())
                print("time", time)
            if counter == 4:
                natoms = int(line.strip())
                if n != natoms:
                    sys.exit('no of atoms do not match')
            if counter == 6:
                L = line.strip()
                L = float(L.split()[1])
            if counter > 9:  # number of header lines
                line = line.strip()
                columns = line.split()
                ids.append(int(columns[0]))
                type.append(int(columns[1]))
                posits.append(np.array([float(columns[2]), float(columns[3]), float(columns[4])]))
                angles.append(np.array([float(columns[5]), float(columns[6]), float(columns[7]), float(columns[8])]))
                major = float(columns[9])
                minor = float(columns[10])
                minor2 = float(columns[11])

                self.particle_data.append(
                    {"type": "SphereOrientable",
                     "id": counter-10 ,
                     "position": posits[counter-10],
                     "velocity": [0, 0, 0],
                     "angle": angles[counter-10],
                     "angvel": [0.5, 0.5, 0],
                     "major": np.float64(1),
                     "minor": np.float64(1),
                     "minor2": np.float64(1),
                     "pinned": False, 'radius': radius})
        box_size = L
        cubeN = math.ceil(n ** (1. / 3.))
        #print("cubeM", cubeN, "box size", box_size, "volume fraction",
         #     ((n * 4. / 3.) * np.pi * major * minor * minor2) / (box_size ** 3), 'particle number', counter-9)
        print('box size', box_size, ', volume fraction after reading file ', (n * 4./3.)* np.pi * float(radius**3)/ (box_size ** 3) )
        # TODO Sophia


        #Uncomment the following lines if you need to check the read positionss
        """
        Uncomment the following lines if you need to check the read positionss
        OvitoFile = open('./data/ellipsoids.dump', 'w')
        OvitoFile.write('ITEM: TIMESTEP \n')
        OvitoFile.write(str(time) + '\n')
        OvitoFile.write('ITEM: NUMBER OF ATOMS \n')
        OvitoFile.write(str(n) + '\n')
        OvitoFile.write('ITEM: BOX BOUNDS pp pp pp \n')
        OvitoFile.write(str(0) + '\t' + str(box_size) + '\n')
        OvitoFile.write(str(0) + '\t' + str(box_size) + '\n')
        OvitoFile.write(str(0) + '\t' + str(box_size) + '\n')
        OvitoFile.write('ITEM: ATOMS id type x y z c_orient[1] c_orient[2] c_orient[3] c_orient[4] c_shape[1] c_shape[2] c_shape[3] \n')
        """

        # Writing format for quarternion px,py,pz,s

        """
        for i in range(n):
            OvitoFile.write(str(int(ids[i])) + '\t' + str(int(ids[i])) + '\t' + str(format(posits[i][0], '.4f')) + '\t' + str(format(posits[i][1], '.4f')) + '\t' + str(format(posits[i][2], '.4f')) + '\t' + str(format(angles[i][0], '.4f')) + '\t' + str(format(angles[i][1], '.4f')) + '\t' + str(format(angles[i][2], '.4f')) + '\t' + str(format(angles[i][3], '.4f')) + '\t' + str(format(major, '.4f')) + '\t' + str(format(minor, '.4f')) + '\t' + str(format(minor, '.4f')) + '\n')
            self.particle_data.append( {"type": "Ellipse", "id": ids[i] ,"position": posits[i], "velocity": np.array([0.,0.,0.]), "angle": angles[i], "angvel": np.array([0.,0.,0.]), "major": np.float64(major), "minor": np.float64(minor), "minor2": np.float64(minor2), "pinned": False} )
        """


class FCCLatticeGenerator(SystemGenerator): # added by SOPHIA
    def generate(self, n, **kwargs):
        """
        Edited by Sophia
        Generate the system with n particles and 'PeriodicBox' boundary conditions.
        @param n Generate n particles
        @param phi packing fraction of the ellipses
        """

        filename = kwargs.get('filename', None)
        phi = kwargs['phi']
        radius = kwargs['radius']
        simulation_id = kwargs['simulation_id']
        box_size = ((n * 4. / 3.) * np.pi * radius ** 3 / phi) ** (1. / 3.)
        print("FCC GENERATE: box size", box_size, "volume fraction", phi)
        cells_per_row = int(math.floor(float(box_size) / 2 / radius))  # TODO: ursprünglich ceil gewesen (Sophia)
        print("cells_per_row",cells_per_row,box_size/cells_per_row)
        if filename is None:
            self.__generate_particles_FCC_SOPHIA(n, phi, box_size, radius,simulation_id)
        else:
            self.__generate_particles_from_file(n, radius, filename)
        grid_space = np.linspace(start=0, stop=box_size, num=cells_per_row+1,
                                 endpoint=True).tolist()  # num+1 because endpoint counts as well...

        self.grid_data = [[float(x), float(y), float(z)] for x in grid_space for y in grid_space for z in grid_space]

        # boundary conditions
        self.boundary_data = "PeriodicBox"

    def __generate_particles_FCC_SOPHIA(self, n, Phi, box_size, radius, simulation_id):
        # particle_data => Ellipse# => {position: [x,y], velocity=[vx, vy], ...}


        # print("box size",box_size,"volume fraction",((n*4./3.)*np.pi*major**3)/(box_size**3))

        np.random.seed(simulation_id+1)
        self.particle_data = []

        # check if n is a suitable particle number for fcc lattice
        check = round((n / 4) ** (1. / 3), 3)

        if not check ** 3 * 4 == n:
            print('no valid particle number! \n')
            sys.exit('no valid number')

        # fcc setup
        cells_per_row = check
        shift = box_size / (cells_per_row * 2)
        print(">>>",box_size,cells_per_row)

        base = np.linspace(start=0, stop=box_size, num=cells_per_row, endpoint=False)
        xx, yy, zz = np.meshgrid(base, base, base)

        fcc_1 = xx + shift
        fcc_2 = yy + shift
        fcc_3 = zz + shift

        # damit das output format fürs file passt:
        fcc_x = np.concatenate([xx, fcc_1, fcc_1, xx], axis=0).flatten()
        fcc_y = np.concatenate([yy, fcc_2, yy, fcc_2], axis=0).flatten()
        fcc_z = np.concatenate([zz, zz, fcc_3, fcc_3], axis=0).flatten()

        a = fcc_x.reshape(1, -1).T
        b = fcc_y.reshape(1, -1).T
        c = fcc_z.reshape(1, -1).T

        #quickens melting process later on :)
        fcc_lattice = np.concatenate((a, b, c), axis=1)

        ## check radius for insurance #TODO kann eigentlich raus (Sophia)
        #r = (Phi * 3 / (4 * np.pi * 4)) ** (1. / 3) * (box_size / cells_per_row)
        #print('radius really is ', r)
        ## check density:
        r_max = np.sqrt(2) * (box_size / cells_per_row) / 4
        if radius > r_max:
            print('ATTENTION: choose lower packing fraction')

        #Write Ovitofile for an outside System Check
        with open('test.txt', 'w') as outfile:
            outfile.write('ITEM: TIMESTEP \n')
            outfile.write(str(1) + '\n')
            outfile.write('ITEM: NUMBER OF ATOMS pp pp pp \n')
            outfile.write(str(n) + '\n')
            outfile.write('ITEM: BOX BOUNDS \n ')
            outfile.write("{a}\t{b}\n".format(a=0, b=box_size) * 3)
            outfile.write('ITEM: ATOMS id type x y z \n')

            factor = r_max - radius
            #factor = 1

            print(factor)
            for i in range(len(fcc_x - 1)):
                #fcc_lattice[i] = fcc_lattice[i].astype(float) + np.random.rand(3) * (1./(400 * Phi) + 0.1)
                fcc_lattice[i] = fcc_lattice[i].astype('f') + np.random.uniform(0.0, float(factor), 3)
                outfile.write(
                    '{w:} \t {w:}\t {x: .4f} \t {y: .4f} \t {z: .4f} \n'.format(w=i, x=fcc_x[i], y=fcc_y[i],

                                                                              z=fcc_z[i]))

                #print(fcc_lattice[i])
                self.particle_data.append(
                    {"type": "SphereOrientable", "id": i, "position": np.float64(fcc_lattice[i]), "velocity": np.float64([0, 0, 0]),
                     "angle": np.float64([0, 0, 0, 0]),
                     "angvel": np.float64([0.5, 0.5, 0]), "major": np.float64(1), "minor": np.float64(1),
                     "minor2": np.float64(1), "pinned": False, 'radius': np.float64(radius)})
        outfile.close()


        print(
            '\n --> Sophias FCC lattice is working with Circles right now \n --> and the produced file can be read')

    def __generate_particles_from_file(self, n, radius,  filename):
        """
        Reads the last frame of an ovito file format - A
        """

        f = open(filename, 'r')

        counter = 0
        ids = []
        type = []
        posits = []
        angles = []
        self.particle_data = []

        for line in f:
            counter += 1
            if counter == 2:
                time = float(line.strip())
                print("time", time)
            if counter == 4:
                natoms = int(line.strip())
                if n != natoms:
                    sys.exit('no of atoms do not match')
            if counter == 6:
                L = line.strip()
                L = float(L.split()[1])
            if counter > 9:  # number of header lines
                line = line.strip()
                columns = line.split()
                ids.append(int(columns[0]))
                type.append(int(columns[1]))
                posits.append(np.array([float(columns[2]), float(columns[3]), float(columns[4])]))
                angles.append(
                    np.array([float(columns[5]), float(columns[6]), float(columns[7]), float(columns[8])]))
                major = float(columns[9])
                minor = float(columns[10])
                minor2 = float(columns[11])

                self.particle_data.append(
                    {"type": "SphereOrientable",
                     "id": counter - 10,
                     "position": posits[counter - 10],
                     "velocity": [0, 0, 0],
                     "angle": angles[counter - 10],
                     "angvel": [0.5, 0.5, 0],
                     "major": np.float64(1),
                     "minor": np.float64(1),
                     "minor2": np.float64(1),
                     "pinned": False, 'radius': radius})

        box_size = L
        cubeN = math.ceil(n ** (1. / 3.))
        #print("cubeM", cubeN, "box size", box_size, "volume fraction",
         #     ((n * 4. / 3.) * np.pi * major * minor * minor2) / (box_size ** 3))

        """
        Uncomment the following lines if you need to check the read positionss
        OvitoFile = open('./data/ellipsoids.dump', 'w')
        OvitoFile.write('ITEM: TIMESTEP \n')
        OvitoFile.write(str(time) + '\n')
        OvitoFile.write('ITEM: NUMBER OF ATOMS \n')
        OvitoFile.write(str(n) + '\n')
        OvitoFile.write('ITEM: BOX BOUNDS pp pp pp \n')
        OvitoFile.write(str(0) + '\t' + str(box_size) + '\n')
        OvitoFile.write(str(0) + '\t' + str(box_size) + '\n')
        OvitoFile.write(str(0) + '\t' + str(box_size) + '\n')
        OvitoFile.write('ITEM: ATOMS id type x y z c_orient[1] c_orient[2] c_orient[3] c_orient[4] c_shape[1] c_shape[2] c_shape[3] \n')
        """

        # Writing format for quarternion px,py,pz,s

        """
        for i in range(n):
            OvitoFile.write(str(int(ids[i])) + '\t' + str(int(ids[i])) + '\t' + str(format(posits[i][0], '.4f')) + '\t' + str(format(posits[i][1], '.4f')) + '\t' + str(format(posits[i][2], '.4f')) + '\t' + str(format(angles[i][0], '.4f')) + '\t' + str(format(angles[i][1], '.4f')) + '\t' + str(format(angles[i][2], '.4f')) + '\t' + str(format(angles[i][3], '.4f')) + '\t' + str(format(major, '.4f')) + '\t' + str(format(minor, '.4f')) + '\t' + str(format(minor, '.4f')) + '\n')
            self.particle_data.append( {"type": "Ellipse", "id": ids[i] ,"position": posits[i], "velocity": np.array([0.,0.,0.]), "angle": angles[i], "angvel": np.array([0.,0.,0.]), "major": np.float64(major), "minor": np.float64(minor), "minor2": np.float64(minor2), "pinned": False} )
        """

        # # generate the cellspace grid
        # grid_space = np.linspace(start=0, stop=box_size, num=cells_per_row + 1,
        #                          endpoint=True).tolist()  # num+1 because endpoint counts as well...
        # self.grid_data = [[float(x), float(y), float(z)] for x in grid_space for y in grid_space for z in grid_space]
        #
        # # boundary conditions
        # self.boundary_data = "PeriodicBox"


class TestSystem(SystemGenerator):

    def __init__(self,box_size,particles):
        SystemGenerator.__init__(self)
        self.box_size = box_size
        self.particles = particles
        pass

    def generate(self, **kwargs):

        grid_space = np.linspace(start=0, stop=self.box_size/2, num=2,
                                 endpoint=True).tolist()  # num+1 because endpoint counts as well...
        self.grid_data = [[float(x), float(y), float(z)] for x in grid_space for y in grid_space for z in grid_space]
        # boundary conditions
        self.boundary_data = "PeriodicBox"

        self.particle_data = []
        for particle in self.particles:
            self.particle_data.append(particle)
