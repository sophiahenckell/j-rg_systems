import math

from bedsim.files import BedfileWriter

import numpy as np

class SystemGenerator(object):
    def generate(self, n, **kwargs):
        raise NotImplementedError()

    def save_to_file(self, filename): # FIXME: move this to files.py with metaformat
        f = BedfileWriter(filename, 'LongFormat', None)
        f.system.system_properties = {'localtime': 0}
        f.system.boundary = np.string_(self.boundary_data)
        f.system.grid = self.grid_data
        f.particles.save_statics(self.particle_data)


    def __init__(self):
        self.particle_data = None
        self.grid_data = []
        self.boundary_data = None



class CircleBenchmark(SystemGenerator):

    def generate(self, n, **kwargs):
        """
        Generate the system with n particles and 'PeriodicBox' boundary conditions.
        @param n Generate n particles
        @param phi packing fraction of the circles
        """
        phi = kwargs['phi']
        box_size = math.sqrt(math.pi * n/phi)
        cells_per_row = math.floor(box_size / (n/8.0)) # 4
        cell_size = box_size / cells_per_row

        self.__generate_particles_verbose(n, phi, box_size)

        # generate the cellspace grid
        grid_space = np.linspace(start=0, stop=box_size, num=cells_per_row, endpoint=True).tolist()
        self.grid_data = [[float(x),float(y)] for x in grid_space for y in grid_space]

        # boundary conditions
        self.boundary_data = "PeriodicBox"


    def __generate_particles_verbose(self, n, phi, box_size):
        # particle_data => Circle# => {position: [x,y], velocity=[vx, vy], ...}
        radius = math.sqrt(phi * box_size**2 / (n * math.pi))
        sqrtN = math.sqrt(n)
        self.particle_data = []
        for i in range(n):
            #partname = "Circle_%s" % i
            x = np.fmod(((i + 0.5)/sqrtN), 1.0) * box_size
            y = (math.floor((i+0.5)/sqrtN) + 0.5) * box_size / sqrtN
            #self.particle_data[partname] = {"position": np.array([x,y]), "velocity": np.array([0,0]), "radius": radius}
            #self.particle_data[partname] = {"type": "Circle", "id": i, "position": np.array([x,y], dtype=np.float64), "velocity": np.array([0.0,0.0], dtype=np.float64), "radius": radius, "major": radius, "minor": radius, "angle": 0.0} # FIXME
            self.particle_data.append( {"type": "Circle", "id": i, "position": np.array([x,y], dtype=np.float64), "velocity": np.array([0.0,0.0], dtype=np.float64), "radius": radius, "pinned": False} )