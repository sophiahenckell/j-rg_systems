# Checked -A
"""
Created on 13.01.2015

@author: mahe
"""

import sys
from math import pi, log10, floor

import numpy as np

import bedsim.boundary
import bedsim.particle
from bedsim.events import EventManager
from bedsim.cell import DelaunayCellspace
# from boundary import PeriodicBox

from bedsim.files import BedfileReader, BedfileWriter


# from bedsim.plot.sysplot import Plot


class SystemProperties(object):  # FIXME: das geht schoener ;)
    """
    Store properties of a simulated system in a unified way
    """

    @property
    def summary(self):
        res = {}
        """
        res['localtime']         = (self.localtime_round(), 'int64')
        res['brownian_timestep'] = (self.brownian_timestep, 'float64')
        res['summary_timestep']  = (self.summary_timestep,  'float64')
        """
        res['localtime'] = self.localtime_round()
        res['brownian_timestep'] = self.brownian_timestep
        res['summary_timestep'] = self.summary_timestep
        res['swelling_rate'] = self.swelling_rate
        res['lifetime'] = self.lifetime
        res['summary_timestep_number'] = self.summary_timestep_number
        res['brownian_timestep_number'] = self.brownian_timestep_number
        res['swelling_rate_number'] = self.swelling_rate_number
        return res

    @summary.setter
    def summary(self, bmf_data):
        # self.brownian_timestep = bmf_data.pop('brownian_timestep')
        ##self.localtime = bmf_data.pop('localtime') # FIXME: use other timing system
        # self.summary_timestep = bmf_data.pop('summary_timestep')
        pass

    # FIXME: maybe use decorator
    def localtime(self):  # FIXME: generalize this # FIXME: set this to the proper value from the h5 config file!
        localtime = 0.0
        try:  # FIXME: this is crap...
            localtime = self.system.event_manager.env.now
        except NameError:
            localtime = 0.0
        return localtime

    def localtime_round(self):
        """
        Rounds the localtime result to a proper value determined by the magnitude of
        the system summary time ts. Use this method for output values only not for
        simulation calculations. 
        """
        ts = self.summary_timestep
        if ts <= 1 and ts != 0:
            decimals = -floor(log10(ts))
        else:
            decimals = 0
        return round(self.localtime(), int(decimals))  ### NOTE: int() cast needed for python2 (maybe remove later)

    def __init__(self, system):
        self.system = system

        """ Parameters provided by argparse """
        self.lifetime = 0.0  # lifetime of the system, i.e. simulation duration
        self.brownian_timestep = 0.0  # brownian timestep
        self.summary_timestep = 0.0  # system summary timestep, i.e. system saving time interval
        self.swelling_rate = 0.0
        self.verbose = False  # Print simulation progress info during runtime

        self.summary_timestep_number = int(0)
        self.brownian_timestep_number = int(0)
        self.swelling_rate_number = int(0)


class System(object):
    """
    System base class.
    """

    def simulate(self):
        """
        Start the system dynamics after initialization.
        Activate event prediction + Event processing
        """
        self.event_manager.start()  # FIXME: generalize this #####
        if self.cnit != 0:
            print("mean it: ", self.nit / self.cnit)  ### PROFILING INFO!!

    # FIXME: obsolete... will be calculated in btools package
    # NOTE: this function here was for testing purposes only during the 'Antrag series'
    def get_packing_fraction(self):  # FIXME: maybe use decorator instead
        """
        FIXME: generalize this!
        use radius = systemsize/4
        NOTE: consider system center
        I still dont know what this is for. Check the equations :P -A
        """
        center = np.array([1, 1, 1])  # FIXME! => calculate system center from grid!
        radius = 0.8  # FIXME: use systemsize/4 or systemsize/2 * 0.8 or similar
        total_volume = 4. * pi * radius ** 3 / 3.
        covered_volume = 0
        for particle in self._particles:
            if np.linalg.norm(particle.position - center) <= radius:
                covered_volume += particle.volume
        return covered_volume / total_volume

    def compute_volume_fraction(self):
        [x, y, z] = self.cellspace._grid_points.transpose()
        [xmin, xmax, ymin, ymax, zmin, zmax] = [np.amin(x), np.amax(x), np.amin(y), np.amax(y), np.amin(z), np.amax(z)]
        [width, height, depth] = [xmax - xmin, ymax - ymin, zmax - zmin]
        total_volume = width * height * depth
        covered_volume = 0
        for particle in self._particles:
            covered_volume += particle.volume()  # particle.major*particle.minor*particle.minor2 #
        return covered_volume / total_volume  # 4.*np.pi*covered_volume/(total_volume*3)

    def box_volume(self):
        [x, y, z] = self.cellspace._grid_points.transpose()
        [xmin, xmax, ymin, ymax, zmin, zmax] = [np.amin(x), np.amax(x), np.amin(y), np.amax(y), np.amin(z), np.amax(z)]
        [width, height, depth] = [xmax - xmin, ymax - ymin, zmax - zmin]
        print("inside box volume")
        return width * height * depth

    def __get_simbox_packing_fraction(self):
        [x, y, z] = self.cellspace._grid_points.transpose()
        [xmin, xmax, ymin, ymax, zmin, zmax] = [np.amin(x), np.amax(x), np.amin(y), np.amax(y), np.amin(z), np.amax(z)]
        [width, height, depth] = [xmax - xmin, ymax - ymin, zmax - zmin]
        total_volume = width * height * depth
        covered_volume = 0
        for particle in self._particles:
            covered_volume += particle.volume
        return covered_volume / total_volume

    """
    IO related methods
    ==================
    """

    def load_from_file(self):
        """ 0. Load file """
        bf = BedfileReader(self.config_filename)

        """ 1. Load system properties """
        self.system_properties.summary = bf.system.system_properties

        """ 2. Load Boundary and cell grid """
        self.cellspace = DelaunayCellspace(system=self, grid_points=bf.system.grid)
        print('hey?')
        self.boundary = getattr(bedsim.boundary, bf.system.boundary)(system=self)

        """ 3. Load Particles """
        particle_list = next(bf.particles.load_statics())
        self._particles = [getattr(bedsim.particle, p['type'])(**p) for p in particle_list]

        """ 4. Update Particles to latest configuration if available """
        particles_last = bf.particles.load_last_dynamics()
        if particles_last is not None:
            for particle in self._particles:
                for pl in particles_last:
                    if int(pl['id']) == particle._id:
                        particle.position = pl['position']
                        particle.angle = pl['angle']
                        particle.major = pl['major']
                        particle.minor = pl['minor']
                        particle.minor2 = pl['minor2']
                        particle.cumulative_position = pl['cumulative_position'] # LOAD THESE ONLY IF YOU WANT TO RESUME A SIMULATION!! ELSE IT WILL MESS UP THE MSD!
                        # particle.cumulative_angle = pl['cumulative_angle']
                        # particle.velocity = pl['velocity']
                        # particle.angvel = pl['angvel']

        """ 5. If output filename is not input filename then save initial configuration as well """
        if self.config_filename != self.output_filename:
            self.save_to_file_initial()

        """ 6. Assign particle to a cell """

        [self.cellspace.assign_particle(particle) for particle in self._particles]

        # ENABLE IF YOU WANT TO PLOT SYSTEM DURING SIMULATION
        # self.__plotter.system_sync() # synchronize plotter with system # FIXME: geht schoener... hoffentlich :P

    def save_to_file_initial(self):  # FIXME: change in all clients needed...
        summary = [particle.get_summary() for particle in
                   self._particles]  # save complete initial state to the statics table
        # num_timesteps = self.system_properties.lifetime / self.system_properties.summary_timestep + 1

        # f = BedfileWriter(filename=self.output_filename, particle_format='ShortFormat', num_timesteps=num_timesteps)
        f = BedfileWriter(filename=self.output_filename, particle_format='ShortFormat',
                          num_timesteps=None)  # only use static methods
        f.particles.save_statics(data=summary)
        f.system.system_properties = self.system_properties.summary
        f.system.boundary = self.boundary.__class__.__name__
        f.system.grid = self.cellspace._grid_points

    # FIXME: save statics even if simulation was not created by a h5 file!
    def save_to_file(self):
        """
        Saves the current system state to the file provided by 'handle'.
        @param handle Handle to an hdf file.
        """
        # FIXME: save writer handle on first call
        timevar = self.system_properties.localtime_round()
        timestep = round(timevar / self.system_properties.summary_timestep)

        """ Create summary and save (NEW METHOD) """
        summary = [particle.get_summary() for particle in
                   self._particles]  # FIXME: only use dynamics here when changing to full API
        num_timesteps = self.system_properties.lifetime / self.system_properties.summary_timestep + 1

        f = BedfileWriter(filename=self.output_filename, particle_format='ShortFormat', num_timesteps=num_timesteps)
        #print("where is the file", self.output_filename) #SOPHIA
        # f.particles.save_dynamics(data=summary, time=timevar)
        f.particles.save_dynamics(data=summary, time=timestep)
        f.system.system_properties = self.system_properties.summary

        # self.__plotter.plot_system("foo_%s.png" % timestep)
        if self.system_properties.verbose:  ## FIXME: maybe remove from save_file and update only after e.g. 2seconds elsewhere
            self.progress_info()
            self.write_trajectories()





        #     def __mean_squared_velocity(self):
        #         # check mean squared velocity -> save this to file and don't print to cli
        #         # could be added to saving routine or wherever you want
        #         vs = 0
        #         for particle in self._particles:
        #             vs += np.dot(particle.velocity,particle.velocity)
        #         print("Deviation ", vs/len(self._particles))

    def progress_info(self):
        """
        Prints progress status of the current simulation to stdout.
        """
        progress = round(self.system_properties.localtime_round() / self.system_properties.lifetime * 100, 2)
        sys.stdout.write("\rSimulation status: %.2f%%" % progress)
        sys.stdout.flush()

    def write_trajectories(self):

        n = len(self._particles)

        [x, y, z] = self.cellspace._grid_points.transpose()
        [xmin, xmax, ymin, ymax, zmin, zmax] = [np.amin(x), np.amax(x), np.amin(y), np.amax(y), np.amin(z), np.amax(z)]
        [width, height, depth] = [xmax - xmin, ymax - ymin, zmax - zmin]

        now = self.system_properties.localtime()  # or put localtime_round() instead, but for checking this does not help

        if now == 0:
            writemode = "w"
        else:
            writemode = "a+"

        OvitoFile = open(self.output_filename + '.dump', writemode)
        OvitoFile.write('ITEM: TIMESTEP \n')
        OvitoFile.write(str(now) + '\n')
        OvitoFile.write('ITEM: NUMBER OF ATOMS \n')
        OvitoFile.write(str(n) + '\n')
        OvitoFile.write('ITEM: BOX BOUNDS pp pp pp \n')
        OvitoFile.write(str(0) + '\t' + str(width) + '\n')
        OvitoFile.write(str(0) + '\t' + str(height) + '\n')
        OvitoFile.write(str(0) + '\t' + str(depth) + '\n')
        OvitoFile.write(
            'ITEM: ATOMS id type x y z c_orient[1] c_orient[2] c_orient[3] c_orient[4] c_shape[1] c_shape[2] c_shape[3] \n')
        """
        Writing format for quarternion px,py,pz,s
        """
        print("saving file *** time", now)
        for particle in self._particles:
            OvitoFile.write(str(int(particle._id)) + '\t' + str(int(particle._id)) + '\t' + str(
                format(particle.position[0], '.4f')) + '\t' + str(format(particle.position[1], '.4f')) + '\t' + str(
                format(particle.position[2], '.4f')) + '\t' + str(format(particle.angle[0], '.4f')) + '\t' + str(
                format(particle.angle[1], '.4f')) + '\t' + str(format(particle.angle[2], '.4f')) + '\t' + str(
                format(particle.angle[3], '.4f')) + '\t' + str(format(particle.major, '.4f')) + '\t' + str(
                format(particle.minor, '.4f')) + '\t' + str(format(particle.minor2, '.4f')) + '\n')
        """
        Writing to extra file for MSD calculations SOPHIA
        """
        OvitoFile1 = open(self.output_filename + '.comulative.dump', writemode)
        OvitoFile1.write('ITEM: TIMESTEP \n')
        OvitoFile1.write(str(now) + '\n')
        OvitoFile1.write('ITEM: ATOMS id type x y z \n')
        for particle in self._particles:
            OvitoFile1.write(str(int(particle._id)) + '\t' + str(format(particle.cumulative_position[0], '.4f')) + '\t' + str(format(particle.cumulative_position[1], '.4f')) + '\t' + str(
                format(particle.cumulative_position[2], '.4f')) + '\n')

    """
    Constructor
    ===========
    """

    def __init__(self, particles=[], system_properties=None):
        self.nit = 0
        self.cnit = 0

        self.config_filename = None
        self.output_filename = None
        self._particles = particles  # save for quick reference, instead of traversing self.cellspace.cells[].particles

        if system_properties is not None:
            self.system_properties = system_properties
        else:
            self.system_properties = SystemProperties(self)
        self.event_manager = EventManager(self)  # FIXME: generalize this! (Goal 2)

        self.cellspace = None
        self.boundary = None
        # self.__plotter = Plot(self)
        # FIXME: correct variance?

        # np.random.seed(0) ## FIXED SEED FOR DEBUGGING
