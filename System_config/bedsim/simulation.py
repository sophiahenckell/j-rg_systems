# Checked -A
"""
To do:
1. Add minor2 in the list of main variables 

"""
"""@package bedsim
Simulate brownian dynamics of particles with an event driven approach.

Created on 13.01.2015

@author: mahe
"""

import argparse

from bedsim.system import System


class Simulation(object):
    """
    Edited by Aiyin
    Simulation class. This is the main class of the simulation framework. 
    """

    def start(self):
        self.system.simulate()

    def __init__(self, config_filename, output_filename, system_lifetime, brownian_timestep, saving_timestep,
                 verbose, **kwargs):

        self.system = System()
        self.system.config_filename = config_filename
        self.system.output_filename = output_filename
        self.system.load_from_file()
        # need to print some ovito writing procedure here
        self.system.system_properties.lifetime = system_lifetime
        self.system.system_properties.brownian_timestep = brownian_timestep
        self.system.system_properties.summary_timestep = saving_timestep
        for k in ["swelling_rate", "packing_fraction", "aspect_ratio"]:
            if k in kwargs:
                setattr(self.system.system_properties, k, kwargs[k])
        # self.system.system_properties.swelling_rate = swelling_rate
        # self.system.system_properties.packing_fraction = packing_fraction
        # self.system.system_properties.aspect_ratio = aspect_ratio
        self.system.system_properties.verbose = verbose

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform a colloid simulation.')
    simgroup = parser.add_argument_group('Simulation settings')
    # parser.add_argument('--config-filename', metavar='cfg.h5', type=argparse.FileType('r'), help='Filename of the hdf5 system configuration file.', nargs=1, required=True)
    simgroup.add_argument('--config-filename', metavar='cfg.h5', help='Filename of the hdf5 system configuration file.',
                          required=True)
    simgroup.add_argument('--output-filename', metavar='out.h5',
                          help='Filename of the hdf5 simulation file. If omitted, simulation data is saved to the file specified by --config-filename')
    simgroup.add_argument('--brownian-timestep', metavar="t_B", type=float,
                          help='Brownian time step. If not set use Newtonian dynamics.')
    simgroup.add_argument('--saving-timestep', metavar="t_s", type=float,
                          help='Timestep for saving the system summary.')
    simgroup.add_argument('--swelling-rate', metavar="gammat", type=float, help='Rate for particle swelling.')
    simgroup.add_argument('--system-lifetime', metavar="T", type=float, help='Simulation duration.')
    simgroup.add_argument('--packing-fraction', metavar="phi", type=float, help='Packing fraction of the system.')
    simgroup.add_argument('--aspect-ratio', metavar="k", type=float, help='Ratio of major to minor axis')
    parser.add_argument('--verbose', help='Print simulation status to stdout.', action='store_true')
    args = parser.parse_args()
    # print(args.config_filename)
    # print(args)

    sim = Simulation(args.config_filename, args.system_lifetime, args.brownian_timestep, args.saving_timestep,
                     args.swelling_rate, args.packing_fraction, args.aspect_ratio, args.verbose)
    sim.start()
