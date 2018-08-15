# Checked -A
"""
To do:
1. Add minor2 in the main input
2. see notes on how to correct overlaps
"""
import argparse
import os
import sys
import numpy as np

from bedsim.simulation import Simulation
from bedsim.sysgen.ellipse_grid import EllipseBenchmark, EllipseBenchmarkRandom, FCCLatticeGenerator
from bedsim.btools.statprop import Statprop

def start(path, brownian_timestep, saving_timestep, swelling_rate, system_lifetime, packing_fraction, particle_number, aspect_ratio, radius, simulation_id,verbose):

    np.random.seed(simulation_id) # set seeding here so that it is easier to see

    ## choose system type
    gen = EllipseBenchmarkRandom()
    #gen = FCCLatticeGenerator()

    ########################
    #    choose process    #
    ########################

    # # A) if system already exists and only simluation of an equilibrated system is wanted:
    # oldfile = '/data/scc/sophia/FreeBROWNIANParticle1-n1372-t1000-phi0.05-k0-tb0.01-ts0.001/Tails/Tail{}.dump'.format(simulation_id)
    # #oldfile = '/home/newton/sophia/Desktop/curr_MASTER/Project/data/Tail.dump'
    #
    # gen.generate(n=particle_number, k=aspect_ratio, phi=packing_fraction, radius = radius, filename=oldfile, simulation_id = simulation_id)
    #
    # config_filename  = "%s/Circle-BrownianNotRescaled-n%s-k%s-phi%s-id%s.h5" % (path, particle_number, aspect_ratio, packing_fraction, simulation_id)
    # try:
    #     os.remove(config_filename) # cleanup
    # except FileNotFoundError:
    #     pass
    # gen.save_to_file(config_filename)
    #
    # sim2 = Simulation(config_filename=config_filename, output_filename=config_filename, system_lifetime=system_lifetime, brownian_timestep=brownian_timestep,  saving_timestep=saving_timestep, swelling_rate = None, packing_fraction = packing_fraction, aspect_ratio = aspect_ratio, verbose=verbose)
    # sim2.start()
    # sys.exit()

    ## B) if system was just generated and has to be equilibrated
    oldfile = None
    gen.generate(n=particle_number, k=aspect_ratio, phi=packing_fraction, radius = radius, filename=oldfile, simulation_id = simulation_id)
    config_filename = "%s/Ellipse-Benchmark-n%s-k%s-phi%s-id%s-equilibrate.h5" % (path, particle_number, aspect_ratio, packing_fraction, simulation_id)
    try:
        os.remove(config_filename) # cleanup
    except FileNotFoundError:
        pass
    gen.save_to_file(config_filename)

    ## START NEWTONIAN EQUILIBRATION
    sim = Simulation(config_filename=config_filename, output_filename=config_filename, system_lifetime=system_lifetime, brownian_timestep=brownian_timestep,  saving_timestep=saving_timestep, swelling_rate = None, packing_fraction = packing_fraction, aspect_ratio = aspect_ratio, verbose=verbose)
    sim.start()
    sys.exit()

    # --> THEN:  START THE BROWNIAN DYNAMICS
    # --> THEN: Calculate statistics etc




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform a Brownian simulation of the Wötzel benchmark ellipse system.')
    simgroup = parser.add_argument_group('Simulation settings')
    simgroup.add_argument('--path', metavar='path', help='Path to save simulation output to.', required=True)
    simgroup.add_argument('--brownian-timestep', metavar="t_B", type=float, help='Brownian time step. If not set use Newtonian dynamics.', required=True)
    simgroup.add_argument('--saving-timestep', metavar="t_s", type=float, help='Timestep for saving the system summary.', required=True)
    simgroup.add_argument('--swelling-rate', metavar="t_s", type=float, help='Timestep for saving the system summary.', required=True)
    simgroup.add_argument('--system-lifetime', metavar="T", type=float, help='Simulation duration.', required=True)
    simgroup.add_argument('--packing-fraction', metavar="phi", type=float, help='Packing fraction of the system between 0 and pi/4.', required=True)
    simgroup.add_argument('--particle-number', metavar="n", type=int, help='Number of particles in the system. Must be a square number.', required=True)
    simgroup.add_argument('--aspect-ratio', metavar="k", type=float, help='Aspect ratio of the ellipses', required=True)
    simgroup.add_argument('--radius', metavar="r", type=float, help='Particle radius. Default is 1.', default=1)
    simgroup.add_argument('--simulation_id', metavar="id", type=int, help='Simulation number, useful for array jobs. Default is 0.', default=0)
    parser.add_argument('--debug', action='store_true', help='print debug messages to stderr')
    parser.add_argument('--verbose', help='Print simulation status to stdout.', action='store_true')
    args = parser.parse_args()

    start(args.path, args.brownian_timestep, args.saving_timestep, args.swelling_rate, args.system_lifetime, args.packing_fraction, args.particle_number, args.aspect_ratio, np.float64(args.radius), args.simulation_id,  args.verbose)

    #start("./data", 0.01, 0.001, 0, 1000, 0.03, 32, 0, np.float64(1), 1, True)

