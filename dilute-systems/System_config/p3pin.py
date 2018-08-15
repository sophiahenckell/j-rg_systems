import argparse
import os

from bedsim.simulation import Simulation
from bedsim.sysgen.ellipse_grid import P3EllipseGridPin
#from bedsim.btools.animate import Animate
#from bedsim.btools.statprop import Statprop

def start(path, brownian_timestep, saving_timestep, system_lifetime, packing_fraction, p3_iterations, simulation_id, verbose):
    # 1. generate system
    gen = P3EllipseGridPin()
    gen.generate(n=p3_iterations,phi=packing_fraction)
    config_filename = "%s/Ellipse-Penrose-n%s-phi%s-id%s.h5" % (path, p3_iterations, packing_fraction, simulation_id)
    try:
        os.remove(config_filename) # cleanup
    except FileNotFoundError:
        pass
    #gen.particle_data[0]['velocity']=[3.,0.5]
    gen.save_to_file(config_filename)
    
    # 2. simulate system
    sim = Simulation(config_filename=config_filename, output_filename=config_filename, system_lifetime=system_lifetime, brownian_timestep=brownian_timestep,  saving_timestep=saving_timestep, verbose=verbose) # 'production' lifetime approx 1000
    sim.start()
    
    # 3. animate
    #Animate(input_filename=config_filename, output_filename=config_filename+".mp4")

    # 4. determine statistic data
    #Statprop(input_filename=config_filename)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform a Brownian simulation of a Penrose P3 tiling with pinned boundary particles.')
    simgroup = parser.add_argument_group('Simulation settings')
    simgroup.add_argument('--path', metavar='path', help='Path to save simulation output to.', required=True)
    simgroup.add_argument('--brownian-timestep', metavar="t_B", type=float, help='Brownian time step. If not set use Newtonian dynamics.', required=True)
    simgroup.add_argument('--saving-timestep', metavar="t_s", type=float, help='Timestep for saving the system summary.', required=True)
    simgroup.add_argument('--system-lifetime', metavar="T", type=float, help='Simulation duration.', required=True)
    simgroup.add_argument('--packing-fraction', metavar="phi", type=float, help='Packing fraction of the system between 0 and pi/4.', required=True)
    simgroup.add_argument('--p3-iterations', metavar="n", type=int, help='Number of P3 subdivisions to perform.', required=True)
    simgroup.add_argument('--simulation-id', metavar="id", type=int, help='Simulation number, useful for array jobs. Default is 0.', default=0)
    parser.add_argument('--verbose', help='Print simulation status to stdout.', action='store_true')
    args = parser.parse_args()

    start(args.path, args.brownian_timestep, args.saving_timestep, args.system_lifetime, args.packing_fraction, args.p3_iterations, args.simulation_id, args.verbose)
