import os
from math import sqrt, pi

#from ..bedsim import simulation
from bedsim.simulation import Simulation
from bedsim.sysgen.circle_grid import CircleBenchmark
from bedsim.btools.animate import Animate
from bedsim.btools.statprop import Statprop


def start():
    # 1. generate system
    gen = CircleBenchmark()
    (n,phi) = (16, 0.5) # phi=0.5 get radius approx 1
    gen.generate(n=n, phi=phi)
    config_filename = "/home/mahe/Dokumente/Eclipse/pyBDsim/pyBDsim/bedsim/data/test/CircleBenchmark-NewtonsCradle-n%s-phi%s.h5" % (n,phi)
    try:
        os.remove(config_filename) # cleanup
    except FileNotFoundError:
        pass
    #config_filename = "bedsim/data/test/CircleBenchmark-n%s-phi%s.h5" % (n,phi)
    gen.particle_data[0]['velocity']=[2.0, 0.0] ###
    distance = sqrt(2*pi) - 2
    #tol = 0.0001
    tol = 0.0
    gen.particle_data[1]['position'] += [distance-tol, 0.0] ###
    gen.particle_data[3]['position'] -= [distance-tol, 0.0] ###
    #gen.particle_data[3]['pinned'] = True ###
    #gen.particle_data['Circle_0']['velocity']=[4.,1.] ###
    gen.save_to_file(config_filename)
    
    # 2. simulate system
    #sim = Simulation(config_filename=config_filename, system_lifetime=10, brownian_timestep=0.1, saving_timestep=0.02, verbose=True)
    sim = Simulation(config_filename=config_filename, output_filename=config_filename, system_lifetime=10, brownian_timestep=None, saving_timestep=0.02, verbose=True)
    sim.start()
    
    # 3. animate
    Animate(input_filename=config_filename, output_filename=config_filename+".mp4")

    # 4. determine statistic data
    #Statprop(input_filename=config_filename)


if __name__ == '__main__':
    start()
