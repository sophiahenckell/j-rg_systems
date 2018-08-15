import os

from bedsim.simulation import Simulation
from bedsim.sysgen.ellipse_grid import EllipseBenchmark, EllipseBenchmarkRandom
from bedsim.btools.animate import Animate 
from bedsim.btools.statprop import Statprop

def start():
    # 1. generate system
    gen = EllipseBenchmark()
    #gen = EllipseBenchmarkRandom()
    (n,k,phi) = (25, 2, 0.7) # 32^2=1024
    gen.generate(n=n, k=k, phi=phi)
    config_filename = "bedsim/data/test/EllipseBenchmark-n%s-k%s-phi%s.h5" % (n,k,phi)
    
    try:
        os.remove(config_filename) # cleanup
    except FileNotFoundError:
       pass
    
    #gen.particle_data[0]['velocity']=[3.,0.5]
    #gen.particle_data[0]['pinned'] = True
    gen.save_to_file(config_filename)
    
    # 2. simulate system
    #sim = Simulation(config_filename=config_filename, output_filename=config_filename, system_lifetime=10, brownian_timestep=None, saving_timestep=0.02, verbose=True)
    #sim = Simulation(config_filename=config_filename, output_filename=config_filename, system_lifetime=10, brownian_timestep=0.1,  saving_timestep=0.02, verbose=True)
    #sim = Simulation(config_filename=config_filename, output_filename=config_filename, system_lifetime=1, brownian_timestep=2,  saving_timestep=0.02, verbose=True)
    sim = Simulation(config_filename=config_filename, output_filename=config_filename, system_lifetime=10, brownian_timestep=0.1, saving_timestep=0.02, verbose=True)
    #sim = Simulation(config_filename=config_filename, output_filename=config_filename, system_lifetime=0.2, brownian_timestep=0.1, saving_timestep=0.001, verbose=True)
    sim.start()
    
    # 3. animate
    Animate(input_filename=config_filename, output_filename=config_filename+".mp4")

    # 4. determine statistic data
    Statprop(input_filename=config_filename)


    '''
    config_filename2 = "bedsim/data/test/EllipseBenchmark-n%s-k%s-phi%s-2.h5" % (n,k,phi)
    try:
        os.remove(config_filename2) # cleanup
    except FileNotFoundError:
       pass
    sim2 = Simulation(config_filename=config_filename, output_filename=config_filename2, system_lifetime=20, brownian_timestep=0.1,  saving_timestep=0.02, verbose=True)
    sim2.start()
    Animate(input_filename=config_filename2, output_filename=config_filename2+".mp4")
    Statprop(input_filename=config_filename2)
    '''


if __name__ == '__main__':
    start()