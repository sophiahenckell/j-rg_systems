import os

#from ..bedsim import simulation
from bedsim.simulation import Simulation
from bedsim.sysgen.ellipse_grid import P3EllipseGridPin
from bedsim.btools.animate import Animate 
#from bedsim.btools.statprop import Statprop

def start():
    # 1. generate system
    gen = P3EllipseGridPin()
    #(n,k,phi) = (16, 2, 0.6) # 32^2=1024
    n = 5
    phi = 0.70
    gen.generate(n=n,phi=phi)
    config_filename = "bedsim/data/test/EllipseBenchmark-Penrose-n%s-phi%s.h5" % (n, phi)
    try:
        os.remove(config_filename) # cleanup
    except FileNotFoundError:
        pass
    #gen.particle_data[0]['velocity']=[3.,0.5]
    gen.save_to_file(config_filename)
    
    # 2. simulate system
    #sim = Simulation(config_filename=config_filename, output_filename=config_filename, system_lifetime=10, brownian_timestep=None, saving_timestep=0.02, verbose=True)
    sim = Simulation(config_filename=config_filename, output_filename=config_filename, system_lifetime=10, brownian_timestep=0.1,  saving_timestep=0.5, verbose=True) # 'production' lifetime approx 1000
    #sim = Simulation(config_filename=config_filename, output_filename=config_filename, system_lifetime=10, brownian_timestep=None, saving_timestep=0.02, verbose=True)
    sim.start()
    
    # 3. animate
    Animate(input_filename=config_filename, output_filename=config_filename+".mp4")

    # 4. determine statistic data
    #Statprop(input_filename=config_filename)

if __name__ == '__main__':
    start()