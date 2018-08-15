#!/usr/bin/env python3
"""
Created on Fri Feb 23 15:09:39 2018

@author: sophia
"""

from bedsim.sysgen import ellipse_grid

def test_grid_generation():
    grid_gen = ellipse_grid.EllipseBenchmarkRandom()
    grid_gen.generate(n=500,k=1.,phi=0.2)
    print(grid_gen.particle_data)

    grid_gen = ellipse_grid.FCCLatticeGenerator()
    grid_gen.generate(n=500,k=1.,phi=0.2)
    print(grid_gen.particle_data)
    print(getattr(grid_gen,"particle_data"))

    pass

if __name__ == "__main__":
    test_grid_generation()
    pass