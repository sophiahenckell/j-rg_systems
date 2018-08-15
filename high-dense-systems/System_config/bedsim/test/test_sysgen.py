import os

import unittest
#import numpy as np

from bedsim.sysgen.ellipse_grid import EllipseGrid, P3EllipseGrid, EllipseBenchmark



class arb(object):
    """
    Arbitrary class for testing purposes.
    """
    pass


def _cleanup_and_save(filename, gen_handle):
        try:
            os.remove(filename) # cleanup
        except FileNotFoundError:
            pass
        gen_handle.save_to_file(filename)



class TestEllipseGrid(unittest.TestCase):

    def setUp(self):
        self.gen = EllipseGrid()

    def test_2x2(self):
        self.gen.generate(2)
        _cleanup_and_save("bedsim/data/test/EllipseGrid-2x2.h5", self.gen)

    def test_5x5(self):
        self.gen.generate(5)
        _cleanup_and_save("bedsim/data/test/EllipseGrid-5x5.h5", self.gen)


class TestEllipseBenchmark(unittest.TestCase):
    """
    Generate sample systems for later benchmarking
    """
    def setUp(self):
        self.gen = EllipseBenchmark()

    def test_gen(self):
        self.gen.generate(n=49, k=2, phi=0.5)
        _cleanup_and_save("bedsim/data/test/EllipseBenchmark-n7x7-phi0.5-k2.h5", self.gen)



class TestP3Grid(unittest.TestCase):
    
    def setUp(self):
        self.gen = P3EllipseGrid()

    def test_2x2(self):
        self.gen.generate(3)
        _cleanup_and_save("bedsim/data/test/P3EllipseGrid-2x2.h5", self.gen)

    """
    def test_P3_standard(self):
        self.gen.generate_standard_P3(2)
    """
