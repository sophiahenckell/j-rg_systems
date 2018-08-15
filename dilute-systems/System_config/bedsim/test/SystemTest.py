#!/usr/bin/python

import h5py
import System




handle = h5py.File("ellipsesconfig.h5", "r")
ms = System.EllipseSystem()
ms.load_from_file(handle)
print(ms.particles)


def basicHdf5tests():
    '''
    this code is just testing stuff... fix later
    '''
    partnum = 10
    particles = [Ellipse(position=[x,0], velocity=[0,0], angle=2, angvel=0, major=1, minor=1) for x in range(partnum)]
    data = np.array([particle.position+[particle.angle] for particle in particles]) # FIXME: only works if position is already a list!
    
    handle = h5py.File("test.hdf5", "w")
    handle.create_group('particles')
    handle.create_dataset('/particles/ParticleData', data=data)




'''
Some inheritance tests
'''
def test():
    A = bar()
    A.a()

class foo(object):
    def a(self):
        print(type(self).__name__)

class bar(foo):
    def a(self):
        foo.a(self)
