import unittest
import numpy as np

import bedsim
from bedsim.system import System
from bedsim.particle import Ellipse
from bedsim.cell import DelaunayCellspace
from bedsim.plot.sysplot import Plot
from numpy import dtype


def numeq(x, y):
    """
    Function to determine approx equality of two numeric variables
    """
    epsilon = 10**-3
    return np.abs(x-y) < epsilon



class TestCollision(unittest.TestCase):
    def setUp(self):
        #grid_space = np.linspace(start=-3, stop=3, num=2, endpoint=True).tolist() # num+1 because endpoint counts as well...
        grid_space = np.linspace(start=-15, stop=15, num=2, endpoint=True).tolist() # num+1 because endpoint counts as well...
        grid_data = np.array([[x,y] for x in grid_space for y in grid_space])

        self.system = System()
        self.system.cellspace = DelaunayCellspace(system=self.system, grid_points=grid_data)
        self.system.boundary  = getattr(bedsim.boundary, 'PeriodicBox')(system=self.system)

        e1 = Ellipse(position=[0,0], velocity=[0,0], angle=0, angvel=0, major = 1, minor = 0.5)
        e2 = Ellipse(position=[1,1], velocity=[0,0], angle=0, angvel=0, major = 1, minor = 0.5)
        
        #e1 = Ellipse(position=[0,0], velocity=[0,0], angle=np.pi/4, angvel=0, major = 1, minor = 0.5)
        #e2 = Ellipse(position=[1.6,0.3], velocity=[0.0000001,0], angle=np.pi/5, angvel=-1*np.pi, major = 1, minor = 0.5)
        e1._id = 1
        e2._id = 2
        
        particles = [e1,e2]
        self.system._particles = particles ## check if this works
        [self.system.cellspace.assign_particle(particle) for particle in self.system._particles]
    

    def test_collision_case1(self):
        [e1, e2] = self.system._particles
        
        e1.angle = np.pi/4
        
        e2.position = np.array([1.6,0.3])
        e2.velocity = np.array([0.0000001,0.0])
        e2.angle = np.pi/5
        e2.angvel = -np.pi
        
        ct = e1.collision_time(e2)
        #print("ct=",ct)
        self.assertTrue( numeq(ct, np.array([0.05743945, 0.71263422, 0.15376019])).all() )
        ### 
        #print("\n\n>>>\nStandard method: ", ct)
        #print("tqw method: ", e1.tpw_collision_time(e2, T=2))
        tpw_ct = e1.tpw_collision_time3(e2, T=2)
        self.assertTrue( numeq(tpw_ct, np.array([0.05743945, 0.71263422, 0.15376019])).all() )
        ###

        p = Plot(self.system)
        p.system_sync()
        p.plot_system('bedsim/data/test/EllipseCollisionTest.eps')


    #@unittest.skip("Can be fixed with correction move.")
    def test_collision_case2(self):
        [e1, e2] = self.system._particles
        
        """
        # spätester Zeitpunkt
        #e1: angle = 6.2046, pos=[3.848,5.971], vel=[0.196, -0.1] ## major=1.41 => scale here!
        #e2: angle = 5.6671, pos=[3.188,7.788], vel=[0.067,0.069]        
        e1.velocity = np.array([0.196,-0.1])
        e1.angle = 6.2046
        e1.angvel = 0.01 * np.pi # ?

        e2.position = np.array([-0.66,1.871])/1.41
        e2.velocity = np.array([0.067, 0.069])
        e2.angle = 5.6671
        e2.angvel = -np.pi # ?
        
        
        # earlier point of time
        # e1: angle = 0.2667, pos=[3.141, 7.740], vel=[0.06735, 0.06891]
        # e2: angle = 6.1283, pos=[3.710, 6.041], vel=[0.19633,-0.09960]
        e1.velocity = np.array([0.06735, 0.06891])
        e1.angle = 0.2667
        e1.angvel = 0.01 * np.pi # ?

        e2.position = np.array([-0.569, 1.699])/1.41
        e2.velocity = np.array([0.19633,-0.09960])
        e2.angle = 6.1283
        e2.angvel = -np.pi # ?
        """
        
        # even earlier
        # e1: angle = 6.253, pos=[3.582, 5.990], vel=[0.4876, 0.5994] 
        # e2: angle = 0.424, pos=[3.070, 7.613], vel=[0.3390,-0.0162]
        e1.velocity = np.array([0.4876, 0.5994])
        e1.angle = 6.253
        e1.angvel = 0.01 * np.pi # ?

        e2.position = np.array([-0.512, 1.623])/1.41
        e2.velocity = np.array([0.3390,-0.0162])
        e2.angle = 0.424
        e2.angvel = -np.pi # ?
        
        
        self.system.particles = [e1,e2]
        
        #ct = e1.collision_time(e2)
        #print("ct=",ct)
        #self.assertTrue( numeq(ct, np.array([ 0.19892904, 0.01777163, 0.61919515 ])).all() )  # not working with SLSQP
        
        tpw_ct = e1.tpw_collision_time3(e2, T=2)
        #print(">> tpw=",tpw_ct)
        self.assertTrue( numeq(tpw_ct, np.array([ 0.19892904, 0.01777163, 0.61919515 ])).all() )

        #ct = [0.3]
        #self.system.system_properties.localtime = lambda: ct[0]
        #[p.update() for p in self.system._particles]

        p = Plot(self.system)
        p.system_sync()
        p.plot_system('bedsim/data/test/EllipseCollisionTest2.eps')


    #@unittest.skip("Can be fixed with correction move.")
    def test_collision_case3(self):
        [e1, e2] = self.system._particles
        
        # even earlier
        # e1: angle = 2.007, pos=[6.464,1.847], vel=[-0.0512,0.64693] 
        # e2: angle = 0.146, pos=[8.834,3.417], vel=[0.2624,-0.7819]
        e1.velocity = np.array([-0.0512,0.64693])
        e1.angle = 2.007
        e1.angvel = 1 * np.pi # ?

        e2.position = np.array([2.37, 1.57])/1.41
        e2.velocity = np.array([0.2624,-0.7819])
        e2.angle = 0.146
        e2.angvel = 0.01 * np.pi # ?
        
        self.system.particles = [e1,e2]
        
        #ct = e1.collision_time(e2)
        #print("ct=",ct)
        #self.assertTrue( numeq(ct, np.array([0.37866501, 0.87956003, 0.50029793])).all() ) # not working with SLSQP 

        tpw_ct = e1.tpw_collision_time3(e2, T=2)
        #print(">> tpw=",tpw_ct)
        self.assertTrue( numeq(tpw_ct, np.array([0.37866501, 0.87956003, 0.50029793])).all() )

        #ct = [0.4]
        #self.system.system_properties.localtime = lambda: ct[0]
        #[p.update() for p in self.system._particles]

        p = Plot(self.system)
        p.system_sync()
        p.plot_system('bedsim/data/test/EllipseCollisionTest3.eps')
    
    
    #@unittest.skip("Can be fixed with correction move.")
    def test_collision_case4(self):
        [e1, e2] = self.system._particles
        
        # even earlier
        # e1: angle=2.429, pos=[6.286,1.211], vel=[-1.111,0.818] 
        # e2: angle=0.111, pos=[6.004,3.398], vel=[0.584,-0.869]
        e1.velocity = np.array([-1.111,0.818])
        e1.angle = 2.429
        e1.angvel = 1 * np.pi # ?

        e2.position = np.array([-0.282, 2.187])/1.41
        e2.velocity = np.array([0.584,-0.869])
        e2.angle = 0.111
        e2.angvel = 0.1 * np.pi # ?
        
        self.system.particles = [e1,e2]
        
        #ct = e1.collision_time(e2)
        #print("ct=",ct) 
        #self.assertTrue( numeq(ct, np.array([0.29561521, -0.16980752, 0.76979694])).all() )  # not working with SLSQP

        tpw_ct = e1.tpw_collision_time3(e2, T=2)
        #print(">> tpw=",tpw_ct)
        self.assertTrue( numeq(tpw_ct, np.array([0.29561521, -0.16980752, 0.76979694])).all() )  # not working with SLSQP

        #ct = [0.4]
        #self.system.system_properties.localtime = lambda: ct[0]
        #[p.update() for p in self.system._particles]

        p = Plot(self.system)
        p.system_sync()
        p.plot_system('bedsim/data/test/EllipseCollisionTest4.eps')
        
    
    #@unittest.skip("Can be fixed with correction move.")
    def test_collision_case5(self): # particles touch at t=1.687 (slight overlap for delta t_c approx 0.02
        #print("\n>>> C5")
        [e1, e2] = self.system._particles

        #e1: pos= [ 5.3607846   1.03103814] ; vel= [ 0.32634922  0.1371882 ] ; angle= 0.724034954957 ; angvel= 0.641544379574 ; a= 1.41421356237 ; b= 0.707106781187
        #e2: pos= [ 5.50824035 -1.01659235] ; vel= [ 0.36867608 -0.23015824] ; angle= 0.602571218602 ; angvel= -1.09635477791 ; a= 1.41421356237 ; b= 0.707106781187

        e1.position = np.array([5.3607846, 1.03103814])
        e1.velocity = np.array([0.32634922, 0.1371882])
        e1.angle = 0.724034954957
        e1.angvel = 0.641544379574
        e1.axes = (1.41421356237, 0.707106781187)

        #e2.position = np.array([5.50824035, -1.01659235]) - np.array([5.3607846, 1.03103814])
        e2.position = np.array([5.50824035, -1.01659235])
        e2.velocity = np.array([0.36867608, -0.23015824])
        e2.angle = 0.602571218602
        e2.angvel = -1.09635477791
        e2.axes = (1.41421356237, 0.707106781187)
        
        self.system.particles = [e1,e2]
        
        #ct = e1.collision_time(e2) 
        #self.assertTrue( numeq(ct, np.array([0.29561521, -0.16980752, 0.76979694])).all() )  # not working with SLSQP
        #self.assertIsNone(ct)

        tpw_ct = e1.tpw_collision_time3(e2, T=2)
        #print(">> tpw=",tpw_ct)
        self.assertIsNone(tpw_ct)
        ##self.assertTrue( numeq(tpw_ct, np.array([1.68801957,  5.98710966, -0.08586885])).all() )

        #ct = [0.4]
        #self.system.system_properties.localtime = lambda: ct[0]
        #[p.update() for p in self.system._particles]

        tc = 0
        e1.position += e1.velocity * tc
        e1.angle += e1.angvel * tc
        e2.position += e2.velocity * tc
        e2.angle += e2.angvel * tc

        p = Plot(self.system)
        p.system_sync()
        p.plot_system('bedsim/data/test/EllipseCollisionTest5.eps')
        
        
        
    #@unittest.skip("Can be fixed with correction move.")
    def test_collision_case6(self): # particles slightly overlap at t approx 0.3-0.4, but collide heavily at t approx 1.95
        #print("\n\n>>> C6")
        [e1, e2] = self.system._particles

        e1.position = np.array([ 9.52225511, -3.20048553])
        e1.velocity = np.array([ 0.64510243,  0.12668314])
        e1.angle = 0.786639981339
        e1.angvel = 0.867581306646
        e1.axes = (1.41421356237, 0.707106781187)

        e2.position = np.array([ 11.68788583,  -1.0466838 ])
        e2.velocity = np.array([-0.50144804, -0.04025363])
        e2.angle = 0.81061092758
        e2.angvel = -0.8609267932
        e2.axes = (1.41421356237, 0.707106781187)
        
        self.system.particles = [e1,e2]
        
        #ct = e1.collision_time(e2)
        #print("ct=",ct) 
        #self.assertTrue( numeq(ct, np.array([0.29561521, -0.16980752, 0.76979694])).all() )  # not working with SLSQP


        from bedsim.particle import ConvergenceError
        try:
            self.assertRaises(ConvergenceError, e1.tpw_collision_time3(e2, T=2))
        except ConvergenceError as e:
            ((ce,),) = e.args # FIXME: this needs apaption once the arg passing is fixed!
            self.assertTrue( ( ce-0.3587230485092775) < 1e-4 )
        else:
            self.assertTrue(False)


        #tpw_ct = e1.tpw_collision_time3(e2, T=2)
        #print(">> tpw=",tpw_ct)
        #self.assertTrue( numeq(tpw_ct, np.array([1.95323102,  10.84500997,  -2.15133253])).all() )  # not working with SLSQP
        #self.assertTrue( numeq(tpw_ct, np.array([1.94231489,  10.84123835,  -2.24599371])).all() )  # not working with SLSQP
        #self.assertTrue( numeq(tpw_ct, np.array([0.35872289,  10.49897711,  -1.79823235])).all() )  # not working with SLSQP
        
        #ct = [0.4]
        #self.system.system_properties.localtime = lambda: ct[0]
        #[p.update() for p in self.system._particles]

        # manual testing
        tc = 0
        e1.position += e1.velocity * tc
        e1.angle += e1.angvel * tc
        e2.position += e2.velocity * tc
        e2.angle += e2.angvel * tc
        
        p = Plot(self.system)
        p.system_sync()
        p.plot_system('bedsim/data/test/EllipseCollisionTest6.eps')
        
    
    #@unittest.skip("Can be fixed with correction move.")
    def test_collision_case7(self):
        [e1, e2] = self.system._particles

        e1.angle = 0.0
        e1.angvel = 0
        e1.axes = (1.41421356237, 1.0)

        e2.position = np.array([ 1.41421356237*5, 0 ])
        e2.velocity = np.array([-1.5, -0.04025363])
        e2.angle = 0
        e2.angvel = -0.8609267932
        e2.axes = (1.41421356237*3, 0.707106781187)
        
        self.system.particles = [e1,e2]
        
        #ct = e1.collision_time(e2)
        #print("ct=",ct) # yields a wrong result!
        #self.assertTrue( numeq(ct, np.array([0.29561521, -0.16980752, 0.76979694])).all() )  # not working with SLSQP

        tpw_ct  = e1.tpw_collision_time3(e2, T=2)
        tpw_ct2 = e1.tpw_collision_time3(e2, T=10)
        #print(">> tpw=",tpw_ct) # yields None in interval [0,2] which is correct
        #print(">> tpw2=",tpw_ct2)
        self.assertIsNone(tpw_ct)
        self.assertTrue( numeq(tpw_ct2, np.array([ 2.8615022,   1.12294173, -0.60786588])).all() )
        #self.assertTrue( numeq(tpw_ct, np.array([  1.95323102,  10.84500997,  -2.15133253])).all() )  # not working with SLSQP
        
        #ct = [0.4]
        #self.system.system_properties.localtime = lambda: ct[0]
        #[p.update() for p in self.system._particles]

        # manual testing
        #tc = tpw_ct2[0]
        tc = 0
        e1.position += e1.velocity * tc
        e1.angle += e1.angvel * tc
        e2.position += e2.velocity * tc
        e2.angle += e2.angvel * tc
        
        p = Plot(self.system)
        p.system_sync()
        p.plot_system('bedsim/data/test/EllipseCollisionTest7.eps')
        
        
        
    def test_collision_case8(self):
        [e1, e2] = self.system._particles

        
        e1.position = np.array([ 1.07613011,  7.53117917])
        e1.velocity = np.array([ 1.03287279, -0.07172733])
        e1.angle = 1.17372250848
        e1.angvel = -0.748346176117
        e1.axes = (1.41421356237, 0.707106781187)
         
        e2.position = np.array([ 3.10022877,  9.30370458])
        e2.velocity = np.array([-0.5131618,  -0.12945531])
        e2.angle = 0.609417513578
        e2.angvel = 0.147163595785
        e2.axes = (1.41421356237, 0.707106781187)
        
        self.system.particles = [e1,e2]
        #ct = e1.collision_time(e2)
        #print("ct=",ct) # yields a wrong result!
        #self.assertTrue( numeq(ct, np.array([0.29561521, -0.16980752, 0.76979694])).all() )  # not working with SLSQP

        tpw_ct = e1.tpw_collision_time3(e2, T=2)
        #print(">> tpw=",tpw_ct) # yields None in interval [0,2] which is correct
        #self.assertTrue( numeq(tpw_ct, np.array([  1.95323102,  10.84500997,  -2.15133253])).all() )  # not working with SLSQP
        self.assertTrue( numeq(tpw_ct, np.array([ 0.0198729,   1.89241964,  8.5635192 ])).all() )
        
        
        
        #ct = [0.4]
        #self.system.system_properties.localtime = lambda: ct[0]
        #[p.update() for p in self.system._particles]

        # manual testing
        tc = 0
        e1.position += e1.velocity * tc
        e1.angle += e1.angvel * tc
        e2.position += e2.velocity * tc
        e2.angle += e2.angvel * tc
        
        p = Plot(self.system)
        p.system_sync()
        p.plot_system('bedsim/data/test/EllipseCollisionTest8.eps')



    def test_collision_case9(self):
        [e1, e2] = self.system._particles

         
        e1.position = np.array([ 0.89375006, 5.1593676 ])
        e1.velocity = np.array([ 1.06831077, 1.18676223])
        e1.angle = 0.61559279607
        e1.angvel = 1.37448242157
        e1.axes = (1.41421356237, 0.707106781187)
         
        e2.position = np.array([ 3.36288631, 5.39311034])
        e2.velocity = np.array([ 1.04045937, 2.34821381])
        e2.angle = 0.828666163424
        e2.angvel = 0.132454699728
        e2.axes = (1.41421356237, 0.707106781187)
        
        self.system.particles = [e1,e2]

        tpw_ct = e1.tpw_collision_time3(e2, T=2)
        #print(">> tpw=",tpw_ct) # yields None in interval [0,2] which is correct
        self.assertIsNone(tpw_ct) # FIXME: this is correct
        #from bedsim.particle import ConvergenceError
        #self.assertRaises(ConvergenceError, e1.tpw_collision_time, e2, 2) # should be None... but better than nothing :P

        # manual testing
        tc = 0
        e1.position += e1.velocity * tc
        e1.angle += e1.angvel * tc
        e2.position += e2.velocity * tc
        e2.angle += e2.angvel * tc
        
        p = Plot(self.system)
        p.system_sync()
        p.plot_system('bedsim/data/test/EllipseCollisionTest9.eps')


    def test_collision_case10(self):
        [e1, e2] = self.system._particles


        e1.position = np.array([ 0.19133246, -0.31316855])
        e1.velocity = np.array([-2.05093999, -2.05897476])
        e1.angle = 0.872930657827
        e1.angvel = -1.59162222286
        e1.axes = (1.8, 1)
         
        e2.position = np.array([-1.09235624, -3.60720701])
        e2.velocity = np.array([ 1.76458492,  0.33911143])
        e2.angle = 6.00601742514
        e2.angvel =  0.747276535373
        e2.axes = (3, 0.5)
        
        self.system.particles = [e1,e2]

        tpw_ct = e1.tpw_collision_time3(e2, T=2)
        #print(">> tpw=",tpw_ct) # yields None in interval [0,2] which is correct
        self.assertTrue( numeq(tpw_ct, np.array([ 0.78107611, -0.21226556, -2.97946729])).all() )
        #self.assertIsNone(tpw_ct) # FIXME: this is correct
        #from bedsim.particle import ConvergenceError
        #self.assertRaises(ConvergenceError, e1.tpw_collision_time, e2, 2) # should be None... but better than nothing :P
        
        # manual testing
        tc = 0
        e1.position += e1.velocity * tc
        e1.angle += e1.angvel * tc
        e2.position += e2.velocity * tc
        e2.angle += e2.angvel * tc
        
        p = Plot(self.system)
        p.system_sync()
        p.plot_system('bedsim/data/test/EllipseCollisionTest10.eps')


         
    def test_collision_case11(self):
        [e1, e2] = self.system._particles
 
        e1.position = np.array([ 1.23602047,  5.47299557])
        e1.velocity = np.array([ 0.94279977,  0.56015277])
        e1.angle = 0.785398163397
        e1.angvel = -2.44546089624
        e1.axes = (1.41421356237, 0.707106781187)
         
        e2.position = np.array([-0.88246708,  5.47299557])
        e2.velocity = np.array([ 2.31053556, -2.21321045])
        e2.angle = 0.785398163397
        e2.angvel =   0.0320244428061
        e2.axes = (1.41421356237, 0.707106781187)
        
        self.system.particles = [e1,e2]

        tpw_ct = e1.tpw_collision_time3(e2, T=2)
        #print(">> tpw=",tpw_ct) # yields None in interval [0,2] which is correct
        self.assertTrue( numeq(tpw_ct, np.array([ 0.07169157,  0.2070478,   5.3672566 ])).all() )
        #self.assertIsNone(tpw_ct) # FIXME: this is correct
        #from bedsim.particle import ConvergenceError
        #self.assertRaises(ConvergenceError, e1.tpw_collision_time, e2, 2) # should be None... but better than nothing :P
        
        # manual testing
        tc = 0
        e1.position += e1.velocity * tc
        e1.angle += e1.angvel * tc
        e2.position += e2.velocity * tc
        e2.angle += e2.angvel * tc
        
        p = Plot(self.system)
        p.system_sync()
        p.plot_system('bedsim/data/test/EllipseCollisionTest11.eps')


    def test_collision_case12(self):
        [e1, e2] = self.system._particles
        
        e1.position = np.array([ 3.49566274,  0.99799715])
        e1.velocity = np.array([ 0.01297014, -0.09947678])
        e1.angle = 0.890346390659
        e1.angvel = 2.00353584
        e1.axes = (1.41421356237, 0.707106781187)
         
        e2.position = np.array([ 3.65316448, -0.69471609])
        e2.velocity = np.array([ 0.32197848, -1.54493718])
        e2.angle = 0.58904859881
        e2.angvel = -2.84414048618
        e2.axes = (1.41421356237, 0.707106781187)
        
        self.system.particles = [e1,e2]


        #tpw_ct = e1.tpw_collision_time3(e2, T=2)
        #print(">> tpw=",tpw_ct)
        
        from bedsim.particle import ConvergenceError
        try:
            self.assertRaises(ConvergenceError, e1.tpw_collision_time3(e2, T=2))
        except ConvergenceError as e:
            ((ce,),) = e.args # FIXME: this needs apaption once the arg passing is fixed!
            self.assertTrue( ( ce-0.008919317881258923) < 1e-4 )
        else:
            self.assertTrue(False) 
        
        #print(">> tpw=",tpw_ct) # yields None in interval [0,2] which is correct
        #self.assertTrue( numeq(tpw_ct, np.array([ 0.07169157,  0.2070478,   5.3672566 ])).all() )
        #self.assertIsNone(tpw_ct) # FIXME: this is correct
        #from bedsim.particle import ConvergenceError
        #self.assertRaises(ConvergenceError, e1.tpw_collision_time, e2, 2) # should be None... but better than nothing :P
        
        # manual testing
        tc = 0
        e1.position += e1.velocity * tc
        e1.angle += e1.angvel * tc
        e2.position += e2.velocity * tc
        e2.angle += e2.angvel * tc
        
        p = Plot(self.system)
        p.system_sync()
        p.plot_system('bedsim/data/test/EllipseCollisionTest12.eps')



    @unittest.skip("circle_approximation returns None, because particles collision happens for large t over boundary")
    def test_collision_case13(self):
        [e1, e2] = self.system._particles

        e1.position = np.array([-3.10325725,  2.99446503])
        e1.velocity = np.array([ 3.26104863, -0.2039209 ])
        e1.angle = 2.62376628575
        e1.angvel = 1.44472197702
        e1.axes     = (1.8, 1)
         
        e2.position = np.array([-0.99901586, -3.20530504])
        e2.velocity = np.array([-1.73391993,  0.47744461])
        e2.angle = 0.67156889502
        e2.angvel = 0.189738793586
        e2.axes     = (3, 0.5)
        
        self.system.particles = [e1,e2]

        tpw_ct = e1.tpw_collision_time3(e2, T=10) # ct \approx 6.04
        """
        from bedsim.particle import ConvergenceError
        try:
            self.assertRaises(ConvergenceError, e1.tpw_collision_time3(e2, T=2))
        except ConvergenceError as e:
            ((ce,),) = e.args # FIXME: this needs apaption once the arg passing is fixed!
            self.assertTrue( ( ce-0.008919317881258923) < 1e-4 )
        else:
            self.assertTrue(False)
        """
        print(">> tpw=",tpw_ct) # yields None in interval [0,2] which is correct
        #self.assertTrue( numeq(tpw_ct, np.array([ 0.07169157,  0.2070478,   5.3672566 ])).all() )
        #self.assertIsNone(tpw_ct) # FIXME: this is correct
        #from bedsim.particle import ConvergenceError
        #self.assertRaises(ConvergenceError, e1.tpw_collision_time, e2, 2) # should be None... but better than nothing :P
        
        # manual testing
        tc = 0
        e1.position += e1.velocity * tc
        e1.position = self.system.boundary.unwrap(e1.position)
        e1.angle += e1.angvel * tc
        e2.position += e2.velocity * tc
        e2.position = self.system.boundary.unwrap(e2.position)
        e2.angle += e2.angvel * tc
        
        p = Plot(self.system)
        p.system_sync()
        p.plot_system('bedsim/data/test/EllipseCollisionTest13.eps')
    
    
    def test_collision_case14(self):
        [e1, e2] = self.system._particles

        e1.position = np.array([ 9.43819383,  5.32214752])
        e1.velocity = np.array([-1.310889,    0.58227123])
        e1.angle = 0.650894016558
        e1.angvel = -2.73007848412
        e1.axes     = (1.41421356237, 0.707106781187)
         
        e2.position = np.array([ 7.85125889,  7.49351849])
        e2.velocity = np.array([ 0.63353132,  0.63763138])
        e2.angle = 0.637435228025
        e2.angvel = 0.227881904481
        e2.axes     = (1.41421356237, 0.707106781187)
        
        self.system.particles = [e1,e2]

        tpw_ct = e1.tpw_collision_time3(e2, T=10) # ct \approx 6.04
        """
        from bedsim.particle import ConvergenceError
        try:
            self.assertRaises(ConvergenceError, e1.tpw_collision_time3(e2, T=2))
        except ConvergenceError as e:
            ((ce,),) = e.args # FIXME: this needs apaption once the arg passing is fixed!
            self.assertTrue( ( ce-0.008919317881258923) < 1e-4 )
        else:
            self.assertTrue(False)
        """
        #print(">> tpw=",tpw_ct) # yields None in interval [0,2] which is correct
        self.assertTrue( numeq(tpw_ct, np.array([ 0.5745462,   7.9195222,   6.83157222])).all() )
        #self.assertIsNone(tpw_ct) # FIXME: this is correct
        #from bedsim.particle import ConvergenceError
        #self.assertRaises(ConvergenceError, e1.tpw_collision_time, e2, 2) # should be None... but better than nothing :P
        
        # manual testing
        tc = 0
        e1.position += e1.velocity * tc
        e1.position = self.system.boundary.unwrap(e1.position)
        e1.angle += e1.angvel * tc
        e2.position += e2.velocity * tc
        e2.position = self.system.boundary.unwrap(e2.position)
        e2.angle += e2.angvel * tc
        
        p = Plot(self.system)
        p.system_sync()
        p.plot_system('bedsim/data/test/EllipseCollisionTest14.eps')



    def test_collision_case15(self):
        [e1, e2] = self.system._particles

        e1.position = np.array([-3.1960663,  -0.53429538])
        e1.velocity = np.array([-0.64493742, -1.54243142])
        e1.angle = 0.815759462389
        e1.angvel = -1.40452331519
        e1.axes     = (1.8, 1)
         
        e2.position = np.array([ 2.92504462,  1.43799687])
        e2.velocity = np.array([-1.58707283, -2.27222112])
        e2.angle = 2.39079213389
        e2.angvel = 0.845875132896
        e2.axes     = (3, 0.5)
        
        self.system.particles = [e1,e2]

        tpw_ct = e1.tpw_collision_time3(e2, T=10) # ct \approx 6.04
        """
        from bedsim.particle import ConvergenceError
        try:
            self.assertRaises(ConvergenceError, e1.tpw_collision_time3(e2, T=2))
        except ConvergenceError as e:
            ((ce,),) = e.args # FIXME: this needs apaption once the arg passing is fixed!
            self.assertTrue( ( ce-0.008919317881258923) < 1e-4 )
        else:
            self.assertTrue(False)
        """
        #print(">> tpw=",tpw_ct) # yields None in interval [0,2] which is correct
        #self.assertTrue( numeq(tpw_ct, np.array([ 0.5745462,   7.9195222,   6.83157222])).all() ) # [ 4.67647932 -6.55797273 -8.95166136] vs 3.68 ### TPW FALSCH!
        self.assertTrue( numeq(tpw_ct, np.array([ 3.67595465, -4.87485368, -5.34877789])).all() )
        #self.assertIsNone(tpw_ct) # FIXME: this is correct
        #from bedsim.particle import ConvergenceError
        #self.assertRaises(ConvergenceError, e1.tpw_collision_time, e2, 2) # should be None... but better than nothing :P
        
        # manual testing
        tc = 0
        e1.position += e1.velocity * tc
        e1.position = self.system.boundary.unwrap(e1.position)
        e1.angle += e1.angvel * tc
        e2.position += e2.velocity * tc
        e2.position = self.system.boundary.unwrap(e2.position)
        e2.angle += e2.angvel * tc
        
        p = Plot(self.system)
        p.system_sync()
        p.plot_system('bedsim/data/test/EllipseCollisionTest15.eps')


    
    #@unittest.skip("Can be fixed with correction move.")
    def test_random_configuration(self):
        [e1, e2] = self.system._particles

        e1.position = np.random.uniform(-4,4,2)
        e1.velocity = np.random.normal(0,2,2)
        e1.angle    = np.random.uniform(0,2*np.pi,1)[0]
        e1.angvel   = np.random.normal(0,1,1)[0]
        e1.axes     = (1.8, 1)
        
        e2.position = np.random.uniform(-4,4,2)
        e2.velocity = np.random.normal(0,1,2)
        e2.angle    = np.random.uniform(0,2*np.pi,1)[0]
        e2.angvel   = np.random.normal(0,1,1)[0]
        e2.axes     = (3, 0.5)
        
        self.system.particles = [e1,e2]

        if e1.overlap(e2, 0):
            ref_tc = -1 # which means correction move ;)
        else: # get collision time with brute force min timesteps
            t0 = 0
            found = False
            while t0 < 10 and not found:
                t0 += 0.001
                if e1.overlap(e2, t0):
                    ref_tc = t0
                    ref_tc = round(ref_tc,2)
                    found = True
            if not found:
                ref_tc = None

        try:
            tpw_ct = e1.tpw_collision_time3(e2, T=10)
        except ValueError:
            tpw_ct = -1 # which means correction
        
        if tpw_ct is not None and tpw_ct is not -1:
            tpw_ct = tpw_ct[0] # just take time var
            tpw_ct = round(tpw_ct, 2)
        
        print("\n>>> Random test result")
        print("tpw: ", tpw_ct)
        print("ref: ", ref_tc)
        if tpw_ct == ref_tc:
            print("MATCH!")
        else:
            print("No Match!")
            print("p1: pos:",e1.position, " vel:",e1.velocity, " ang:",e1.angle, " angvel:",e1.angvel)
            print("p2: pos:",e2.position, " vel:",e2.velocity, " ang:",e2.angle, " angvel:",e2.angvel)

        if tpw_ct is None and ref_tc is not None and ref_tc > 0:
            tc = ref_tc*0.95
            e1.position += e1.velocity * tc
            e1.angle += e1.angvel * tc
            e2.position += e2.velocity * tc
            e2.angle += e2.angvel * tc
        elif ref_tc is not None and ref_tc > 0:
            tc = ref_tc
            e1.position += e1.velocity * tc
            e1.angle += e1.angvel * tc
            e2.position += e2.velocity * tc
            e2.angle += e2.angvel * tc


        p = Plot(self.system)
        p.system_sync()
        p.plot_system('bedsim/data/test/EllipseCollisionTestX.eps')