import unittest
import numpy as np
from collections import deque

from bedsim.boundary import PeriodicBox, has_identical_coords, list_duplicates, has_shifted_coords, equivalent_indices, list_shifted_duplicates
from bedsim.cell import DelaunayCellspace


def numeq(x, y):
    """
    Function to determine approx equality of two numeric variables
    """
    epsilon = 10**-10
    return np.abs(x-y) < epsilon


class arb(object):
    """
    Arbitrary class for testing purposes.
    """
    pass


class TestObsoleteHelers(unittest.TestCase):
    """
    Test class for unused helpers.
    """
    def setUp(self):
        pass

    def test_id_coords(self):
        self.assertEqual(has_identical_coords([[1,0], [0,0]], [[1,0], [0,0]]), True)
        self.assertEqual(has_identical_coords([[1,0], [0,0]], [[0,0], [1,0]]), True)
        self.assertEqual(has_identical_coords([[1,0], [0,0]], [[0,1], [0,0]]), False)


class TestHelpers(unittest.TestCase):
    """
    Test class for helpers.
    """
    def setUp(self):
        pass

    def test_dups(self):
        self.assertEqual(list_duplicates([1,2,3,1,5,2,1,3,4,8,3,1], 1), [0,3,6,11])
        self.assertEqual(list_duplicates([1,2,3], 4), [])
    
    def test_shifted_coords(self):
        self.assertEqual(has_shifted_coords([[0,0],[1,0]], [[0,2],[1,2]]), True)  # test basic shifting
        self.assertEqual(has_shifted_coords([[1,0],[0,0]], [[0,2],[1,2]]), True)  # test if function is symmetric as desired
        self.assertEqual(has_shifted_coords([[0,0],[1,0]], [[0,5],[1,5]]), True)  # basic shifting over 5 units
        self.assertEqual(has_shifted_coords([[1,0],[0,0]], [[0,1],[0,0]]), False) # should fail because vertices are orthogonal 
        self.assertEqual(has_shifted_coords([[0,0],[1,0]], [[0,0],[0,1]]), False) # should fail because vertices are orthogonal
        self.assertEqual(has_shifted_coords([[1,0],[0,0]], [[1,2],[2,2]]), False) # should fail because of shift in 2 directions 
        self.assertEqual(has_shifted_coords([[1,0],[2,0]], [[1,2],[3,2]]), False) # should fail because of different lengths
        self.assertEqual(has_shifted_coords([[0,0],[1,0]], [[2,0],[3,0]]), False) # should fail because only shifted along its axis
        self.assertEqual(has_shifted_coords([[0,0],[1,0]], [[3,0],[4,0]]), False) # should fail because only shifted along its axis

    def test_list_shifted_duplicates(self):
        haystack = [[x,0] for x in range(6)] + [[5,y+1] for y in range(5)] + [[4-x,5] for x in range(5)] + [[0,4-y] for y in range(4)] # line segments around a 5x5 cellspace
        
        h2 = deque(haystack)
        h2.rotate(1)
        h2 = list(h2)
        haystack = [[h1, h2] for (h1,h2) in zip(haystack, h2)]
        #print("HAYSTACK=", haystack)
        
        needle = [[0,0],[1,0]]
        self.assertEqual(list_shifted_duplicates(haystack, needle), [1,15])
        #print("SD: ",list_shifted_duplicates(haystack, needle))
        #print("MANUAL MODE")
        #w = [has_shifted_coords(e,needle) for e in haystack]
        #print("w=", w)

    def test_equivalent_indices(self):
        pass



class TestPeriodicBox(unittest.TestCase):
    """
    Test class for the periodic box.
    """
    def setUp(self):
        #print("\n\n>>> TPB: ")
        system = arb()
        system.system_properties = arb()
        #system.system_properties.boxSize = np.array([2,2])
        
        grid_points = np.array([[0,0], [0,1], [1,1], [1,0], [0,2], [1,2], [2,2], [2,1], [2,0]])
        system.cellspace = DelaunayCellspace(system, grid_points)
        
        self.system = system
        self.boundary = PeriodicBox(system)


    def test_unwrap(self):
        self.assertEqual((self.boundary.unwrap([3,1]) == np.array([1,1])).any(), True)
        self.assertEqual((self.boundary.unwrap([1,3]) == np.array([1,1])).any(), True)
        self.assertEqual((self.boundary.unwrap([-1,1]) == np.array([1,1])).any(), True)
        self.assertEqual((self.boundary.unwrap([0,0]) == np.array([0,0])).any(), True)
        self.assertEqual((self.boundary.unwrap([0,1]) == np.array([0,1])).any(), True)

    def test_distance(self):
        self.assertEqual((self.boundary.delta(np.array([0,0]), np.array([0.1, 0.1])) == np.array([0.1, 0.1])).all(), True)
        self.assertEqual(numeq(self.boundary.delta(np.array([0.2,1]), np.array([1.8, 1])), np.array([0.4, 0.])).all(), True)
        self.assertEqual(numeq(self.boundary.delta(np.array([1.8,0.3]), np.array([0.2, 0.2])), np.array([0.4, 0.1])).all(), True)
        self.assertEqual(numeq(self.boundary.delta(np.array([0.2,0.2]), np.array([1.8, 0.3])), np.array([0.4, 0.1])).all(), True)
    
    def test_direction(self):
        self.assertEqual(numeq(self.boundary.delta_dir(np.array([1.8,0.3]), np.array([0.2, 0.2])), np.array([0.4, -0.1])).all(), True)
        self.assertEqual(numeq(self.boundary.delta_dir(np.array([0.2,0.2]), np.array([1.8, 0.3])), np.array([-0.4, 0.1])).all(), True)
        self.assertEqual(numeq(self.boundary.delta_dir(np.array([0.2,0.2]), np.array([0.1, 0.1])), np.array([-0.1, -0.1])).all(), True)
        self.assertEqual(numeq(self.boundary.delta_dir(np.array([0.1,0.1]), np.array([0.2, 0.2])), np.array([0.1, 0.1])).all(), True)
        #self.assertEqual(numeq(self.boundary.delta_dir(np.array([0.2,0.4]), np.array([0.2, 1.4])), np.array([0, 1])).all(), True) # question of definition on a 2x2 grid => may be true or false
        #self.assertEqual(numeq(self.boundary.delta_dir(np.array([0.2,0.4]), np.array([1.2, 0.4])), np.array([-1, 0])).all(), True) # question of definition on a 2x2 grid
        self.assertEqual(numeq(self.boundary.delta_dir(np.array([0.2,0.5]), np.array([0.2, 1.4])), np.array([0, 0.9])).all(), True) # this is unique

    def test_delta_wrap_identity(self):
        """
        Test if unwrap of accumulated "delta_dir"s and single "unwrap"s yield the same result.
        """
        steps = 1000
        moves = np.random.uniform(-1,1,(steps,2)) # perform random moves inside the 2x2 box
        r0 = np.array([1,1]) # start in the middle of the box
        #acc = np.array([0,0]) # accumulated movement
        acc = np.array([1,1]) # is not anymore the accumulated movement but also considers starting point which makes comparison at the end easier (+[1,1] not needed)
        
        for move in moves:
            r1 = self.boundary.unwrap(r0 + move) # update the position
            acc = acc + self.boundary.delta_dir(r0, r1) # recalculate distance from r0 to r1
            r0 = r1

        self.assertTrue(numeq(r0, self.boundary.unwrap(acc)).all()) # check if single unwraps did the same thing as unwrapping the accumulator


    def test_edge_assignment(self):
        # TODO!
        for c in self.system.cellspace.cells:
            #print(len(c.edges))
            for e in c.edges:
                (c1, c2) = e.cells
                (p1c1, p1c2) = e.point1 # c1 points are points of this cell
                (p2c1, p2c2) = e.point2
                #print("From cell %s to cell %s" % (c1, c2))
                #print("Cell 1 Edge: (%s, %s)" % (p1c1, p2c1))
                #print("Cell 2 Edge: (%s, %s)" % (p2c2, p2c2))
                pass


class TestBigPeriodicBox(unittest.TestCase):
    """
    Test class for a bigger periodic box.
    """
    def setUp(self):
        system = arb()
        system.system_properties = arb()
        
        n = 6 # issue with range(n) => range(6) means 0,1,2,3,4,5
        grid_points = np.array([np.array([x,y]) for x in range(n) for y in range(n)])
        system.cellspace = DelaunayCellspace(system, grid_points)
        
        self.system = system
        self.boundary = PeriodicBox(system)

    def test_edge_assignment(self):
        pid_to_coord = lambda x: self.system.cellspace._grid_points[x]
        flatten = lambda list_of_lists: [val for sublist in list_of_lists for val in sublist]
        
        ###print("#Border egdes: ", len(self.boundary._bordervertices)) # should be 20
        
        #print("BV coord::")
        bordervertices_coord = [[pid_to_coord(x).tolist() for x in bv] for bv in self.boundary._bordervertices]
        ###print(bordervertices_coord)
        
        ###print("eq idx: ", equivalent_indices(bordervertices_coord))
        ###print(self.boundary._eq_vertices)
        
        ###print("# Edges per cell")
        for c in self.system.cellspace.cells:
            ###print(len(c.edges))
            pass