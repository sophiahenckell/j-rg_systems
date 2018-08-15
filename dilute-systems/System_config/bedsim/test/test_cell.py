import unittest
import numpy as np
import matplotlib.pyplot as plt

#import bedsim.cell
from bedsim import cell


class arb(object):
    """
    Arbitrary class for testing purposes.
    """
    pass


class TestDelaunayCellspace(unittest.TestCase):
    
    def setUp(self):
        self.points = np.array([[0,0], [0,1], [1,1], [1,0], [0,2], [1,2], [2,2], [2,1], [2,0]])
        self.dcs = cell.DelaunayCellspace(None, self.points)
        
        """
        for c in self.dcs.cells:
            print(c._simplex)
            print("has neighbours")
            print([list(x._simplex) for x in c.neighbours])
            #print(list(c._simplex) + " has neighbours: " + list(c.neighbours._simplex))
        """

    @unittest.skip("Activate if you want to see an example plot.")
    def test_plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.triplot(self.points[:,0], self.points[:,1], self.dcs.triang.simplices.copy())
        for xy in list(self.points):
            xy2 = tuple(list(xy))
            ax.annotate('(%s, %s)' % xy2, xy=xy2, textcoords='offset points')
        ax.set_ylim(-1,3)
        ax.set_xlim(-1,3)
        plt.grid()
        plt.savefig('bedsim/data/test/dcs.eps')

    def test_triangulation(self):
        # Triangulation should yield the following vertex configuration
        eqvertices = [[5, 1, 2], [1, 5, 4], [1, 3, 2], [3, 1, 0], [7, 5, 2], [5, 7, 6], [3, 7, 2], [7, 3, 8]]
        self.assertEqual(self.dcs.triang.vertices.tolist(), eqvertices)
    
    def test_assign_particle(self):
        particle = arb()
        particle.position = np.array([0.2, 0.2])
        self.dcs.assign_particle(particle)
        # FIXME: test which vertex, etc

    def test_edge_assignment(self):
        # TODO
        pass


    @unittest.skip("Activate if you want to see cell assignment with a nice plot.")
    def test_assign_particle_plt(self):
        pp = np.array([1.2,1.2]) # particle position
        ps = self.dcs.triang.find_simplex(pp)
        #print(ps)

        print(self.dcs.triang.find_simplex([2,2]))
        print(self.dcs.triang.find_simplex([1,2]))

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax2 = fig.add_subplot(111)
        ax.triplot(self.points[:,0], self.points[:,1], self.dcs.triang.simplices.copy())
        ax2.plot(self.points[self.dcs.triang.simplices[ps]][:,0],self.points[self.dcs.triang.simplices[ps]][:,1])
        #print(self.points[self.dcs.triang.simplices[ps]])
        for xy in list(self.points):
            xy2 = tuple(list(xy))
            ax.annotate('(%s, %s)' % xy2, xy=xy2, textcoords='offset points')
        ax.set_ylim(-1,3)
        ax.set_xlim(-1,3)
        plt.grid()
        plt.savefig('bedsim/data/test/dcs2.eps')
    
if __name__ == '__main__':
    unittest.main()
