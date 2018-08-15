import numpy as  np
from itertools import product

"""
Attempt to imitate the Delaunay tesselation done in the code -A
"""


class CubicLattice(object):
    """ 
                           Cell indices are given by: (r,c,d) which are row, columns and depth indices in x,y,z respectively
    z(top)                 i.e. if box is divided into 8 smaller cubes, the small cubes are labeled as follows:
    |   y(front)           left-back-bottom 0    left-front-bottom 2   right-back-bottom 4     right-front-bottom 6     
    | /                    left-back-top 1       left-front-top 3      right-back-top 5        right-front-top 7
    |/______x(right)       Each smaller box has 8 vertices/simplices 6 faces
                           
    """

    def __init__(self, gridpoints):
        self.grid_points = gridpoints
        self.lattice_length = np.linalg.norm(self.grid_points[0] - self.grid_points[1])
        self.box_length = np.amax(self.grid_points)
        self.nrows = int(
            self.box_length / self.lattice_length)  # number of rows in x y and z # change this to a non cubic box
        print("number of rows",self.nrows,self.lattice_length)
        self.latticepoints = self.all_gridpoints()  # change this to cell indices ie compute_indices
        #print("inside cubiclattice self.latticepoints",self.latticepoints) ## SOPHIA Juni
        # self.neigbours = self.compute_neighbour_cells()

    def cell_coordinate(self, position):

        for j in [0, 1,
                  2]:  # use boundary file maybe, find a way to remove this like unwrap should be implemented already
            if position[j] >= self.box_length:
                position[j] = position[j] - int(position[j] / self.box_length) * self.box_length
            else:
                position[j] = position[j] + int(position[j] / self.box_length) * self.box_length

        gridlabel = np.array(list(map(int, position / self.lattice_length)))

        return position - gridlabel * self.lattice_length

    def all_gridpoints(self):
        """
        This will return a list of indices of each lattice cubes found in the simulation box
        """
        grid_space = np.linspace(start=0, stop=self.nrows, num=self.nrows, endpoint=False).tolist()
        grid_points = [[int(i), int(j), int(k)] for i in grid_space for j in grid_space for k in grid_space]

        return grid_points

    def find_gridpoint(self, position):

        lat_len = self.lattice_length
        box_len = self.box_length

        # find a way to remove this
        for j in [0, 1, 2]:
            if position[j] >= box_len:
                position[j] = position[j] - int(position[j] / box_len) * box_len
            else:
                position[j] = position[j] + int(position[j] / box_len) * box_len

        grid_point = np.array(list(map(int, position / lat_len)))
        return grid_point
