# -checked -A (should be fine but still need testing)
"""
Created on 22.01.2015

@author: mahe
"""

from collections import deque

import numpy as np 
#import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from bedsim.cubiclattice import CubicLattice #SquareLattice
from itertools import product
from decimal import Decimal

"""
Cell objects
============
"""

class Cell(object):
    """
    An abstract 'Cell' which lives in a 'Cellspace'.
    """
    
    def crossing_time(self, position, velocity):
        """
        Calculate when an object at position 'position' with velocity 'velocity' leaves the Cell.
        @param position Position of the object inside the cell
        @param velocity Velocity of the object inside the cell
        """
        raise NotImplementedError()
    
    def _find_neighbours(self):
        """
        Find all neighbours of current cell and save to self.neighbours.
        This methods needs to be called from the 'Cellspace' and NOT from
        the 'Cell' constructor, because all cells should be created before
        search for neighbours begins.
        """
        raise NotImplementedError()
    
    def __init__(self, cellspace):
        """
        Constructor for the Cell.
        @param cellspace Reference to the 'Cellspace' object in which this cell lives.
        """
        self.cellspace = cellspace
        self.neighbours = []
        self.particles = []


class DelaunayCell(Cell):
    """
    A 'Cell' object identified by one 'lattice point' which is an array of three point ids
    """

    
    def crossing_time(self, position, velocity):
        """
        Plane equation: n.(rp - r) = 0 where n is the normal vector of the plane, rp - r is a vector lying on the plane,
        and r = r0 + vt, r0 and v are the position and velocity of the particle respectively.
        There are 6 faces of the cube, the sign of the velocity will give the three faces the particle will cross.
        Cube faces are defined by the normal vectors and a point on the plane.
        Velocity signs are stored in an array velocity signs.
        Velocity faces are the equivalent faces for each velocity sign combinations
        """

        if not velocity.any(): # all components of the velocity is zero so the particle will not cross any cell
            return None

#         print("from cell", self._latticepoint,"init position",position,"init velocity",velocity)
        gp = self.cellspace._grid_points
        
        lat_len = np.linalg.norm(np.array(gp[0])-np.array(gp[1])) # smallest grid point separation
        box_len = np.amax(gp) # maximum length available if the array if flattened
        
        nn = int(box_len/lat_len) # number of sublattices axes
        
        n = np.array([[1,0,0],[1,0,0],[0,1,0],[0,1,0],[0,0,1],[0,0,1]]) # all possible normal vectors
        rp = lat_len*np.array([[0,0,0],[1,0,0],[0,0,0],[0,1,0],[0,0,0],[0,0,1]]) # all possible points in each planes of the box
            
        # For cubic lattice, here's the difference in indices corresponding to cell at the [bottom, top, left, right, back, front]
        #neighborcells = np.array([-1,+1,-Nz*Ny,Nz*Ny,-Nz,+Nz])
        
        to_cell = None
        cell_vertices = []
        
        for point in gp:
            if not any(point >= box_len): # this removes the end points so that it will satisfy periodic boundary conditions
                cell_vertices.append(point)  # contains the cell index and the coordinate of the cell origin
    
#         for j in range(3):
#             if position[j] >= box_len:
#                 position[j] = position[j] - int(position[j]/box_len)*box_len
#             else:
#                 position[j] = position[j] + int(position[j]/box_len)*box_len

        for j in range(3):
            if position[j] >= box_len:
                position[j] = position[j] - box_len
            elif position[j] < 0:
                position[j] = position[j] + box_len
    
#         gridlabel = np.array(list(map(int,position/lat_len))) # find the lattice point where the particle lies
#         print("before gridlabel",gridlabel,self._latticepoint,position,velocity)
#         if not (np.array_equal(gridlabel, self._latticepoint)):
#             print("WTF",gridlabel,self._latticepoint,position,"how to update particle cells at the boundaries?")

        pos_cell = position - np.array(self._latticepoint)*lat_len # position wrt to the cell's origin
        
        
        
        planes = []
        
        velocitysigns = np.sign(velocity)
        if velocitysigns[0] < 0.:
            planes.append(0) #left plane (i.e. the face in the -x direction)
        elif velocitysigns[0] > 0.:
            planes.append(1) #right plane (+x)
        if velocitysigns[1] < 0.:
            planes.append(2) #back plane (-y)
        elif velocitysigns[1] > 0.:
            planes.append(3) #front plane (+y)
        if velocitysigns[2] < 0.:
            planes.append(4) # bottom plane (-z)
        elif velocitysigns[2] > 0.:
            planes.append(5) # top plane (+z)
        
        crosstime = []
        crossface = []
        togridlabels = []
        
        #if planes: # if the list is not empty
        for faces in planes: 
            tt = np.dot(n[faces],rp[faces] - pos_cell) / np.dot(n[faces],velocity)
            crosstime.append(tt)
            crossface.append(faces)
            if (faces % 2 == 0): # even means to the negative direction, check velocitysigns array above
                togridlabel = np.array(self._latticepoint) - n[faces]
            else:                # odd means to the positive direction
                togridlabel = np.array(self._latticepoint) + n[faces]
                
            togridlabel[togridlabel == -1] = nn-1 # if one of the indices is negative, it means it is outside the boundary so we replace this negative value with the maximum row/column/depth
            togridlabel[togridlabel == nn] = 0    # similar concept, if one of the indices exceeds the max row/column/depth we revert them to zero
            togridlabels.append(togridlabel)
        
        
#         print("all crosstime",crosstime,"all planes",planes,n[faces],rp[faces],pos_cell)
        t = min(crosstime)
        face = crosstime.index(t)
        crossing_point = position + velocity * t
        tocell = togridlabels[face]

#         for j in range(3):
#             if crossing_point[j] >= box_len:
#                 crossing_point[j] = crossing_point[j] - int(crossing_point[j]/box_len)*box_len
#             else:
#                 crossing_point[j] = crossing_point[j] + int(crossing_point[j]/box_len)*box_len

        for j in range(3):
            if crossing_point[j] >= box_len:
                crossing_point[j] = crossing_point[j] - box_len
            elif crossing_point[j] < 0.:
                crossing_point[j] = crossing_point[j] + box_len


        
        for n in self.neighbours:
            if np.array_equal(n._latticepoint, tocell): #all(n._latticepoint == gridlabel):
                to_cell = n

        return (t,to_cell,crossing_point)

    
    def _find_neighbours(self):

        """
        Each lattice cell have 26 neighbours in total, which was obtained by adding the combinations of [-1,0,1] to the lattice points index
        """
        
        gp = self.cellspace._grid_points
        
        lattice_length = np.linalg.norm(gp[0] - gp[1])
        box_length = np.amax(gp)
        nrows = int(box_length/lattice_length)
        
        points = [-1,0,1]
        neighbours_cells = list(product(points,repeat=3)) # this will return an array, with a length of 27 including the point 000 
                                                          # which when added to the cell indices will give the neighbours 
        neighbour_index = []                                     
        for indices in neighbours_cells:
            nn = self._latticepoint + np.array(indices)
            if not np.array_equal(self._latticepoint, nn):
                nn[nn == -1] = nrows - 1
                nn[nn == nrows] = 0
                neighbour_index.append(nn)
        
#         print("inside _find_neighbours",self._latticepoint,len(self.cellspace.cells),len(neighbour_index))        
        for cell in self.cellspace.cells:
            for neigh in neighbour_index:
#                 print("find n",neigh,cell._latticepoint)
                if np.array_equal(cell._latticepoint, neigh):
                    self.neighbours.append(cell)
                    
        self.neighbours = list(set(self.neighbours)) #remove duplicates, not needed if system is big
#         print("len self.neighbours",len(self.neighbours))
    
    def __init__(self, cellspace, latticepoint):
        Cell.__init__(self, cellspace)
        self._latticepoint = latticepoint


"""
Cellspaces
==========
"""

class Cellspace(object):

    def assign_particle(self, particle):
        """
        Calculate to which "Cell" a particle should belong.
        @param particle Particle which should be assigned to a cell
        """
        raise NotImplementedError()
                  
    def _create_cells(self):
        """
        Create cells in the current 'Cellspace'.
        """
        raise NotImplementedError()
        
    def _calculate_neighbours(self):
        """
        Calculate neighbours of each cell in the 'Cellspace' and save neighbour lists to each 'Cell'.
        """
        for cell in self.cells:
            cell._find_neighbours()

    def __init__(self, system):
        self.system = system
        self.cells = []
        self._create_cells()
        self._calculate_neighbours()
        print('worked')
        


class DelaunayCellspace(Cellspace):
    
    def assign_particle(self, particle):
        pcell = self._locate_cell(particle.position)
        pcell.particles.append(particle)
        particle.cell = pcell
#        print(particle._id,pcell._latticepoint,particle.position)
    
    def reassign_particle(self, particle):
        pcell = self._locate_cell(particle.position)
        if pcell is not particle.cell: # only do this if we are not in the correct cell
            particle.cell.particles = list(filter(lambda x: x is not particle, particle.cell.particles)) # remove particle from old cell
            pcell.particles.append(particle)
            particle.cell = pcell
    
    def _create_cells(self):
        self.__grid_to_triang()
        self.__triang_to_cells()

    def _locate_cell(self, position):
        """
        Determine in which cell the coordinates 'position' reside in and return the reference to that 'Cell'.
        This Method doesn't use the reference lists of 'Particle' and 'Cell'!
        """

        ps = self.triang.find_gridpoint(position) # just get the index of the cell to find the equivalent simplex
        for c in self.cells:
            if np.array_equal(c._latticepoint,ps):
#             if (c._latticepoint == ps).all(): # if they have the same grid points
#             if (c._simplex == self.triang.simplices[ps]).all():

                return c
    
    def __grid_to_triang(self):
        """
        This makes the cubic Lattice
        Gives indices of each point and the 8 coordinates of each vertices in the cube
        """
        self.triang = CubicLattice(self._grid_points) # generates an array of triangulation objects
        #self.triang = SquareLattice(self._grid_points)
        
    def __triang_to_cells(self):
        """
        Create cells from triangulation self.triang
        """
        for s in self.triang.latticepoints:
            self.cells.append(DelaunayCell(self, s))
    
    
    def __init__(self, system, grid_points):
        """
        The 'Cellspace' constructor already calls _create_cells(). 
        @param system Reference to a System object.
        @param grid_points List of points which should be used to construct the Cellspace. IMPORTANT: The grid needs translational symmetry to fulfill periodic boundary conditions!
        """

        self._grid_points = grid_points
        self.triang = None
        Cellspace.__init__(self, system)
        print('delaunay')