# Checked -A
import numpy as np

# from bedsim.cell import Edge

"""
boundary module
===============
This module should handle boundary conditions of the cellspace transparently.

IMPORTANT INFORMATION:
Throughout this module 'simplex' refers to just a line segment consisting of exactly two points!
"""


"""
Obsolete functions... but may be useful later :P
"""

#boundary.eqassign([1,2,3,4,5,6,7,8], lambda x,y: x % 5 == y)
# def eqassign(elements, eqfun, assign = {}):
#     """
#     Assign 'elements' to each other by an equality function 'eqfun'.
#     @param elements List of elements. Must be hashable.
#     @param eqfun Function to determine equality of the elements.
#     @param assign Dict to store results for iteration. Don't touch! 
#     
#     If not all elements can be assigned no error is raised.
#     Function just terminates and returns assigned elements.
#     """
#     x = elements.pop()
#     for y in elements:
#         if eqfun(x,y):
#             elements.remove(y)
#             assign[x] = y
#             assign[y] = x
#     
#     if elements != []:
#         return eqassign(elements, eqfun, assign)
#     else:
#         return assign
# 
# 
# def has_identical_coords(vertex1, vertex2):
#     """
#     Tells whether two vertices have the same coordinates.
#     BUT does not say anything about having the same line order!
#     @param vertex1,vertex2 vertex = [[x1,y2], [x2,y2], ...], where x_i, y_i are real coordinates
#     """
#     return all([(v1 in vertex2) for v1 in vertex1])
# 
# 
# """
# Helper functions
# ================
# """
# 
# def list_duplicates(elements, item):
#     """
#     Return a list of indices with all occurrences of 'item' in elements
#     @param elements List of elements, i.e. haystack
#     @param item Item to search for in elements, i.e. needle
#     """
#     start_at = -1
#     locations = []
#     while True:
#         try:
#             location = elements.index(item,start_at+1)
#         except ValueError:
#             break
#         else:
#             locations.append(location)
#             start_at = location
#     return locations
# 
# 
# def has_shifted_coords(vertex1, vertex2):
#     """
#     Check if one vertex is a shifted version of the other vertex.
#     In this context 'shifted' means, that vertex1 and vertex2 have the same coordinates along one axis but
#     distinct coordinates along the other axis. If all coordinates are the same, then 'vertex1==vertex2' which
#     is also considered a valid shift (with offset 0). 
#     @param vertex1,vertex2 vertex = [[x1,y2], [x2,y2], ...], where x_i, y_i are real coordinates
#     """
#     (vertex1, vertex2) = (np.array(v) for v in [vertex1, vertex2]) # make sure that vertices are represented as numpy arrays
#     
#     """ 0. calculate the segment vectors of the vertices """ 
#     aprime = vertex1[1] - vertex1[0]
#     bprime = vertex2[1] - vertex2[0]
# 
#     """ 1. Check that segments have the same length. If not they are not just shifted. """
#     if np.linalg.norm(aprime) != np.linalg.norm(bprime): # i.e. segments of different lengths
#         return False
#     
#     """
#     2. If aprime and bprime are parallel they need to have either the same signs or -1 times the other's signs.
#     NOTE: This does NOT guarantee that aprime and bprime are parallel, but we will find out later if the scalar
#     product between difference vectors and a prime vectors is not exactly 0.
#     """
#     if (np.sign(aprime) - np.sign(bprime)).any() != 0: 
#         vertex2 = np.array([vertex2[1], vertex2[0]]) # i.e. rotate
#         bprime = np.sign(aprime) * bprime            # i.e. switch signs according to the other segment vector 
# 
#     """
#     3. Check if vertex1 is parallel to vertex2.
#     For most cases 2. + 4. are sufficient but there is a special case ([[0,0],[1,0]] and [[0,0],[0,1]], see testcase)
#     where manual checking is needed.
#     """
#     if np.dot(aprime, bprime) != np.linalg.norm(aprime)**2: # **2 because norms were already checked for equality
#         return False
# 
#     """
#     4. vertices are properly shifted if difference vectors of the vertices and appropriate segment vectors are orthogonal.
#     """
#     [d1, d2] = vertex2 - vertex1 # difference vectors
#     if np.dot(aprime, d1) == 0 and np.dot(bprime, d2) == 0: # vertices are shifted if appropriate segment vectors and difference vectors are orthogonal
#         return True
#     else:
#         return False
#     
# 
# # UNUSED: maybe remove later
# def has_shifted_coords2(vertex1, vertex2):
#     """
#     Check if one vertex is a shifted version of the other vertex.
#     In this context 'shifted' means, that vertex1 and vertex2 have the same coordinates along one axis but
#     distinct coordinates along the other axis. If all coordinates are the same, then 'vertex1==vertex2' which
#     is also considered a valid shift (with offset 0). 
#     @param vertex1,vertex2 vertex = [[x1,y2], [x2,y2], ...], where x_i, y_i are real coordinates
#     """
#     [v1x, v1y] = [list(x) for x in zip(*vertex1)] # Transpose
#     [v2x, v2y] = [list(x) for x in zip(*vertex2)] # Transpose
#     if (set(v1x) == set(v2x) and set(v1y).isdisjoint(set(v2y))) or (set(v1y) == set(v2y) and set(v1x).isdisjoint(set(v2x))):
#         return True
#     elif set(v1x) == set(v2x) and set(v1y) == set(v2y): # "vertex1 == vertex2"
#         return True # IS THIS OK?
#     else:
#         return False
#  
# 
# def list_shifted_duplicates(elements, item):
#     """
#     Determine all shifted duplicates of 'item' in 'elements'
#     @param elements Haystack of border line segments
#     @param item Needle line segment whose shifted duplicates should be searched for.
#     """
#     w = [has_shifted_coords(e,item) for e in elements]
#     return list_duplicates(w, True)
# 
# 
# def equivalent_indices(elements):
#     """
#     Assume elements is a list with an even amount of items.
#     Further each item  is assumed to be equivalent to exactly one other item (no duplicates!).
#     @param elements List of items where each item is exactly equivalent to another item.
#     """
#     idx=[]
#     for e in elements:
#         idx.append(list_shifted_duplicates(elements, e))
#     return idx




"""
Boundary classes
================
"""

class Boundary(object):
    """
    General Boundary class
    """
    
    def unwrap(self, x):
        raise NotImplementedError()
    
    def delta(self, x, y):
        raise NotImplementedError()
    
    def delta_dir(self, x, y):
        raise NotImplementedError()
    
#     def to_space(self, x, y):
#         raise NotImplementedError()
#     
#     def to_body(self, x, y):
#         raise NotImplementedError()

    def __init__(self, system):
        self.system = system


class PeriodicBox(Boundary):
    """
    Boundary conditions for a periodic box.
    """
    
    def unwrap(self, x):
        """
        Unwrap a vector x inside the periodic box.
        @param x Position vector.
        """
        (xmin, xmax, ymin, ymax, zmin, zmax) = self._box_corners
        # Perform transformations, so that unwrap also works for boxes where (xmin, ymin) != (0,0) 
#         return np.array([xmin, ymin, zmin]) + np.remainder(np.array(x)-np.array([xmin, ymin, zmin]), self.__box_size)
        size = np.absolute([xmax - xmin, ymax - ymin, zmax - zmin])

        for j in range(3):
            if x[j] >= size[j]:
                x[j] = x[j] - size[j]
            elif x[j] < 0.:
                x[j] = x[j] + size[j]

        return x
                
#         newx = x - np.sign(x)*size
#         
        #return np.remainder(x, self.__box_size) # Non transformed version

    def delta(self, x, y):
        """
        Absolute difference vector between the vectors x and y considering boundary conditions.
        @param x,y Position vector. 
        """
        (xmin, xmax, ymin, ymax, zmin, zmax) = self._box_corners
        size = np.absolute([xmax - xmin, ymax - ymin, zmax - zmin])
        (x, y) = (np.array(x), np.array(y))

        delta = np.absolute(x - y)
        delta = np.where(delta > 0.5*size, size - delta, delta)
        return delta


    def delta_dir(self, x, y):

        """
        Difference vector between the vectors x and y considering boundary conditions.
        The direction of the vector points from x to y, i.e. signs of the components are preserved.   
        Origin is at one of the vertices
        @param x,y Position vector.
        """ 
        (xmin, xmax, ymin, ymax, zmin, zmax) = self._box_corners
        size = np.absolute([xmax - xmin, ymax - ymin, zmax - zmin])
        z = y-x
        
        delta = np.float64(np.where(abs(z) > 0.5*size, z-np.sign(z)*size, z)) #SOPHIA: np.where(if cond true, do x, otherwise do y) element wise

        return delta

    def delta_norm(self, x, y):
        """
        Absolute difference vector between the vectors x and y considering boundary conditions.
        @param x,y Position vector.
        """
        (xmin, xmax, ymin, ymax, zmin, zmax) = self._box_corners
        size = np.float64(np.absolute([xmax - xmin, ymax - ymin, zmax - zmin]))
        (x, y) = (np.array(x), np.array(y))

        delta = np.float64(np.absolute(x - y))  # component wise distances from two considered particles x,y
        delta = np.where(delta > 0.5 * size, size - delta, delta)  # set to minimum distance

        return np.float64(np.linalg.norm(delta))
    
#     def to_space(self, x, y):
#         """
#         Transform from vector vec to rotated coordinates in quaternion form q, from body 
#         coordinates to space coordinates 
#         """
#         
#         x1, y1, z1, w1 = y
#       
#         p = np.array([[x1,y1,z1]])
#         P = np.array([[0,-z1,y1],[z1,0,-x1],[-y1,x1,0]],dtype=np.float64)
#         Q= 2.0*(np.matmul(p.transpose(),p) - w1*P + (w1*w1 - 0.5)*np.identity(3, dtype=np.float64))
#       
#         vecspace = np.matmul(Q.transpose(),x)
#         
#         return vecspace
#         
#     def to_body(self, x, y):
#         """
#         Transform from vector vec to rotated coordinates in quaternion form q, from space 
#         coordinates to body coordinates
#         """
#         
#         x1, y1, z1, w1 = y
#       
#         p = np.array([[x1,y1,z1]])
#         P = np.array([[0,-z1,y1],[z1,0,-x1],[-y1,x1,0]],dtype=np.float64)
#         Q= 2.0*(np.matmul(p.transpose(),p) - w1*P + (w1*w1 - 0.5)*np.identity(3, dtype=np.float64))
#       
#         vecbody = np.matmul(Q,x)
# 
#         return vecbody

#     def __boundary_simplex_to_cell(self, boundary_simplex):
#         """
#         Search in the cellspace to which cell the 'boundary_simplex' line
#         segment belongs. Only works if boundary_simplex is on the cellspace
#         boundary, else the assignment would not be unique.
#         @param boundary_simplex A line segment one the cellspace boundary in list format.
#         """
#         for cell in self.system.cellspace.cells:
#             if set(boundary_simplex).issubset(set(cell._simplex.tolist())):
#                 return cell
#     
#     
#     def __find_boundary_neighbours(self):
#         """
#         Idea:
#         A normal neighbour cell has at least one point in common.
#         In case of a boundary neighbour cell those points are not identical but have delta_dir = 0.
#         NOTE: this algorithm is not optimal and not elegant, but very straight forward and sufficient for now.
#         """
#         flatten = lambda list_of_lists: [val for sublist in list_of_lists for val in sublist]
#         
#         # get point ids of cells which have at least one point lying on the boundary
#         border_point_ids = flatten(self.system.cellspace.triang.convex_hull.tolist())
#         
#         # get list of cells which have at least one point lying on the boundary
#         bordercells = []
#         for point_id in border_point_ids:
#             for cell in self.system.cellspace.cells:
#                 if point_id in cell._simplex:
#                     bordercells.append(cell)
#         bordercells = list(set(bordercells)) 
#         
#         # check each border cell with all other border cells 
#         for bordercell1 in bordercells:
#             for bordercell2 in bordercells:
#                 if bordercell1 is not bordercell2:
#                     for point_id1 in bordercell1._simplex: # check all points in the cell simplex with all points of the other cell's simplex
#                         for point_id2 in bordercell2._simplex:
#                             x = self.system.cellspace._grid_points[point_id1]
#                             y = self.system.cellspace._grid_points[point_id2]
#                             # if distance is 0 (i.e. numerically small), the cells are neighbours
#                             if np.linalg.norm(self.delta_dir(x, y)) < 1e-4: # FIXME: calc epsilon from cell size!! 
#                                 bordercell1.neighbours.append(bordercell2)
#                                 bordercell2.neighbours.append(bordercell1)
#         
#         for bordercell in bordercells: # remove duplicates
#             bordercell.neighbours = list(set(bordercell.neighbours))
#             
#     
#     def __boundary_edge_assignment(self): # Tested and working
#         """
#         Assign to each boundary cell a corresponding boundary neighbour cell.
#         How the algorithm works:
#         1. Determine boundary segments
#         2. Check which segments belong together, i.e. those on opposite sites of the box at same height (roughly speaking)
#         3. Take equivalent segments and determine corresponding cells
#         4. Write the additional neighbour information to the cells 
#         5. Save border vertices and which of them are equivalent for further use
#         """
#         # Calculate point coordinates from point ID
#         #pid_to_coord = lambda x: self.system.cellspace._grid_points[x]
#         pid_to_coord = self.system.cellspace._grid_points
#         
#         # 1. Find all border vertices (NOTE: here vertex means line segment consisting of exactly 2 points)
#         bordervertices = self.system.cellspace.triang.convex_hull.tolist()
#         
#         # 2. convert bordervertices to coordinates
#         #bordervertices_coord = [[pid_to_coord(x).tolist() for x in bv] for bv in bordervertices]
#         bordervertices_coord = [[pid_to_coord[x].tolist() for x in bv] for bv in bordervertices]
#         
#         # 3. Check which vertices are equal. here 'equal' means that either every x or y coords must be pairwise identical 
#         eq_indices = equivalent_indices(bordervertices_coord)
# 
#         # 4. construct associative arrays of equivalent vertices
#         eq_vertices = {}
#         for i in eq_indices: # REMARK: traverses everything twice due to eq_indices has already considered this but without proper element order...
#             if len(i) == 2: # FIXME
#                 [i1, i2] = i # NOTE: i1 and i2 are vertex indices, not vertex IDs!
#                 # save eq_vertices to an internal index for later use 
#                 eq_vertices[i1] = i2
#                 eq_vertices[i2] = i1
#                 
#                 # determine cells from the border vertices
#                 c1 = self.__boundary_simplex_to_cell(bordervertices[i1])
#                 c2 = self.__boundary_simplex_to_cell(bordervertices[i2])
#     
#                 # save equivalent edges to cells' edges list
#                 [point1, point2] = np.transpose([bordervertices[i1], bordervertices[i2]]).tolist() #[point1, point2] = [(p1c1, p1c2), (p2c1, p2c2)]
#                 e1 = Edge((c1, c2), tuple(point1), tuple(point2))
#                 e2 = Edge((c2, c1), tuple(point1)[::-1], tuple(point2)[::-1])
#                 #e1 = Edge((c1, c2), tuple(bordervertices[i1]), tuple(bordervertices[i2])) # FIXME: TRANSPOSE!!!
#                 #e2 = Edge((c2, c1), tuple(bordervertices[i2]), tuple(bordervertices[i1])) # FIXME: TRANSPOSE!!!
#                 c1.edges.append(e1)
#                 c2.edges.append(e2)
# 
#         # save raw data for possible later use
#         self._bordervertices = bordervertices
#         self._bordervertices_coord = bordervertices_coord
#         self._eq_vertices = eq_vertices
#         
#         """       
#         print("\nBOUNDARY EDGE ASSIGNMENT\n========================")
#         print("bordervertices: ", bordervertices)
#         print("bordervertices_coord: ", bordervertices_coord)
#         print("eq idx: ", equivalent_indices(bordervertices_coord))
#         print("eq vertices: ", eq_vertices)
#         """

    
    def __grid_to_box_corners(self):
        """
        Calculate the corners of the simulation box from the cell grid.
        """
        x,y,z = self.system.cellspace._grid_points.transpose()
        [xmin, xmax, ymin, ymax, zmin, zmax] = [np.amin(x), np.amax(x), np.amin(y), np.amax(y), np.amin(z), np.amax(z)]
        [width, height, depth] = [xmax-xmin, ymax-ymin, zmax - zmin]
        self.__box_size = np.array([width, height, depth])
        self.__box_center = np.array([(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2])
        self._box_corners = (xmin, xmax, ymin, ymax, zmin, zmax)
        print('box corners', self._box_corners)
        self.extent = np.fabs(np.array([xmax, ymax, zmax]) - np.array([xmin, ymin, zmin])) 
    
    def __init__(self, system):
        Boundary.__init__(self, system)
        self.__grid_to_box_corners()
#         self.__find_boundary_neighbours() # assign bounray neighbours
#         self.__boundary_edge_assignment() # assign boundary edges